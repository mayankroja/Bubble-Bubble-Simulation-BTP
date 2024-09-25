#include "Interface.H"
#include <AdvectLS_helper.H>
#include <WeightedENO.h>
#include <AMReX_ParmParse.H>
#include <algoim_hocp.hpp>

namespace mycode
{


//* NOTE : as the psi is computed on grown box
//*      : do not do phaseFieldBoundary here, this was the difference in 
//*      : earlier multiphase code(Algoim::Reinit in interface.cpp)
//*      : The phaseFieldBoundary call from this was causing the code to fail
//*      : due to overwriting the computed values with communicated values
void Interface::AlgoimReinit()
{
    BL_PROFILE("Interface::AlgoimReinit()");
    //PhaseFieldBC();    
    const amrex::Real *dx = mesh_.geometry().CellSize();
    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const &f = psi.array(mfi);

        blitz::Array<amrex::Real, AMREX_SPACEDIM> phi(AMREX_D_DECL(bx.length(0), bx.length(1), bx.length(2)));

        /// copy Psi in phi
        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
        {
            int ilocal = i - bx.smallEnd(0);
            int jlocal = j - bx.smallEnd(1);
            phi(ilocal, jlocal) = f(i, j, k);
        });

        /// reinit
        Algoim::reinit<AMREX_SPACEDIM, 3>(phi, dx[0], 1.8 * LAYERS * dx[0]);

        /// copy phi in psi
        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
        {
            int ilocal = i - bx.smallEnd(0);
            int jlocal = j - bx.smallEnd(1);
            f(i, j, k) = phi(ilocal, jlocal);
        });
    }
}

void Interface::Reinit()
{
    BL_PROFILE("Interface::Reinit()");

    //interface_ = IF;
    if (use_FMM)
    {
        AlgoimReinit();
        return;
    }

    const amrex::Real *dx = mesh_.geometry().CellSize();
    amrex::Real Tol, tau, Residue;

    tau = 0.25 * Minimum2(dx[0], dx[1]);
    Tol = tau * tau;
    Tol = 1.0E-6;
    LS_Iter = 0;

    LSGammaIdentification();
    LSSIdentification();

    for (amrex::MFIter mfi(LSReinit_Source); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        auto source = LSReinit_Source.array(mfi);
        //auto psi = Psi().const_array(mfi);
        amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);

        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            source(i, j, k, 1) = f(i, j, k);
            source(i, j, k, 2) = f(i, j, k) / std::sqrt(f(i, j, k) * f(i, j, k) + dx[0] * dx[0]);
        });
    }

    do
    {
        for (amrex::MFIter mfi(LSReinit_Source); mfi.isValid(); ++mfi)
        {
            auto source = LSReinit_Source.array(mfi);
            //auto psi = Psi().const_array(mfi);
            amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);
            auto index = Index()[mfi];

            for (auto &&cell : index)
            {
                source(cell, 3) = f(cell);
            }
        }

        TubeComputeLevelSetRHS();

        for (amrex::MFIter mfi(LSReinit_Source); mfi.isValid(); ++mfi)
        {
            auto source = LSReinit_Source.array(mfi);
            //auto psi = Psi().array(mfi);
            amrex::Array4<amrex::Real> const &f = psi.array(mfi);
            auto index = Index()[mfi];

            for (auto &&cell : index)
            {
                f(cell) = source(cell, 3) + tau * source(cell, 0);
            }
        }

        PhaseFieldBC();
        TubeComputeLevelSetRHS();

        for (amrex::MFIter mfi(LSReinit_Source); mfi.isValid(); ++mfi)
        {
            auto source = LSReinit_Source.array(mfi);
            //auto psi = Psi().array(mfi);
            amrex::Array4<amrex::Real> const &f = psi.array(mfi);
            auto index = Index()[mfi];

            for (auto &&cell : index)
            {
                f(cell) = 0.75 * source(cell, 3) + 0.25 * f(cell) + 0.25 * tau * source(cell, 0);
            }
        }

        PhaseFieldBC();
        TubeComputeLevelSetRHS();

        for (amrex::MFIter mfi(LSReinit_Source); mfi.isValid(); ++mfi)
        {
            auto source = LSReinit_Source.array(mfi);
            //auto psi = Psi().array(mfi);
            amrex::Array4<amrex::Real> const &f = psi.array(mfi);
            auto index = Index()[mfi];

            for (auto &&cell : index)
            {
                f(cell) = (1.0/3.0) * source(cell, 3) + (2.0/3.0) * f(cell) + (2.0/3.0) * tau * source(cell, 0);
            }
        }

        PhaseFieldBC();
        Residue = 0.0;
        for (amrex::MFIter mfi(LSReinit_Source); mfi.isValid(); ++mfi)
        {
            auto source = LSReinit_Source.array(mfi);
            //auto psi = Psi().array(mfi);
            amrex::Array4<amrex::Real> const &f = psi.array(mfi);
            auto index = Index()[mfi];

            for (auto &&cell : index)
            {
                if (Residue < fabs(source(cell, 3) - f(cell)))
                    Residue = fabs(source(cell, 3) - f(cell));
            }
        }

        amrex::ParallelDescriptor::ReduceRealMax(Residue);
        LS_REINIT_TOL = Residue;
        LS_Iter++;

    } while ((Residue > Tol) && (LS_Iter < LS_MAXITER));
}

void Interface::LSGammaIdentification()
{
    LS_Gamma.setVal(-1);
    amrex::Real Xp, Xm, Yp, Ym;
    for (amrex::MFIter mfi(LS_Gamma); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        //auto psi = Psi().const_array(mfi);
        amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);
        auto ls_gamma = LS_Gamma.array(mfi);

        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Xp = f(i, j, k) * f(i + 1, j, k);
            Xm = f(i, j, k) * f(i - 1, j, k);
            Yp = f(i, j, k) * f(i, j + 1, k);
            Ym = f(i, j, k) * f(i, j - 1, k);

            if (fabs(Xp) < 1.0E-12 || fabs(Xm) < 1.0E-12 || fabs(Yp) < 1.0E-12 || fabs(Ym) < 1.0E-12)
            {
                ls_gamma(i, j, k) = 1;
            }
            if (Xp < 0.0 || Xm < 0.0 || Yp < 0.0 || Ym < 0.0)
            {
                ls_gamma(i, j, k) = 1;
            }
        });
    }
}

void Interface::LSSIdentification()
{
    amrex::Real Xp, Xm, Yp, Ym, sum;
    LS_S.setVal(-1);
    LS_R.setVal(0.0);
    for (amrex::MFIter mfi(LS_S); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        auto ls_s = LS_S.array(mfi);
        auto ls_r = LS_R.array(mfi);
        auto ls_gamma = LS_Gamma.const_array(mfi);
        //auto psi = Psi().const_array(mfi);
        amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);

        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (ls_gamma(i,j,k) > 0)
            {
                Xp = f(i, j, k) * f(i + 1, j, k);
                Xm = f(i, j, k) * f(i - 1, j, k);
                Yp = f(i, j, k) * f(i, j + 1, k);
                Ym = f(i, j, k) * f(i, j - 1, k);

                sum = 0.0 ;
                if( Xp < 0.0 ) { ls_s(i, j, k, 0) = 1 ; sum += f(i+1,j,k) ; }
                if( Xm < 0.0 ) { ls_s(i, j, k, 1) = 1 ; sum += f(i-1,j,k) ; }
                if( Yp < 0.0 ) { ls_s(i, j, k, 2) = 1 ; sum += f(i,j+1,k) ; }
                if( Ym < 0.0 ) { ls_s(i, j, k, 3) = 1 ; sum += f(i,j-1,k) ; }
                if
                (
                    ls_s(i, j, k, 0) == 1 || 
                    ls_s(i, j, k, 1) == 1 || 
                    ls_s(i, j, k, 2) == 1 || 
                    ls_s(i, j, k, 3) == 1
                )
                {
                    ls_r(i, j, k) = f(i, j, k) / sum;
                }
            }            
        });
    }
}

void Interface::TubeComputeLevelSetRHS()
{
    amrex::Real Psix_L, Psix_R, Psiy_L, Psiy_R, SignPsi;
    amrex::Real Psix_L_Plus, Psix_R_Plus, Psiy_L_Plus, Psiy_R_Plus;
    amrex::Real Psix_L_Minus, Psix_R_Minus, Psiy_L_Minus, Psiy_R_Minus;
    amrex::Real A, B;

    LSCIdentification();
    LSSum();

    const amrex::Real *dx = mesh_.geometry().CellSize();

    for (amrex::MFIter mfi(LSReinit_Source); mfi.isValid(); ++mfi)
    {
        auto source = LSReinit_Source.array(mfi);
        auto mask = Mask().const_array(mfi);
        //auto psi = Psi().const_array(mfi);
        amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);
        auto ls_r = LS_R.const_array(mfi);
        auto ls_sum = LS_Sum.const_array(mfi);
        auto ls_delta = LS_Dealta.const_array(mfi);

        auto index = Index()[mfi];
        for (auto &&cell : index)
        {
            int i = cell[0], j = cell[1], k = 0;
            if (mask(cell) != 3 && mask(cell) != 2 && mask(cell) != 1)
            {
                amrex::PrintToFile("log") << "Trouble in defining mask for tube T \n";
                amrex::PrintToFile("log") << mask(cell) << "\t" << cell << "\n";
                exit(1);
            }

            if (mask(cell) == 3)
            {
                WENO5_LS(Psix_L, Psix_R, Psiy_L, Psiy_R, i, j, dx, f);
            }
            else if (mask(cell) == 2 || mask(cell) == 1)
            {
                Psix_R = (f(i + 1, j, k) - f(i, j, k)) / dx[0];
                Psix_L = (f(i, j, k) - f(i - 1, j, k)) / dx[0];
                Psiy_R = (f(i, j + 1, k) - f(i, j, k)) / dx[1];
                Psiy_L = (f(i, j, k) - f(i, j - 1, k)) / dx[1];
            }
            SignPsi = source(i, j, k, 2);
            Psix_L_Plus = Maximum2(Psix_L,0.0) ; Psix_L_Minus = Minimum2(Psix_L,0.0) ;
            Psix_R_Plus = Maximum2(Psix_R,0.0) ; Psix_R_Minus = Minimum2(Psix_R,0.0) ;
            Psiy_L_Plus = Maximum2(Psiy_L,0.0) ; Psiy_L_Minus = Minimum2(Psiy_L,0.0) ;
            Psiy_R_Plus = Maximum2(Psiy_R,0.0) ; Psiy_R_Minus = Minimum2(Psiy_R,0.0) ;

            if( source(i,j,k,1) > 0.0 || fabs(source(i,j,k,1)) < 1.0E-12 ) 
            {
                A = Maximum2(Psix_L_Plus * Psix_L_Plus, Psix_R_Minus * Psix_R_Minus);
                B = Maximum2(Psiy_L_Plus * Psiy_L_Plus, Psiy_R_Minus * Psiy_R_Minus);
            }
            else
            {
                A = Maximum2(Psix_L_Minus * Psix_L_Minus, Psix_R_Plus * Psix_R_Plus);
                B = Maximum2(Psiy_L_Minus * Psiy_L_Minus, Psiy_R_Plus * Psiy_R_Plus);
            }

            source(i, j, k, 0) = SignPsi * (1.0 - sqrt(A + B));
            source(i, j, k, 0) += 0.5 * ls_delta(i, j, k) * (ls_r(i, j, k) * ls_sum(i, j, k) - f(i, j, k)) / dx[0];
        }
    }
}

//* NOTE : LS_C in rks 2d own-ls codes which is not used in that code
//*        So here it is skipped.. CONFIRM WITH RKS
void Interface::LSCIdentification()
{
    amrex::Real Xp, Xm, Yp, Ym;
    LS_Dealta.setVal(0.0);

    for (amrex::MFIter mfi(LS_Gamma); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();

        auto ls_delta = LS_Dealta.array(mfi);
        auto ls_gamma = LS_Gamma.const_array(mfi);
        auto ls_s = LS_S.const_array(mfi);
        //auto psi = Psi().const_array(mfi);
        amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);

        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if(ls_gamma(i,j,k) > 0)
            {
                Xp = f(i, j, k) * f(i + 1, j, k);
                Yp = f(i, j, k) * f(i, j + 1, k);
                Xm = f(i, j, k) * f(i - 1, j, k);
                Ym = f(i, j, k) * f(i, j - 1, k);

                if (Xp < 0.0 && ls_s(i, j, k, 0) > 0)
                    ls_delta(i, j, k) = 1.0;
                if (Xm < 0.0 && ls_s(i, j, k, 1) > 0)
                    ls_delta(i, j, k) = 1.0;
                if (Yp < 0.0 && ls_s(i, j, k, 2) > 0)
                    ls_delta(i, j, k) = 1.0;
                if (Ym < 0.0 && ls_s(i, j, k, 3) > 0)
                    ls_delta(i, j, k) = 1.0;
            }
        });
    }
}

void Interface::LSSum()
{
    amrex::Real sum;
    for (amrex::MFIter mfi(LS_Gamma); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        auto ls_sum = LS_Sum.array(mfi);
        auto ls_s = LS_S.const_array(mfi);
        auto ls_gamma = LS_Gamma.const_array(mfi);
        //auto psi = Psi().const_array(mfi);
        amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);

        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (ls_gamma(i,j,k) > 0)
            {
                sum = 0.0;
                if( ls_s(i, j, k, 0) == 1 ) { sum += f(i+1,j,k) ; }
                if( ls_s(i, j, k, 1) == 1 ) { sum += f(i-1,j,k) ; }
                if( ls_s(i, j, k, 2) == 1 ) { sum += f(i,j+1,k) ; }
                if( ls_s(i, j, k, 3) == 1 ) { sum += f(i,j-1,k) ; }
                
                ls_sum(i, j, k) = sum;
            }
        });
    }
}

} /*End namespace mycode */

