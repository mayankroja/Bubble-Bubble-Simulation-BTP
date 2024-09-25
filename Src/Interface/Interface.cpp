#include "Interface.H"
#include <AMReX_ParmParse.H>
#include <algoim_hocp.hpp>
#include <SolveCubicEqn.h>

namespace mycode
{

Interface::Interface(const amrexMesh &mesh, const std::string& name)
: 
mesh_(mesh),
Name(name)
{
    amrex::ParmParse pp;

    pp.get("Nghost", nghost);
    pp.query("DOF", DOF);
    pp.query("Axisymmetric", isAxisymmetric);
    pp.query("RKOrder", RKOrder);
    pp.get("LAYERS",LAYERS);
    
    psi.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
    psi.setVal(0.0);

    psi_old.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
    psi_old.setVal(0.0);

    psi_prev.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
    psi_prev.setVal(0.0);

    source_reinit.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
    source_reinit.setVal(0.0);

    kappa_.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
    kappa_.setVal(0.0);

    normal.define(mesh_.grid(), mesh_.dmap(), AMREX_SPACEDIM, nghost);
    normal.setVal(0.0);

    intercept_data.define(mesh_.grid(), mesh_.dmap());
    //For Reinit
    amrex::ParmParse pp1("LS_Reinit");
    LS_MAXITER = 20;
    use_FMM = false;
    LAYERS = 15;
    pp1.query("maxIter", LS_MAXITER);
    pp1.query("use_FMM", use_FMM);

    if (!use_FMM)
    {
        LS_Gamma.define(mesh_.grid(), mesh_.dmap(), 1, 0);
        LS_S.define(mesh_.grid(), mesh_.dmap(), 4, 0);
        // LS_C.define(mesh_.grid(), mesh_.dmap(), 4, 0);
        LS_R.define(mesh_.grid(), mesh_.dmap(), 1, 0);
        LS_Dealta.define(mesh_.grid(), mesh_.dmap(), 1, 0);
        LS_Sum.define(mesh_.grid(), mesh_.dmap(), 1, 0);
        LSReinit_Source.define(mesh_.grid(), mesh_.dmap(), 4, 0);
    }
    if(RKOrder >= 2)
    {
        FRK_psi.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
        FRK_psi.setVal(0.0);

        RK1_psi.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
        RK1_psi.setVal(0.0);

        RK2_psi.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
        RK2_psi.setVal(0.0);

        RK3_psi.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
        RK3_psi.setVal(0.0);
    }


}

Interface::~Interface() {}

void Interface::Reinit_algoim()
{
    //amrex::Print()<<"Algoim"<<'\n';
    PhaseFieldBC();
    const amrex::Real *dx = mesh_.geometry().CellSize();
    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.growntilebox();
        const amrex::Box &bx_valid = mfi.validbox();
        amrex::Array4<amrex::Real> const &f = psi.array(mfi);
        //amrex::Print(-1)<<"bx = "<<bx<<'\n';
        //amrex::Print(-1)<<"bx_valid = "<<bx_valid<<'\n';
        //amrex::Print(-1)<<"dx = "<<dx[0]<<'\n';

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
        Algoim::reinit<AMREX_SPACEDIM, 3>(phi, dx[0], 1.5 * LAYERS * dx[0]);

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


void Interface::Reinit2()
{
    int MAX_ITER = 50;
    amrex::Real C1, C2, C3;

    const amrex::Real *dx = mesh_.geometry().CellSize();

    int Iter = 0;

    amrex::Real reinit_tol = 0.0;
    amrex::MultiFab::Copy(psi_prev, psi, 0, 0, 1, nghost);
    do
    {
        //amrex::Print()<<"Reinit iter = "<<Iter<<"\n";

        const amrex::Real *prob_lo = mesh_.geometry().ProbLo();
        const amrex::Real *prob_hi = mesh_.geometry().ProbHi();

        amrex::MultiFab::Copy(psi_old, psi, 0, 0, 1, nghost);
        for(int istep = 0; istep < 2; istep++)
        {
            //amrex::Print()<<"step = "<<istep + 1<<'\n';
            if(istep == 0)
            {
                C1 = 1.0;
                C2 = 0.0;
                C3 = 1.0;
            }
            else if(istep == 1)
            {
                C1 = 0.5;
                C2 = 0.5;
                C3 = 0.5;
            }
            PhaseFieldBC();

            compute_source_reinit();
            
            amrex::Real dtau = 0.4*dx[0]; 
            for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();
                amrex::Array4<amrex::Real> const &f = psi.array(mfi);
                amrex::Array4<amrex::Real> const &fold = psi_old.array(mfi);
                amrex::Array4<amrex::Real> const &fprev = psi_prev.array(mfi);
                amrex::Array4<amrex::Real> const &L = source_reinit.array(mfi);
                amrex::Array4<int> const &mask = Mask_.array(mfi);
                /// copy phi in psi
                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if( f(i,j,k)*f(i + 1,j,k) > 0 &&
                        f(i,j,k)*f(i - 1,j,k) > 0 &&
                        f(i,j,k)*f(i,j + 1,k) > 0 &&
                        f(i,j,k)*f(i,j - 1,k) > 0 )
                    {
                        f(i , j , k) = C1*fold(i, j, k) + 
                                       C2*f(i, j, k) -
                                       C3*dtau*L(i, j, k);
                        f(i, j, k) = std::copysign(f(i, j, k),fprev(i, j, k));
                        //if(mask(i ,j ,k) == 1)
                        {
                            //amrex::Print()<<"L(i, j, k) = "<<std::fabs(L(i, j, k))<<'\n';
                            //reinit_tol = std::max(reinit_tol, std::fabs(L(i ,j ,k)));
                        }
                        //if(i == 507 && j == 531)
                        //{
                        //    amrex::Print()<<" i = "<<i<<" , j = "<<j<<'\n';
                        //    amrex::Print()<<" L( i, j, k) = "<<L(i,j,k)<<'\n';      
                        //    amrex::Print()<<"f ="<<f(i,j,k)<<" , fprev = "<<fprev(i,j,k)<<'\n';
                        //}
                    }
                    else if(mask(i,j,k) < 3)
                    {
                        f(i , j , k) = C1*fold(i, j, k) +    
                                       C2*f(i, j, k) -
                                       C3*dtau*L(i, j, k);
                    }
                     
                            //if(i == 507 && j == 531)
                            //{
                            //    amrex::Print()<<" i = "<<i<<" , j = "<<j<<'\n';
                            //    amrex::Print()<<" L( i, j, k) = "<<L(i,j,k)<<'\n';       
                            //    amrex::Print()<<"f ="<<f(i,j,k)<<" , fprev = "<<fprev(i,j,k)<<'\n';
                            //}

                });
                //amrex::Print()<<"Max L = "<<std::max(std::fabs(L))<<"\n";
            }
            amrex::ParallelDescriptor::ReduceRealMax(reinit_tol);

        }
        Iter++;


        reinit_tol = 0.0;
        for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real> const &f = psi.array(mfi);
            amrex::Array4<amrex::Real> const &fold = psi_old.array(mfi);
            amrex::Array4<int> const &mask = Mask_.array(mfi);
            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if(mask(i ,j ,k) > 0)
                {   
                    reinit_tol = std::max(reinit_tol, fabs((f(i,j,k)-fold(i,j,k))/f(i,j,k)));
                }
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(reinit_tol);

        //amrex::Print()<<"reinit_tol = "<<reinit_tol<<'\n';
    }
    while(Iter <= MAX_ITER && reinit_tol > max_reinit_tol);

}

void Interface::compute_source_reinit()
{
    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        amrex::Array4<amrex::Real> const &f = psi.array(mfi);
        amrex::Array4<amrex::Real> const &fprev = psi_prev.array(mfi);
        amrex::Array4<amrex::Real> const &L = source_reinit.array(mfi);
        const amrex::Real *dx = mesh_.geometry().CellSize();
        /// copy phi in psi
        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            amrex::Real sgn_f = std::copysign(1.0,fprev(i,j,k));

            amrex::Real f_xx_p = (f(i + 1, j, k) - 2.0*f(i, j, k) + f(i - 1, j, k))/(dx[0]*dx[0]);
            amrex::Real f_xx_e = (f(i + 2, j, k) - 2.0*f(i + 1, j, k) + f(i, j, k))/(dx[0]*dx[0]);
            amrex::Real f_xx_w = (f(i, j, k) - 2.0*f(i - 1, j, k) + f(i - 2, j, k))/(dx[0]*dx[0]);
            amrex::Real f_x_p = (f(i + 1, j, k) - f(i, j, k))/dx[0] - 0.5*dx[0]*minmod(f_xx_p, f_xx_e);
            amrex::Real f_x_m = (f(i, j, k) - f(i - 1, j, k))/dx[0] + 0.5*dx[0]*minmod(f_xx_p, f_xx_w);

            amrex::Real f_yy_p = (f(i, j + 1, k) - 2.0*f(i, j, k) + f(i, j - 1, k))/(dx[1]*dx[1]);
            amrex::Real f_yy_n = (f(i, j + 2, k) - 2.0*f(i, j + 1, k) + f(i, j, k))/(dx[1]*dx[1]);
            amrex::Real f_yy_s = (f(i, j, k) - 2.0*f(i, j - 1, k) + f(i, j - 2, k))/(dx[1]*dx[1]);
            amrex::Real f_y_p = (f(i , j + 1, k) - f(i, j, k))/dx[1] - 0.5*dx[1]*minmod(f_yy_p, f_yy_n);
            amrex::Real f_y_m = (f(i , j, k) - f(i, j - 1, k))/dx[1] + 0.5*dx[1]*minmod(f_yy_p, f_yy_s);

            amrex::Real H_ijk = 0.0;
            //amrex::Print()<<"f("<<i<<" , "<<j<<" , "<<k << ") = "<<f(i, j, k)<<" , sign = "<<sgn_f<<'\n';
            amrex::Real a_p = std::max(f_x_p,0.0);
            amrex::Real a_m = std::min(f_x_p,0.0); 
            amrex::Real b_p = std::max(f_x_m,0.0);
            amrex::Real b_m = std::min(f_x_m,0.0);
            amrex::Real c_p = std::max(f_y_p,0.0);
            amrex::Real c_m = std::min(f_y_p,0.0);
            amrex::Real d_p = std::max(f_y_m,0.0);
            amrex::Real d_m = std::min(f_y_m,0.0);
            if(sgn_f >= 0 )
                H_ijk = std::sqrt(
                           std::max(a_m*a_m,b_p*b_p) +
                           std::max(c_m*c_m,d_p*d_p));
            else
                H_ijk = std::sqrt(
                           std::max(a_p*a_p,b_m*b_m) +
                           std::max(c_p*c_p,d_m*d_m));
            L(i,j,k) = sgn_f*(H_ijk - 1.0); 
            //if(i == 127 && j == 141)
            //{
            //    amrex::Print()<<" i = "<<i<<" , j = "<<j<<'\n';
            //    amrex::Print()<<" f = "<<f(i ,j ,k)<<" ,f_e = "<<f(i+1,j,k)<<" , f_w = "<<f(i-1,j,k)<<" , f_n = "<<f(i, j+1, k)<<" , f_s = "<<f(i , j-1,k)<<'\n';
            //    amrex::Print()<<" L( i, j, k) = "<<L(i,j,k)<<" , sgn_f = "<<sgn_f<<'\n'; 
            //}

        });
    }
}


void Interface::PhaseFieldBC()
{
    psi.FillBoundary();
    const amrex::Box &domain = mesh_.geometry().Domain();
    int ngrow = psi.nGrow();
    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        amrex::Array4<amrex::Real> const &f = psi.array(mfi);

        int k = 0; // as 2d

        amrex::Box gbx(bx);
        gbx.grow(ngrow);

        /// left
        if (bx.smallEnd(0) == domain.smallEnd(0))
        {
            //int ig = bx.smallEnd(0) - ngrow;
            //int ip = bx.smallEnd(0) + ngrow - 1;

            int ig = bx.smallEnd(0) - 1;
            for (int i = bx.smallEnd(0); i < bx.smallEnd(0) + ngrow; i++)
            {
                for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++)
                {
                    //f(i-ngrow, j, k) = f(2 * ngrow - 1 - i, j, k);
                    f(ig, j, k) = f(i, j, k);
                }
                ig--;
            }
        }
        /// bottom
        if (bx.smallEnd(1) == domain.smallEnd(1))
        {
            //int jm = bx.smallEnd(1) - ngrow;
            //int jp = bx.smallEnd(1) + ngrow - 1;

            int jg = bx.smallEnd(1) - 1;
            for (int j = bx.smallEnd(1); j < bx.smallEnd(1) + ngrow; j++)
            {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++)
                {
                    f(i, jg, k) = f(i, j, k);
                }
                jg--;
            }
        }
        /// right
        if (bx.bigEnd(0) == domain.bigEnd(0))
        {
            int ig = bx.bigEnd(0) + 1;
            for (int i = bx.bigEnd(0); i > bx.bigEnd(0) - ngrow; i--)
            {
                for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++)
                {
                    f(ig, j, k) = f(i, j, k);
                }
                ig++;
            }
        }
        /// top
        if (bx.bigEnd(1) == domain.bigEnd(1))
        {
            int jg = bx.bigEnd(1) + 1;
            for (int j = bx.bigEnd(1); j > bx.bigEnd(1) - ngrow; j--)
            {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++)
                {
                    f(i, jg, k) = f(i, j, k);
                }
                jg++;
            }
        }               
    }
}

amrex::Real Interface::ComputePhiFromPsi(const amrex::Real &Ps)
{
    amrex::Real H, Ph;
    const amrex::Real *dx = mesh_.geometry().CellSize();
    amrex::Real EPS_LS = 0.75 * std::min(dx[0], dx[1]);
    Ph = 0.5 + 0.5 * tanh(Ps / (2.0 * EPS_LS)); // this should read 0.5 + 0.5*tanh(Ps/(2.0*EPS_LS))
    return Ph;
}

amrex::Real Interface::Function2(amrex::Real Value)
{
    if(Value < 0.0) return 0.0;
    else return Value*Value;
}
//TODO ... for a case where there is fluid inside interface
//TODO ... this will give vol outside the interface
//* Note : for Volume_ computation move over the nodal box.. 
//*        the higher end of nodal box is not considered
//*        so owner mask is not required
//TODO ... for a case where there is fluid inside interface
//TODO ... this will give vol outside the interface
//* Note : for Volume_ computation move over the nodal box.. 
//*        the higher end of nodal box is not considered
//*        so owner mask is not required

//* NOTE : the iteration is over lower nodal boxes but the actual 
//*        computation is over higher nodal box
//*        as nbrs are ij, i+1 j, i+1 j+1, i j+1 
//*        This will also give correct result
void Interface::ComputeVolume(bool init)
{
    if (isAxisymmetric)
    {
        //amrex::Print()<<"Computing Axisymmetric Vol :"<<'\n';
        computeVolumeAxisymmetric(init);
    }
    else
    {
        amrex::Real dr_, xloc, yloc, lx_, ly_;
        amrex::Real Vol = 0.0;
        const amrex::Real *dx = mesh_.geometry().CellSize();
        // const amrex::Box domain = mesh_.geometry().Domain();
        // amrex::Box nddomain = amrex::convert(domain, amrex::IntVect::TheNodeVector());

        Vol_Diff_ = Volume_ = 0.0;

        for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                Vol_Diff_ += dx[0] * dx[1] * ComputePhiFromPsi(f(i, j, k));
            });
        }

        amrex::ParallelDescriptor::ReduceRealSum(Vol_Diff_);

        for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);

            //! Note : need to iterate on nodal box..
            //!        do not consider higher end nodes

            amrex::Box ndbx = amrex::convert(bx, amrex::IntVect::TheNodeVector());
            //! shrink the higher end by 1 cell in both x and y dir
            ndbx.growHi(0, -1);
            ndbx.growHi(1, -1);

            amrex::ParallelFor(ndbx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                {
                    // all gas, CASE 1
                    Vol += dx[0] * dx[1];
                }
                else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                {
                    // only (i,j) is in liquid, CASE 2
                    Vol += dx[0] * dx[1];
                    lx_ = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                    ly_ = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                    if (lx_ < 0.0 || ly_ < 0.0)
                    {
                        amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                        exit(1);
                    }
                    else
                        Vol -= 0.5 * lx_ * ly_;
                }
                else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                {
                    // only (i+1,j) is in liquid, CASE 3
                    Vol += dx[0] * dx[1];
                    lx_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                    ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                    if (lx_ < 0.0 || ly_ < 0.0)
                    {
                        amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                        exit(1);
                    }
                    else
                        Vol -= 0.5 * lx_ * ly_;
                }
                else if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                {
                    // only (i,j+1) is in liquida, CASE 4
                    Vol += dx[0] * dx[1];
                    lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                    ly_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                    if (lx_ < 0.0 || ly_ < 0.0)
                    {
                        amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                        exit(1);
                    }
                    else
                        Vol -= 0.5 * lx_ * ly_;
                }
                else if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                {
                    // only (i+1,j+1) is in liquid, CASE 5
                    Vol += dx[0] * dx[1];
                    lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                    ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                    if (lx_ < 0.0 || ly_ < 0.0)
                    {
                        amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                        exit(1);
                    }
                    else
                        Vol -= 0.5 * lx_ * ly_;
                }
                else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                {
                    // nodes (i,j) and (i+1,j) are in liquid, CASE 6
                    yloc = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                    ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                    if (ly_ < 0.0 || yloc < 0.0)
                    {
                        amrex::Print() << " ly_i: " << yloc << "  ly_{i+1}: " << ly_ << "\n";
                        exit(1);
                    }
                    else
                        Vol += dx[0] * (ly_ + yloc) / 2.0;
                }
                else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                {
                    // nodes (i,j) and (i,j+1) are in liquid, CASE 7
                    xloc = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                    lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                    if (lx_ < 0.0 || xloc < 0.0)
                    {
                        amrex::Print() << " lx_j: " << xloc << "  lx_{j+1}: " << lx_ << "\n";
                        exit(1);
                    }
                    else
                        Vol += dx[1] * (xloc + lx_) / 2.0;
                }
                else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                {
                    // nodes (i,j) and (i+1,j+1) are in liquid, CASE 8
                    Vol += dx[0] * dx[1];
                    lx_ = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                    ly_ = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                    if (lx_ < 0.0 || ly_ < 0.0)
                    {
                        amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                        exit(1);
                    }
                    else
                        Vol -= 0.5 * lx_ * ly_;

                    lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                    ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                    if (lx_ < 0.0 || ly_ < 0.0)
                    {
                        amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                        exit(1);
                    }
                    else
                        Vol -= 0.5 * lx_ * ly_;
                }
                else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                {
                    // nodes (i+1,j) and (i,j+1) are in liquid, CASE 9
                    Vol += dx[0] * dx[1];
                    lx_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                    ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                    if (lx_ < 0.0 || ly_ < 0.0)
                    {
                        amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                        exit(1);
                    }
                    else
                        Vol -= 0.5 * lx_ * ly_;

                    lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                    ly_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                    if (lx_ < 0.0 || ly_ < 0.0)
                    {
                        amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                        exit(1);
                    }
                    else
                        Vol -= 0.5 * lx_ * ly_;
                }
                else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                {
                    // nodes (i+1,j) and (i+1,j+1) are in liquid, CASE 10
                    xloc = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                    lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                    Vol += dx[1] * (xloc + lx_) / 2.0;
                }
                else if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                { // nodes (i,j+1) and (i+1,j+1) are in liquid, CASE 11
                    yloc = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                    ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                    Vol += dx[0] * (ly_ + yloc) / 2.0;
                }
                else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                {
                    // node (i,j) is in gas, CASE 12
                    lx_ = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                    ly_ = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                    if (lx_ < 0.0 || ly_ < 0.0)
                    {
                        amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                        exit(1);
                    }
                    else
                        Vol += 0.5 * lx_ * ly_;
                }
                else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                {
                    // node (i+1,j) is in gas, CASE 13
                    lx_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                    ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                    if (lx_ < 0.0 || ly_ < 0.0)
                    {
                        amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                        exit(1);
                    }
                    else
                        Vol += 0.5 * lx_ * ly_;
                }
                else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                {
                    // node (i,j+1) is in gas, CASE 14
                    lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                    ly_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                    if (lx_ < 0.0 || ly_ < 0.0)
                    {
                        amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                        exit(1);
                    }
                    else
                        Vol += 0.5 * lx_ * ly_;
                }
                else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                {
                    // node (i+1,j+1) is in gas, CASE 15
                    lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                    ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                    if (lx_ < 0.0 || ly_ < 0.0)
                    {
                        amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                        exit(1);
                    }
                    else
                        Vol += 0.5 * lx_ * ly_;
                }
                else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                {
                    // all liquid., CASE 16
                }
                else
                {
                    amrex::Print() << "Logical error in phase identification for volume calculation"
                                    << "\n";
                    exit(1);
                }
            });
        }

        Volume_ = Vol;
        amrex::ParallelDescriptor::ReduceRealSum(Volume_);
        if (init)
        {
            Volume0_ = Volume_;
            Vol_Diff0_ = Vol_Diff_;
        }
    }
}

void Interface::computeVolumeAxisymmetric(bool init)
{

    amrex::Real dr_, xloc, yloc, lx_, ly_;
    amrex::Real Vol = 0.0;
    //amrexMesh mesh__ = getMesh();
    const amrex::Real *dx = mesh_.geometry().CellSize();
    const amrex::Real *prob_lo = mesh_.geometry().ProbLo();

    Vol_Diff_ = Volume_ = 0.0;

    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);

        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
        {
            amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
            Vol_Diff_ += M_PI * dx[0] * dx[1] * fabs(y) * ComputePhiFromPsi(f(i, j, k));
        });
    }

    amrex::ParallelDescriptor::ReduceRealSum(Vol_Diff_);

    //! iterate over mid line
    const amrex::Box &domain = mesh_.geometry().Domain();
    amrex::Box upper_half_domain(domain);
    amrex::IntVect dlo = domain.smallEnd();
    dlo.shift(1, domain.length(1) / 2);
    upper_half_domain.setSmall(dlo);

    //! nodal upper half domain box
    upper_half_domain.convert(amrex::IntVect::TheNodeVector());

    //! create the centre(actually 1 below centre) line box
    amrex::IntVect midlo = upper_half_domain.smallEnd();
    amrex::IntVect midhi = upper_half_domain.bigEnd();
    midhi.setVal(1, midlo[1] - 1);
    midlo.setVal(1, midlo[1] - 1);
    amrex::Box mid(midlo, midhi);

    //! first loop over mid line
    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);
        amrex::Box ndbx = convert(bx, amrex::IntVect::TheNodeVector());
        ndbx.growHi(0, -1);
        ndbx.growHi(1, -1);
        amrex::Box midisect = ndbx & mid;
        //amrex::Print()<<"midisect = "<<midisect<<'\n';
        //amrex::Print()<<"ndbx.SmallEnd(0) = "<<ndbx.smallEnd(0)<<" , ndbx.SmallEnd(1) = "<<ndbx.smallEnd(1)<<'\n';
        //amrex::Print()<<"ndbx.BigEnd(0) = "<<ndbx.bigEnd(0)<<" , ndbx.bigEnd(1) = "<<ndbx.bigEnd(1)<<'\n';
        if (midisect.ok())
        {
            //! iterate only on mid line
            amrex::ParallelFor(midisect,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                {
                    //amrex::Print()<<"j = "<<j<<" , f(i , j , k ) = "<<f(i,j,k)<<'\n';
                    amrex::Real y = prob_lo[1] + dx[1] * (j + 1 + 0.5);
                    if (f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                    {
                        Vol += M_PI * dx[0] * dx[1] * y;
                    }
                    else if (f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                    {   // nodes (i,j) and (i,j+1) are in liquid, CASE 7
                        //	xloc = (Psi[i+1][j][p]/(Psi[i+1][j][p] - f(i, j, k)))*(x[i+1] - x[i]) ;
                        lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                        if (lx_ < 0.0)
                        {
                            amrex::Print() << "  lx_{j+1}: " << lx_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * dx[1] * (y * lx_); // assume lx_ = xloc
                    }
                    else if (f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    {
                        //	xloc = (f(i, j, k)/(f(i, j, k) - Psi[i+1][j][p]))*(x[i+1] - x[i]) ;
                        lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                        if (lx_ < 0.0)
                        {
                            amrex::Print() << "  lx_{j+1}: " << lx_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * dx[1] * y * lx_; // assume lx_ = xloc
                    }
                    else if (f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    {
                    }
                    else
                    {
                        amrex::Print() << "Logical error in phase identification for volume calculation" << "\n";
                        amrex::Print() << f(i, j, k) << "\t" << f(i + 1, j, k) << "\t" << f(i, j + 1, k) << "\t" << f(i + 1, j + 1, k) << "\n";
                        // amrex::Print() << "Location:" << x[i] << "\t" << y1 << "\n";
                        // WriteFile();
                        // Write_Interface();
                        exit(1);
                    }
                }
            });
        }
    }

    //! loop over upper half nodal box
    //TODO : if shrink higher end by 1, extra volume on outer region will be counted
    //TODO : correct way would be to add 1/2 the contro at nodes near higher domain bndry
    //TODO : for now outer volm is not used so correct later
    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);

        //! Note : need to iterate on nodal box..
        //!        do not consider higher end nodes

        amrex::Box ndbx = amrex::convert(bx, amrex::IntVect::TheNodeVector());
        //! shrink the higher end by 1 cell in both x and y dir
        ndbx.growHi(0, -1);
        ndbx.growHi(1, -1);
        amrex::Box isect = ndbx & upper_half_domain;
        if (isect.ok())
        {
            amrex::ParallelFor(isect,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                //amrex::Print()<<"i = "<<i<<" , j ="<<j<<'\n';
                //if(j > (ndbx.bigEnd(1) - ndbx.smallEnd(1) - 1) / 2)
                {
                    //amrex::Print()<<"i = "<<i<<" , j ="<<j<<" , f(i , j ,k) = "<<f(i,j,k)<<'\n';
                    amrex::Real y1 = prob_lo[1] + dx[1] * (j + 0.5);
                    amrex::Real y2 = prob_lo[1] + dx[1] * (j + 1 + 0.5);
                    if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // all gas, CASE 1
                        Vol += M_PI * dx[0] * dx[1] * (y1 + y2);
                        //amrex::Print()<<"Vol = "<<Vol<<'\n';
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // only (i,j) is in liquid, CASE 2
                        Vol += M_PI * dx[0] * dx[1] * (y1 + y2);
                        lx_ = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                        ly_ = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y1 + ly_ / 3.0);
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // only (i+1,j) is in liquid, CASE 3
                        Vol += M_PI * dx[0] * dx[1] * (y1 + y2);
                        lx_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                        ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y1 + ly_ / 3.0);
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // only (i,j+1) is in liquida, CASE 4
                        Vol += M_PI * dx[0] * dx[1] * (y1 + y2);
                        lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                        ly_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y2 - ly_ / 3.0);
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // only (i+1,j+1) is in liquid, CASE 5
                        Vol += M_PI * dx[0] * dx[1] * (y1 + y2);
                        lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                        ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y2 - ly_ / 3.0);
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // nodes (i,j) and (i+1,j) are in liquid, CASE 6
                        yloc = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                        ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                        if (ly_ < 0.0 || yloc < 0.0)
                        {
                            amrex::Print() << " ly_i: " << yloc << "  ly_{i+1}: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * dx[0] * (y2 * (ly_ + yloc) - (ly_ * ly_ + ly_ * yloc + yloc * yloc) / 3.0);
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // nodes (i,j) and (i,j+1) are in liquid, CASE 7
                        xloc = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                        lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                        if (lx_ < 0.0 || xloc < 0.0)
                        {
                            amrex::Print() << " lx_j: " << xloc << "  lx_{j+1}: " << lx_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * dx[1] * (y1 * xloc + y2 * lx_ + dx[1] * (xloc - lx_) / 3.0);
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // nodes (i,j) and (i+1,j+1) are in liquid, CASE 8
                        Vol += M_PI * dx[0] * dx[1] * (y1 + y2);
                        lx_ = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                        ly_ = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y1 + ly_ / 3.0);
                    
                        lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                        ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y2 - ly_ / 3.0);
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // nodes (i+1,j) and (i,j+1) are in liquid, CASE 9
                        Vol += M_PI * dx[0] * dx[1] * (y1 + y2);
                        lx_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                        ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y1 + ly_ / 3.0);
                    
                        lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                        ly_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y2 - ly_ / 3.0);
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // nodes (i+1,j) and (i+1,j+1) are in liquid, CASE 10
                        xloc = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                        lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                        Vol += M_PI * dx[1] * (y1 * xloc + y2 * lx_ + dx[1] * (xloc - lx_) / 3.0);
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // nodes (i,j+1) and (i+1,j+1) are in liquid, CASE 11
                        yloc = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                        ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                        Vol += M_PI * dx[0] * (y1 * (ly_ + yloc) + (ly_ * ly_ + ly_ * yloc + yloc * yloc) / 3.0);
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // node (i,j) is in gas, CASE 12
                        lx_ = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                        ly_ = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * lx_ * ly_ * (y1 + ly_ / 3.0);
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // node (i+1,j) is in gas, CASE 13
                        lx_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                        ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * lx_ * ly_ * (y1 + ly_ / 3.0);
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // node (i,j+1) is in gas, CASE 14
                        lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                        ly_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * lx_ * ly_ * (y2 - ly_ / 3.0);
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // node (i+1,j+1) is in gas, CASE 15
                        lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                        ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * lx_ * ly_ * (y2 - ly_ / 3.0);
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // all liquid., CASE 16
                    }
                    else
                    {
                        amrex::Real x1 = prob_lo[0] + dx[0] * (i + 0.5);
                        amrex::Print() << "Logical error in phase identification for volume calculation"
                                       << "\n";
                        amrex::Print() << f(i, j, k) << "\t" << f(i + 1, j, k) << "\t" << f(i, j + 1, k) << "\t" << f(i + 1, j + 1, k) << "\n";
                        amrex::Print() << "Location:" << x1 << "\t" << y1 << "\n";
                        // WriteFile();
                        // Write_Interface();
                        exit(1);
                    }
                }
            });
        }
    }

    Volume_ = Vol;
    amrex::ParallelDescriptor::ReduceRealSum(Volume_);
    //amrex::Print()<<"Vol ="<<Vol <<'\n';
    if (init)
    {
        Volume0_ = Volume_;
        Vol_Diff0_ = Vol_Diff_;
    }

}

void Interface::ComputeVolume(bool init, const amrex::iMultiFab &cfmask_)
{
    if (isAxisymmetric)
    {
        //amrex::Print()<<"Computing Axisymmetric Vol :"<<'\n';
        computeVolumeAxisymmetric(init,cfmask_);
    }
    else
    {
        amrex::Real dr_, xloc, yloc, lx_, ly_;
        amrex::Real Vol = 0.0;
        const amrex::Real *dx = mesh_.geometry().CellSize();
        const amrex::Real *prob_lo = mesh_.geometry().ProbLo();
        amrex::Real vol_moment_diff = 0;
        // const amrex::Box domain = mesh_.geometry().Domain();
        // amrex::Box nddomain = amrex::convert(domain, amrex::IntVect::TheNodeVector());
        //amrex::Print()<<"M_PI = "<<M_PI<<'\n';

        Part_Vol_Diff_ = Part_Volume_ = 0.0;

        for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);
            amrex::Array4<int const> const &cfmask = cfmask_.const_array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if(!cfmask(i,j,k) == CFMask::covered)
		        {
		            amrex::Real xc = prob_lo[0] + dx[0] * (i + 0.5);
                    Part_Vol_Diff_ += dx[0] * dx[1] * ComputePhiFromPsi(f(i, j, k));
		            vol_moment_diff += xc * dx[0] * dx[1] * ComputePhiFromPsi(f(i, j, k));
		        }
            });
        }

        amrex::ParallelDescriptor::ReduceRealSum(Part_Vol_Diff_);
	amrex::ParallelDescriptor::ReduceRealSum(vol_moment_diff);

	xcp = vol_moment_diff/Part_Vol_Diff_;

        for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);
            amrex::Array4<int const> const &cfmask = cfmask_.const_array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if(!cfmask(i,j,k) == CFMask::covered)
                {
                    //Part_Vol_Diff_ += dx[0] * dx[1] * ComputePhiFromPsi(f(i, j, k));
                    amrex::Real xc = prob_lo[0] + dx[0] * (i + 0.5);
                    amrex::Real yc = prob_lo[1] + dx[1] * (j + 0.5);
                    amrex::Real zc = 0.0;

                    amrex::Real fp = f(i ,j, k);
                    amrex::Real fe = f(i + 1, j, k);
                    amrex::Real fw = f(i - 1, j, k);
                    amrex::Real fn = f(i, j + 1, k);
                    amrex::Real fs = f(i, j - 1, k);
                    amrex::Real ft = 0.0;//For 3D
                    amrex::Real fb = 0.0;//For 3D
             
                    amrex::Real nx = (fe - fw)/(2.0*dx[0]);
                    amrex::Real ny = (fn - fs)/(2.0*dx[1]); 
                    amrex::Real nz = 0.0;//For 3d, NEEDS TO BE FIXED LATER
                    amrex::Real norm = std::sqrt(nx*nx + ny*ny + nz*nz) + 1e-20;
                    amrex::Real xfoot = xc - fp*nx;
                    amrex::Real yfoot = yc - fp*ny;
                    amrex::Real zfoot = zc - fp*nz; 
 
                    nx = nx/norm + 1e-20;
                    ny = ny/norm + 1e-20;
                    nz = nz/norm + 1e-20;

                    // Equation of line tangent to the 0 ls contout at the foot of normal is
                    // ax+by=c
                    amrex::Real a = nx;
                    amrex::Real b = ny;
                    amrex::Real c = nx*xfoot + ny*yfoot;

                    //Coordiantes of the vertices of the rectangular cell 
                    amrex::Real x1 = xc - 0.5*dx[0];
                    amrex::Real y1 = yc - 0.5*dx[1];
                    amrex::Real sgn1 = copysign(1.0, a*x1 + b*y1 - c);
                 
                    amrex::Real x2 = xc + 0.5*dx[0];
                    amrex::Real y2 = yc - 0.5*dx[1];
                    amrex::Real sgn2 = copysign(1.0, a*x2 + b*y2 - c);

                    amrex::Real x3 = xc + 0.5*dx[0];
                    amrex::Real y3 = yc + 0.5*dx[1];
                    amrex::Real sgn3 = copysign(1.0, a*x3 + b*y3 - c);

                    amrex::Real x4 = xc - 0.5*dx[0];
                    amrex::Real y4 = yc + 0.5*dx[1];
                    amrex::Real sgn4 = copysign(1.0, a*x4 + b*y4 - c);

                    amrex::Real sgn5 = copysign(1.0, a*xc + b*yc - c);

                    /*if((i == 20 && j == 5) || (i == 13 && j == 6))
                    {
                        amrex::Print()<<"i = "<<i<<", j = "<<j<<'\n';
                        amrex::Print()<<"xc = "<<xc<<" , yc = "<<yc<<'\n';
                        amrex::Print()<<"xfoot = "<<xfoot<<" , yfoot = "<<yfoot<<"\n";
                        amrex::Print()<<"a = "<<a<<" , b = "<<b<<", c = "<<c<<"\n";
                        amrex::Print()<<"sgn1 = "<<sgn1<<" , sgn2 = "<<sgn2<<" , sgn3 = "<<sgn3<<" , sgn4 = "<<sgn4<<" , sgn5 = "<<sgn5<<'\n';
                    }*/


                    amrex::Real vol_frac = 0.0;

                    bool interfacial_point = false;

                    if((fp*fe < 0) || (fp*fw < 0) || (fp*fn < 0) || (fp*fs < 0))
                        interfacial_point = true;

                    if(fp > 0 && !interfacial_point) vol_frac = 1.0;
                    if(interfacial_point)
                    {
                        if( fp > 0 && std::abs(sgn1 + sgn2 + sgn3 + sgn4 + sgn5) == 5 )
                        {
                            vol_frac = 1.0;
                            //amrex::Print()<<"gas points at the interface: "<<i<<" , "<<j<<"\n";
                        }
                        else if(std::abs(sgn1 + sgn2 + sgn3 + sgn4 + sgn5) < 5)
                        {
                            //amrex::Print()<<"Fraction points at the interface: "<<i<<" , "<<j<<'\n';
                            if(std::abs(a/b) < 0.0001)
                            {
                                if(fp >= 0)
                                    vol_frac = 0.5 + std::abs(fp)/dx[1];
                                else if(fp < 0)
                                    vol_frac = 0.5 - std::abs(fp)/dx[1];
                            }
                            else if(std::abs(a/b) > 1000.0)
                            {
                                if(fp >= 0)
                                    vol_frac = 0.5 + std::abs(fp)/dx[0];
                                else if(fp < 0)
                                    vol_frac = 0.5 - std::abs(fp)/dx[0];

                            }
                            else if(a/b >= 0.0001)
                            {
                                amrex::Real a_prime = a;
                                amrex::Real b_prime = b;
                                amrex::Real c_prime = c - a*x1 -b*y1;
                                
                                vol_frac = ((c_prime*c_prime)/(2.0*a_prime*b_prime)) - 
                                           0.5 * Function2(c_prime/a_prime - dx[0])*(a_prime/b_prime) - 
                                           0.5 * Function2(c_prime/b_prime - dx[1])*(b_prime/a_prime); 
                                
                                vol_frac = vol_frac/(dx[0]*dx[1]);
                                if(fp <= 0.0 && (sgn1 == sgn5))
                                    vol_frac = 1 - vol_frac;
                                else if(fp >= 0.0 && (sgn1 != sgn5))
                                    vol_frac = 1 - vol_frac;  
                                
                            }
                            else if(a/b <= -0.0001)
                            {
                                amrex::Real a_prime = -1.0*a;
                                amrex::Real b_prime = b;
                                amrex::Real c_prime = c - a*x2 -b*y2;

                                vol_frac = ((c_prime*c_prime)/(2.0*a_prime*b_prime)) -  
                                           0.5 * Function2(c_prime/a_prime - dx[0])*(a_prime/b_prime) -  
                                           0.5 * Function2(c_prime/b_prime - dx[1])*(b_prime/a_prime);  

                                vol_frac = vol_frac/(dx[0]*dx[1]);
                                if(fp <= 0.0 && (sgn2 == sgn5))
                                    vol_frac = 1 - vol_frac;
                                else if(fp >= 0.0 && (sgn2 != sgn5))
                                    vol_frac = 1 - vol_frac;
                                 
                            }
                        }
                        //amrex::Print()<<"vol_frac = "<<vol_frac<<'\n';
                    }
                    Vol += dx[0] * dx[1] * vol_frac; 
                }
            });
        }

        /*
        for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);
            amrex::Array4<int const> const &cfmask = cfmask_.const_array(mfi);

            //! Note : need to iterate on nodal box..
            //!        do not consider higher end nodes

            amrex::Box ndbx = amrex::convert(bx, amrex::IntVect::TheNodeVector());
            //! shrink the higher end by 1 cell in both x and y dir
            ndbx.growHi(0, -1);
            ndbx.growHi(1, -1);
  
            amrex::ParallelFor(ndbx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                if(!cfmask(i,j,k) == CFMask::covered)
                {
                    if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                    {
                        // all gas, CASE 1
                        Vol += dx[0] * dx[1];
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                    {
                        // only (i,j) is in liquid, CASE 2
                        Vol += dx[0] * dx[1];
                        lx_ = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                        ly_ = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= 0.5 * lx_ * ly_;
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                    {
                        // only (i+1,j) is in liquid, CASE 3
                        Vol += dx[0] * dx[1];
                        lx_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                        ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= 0.5 * lx_ * ly_;
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                    {
                        // only (i,j+1) is in liquida, CASE 4
                        Vol += dx[0] * dx[1];
                        lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                        ly_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= 0.5 * lx_ * ly_;
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    {
                        // only (i+1,j+1) is in liquid, CASE 5
                        Vol += dx[0] * dx[1];
                        lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                        ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= 0.5 * lx_ * ly_;
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                    {
                        // nodes (i,j) and (i+1,j) are in liquid, CASE 6
                        yloc = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                        ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                        if (ly_ < 0.0 || yloc < 0.0)
                        {
                            amrex::Print() << " ly_i: " << yloc << "  ly_{i+1}: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += dx[0] * (ly_ + yloc) / 2.0;
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                    {
                        // nodes (i,j) and (i,j+1) are in liquid, CASE 7
                        xloc = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                        lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                        if (lx_ < 0.0 || xloc < 0.0)
                        {
                            amrex::Print() << " lx_j: " << xloc << "  lx_{j+1}: " << lx_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += dx[1] * (xloc + lx_) / 2.0;
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    {
                        // nodes (i,j) and (i+1,j+1) are in liquid, CASE 8
                        Vol += dx[0] * dx[1];
                        lx_ = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                        ly_ = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= 0.5 * lx_ * ly_;
                    
                        lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                        ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= 0.5 * lx_ * ly_;
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                    {
                        // nodes (i+1,j) and (i,j+1) are in liquid, CASE 9
                        Vol += dx[0] * dx[1];
                        lx_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                        ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= 0.5 * lx_ * ly_;
                    
                        lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                        ly_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= 0.5 * lx_ * ly_;
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    {
                        // nodes (i+1,j) and (i+1,j+1) are in liquid, CASE 10
                        xloc = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                        lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                        Vol += dx[1] * (xloc + lx_) / 2.0;
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // nodes (i,j+1) and (i+1,j+1) are in liquid, CASE 11
                        yloc = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                        ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                        Vol += dx[0] * (ly_ + yloc) / 2.0;
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    {
                        // node (i,j) is in gas, CASE 12
                        lx_ = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                        ly_ = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += 0.5 * lx_ * ly_;
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    {
                        // node (i+1,j) is in gas, CASE 13
                        lx_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                        ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += 0.5 * lx_ * ly_;
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    {
                        // node (i,j+1) is in gas, CASE 14
                        lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                        ly_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += 0.5 * lx_ * ly_;
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                    {
                        // node (i+1,j+1) is in gas, CASE 15
                        lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                        ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += 0.5 * lx_ * ly_;
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    {
                        // all liquid., CASE 16
                    }
                    else
                    {
                        amrex::Print() << "Logical error in phase identification for volume calculation"
                                        << "\n";
                        exit(1);
                    }
                }
            });
        }
        */

        Part_Volume_ = Vol;
        amrex::ParallelDescriptor::ReduceRealSum(Part_Volume_);
        if (init)
        {
            Part_Volume0_ = Part_Volume_;
            Part_Vol_Diff0_ = Part_Vol_Diff_;
        }
    }
}

void Interface::computeVolumeAxisymmetric(bool init, const amrex::iMultiFab &cfmask_)
{

    amrex::Real dr_, xloc, yloc, lx_, ly_;
    amrex::Real Vol = 0.0;
    amrex::Real cut_surface_area = 0.0;
    //amrexMesh mesh__ = getMesh();
    const amrex::Real *dx = mesh_.geometry().CellSize();
    const amrex::Real *prob_lo = mesh_.geometry().ProbLo();
    amrex::Real vol_moment_diff = 0;

    Part_Vol_Diff_ = Part_Volume_ = 0.0;
    part_surf_area = 0.0;

    Part_resolution = 0;

    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);
        amrex::Array4<int const> const &cfmask = cfmask_.const_array(mfi);

        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
        {
	    amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
            amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
            if(!cfmask(i,j,k) == CFMask::covered)
	    {
                Part_Vol_Diff_ += 2.0 * M_PI * dx[0] * dx[1] * fabs(y) * ComputePhiFromPsi(f(i, j, k));
                vol_moment_diff += 2.0 * M_PI * x * dx[0] * dx[1] * fabs(y) * ComputePhiFromPsi(f(i, j, k));
                if(f(i,j,k) >= 0.0)
                    Part_resolution++;
            }	
        });
    }
    amrex::ParallelDescriptor::ReduceRealSum(Part_Vol_Diff_);
    amrex::ParallelDescriptor::ReduceRealSum(vol_moment_diff);
    amrex::ParallelDescriptor::ReduceIntSum(Part_resolution);
    xcp = vol_moment_diff/Part_Vol_Diff_ + 1e-20;
    //amrex::Print()<<"xcp = "<<xcp<<'\n';

    
    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);
        amrex::Array4<int const> const &cfmask = cfmask_.const_array(mfi);
    
        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if(!cfmask(i,j,k) == CFMask::covered)
            {
                //Part_Vol_Diff_ += dx[0] * dx[1] * ComputePhiFromPsi(f(i, j, k));
                amrex::Real xc = prob_lo[0] + dx[0] * (i + 0.5);
                amrex::Real yc = prob_lo[1] + dx[1] * (j + 0.5);
                amrex::Real zc = 0.0;
    
                amrex::Real fp = f(i ,j, k);
                amrex::Real fe = f(i + 1, j, k);
                amrex::Real fw = f(i - 1, j, k);
                amrex::Real fn = f(i, j + 1, k);
                amrex::Real fs = f(i, j - 1, k);
                amrex::Real ft = 0.0;//For 3D
                amrex::Real fb = 0.0;//For 3D
            
                amrex::Real nx = (fe - fw)/(2.0*dx[0]);
                amrex::Real ny = (fn - fs)/(2.0*dx[1]); 
                amrex::Real nz = 0.0;//For 3d, NEEDS TO BE FIXED LATER
                amrex::Real norm = std::sqrt(nx*nx + ny*ny + nz*nz) + 1e-20;
                amrex::Real xfoot = xc - fp*nx;
                amrex::Real yfoot = yc - fp*ny;
                amrex::Real zfoot = zc - fp*nz; 
     
                nx = nx/norm + 1e-20;
                ny = ny/norm + 1e-20;
                nz = nz/norm + 1e-20;
    
                // Equation of line tangent to the 0 ls contout at the foot of normal is
                // ax+by=c
                amrex::Real a = nx;
                amrex::Real b = ny;
                amrex::Real c = nx*xfoot + ny*yfoot;
    
                //Coordiantes of the vertices of the rectangular cell 
                amrex::Real x1 = xc - 0.5*dx[0];
                amrex::Real y1 = yc - 0.5*dx[1];
                amrex::Real sgn1 = copysign(1.0, a*x1 + b*y1 - c);
                
                amrex::Real x2 = xc + 0.5*dx[0];
                amrex::Real y2 = yc - 0.5*dx[1];
                amrex::Real sgn2 = copysign(1.0, a*x2 + b*y2 - c);
    
                amrex::Real x3 = xc + 0.5*dx[0];
                amrex::Real y3 = yc + 0.5*dx[1];
                amrex::Real sgn3 = copysign(1.0, a*x3 + b*y3 - c);
    
                amrex::Real x4 = xc - 0.5*dx[0];
                amrex::Real y4 = yc + 0.5*dx[1];
                amrex::Real sgn4 = copysign(1.0, a*x4 + b*y4 - c);
    
                amrex::Real sgn5 = copysign(1.0, a*xc + b*yc - c);
    
                amrex::Real vol_frac = 0.0;

                amrex::Real local_cut_surface_area = 0.0;
    
                bool interfacial_point = false;
    
                if((fp*fe < 0) || (fp*fw < 0) || (fp*fn < 0) || (fp*fs < 0))
                    interfacial_point = true;
    
                if(fp > 0 && !interfacial_point) vol_frac = 1.0;
                if(interfacial_point)
                {
                    if( fp > 0 && std::abs(sgn1 + sgn2 + sgn3 + sgn4 + sgn5) == 5 )
                    {
                        vol_frac = 1.0;
                        //amrex::Print()<<"gas points at the interface: "<<i<<" , "<<j<<"\n";
                    }
                    else if(std::abs(sgn1 + sgn2 + sgn3 + sgn4 + sgn5) < 5)
                    {
                        //amrex::Print()<<"Fraction points at the interface: "<<i<<" , "<<j<<'\n';
                        if(std::abs(a/b) < 0.0001)
                        {
                            if(fp >= 0)
                                vol_frac = 0.5 + std::abs(fp)/dx[1];
                            else if(fp < 0)
                                vol_frac = 0.5 - std::abs(fp)/dx[1];
                        }
                        else if(std::abs(a/b) > 1000.0)
                        {
                            if(fp >= 0)
                                vol_frac = 0.5 + std::abs(fp)/dx[0];
                            else if(fp < 0)
                                vol_frac = 0.5 - std::abs(fp)/dx[0];
    
                        }
                        else if(a/b >= 0.0001)
                        {
                            amrex::Real a_prime = a;
                            amrex::Real b_prime = b;
                            amrex::Real c_prime = c - a*x1 -b*y1;
                            
                            //vol_frac = (M_PI*(y1 + (c_prime/b_prime)/3.0) * (c_prime*c_prime)/(2.0*a_prime*b_prime)) - 
                            //            M_PI*(y1 + (c_prime/a_prime - dx[0])*(a_prime/b_prime)/3.0) * 0.5 * Function2(c_prime/a_prime - dx[0])*(a_prime/b_prime) - 
                            //            M_PI*(y1 + (c_prime/b_prime - dx[1])/3.0) * 0.5 * Function2(c_prime/b_prime - dx[1])*(b_prime/a_prime); 
                            
                            //vol_frac = vol_frac/(M_PI*fabs(yc)*dx[0]*dx[1]);
                            vol_frac = ( (c_prime*c_prime)/(2.0*a_prime*b_prime)) - 
                                         0.5 * Function2(c_prime/a_prime - dx[0])*(a_prime/b_prime) -
                                         0.5 * Function2(c_prime/b_prime - dx[1])*(b_prime/a_prime);

                            vol_frac = vol_frac/(dx[0]*dx[1]);

                            if(fp <= 0.0 && (sgn1 == sgn5))
                                vol_frac = 1 - vol_frac;
                            else if(fp >= 0.0 && (sgn1 != sgn5))
                                vol_frac = 1 - vol_frac;  
                            
                        }
                        else if(a/b <= -0.0001)
                        {
                            amrex::Real a_prime = -1.0*a;
                            amrex::Real b_prime = b;
                            amrex::Real c_prime = c - a*x2 -b*y2;
    
                            vol_frac = ((c_prime*c_prime)/(2.0*a_prime*b_prime)) -  
                                        0.5 * Function2(c_prime/a_prime - dx[0])*(a_prime/b_prime) -  
                                        0.5 * Function2(c_prime/b_prime - dx[1])*(b_prime/a_prime);  
    
                            vol_frac = vol_frac/(dx[0]*dx[1]);
                            if(fp <= 0.0 && (sgn2 == sgn5))
                                vol_frac = 1 - vol_frac;
                            else if(fp >= 0.0 && (sgn2 != sgn5))
                                vol_frac = 1 - vol_frac;
                                
                        }

                        //Compute the area of cut surface
                        double irdist = 0.5*std::sqrt(std::pow(dx[0],2) + std::pow(dx[1],2));
                        std::vector<double> ix;
                        std::vector<double> iy;

                        double ix1 = (c - b*y1)/a;
                        double iy1 = y1;
                        double ir1 = std::sqrt(std::pow(ix1 - xc,2) + std::pow(iy1 - yc,2));
                        if(ir1 <= irdist)
                        {
                            ix.push_back(ix1);
                            iy.push_back(iy1);
                        }

                        double ix2 = x2;
                        double iy2 = (c - a*x2)/b;
                        double ir2 = std::sqrt(std::pow(ix2 - xc,2) + std::pow(iy2 - yc,2));
                        if(ir2 <= irdist)
                        {
                            ix.push_back(ix2);
                            iy.push_back(iy2);
                        }

                        double ix3 = (c - b*y3)/a;
                        double iy3 = y3;
                        double ir3 = std::sqrt(std::pow(ix3 - xc,2) + std::pow(iy3 - yc,2));
                        if(ir3 <= irdist)
                        {
                            ix.push_back(ix3);
                            iy.push_back(iy3);
                        }

                        double ix4 = x4;
                        double iy4 = (c - a*x4)/b;
                        double ir4 = std::sqrt(std::pow(ix4 - xc,2) + std::pow(iy4 - yc,2));
                        if(ir4 <= irdist)
                        {
                            ix.push_back(ix4);
                            iy.push_back(iy4);
                        }

                        if(ix.size() > 2 || iy.size() > 2)
                        {
                            amrex::Print()<<"Error in surface area calculation: more than two cut points found"<<'\n';
                            std::exit(9);      
                        }
                        else if(ix.size() == 2 && iy.size() == 2)
                            local_cut_surface_area = std::sqrt(std::pow(ix[0] - ix[1],2) + std::pow(iy[0] - iy[1],2));

                    }
                    //amrex::Print()<<"vol_frac = "<<vol_frac<<'\n';
                }
                Vol += 2.0 * M_PI * dx[0] * dx[1] * fabs(yc) * vol_frac;
                cut_surface_area += 2.0 * M_PI * yfoot * local_cut_surface_area; 
            }
        });
    }
    
    /*
    //! iterate over mid line
    const amrex::Box &domain = mesh_.geometry().Domain();
    amrex::Box upper_half_domain(domain);
    amrex::IntVect dlo = domain.smallEnd();
    dlo.shift(1, domain.length(1) / 2);
    upper_half_domain.setSmall(dlo);
   
    //amrex::Print()<<"upper_half_domain = "<<upper_half_domain<<'\n';

    //! nodal upper half domain box
    upper_half_domain.convert(amrex::IntVect::TheNodeVector());

    //amrex::Print()<<"Nodal upper_half_domain = "<<upper_half_domain<<'\n'; 

    //! create the centre(actually 1 below centre) line box
    amrex::IntVect midlo = upper_half_domain.smallEnd();
    amrex::IntVect midhi = upper_half_domain.bigEnd();
    midhi.setVal(1, midlo[1] - 1);
    midlo.setVal(1, midlo[1] - 1);
    amrex::Box mid(midlo, midhi);

    //amrex::Print()<<"midlo = "<<midlo<<'\n';
    //amrex::Print()<<"midhi = "<<midhi<<'\n';
    //amrex::Print()<<"mid = "<<mid<<'\n';

    //! first loop over mid line
    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);
        amrex::Array4<int const> const &cfmask = cfmask_.const_array(mfi);
        amrex::Box ndbx = convert(bx, amrex::IntVect::TheNodeVector());
        ndbx.growHi(0, -1);
        ndbx.growHi(1, -1);
        amrex::Box midisect = ndbx & mid;
        //amrex::Print()<<"midisect = "<<midisect<<" , midisect.ok() = "<<midisect.ok()<<'\n';
        //amrex::Print()<<"ndbx.SmallEnd(0) = "<<ndbx.smallEnd(0)<<" , ndbx.SmallEnd(1) = "<<ndbx.smallEnd(1)<<'\n';
        //amrex::Print()<<"ndbx.BigEnd(0) = "<<ndbx.bigEnd(0)<<" , ndbx.bigEnd(1) = "<<ndbx.bigEnd(1)<<'\n';
        //amrex::Print()<<" dx = "<<dx[0]<<" , "<<dx[1]<<'\n';
        if (midisect.ok())
        {
            //! iterate only on mid line
            amrex::ParallelFor(midisect,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                if(!cfmask(i,j,k) == CFMask::covered)
                {
                    amrex::Real y = prob_lo[1] + dx[1] * (j + 1 + 0.5);
                    //amrex::Print()<<"j = "<<j<<" , f("<<i<<" , "<<j<<" , "<<k<<" ) = "<<f(i,j,k)<<" ,y = "<<y<<'\n';
                    if (f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                    {
                        Vol += M_PI * dx[0] * dx[1] * y;
                    }
                    else if (f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                    {   // nodes (i,j) and (i,j+1) are in liquid, CASE 7
                        //	xloc = (Psi[i+1][j][p]/(Psi[i+1][j][p] - f(i, j, k)))*(x[i+1] - x[i]) ;
                        lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                        if (lx_ < 0.0)
                        {
                            amrex::Print() << "  lx_{j+1}: " << lx_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * dx[1] * (y * lx_); // assume lx_ = xloc
                    }
                    else if (f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    {
                        //	xloc = (f(i, j, k)/(f(i, j, k) - Psi[i+1][j][p]))*(x[i+1] - x[i]) ;
                        lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                        if (lx_ < 0.0)
                        {
                            amrex::Print() << "  lx_{j+1}: " << lx_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * dx[1] * y * lx_; // assume lx_ = xloc
                    }
                    else if (f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    {
                    }
                    else
                    {
                        amrex::Print() << "Logical error in phase identification for volume calculation" << "\n";
                        amrex::Print() << f(i, j, k) << "\t" << f(i + 1, j, k) << "\t" << f(i, j + 1, k) << "\t" << f(i + 1, j + 1, k) << "\n";
                        // amrex::Print() << "Location:" << x[i] << "\t" << y1 << "\n";
                        // WriteFile();
                        // Write_Interface();
                        exit(1);
                    }
                }
            });
        }
    }

    //! loop over upper half nodal box
    //TODO : if shrink higher end by 1, extra volume on outer region will be counted
    //TODO : correct way would be to add 1/2 the contro at nodes near higher domain bndry
    //TODO : for now outer volm is not used so correct later
    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);
        amrex::Array4<int const> const &cfmask = cfmask_.const_array(mfi);

        //! Note : need to iterate on nodal box..
        //!        do not consider higher end nodes

        amrex::Box ndbx = amrex::convert(bx, amrex::IntVect::TheNodeVector());
        //amrex::Print()<<"bx = "<<bx<<'\n';
        //amrex::Print()<<"ndbx 1= "<<ndbx<<'\n';
        //! shrink the higher end by 1 cell in both x and y dir
        ndbx.growHi(0, -1);
        ndbx.growHi(1, -1);
        //amrex::Print()<<"ndbx 2= "<<ndbx<<'\n';
        amrex::Box isect = ndbx & upper_half_domain;
        //amrex::Print()<<"isect = "<<isect<<"isect.ok() = "<<isect.ok()<<'\n';
        if (isect.ok())
        {
            amrex::ParallelFor(isect,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                //amrex::Print()<<"i = "<<i<<" , j ="<<j<<'\n';
                //if(j > (ndbx.bigEnd(1) - ndbx.smallEnd(1) - 1) / 2)
                if(!cfmask(i,j,k) == CFMask::covered)
                {
                    //amrex::Print()<<"i = "<<i<<" , j ="<<j<<" , f(i , j ,k) = "<<f(i,j,k)<<'\n';
                    amrex::Real y1 = prob_lo[1] + dx[1] * (j + 0.5);
                    amrex::Real y2 = prob_lo[1] + dx[1] * (j + 1 + 0.5);
                    if(i == 24 && j == 43)
                    {
                        amrex::Print()<<"i = "<<i<<" , j = "<<j<<'\n';
                        amrex::Print()<<"y1 = "<<y1<<" , y2 = "<<y2<<'\n';
                    }
                    if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // all gas, CASE 1
                        Vol += M_PI * dx[0] * dx[1] * (y1 + y2);
                        //amrex::Print()<<"Vol = "<<Vol<<'\n';
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // only (i,j) is in liquid, CASE 2
                        Vol += M_PI * dx[0] * dx[1] * (y1 + y2);
                        lx_ = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                        ly_ = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y1 + ly_ / 3.0);
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // only (i+1,j) is in liquid, CASE 3
                        Vol += M_PI * dx[0] * dx[1] * (y1 + y2);
                        lx_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                        ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y1 + ly_ / 3.0);
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // only (i,j+1) is in liquida, CASE 4
                        Vol += M_PI * dx[0] * dx[1] * (y1 + y2);
                        lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                        ly_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y2 - ly_ / 3.0);
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // only (i+1,j+1) is in liquid, CASE 5
                        Vol += M_PI * dx[0] * dx[1] * (y1 + y2);
                        lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                        ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y2 - ly_ / 3.0);
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // nodes (i,j) and (i+1,j) are in liquid, CASE 6
                        yloc = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                        ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                        if (ly_ < 0.0 || yloc < 0.0)
                        {
                            amrex::Print() << " ly_i: " << yloc << "  ly_{i+1}: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * dx[0] * (y2 * (ly_ + yloc) - (ly_ * ly_ + ly_ * yloc + yloc * yloc) / 3.0);
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // nodes (i,j) and (i,j+1) are in liquid, CASE 7
                        xloc = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                        lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                        if (lx_ < 0.0 || xloc < 0.0)
                        {
                            amrex::Print() << " lx_j: " << xloc << "  lx_{j+1}: " << lx_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * dx[1] * (y1 * xloc + y2 * lx_ + dx[1] * (xloc - lx_) / 3.0);
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // nodes (i,j) and (i+1,j+1) are in liquid, CASE 8
                        Vol += M_PI * dx[0] * dx[1] * (y1 + y2);
                        lx_ = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                        ly_ = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y1 + ly_ / 3.0);
                    
                        lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                        ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y2 - ly_ / 3.0);
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // nodes (i+1,j) and (i,j+1) are in liquid, CASE 9
                        Vol += M_PI * dx[0] * dx[1] * (y1 + y2);
                        lx_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                        ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y1 + ly_ / 3.0);
                    
                        lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                        ly_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol -= M_PI * lx_ * ly_ * (y2 - ly_ / 3.0);
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // nodes (i+1,j) and (i+1,j+1) are in liquid, CASE 10
                        xloc = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                        lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                        Vol += M_PI * dx[1] * (y1 * xloc + y2 * lx_ + dx[1] * (xloc - lx_) / 3.0);
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // nodes (i,j+1) and (i+1,j+1) are in liquid, CASE 11
                        yloc = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                        ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                        Vol += M_PI * dx[0] * (y1 * (ly_ + yloc) + (ly_ * ly_ + ly_ * yloc + yloc * yloc) / 3.0);
                    }
                    else if (f(i, j, k) > 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // node (i,j) is in gas, CASE 12
                        lx_ = (f(i, j, k) / (f(i, j, k) - f(i + 1, j, k))) * dx[0];
                        ly_ = (f(i, j, k) / (f(i, j, k) - f(i, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * lx_ * ly_ * (y1 + ly_ / 3.0);
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) > 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // node (i+1,j) is in gas, CASE 13
                        lx_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i, j, k))) * dx[0];
                        ly_ = (f(i + 1, j, k) / (f(i + 1, j, k) - f(i + 1, j + 1, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * lx_ * ly_ * (y1 + ly_ / 3.0);
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) > 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // node (i,j+1) is in gas, CASE 14
                        lx_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i + 1, j + 1, k))) * dx[0];
                        ly_ = (f(i, j + 1, k) / (f(i, j + 1, k) - f(i, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * lx_ * ly_ * (y2 - ly_ / 3.0);
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) > 0.0)
                    { // node (i+1,j+1) is in gas, CASE 15
                        lx_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i, j + 1, k))) * dx[0];
                        ly_ = (f(i + 1, j + 1, k) / (f(i + 1, j + 1, k) - f(i + 1, j, k))) * dx[1];
                        if (lx_ < 0.0 || ly_ < 0.0)
                        {
                            amrex::Print() << " lx: " << lx_ << "  ly: " << ly_ << "\n";
                            exit(1);
                        }
                        else
                            Vol += M_PI * lx_ * ly_ * (y2 - ly_ / 3.0);
                    }
                    else if (f(i, j, k) <= 0.0 && f(i + 1, j, k) <= 0.0 && f(i, j + 1, k) <= 0.0 && f(i + 1, j + 1, k) <= 0.0)
                    { // all liquid., CASE 16
                    }
                    else
                    {
                        amrex::Real x1 = prob_lo[0] + dx[0] * (i + 0.5);
                        amrex::Print() << "Logical error in phase identification for volume calculation"
                                       << "\n";
                        amrex::Print() << f(i, j, k) << "\t" << f(i + 1, j, k) << "\t" << f(i, j + 1, k) << "\t" << f(i + 1, j + 1, k) << "\n";
                        amrex::Print() << "Location:" << x1 << "\t" << y1 << "\n";
                        // WriteFile();
                        // Write_Interface();
                        exit(1);
                    }
                }
            });
        }
    }
   */ 
    Part_Volume_ = Vol;
    amrex::ParallelDescriptor::ReduceRealSum(Part_Volume_);

    part_surf_area = cut_surface_area;
    amrex::ParallelDescriptor::ReduceRealSum(part_surf_area);
    //amrex::Print()<<"Vol ="<<Vol <<'\n';
    if (init)
    {
        Part_Volume0_ = Part_Volume_;
        Part_Vol_Diff0_ = Part_Vol_Diff_;
        Part_resolution0 = Part_resolution;
    }

}


void Interface::ComputeAvgIntVel(const amrex::MultiFab &u, const amrex::MultiFab &v, const amrex::iMultiFab &PMask)
{
    //amrex::Print()<<"ComputeAvgIntVel"<<'\n';
    int myproc = amrex::ParallelDescriptor::MyProc();
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    AvgRadius_ = 0.0;
    AvgIntVel_ = 0.0;
    AvgIntPres_ = 0.0;

    if (isAxisymmetric)
    {
    const amrex::Real *prob_lo = mesh_.geometry().ProbLo();
    const amrex::Real *dx = mesh_.geometry().CellSize();
    amrex::Real xloc, yloc,ravg = 0.0,r = 0.0 ,ur_avg= 0.0,p_I_avg = 0.0 ,icpt = 0.0, e_I_avg = 0.0, icpt_e = 0.0;
    
    //Assuming bubble remains spherical
        for (amrex::MFIter mfi(PMask); mfi.isValid(); ++mfi)
        {
            amrex::Array4<amrex::Real const> const &xvel = u.const_array(mfi);
            amrex::Array4<amrex::Real const> const &yvel = v.const_array(mfi);
                
            const amrex::Box &bx = mfi.validbox();
            auto &icpt_data = getInterceptData()[mfi];

            for (auto &&idt : icpt_data)
            {
                if (bx.contains(idt.cellid_))
                {
                    for (int i = 0; i < idt.n_intercepts; i++)
                    {
                        int i0 = idt.cellid_[0], j0 = idt.cellid_[1], k0_ = 0;
                        icpt++;
                        if (idt.type_[i] == 0)
                        {
                            xloc = prob_lo[0] + dx[0] * (i0 + 0.5 - idt.frac_[i]);
                            yloc = prob_lo[1] + dx[1] * (j0 + 0.5);
                        }
                        else if (idt.type_[i] == 1)
                        {
                            xloc = prob_lo[0] + dx[0] * (i0 + 0.5); 
                            yloc = prob_lo[1] + dx[1] * (j0 + 0.5 - idt.frac_[i]);
                        }

                       
                        r = sqrt((xloc - xcp)*(xloc - xcp) + (yloc-ycp)*(yloc - ycp)); 
                        ravg += r; 

                        ur_avg +=  (xvel(i0,j0,k0_)*(xloc - xcp)/r + yvel(i0,j0,k0_)*(yloc-ycp)/r); 
                        p_I_avg += idt.P[i];
                        amrex::Real ep_max;
                        {
                            amrex::Real C11 = idt.F11_Int[i]*idt.F11_Int[i] + idt.F21_Int[i]*idt.F21_Int[i];
                            amrex::Real C12 = idt.F11_Int[i]*idt.F12_Int[i] + idt.F21_Int[i]*idt.F22_Int[i];
                            amrex::Real C21 = idt.F12_Int[i]*idt.F11_Int[i] + idt.F22_Int[i]*idt.F21_Int[i];
                            amrex::Real C22 = idt.F12_Int[i]*idt.F12_Int[i] + idt.F22_Int[i]*idt.F22_Int[i];
                            amrex::Real C33 = idt.F33_Int[i]*idt.F33_Int[i];

                            amrex::Real first_inv = C11+C22+C33;
                            amrex::Real second_inv = C22*C33 + C11*C33 + C11*C22 - C12*C21;
                            amrex::Real third_inv = C33*(C11*C22 - C12*C21);

                            amrex::Real lambda_1;// = 0.5*(tr_C + std::sqrt(tr_C*tr_C - 4.0*det_C));
                            amrex::Real lambda_2;// = 0.5*(tr_C - std::sqrt(tr_C*tr_C - 4.0*det_C));
                            amrex::Real lambda_3;
                            SolveCubicEqn(lambda_1, lambda_2, lambda_3, first_inv, second_inv, third_inv);
                            //if(std::isnan(lambda_1)) lambda_1 = 1.0;
                            //if(std::isnan(lambda_2)) lambda_2 = 1.0;
                            //if(std::isnan(lambda_3)) lambda_3 = 1.0;
			    if(std::isnan(lambda_1)||std::isnan(lambda_1)||std::isnan(lambda_1))
		            {

			    }
			    else
		            {
				icpt_e++;
                                ep_max = std::max(std::abs(0.5*(lambda_1 - 1.0)),
                                               std::abs(0.5*(lambda_2 - 1.0)));
                                ep_max = std::max(ep_max, std::abs(0.5*(lambda_3 - 1.0)));

			        //amrex::Print()<<"ep_max = "<<ep_max<<", C11 = "<<C11<<" , C22 = "<<C22<<" , idt.P[i] = "<<idt.P[i]<<'\n';

                                e_I_avg += ep_max;
			    }
                        }
			//Compute force on interface
			F_p += idt.P[i]*idt.psix_[i]/std::hypot(idt.psix_[i],idt.psiy_[i]);
			F_visc_nn += idt.norm_shear_[i]*idt.psix_[i]/std::hypot(idt.psix_[i],idt.psiy_[i]);
			F_visc_nt += idt.tan_shear_[i]*idt.psix_[i]/std::hypot(idt.psix_[i],idt.psiy_[i]);
                    }
                }
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(ravg);
        amrex::ParallelDescriptor::ReduceRealSum(ur_avg);
        amrex::ParallelDescriptor::ReduceRealSum(p_I_avg); 
	amrex::ParallelDescriptor::ReduceRealSum(e_I_avg);
        amrex::ParallelDescriptor::ReduceRealSum(icpt);
	amrex::ParallelDescriptor::ReduceRealSum(icpt_e);
	amrex::ParallelDescriptor::ReduceRealSum(F_p);
	amrex::ParallelDescriptor::ReduceRealSum(F_visc_nn);
	amrex::ParallelDescriptor::ReduceRealSum(F_visc_nt);
        AvgRadius_ = ravg/icpt;
        AvgIntVel_ = ur_avg/icpt;
        AvgIntPres_ = p_I_avg/icpt;
	AvgIntStrain_ = e_I_avg/icpt_e;
        //std::cout<<myproc<<"\t"<< AvgRadius_ <<"\t"<< AvgIntVel_ <<"\t"<<AvgIntPres_<<"\t"<<Vol_Diff_<<"\n";

    }
    else
    {
    }

    return ;
}

void Interface::ClearInterfaceData()
{

    psi.clear();
    psi_old.clear();
    psi_prev.clear();
    source_reinit.clear();
    intercept_data.clear();
    Mask_.clear();
    Index_.clear();
    Source.clear();
    kappa_.clear();
    normal.clear();
    if (!use_FMM)
    {
        LS_Gamma.clear();
        LS_S.clear();
        // LS_C.define(mesh_.grid(), mesh_.dmap(), 4, 0);
        LS_R.clear();
        LS_Dealta.clear();
        LS_Sum.clear();
        LSReinit_Source.clear();
    }
    if(RKOrder >= 2)
    {
        FRK_psi.clear();
        RK1_psi.clear();
        RK2_psi.clear();
        RK3_psi.clear();
    }
}

void Interface::RemakeInterface(const amrexMesh &newmesh, const amrex::MultiFab &new_psi)
{
    mesh_ = newmesh;

    psi.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
    psi.setVal(0.0);
    amrex::MultiFab::Copy(psi, new_psi, 0, 0, 1, nghost); 

    psi_old.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
    psi_old.setVal(0.0);

    psi_prev.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
    psi_prev.setVal(0.0);

    source_reinit.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
    source_reinit.setVal(0.0);

    intercept_data.define(mesh_.grid(), mesh_.dmap());

    Mask_.define(mesh_.grid(), mesh_.dmap(), 1, 0);
    Mask_.setVal(0);

    Index_.define(mesh_.grid(), mesh_.dmap());
    Source.define(mesh_.grid(), mesh_.dmap());

    kappa_.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
    kappa_.setVal(0.0);

    normal.define(mesh_.grid(), mesh_.dmap(), AMREX_SPACEDIM, nghost);
    normal.setVal(0.0);

    const amrex::Real *dx = mesh_.geometry().CellSize();
    L_BETA = 4.0 * dx[0];
    L_GAMMA = 9.0 * dx[0];
    if (!use_FMM)
    {
        LS_Gamma.define(mesh_.grid(), mesh_.dmap(), 1, 0);
        LS_S.define(mesh_.grid(), mesh_.dmap(), 4, 0);
        // LS_C.define(mesh_.grid(), mesh_.dmap(), 4, 0);
        LS_R.define(mesh_.grid(), mesh_.dmap(), 1, 0);
        LS_Dealta.define(mesh_.grid(), mesh_.dmap(), 1, 0);
        LS_Sum.define(mesh_.grid(), mesh_.dmap(), 1, 0);
        LSReinit_Source.define(mesh_.grid(), mesh_.dmap(), 4, 0);

        LS_Gamma.setVal(0);
        LS_S.setVal(0);
        // LS_C.define(mesh_.grid(), mesh_.dmap(), 4, 0);
        LS_R.setVal(0);
        LS_Dealta.setVal(0);
        LS_Sum.setVal(0);
        LSReinit_Source.setVal(0);
    }
   if(RKOrder >= 2)
   {
       FRK_psi.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
       FRK_psi.setVal(0.0);
       amrex::MultiFab::Copy(FRK_psi, new_psi, 0, 0, 1, nghost);

       RK1_psi.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
       RK1_psi.setVal(0.0);
       amrex::MultiFab::Copy(RK1_psi, new_psi, 0, 0, 1, nghost);

       RK2_psi.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
       RK2_psi.setVal(0.0);
       amrex::MultiFab::Copy(RK2_psi, new_psi, 0, 0, 1, nghost);

       RK3_psi.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
       RK3_psi.setVal(0.0);
       amrex::MultiFab::Copy(RK3_psi, new_psi, 0, 0, 1, nghost);
   }

}

void Interface::MakeInterfaceFromCoarse(const amrexMesh &newmesh, const amrex::MultiFab &new_psi)
{
    //mesh_ = newmesh;

    //psi.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
    //psi.setVal(0.0);
    amrex::MultiFab::Copy(psi, new_psi, 0, 0, 1, nghost);

    //psi_old.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
    //psi_old.setVal(0.0);

    //psi_prev.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
    //psi_prev.setVal(0.0);

    //source_reinit.define(mesh_.grid(), mesh_.dmap(), 1, nghost);
    //source_reinit.setVal(0.0);

    //intercept_data.define(mesh_.grid(), mesh_.dmap());

    //Mask_.define(mesh_.grid(), mesh_.dmap(), 1, 0);
    //Mask_.setVal(0);

    //Index_.define(mesh_.grid(), mesh_.dmap());
    //Source.define(mesh_.grid(), mesh_.dmap());

    //const amrex::Real *dx = mesh_.geometry().CellSize();
    //L_BETA = 4.0 * dx[0];
    //L_GAMMA = 9.0 * dx[0];
}

void Interface::copyFRK_Psi()
{
   amrex::MultiFab::Copy(FRK_psi, psi, 0, 0, 1, nghost);
}

void Interface::TVDRK2Avg_Psi()
{
    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();

        amrex::Array4<amrex::Real> const &Psi = psi.array(mfi);
        amrex::Array4<amrex::Real const> const &frk_psi = FRK_psi.const_array(mfi);

        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Psi(i, j, k) = 0.5 * (frk_psi(i, j, k) + Psi(i, j, k));
        });
    }
}

void Interface::copyRK1_Psi()
{
   amrex::MultiFab::Copy(RK1_psi, psi, 0, 0, 1, nghost);
}

void Interface::copyRK2_Psi()
{
   amrex::MultiFab::Copy(RK2_psi, psi, 0, 0, 1, nghost);
}

void Interface::copyRK3_Psi()
{
   amrex::MultiFab::Copy(RK3_psi, psi, 0, 0, 1, nghost);
}


void Interface::RK2Avg_Psi()
{
    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();

        amrex::Array4<amrex::Real> const &Psi = psi.array(mfi);
        amrex::Array4<amrex::Real const> const &rk1_psi = RK1_psi.const_array(mfi);

        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Psi(i, j, k) = 0.5 * (rk1_psi(i, j, k) + Psi(i, j, k));
        });
    }
}

void Interface::displayMesh()
{
    //const amrex::Real *prob_lo = mesh_.geometry().ProbLo();
    //const amrex::Real *dx = mesh_.geometry().CellSize();
    amrex::Print()<<"##In displayMesh() ##"<<'\n';
    const amrex::Geometry & mesh_geom = mesh_.geometry();
    const amrex::Real *prob_lo = mesh_geom.ProbLo();
    const amrex::Real *dx = mesh_geom.CellSize();
    amrex::Print()<<"Grid = "<<mesh_.grid()<<'\n';
    amrex::Print()<<" dx = "<<dx[0]<<" , "<<dx[1]<<'\n';
    amrex::Print()<<"prob_lo = "<<prob_lo[0]<<" , "<<prob_lo[1]<<'\n';
}
void Interface::DetectUnderresolvedJetAlongAxis()
{
    normal.FillBoundary();
    const amrex::Geometry & mesh_geom = mesh_.geometry();
    const amrex::Real *prob_lo = mesh_geom.ProbLo();
    const amrex::Real *dx = mesh_geom.CellSize();
    has_jet = 0;

    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
	amrex::Array4<amrex::Real> const &Psi = psi.array(mfi);
        amrex::Array4<amrex::Real> const &Kappa = kappa_.array(mfi);
        amrex::Array4<amrex::Real> const &Normal = normal.array(mfi);
        auto &icpt_data = getInterceptData()[mfi];
        for (auto &&idt : icpt_data)
        {
	    int i = idt.cellid_[0];
	    int j = idt.cellid_[1];
	    int k = 0;
            amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
            amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
	    amrex::Real nx = Normal(i, j, k, 0);
	    amrex::Real ny = Normal(i, j, k, 1);

	    amrex::Real probe_length = 2.0*dx[1];

	    amrex::Real probe_x = x - nx*(probe_length + Psi(i, j, k));
	    amrex::Real probe_y = y - ny*(probe_length + Psi(i, j, k));

	    int i_probe = std::round((probe_x - prob_lo[0])/dx[0] - 0.5); 
	    int j_probe = std::round((probe_y - prob_lo[1])/dx[1] - 0.5);

	    if(probe_y < prob_lo[1])//i.e. probe falls below axis of symmetry
            {
	        idt.r_flag = 1;
		has_jet = 1;
	    }

            if(Psi(i_probe,j_probe,k) >= 0)//i.e. probe falls inside bubble
            {
                idt.r_flag = 2;
                has_jet = 1;
		//amrex::Print()<<" i = "<<i<<" , j = "<<j<<'\n';
                //amrex::Print()<<" x = "<<x<<" , y = "<<y<<'\n';
		//amrex::Print()<<"normal = "<<nx<<" , "<<ny<<'\n';
		//amrex::Print()<<" i_probe = "<<i_probe<<" , j_probe = "<<j_probe<<'\n';
		//amrex::Print()<<" x_probe = "<<probe_x<<" , y_probe = "<<probe_y<<'\n';
            }
        }
    }
    amrex::ParallelDescriptor::ReduceIntMax(has_jet);
}

void Interface::RemoveJetCellByCell()
{
    //amrex::Print()<<"RemoveJetCellByCell"<<'\n';
    const amrex::Geometry & mesh_geom = mesh_.geometry();
    const amrex::Real *prob_lo = mesh_geom.ProbLo();
    const amrex::Real *dx = mesh_geom.CellSize();
    const amrex::Box &domain = mesh_geom.Domain();

    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        amrex::Array4<amrex::Real> const &Psi = psi.array(mfi);
	amrex::Array4<amrex::Real> const &Psiold = psi_old.array(mfi);
        amrex::Array4<int> const &mask = Mask_.array(mfi);
        amrex::Array4<amrex::Real> const &Kappa = kappa_.array(mfi);
        amrex::Array4<amrex::Real> const &Normal = normal.array(mfi);
        auto &icpt_data = getInterceptData()[mfi];
        for (auto &&idt : icpt_data)
        {
	    if(idt.r_flag == 1)
	    {
                int i = idt.cellid_[0];
                int j = idt.cellid_[1];
                int k = 0;
                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                amrex::Real nx = Normal(i, j, k, 0);
                amrex::Real ny = Normal(i, j, k, 1);

	        //Detect the cell where interface needs to be moved
	        // Equation of line tangent to the 0 ls
                // a1*x + b1*y = c1
                amrex::Real a1 = nx;
                amrex::Real b1 = ny;
                amrex::Real c1 = nx*x + ny*y;
                // Equation of line perpendicular to the 0 ls 
                // a2*x + b2*y = c2
                amrex::Real a2 = ny;
                amrex::Real b2 = -1.0*nx;
                amrex::Real c2 = ny*x - nx*y;

	        amrex::Real y_below = prob_lo[1] + dx[1] * (std::max(j - 1,0) + 1.0);
	        amrex::Real x_below = (c2 + nx*y_below)/ny;

		int i_below = std::round((x_below - prob_lo[0])/(dx[0]) - 0.5);
		int j_below = std::max(j - 1,0);
		//amrex::Print()<<"i = "<<i<<" , j = "<<j<<'\n';
		//amrex::Print()<<"i_below = "<<i_below<<" , j_below = "<<j_below<<'\n';
                Psi(i_below, j_below, k ) = dx[1];
		//exit(1);
            }
	    else if(idt.r_flag == 2)
            {
                int i = idt.cellid_[0];
                int j = idt.cellid_[1];
                int k = 0;
                amrex::IntVect cellid_ = amrex::IntVect(AMREX_D_DECL(i, j, k));
                amrex::Box gbx(cellid_, cellid_);
                gbx.grow(4);
                amrex::Box gbx_isect = gbx & domain;
                amrex::Real avg_psi = 0.0;
                amrex::Real ncount = 0.0;
                for (amrex::BoxIterator bit(gbx_isect); bit.ok(); ++bit)
                {
                    const amrex::IntVect &iv = bit();
                    avg_psi += Psiold(iv[0],iv[1],0);
                    ncount++;
                }
                Psi(cellid_) = avg_psi/ncount;
	    }
        }
    }
}

void Interface::Regularization()
{
    const amrex::Real *dx = mesh_.geometry().CellSize();
    const amrex::Real *prob_lo = mesh_.geometry().ProbLo();
    const amrex::Real *prob_hi = mesh_.geometry().ProbHi();
    const amrex::Box &domain = mesh_.geometry().Domain();
    amrex::MultiFab::Copy(psi_old, psi, 0, 0, 1, nghost);

    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        amrex::Array4<amrex::Real> const &f = psi.array(mfi);
        amrex::Array4<amrex::Real> const &fold = psi_old.array(mfi);
        amrex::Array4<int> const &mask = Mask_.array(mfi);
        /// copy phi in psi
        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            /// create grown box around cell
	        amrex::IntVect cellid_ = amrex::IntVect(AMREX_D_DECL(i, j, k));
            amrex::Box gbx(cellid_, cellid_);
            gbx.grow(2);
            amrex::Box gbx_isect = gbx & domain;
            amrex::Real avg_psi = 0.0;
            amrex::Real ncount = 0.0;
            for (amrex::BoxIterator bit(gbx_isect); bit.ok(); ++bit)
    	    {
	            const amrex::IntVect &iv = bit();
		        avg_psi += fold(iv[0],iv[1],0);
		        ncount++;
	        } 
            f(cellid_) = avg_psi/ncount;
	    });
    }
}

void Interface::Regularization2(int MAX_ITER)
{
    MAX_ITER = std::min(100, MAX_ITER);//minimum 100 iterations
    amrex::Real C1, C2, C3, c4;

    const amrex::Real *dx = mesh_.geometry().CellSize();

    int Iter = 0;

    amrex::Real reinit_tol = 0.0;
    amrex::MultiFab::Copy(psi_prev, psi, 0, 0, 1, nghost);
    do
    {
        //amrex::Print()<<"Regilarization iter = "<<Iter<<"\n";

        const amrex::Real *prob_lo = mesh_.geometry().ProbLo();
        const amrex::Real *prob_hi = mesh_.geometry().ProbHi();

        amrex::MultiFab::Copy(psi_old, psi, 0, 0, 1, nghost);
        Compute_Normal_Curvature();
        //amrex::Print()<<"max_kappa = "<<max_kappa<<'\n';
        amrex::Real dtau = 0.05/(max_kappa*(1.0/dx[0] + 1.0/dx[1]));

        //for(int istep = 0; istep < 2; istep++)
	for(int istep = 0; istep < 3; istep++)
        {
            //amrex::Print()<<"step = "<<istep + 1<<'\n';
            //if(istep == 0)
            //{
            //    C1 = 1.0;
            //    C2 = 0.0;
            //    C3 = 1.0;
            //}
            //else if(istep == 1)
            //{
            //    C1 = 0.5;
            //    C2 = 0.5;
            //    C3 = 0.5;
            //}
	    if(istep == 0)
            {
                C1 = 1.0;
                C2 = 0.0;
                C3 = 1.0;
            }
            else if(istep == 1)
            {
                C1 = 3.0/4.0;
                C2 = 1.0/4.0;
                C3 = 1.0/4.0;
            }
	    else if(istep == 2)
            {
                C1 = 1.0/3.0;
                C2 = 2.0/3.0;
                C3 = 2.0/3.0;
            }


            PhaseFieldBC();

            compute_source_regularization();
	    Compute_Normal_Curvature();
            
            for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();
                amrex::Array4<amrex::Real> const &f = psi.array(mfi);
                amrex::Array4<amrex::Real> const &fold = psi_old.array(mfi);
                amrex::Array4<amrex::Real> const &fprev = psi_prev.array(mfi);
                amrex::Array4<amrex::Real> const &L = source_reinit.array(mfi);
		amrex::Array4<amrex::Real> const &Kappa = kappa_.array(mfi);
                amrex::Array4<int> const &mask = Mask_.array(mfi);
                /// copy phi in psi
                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    //if( f(i,j,k)*f(i + 1,j,k) > 0 &&
                    //    f(i,j,k)*f(i - 1,j,k) > 0 &&
                    //    f(i,j,k)*f(i,j + 1,k) > 0 &&
                    //    f(i,j,k)*f(i,j - 1,k) > 0 )
                    {
		        amrex::Real F = 0;

			F = std::max(Kappa(i, j, k),0.0); 
		        
                        f(i , j , k) = C1*fold(i, j, k) + 
                                       C2*f(i, j, k) +
                                       C3*dtau*L(i, j, k)*std::max(Kappa(i, j, k),0.0);
                        if(std::isnan(f(i , j , k)))
			{
			    amrex::Print()<<"Loc 2 "<<'\n';
			    amrex::Print()<<"i = "<<i<<" , j = "<<j<<" is nan "<<'\n';
			    amrex::Print()<<"fold(i, j, k) =  "<<fold(i, j, k)<<'\n';
			    amrex::Print()<<"L(i, j, k) =  "<<L(i, j, k)<<'\n';
			    amrex::Print()<<"Kappa(i, j, k) =  "<<Kappa(i, j, k)<<'\n';
			    exit(9);
			}
                        //f(i, j, k) = std::copysign(f(i, j, k),fprev(i, j, k));
                        //if(mask(i ,j ,k) == 1)
                        {
                            //amrex::Print()<<"L(i, j, k) = "<<std::fabs(L(i, j, k))<<'\n';
                            //reinit_tol = std::max(reinit_tol, std::fabs(L(i ,j ,k)));
                        }
                        //if(i == 507 && j == 531)
                        //{
                        //    amrex::Print()<<" i = "<<i<<" , j = "<<j<<'\n';
                        //    amrex::Print()<<" L( i, j, k) = "<<L(i,j,k)<<'\n';      
                        //    amrex::Print()<<"f ="<<f(i,j,k)<<" , fprev = "<<fprev(i,j,k)<<'\n';
                        //}
                    }
                    //else if(mask(i,j,k) < 3)
                    //{
                    //    f(i , j , k) = C1*fold(i, j, k) +    
                    //                   C2*f(i, j, k) -
                    //                   C3*dtau*L(i, j, k);
                    //}
                     
                            //if(i == 507 && j == 531)
                            //{
                            //    amrex::Print()<<" i = "<<i<<" , j = "<<j<<'\n';
                            //    amrex::Print()<<" L( i, j, k) = "<<L(i,j,k)<<'\n';       
                            //    amrex::Print()<<"f ="<<f(i,j,k)<<" , fprev = "<<fprev(i,j,k)<<'\n';
                            //}

                });
                //amrex::Print()<<"Max L = "<<std::max(std::fabs(L))<<"\n";
            }
        }
        Iter++;
    }
    while(Iter <= MAX_ITER );
}
void Interface::RemoveSpeckles()
{
    const amrex::Real *dx = mesh_.geometry().CellSize();
    has_spekle = 0;

    {
        //Remove speckles
        for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real> const &f = psi.array(mfi);
            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if(f(i ,j ,k) < 0)
                {   
		    if(f(i, j, k)*f(i + 1, j, k) < 0 &&
		       f(i, j, k)*f(i - 1, j, k) < 0 &&
		       f(i, j, k)*f(i, j + 1, k) < 0 &&
		       f(i, j, k)*f(i, j - 1, k) < 0)
		    {
                        f(i, j, k ) = 0.25*(f(i + 1, j, k) + f(i - 1, j, k) + f(i,j + 1, k) + f(i, j - 1, k));
			has_spekle = 1;
		    }
                    else if( f(i - 1, j, k) > 0 && f(i, j + 1, k) > 0 && f(i, j - 1, k) > 0)
		    {
		        f(i, j, k ) = -1.0*f(i + 1, j, k);//(1.0/3.0)*(f(i - 1, j, k) + f(i,j + 1, k) + f(i, j - 1, k));     
			has_spekle = 1;
		    }
		    else if( f(i + 1, j, k) > 0 && f(i, j + 1, k) > 0 && f(i, j - 1, k) > 0)
		    {
                        f(i, j, k ) = -1.0*f(i - 1, j, k);//(1.0/3.0)*(f(i + 1, j, k) + f(i,j + 1, k) + f(i, j - 1, k));
			has_spekle = 1;
		    }
                    else if( f(i + 1, j, k) > 0 && f(i - 1, j, k) > 0 && f(i, j + 1, k) > 0)
		    {
                        f(i, j, k ) = -1.0*f(i, j - 1, k);//(1.0/3.0)*(f(i + 1, j, k) + f(i - 1, j, k) + f(i,j + 1, k) );
			has_spekle = 1;
		    }
                    else if( f(i + 1, j, k) > 0 && f(i - 1, j, k) > 0 && f(i, j - 1, k) > 0)
	            {
                        f(i, j, k ) = -1.0*f(i, j + 1, k);//(1.0/3.0)*(f(i + 1, j, k) + f(i - 1, j, k) + f(i,j - 1, k) );
			has_spekle = 1;
		    }
                }
            });
        }
    }
    amrex::ParallelDescriptor::ReduceIntMax(has_spekle);
}

void Interface::compute_source_regularization()
{
    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        amrex::Array4<amrex::Real> const &f = psi.array(mfi);
        amrex::Array4<amrex::Real> const &fprev = psi_prev.array(mfi);
        amrex::Array4<amrex::Real> const &L = source_reinit.array(mfi);
        const amrex::Real *dx = mesh_.geometry().CellSize();
        /// copy phi in psi
        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            amrex::Real sgn_f = std::copysign(1.0,fprev(i,j,k));

            amrex::Real f_xx_p = (f(i + 1, j, k) - 2.0*f(i, j, k) + f(i - 1, j, k))/(dx[0]*dx[0]);
            amrex::Real f_xx_e = (f(i + 2, j, k) - 2.0*f(i + 1, j, k) + f(i, j, k))/(dx[0]*dx[0]);
            amrex::Real f_xx_w = (f(i, j, k) - 2.0*f(i - 1, j, k) + f(i - 2, j, k))/(dx[0]*dx[0]);
            amrex::Real f_x_p = (f(i + 1, j, k) - f(i, j, k))/dx[0] - 0.5*dx[0]*minmod(f_xx_p, f_xx_e);
            amrex::Real f_x_m = (f(i, j, k) - f(i - 1, j, k))/dx[0] + 0.5*dx[0]*minmod(f_xx_p, f_xx_w);

            amrex::Real f_yy_p = (f(i, j + 1, k) - 2.0*f(i, j, k) + f(i, j - 1, k))/(dx[1]*dx[1]);
            amrex::Real f_yy_n = (f(i, j + 2, k) - 2.0*f(i, j + 1, k) + f(i, j, k))/(dx[1]*dx[1]);
            amrex::Real f_yy_s = (f(i, j, k) - 2.0*f(i, j - 1, k) + f(i, j - 2, k))/(dx[1]*dx[1]);
            amrex::Real f_y_p = (f(i , j + 1, k) - f(i, j, k))/dx[1] - 0.5*dx[1]*minmod(f_yy_p, f_yy_n);
            amrex::Real f_y_m = (f(i , j, k) - f(i, j - 1, k))/dx[1] + 0.5*dx[1]*minmod(f_yy_p, f_yy_s);

            amrex::Real H_ijk = 0.0;
            //amrex::Print()<<"f("<<i<<" , "<<j<<" , "<<k << ") = "<<f(i, j, k)<<" , sign = "<<sgn_f<<'\n';
            amrex::Real a_p = std::max(f_x_p,0.0);
            amrex::Real a_m = std::min(f_x_p,0.0); 
            amrex::Real b_p = std::max(f_x_m,0.0);
            amrex::Real b_m = std::min(f_x_m,0.0);
            amrex::Real c_p = std::max(f_y_p,0.0);
            amrex::Real c_m = std::min(f_y_p,0.0);
            amrex::Real d_p = std::max(f_y_m,0.0);
            amrex::Real d_m = std::min(f_y_m,0.0);
            if(sgn_f >= 0 )
                H_ijk = std::sqrt(
                           std::max(a_m*a_m,b_p*b_p) +
                           std::max(c_m*c_m,d_p*d_p));
            else
                H_ijk = std::sqrt(
                           std::max(a_p*a_p,b_m*b_m) +
                           std::max(c_p*c_p,d_m*d_m));
            L(i,j,k) = H_ijk; 
            //if(i == 127 && j == 141)
            //{
            //    amrex::Print()<<" i = "<<i<<" , j = "<<j<<'\n';
            //    amrex::Print()<<" f = "<<f(i ,j ,k)<<" ,f_e = "<<f(i+1,j,k)<<" , f_w = "<<f(i-1,j,k)<<" , f_n = "<<f(i, j+1, k)<<" , f_s = "<<f(i , j-1,k)<<'\n';
            //    amrex::Print()<<" L( i, j, k) = "<<L(i,j,k)<<" , sgn_f = "<<sgn_f<<'\n'; 
            //}

        });
    }
}

void Interface::Compute_Normal_Curvature()
{
    max_kappa = 0.0;
    max_kappa_interface = 0.0;
    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        const amrex::Real *prob_lo = mesh_.geometry().ProbLo();
        const amrex::Real *dx = mesh_.geometry().CellSize();

        amrex::Array4<amrex::Real> const &Psi = psi.array(mfi);
        amrex::Array4<amrex::Real> const &Kappa = kappa_.array(mfi);
        amrex::Array4<amrex::Real> const &Normal = normal.array(mfi);
        amrex::Array4<int> const &mask = Mask_.array(mfi);
        amrex::Real deltax = dx[1];
        /// copy phi in psi
        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
	    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
	    Kappa(i, j, k) = 0.0;
            if(mask(i, j, k) > 2)
            {
                amrex::Real Psi_xx, Psi_xy, Psi_yy,dPsi_x,dPsi_y ;
                dPsi_x = 0.5 * (Psi(i + 1, j, k) - Psi(i - 1, j, k));
                dPsi_y = 0.5 * (Psi(i, j + 1, k) - Psi(i, j - 1, k));
                Psi_xx = Psi(i + 1, j, k) - 2.0 * Psi(i, j, k) + Psi(i - 1, j, k);
                Psi_yy = Psi(i, j + 1, k) - 2.0 * Psi(i, j, k) + Psi(i, j - 1, k);
                Psi_xy = 0.25 * (Psi(i + 1, j + 1, k) - Psi(i + 1, j - 1, k) - Psi(i - 1, j + 1, k) + Psi(i - 1, j - 1, k));

                Kappa(i, j, k) = (Psi_xx * (dPsi_y * dPsi_y) + Psi_yy * (dPsi_x * dPsi_x) - 2.0 * dPsi_x * dPsi_y * Psi_xy) / (deltax * pow(dPsi_x * dPsi_x + dPsi_y * dPsi_y + 1.0E-14, 1.5));
                if (isAxisymmetric)
                {
                    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                    //Kappa(i, j, k) += (1.0 / y) * dPsi_y / sqrt(dPsi_x * dPsi_x + dPsi_y * dPsi_y + 1.0E-14);
                }

	        	max_kappa = std::max(max_kappa,std::abs(Kappa(i, j , k)));
	        	if((Psi(i + 1, j, k)*Psi(i, j, k) < 0.0 ||
                Psi(i - 1, j, k)*Psi(i, j, k) < 0.0 ||
		        Psi(i, j, k)*Psi(i, j + 1, k) < 0.0 ||
		        Psi(i, j, k)*Psi(i, j - 1, k) < 0.0) && j > 2)
	            {
		            max_kappa_interface = std::max(max_kappa_interface,Kappa(i, j , k));
                    if(Kappa(i, j , k) > 10.0)
                    {
                        //amrex::Print()<<" i = "<<i<<" , j = "<<j<<" , kappa = " << Kappa(i, j , k)<<'\n';
                    }
                }
	            amrex::Real mod_grad_phi = std::hypot(dPsi_x,dPsi_y) + 1.0E-14 ;	
                Normal(i, j, k, 0) = dPsi_x/mod_grad_phi;
                Normal(i, j, k, 1) = dPsi_y/mod_grad_phi;
                if(std::isnan(Kappa(i,j,k)))
                {
                    amrex::Print()<<" i = "<<i<<" , j = "<<j<<'\n';
                    amrex::Print()<<"Kappa(i, j, k) ="<<Kappa(i, j, k)<<'\n';
                    amrex::Print()<<"dPsi_x = "<<dPsi_x<<" , dPsi_y = "<<dPsi_y<<" , Psi_xx = "<<Psi_xx<<" , Psi_yy = "<<Psi_yy<<" , Psi_xy = "<<Psi_xy<<'\n';
                    amrex::Print()<<"Normal = "<<Normal(i, j, k, 0)<<" , "<<Normal(i, j, k, 1)<<'\n'; 
                    Kappa(i, j, k) = 1e9;
// exit(9);
                }
    	    }
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(max_kappa);
    amrex::ParallelDescriptor::ReduceRealMax(max_kappa_interface);
    //amrex::Print()<<"max_kappa_interface = "<<max_kappa_interface<<'\n';
}


} /*End namespace mycode */

