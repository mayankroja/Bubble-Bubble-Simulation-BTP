#include "incFSI.H"
#include <CFMask.H>
#include <AMReX_PlotFileUtil.H>
#include <SolveCubicEqn.h>
#include <AdvectLS_helper.H>

namespace mycode
{

void incFSI::WriteFile()
{
    int ncomp = 6;
    if(!DamageModel)
    {
        plot_strain_error = false;
    }
    if(plot_strain_rate) ncomp++;
    if(N_IF > 0) ncomp = ncomp + N_IF;
    if(TempField) ncomp++;
    if(PhaseField) ncomp++;
    if(DamageModel) ncomp=ncomp+7;
    if(plot_strain_error) ncomp=ncomp+4;
    if(plot_viscosity) ncomp++;

    const std::string& pltfile = amrex::Concatenate("Output/plt",Iter,5);
    
    amrex::Vector<amrex::MultiFab> plotmf(finest_level + 1);
    amrex::Vector<amrex::MultiFab> Pressure_plot(finest_level + 1);
    amrex::Vector<amrex::MultiFab> Temperature_plot(finest_level + 1);
    amrex::Vector<amrex::MultiFab> Phi_plot(finest_level + 1);
    amrex::Vector<amrex::MultiFab> dmg_plot(finest_level + 1);
    amrex::Vector<amrex::MultiFab> strain_rate_plot(finest_level + 1);
    amrex::Vector<amrex::MultiFab> strain_exact(finest_level + 1);//for axisymmetric bubble oscillation only
    amrex::Vector<amrex::MultiFab> strain_error(finest_level + 1);//for axisymmetric bubble oscillation only
    amrex::Vector<amrex::MultiFab> viscosity(finest_level + 1);
    amrex::Vector<amrex::MultiFab> div_F_r(finest_level + 1);
    amrex::Vector<amrex::MultiFab> div_F_z(finest_level + 1);
    for (int lev = 0; lev <= finest_level; lev++)
    {
	const amrex::Real *dx = geom[lev].CellSize();
        const amrex::Real *prob_lo = geom[lev].ProbLo();
        const amrex::Real *prob_hi = geom[lev].ProbHi();
        amrex::iMultiFab *PMask;
        amrex::iMultiFab PMask_NI;//NI stands for No Interface, need a better way of doing this
        if(N_IF > 0)
        {
            auto &&mask__= mask_[lev];
            PMask = &mask__->getPMask();
        }
        else
        {
            PMask_NI.define(grids[lev], dmap[lev], 1, Nghost);
            PMask_NI.setVal(1);
        }

        Pressure_plot[lev].define(grids[lev], dmap[lev], 1 , 0);
        amrex::MultiFab::Copy(Pressure_plot[lev], Pressure[lev],  0, 0, 1, 0);
        for(amrex::MFIter mfi(Pressure_plot[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox(); 
            amrex::Array4<amrex::Real> const  &P = Pressure_plot[lev].array(mfi);
            amrex::Array4<int const> pmask;
            if(N_IF > 0)
                pmask = PMask->const_array(mfi);
            else
                pmask = PMask_NI.const_array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if(N_IF > 0)
                {
                    if(pmask(i , j , k) != 1)
                    {
                        P(i,j,k) = 0.0;
                        if(N_IF == 1)
                        {
                            auto &solid = interfaces[finest_level][0];
                            P(i,j,k) = solid->getP_interface();
                        }
                    }
                }
            });
        }
        if(TempField)
        {
            Temperature_plot[lev].define(grids[lev], dmap[lev], 1 , 0);
            amrex::MultiFab::Copy(Temperature_plot[lev], Theta[lev],  0, 0, 1, 0);
            for(amrex::MFIter mfi(Pressure_plot[lev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<amrex::Real> const  &T = Temperature_plot[lev].array(mfi);
                amrex::Array4<int const> pmask;
                if(N_IF > 0)
                    pmask = PMask->const_array(mfi);
                else
                    pmask = PMask_NI.const_array(mfi);

                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if(N_IF > 0)
                    {
                        if(pmask(i , j , k) != 1)
                            T(i,j,k) = 0.0;
                    }
                });
            }
        }
        if(PhaseField)
        {
            Phi_plot[lev].define(grids[lev], dmap[lev], 1 , 0);
            amrex::MultiFab::Copy(Phi_plot[lev], Phi[lev],  0, 0, 1, 0);
            for(amrex::MFIter mfi(Phi_plot[lev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<amrex::Real> const  &phi = Phi_plot[lev].array(mfi);
                amrex::Array4<int const> pmask;
                if(N_IF > 0)
                    pmask = PMask->const_array(mfi);
                else
                    pmask = PMask_NI.const_array(mfi);

                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if(N_IF > 0)
                    {
                        if(pmask(i , j , k) != 1)
                            phi(i,j,k) = -1.0;
                    }

                });
            }
        }
	if(DamageModel)
        {
            dmg_plot[lev].define(grids[lev], dmap[lev], 1 , 0);
            amrex::MultiFab::Copy(dmg_plot[lev], Scalars[lev][8],  0, 0, 1, 0);
            for(amrex::MFIter mfi(dmg_plot[lev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<amrex::Real> const  &dmg = dmg_plot[lev].array(mfi);
                amrex::Array4<int const> pmask;
                if(N_IF > 0)
                    pmask = PMask->const_array(mfi);
                else
                    pmask = PMask_NI.const_array(mfi);

                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if(N_IF > 0)
                    {
                        if(pmask(i , j , k) != 1)
                            dmg(i,j,k) = 0.0;
                    }
                });
            }
        }
	if(plot_strain_rate)
        {
	    const amrex::Real *dx = geom[lev].CellSize();
            strain_rate_plot[lev].define(grids[lev], dmap[lev], 1 , 0);
            for(amrex::MFIter mfi(U[lev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<amrex::Real> const  &strain_rate = strain_rate_plot[lev].array(mfi);
		amrex::Array4<amrex::Real const> const &u = xvel[lev].const_array(mfi);
                amrex::Array4<amrex::Real const> const &v = yvel[lev].const_array(mfi);
                amrex::Array4<amrex::Real const> e_max;
                if(DamageModel)
                    e_max = Scalars[lev][7].const_array(mfi);

                amrex::Array4<int const> pmask;
                if(N_IF > 0)
                    pmask = PMask->const_array(mfi);
                else
                    pmask = PMask_NI.const_array(mfi);

                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);

                    amrex::Real ux = (u(i + 1, j, k) - u(i, j, k)) / dx[0];
                    amrex::Real uy = (u(i + 1, j + 1, k) - u(i + 1, j - 1, k) + u(i, j + 1, k) - u(i, j - 1, k)) / (4.0 * dx[1]);
                    amrex::Real vx = (v(i + 1, j + 1, k) - v(i - 1, j + 1, k) + v(i + 1, j, k) - v(i - 1, j, k)) / (4.0 * dx[0]);
                    amrex::Real vy = (v(i, j + 1, k) - v(i, j, k)) / dx[1];
		    strain_rate(i, j, k) = visc_.ComputeGammaDot(ux, 0.5 * (uy + vx), vy);
                    if(N_IF > 0)
                    {
                        if(pmask(i , j , k) != 1)
                            strain_rate(i,j,k) = 0.0;
                    }
                });
            }
        }
        if(plot_strain_error)
        {
            auto &solid = interfaces[finest_level][0];
	    amrex::Real bubble_radius = solid->AvgRadius();
            const amrex::Real *dx = geom[lev].CellSize();
            strain_rate_plot[lev].define(grids[lev], dmap[lev], 1 , 0);
	    strain_exact[lev].define(grids[lev], dmap[lev], 1 , 0);
	    strain_error[lev].define(grids[lev], dmap[lev], 1 , 0);
            div_F_r[lev].define(grids[lev], dmap[lev], 1 , 0);
            div_F_z[lev].define(grids[lev], dmap[lev], 1 , 0);
            for(amrex::MFIter mfi(U[lev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<amrex::Real> const  &Strain_exact = strain_exact[lev].array(mfi);
		amrex::Array4<amrex::Real> const  &Strain_error = strain_error[lev].array(mfi);
                amrex::Array4<amrex::Real> const  &divFr = div_F_r[lev].array(mfi);
                amrex::Array4<amrex::Real> const  &divFz = div_F_z[lev].array(mfi); 
                amrex::Array4<amrex::Real const> e_max, F11,F12,F21,F22,F33;
                if(DamageModel)
                {
                    F11 = Scalars[lev][2].const_array(mfi);
                    F12 = Scalars[lev][3].const_array(mfi);
                    F22 = Scalars[lev][4].const_array(mfi);
                    F21 = Scalars[lev][5].const_array(mfi);
                    F33 = Scalars[lev][6].const_array(mfi);
                    e_max = Scalars[lev][7].const_array(mfi);
                }

                amrex::Array4<int const> pmask;
                if(N_IF > 0)
                    pmask = PMask->const_array(mfi);
                else
                    pmask = PMask_NI.const_array(mfi);

                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
		    Strain_exact(i ,j ,k) = 0.0;
		    Strain_error(i, j, k) = 0.0;
                    amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                    amrex::Real r = std::hypot(x-solid->Xcp(), y-solid->Ycp());
                    amrex::Real r_0 = std::pow(r*r*r + (0.5*0.5*0.5 - bubble_radius*bubble_radius*bubble_radius),1.0/3.0);
                    if(pmask(i, j, k) == 1)
                        Strain_exact(i ,j ,k) = (r*r - r_0*r_0)/(2.0*r_0*r_0);

                    if(Strain_exact(i ,j ,k) != 0.0 && pmask(i, j, k) == 1)
                        Strain_error(i ,j ,k) = std::abs((Strain_exact(i ,j ,k) - e_max(i, j, k)));

                    if(DamageModel)
                    {
                        divFr(i, j, k) = (F22(i, j + 1, k) - F22(i, j - 1, k))/(2.0 * dx[1]) + 
                                         (F12(i + 1, j, k) - F21(i - 1, j, k))/(2.0 * dx[0]) + 
                                          F22(i, j, k)/y - 3.0 * F33(i, j, k)/y;
                        divFz(i, j, k) = (F21(i, j + 1, k) - F21(i, j - 1, k))/(2.0 * dx[1]) +
                                         (F11(i + 1, j, k) - F11(i - 1, j, k))/(2.0 * dx[0]) +  
                                          F21(i, j, k)/y;
         
                    } 
                });
            }
        }
	if(plot_viscosity)
	{
            const amrex::Real *dx = geom[lev].CellSize();
            viscosity[lev].define(grids[lev], dmap[lev], 1 , 0);
            for(amrex::MFIter mfi(U[lev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<amrex::Real const> const &u = xvel[lev].const_array(mfi);
                amrex::Array4<amrex::Real const> const &v = yvel[lev].const_array(mfi);
                amrex::Array4<amrex::Real> const  &viscosity_ = viscosity[lev].array(mfi);
                amrex::Array4<amrex::Real const> phi, dmg;
                if(PhaseField)
                    phi = Phi[lev].array(mfi);
                if(DamageModel)
                    dmg = Scalars[lev][8].array(mfi);

                amrex::Array4<int const> pmask;
                if(N_IF > 0)
                    pmask = PMask->const_array(mfi);
                else
                    pmask = PMask_NI.const_array(mfi);

                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);

                    amrex::Real ux = (u(i + 1, j, k) - u(i, j, k)) / dx[0];
                    amrex::Real uy = (u(i + 1, j + 1, k) - u(i + 1, j - 1, k) + u(i, j + 1, k) - u(i, j - 1, k)) / (4.0 * dx[1]);
                    amrex::Real vx = (v(i + 1, j + 1, k) - v(i - 1, j + 1, k) + v(i + 1, j, k) - v(i - 1, j, k)) / (4.0 * dx[0]);
                    amrex::Real vy = (v(i, j + 1, k) - v(i, j, k)) / dx[1];

                    amrex::Real S11_ = ux;
                    amrex::Real S12_ = 0.5*(uy + vx);
                    amrex::Real S22_ = vy;

                    amrex::Real phi_, dmg_;
                    if(PhaseField)
                        phi_ = phi(i, j, k);
                    else
                        phi_ = 0.0;
                    if(DamageModel)
                        dmg_ = dmg(i, j, k);
                    else
                        dmg_ = 0.0;
                    
		    viscosity_(i, j, k) = visc_.GetViscosity(S11_, S12_, S22_, phi_, dmg_);
                    if(N_IF > 0)
                    {
                        if(pmask(i , j , k) != 1)
                            viscosity_(i, j, k) = 0.0;
                    }
                });
            }
	}



        int nvar = 0;
        plotmf[lev].define(grids[lev], dmap[lev], ncomp, 0);
        amrex::MultiFab::Copy(plotmf[lev], U[lev],         0, nvar, 1, 0);
        amrex::MultiFab::Copy(plotmf[lev], U[lev],         1, ++nvar, 1, 0);
        amrex::MultiFab::Copy(plotmf[lev], Pressure_plot[lev],  0, ++nvar, 1, 0);
        amrex::MultiFab::Copy(plotmf[lev], Umag[lev],      0, ++nvar, 1, 0);
        amrex::MultiFab::Copy(plotmf[lev], Vorticity[lev], 0, ++nvar, 1, 0);

        const amrex::Real *dxinv = geom[lev].InvCellSize();
        for (amrex::MFIter mfi(plotmf[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real> const &div = plotmf[lev].array(mfi);
            amrex::Array4<amrex::Real const> const &u = xvel[lev].const_array(mfi);
            amrex::Array4<amrex::Real const> const &v = yvel[lev].const_array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
	        amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
		amrex::Real ys = prob_lo[1] + dx[1] * (j );
		amrex::Real yn = prob_lo[1] + dx[1] * (j + 1.0);

	        if(isAxisymmetric)
		{
		    div(i, j, k, 5) = std::abs(dxinv[0] * (u(i + 1, j, k) - u(i, j, k)) + (1.0/y)*dxinv[1] * (yn*v(i, j + 1, k) - ys*v(i, j, k)));
		}
		else
		{
                    div(i, j, k, 5) = dxinv[0] * (u(i + 1, j, k) - u(i, j, k)) + dxinv[1] * (v(i, j + 1, k) - v(i, j, k));
		}
            });
        }
	++nvar;

	if(plot_strain_rate) amrex::MultiFab::Copy(plotmf[lev], strain_rate_plot[lev], 0, ++nvar, 1, 0);

        int nsolid = 0;        
        for (auto &&solid : interfaces[lev])
        {
            nsolid++;
            amrex::MultiFab::Copy(plotmf[lev], solid->Psi(), 0, ++nvar, 1, 0);
        }
        if(TempField)amrex::MultiFab::Copy(plotmf[lev], Temperature_plot[lev], 0, ++nvar, 1, 0);
        if(PhaseField)amrex::MultiFab::Copy(plotmf[lev], Phi_plot[lev], 0, ++nvar, 1, 0);
	if(DamageModel)
        {
	    for(int iscalar = 2; iscalar <= 6; iscalar++)
            {
	        amrex::MultiFab::Copy(plotmf[lev], Scalars[lev][iscalar], 0, ++nvar, 1, 0);
	    }
	    amrex::MultiFab::Copy(plotmf[lev], Scalars[lev][7], 0, ++nvar, 1, 0);
	    amrex::MultiFab::Copy(plotmf[lev], dmg_plot[lev], 0, ++nvar, 1, 0);
	}
        if(plot_strain_error)
        {
            amrex::MultiFab::Copy(plotmf[lev], strain_exact[lev], 0, ++nvar, 1, 0);
            amrex::MultiFab::Copy(plotmf[lev], strain_error[lev], 0, ++nvar, 1, 0);
            amrex::MultiFab::Copy(plotmf[lev], div_F_r[lev], 0, ++nvar, 1, 0);
            amrex::MultiFab::Copy(plotmf[lev], div_F_z[lev], 0, ++nvar, 1, 0);
        }
	if(plot_viscosity)
	    amrex::MultiFab::Copy(plotmf[lev], viscosity[lev], 0, ++nvar, 1, 0);


    }
            
    amrex::Vector<std::string> varname = {"u", "v", "P", "velmag", "omega", "divU"};   
    if(plot_strain_rate)
        varname.push_back("strain rate");
    if(N_IF > 0)
    {
        for (auto &&solid : interfaces[finest_level])
        {
            varname.push_back("Psi");
        }
    }
    if(TempField)
        varname.push_back("Temperature");
    if(PhaseField)
        varname.push_back("Phi");
    if(DamageModel)
    {
        //varname.push_back("X");
        //varname.push_back("Y");	
	varname.push_back("F11");
	varname.push_back("F12");
        varname.push_back("F21");
	varname.push_back("F22");
	varname.push_back("F33");
	//varname.push_back("eps_p_11");
	varname.push_back("max_e");
	varname.push_back("Damage");
    }
    if(plot_strain_error)
    {
        varname.push_back("e_exact");
        varname.push_back("e_error"); 
        varname.push_back("div_F_r");
        varname.push_back("div_F_z");
    }
    if(plot_viscosity)
        varname.push_back("viscosity");

    if(DamageModel) 
        amrex::VisMF::SetNOutFiles(8);
    else
       amrex::VisMF::SetNOutFiles(4);

    amrex::WriteMultiLevelPlotfile(pltfile, finest_level+1, amrex::GetVecOfConstPtrs(plotmf),
                                   varname, geom, Time, amrex::Vector<int>(finest_level+1, Iter),
                                   ref_ratio);
}

void incFSI::WriteFileParaview()
{
    int ncomp = 3;

    const std::string& pltfile = amrex::Concatenate("Output/plt",Iter,5);
    
    amrex::Vector<amrex::MultiFab> plotmf(finest_level + 1);
    for (int lev = 0; lev <= finest_level; lev++)
    {
        const amrex::Real *dx = geom[lev].CellSize();
        const amrex::Real *prob_lo = geom[lev].ProbLo();
        const amrex::Real *prob_hi = geom[lev].ProbHi();
        

        int nvar = 0;
        plotmf[lev].define(grids[lev], dmap[lev], ncomp, 0);
        amrex::MultiFab::Copy(plotmf[lev], U[lev],         0, nvar, 1, 0);
        amrex::MultiFab::Copy(plotmf[lev], U[lev],         1, ++nvar, 1, 0);
        amrex::MultiFab::Copy(plotmf[lev], Pressure[lev],  0, ++nvar, 1, 0);
    }        
    amrex::Vector<std::string> varname = {"u", "v", "P"};   
 

    amrex::WriteMultiLevelPlotfile(pltfile, finest_level+1, amrex::GetVecOfConstPtrs(plotmf),
                                   varname, geom, Time, amrex::Vector<int>(finest_level+1, Iter),
                                   ref_ratio);
}

void incFSI::Derive()
{
    //! on all levels compute vorticity and vel mag
    MaxU = 0.0;
    for (int lev = 0; lev <= finest_level; lev++)
    {
        const amrex::Real *dx = geom[lev].CellSize();
        Vorticity[lev].setVal(0.0);

        amrex::iMultiFab *PMask;
        amrex::iMultiFab PMask_NI;//NI stands for No Interface, need a better way of doing this
        if(N_IF > 0)
        {
            auto &&mask__= mask_[lev];
            PMask = &mask__->getPMask();
        }
        else
        {
            PMask_NI.define(grids[lev], dmap[lev], 1, Nghost);
            PMask_NI.setVal(1);
        }

        for (amrex::MFIter mfi(Umag[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real const> const &vel = U[lev].const_array(mfi);
            amrex::Array4<amrex::Real const> const &u = xvel[lev].const_array(mfi);
            amrex::Array4<amrex::Real const> const &v = yvel[lev].const_array(mfi);
            amrex::Array4<amrex::Real> const &vel_mag = Umag[lev].array(mfi);
            amrex::Array4<amrex::Real> const &omega = Vorticity[lev].array(mfi);

            amrex::Array4<int const> pmask;
            if(N_IF > 0)
                pmask = PMask->const_array(mfi);
            else
                pmask = PMask_NI.const_array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                vel_mag(i, j, k) = std::hypot(vel(i, j, k, 0), vel(i, j, k, 1));
                MaxU = std::max(MaxU, vel_mag(i, j, k));

                // vorticity computation
                //ux = (u(i + 1, j, k) - u(i, j, k)) / dx[0];
                amrex::Real uy = 0.25 * (u(i, j + 1, k) - u(i, j - 1, k)) / dx[1] + 0.25 * (u(i + 1, j + 1, k) - u(i, j - 1, k)) / dx[1];
                amrex::Real vx = 0.25 * (v(i + 1, j + 1, k) - v(i - 1, j + 1, k)) / dx[0] + 0.25 * (v(i + 1, j, k) - v(i - 1, j, k)) / dx[0];
                //amrex::Real uy = (vel(i,j + 1, k,0) - vel(i,j-1,k,0))/(2.0*dx[1]);
                //amrex::Real vx = (vel(i + 1,j , k, 1) - vel(i - 1,j,k, 1 ))/(2.0*dx[0]);
                //vy = (v(i, j + 1, k) - v(i, j, k)) / dx[1];
                omega(i, j, k) = vx - uy;
                if(std::isnan(omega(i, j, k)))
                {
                    //amrex::Print()<<"lev = "<<lev<<", Location :"<<i <<" , "<<j<<" , "<<k<<'\n';
                    //amrex::Print()<<"vx = "<<vx<<" , uy = "<<uy<<'\n';
                    
                }
		        if(N_IF > 0 && pmask(i,j,k) != 1)
		    	omega(i, j, k) = 0.0;
            });
        }
    }

    amrex::ParallelDescriptor::ReduceRealMax(MaxU);
}

void incFSI::CheckConservationPhaseField()
{
    //! on all levels compute vorticity and vel mag
    amrex::Real tot_mass = 0.0;
    amrex::Real tot_x_mom = 0.0;
    amrex::Real tot_y_mom = 0.0;
    amrex::Real max_abs_u = 0.0;
    amrex::Real max_abs_v = 0.0;
    amrex::Real damage_front_left = 1e9;
    amrex::Real damage_front_right = -1e9;
    for (int lev = 0; lev <= finest_level; lev++)
    {
        amrex::Real tot_mass_lev = 0.0;
        amrex::Real tot_x_mom_lev = 0.0;
        amrex::Real tot_y_mom_lev = 0.0;
        amrex::Real max_abs_u_lev = 0.0;
        amrex::Real max_abs_v_lev = 0.0;
        amrex::Real damage_front_left_lev = 1e9;
        amrex::Real damage_front_right_lev = -1e9;

        CFMask cfmask_p_;
        cfmask_p_.define(lev, geom, grids, dmap);

	    amrex::iMultiFab *PMask;
        amrex::iMultiFab PMask_NI;//NI stands for No Interface, need a better way of doing this
        if(N_IF > 0)
        {
            auto &&mask__= mask_[lev];
            PMask = &mask__->getPMask();
        }
        else
        {
            PMask_NI.define(grids[lev], dmap[lev], 1, Nghost);
            PMask_NI.setVal(1);
        }

	    auto &solid = interfaces[finest_level][0];
	    amrex::MultiFab *psi;
	    if(N_IF > 0)psi = &solid->Psi();


        const amrex::Real *dx = geom[lev].CellSize();
	    const amrex::Real *prob_lo = geom[lev].ProbLo();
        amrex::Real grid_const = std::sqrt(2.0*4*dx[0]*dx[1]);

        int found_left_edge = 0;
        int found_right_edge = 0;

	
        for (amrex::MFIter mfi(Phi[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real> const &phi = Phi[lev].array(mfi);
            amrex::Array4<amrex::Real const> const &Vel = U[lev].const_array(mfi);
            amrex::Array4<int const> const &cfmask = cfmask_p_.Mask().const_array(mfi);
            amrex::Array4<amrex::Real>  Psi;
            if(lev == finest_level && N_IF > 0)Psi = psi->array(mfi);

            amrex::Array4<int const> pmask;
            if(N_IF > 0)
                pmask = PMask->const_array(mfi);
            else
                pmask = PMask_NI.const_array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
	            amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                bool interfacial_pt = false;
                 //if(pmask(i, j, k) == 1 && cfmask(i,j,k) != CFMask::covered)
                //if(cfmask(i,j,k) != CFMask::covered)
                if(pmask(i, j, k) == 1 && cfmask(i,j,k) != CFMask::covered)
                {
		            /*
		            if(lev == finest_level)
		            {
		                if((Psi(i, j, k) * Psi(i + 1, j, k) < 0 ||
	                            Psi(i, j, k) * Psi(i - 1, j, k) < 0 ||
		                    Psi(i, j, k) * Psi(i, j + 1, k) < 0 ||
		                    Psi(i, j, k) * Psi(i, j - 1, k) < 0 ))
		                    interfacial_pt  = true;
		            }
		            if(!interfacial_pt)
		                if(lev != finest_level || (lev == finest_level && pmask(i, j, k) == 1))
		                    tot_mass_lev += (1.0 - phi(i,j,k))*dx[0]*dx[1]*2*M_PI*std::abs(y);
	                */
                    if (isAxisymmetric)
                    {
                        tot_mass_lev += 0.5*(1.0 - phi(i,j,k))*dx[0]*dx[1]*2*M_PI*std::abs(y);
                        tot_x_mom_lev += 0.5*(1.0 - phi(i,j,k))*dx[0]*dx[1]*2*M_PI*std::abs(y)*Vel(i,j,k,0);
                    }
		            else
                    {
	                    tot_mass_lev += 0.5*(1.0 - phi(i,j,k))*dx[0]*dx[1];
                        tot_x_mom_lev += 0.5*(1.0 - phi(i,j,k))*dx[0]*dx[1]*Vel(i,j,k,0);
                    }
                    if(std::abs(Vel(i,j,k,0)) > std::abs(max_abs_u_lev))
                    {
                        max_abs_u_lev = std::abs(Vel(i,j,k,0));
                    }
		        }
		        if(pmask(i, j, k) == 1 && lev == finest_level)
		        {	
	                //int found_left_edge = 0;
		            //int found_right_edge = 0;
		            if(phi(i, j, k) >= 0)
	                {
		        	    if(Psi(i,j,k)*Psi(i + 1,j,k) <= 0)
	                        damage_front_left_lev = std::min(damage_front_left_lev, x - Psi(i,j,k));
		            }
                    else
		            {
	                    if(phi(i - 1,j ,k) >=0 && phi(i, j, k) < 0)
	                    {
	                        found_left_edge = 1;
		                    damage_front_left_lev = std::min(damage_front_left_lev, x + grid_const * std::atanh(phi(i, j, k)));
		                }
		            }
        
		            if(phi(i, j, k) >= 0)
		            {
		        	    if(Psi(i - 1,j,k)*Psi(i,j,k) <= 0)
                            damage_front_right_lev = std::max(damage_front_right_lev, x - Psi(i,j,k));
                    }
                    else
		            {
		               if(phi(i + 1,j ,k) >=0 && phi(i, j, k) < 0)
                       {
		                    found_right_edge = 1;
                            damage_front_right_lev = std::max(damage_front_right_lev, x + grid_const * std::atanh(phi(i, j, k)));
                       }
		            }
		        }
            });
        }

        tot_mass += tot_mass_lev;
        tot_x_mom += tot_x_mom_lev;
        if(std::abs(max_abs_u_lev) > std::abs(max_abs_u))
            max_abs_u = max_abs_u_lev;
	    if(lev == finest_level)
	    {
	        damage_front_left = damage_front_left_lev;
            damage_front_right = damage_front_right_lev;
	    }
    }

    amrex::ParallelDescriptor::ReduceRealSum(tot_mass);
    amrex::ParallelDescriptor::ReduceRealSum(tot_x_mom);
    amrex::ParallelDescriptor::ReduceRealMax(max_abs_u);
    amrex::ParallelDescriptor::ReduceRealMax(damage_front_right);
    amrex::ParallelDescriptor::ReduceRealMin(damage_front_left);
    //amrex::Print()<<"tot_mass = "<<tot_mass<<'\n';
    //amrex::Print()<<"damage_front_right = "<<damage_front_right<<'\n';
    //amrex::Print()<<"damage_front_left = "<<damage_front_left<<'\n';

    //For now, hard coded for the first bubble
    {
        char *s, *s1, *s2;
        s = new char[50];
        s1 = new char[50];
        s2 = new char[50];
        sprintf(s, "%2.10f", float(Time));
        strcpy(s1, "Output/Mass.dat");
        //strcat(s1, s);
        // std::string filename(s1);
        std::ofstream ofs_(s1, std::ofstream::out | std::ofstream::app);
        ofs_.flags(std::ios::dec | std::ios::scientific);
        ofs_.precision(7);
        //for (auto &&solid : interfaces)
        //{
            auto &solid = interfaces[finest_level][0];
            amrex::Print(ofs_) << Time <<"\t"<< tot_mass <<"\t"<< damage_front_left<<"\t"<< damage_front_right<<"\t"<<tot_x_mom<<"\t"<<max_abs_u<<"\n";
        //}
        ofs_.close();
    }
    //exit(7);
}

void incFSI::writeDataFile()
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        std::ofstream iofile("Data.dat", std::ofstream::out | std::ofstream::app);
        iofile.flags(std::ofstream::dec | std::ofstream::scientific);
        iofile.precision(8);
        if (!iofile)
        {
            std::cerr << "Error: Output file for time trace could not be opened.\n";
            exit(1);
        }

        iofile << std::endl
               << Time << "\t" << poisson_iter << "\t" << system_size << "\t" << poisson_tol << "\t"
               << dt << "\t";
               
        if (StopAtSteady)
        	iofile << "\t" << SteadyDeviation << "\n";

        iofile.close();
    }
}

void incFSI::WriteFileTecplot() {
	
    //amrex::Print()<<"Writing data in tecplot format"<<'\n';
    char *s, *s1, *s2;
    s = new char[50];
    s1 = new char[50];
 
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    int myproc = amrex::ParallelDescriptor::MyProc();
    sprintf(s,"%d",Iter) ; 
    strcpy(s1,"Output/Out") ; strcat(s1,s) ; 
    sprintf(s,"-") ; strcat(s1,s) ;
    sprintf(s,"%d",myproc) ;
    strcat(s1,s) ;
    strcat(s1,".dat") ;
    std::ofstream File(s1, std::ios::out);
    File.flags( std::ios::dec | std::ios::scientific );
    File.precision(5);


    if(!File) {
            std::cout << "Error: Output file for field data couldnot be opened.\n";
            return;
    }
    int fcount = 0;

    amrex::Vector<CFMask> cfmask_(finest_level + 1);
    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        //amrex::Print()<<"ilev = "<<ilev<<'\n';
        cfmask_[ilev].define(ilev, geom, grids, dmap);
    }

    amrex::Real F11_L1_err = 0.0;
    amrex::Real F12_L1_err = 0.0;
    amrex::Real F21_L1_err = 0.0;
    amrex::Real F22_L1_err = 0.0;
    amrex::Real F33_L1_err = 0.0;

    for (int lev = 0; lev <= finest_level ; lev++)
    {
        const amrex::Real *dx = geom[lev].CellSize();
        const amrex::Real *prob_lo = geom[lev].ProbLo();
        const amrex::Real *prob_hi = geom[lev].ProbHi();
        const amrex::Box &domain = geom[lev].Domain();
        amrex::BoxArray xba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(0));
        amrex::BoxArray yba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(1));

        amrex::MultiFab *psi; 
        amrex::MultiFab *kappa; 
        amrex::MultiFab *norm;
        amrex::iMultiFab *solidMask;
        amrex::iMultiFab *PMask;
        amrex::iMultiFab *UMask;
        amrex::iMultiFab *VMask;

        amrex::iMultiFab UMask_NI;//NI stands for No Interface, need a better way of doing this
        amrex::iMultiFab VMask_NI;//NI stands for No Interface, need a better way of doing this
        amrex::iMultiFab PMask_NI;//NI stands for No Interface, need a better way of doing this

        if(N_IF > 0)
        {
            for(int IF = 0;IF<N_IF;IF++)
            {
                auto &&solid = interfaces[lev][IF];
                psi = &solid->Psi();
		kappa = &solid->Kappa_();
		norm = &solid->Normal_(); 
                solidMask = &solid->Mask();
            }
            auto &&mask__= mask_[lev];
            PMask = &mask__->getPMask(); 
            UMask = &mask__->getUMask();
            VMask = &mask__->getVMask();
        }
        else
        {
            UMask_NI.define(xba, dmap[lev], 1, Nghost);
            VMask_NI.define(yba, dmap[lev], 1, Nghost);
            PMask_NI.define(yba, dmap[lev], 1, Nghost);

            UMask_NI.setVal(1);
            VMask_NI.setVal(1);
            PMask_NI.setVal(1);
        }

        if(N_IF > 0)
            File << "TITLE = Flow" << std::endl << "VARIABLES = X, Y, p, Omega, u, v, psi,kappa,norm_x, norm_y,int_mask, pmask, umask , vmask , i, j, lev, RHS, cover, Temperature, Phi,F11,F12,F21,F22,F33,F11_ana,F12_ana,F21_ana,F22_ana,F33_ana,F11_err,F12_err,F21_err,F22_err,F33_err,eps_max,eps_ana,eps_err,ue,ve,divFr,divFz,J_cc,vel_mag,dmg_cc" << std::endl;
	    else if(PhaseField && DamageModel)
	        File << "TITLE = Flow" << std::endl << "VARIABLES = X, Y, p, Omega, u, v, phi, i, j, lev, RHS, F11,F12, F21,F22,F33,eps" << std::endl;
        else
            File << "TITLE = Flow" << std::endl << "VARIABLES = X, Y, p, Omega, u, v, phi, i, j, lev, RHS" << std::endl;
        for(amrex::MFIter mfi(Pressure[lev]); mfi.isValid(); ++mfi)
    	{
            fcount++;
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<amrex::Real const> const &u = xvel[lev].const_array(mfi);
            amrex::Array4<amrex::Real const> const &v = yvel[lev].const_array(mfi);
            amrex::Array4<amrex::Real const> const &P = Pressure[lev].const_array(mfi);
            amrex::Array4<amrex::Real const> const &Vel = U[lev].const_array(mfi);
            amrex::Array4<amrex::Real> const &omega = Vorticity[lev].array(mfi);
            amrex::Array4<amrex::Real const> const &rhs = RHS[lev].const_array(mfi);
            amrex::Array4<int const> const &cfmask = cfmask_[lev].Mask().const_array(mfi);
            amrex::Array4<int> pmask, umask, vmask;
            amrex::Array4<amrex::Real>  Psi, Kappa,Norm; 
            amrex::Array4<int> solid_mask;
            if(N_IF > 0)
            {
                Psi = psi->array(mfi);
	        	Kappa = kappa->array(mfi);
	        	Norm = norm->array(mfi);
                solid_mask = solidMask->array(mfi);
                pmask = PMask->array(mfi);
                umask = UMask->array(mfi);
                vmask = VMask->array(mfi);
            }
            else
            {
                umask = UMask_NI.array(mfi);
                vmask = VMask_NI.array(mfi);
                pmask = PMask_NI.array(mfi);
                solid_mask = PMask_NI.array(mfi);
            }
            amrex::Array4<amrex::Real const> phi, T;
            if(TempField) T = Theta[lev].const_array(mfi);
            if(PhaseField) phi = Phi[lev].const_array(mfi);
    	    amrex::Array4<amrex::Real const> X_ref,Y_ref,F11,F12,F21,F22,F33,eps_max,dmg;
	        if(DamageModel)
            {
                X_ref = Scalars[lev][0].const_array(mfi);
                Y_ref = Scalars[lev][1].const_array(mfi);
                F11 = Scalars[lev][2].const_array(mfi);
                F12 = Scalars[lev][3].const_array(mfi);
                F21 = Scalars[lev][4].const_array(mfi);
                F22 = Scalars[lev][5].const_array(mfi); 
                F33 = Scalars[lev][6].const_array(mfi);
                eps_max = Scalars[lev][7].const_array(mfi);
                dmg = Scalars[lev][8].const_array(mfi);
            }
            amrex::Real bubble_radius = 0.0;
            amrex::Real bubble_radius_0 = 0.0;
            amrex::Real bubble_xcp = 0.0;
            amrex::Real bubble_ycp = 0.0;
            amrex::Real bubble_strain = 0.0;
            amrex::Real bubble_r_dot = 0.0;
            amrex::Real bubble_pressure = 0.0;
            if(N_IF > 0)
            {
                auto &solid = interfaces[finest_level][0];
                bubble_radius_0 = std::pow(solid->Volume_0()/(4.0/3.0*M_PI),(1.0/3.0));
                bubble_radius = std::pow(solid->Volume()/(4.0/3.0*M_PI),(1.0/3.0));
                bubble_xcp = solid->Xcp();
                bubble_ycp = solid->Ycp();
                bubble_strain = solid->AvgIntStrain();
                bubble_r_dot = solid->AvgIntVel();
                bubble_pressure = solid->getP_interface();
            }

            if(std::isnan(bubble_radius))
            {
                amrex::Print()<<"bubble radius = "<<bubble_radius<<'\n';
                exit(1);
            }

            File << "Zone T = Omega I = " << bx.bigEnd(0) - bx.smallEnd(0) + 1 << " J = " << bx.bigEnd(1) - bx.smallEnd(1)+ 1 << std::endl ;

	        amrex::ParallelFor(bx,
       		[&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        	{
	            amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
        	    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                //amrex::Real u_cc = Vel(i, j, k, 0);//0.5*(u(i, j, k) + u(i+1, j, k));
                //amrex::Real v_cc = Vel(i, j, k, 1);//0.5*(v(i, j, k) + v(i, j+1, k)); 
                amrex::Real r = std::hypot(x - bubble_xcp, y - bubble_ycp);
                amrex::Real theta = std::asin((y - bubble_ycp)/r);
                if( x - bubble_xcp < 0 ) theta = M_PI - theta;

                amrex::Real v_r = cos(theta)*(1.0 - bubble_radius/r);
                amrex::Real v_theta = -1.0*sin(theta)*(1.0 - bubble_radius/(2.0*r));

                v_theta = 0.0;
                v_r = -1.0/(3.0*r*r);
                //v_r = bubble_radius*bubble_radius*bubble_r_dot/(r*r);

                amrex::Real u_e = v_r*cos(theta) - v_theta*sin(theta);
                amrex::Real v_e = v_r*sin(theta) + v_theta*cos(theta);
                                                                       
                amrex::Real u_cc = 0.5*(u(i, j, k) + u(i+1, j, k));
                amrex::Real v_cc = 0.5*(v(i, j, k) + v(i, j+1, k)); 

                amrex::Real T_cc = 0.0;
                if(TempField) T_cc = T(i, j, k);
    		    if(std::isnan(T_cc)) T_cc = 100.0;
                amrex::Real phi_cc = 0.0;
                if(PhaseField) phi_cc = phi(i, j, k);
                amrex::Real X_ref_cc = 0.0;
                amrex::Real Y_ref_cc = 0.0;
                amrex::Real F11_cc = 0.0;
                amrex::Real F12_cc = 0.0;
                amrex::Real F21_cc = 0.0;
                amrex::Real F22_cc = 0.0;
                amrex::Real F33_cc = 0.0;
                amrex::Real eps_max_cc = 0.0;
                amrex::Real dmg_cc = 0.0;
                amrex::Real eps_ana = 0.0;
                amrex::Real eps_err = 0.0;
                amrex::Real divFr = 0.0;
                amrex::Real divFz = 0.0;
                amrex::Real F11_cc_ana = 0.0;
                amrex::Real F12_cc_ana = 0.0;
                amrex::Real F21_cc_ana = 0.0;
                amrex::Real F22_cc_ana = 0.0;
                amrex::Real F33_cc_ana = 0.0;
                amrex::Real F11_cc_err = 0.0;
                amrex::Real F12_cc_err = 0.0;
                amrex::Real F21_cc_err = 0.0;
                amrex::Real F22_cc_err = 0.0;
                amrex::Real F33_cc_err = 0.0;
                amrex::Real J_cc = 0.0;
                auto &solid = interfaces[finest_level][0];
                if(DamageModel)
                {
                    X_ref_cc = X_ref(i, j, k);
                    Y_ref_cc = Y_ref(i, j, k);
                    F11_cc = F11(i, j, k);
                    F12_cc = F12(i, j, k);
                    F21_cc = F21(i, j, k);
                    F22_cc = F22(i, j, k);
                    F33_cc = F33(i, j, k);
                    eps_max_cc = eps_max(i, j, k);
                    dmg_cc = dmg(i, j, k);

                    //amrex::Print()<<"bubble radius = "<<bubble_radius<<" , bubble_radius_0 = "<<bubble_radius_0<<'\n';
                    //amrex::Print()<<"bubble center = "<<bubble_xcp<<","<<bubble_ycp<<'\n';
                    amrex::Real r = std::hypot(x-bubble_xcp, y-bubble_ycp);
                    amrex::Real r_0 = std::pow(r*r*r + (std::pow(bubble_radius_0,3.0) - std::pow(bubble_radius,3.0)),1.0/3.0);
                    if(std::isnan(r_0)) r_0 = dx[1];
                    amrex::Real theta = std::acos((x - bubble_xcp)/r);

                    if(pmask(i, j, k) == 1)
                    {
                        amrex::Real ep_1_ana = 0.5*(std::pow(r_0/r,4.0) - 1.0);
                        amrex::Real ep_2_ana = (r*r - r_0*r_0)/(2.0*r_0*r_0);
                        amrex::Real ep_3_ana = (r*r - r_0*r_0)/(2.0*r_0*r_0);
                        eps_ana = std::max(ep_1_ana,ep_2_ana);
                        eps_ana = std::max(eps_ana,ep_3_ana);
                    }
                    if(std::isnan(eps_ana))
                    {   
                        eps_ana = 0.0;            
                    }
                    if(eps_ana != 0.0 && pmask(i, j, k) == 1)
                        eps_err = std::abs((eps_ana - eps_max(i, j, k)))/bubble_strain;

                    divFr = (F22(i, j + 1, k) - F22(i, j - 1, k))/(2.0 * dx[1]) +
                                        (F12(i + 1, j, k) - F12(i - 1, j, k))/(2.0 * dx[0]) +
                                        F22(i, j, k)/y - F33(i, j, k)/y;
                    divFz = (F21(i, j + 1, k) - F21(i, j - 1, k))/(2.0 * dx[1]) +
                                        (F11(i + 1, j, k) - F11(i - 1, j, k))/(2.0 * dx[0]) +
                                        F21(i, j, k)/y;

                    J_cc = (F11_cc*F22_cc - F12_cc*F21_cc);
                    if(isAxisymmetric)
                        J_cc *= F33_cc;
                    if(std::isnan(divFz)) divFz = 0.0;
                                if(std::isnan(divFr)) divFr = 0.0;
                    if(std::isnan(J_cc)) J_cc = 0.0;
                    //Compute Fij from analytical solution
                    amrex::Real R_by_r = r_0/r;
                    //amrex::Print()<<"r_gc = "<<r_gc<<" , theta = "<<theta<<'\n';
                    if(pmask(i, j, k) == 1)
                    {	
                        F11_cc_ana = R_by_r*R_by_r*std::pow(std::cos(theta),2.0) + 1.0/R_by_r * std::pow(std::sin(theta),2.0);
                        F12_cc_ana = (R_by_r*R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);
                        F21_cc_ana = (R_by_r*R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);
                        F22_cc_ana = R_by_r*R_by_r*std::pow(std::sin(theta),2.0) + 1.0/R_by_r * std::pow(std::cos(theta),2.0);
                        F33_cc_ana = 1.0/R_by_r;

                        F11_cc_err = (F11_cc_ana - F11_cc)/F11_cc_ana;
                        F12_cc_err = (F12_cc_ana - F12_cc)/F12_cc_ana;
                        F21_cc_err = (F21_cc_ana - F21_cc)/F21_cc_ana;
                        F22_cc_err = (F22_cc_ana - F22_cc)/F22_cc_ana;
                        F33_cc_err = (F33_cc_ana - F33_cc)/F33_cc_ana;
                                    
                        if(lev == finest_level)
                        {
                            F11_L1_err = std::max(F11_L1_err,std::abs(F11_cc_err));
                            F12_L1_err = std::max(F12_L1_err,std::abs(F12_cc_err));
                            F21_L1_err = std::max(F21_L1_err,std::abs(F21_cc_err));
                            F22_L1_err = std::max(F22_L1_err,std::abs(F22_cc_err));
                            F33_L1_err = std::max(F33_L1_err,std::abs(F33_cc_err));
                        }
                    }
                    if(std::isnan(F11_cc_ana))
                    {
                        amrex::Print()<<"i = "<<i<<" , j = "<<j<<" , pmask  = "<<pmask(i, j, k)<<'\n';
                        amrex::Print()<<"r = "<<r<<" , r_0 = "<<r_0<<" , theta = "<<theta<<"\n";
                        //exit(9);
                    }
                        
                    amrex::Real C11 = F11_cc*F11_cc + F21_cc*F21_cc;
                    amrex::Real C12 = F11_cc*F12_cc + F21_cc*F22_cc;
                    amrex::Real C21 = F12_cc*F11_cc + F22_cc*F21_cc;
                    amrex::Real C22 = F12_cc*F12_cc + F22_cc*F22_cc;
                    amrex::Real C33 = F33_cc*F33_cc;
                    /* 
                    amrex::Real C11 = F11_cc_ana*F11_cc_ana + F21_cc_ana*F21_cc_ana;
                    amrex::Real C12 = F11_cc_ana*F12_cc_ana + F21_cc_ana*F22_cc_ana;
                    amrex::Real C21 = F12_cc_ana*F11_cc_ana + F22_cc_ana*F21_cc_ana;
                    amrex::Real C22 = F12_cc_ana*F12_cc_ana + F22_cc_ana*F22_cc_ana;
                    amrex::Real C33 = F33_cc_ana*F33_cc_ana;
                    */
                    amrex::Real tr_C = C11+C22;
                    amrex::Real det_C = C11*C22 - C21*C12;
                    amrex::Real first_inv = C11+C22+C33;
                    amrex::Real second_inv = C22*C33 + C11*C33 + C11*C22 - C12*C21;
                    amrex::Real third_inv = C33*(C11*C22 - C12*C21);

                    amrex::Real lambda_1;// = 0.5*(tr_C + std::sqrt(tr_C*tr_C - 4.0*det_C));
                    amrex::Real lambda_2;// = 0.5*(tr_C - std::sqrt(tr_C*tr_C - 4.0*det_C));
                    amrex::Real lambda_3;
                    //amrex::Print()<<"first_inv = "<<first_inv<<" , second_inv = "<<second_inv<<" , third_inv = "<<third_inv<<'\n';
                    SolveQuadraticEqn(lambda_1,lambda_2, tr_C, det_C);
                    lambda_3 = C33;
                    //SolveCubicEqn(lambda_1, lambda_2, lambda_3, first_inv, second_inv, third_inv);
                    //if(std::isnan(lambda_1)) lambda_1 = 1.0;
                    //if(std::isnan(lambda_2)) lambda_2 = 1.0;
                    //if(std::isnan(lambda_3)) lambda_3 = 1.0;
                    eps_max_cc = std::max((0.5*(lambda_1 - 1.0)),
                                            (0.5*(lambda_2 - 1.0)));
                    eps_max_cc = std::max(eps_max_cc,(0.5*(lambda_3 - 1.0)));

                    if(eps_ana != 0.0 && pmask(i, j, k) == 1)
                                    eps_err = std::abs((eps_ana - eps_max_cc))/bubble_strain;
                    if(std::isnan(eps_err)) eps_err = 0.0;
                    if(std::isnan(eps_max_cc)) eps_max_cc = 0.0;
                    if(Psi(i, j ,k)>0)eps_max_cc = 0.0;
                    if(std::isnan(eps_ana)) eps_ana = 0.0;
                    if(std::isnan(omega(i, j, k))) omega(i, j, k) = 0.0;
                    /*
                    if(i == 16383 && j == 40)
                        {
                        amrex::Print()<<" i = "<<i<<" , j = "<<j<<'\n';
                        amrex::Print()<<"x = "<<x<<"y = "<<y<<'\n';
                        amrex::Print()<<"bubble radius = "<<bubble_radius<<" , bubble_radius_0 = "<<bubble_radius_0<<'\n';
                                    amrex::Print()<<"r = "<<r<<" , r_0 = "<<r_0<<" , theta = "<<theta<<'\n';
                        amrex::Print()<<"first_inv = "<<first_inv<<" , "<<std::pow(r_0/r,4.0) + 2.0*std::pow(r/r_0,2)<<"\n";
                        amrex::Print()<<"second_inv = "<<second_inv<<" , "<<std::pow(r/r_0,4.0) + 2.0*std::pow(r_0/r,2)<<"\n";
                        amrex::Print()<<"third_inv = "<<third_inv<<"\n";
                        amrex::Print()<<"lambda_1 = "<<lambda_1<<" , lambda_2 = "<<lambda_2<<" , lambda_3 = "<<lambda_3<<"\n";
                        amrex::Print()<<"ep-1 = "<<0.5*(lambda_1 - 1.0)<<" , ep_2 = "<<0.5*(lambda_2 - 1.0)<<" , ep_3 = "<<0.5*(lambda_3 - 1.0)<<"\n";
                        amrex::Print()<<"lambda_1 = "<<std::pow(r_0/r,4.0)<<" , lambda_2 = "<<std::pow(r/r_0,2.0)<<" , lambda_3 = "<<std::pow(r/r_0,2.0)<<"\n";
                        amrex::Print()<<" eps_ana = "<<eps_ana<<", eps_max_cc = "<<eps_max_cc<<" , eps_err = "<<eps_err<<'\n';
                    }
                    */
                    if(std::isnan(F11_cc))
                    {
                        amrex::Print()<<"lev = "<<lev<<" , i = "<<i<<" , j = "<<j<<'\n';
                        amrex::Print()<<"F11_cc = "<<F11_cc<<'\n';
                        //exit(2);
                    }
                    if(std::isnan(F12_cc))
                    {
                        amrex::Print()<<"lev = "<<lev<<" , i = "<<i<<" , j = "<<j<<'\n';
                        amrex::Print()<<"F12_cc = "<<F12_cc<<'\n';
                        // exit(3);
                    }
                    if(std::isnan(F21_cc))
                    {
                        amrex::Print()<<"lev = "<<lev<<" , i = "<<i<<" , j = "<<j<<'\n';
                        amrex::Print()<<"F21_cc = "<<F21_cc<<'\n';
                        //exit(4);
                    }
                    if(std::isnan(F22_cc))
                    {
                        amrex::Print()<<"lev = "<<lev<<" , i = "<<i<<" , j = "<<j<<'\n';
                        amrex::Print()<<"F22_cc = "<<F22_cc<<'\n';
                        //exit(5);
                    }
                    if(std::isnan(F33_cc))
                    {
                        amrex::Print()<<"lev = "<<lev<<" , i = "<<i<<" , j = "<<j<<'\n';
                        amrex::Print()<<"F33_cc = "<<F33_cc<<'\n';
                        //exit(6);
                    }
                }

                amrex::Real cover = 0.0;
                if(cfmask(i,j,k) == CFMask::covered)
                    cover = 1.0;

                if(lev == finest_level && i == 32868 && j == 0 )
                {
                    //std::cout.precision(17);
                    std::cout<<"Cell : "<<i <<" , "<<j<<'\n';
                    //std::cout<<"P dn = "<<std::fixed<<P(i,j,k)<<'\n';
                    //std::cout<<"P up = "<<std::fixed<<P(i,j+1,k)<<'\n';
                    //std::cout<<"u dn = "<<std::fixed<<u(i,j,k)<<'\n';
                    //std::cout<<"u up = "<<std::fixed<<u(i,j+1,k)<<'\n';
                    //std::cout<<"y dn = "<<std::fixed<<v(i,j,k)<<'\n';
                    //std::cout<<"y 0  = "<<std::fixed<<v(i,j+1,k)<<'\n';
                    //std::cout<<"y up = "<<std::fixed<<v(i,j+2,k)<<'\n';
                    amrex::Print()<<"pmask = "<<pmask(i,j,k)<<'\n';
                }


                if(N_IF > 0)
                {
                    double P_cc = P(i,j,k);
                    double omega_cc = omega(i,j,k); 
                    if(Psi(i, j ,k)>0)
                    {

                        P_cc = bubble_pressure;
                        u_cc = 0.0;
                        v_cc = 0.0;
                        omega_cc = 0.0;
                        dmg_cc = 0.0;
			phi_cc = -1.0;
                    } 
                    if(phi_cc < -.98)
                    {
                        dmg_cc = 0.0;
                    }
                    double vel_mag = std::sqrt(u_cc*u_cc + v_cc*v_cc);
                    File << x <<"\t"<< y <<"\t"
                            << P_cc <<"\t"<< omega_cc <<"\t"
                            << u_cc <<"\t"<< v_cc  <<"\t"<<Psi(i, j ,k)<<"\t"<<Kappa(i, j ,k)<<"\t"
                            << Norm(i, j ,k,0)<<"\t"<<Norm(i, j ,k,1)<<"\t"
                            << solid_mask(i, j, k)<<"\t"<<pmask(i, j, k)<<"\t"
                            << umask(i, j, k)<<"\t"<<vmask(i, j, k)<<"\t"
                            << i <<"\t"<< j <<"\t"
                            << lev <<"\t"<<rhs(i,j,k)<<'\t'<<cover<<'\t'
                            << T_cc <<"\t"<< phi_cc <<"\t"
                            <<F11_cc<<"\t"<<F12_cc<<"\t"
                            <<F21_cc<<"\t"<<F22_cc<<"\t" 
                            <<F33_cc<<"\t"
                            <<F11_cc_ana<<"\t"<<F12_cc_ana<<"\t"
                            <<F21_cc_ana<<"\t"<<F22_cc_ana<<"\t"
                            <<F33_cc_ana<<"\t"
                            <<F11_cc_err<<"\t"<<F12_cc_err<<"\t"
                            <<F21_cc_err<<"\t"<<F22_cc_err<<"\t"
                            <<F33_cc_err<<"\t"
                            <<eps_max_cc<<"\t"<<eps_ana<<'\t'<<eps_err<<'\t'
                            <<u_e<<'\t'<<v_e<<'\t'
                            <<divFr<<'\t'<<divFz<<'\t'<<J_cc<<"\t"<<vel_mag<<'\t'<<dmg_cc<<'\n';
                }
                else
                {
                    if(PhaseField && DamageModel)
                            File << x <<"\t"<< y <<"\t"
                            << P(i,j,k) <<"\t"<< omega(i,j,k) <<"\t"
                            << u_cc <<"\t"<< v_cc <<"\t"<< phi(i,j,k) <<"\t"
                            << i <<"\t"<< j <<"\t"
                            << lev <<"\t"<<rhs(i,j,k)<<'\t'<<F11_cc<<"\t"<<F12_cc<<"\t"
                            <<F21_cc<<"\t"<<F22_cc<<"\t"
                            <<F33_cc<<"\t"<<eps_max_cc<<'\n';
                    else if(PhaseField)
                            File << x <<"\t"<< y <<"\t"
                            << P(i,j,k) <<"\t"<< omega(i,j,k) <<"\t"
                            << u_cc <<"\t"<< v_cc <<"\t"<< phi(i,j,k) <<"\t"
                            << i <<"\t"<< j <<"\t"
                            << lev <<"\t"<<rhs(i,j,k)<<'\t'<<'\n';


                    else
                            File << x <<"\t"<< y <<"\t"
                            << P(i,j,k) <<"\t"<< omega(i,j,k) <<"\t"
                            << u_cc <<"\t"<< v_cc <<"\t"<< "0.0" <<"\t"
                            << i <<"\t"<< j <<"\t"
                            << lev <<"\t"<<rhs(i,j,k)<<'\t'<<'\n';
                }
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(F11_L1_err);
        amrex::ParallelDescriptor::ReduceRealMax(F21_L1_err);
        amrex::ParallelDescriptor::ReduceRealMax(F12_L1_err);
        amrex::ParallelDescriptor::ReduceRealMax(F22_L1_err);
        amrex::ParallelDescriptor::ReduceRealMax(F33_L1_err);

        if(lev == finest_level)
        {
            amrex::Print()<<"F11_L1_err = "<<F11_L1_err<<'\n';
            char *s, *s1, *s2;
            s = new char[50];
            s1 = new char[50];
            s2 = new char[50];
            strcpy(s1, "Output/FL1_err");
            strcpy(s2, ".dat");
            strcat(s1, s2);
            // std::string filename(s1);
            std::ofstream ofs_(s1, std::ofstream::out | std::ofstream::app);
            ofs_.flags(std::ios::dec | std::ios::scientific);
            ofs_.precision(7);
            //for (auto &&solid : interfaces)
            //{
            amrex::Print(ofs_) <<Time<<"\t"<< F11_L1_err<<"\t"<<F12_L1_err
                                  <<"\t"<< F21_L1_err<<"\t"<<F22_L1_err
                                  <<"\t"<< F33_L1_err
              <<'\n';
            //} 
            ofs_.close();
         }
    }
    if (fcount == 0)
    {
                amrex::Print()<<"Writing empty file"<<'\n';
                File << "Zone T = Omega I = " << " 1 " << " J = " << " 1 " << std::endl ;
                if(N_IF > 0)
                    File << "0.0" <<"\t"<< "0.0" <<"\t"<<"0.0" <<"\t"<<"0.0" <<"\t"
                         <<"\t"<<"0.0" <<"\t"
                         << "0.0" <<"\t"<< "0.0" <<"\t"<<"0.0" <<"\t"
			             << "0.0" <<"\t"<< "0.0" <<"\t"<<"0.0" <<"\t"
                         << "0.0" <<"\t"<< "0.0" <<"\t"<<"0.0" <<"\t"<< "0.0" <<"\t"<<"0.0" <<"\t"
                         << "0.0" <<"\t"<< "0.0" <<"\t"<< "0.0" <<"\t"<< "0.0" <<"\t"<< "0.0" <<"\t"
                         << "0" <<"\t"<< "0" << '\t' << "0"<<"\t"<<"0.0"<<'\t'<<"0.0"<<'\t'
			             << "0"<<"\t"<<"0.0"<<'\t'<<"0.0"<<'\t'<<"0.0"<<'\t'<<"0.0"<<'\t'<<"0.0"<<'\t'
			             << "0"<<"\t"<<"0.0"<<'\t'<<"0.0"<<'\t'<<"0.0"<<'\t'<<"0.0"<<'\t'<<"0.0"<<'\t'
                         <<"0.0"<<'\t'<<"0.0"<<'\t'<<"0.0"<<'\t'<<"0.0"<<'\t'<<"0.0"<<'\t'
                         <<"0.0"<<'\t'<<"0.0"<<'\t'<<"0.0"<<'\t'<<"0.0"<<'\t'<<"0.0"<<'\n';
                else
                    File << "0.0" <<"\t"<< "0.0" <<"\t"<<"0.0" <<"\t"
                         << "0.0" <<"\t"<< "0.0" <<"\t"<<"0.0" <<"\t"<<"0.0"<<"\t"
                         << "0" <<"\t"<< "0" <<'\t' <<"0"<<"\t"<<"0.0"<<'\t'
                         <<'\n';

                
     }
     File.close();
     //if(N_IF > 0)WriteInterface();
}

void incFSI::WriteFileTecplotW_Ghost() {
	
    //amrex::Print()<<"Writing data in tecplot format; N_IF = "<< N_IF <<'\n';
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    int myproc = amrex::ParallelDescriptor::MyProc();

    bool exit  = false;

    for (int lev = 0; lev <= finest_level ; lev++)
    {
        char *s, *s1, *s2;
        s = new char[50];
        s1 = new char[50];
        s2 = new char[50];
        int fcount = 0;
        const amrex::Real *dx = geom[lev].CellSize();
        const amrex::Real *prob_lo = geom[lev].ProbLo();
        const amrex::Real *prob_hi = geom[lev].ProbHi();
        const amrex::Box &domain = geom[lev].Domain();
        amrex::MultiFab *psi; 
        amrex::iMultiFab *solidMask;
        amrex::iMultiFab *PMask;
        amrex::iMultiFab *UMask;
        amrex::iMultiFab *VMask;
        if(N_IF > 0)
        {
            for(int IF = 0;IF<N_IF;IF++)
            {
                auto &&solid = interfaces[lev][0];
                psi = &solid->Psi();
                solidMask = &solid->Mask();
            }
            auto &&mask__= mask_[lev];
            PMask = &mask__->getPMask(); 
            UMask = &mask__->getUMask();
            VMask = &mask__->getVMask();
        }        

        sprintf(s,"%d",Iter) ; 
        strcpy(s1,"Output/Out_w_patchu_Iter") ; strcat(s1,s) ; 
        sprintf(s,"-proc") ; strcat(s1,s) ;
        sprintf(s,"%d",myproc) ; strcat(s1,s) ;
        strcpy(s2,s1);

        for(amrex::MFIter mfi(Pressure[lev]); mfi.isValid(); ++mfi)
    	{
            fcount++;
            strcpy(s2,s1);
            sprintf(s,"-") ; strcat(s2,s) ;
	        sprintf(s,"%d",lev) ; strcat(s2,s) ;
            sprintf(s,"-") ; strcat(s2,s) ;
	        sprintf(s,"%d",fcount) ; strcat(s2,s) ;
            strcat(s2,".dat") ;
            std::ofstream File(s2, std::ios::out);
            if(!File) 
            {
                std::cout << "Error: Output file for field data couldnot be opened.\n";
                return;
            }
	        File.flags( std::ios::dec | std::ios::scientific );
	        File.precision(5);

            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<amrex::Real > const &u = xvel[lev].array(mfi);
            amrex::Array4<amrex::Real > const &v = yvel[lev].array(mfi);
            amrex::Array4<amrex::Real > const &P = Pressure[lev].array(mfi);
            amrex::Array4<amrex::Real> const &Vel = U[lev].array(mfi);
            amrex::Array4<int> pmask, umask, vmask;
            amrex::Array4<amrex::Real>  Psi; 
            amrex::Array4<int> solid_mask;
            amrex::Array4<amrex::Real const> phi, T;
            if(TempField) T = Theta[lev].const_array(mfi);
            if(PhaseField) phi = Phi[lev].const_array(mfi);

            if(N_IF > 0)
            {
                Psi = psi->array(mfi);
                solid_mask = solidMask->array(mfi);
                pmask = PMask->array(mfi);
                umask = UMask->array(mfi);
                vmask = VMask->array(mfi);
            }            
            amrex::Array4<amrex::Real const> X_ref,Y_ref,F11,F12,F21,F22;
            if(DamageModel)
            {
                 X_ref = Scalars[lev][0].const_array(mfi);
                 Y_ref = Scalars[lev][1].const_array(mfi);
                 F11 = Scalars[lev][2].const_array(mfi);
                 F12 = Scalars[lev][3].const_array(mfi);
                 F21 = Scalars[lev][4].const_array(mfi);
                 F22 = Scalars[lev][5].const_array(mfi);
            }

            if(N_IF > 0 )File << "TITLE = Flow" << std::endl << "VARIABLES = X, Y, p, u, v, i, j, lev, ghost,int_mask,pmask,umask,vmask, Psi, Temperature, Phi,F11,F12,F21,F22" << std::endl;
            else File << "TITLE = Flow" << std::endl << "VARIABLES = X, Y, p, u, v, phi, i, j, lev, ghost" << std::endl;
            File << "Zone T = Omega I = " << bx.bigEnd(0) + Nghost - (bx.smallEnd(0) - Nghost) + 1 << " J = " << bx.bigEnd(1) + Nghost - (bx.smallEnd(1) - Nghost) + 1 << std::endl ;

			for(int j = bx.smallEnd(1) - Nghost; j <= bx.bigEnd(1) + Nghost ; j++)
			{
			    for(int i = bx.smallEnd(0) - Nghost; i <= bx.bigEnd(0) + Nghost ; i++)
			    {
				int k = 0;
				amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
        	        	amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                    		amrex::Real u_cc = Vel(i ,j ,k ,0);//0.5*(u(i, j, k) + u(i+1, j, k));
                    		amrex::Real v_cc = Vel(i ,j ,k ,1);//0.5*(v(i, j, k) + v(i, j+1, k));
				//if(i == 16463 && j == 0)
			        //{
				//    amrex::Print(-1)<<" i = "<<i<<" , j = "<<j<<'\n';
				//    amrex::Print(-1)<<"u_cc = "<<u_cc<<" , v_cc = "<<v_cc<<'\n';
				//}
                                //amrex::Real u_cc = 0.5*(u(i, j, k) + u(i+1, j, k));
                                //amrex::Real v_cc = 0.5*(v(i, j, k) + v(i, j+1, k));
                                amrex::Real T_cc = 0.0;
                                if(TempField) T_cc = T(i, j, k);
                                amrex::Real phi_cc = 0.0;
                                if(PhaseField) phi_cc = phi(i, j, k);
                                amrex::Real X_ref_cc = 0.0;
                                amrex::Real Y_ref_cc = 0.0;
                                amrex::Real F11_cc = 0.0;
                                amrex::Real F12_cc = 0.0;
                                amrex::Real F21_cc = 0.0;
                                amrex::Real F22_cc = 0.0;
                                if(DamageModel)
                                {
                                    X_ref_cc = X_ref(i, j, k);
                                    Y_ref_cc = Y_ref(i, j, k);
                                    F11_cc = F11(i, j, k);
                                    F12_cc = F12(i, j, k);
                                    F21_cc = F21(i, j, k);
                                    F22_cc = F22(i, j, k);
                               }




                    		int isghost = 1; 
                    		if(i>= bx.smallEnd(0) && i <= bx.bigEnd(0) && j>= bx.smallEnd(1) && j <= bx.bigEnd(1) )
                       		isghost = 0;
				
                                if(std::isnan(P(i,j,k)))
                                {
                                   //amrex::Print()<<"P nan at "<<i<<" , "<<j<<" lev = "<<lev<<'\n';
                                   exit = true;
                                   P(i,j,k) = 999.0;
                                }
                                if(std::isnan(u_cc))
                                {
                                   //amrex::Print()<<"u_cc nan at "<<i<<" , "<<j<<" lev = "<<lev<<'\n';
                                   exit = true;
                                   u_cc = 999.0;
                                }
                                if(std::isnan(v_cc))
                                {
                                   //amrex::Print()<<"v_cc nan at "<<i<<" , "<<j<<" lev = "<<lev<<'\n';
                                   exit = true;
                                   v_cc = 999.0;
                                }
                                if(std::isnan(u(i, j, k)))
                                {
                                   //amrex::Print()<<"u nan at "<<i<<" , "<<j<<" lev = "<<lev<<'\n';
                                   exit = true;
                                   //u_cc = 999.0;
                                }
                                if(std::isnan(v(i, j, k)))
                                {
                                   //amrex::Print()<<"v nan at "<<i<<" , "<<j<<" lev = "<<lev<<'\n';
                                   exit = true;
                                   //u_cc = 999.0;
                                }
                                if(std::isnan(F11_cc))
                                {
                                   //amrex::Print()<<"F11_cc nan at "<<i<<" , "<<j<<" lev = "<<lev<<'\n';
                                   exit = true;
                                   F11_cc = 999.0;
                                }
                                if(std::isnan(F12_cc))
                                {
                                   //amrex::Print()<<"F12_cc nan at "<<i<<" , "<<j<<" lev = "<<lev<<'\n';
                                   exit = true;
                                   F12_cc = 999.0;
                                }
                                if(std::isnan(F21_cc))
                                {
                                   //amrex::Print()<<"F21_cc nan at "<<i<<" , "<<j<<" lev = "<<lev<<'\n';
                                   exit = true;
                                   F21_cc = 999.0;
                                }
                                if(std::isnan(F22_cc))
                                {
                                   //amrex::Print()<<"F22_cc nan at "<<i<<" , "<<j<<" lev = "<<lev<<'\n';
                                   exit = true;
                                   F22_cc = 999.0;
                                }
				amrex::Real Psi_cc = Psi(i,j,k);
				if(std::isnan(Psi(i,j,k)))
			            Psi_cc = 999.0;
					




                               
                                if(N_IF > 0)
                                {
                    	 	    File << x <<"\t"<< y <<"\t"
                             	    << P(i,j,k) <<"\t"
                             	    << u_cc <<"\t"<< v_cc <<"\t"
                             	    << i <<"\t"<< j <<"\t"
                             	    << lev <<"\t"<<isghost<<"\t"
                                    << solid_mask(i, j, k)<<"\t"<<pmask(i, j, k)<<"\t"
                                    << umask(i, j, k)<<"\t"<<vmask(i, j, k)<<"\t"<<Psi_cc<<"\t"
                                    << T_cc <<"\t"<< phi_cc <<"\t"
                                    <<F11_cc<<"\t"<<F12_cc<<"\t"
                                    <<F21_cc<<"\t"<<F22_cc<<"\t"
				    <<'\n';
                                }
                                else
                                {
                                    if(PhaseField)
                                        File << x <<"\t"<< y <<"\t"
                                        << P(i,j,k) <<"\t"
                                        << u_cc <<"\t"<< v_cc <<"\t"<< phi(i,j,k) <<'\t'
                                        << i <<"\t"<< j <<"\t"
                                        << lev <<"\t"<<isghost<<"\t"
                                        <<'\n';

                                    else
                                        File << x <<"\t"<< y <<"\t"
                                        << P(i,j,k) <<"\t"
                                        << u_cc <<"\t"<< v_cc <<"\t"<<0.0<<'\t'
                                        << i <<"\t"<< j <<"\t"
                                        << lev <<"\t"<<isghost<<"\t"
                                        <<'\n';
                                }
			    }
			}
			
			File.close();
        }
    }
    if(exit) 
    {
        amrex::Print()<<"Found nan"<<'\n';
        //std::exit(9);
    } 
}

void incFSI::WriteInterface()
{
    //int r_end[2] = {0,0};
    //int l_end[2] = {100000000,100000000};
    //int r_end_axis[2] = {0,0};
    //int l_end_axis[2] = {100000000,100000000};
    //Jet volume calculation implemented assuming single bubble
    amrex::Real r_end = -1e9;
    amrex::Real l_end = 1e9;
    amrex::Real r_end_axis = -1e9;
    amrex::Real l_end_axis = 1e9;
    amrex::Real ur_end = 0.0;
    amrex::Real ul_end = 0.0;

    //for(int lev = 0; lev <= max_level;lev++)
    {
	int lev = finest_level;
        int nsolid = 0;
        const amrex::Real *prob_lo = geom[lev].ProbLo();
        const amrex::Real *dx = geom[lev].CellSize();
        int myproc = amrex::ParallelDescriptor::MyProc();

        for (auto &&solid : interfaces[lev])
        {
            char *s, *s1, *s2;
            s = new char[50];
            s1 = new char[50];
            strcpy(s1, "Output/Interface-");
            strcat(s1, "lev"); sprintf(s, "%d", lev); strcat(s1, s); strcat(s1, "-");
            strcat(s1, "obj"); sprintf(s, "%d", nsolid); strcat(s1, s); strcat(s1, "-");
            strcat(s1, "Iter"); sprintf(s, "%d", Iter); strcat(s1, s); strcat(s1, "-");
            strcat(s1, "proc"); sprintf(s, "%d", myproc); strcat(s1, s); strcat(s1, ".dat");
            std::ofstream ofs(s1, std::ofstream::out | std::ofstream::app);
            ofs.flags(std::ios::dec | std::ios::scientific);
            ofs.precision(5);
            //ofs  << "VARIABLES = X, Y, p_I,T_I ,Norm_shear,type,frac,Interpolation_weights,i_intercepts,stencil_" << std::endl;
            ofs  << "VARIABLES = X, Y,Xc,Yc,type,frac,weights, i_intercepts_i,P_int,u_int,v_int, ep_max, dmg, phi, reg_stat, kappa" << std::endl;

            char *s3, *s4, *s5;
            s3 = new char[50];
            s4 = new char[50];
            strcpy(s4, "Output/Interface_nbr-");
            strcat(s4, "lev"); sprintf(s3, "%d", lev); strcat(s4, s3); strcat(s4, "-");
            strcat(s4, "obj"); sprintf(s3, "%d", nsolid); strcat(s4, s3); strcat(s4, "-");
            strcat(s4, "Iter"); sprintf(s, "%d", Iter); strcat(s4, s3); strcat(s4, "-");
            strcat(s4, "proc"); sprintf(s3, "%d", myproc); strcat(s4, s3); strcat(s4, ".dat");
            std::ofstream ofs1(s4, std::ofstream::out | std::ofstream::app);
            ofs1.flags(std::ios::dec | std::ios::scientific);
            ofs1.precision(5);
            //ofs  << "VARIABLES = X, Y, p_I,T_I ,Norm_shear,type,frac,Interpolation_weights,i_intercepts,stencil_" << std::endl;
            ofs1  << "VARIABLES = X, Y" << std::endl;



            int count_pt = 0;
            int avg_int_ep = 0;

            amrex::MultiFab &psi = solid->Psi();
            for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();
            	auto &icpt_data = solid->getInterceptData()[mfi];
                for (auto &&idt : icpt_data)
                {
                	if (bx.contains(idt.cellid_))
                	{
                            count_pt++;
                            for (int i = 0; i < idt.n_intercepts; i++)
                            {

                                amrex::Real xloc = prob_lo[0] + dx[0] * (idt.cellid_[0] + 0.5);
                                amrex::Real yloc = prob_lo[1] + dx[1] * (idt.cellid_[1] + 0.5);
                                amrex::Real xc = prob_lo[0] + dx[0] * (idt.cellid_[0] + 0.5);
                                amrex::Real yc = prob_lo[1] + dx[1] * (idt.cellid_[1] + 0.5);
                                if (idt.type_[i] == 0)
                   	            {   
                                   xloc += -idt.frac_[i] * dx[0];
                    	        }
                   	            else if (idt.type_[i] == 1)
                    	        {   
                                   yloc += -idt.frac_[i] * dx[1];
                    	        }
			                    amrex::Real ep_max = 0.0;

                               /*
                               if(idt.cellid_[0] <= l_end[0])
			       {
			           l_end[0] = idt.cellid_[0];
				   l_end[1] = idt.cellid_[1];
			       }
                               if(idt.cellid_[0] >= r_end[0])
                               {
                                   r_end[0] = idt.cellid_[0];
                                   r_end[1] = idt.cellid_[1];
                               }
			       if(idt.cellid_[1] == 0)
                               {
                                   if(idt.cellid_[0] <= l_end_axis[0])
                                   {
                                       l_end_axis[0] = idt.cellid_[0];
                                       l_end_axis[1] = idt.cellid_[1];
                                   }
                                   if(idt.cellid_[0] >= r_end_axis[0])
                                   {
                                       r_end_axis[0] = idt.cellid_[0];
                                       r_end_axis[1] = idt.cellid_[1];
                                   }
			       }
			       */
                               if(xloc <= l_end)
                               {
                                   l_end = xloc;
                               }
                               if(xloc >= r_end)
                               {
                                   r_end = xloc;
                               }
                               if(yloc == prob_lo[1] + 0.5*dx[1])
                               {
                                   if(xloc <= l_end_axis)
                                   {
                                       l_end_axis = xloc;
                                   }
                                   if(xloc >= r_end_axis)
                                   {
                                       r_end_axis = xloc;
                                   }
                               }


			       if(DamageModel)
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
                                    if(std::isnan(lambda_1)) lambda_1 = 1.0;
				    if(std::isnan(lambda_2)) lambda_2 = 1.0;
				    if(std::isnan(lambda_3)) lambda_3 = 1.0;
                                    ep_max = std::max(std::abs(0.5*(lambda_1 - 1.0)),
                                                      std::abs(0.5*(lambda_2 - 1.0)));
                                    ep_max = std::max(ep_max, std::abs(0.5*(lambda_3 - 1.0)));

				    avg_int_ep = avg_int_ep + ep_max;
			       }

				ofs << xloc << "\t"
                    		    << yloc << "\t"
				    << xc << "\t"
                                    << yc << "\t"
                                    << idt.type_[i]<<"\t"
                                    << idt.frac_[i]<<'\t'
			            << idt.Interpolation_weights[i] <<'\t'
                                    << i <<'\t'
                                    << idt.P[i]<<'\t'<< idt.u[i]<<'\t'<< idt.v[i]<<'\t'<<ep_max<<'\t'
                                    << idt.Dmg_Int[i]<<'\t'<<idt.phi_Int[i]<<'\t'
				    << idt.r_flag<<'\t'<< idt.kappa_[i]<<"\n";
				    //<< idt.r_flag<<'\t'<< idt.s_flag<<'\t'<<idt.kappa_[i]<<"\n";
				 /*
                                 if( (idt.cellid_[0] == 127 && idt.cellid_[1] == 18)   )
                                 {
                                     std::cout<<"i =  "<<idt.cellid_[0]<<" , j = "<<idt.cellid_[1]<<'\n';
                                     std::cout<<" intercept = "<<i<<'\n';
                                     std::cout<<"xloc =  "<<xloc<<" , yloc = "<<yloc<<'\n';
                                     //std::cout<<" processor = "<< amrex::ParallelDescriptor::MyProc()<<'\n';
				     std::cout<<"FII_int = "<<idt.F11_int[i]<<'\n';
                                 }
				 */
                            }
			    for(int i = 0;i < idt.n_nbr ;i++)
		            {

                                amrex::Real xc = prob_lo[0] + dx[0] * (idt.cellid_[0] + 0.5);
                                amrex::Real yc = prob_lo[1] + dx[1] * (idt.cellid_[1] + 0.5);
				//amrex::Print()<<"idt.n_nbr = "<<idt.n_nbr<<"\n";
				//amrex::Print()<<"x,y = "<<xc<<" , "<<yc<<"\n";
                                amrex::Real xloc = prob_lo[0] + dx[0] * (idt.nbr_cellid_[i][0] + 0.5);
                                amrex::Real yloc = prob_lo[1] + dx[1] * (idt.nbr_cellid_[i][1] + 0.5);
                                ofs1 << xloc << "\t"
                                    << yloc << "\t"
                                    <<"\n";
			    }    
		        }
                }
            }
	    if(count_pt == 0) ofs<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<'\t'<<0.0<<"\t"<<0.0<<'\t'<<0.0<<"\t"<<0.0<<'\t'<<0.0<<'\t'<<0.0<<'\t'<<"\n";
            //if(count_pt == 0) ofs<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<'\t'<<0.0<<"\t"<<0.0<<'\t'<<0.0<<"\t"<<0.0<<'\t'<<0.0<<'\t'<<0.0<<'\t'<<0.0<<"\n";
            
            ofs.close();
	    ofs1.close();
            
            //amrex::ParallelDescriptor::ReduceIntMax(r_end[0]);
	    //amrex::ParallelDescriptor::ReduceIntMin(l_end[0]);
            //amrex::ParallelDescriptor::ReduceIntMax(r_end_axis[0]);
            //amrex::ParallelDescriptor::ReduceIntMin(l_end_axis[0]);
            amrex::ParallelDescriptor::ReduceRealMax(r_end);
            amrex::ParallelDescriptor::ReduceRealMin(l_end);
            amrex::ParallelDescriptor::ReduceRealMax(r_end_axis);
            amrex::ParallelDescriptor::ReduceRealMin(l_end_axis);

	    //amrex::Print(-1)<<"r end = "<<r_end<<" , l end = "<<l_end<<'\n';
            //amrex::Print(-1)<<"r end axis = "<<r_end_axis<<" , l end axis = "<<l_end_axis<<'\n';

            for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
		auto &icpt_data = solid->getInterceptData()[mfi];
		for (auto &&idt : icpt_data)
                {
                    if (bx.contains(idt.cellid_))
                    {
	                for (int i = 0; i < idt.n_intercepts; i++)
		        {
		               amrex::Real xloc = prob_lo[0] + dx[0] * (idt.cellid_[0] + 0.5);
                               amrex::Real yloc = prob_lo[1] + dx[1] * (idt.cellid_[1] + 0.5);
                               if (idt.type_[i] == 0)
                               {
                                   xloc += -idt.frac_[i] * dx[0];
                               }
                               else if (idt.type_[i] == 1)
                               {
                                   yloc += -idt.frac_[i] * dx[1];
                               }
			       if(std::abs(xloc - r_end) <= 0.01*dx[0])
			           ur_end = yloc;
			       if(std::abs(xloc - l_end) <= 0.01*dx[0])
			           ul_end = yloc;
			}
                    }
                }	
            }
            amrex::ParallelDescriptor::ReduceRealMax(ur_end);
            amrex::ParallelDescriptor::ReduceRealMax(ul_end);
	    //amrex::Print(-1)<<"ur_end = "<<ur_end<<" , ul_end = "<<ul_end<<"\n";
            nsolid++;
        }
    }
    //amrex::Print()<<"r end = "<<r_end<<" , l end = "<<l_end<<'\n';
    //amrex::Print()<<"r end axis = "<<r_end_axis<<" , l end axis = "<<l_end_axis<<'\n';
    //amrex::Print()<<"ur_end = "<<ur_end<<" , ul_end = "<<ul_end<<"\n";
    //Jet detection criteria
    int jet = 0;
    if(r_end > r_end_axis || l_end < l_end_axis)
        jet = 1;

    amrex::ParallelDescriptor::ReduceIntMax(jet);
    //amrex::Print()<<"jet = "<<jet<<'\n';
    //Compute jet volume and average jet velocity
    amrex::Vector<CFMask> cfmask_(finest_level + 1);
    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        cfmask_[ilev].define(ilev, geom, grids, dmap);
    }
    amrex::Real jet_vol = 0.0;
    amrex::Real jet_x_mom = 0.0;
    amrex::Real jet_y_mom = 0.0;
    amrex::Real jet_x_vel = 0.0;
    amrex::Real jet_y_vel = 0.0;
    if(jet == 1)
    {
        for (int i = 0; i < N_IF; i++)
        {
            for(int ilev = 0; ilev <= finest_level; ilev++)
            {
                amrex::Real jet_vol_part = 0.0;
                const amrex::Real *dx = geom[ilev].CellSize(); 
                const amrex::Real *prob_lo = geom[ilev].ProbLo();
                auto &&solid = interfaces[ilev][i];

                amrex::MultiFab &psi = solid->Psi();
                for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
                {
                    const amrex::Box &bx = mfi.validbox();
                    amrex::Array4<amrex::Real const> const &f = psi.const_array(mfi);
                    amrex::Array4<int const> const &cfmask = cfmask_[ilev].Mask().const_array(mfi);
                    amrex::Array4<amrex::Real> const &Vel = U[ilev].array(mfi);
                    amrex::ParallelFor(bx,
                    [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                        amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
		        amrex::Real EPS_LS = 0.75 * std::min(dx[0], dx[1]);
			amrex::Real Phi = 0.5 + 0.5 * tanh( -1.0*f(i,j,k) / (2.0 * EPS_LS));
                        if(!cfmask(i,j,k) == CFMask::covered)
                        {
			    if((y < ur_end || y < ul_end) && (x > l_end && x < r_end))
			    {
				jet_vol += 2.0 * M_PI * dx[0] * dx[1] * fabs(y) * Phi;
				jet_x_mom += 2.0 * M_PI * dx[0] * dx[1] * fabs(y) * Phi * Vel(i,j,k,0);
				jet_y_mom += 2.0 * M_PI * dx[0] * dx[1] * fabs(y) * Phi * Vel(i,j,k,1);
			    }
                        }
                     });
                }
            }
	}
	amrex::ParallelDescriptor::ReduceRealSum(jet_vol);
	amrex::ParallelDescriptor::ReduceRealSum(jet_x_mom);
	amrex::ParallelDescriptor::ReduceRealSum(jet_y_mom);
	jet_x_vel = jet_x_mom/jet_vol;
	jet_y_vel = jet_y_mom/jet_vol;
	//amrex::Print(-1)<<"jet_vol = "<<jet_vol<<'\n';
	//amrex::Print(-1)<<"jet_x_vel = "<<jet_x_vel<<'\n';
	//amrex::Print(-1)<<"jet_y_vel = "<<jet_y_vel<<'\n';
    }


    //For now, hard coded for the first bubble
    //amrex::Print()<<"N_IF = "<<N_IF<<'\n';
    for(int i = 0; i < N_IF; i++)
    {
        char *s, *s1, *s2;
        s = new char[50];
        s1 = new char[50];
        s2 = new char[50];
        sprintf(s, "%d", i);
        strcpy(s1, "Output/Radius");
        strcpy(s2, ".dat");
        strcat(s1, s);
        strcat(s1, s2);
        // std::string filename(s1);
        std::ofstream ofs_(s1, std::ofstream::out | std::ofstream::app);
        ofs_.flags(std::ios::dec | std::ios::scientific);
        ofs_.precision(7);
        //for (auto &&solid : interfaces)
        //{
        auto &solid = interfaces[finest_level][i];
	    amrex::Real bubble_radius = std::pow(solid->Volume()/(4.0/3.0*M_PI),(1.0/3.0));
            amrex::Print(ofs_) <<Time<<"\t"<< bubble_radius<<"\t"<<solid->AvgIntVel()
                 <<"\t"<<solid->AvgIntPres()<<"\t"<<solid->Volume()<<"\t"<<solid->AvgIntStrain()
		 <<"\t"<<solid->Xcp()<<"\t"<<jet_vol<<"\t"<<jet_x_vel<<"\t"<<jet_y_vel<<"\t"<<solid->getSphericity()<<'\n';
        //} 
        ofs_.close();
        //amrex::Print()<<"bubble Resolution = "<<solid->getResolution()<<'\n';
    }
}

void incFSI::ComputeMomentumBalance() {
	
    amrex::Print()<<"ComputeMomentumBalance"<<'\n';
    amrex::Real mass_flux = 0.0;
    amrex::Real x_mom_flux_local = 0.0;
    amrex::Real x_mom_flux_global = 0.0;
    amrex::Real y_mom_flux_local = 0.0;
    amrex::Real y_mom_flux_global = 0.0;
    amrex::Real net_traction = 0.0;    

    mass_cv = 0.0;
//Shape of the control volume
//Assuming a rectangular region ( in 2d, a cylindrical region in the axisymmetric geometry)
    double* cv_lo = new double[2]{ -3.0, 0.0 };//x and y coordinate of the lowest corner
    double* cv_hi = new double[2]{ 3.0, 6.0 };//x and y coordinate of the lowest corner


    char *s, *s1, *s2;
    s = new char[50];
    s1 = new char[50];
 
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    int myproc = amrex::ParallelDescriptor::MyProc();


    int fcount = 0;

    amrex::Vector<CFMask> cfmask_(finest_level + 1);
    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        cfmask_[ilev].define(ilev, geom, grids, dmap);
    }

    for (int lev = 0; lev <= finest_level ; lev++)
    {
        const amrex::Real *dx = geom[lev].CellSize();
        const amrex::Real *prob_lo = geom[lev].ProbLo();
        const amrex::Real *prob_hi = geom[lev].ProbHi();
        const amrex::Box &domain = geom[lev].Domain();
        amrex::BoxArray xba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(0));
        amrex::BoxArray yba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(1));

        amrex::iMultiFab *PMask;
        amrex::iMultiFab PMask_NI;//NI stands for No Interface, need a better way of doing this
        if(N_IF > 0)
        {
            auto &&mask__= mask_[lev];
            PMask = &mask__->getPMask();
        }
        else
        {
            PMask_NI.define(grids[lev], dmap[lev], 1, Nghost);
            PMask_NI.setVal(1);
        }        

        int * cv_lo_grid = new int[2]{0,0};
        int * cv_hi_grid = new int[2]{0,0};
        
        //amrex::Print()<<"lev = "<<lev<<" , prob_lo = "<<prob_lo[0]<<" , "<<prob_lo[1]<<" , prob_hi = "<<prob_hi[0]<<" , "<<prob_hi[1]<<" , dx = "<<dx[0]<<" , "<<dx[1]<<"\n";
        cv_lo_grid[0] = convertToGridCoordinate(cv_lo[0], cv_lo[1], lev).first;//gets the integer smaller than the fraction
        cv_lo_grid[1] = convertToGridCoordinate(cv_lo[0], cv_lo[1], lev).second;
        cv_hi_grid[0] = convertToGridCoordinate(cv_hi[0], cv_hi[1], lev).first;
        cv_hi_grid[1] = convertToGridCoordinate(cv_hi[0], cv_hi[1], lev).second; 


        amrex::IntVect cv1,cv2,cv3,cv4;
        cv1.setVal(0,cv_lo_grid[0]);
        cv1.setVal(1,cv_lo_grid[1]);   
        cv2.setVal(0,cv_hi_grid[0]);
        cv2.setVal(1,cv_lo_grid[1]);   
        cv3.setVal(0,cv_hi_grid[0]);
        cv3.setVal(1,cv_hi_grid[1]);   
        cv4.setVal(0,cv_lo_grid[0]);
        cv4.setVal(1,cv_hi_grid[1]);   
        //amrex::Print()<<"cv1 = "<<cv1<<'\n';
        //amrex::Print()<<"cv2 = "<<cv2<<'\n';
        //amrex::Print()<<"cv3 = "<<cv3<<'\n';
        //amrex::Print()<<"cv4 = "<<cv4<<'\n';

        amrex::Box cv(cv1,cv3);
        amrex::Box cv_arm_e(cv2, cv3);
        amrex::Box cv_arm_n(cv4, cv3);
        amrex::Box cv_arm_w(cv1, cv4);
        amrex::Box cv_arm_s(cv1, cv2);   

        //amrex::Print()<<"cv = "<<cv<<'\n';
        //amrex::Print()<<"cv_arm_e = "<<cv_arm_e<<'\n';
        //amrex::Print()<<"cv_arm_w = "<<cv_arm_w<<'\n';
        //amrex::Print()<<"cv_arm_n = "<<cv_arm_n<<'\n';
        //amrex::Print()<<"cv_arm_s = "<<cv_arm_s<<'\n';    

        CFMask cfmask_;
        cfmask_.define(lev, geom, grids, dmap);

        for(amrex::MFIter mfi(Pressure[lev]); mfi.isValid(); ++mfi)
    	{
            fcount++;
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<int const> const &cfmask = cfmask_.Mask().const_array(mfi);
            amrex::Box cv_isect = cv & bx;
            amrex::Box cv_e_isect = cv_arm_e & bx;
            amrex::Box cv_w_isect = cv_arm_w & bx;
            amrex::Box cv_n_isect = cv_arm_n & bx;
            amrex::Box cv_s_isect = cv_arm_s & bx;

            amrex::Array4<amrex::Real const> const &u = xvel[lev].array(mfi);
            amrex::Array4<amrex::Real const> const &v = yvel[lev].array(mfi);
            amrex::Array4<amrex::Real const> const &P = Pressure[lev].const_array(mfi);   
            amrex::Array4<amrex::Real const> phi;         
            if(PhaseField)
                phi = Phi[lev].array(mfi);
            //amrex::Print()<<"cv_arm_e = "<<cv_arm_e<<'\n';
            //amrex::Print()<<"cv_arm_w = "<<cv_arm_w<<'\n';
            //if(cv_e_isect.ok())amrex::Print()<<"cv_e_isect = "<<cv_e_isect<<" ok = "<<cv_e_isect.ok()<<'\n';
            //if(cv_w_isect.ok())amrex::Print()<<"cv_w_isect = "<<cv_w_isect<<" ok = "<<cv_w_isect.ok()<<'\n';
            
            amrex::Array4<int const> pmask;
            if(N_IF > 0)
                pmask = PMask->const_array(mfi);
            else
                pmask = PMask_NI.const_array(mfi);

            if (cv_isect.ok())
            {
                amrex::ParallelFor(cv_isect,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
                {
                    if(cfmask(i,j,k) != CFMask::covered)
                    {
                        amrex::Real u_cc = 0.5*(u(i,j,k) + u(i+1,j,k));
                        amrex::Real v_cc = 0.5*(v(i,j,k) + v(i,j+1,k));
                        amrex::Real area = 2*(22/7)*convertToPhysicalCoordinate(i, j,lev).second*dx[1];

                        if(pmask(i,j,k) == 1)
                        {

                            mass_cv += dx[0]*dx[1]*area;
                        }

                    }
                });
            }
            if (cv_e_isect.ok())
            {
                //! iterate only on mid line
                amrex::ParallelFor(cv_e_isect,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
                {
                    if(cfmask(i,j,k) != CFMask::covered)
                    {
                        amrex::Real alpha = (cv_hi[0] - convertToPhysicalCoordinate(i, j,lev).first)/dx[0];
                        amrex::Real u_face = 0.0;
                        amrex::Real v_face = 0.0;
                        amrex::Real ux = (u(i+1,j,k) - u(i,j,k))/dx[0];
                        amrex::Real uy = (u(i,j + 1,k) - u(i,j,k))/dx[1];
                        amrex::Real vx = (v(i+1,j,k) - v(i,j,k))/dx[0];
                        amrex::Real vy = (v(i,j + 1,k) - v(i,j,k))/dx[1];  
                        amrex::Real Sxx_face = ux;
                        amrex::Real Sxy_face = 0.5*(uy + vx);
                        amrex::Real Syy_face = vy;
                        amrex::Real Mu_face = Mu;                                              
                        amrex::Real phi_face = 0.0;
                        amrex::Real area = dx[1];
                        if(isAxisymmetric)
                            area *= 2*(22/7)*convertToPhysicalCoordinate(i, j,lev).second;
                        if(PhaseField)
                        {
                            phi_face = (1.0 - alpha)*phi(i,j,k) + alpha*phi(i + 1,j,k);
                            Mu_face = visc_.GetViscosity(Sxx_face, Sxy_face, Syy_face, phi_face);
                        }                          

                        if(alpha <= 0.5)
                            u_face = (0.5-alpha)*u(i,j,k)+(0.5+alpha)*u(i+1,j,k);
                        else
                            u_face = (alpha-0.5)*u(i+2,j,k)+(1.5-alpha)*u(i+1,j,k);                        

                        v_face = 0.5*(1.0 - alpha)*(v(i,j,k)+v(i,j+1,k)) + 0.5*alpha*(v(i+1,j,k)+v(i+1,j+1,k));

                        
                        
                        amrex::Real p_face = (1.0 - alpha)*P(i,j,k) + alpha*P(i + 1,j,k);
                        amrex::Real tau_xx_face  = -1.0*p_face + 2.0*Mu_face*Sxx_face;
                        amrex::Real tau_xy_face  =  2.0*Mu_face*Sxy_face;
                        amrex::Real tau_yy_face  = -1.0*p_face + 2.0*Mu_face*Syy_face;                  

                        mass_flux += u_face*area;
                        x_mom_flux_local += u_face*u_face*area;
                        y_mom_flux_local += v_face*u_face*area;
                        //amrex::Print()<<"East Boundary points: "<< i<<" , "<< j<< " , "<<lev<<" , dx = "<<dx[0]<<"\n";
                        //amrex::Print()<<"East Boundary points: "<< convertToPhysicalCoordinate(i, j,lev)<<" , alpha = "<<alpha<<"\n";

                    }
                });
            }
            if (cv_w_isect.ok())
            {
                //! iterate only on mid line
                amrex::ParallelFor(cv_w_isect,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
                {
                    if(cfmask(i,j,k) != CFMask::covered)
                    {
                        amrex::Real alpha = (cv_lo[0] - convertToPhysicalCoordinate(i, j,lev).first)/dx[0];
                        amrex::Real u_face = 0.0;
                        amrex::Real v_face = 0.0;
                        amrex::Real phi_face = 0.0;
                        amrex::Real ux = (u(i+1,j,k) - u(i,j,k))/dx[0];
                        amrex::Real uy = (u(i,j + 1,k) - u(i,j,k))/dx[1];
                        amrex::Real vx = (v(i+1,j,k) - v(i,j,k))/dx[0];
                        amrex::Real vy = (v(i,j + 1,k) - v(i,j,k))/dx[1];
                        amrex::Real Sxx_face = ux;
                        amrex::Real Sxy_face = 0.5*(uy + vx);
                        amrex::Real Syy_face = vy;
                        amrex::Real Mu_face = Mu;
                        amrex::Real area = dx[1];
                        if(isAxisymmetric)
                            area *= 2*(22/7)*convertToPhysicalCoordinate(i, j,lev).second;
                        if(PhaseField)
                        {
                            phi_face = (1.0 - alpha)*phi(i,j,k) + alpha*phi(i + 1,j,k);
                            Mu_face = visc_.GetViscosity(Sxx_face, Sxy_face, Syy_face, phi_face); 
                        }

                        if(alpha <= 0.5)
                            u_face = (0.5-alpha)*u(i,j,k)+(0.5+alpha)*u(i+1,j,k);
                        else
                            u_face = (alpha-0.5)*u(i+2,j,k)+(1.5-alpha)*u(i+1,j,k);                        

                        v_face = 0.5*(1 - alpha)*(v(i,j,k)+v(i,j+1,k)) + 0.5*(alpha)*(v(i+1,j,k)+v(i+1,j+1,k));
                        amrex::Real p_face = (1.0 - alpha)*P(i,j,k) + alpha*P(i + 1,j,k);
                        amrex::Real tau_xx_face  = -1.0*p_face + 2.0*Mu_face*Sxx_face;
                        amrex::Real tau_xy_face  =  2.0*Mu_face*Sxy_face;
                        amrex::Real tau_yy_face  = -1.0*p_face + 2.0*Mu_face*Syy_face; 

                        mass_flux -= u_face*area;
                        x_mom_flux_local -= u_face*u_face*area;
                        y_mom_flux_local -= v_face*u_face*area;

                        //amrex::Print()<<"West Boundary points: "<< i<<" , "<< j<< " , "<<lev<<" , dx = "<<dx[0]<<"\n";
                        //amrex::Print()<<"West Boundary points: "<< convertToPhysicalCoordinate(i, j,lev)<<" , alpha = "<<alpha<<"\n";

                    }
                });
            }

            if (cv_n_isect.ok())
            {
                //! iterate only on mid line
                amrex::ParallelFor(cv_n_isect,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
                {
                    if(cfmask(i,j,k) != CFMask::covered)
                    {
                        amrex::Real alpha = (cv_hi[1] - convertToPhysicalCoordinate(i, j,lev).second)/dx[1];
                        amrex::Real u_face = 0.0;
                        amrex::Real v_face = 0.0;
                        amrex::Real ux = (u(i+1,j,k) - u(i,j,k))/dx[0];
                        amrex::Real uy = (u(i,j + 1,k) - u(i,j,k))/dx[1];
                        amrex::Real vx = (v(i+1,j,k) - v(i,j,k))/dx[0];
                        amrex::Real vy = (v(i,j + 1,k) - v(i,j,k))/dx[1];
                        amrex::Real Sxx_face = ux;
                        amrex::Real Sxy_face = 0.5*(uy + vx);
                        amrex::Real Syy_face = vy;
                        amrex::Real Mu_face = Mu;                                              
                        amrex::Real phi_face = 0.0;
                        amrex::Real area = dx[0];
                        if(isAxisymmetric)
                            area *= 2*(22/7)*convertToPhysicalCoordinate(i, j,lev).second;                        
                        if(PhaseField)
                        {
                            phi_face = (1.0 - alpha)*phi(i,j,k) + alpha*phi(i + 1,j,k);
                            Mu_face = visc_.GetViscosity(Sxx_face, Sxy_face, Syy_face, phi_face);
                        }

                        if(alpha <= 0.5)
                            v_face = (0.5-alpha)*v(i,j,k)+(0.5+alpha)*v(i,j+1,k);
                        else
                            v_face = (alpha-0.5)*v(i,j+2,k)+(1.5-alpha)*v(i,j+1,k);                        

                        u_face = 0.5*(1.0 - alpha)*(u(i,j,k)+u(i+1,j,k)) + 0.5*alpha*(u(i,j+1,k)+u(i+1,j+1,k));
                        amrex::Real p_face = (1.0 - alpha)*P(i,j,k) + alpha*P(i,j+1,k);
                        amrex::Real tau_xx_face  = -1.0*p_face + 2.0*Mu_face*Sxx_face;
                        amrex::Real tau_xy_face  =  2.0*Mu_face*Sxy_face;
                        amrex::Real tau_yy_face  = -1.0*p_face + 2.0*Mu_face*Syy_face;                        

                        mass_flux += v_face*area;
                        x_mom_flux_local += u_face*v_face*area;
                        y_mom_flux_local += v_face*v_face*area;

                         


                        //amrex::Print()<<"north Boundary points: "<< i<<" , "<< j<< " , "<<lev<<" , dx = "<<dx[0]<<"\n";
                        //amrex::Print()<<"north Boundary points: "<< convertToPhysicalCoordinate(i, j,lev)<<" , alpha = "<<alpha<<"\n";

                    }
                });
            }            

            if (cv_s_isect.ok())
            {
                //! iterate only on mid line
                amrex::ParallelFor(cv_s_isect,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
                {
                    if(cfmask(i,j,k) != CFMask::covered)
                    {
                        amrex::Real alpha = (cv_lo[1] - convertToPhysicalCoordinate(i, j,lev).second)/dx[1];
                        amrex::Real u_face = 0.0;
                        amrex::Real v_face = 0.0;
                        amrex::Real ux = (u(i+1,j,k) - u(i,j,k))/dx[0];
                        amrex::Real uy = (u(i,j + 1,k) - u(i,j,k))/dx[1];
                        amrex::Real vx = (v(i+1,j,k) - v(i,j,k))/dx[0];
                        amrex::Real vy = (v(i,j + 1,k) - v(i,j,k))/dx[1];
                        amrex::Real Sxx_face = ux;
                        amrex::Real Sxy_face = 0.5*(uy + vx);
                        amrex::Real Syy_face = vy;
                        amrex::Real Mu_face = Mu;                                              
                        amrex::Real phi_face = 0.0;
                        amrex::Real area = dx[0];
                        if(isAxisymmetric)
                            area *= 2*(22/7)*convertToPhysicalCoordinate(i, j,lev).second;
                        if(PhaseField)
                        {
                            phi_face = (1.0 - alpha)*phi(i,j,k) + alpha*phi(i + 1,j,k);
                            Mu_face = visc_.GetViscosity(Sxx_face, Sxy_face, Syy_face, phi_face);
                        }

                        if(alpha <= 0.5)
                            v_face = (0.5-alpha)*v(i,j,k)+(0.5+alpha)*v(i,j+1,k);
                        else
                            v_face = (alpha-0.5)*v(i,j+2,k)+(1.5-alpha)*v(i,j+1,k);                        

                        u_face = 0.5*(1.0 - alpha)*(u(i,j,k)+u(i+1,j,k)) + 0.5*alpha*(u(i,j+1,k)+u(i+1,j+1,k));
                        amrex::Real p_face = (1.0 - alpha)*P(i,j,k) + alpha*P(i,j+1,k);
                        amrex::Real tau_xx_face  = -1.0*p_face + 2.0*Mu_face*Sxx_face;
                        amrex::Real tau_xy_face  =  2.0*Mu_face*Sxy_face;
                        amrex::Real tau_yy_face  = -1.0*p_face + 2.0*Mu_face*Syy_face; 

                        mass_flux -= v_face*area;
                        x_mom_flux_local -= u_face*v_face*area;
                        y_mom_flux_local -= v_face*v_face*area;

                         


                        //amrex::Print()<<"south Boundary points: "<< i<<" , "<< j<< " , "<<lev<<" , dx = "<<dx[0]<<"\n";
                        //amrex::Print()<<"south Boundary points: "<< convertToPhysicalCoordinate(i, j,lev)<<" , alpha = "<<alpha<<"v_face = "<<v_face<<"\n";

                    }
                });
            }            
        }
    }
    amrex::ParallelDescriptor::ReduceRealSum(mass_flux);
    amrex::ParallelDescriptor::ReduceRealSum(mass_cv);
    amrex::ParallelDescriptor::ReduceRealSum(x_mom_flux_local);
    amrex::ParallelDescriptor::ReduceRealSum(y_mom_flux_local);
    amrex::Real dmdt = (mass_cv - mass_cv_old)/dt;

    auto &bubble = interfaces[finest_level][0];

    amrex::Real dmdt_2 = (bubble->Volume() - bubble->Volume_prev())/dt ;


    //amrex::Print()<<"mass_cv = "<<mass_cv<<" , mass_cv_old = "<<mass_cv_old<<", dmdt = "<<(mass_cv - mass_cv_old)/dt<<'\n';
    //amrex::Print()<<"Mass flux = "<<(mass_cv - mass_cv_old)/dt + mass_flux<<"\n";
    
    {
        char *s, *s1, *s2;
        s = new char[50];
        s1 = new char[50];
        s2 = new char[50];
        sprintf(s, "%2.10f", float(Time));
        strcpy(s1, "Output/CV_analysis.dat");
        //strcat(s1, s);
        // std::string filename(s1);
        std::ofstream ofs_(s1, std::ofstream::out | std::ofstream::app);
        ofs_.flags(std::ios::dec | std::ios::scientific);
        ofs_.precision(7);
        //for (auto &&solid : interfaces)
        //{
            auto &solid = interfaces[finest_level][0];
            amrex::Print(ofs_) << Time <<"\t"<< dt<<"\t"<<mass_cv<<"\t"<<mass_cv_old<<"\t"<<dmdt<<"\t"<<dmdt_2<<"\t"<<mass_flux<<"\t"<<-1.0*dmdt_2 + mass_flux<<"\n";
        //}
        ofs_.close();
    }
    mass_cv_old = mass_cv;
     //if(N_IF > 0)WriteInterface();
     //exit(9);
}

std::pair<double,double> incFSI::convertToPhysicalCoordinate(int i, int j, int lev)
{
    const amrex::Real *dx = geom[lev].CellSize();
    const amrex::Real *prob_lo = geom[lev].ProbLo();
    const amrex::Real *prob_hi = geom[lev].ProbHi();
    return std::make_pair(prob_lo[0] + dx[0] * (i + 0.5),prob_lo[1] + dx[1] * (j + 0.5));
}

std::pair<int,int> incFSI::convertToGridCoordinate(double x, double y, int lev)
{
    const amrex::Real *dx = geom[lev].CellSize();
    const amrex::Real *prob_lo = geom[lev].ProbLo();
    const amrex::Real *prob_hi = geom[lev].ProbHi();
    return std::make_pair(static_cast<int>((x - prob_lo[0])/dx[0] - 0.5 ),
                            static_cast<int>((y - prob_lo[1])/dx[1] - 0.5 ));
}

} // namespace mycode
