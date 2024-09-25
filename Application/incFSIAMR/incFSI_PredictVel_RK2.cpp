#include "incFSI.H"
#include <WeightedENO.h>
#include <AdvectLS_helper.H>

namespace mycode
{

    void incFSI::ComputeIntermediateVelocity_RK2(int lev)
    {

        CFMask cfmask_p_, cfmask_xvel_, cfmask_yvel_;
        //for (int ilev = 0; ilev <= finest_level; ilev++)
        //{
            //amrex::Print()<<"ilev = "<<ilev<<'\n';
        cfmask_p_.define(lev, geom, grids, dmap);
	cfmask_xvel_.define(lev, geom, grids_xvel, dmap);
	cfmask_yvel_.define(lev, geom, grids_yvel, dmap);
        //}


        /// fill internal ghost values
        xvel_old[lev].FillBoundary();
        yvel_old[lev].FillBoundary();
        Pressure[lev].FillBoundary();
        if(TempField)Theta_old[lev].FillBoundary();
        if(PhaseField)Phi_old[lev].FillBoundary();
    
        //! create temporary field data
        
        amrex::BoxArray xba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(0));
        amrex::BoxArray yba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(1));
        amrex::MultiFab tmp_xvel_old(xba, dmap[lev], 1, Nghost);
        amrex::MultiFab tmp_yvel_old(yba, dmap[lev], 1, Nghost);
        amrex::MultiFab tmp_Pressure(grids[lev], dmap[lev], 1, Nghost);
        amrex::MultiFab tmp_Theta_old(grids[lev], dmap[lev], 1, Nghost);
        amrex::MultiFab tmp_Phi_old(grids[lev], dmap[lev], 1, Nghost);
        if (lev > 0)
        {
            //! fill from old current level and coarse level
            //! need to pass the current data also
            amrex::Vector<amrex::MultiFab *> cP(1), cTheta(1), cPhi(1);
            
            amrex::Vector<amrex::Real> ctime(2);
            amrex::Vector<amrex::Real> ftime(2);
    
            amrex::Vector<amrex::MultiFab *> cxvel(1);
            amrex::Vector<amrex::MultiFab *> cyvel(1);
    
            //if (lev > 0)
            {
                ctime[0] = t_old;
                ctime[1] = t_new;
    
                cP[0] = &Pressure[lev - 1];
                if(TempField)cTheta[0] = &Theta_old[lev - 1];
                if(PhaseField)cPhi[0] = &Phi_old[lev - 1];
    
                cxvel[0] = &xvel_old[lev - 1];
                cyvel[0] = &yvel_old[lev - 1];
            }
            ftime[0] = t_old;
            ftime[1] = t_new;
    
            if(divergence_free_interpolation)
            {
                amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> tmp_vel_old = {&tmp_xvel_old, &tmp_yvel_old};
                amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> cvel;
                amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> fvel;
    
                cvel.resize(1);
                fvel.resize(1);
                cvel[0] = {&xvel_old[lev - 1], &yvel_old[lev - 1]};
                fvel[0] = {&xvel_old[lev], &yvel_old[lev]};
    
                FillPatch(lev, Time, tmp_vel_old, cvel, {Time}, fvel, {Time}, 0, 1);
    
                amrex::MultiFab::Copy(tmp_xvel_old, *tmp_vel_old[0], 0, 0, 1, Nghost);
                amrex::MultiFab::Copy(tmp_yvel_old, *tmp_vel_old[1], 0, 0, 1, Nghost);
            }
            else
            {
                FillPatch(lev, Time, tmp_xvel_old, {cxvel[0]}, {ctime[0]}, {&xvel_old[lev]}, {ftime[0]}, 0, 1);
                FillPatch(lev, Time, tmp_yvel_old, {cyvel[0]}, {ctime[0]}, {&yvel_old[lev]}, {ftime[0]}, 0, 1);
            }
    
            //! no need of time interpolation of following
            FillPatch(lev, Time, tmp_Pressure, cP, {ctime[0]}, {&Pressure[lev]}, {ftime[0]}, 0, 1);
            if(TempField) FillPatch(lev, Time, tmp_Theta_old, cTheta, {ctime[0]}, {&Theta_old[lev]}, {ftime[0]}, 0, 1);
            if(PhaseField) FillPatch(lev, Time, tmp_Phi_old, cPhi, {ctime[0]}, {&Phi_old[lev]}, {ftime[0]}, 0, 1);
    
    
            //! need to fill physical bcs
            XVelBoundaryConditions(lev, Time, tmp_xvel_old);
            YVelBoundaryConditions(lev, Time, tmp_yvel_old);
            PressureBoundaryConditions(lev, Time, tmp_Pressure);
            if(TempField) TemperatureBoundaryConditions(lev, Time, tmp_Theta_old);
            if(PhaseField) PhiBoundaryConditions(lev, Time, tmp_Phi_old);
        }
        else
        {
            amrex::MultiFab::Copy(tmp_xvel_old, xvel_old[lev], 0, 0, 1, Nghost);
            amrex::MultiFab::Copy(tmp_yvel_old, yvel_old[lev], 0, 0, 1, Nghost);
            amrex::MultiFab::Copy(tmp_Pressure, Pressure[lev], 0, 0, 1, Nghost);
            if(TempField) amrex::MultiFab::Copy(tmp_Theta_old, Theta_old[lev], 0, 0, 1, Nghost);
            if(PhaseField) amrex::MultiFab::Copy(tmp_Phi_old, Phi_old[lev], 0, 0, 1, Nghost);
        }
        
        const amrex::Real *dx = geom[lev].CellSize();
        const amrex::Real *prob_lo = geom[lev].ProbLo();
    
        amrex::iMultiFab *UMask;
        amrex::iMultiFab *VMask;
        amrex::iMultiFab *PMask;
        amrex::iMultiFab *Mask;

        amrex::iMultiFab UMask_NI;//NI stands for No Interface, need a better way of doing this
        amrex::iMultiFab VMask_NI;//NI stands for No Interface, need a better way of doing this
        amrex::iMultiFab PMask_NI;//NI stands for No Interface, need a better way of doing this
        amrex::iMultiFab Mask_NI;

        if(N_IF > 0 && lev == finest_level)
        {
            auto &&mask__= mask_[lev];
            UMask = &mask__->getUMask();
            VMask = &mask__->getVMask();
            PMask = &mask__->getPMask();
            Mask = &mask__->getMask();
        }
        else
        {
            UMask_NI.define(xba, dmap[lev], 1, Nghost);
            VMask_NI.define(yba, dmap[lev], 1, Nghost);
            PMask_NI.define(grids[lev], dmap[lev], 1, Nghost);
            Mask_NI.define(grids[lev], dmap[lev], 1, Nghost);

            UMask_NI.setVal(1);
            VMask_NI.setVal(1);
            PMask_NI.setVal(1);
            Mask_NI.setVal(1);
        }
    
        Mu_max = 0.0;    
        /// compute x-component
        for(amrex::MFIter mfi(xvel[lev]); mfi.isValid(); ++mfi)
        {
	    amrex::Array4<int const> const &cfmask = cfmask_xvel_.Mask().const_array(mfi);
            const amrex::Box& bx = mfi.validbox();
    
            amrex::Array4<amrex::Real const> const &u = tmp_xvel_old.const_array(mfi);
            amrex::Array4<amrex::Real const> const &v = tmp_yvel_old.const_array(mfi);
            amrex::Array4<amrex::Real const> const &P = tmp_Pressure.const_array(mfi);
            amrex::Array4<amrex::Real> const &ustar = xvel[lev].array(mfi);
            amrex::Array4<amrex::Real> const &uold = FRK_xvel[lev].array(mfi);
            amrex::Box box = get_valid_face_box(lev, bx, 0);
    
            amrex::Array4<int> umask;
    
            if(N_IF == 0 || lev != finest_level)
                umask = UMask_NI.array(mfi);
            else
                umask = UMask->array(mfi);

            amrex::Array4<amrex::Real const> phi, dmg;
            if(PhaseField)
                phi = tmp_Phi_old.array(mfi);
            if(DamageModel)
                dmg = Scalars_old[lev][8].array(mfi);

    
            amrex::ParallelFor(box,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                if(umask(i, j, k) == 1 && cfmask(i,j,k) != CFMask::covered)
                {
                    amrex::Real ue = 0.5 * (u(i, j, k) + u(i + 1, j, k));
                    amrex::Real ux_e = (u(i + 1, j, k) - u(i, j, k)) / dx[0];
                    amrex::Real uy_e = 0.25 * ((u(i, j + 1, k) + u(i + 1, j + 1, k))
                                              -(u(i, j - 1, k) + u(i + 1, j - 1, k)))/ dx[1];
                    amrex::Real vx_e = 0.25*( (v(i + 1, j, k) + v( i + 1, j + 1, k))
                                            -(v(i - 1, j, k) + v( i - 1, j + 1, k)))/ dx[0];
                    amrex::Real vy_e = (v(i, j + 1, k) - v(i, j, k))/dx[1];
                    amrex::Real S11_e = ux_e;
                    amrex::Real S12_e = 0.5*(uy_e + vx_e);
                    amrex::Real S22_e = vy_e;
                    amrex::Real phi_e,dmg_e;
                    if(PhaseField)
                        phi_e = phi(i, j, k);
                    else
                        phi_e = 0.0;
                    if(DamageModel)
                        dmg_e = dmg(i, j, k);
                    else
                        dmg_e = 0.0;
                    amrex::Real Mu_e = Mu;
                    if(PhaseField)
                        Mu_e = visc_.GetViscosity(S11_e, S12_e, S22_e, phi_e, dmg_e);
                    amrex::Real fe = ue * ue - 2.0 * Mu_e * ux_e;

                    
                    amrex::Real uw = 0.5 * (u(i, j, k) + u(i - 1, j, k));
                    amrex::Real ux_w = (u(i, j, k) - u(i - 1, j, k)) / dx[0];
                    amrex::Real uy_w = 0.25 * ((u(i - 1, j + 1, k) + u(i , j + 1, k))
                                             -(u(i - 1, j - 1, k) + u(i , j - 1, k)))/ dx[1];
                    amrex::Real vx_w = 0.25*( (v(i , j, k) + v(i , j + 1, k))
                                           -(v(i - 2, j, k) + v( i - 2, j + 1, k)))/ dx[0];
                    amrex::Real vy_w = (v(i - 1, j + 1, k) - v(i - 1, j, k))/dx[1];
                    amrex::Real S11_w = ux_w;
                    amrex::Real S12_w = 0.5*(uy_w + vx_w);
                    amrex::Real S22_w = vy_w;
                    amrex::Real phi_w, dmg_w;
                    if(PhaseField)
                        phi_w = phi(i - 1, j, k);
                    else
                        phi_w = 0.0;
                    if(DamageModel)
                        dmg_w = dmg(i - 1, j, k);
                    else
                        dmg_w = 0.0;
                    amrex::Real Mu_w = Mu;
                    if(PhaseField)
                        Mu_w = visc_.GetViscosity(S11_w, S12_w, S22_w, phi_w, dmg_w);
                    amrex::Real fw = uw * uw - 2.0 * Mu_w * ux_w;

                    
                    amrex::Real un = 0.5 * (u(i, j, k) + u(i, j + 1, k));
                    amrex::Real vn = 0.5 * (v(i - 1, j + 1, k) + v(i, j + 1, k));
                    amrex::Real ux_n = 0.25 * ( (u(i+1 ,j ,k) + u(i+1 ,j+1 ,k))
                                             -(u(i-1 ,j ,k) + u(i-1 ,j+1 ,k)) )/dx[0];
                    amrex::Real uy_n = (u(i, j + 1, k) - u(i, j, k)) / dx[1];
                    amrex::Real vx_n = (v(i, j + 1, k) - v(i - 1, j + 1, k)) / dx[0];
                    amrex::Real vy_n = 0.25 * ( ((v(i ,j + 2 ,k) + v(i - 1 ,j + 2 ,k))
                                             -(v(i ,j ,k ) + v( i - 1 , j, k))) )/dx[1];
                    amrex::Real S11_n = ux_n;
                    amrex::Real S12_n = 0.5*(uy_n + vx_n);
                    amrex::Real S22_n = vy_n;
                    amrex::Real phi_n,dmg_n;
                    if(PhaseField)
                        phi_n = 0.25*(phi(i, j, k) + phi(i, j+1, k) + phi(i+1, j, k) + phi(i+1, j+1, k));
                    else
                        phi_n = 0.0;
                    if(DamageModel)
                        dmg_n = 0.25*(dmg(i, j, k) + dmg(i, j+1, k) + dmg(i+1, j, k) + dmg(i+1, j+1, k));
                    else
                        dmg_n = 0.0;
                    amrex::Real Mu_n = Mu;
                    if(PhaseField)
                        Mu_n = visc_.GetViscosity(S11_n, S12_n, S22_n, phi_n, dmg_n);
                    amrex::Real fn = un * vn - Mu_n * (uy_n + vx_n);

                    
                    amrex::Real us = 0.5 * (u(i, j, k) + u(i, j - 1, k));
                    amrex::Real vs = 0.5 * (v(i - 1, j, k) + v(i, j, k));
                    amrex::Real ux_s = 0.25 * ( (u(i + 1 , j, k) + u(i + 1 , j - 1, k))
                                             -(u(i - 1 , j, k) + u(i - 1 , j - 1, k)))/dx[0];
                    amrex::Real uy_s = (u(i, j, k) - u(i, j - 1, k)) / dx[1];
                    amrex::Real vx_s = (v(i, j, k) - v(i - 1, j, k)) / dx[0];
                    amrex::Real vy_s = 0.25 * ( (v(i , j + 1 ,k) + v(i - 1, j + 1 , k))
                                             -(v(i , j - 1 ,k) + v(i - 1, j - 1 , k)))/dx[1];
                    amrex::Real S11_s = ux_s;
                    amrex::Real S12_s = 0.5*(uy_s + vx_s);
                    amrex::Real S22_s = vy_s;
                    amrex::Real phi_s, dmg_s;
                    if(PhaseField)
                        phi_s = 0.25*(phi(i, j, k) + phi(i, j - 1, k) + phi(i - 1, j, k) + phi(i - 1, j - 1, k));
                    else
                        phi_s = 0.0;
                    if(DamageModel)
                        dmg_s = 0.25*(dmg(i, j, k) + dmg(i, j - 1, k) + dmg(i - 1, j, k) + dmg(i - 1, j - 1, k));
                    else
                        dmg_s = 0.0;
                    amrex::Real Mu_s = Mu;
                    if(PhaseField)
                        Mu_s = visc_.GetViscosity(S11_s, S12_s, S22_s, phi_s, dmg_s); 
                    amrex::Real fs = us * vs - Mu_s * (uy_s + vx_s);

                    Mu_max = std::max(0.25*(Mu_e+Mu_w+Mu_n+Mu_s), Mu_max);
		    //if(i == 2432 && j == 95 && lev == 6)
		    /*if(Mu_max > 2.0)
	            {
		        amrex::Print()<<"lev = "<<lev<<" , i  = "<<i<<" , j = "<<j<<'\n';
			amrex::Print()<<"Mu = "<<Mu_max<<'\n';
			amrex::Print()<<" Mu_e = "<<Mu_e<<" Mu_w = "<<Mu_w<<" Mu_n = "<<Mu_n<<" Mu_s = "<<Mu_s<<'\n';
			amrex::Print()<<"umask(i, j, k) = "<<umask(i, j, k)<<'\n';
			exit(1);
		    }*/
                    if (isAxisymmetric)
                    {
                        //TODO how to avoid singularity? i.e. what if y = 0
                        //! box is centered in y
                        amrex::Real yn = prob_lo[1] + dx[1] * (j + 1);
                        amrex::Real y  = prob_lo[1] + dx[1] * (j + 0.5);
                        amrex::Real ys = prob_lo[1] + dx[1] * j;
                    
                        fn *= yn;
                        fs *= ys;
                        ustar(i, j, k) = uold(i, j, k) +  dt * ( (fw - fe) / dx[0] + (fs - fn) / (y * dx[1]) + gravity);                    
                    }
                    else
                    {
                        ustar(i, j, k) = uold(i, j, k) +  dt * ( (fw - fe) / dx[0] + (fs - fn) / dx[1] + gravity);
                    }
                    
                    if (Projection)
                    {
                        amrex::Real dpdx = (P(i, j, k) - P(i - 1, j, k)) / dx[0];
                        ustar(i, j, k) -= dt * dpdx;
                    }
                }
            });
        }
    
        /// compute y-component
        for(amrex::MFIter mfi(yvel[lev]); mfi.isValid(); ++mfi)
        {
	    amrex::Array4<int const> const &cfmask = cfmask_yvel_.Mask().const_array(mfi);
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<amrex::Real const> const &u = tmp_xvel_old.const_array(mfi);
            amrex::Array4<amrex::Real const> const &v = tmp_yvel_old.const_array(mfi);
            amrex::Array4<amrex::Real const> const &P = tmp_Pressure.const_array(mfi);
            amrex::Array4<amrex::Real> const &vstar = yvel[lev].array(mfi);
            amrex::Array4<amrex::Real> const &vold = FRK_yvel[lev].array(mfi);
            amrex::Box box = get_valid_face_box(lev, bx, 1);
    
            amrex::Array4<int> vmask;
            if(N_IF == 0 || lev != finest_level)
                vmask = VMask_NI.array(mfi);
            else
                vmask = VMask->array(mfi);
    
            amrex::Array4<amrex::Real const> phi, dmg;
            if(PhaseField)
                phi = tmp_Phi_old.array(mfi);
            if(DamageModel)
                dmg = Scalars_old[lev][8].array(mfi);
    
            amrex::ParallelFor(box,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                if(vmask(i ,j ,k) == 1 && cfmask(i,j,k) != CFMask::covered)
                {
                    amrex::Real ue = 0.5 * (u(i + 1, j, k) + u(i + 1, j - 1, k));
                    amrex::Real ve = 0.5 * (v(i, j, k) + v(i + 1, j, k));
                    amrex::Real ux_e = 0.25 * ( (u(i + 2, j, k) + u(i + 2, j - 1, k))
                                               -(u(i, j, k) + u(i ,j - 1, k)))/dx[0];
                    amrex::Real uy_e = (u(i + 1, j, k) - u(i + 1, j - 1, k)) / dx[1];
                    amrex::Real vx_e = (v(i + 1, j, k) - v(i, j, k)) / dx[0];
                    amrex::Real vy_e = 0.25 * ( (v(i , j+1 ,k) + v(i+1 ,j+1 ,k))
                                               -(v(i, j-1, k) + v(i+1 ,j-1 ,k)))/dx[1];
                    amrex::Real S11_e = ux_e;
                    amrex::Real S12_e = 0.5*(uy_e + vx_e);
                    amrex::Real S22_e = vy_e;
                    amrex::Real phi_e, dmg_e;
                    if(PhaseField)
                        phi_e = 0.25 * (phi(i, j, k) + phi(i + 1, j, k) + phi(i, j - 1, k) + phi(i + 1, j - 1, k));
                    else
                        phi_e = 0.0;
                    if(DamageModel)
                        dmg_e = 0.25 * (dmg(i, j, k) + dmg(i + 1, j, k) + dmg(i, j - 1, k) + dmg(i + 1, j - 1, k));
                    else
                        dmg_e = 0.0;
                    amrex::Real Mu_e = Mu;
                    if(PhaseField)
                        Mu_e = visc_.GetViscosity(S11_e, S12_e, S22_e, phi_e, dmg_e);
                    amrex::Real fe = ue * ve - Mu_e * (uy_e + vx_e);
                    
                    amrex::Real uw = 0.5 * (u(i, j, k) + u(i, j - 1, k));
                    amrex::Real vw = 0.5 * (v(i, j, k) + v(i - 1, j, k));
                    amrex::Real ux_w = 0.25 * ( (u(i + 1, j, k) + u(i + 1, j - 1, k))
                                               -(u(i - 1, j, k) + u(i - 1,j - 1, k)))/dx[0];
                    amrex::Real uy_w = (u(i, j, k) - u(i, j - 1, k)) / dx[1];
                    amrex::Real vx_w = (v(i, j, k) - v(i - 1, j, k)) / dx[0];
                    amrex::Real vy_w = 0.25 * ( (v(i, j + 1, k) + v(i - 1, j + 1 ,k))
                                               -(v(i, j - 1, k) + v(i - 1, j - 1 ,k)))/dx[1];
                    amrex::Real S11_w = ux_w;
                    amrex::Real S12_w = 0.5*(uy_w + vx_w);
                    amrex::Real S22_w = vy_w;
                    amrex::Real phi_w, dmg_w;
                    if(PhaseField)
                        phi_w = 0.25 * (phi(i - 1, j, k) + phi(i, j, k) + phi(i - 1, j - 1, k) + phi(i, j - 1, k));
                    else
                        phi_w = 0.0;
                    if(DamageModel)
                        dmg_w = 0.25 * (dmg(i - 1, j, k) + dmg(i, j, k) + dmg(i - 1, j - 1, k) + dmg(i, j - 1, k));
                    else
                        dmg_w = 0.0;
                    amrex::Real Mu_w = Mu;
                    if(PhaseField)
                        Mu_w = visc_.GetViscosity(S11_w, S12_w, S22_w, phi_w, dmg_w);
                    amrex::Real fw = uw * vw - Mu_w * (uy_w + vx_w);


                    amrex::Real vn = 0.5 * (v(i, j, k) + v(i, j + 1, k));
                    amrex::Real ux_n = (u(i + 1, j, k) - u(i, j, k))/dx[0];
                    amrex::Real uy_n = 0.25 * ( (u(i + 1, j + 1, k) + u(i, j + 1, k))
                                               -(u(i + 1, j - 1, k) + u(i, j - 1, k)))/dx[1];
                    amrex::Real vx_n = 0.25 * ( (v(i + 1, j + 1, k) + v(i + 1, j, k))
                                               -(v(i - 1, j + 1, k) + v(i - 1, j, k)))/dx[0];
                    amrex::Real vy_n = (v(i, j + 1, k) - v(i, j, k)) / dx[1];
                    amrex::Real S11_n = ux_n;
                    amrex::Real S12_n = 0.5*(uy_n + vx_n);
                    amrex::Real S22_n = vy_n;
                    amrex::Real phi_n, dmg_n;
                    if(PhaseField)
                        phi_n = phi(i, j, k);
                    else
                        phi_n = 0.0;
                    if(DamageModel)
                        dmg_n = dmg(i, j, k);
                    else
                        dmg_n = 0.0;
                    amrex::Real Mu_n = Mu;
                    if(PhaseField)
                        Mu_n = visc_.GetViscosity(S11_n, S12_n, S22_n, phi_n, dmg_n);
                    amrex::Real fn = vn * vn - 2.0 * Mu_n * vy_n;
                    
                    amrex::Real vs = 0.5 * (v(i, j, k) + v(i, j - 1, k));
                    amrex::Real ux_s = (u(i + 1, j - 1, k) - u(i, j - 1, k))/dx[0];
                    amrex::Real uy_s = 0.25 * ( (u(i + 1, j, k) + u(i, j, k))
                                               -(u(i + 1, j - 2, k) + u(i, j - 2, k)))/dx[1];
                    amrex::Real vx_s = 0.25 * ( (v(i + 1, j, k) + v(i + 1, j - 1, k))
                                               -(v(i - 1, j, k) + v(i - 1, j - 1, k)))/dx[0];
                    amrex::Real vy_s = (v(i, j, k) - v(i, j - 1, k)) / dx[1];
                    amrex::Real S11_s = ux_s;
                    amrex::Real S12_s = 0.5*(uy_s + vx_s);
                    amrex::Real S22_s = vy_s;
                    amrex::Real phi_s, dmg_s;
                    if(PhaseField)
                        phi_s = phi(i, j - 1, k);
                    else
                        phi_s = 0.0;
                    if(DamageModel)
                        dmg_s = dmg(i, j - 1, k);
                    else
                        dmg_s = 0.0;
                    amrex::Real Mu_s = Mu;
                    if(PhaseField)
                        Mu_s = visc_.GetViscosity(S11_s, S12_s, S22_s, phi_s, dmg_s);
                    amrex::Real fs = vs * vs - 2.0 * Mu_s * vy_s;

                    Mu_max = std::max(0.25*(Mu_e+Mu_w+Mu_n+Mu_s), Mu_max);
                    if (isAxisymmetric)
                    {
                        //! box is nodal in y
                        amrex::Real yn = prob_lo[1] + dx[1] * (j + 0.5);
                        amrex::Real y = prob_lo[1] + dx[1] * j;
                        amrex::Real ys = prob_lo[1] + dx[1] * (j - 0.5);
                        amrex::Real Mu_ = 0.25*(Mu_e+Mu_w+Mu_n+Mu_s);
                    
                        fn *= yn;
                        fs *= ys;
                    
                        //! avoid singularity
                        if (fabs(y) <  1e-10)
                        {
                            vstar(i, j, k) = 0.0;
                        }
                        else
                        {
                            vstar(i, j, k) = vold(i ,j ,k) + dt * ((fw - fe) / dx[0] + (fs - fn) / (y * dx[1]) - 2.0 * Mu_ * v(i, j, k) / (y * y));
                        }
                    }
                    else
                    {
                        vstar(i, j, k) = vold(i ,j ,k) +  dt * ((fw - fe) / dx[0] + (fs - fn) / dx[1]);
                    }
                    
                    if (Projection)
                    {
                        amrex::Real dpdy = (P(i, j, k) - P(i, j - 1, k)) / dx[1];
                        vstar(i, j, k) -= dt * dpdy;
                    }
                }
            });
        }
    
    
        if(TempField)
        {
            /// compute Theta
            amrex::Real k_;
            k_max = 0.0;
            
            /// compute Temperature
            for(amrex::MFIter mfi(Theta[lev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<amrex::Real const> const &xvel = tmp_xvel_old.const_array(mfi);
                amrex::Array4<amrex::Real const> const &yvel = tmp_yvel_old.const_array(mfi);
                amrex::Array4<amrex::Real const> const &P = tmp_Pressure.const_array(mfi);
                amrex::Array4<amrex::Real const> const &Told = tmp_Theta_old.const_array(mfi);
                amrex::Array4<amrex::Real const> const &T_frk = FRK_Theta[lev].const_array(mfi);
                amrex::Array4<amrex::Real> const &T = Theta[lev].array(mfi);
        
                amrex::Array4<int> mask;
                if(N_IF == 0 || lev != finest_level)
                    mask = PMask_NI.array(mfi);
                else
                    mask = PMask->array(mfi);

		amrex::Array4<int const> const &cfmask = cfmask_p_.Mask().const_array(mfi);

                amrex::Array4<amrex::Real const> phi, dmg;
                if(PhaseField)
                    phi = tmp_Phi_old.array(mfi);
                if(DamageModel)
                    dmg = Scalars_old[lev][8].array(mfi);

        
        
                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
                {
                    if(mask(i ,j ,k) == 1 && cfmask(i,j,k) != CFMask::covered)
                    {

                        amrex::Real ux = (xvel(i + 1, j, k) - xvel(i, j, k)) / dx[0];
                        amrex::Real uy = (xvel(i + 1, j + 1, k) - xvel(i + 1, j - 1, k) + xvel(i, j + 1, k) - xvel(i, j - 1, k)) / (4.0 * dx[1]);
                        amrex::Real vx = (yvel(i + 1, j + 1, k) - yvel(i - 1, j + 1, k) + yvel(i + 1, j, k) - yvel(i - 1, j, k)) / (4.0 * dx[0]);
                        amrex::Real vy = (yvel(i, j + 1, k) - yvel(i, j, k)) / dx[1];
                        amrex::Real S11_ = ux;
                        amrex::Real S12_ = 0.5*(uy + vx);
                        amrex::Real S22_ = vy;

                        amrex::Real dGamma_dt = visc_.ComputeGammaDot(ux, 0.5 * (uy + vx), vy);
                        amrex::Real phi_, dmg_;
                        if(PhaseField)
                            phi_ = phi(i, j, k);
                        else
                            phi_ = 0.0;
                        if(DamageModel)
                            dmg_ = dmg(i, j, k);
                        else
                           dmg_ = 0.0;
                        amrex::Real Mu_ = visc_.GetViscosity(S11_, S12_, S22_, phi_, dmg_);

                        amrex::Real theta_e = 0.5 * (Told(i, j, k) + Told(i + 1, j, k));
                        amrex::Real thetax_e = (Told(i + 1, j, k) - Told(i, j, k)) / dx[0];
                        k_ = visc_.GetConductivity(Mu_, Prandtl_no);
                        if (mask(i,j,k) == 1 && mask(i+1,j,k) == 1)
                        {
                            k_max = std::max(k_, k_max);
                        }
                        amrex::Real fe = theta_e * xvel(i + 1, j, k) - k_ * thetax_e;
        
                        amrex::Real theta_w = 0.5 * (Told(i, j, k) + Told(i - 1, j, k));
                        amrex::Real thetax_w = (Told(i, j, k) - Told(i - 1, j, k)) / dx[0];
                        k_ = visc_.GetConductivity(Mu_, Prandtl_no);
                        if (mask(i, j, k) == 1 && mask(i - 1, j, k) == 1)
                        {
                            k_max = std::max(k_, k_max);
                        }
                        amrex::Real fw = theta_w * xvel(i, j, k) - k_ * thetax_w;
                        
                        amrex::Real theta_n = 0.5 * (Told(i, j, k) + Told(i, j + 1, k));
                        amrex::Real thetay_n = (Told(i, j + 1, k) - Told(i, j, k)) / dx[1];
                        k_ = visc_.GetConductivity(Mu_, Prandtl_no);
                        if (mask(i, j, k) == 1 && mask(i, j + 1, k) == 1)
                        {
                            k_max = std::max(k_, k_max);
                        }
                        amrex::Real fn = theta_n * yvel(i, j + 1, k) - k_ * thetay_n;
        
                        amrex::Real theta_s = 0.5 * (Told(i, j, k) + Told(i, j - 1, k));
                        amrex::Real thetay_s = (Told(i, j, k) - Told(i, j - 1, k)) / dx[1];
                        k_ = visc_.GetConductivity(Mu_, Prandtl_no);
                        if (mask(i, j, k) == 1 && mask(i, j - 1, k) == 1)
                        {
                            k_max = std::max(k_, k_max);
                        }
                        amrex::Real fs = theta_s * yvel(i, j, k) - k_ * thetay_s;
        
                        if (isAxisymmetric)
                        {
                            //! box is cell centered
                            amrex::Real yn = prob_lo[1] + dx[1] * (j + 1);
                            amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                            amrex::Real ys = prob_lo[1] + dx[1] * j;
        
                            fn *= yn;
                            fs *= ys;
                            T(i, j, k) = T_frk(i, j, k) + dt * ((fw - fe) / dx[0] + (fs - fn) / (y * dx[1]) +  Mu_ * (u_ref*u_ref/cp_/T_ref) * dGamma_dt* dGamma_dt);
                        }
                        else
                        {
                            T(i, j, k) = T_frk(i, j, k) +  dt * ((fw - fe) / dx[0] + (fs - fn) / dx[1] +  Mu_ * (u_ref*u_ref/cp_/T_ref) * dGamma_dt * dGamma_dt);
                        }
                    }
                    if(mask(i,j,k) != 1 && mask(i,j,k) != 2)
                    {
                        //auto &solid = interfaces[0];
                        T(i, j, k) = 0.0;//P_interface[0] * solid->Volume() / (P_interface0[0] * solid->Volume_0());
                        //amrex::Print()<<"T(i, j, k) = "<<T(i, j, k)<<'\n';
                    }
                });
            }    
        }
        
        /// Advect Phasefield
        if(PhaseField)
        {
            amrex::Real alpha_ = 1.0;
            for(amrex::MFIter mfi(Phi[lev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<amrex::Real const> const &xvel = tmp_xvel_old.const_array(mfi);
                amrex::Array4<amrex::Real const> const &yvel = tmp_yvel_old.const_array(mfi);
                amrex::Array4<amrex::Real const> const &phiold = tmp_Phi_old.const_array(mfi);
                amrex::Array4<amrex::Real const> const &phi_frk = FRK_Phi[lev].const_array(mfi);
                amrex::Array4<amrex::Real> const &phi = Phi[lev].array(mfi);
            
                amrex::Array4<int> pmask, mask;
                if(N_IF == 0 || lev != finest_level)
                    pmask = PMask_NI.array(mfi);
                else
                    pmask = PMask->array(mfi);

                if(N_IF == 0 || lev != finest_level)
                    mask = Mask_NI.array(mfi);
                else
                    mask = Mask->array(mfi);

                amrex::Array4<amrex::Real const> dmg;
                if(DamageModel)
                    dmg = Scalars_old[lev][dmg_num].array(mfi);

                amrex::Array4<int const> const &cfmask = cfmask_p_.Mask().const_array(mfi); 
            
                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
                {
                    if(pmask(i ,j ,k) == 1 && cfmask(i,j,k) != CFMask::covered)
                    {

                        amrex::Real Phix_L, Phix_R,Phiy_L, Phiy_R, Phix, Phiy;
                        amrex::Real u_Mid = 0.5 * (xvel(i, j, k) + xvel(i + 1, j, k));
                        amrex::Real v_Mid = 0.5 * (yvel(i, j, k) + yvel(i, j + 1, k));

                        WENO5_LS(Phix_L, Phix_R, Phiy_L, Phiy_R, i, j, dx, phiold);
                        if(mask(i, j, k) == 3)
                        {
                            Phix_R = (phiold(i + 1, j, k) - phiold(i, j, k)) / dx[0];
                            Phix_L = (phiold(i, j, k) - phiold(i - 1, j, k)) / dx[0];
                            Phiy_R = (phiold(i, j + 1, k) - phiold(i, j, k)) / dx[1];
                            Phiy_L = (phiold(i, j, k) - phiold(i, j - 1, k)) / dx[1];
                        }
                        if (u_Mid > 0.0)
                            Phix = Phix_L;
                        else
                            Phix = Phix_R;
                        if (v_Mid > 0.0)
                            Phiy = Phiy_L;
                        else
                            Phiy = Phiy_R;


                        /*
                        amrex::Real phi_e = 0.5 * (phiold(i, j, k) + phiold(i + 1, j, k));
                        //amrex::Real fe = phiold(i, j, k) * xvel(i + 1, j, k); 
                        //amrex::Real fe = phi_e * xvel(i + 1, j, k);
                        amrex::Real fe = phi_e * xvel(i + 1, j, k) 
                                         + alpha_ * 0.5 * std::abs(xvel(i + 1, j, k)) * (phiold(i, j, k) - phiold(i + 1, j, k));
            
            
                        amrex::Real phi_w = 0.5 * (phiold(i, j, k) + phiold(i - 1, j, k));
                        //amrex::Real fw = phi_w * xvel(i, j, k);
                        amrex::Real fw = phi_w * xvel(i, j, k) 
                                         + alpha_ * 0.5 * std::abs(xvel(i, j, k)) * (phiold(i - 1, j, k) - phiold(i, j, k));
            
                        
                        amrex::Real phi_n = 0.5 * (phiold(i, j, k) + phiold(i, j + 1, k));
                        //amrex::Real fn = phi_n * yvel(i, j + 1, k);
                        amrex::Real fn = phi_n * yvel(i, j + 1, k) 
                                         + alpha_ * 0.5 * std::abs(yvel(i, j + 1, k)) * (phiold(i, j, k) - phiold(i, j + 1, k));
            
            
                        amrex::Real phi_s = 0.5 * (phiold(i, j, k) + phiold(i, j - 1, k));
                        //amrex::Real fs = phi_s * yvel(i, j, k);
                        amrex::Real fs = phi_s * yvel(i, j, k) 
                                         + alpha_ * 0.5 * std::abs(yvel(i, j, k)) * (phiold(i, j - 1, k) - phiold(i, j, k));
            
                        amrex::Real ux = (xvel(i + 1, j, k) - xvel(i, j, k)) / dx[0];
                        amrex::Real uy = (xvel(i + 1, j + 1, k) - xvel(i + 1, j - 1, k) + xvel(i, j + 1, k) - xvel(i, j - 1, k)) / (4.0 * dx[1]);
                        amrex::Real vx = (yvel(i + 1, j + 1, k) - yvel(i - 1, j + 1, k) + yvel(i + 1, j, k) - yvel(i - 1, j, k)) / (4.0 * dx[0]);
                        amrex::Real vy = (yvel(i, j + 1, k) - yvel(i, j, k)) / dx[1];
                        amrex::Real dGamma_dt = visc_.ComputeGammaDot(ux, 0.5 * (uy + vx), vy);
                        amrex::Real Mu_ = visc_.GetViscosity(dGamma_dt);
                        //Mu_ = Mu;
                        */
                        //RHS for interface sharpening
                        amrex::Real RHS_phi = 0.0;

                        {
                            amrex::Real b = b_sharp;
                            amrex::Real laplacian_phi = (1.0/(3.0*dx[0]*dx[0]))*
                                                        (2.0*(phiold(i + 1, j, k) + phiold(i, j + 1, k) + phiold(i - 1, j, k) + phiold(i, j - 1, k ) - 4.0*phiold(i, j, k)) +
                                                         0.5*(phiold(i + 1, j + 1, k) + phiold(i + 1, j - 1, k) + phiold(i - 1, j + 1, k) + phiold(i - 1, j - 1, k) - 4.0*phiold(i, j, k)));
                            amrex::Real grad_phi = (1.0/dx[0])*std::hypot(0.5*(phiold(i + 1, j, k) - phiold(i - 1, j, k)), 0.5*(phiold(i, j + 1, k) - phiold(i, j - 1, k)));
                            amrex::Real curvat = 0.0;
                            curvat = (phiold(i+1,j,k) - phiold(i,j,k))/
                                    (std::hypot(phiold(i+1,j,k) - phiold(i,j,k),(phiold(i+1,j+1,k)+phiold(i,j+1,k)-phiold(i+1,j-1,k)-phiold(i,j-1,k))/4.0) + 1e-15);
                            curvat -= (phiold(i,j,k) - phiold(i-1,j,k))/
                                    (std::hypot(phiold(i,j,k) - phiold(i-1,j,k),(phiold(i-1,j+1,k)+phiold(i,j+1,k)-phiold(i-1,j-1,k)-phiold(i,j-1,k))/4.0) + 1e-15);
                            curvat += (phiold(i,j+1,k) - phiold(i,j,k))/
                                    (std::hypot(phiold(i,j+1,k) - phiold(i,j,k),(phiold(i+1,j+1,k)+phiold(i+1,j,k)-phiold(i-1,j+1,k)-phiold(i-1,j,k))/4.0) + 1e-15);
                            curvat -= (phiold(i,j,k) - phiold(i,j-1,k))/
                                    (std::hypot(phiold(i,j,k) - phiold(i,j-1,k),(phiold(i+1,j-1,k)+phiold(i+1,j,k)-phiold(i-1,j-1,k)-phiold(i-1,j,k))/4.0) + 1e-15);
                            curvat = curvat/dx[0];
                            //if(isAxisymmetric)
                            //{
                            //    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                            //    laplacian_phi += (1.0/y)*(phiold(i , j + 1 , k) - phiold(i , j - 1 , k))/(2.0*dx[1]);
                            //    amrex::Real ny = 0.5*(phiold(i, j + 1, k) - phiold(i, j - 1, k))/(grad_phi*dx[1] + 1e-15 );
                            //    curvat += (1.0/y)*ny;
                            //}
                            //amrex::Print()<<"W_phasefield = "<<W_phasefield*W_phasefield<<'\n';
                            RHS_phi = b*(laplacian_phi + (phiold(i,j,k)*(1.0 - phiold(i,j,k)*phiold(i,j,k))/(W_phasefield * W_phasefield * dx[0]*dx[0])) - grad_phi*curvat);
                            //if(std::isnan(RHS_phi))
                            {
                                //amrex::Print()<<"RHS_phi = "<<RHS_phi<<'\n';
                                //amrex::Print()<<"curvat = "<<curvat<<'\n';
                                //exit(3);
                            }
                        }

                        //if (isAxisymmetric)
                        //{
                            //! box is cell centered
                            //amrex::Real yn = prob_lo[1] + dx[1] * (j + 1);
                            //amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                            //amrex::Real ys = prob_lo[1] + dx[1] * j;
            
                            //fn *= yn;
                            //fs *= ys;
                            //phi(i, j, k) = phi_frk(i, j, k) + dt * ((fw - fe) / dx[0] + (fs - fn) / (y * dx[1]));
                        //}
                        //else
                        //{
                            //phi(i, j, k) = phi_frk(i, j, k) + dt * ((fw - fe) / dx[0] + (fs - fn) / dx[1]);
                        //}
			//Add the sharpening term to RHS
			phi(i, j, k) = phi_frk(i, j, k) - dt * (u_Mid * Phix + v_Mid * Phiy);
			if(sharpPhaseField)phi(i, j, k) += dt*RHS_phi;
			//Decay due to damage
			if(DamageModel)phi(i, j, k) += -1.0*dt*dmg(i, j, k)*(damage_coeff/dx[1])*(1.0 + phiold(i, j, k));
			/*if(DamageModel)
		        {
			    if(dmg(i, j, k) == 1) phi(i, j, k) = -1.0;
			}*/
			//Regularization
			//phi(i, j, k) = std::min(1.0, phi(i, j ,k));
			//phi(i, j, k) = std::max(0.0, phi(i, j ,k));
                    }
                    if(mask(i,j,k) != 1 && mask(i,j,k) != 2)
                    {
                        //auto &solid = interfaces[0];
                        //phi(i, j, k) = 0.0;//P_interface[0] * solid->Volume() / (P_interface0[0] * solid->Volume_0());
                        //amrex::Print()<<"T(i, j, k) = "<<T(i, j, k)<<'\n';
                    }
                });
            }       
        }
    }

void incFSI::AdvectScalars_RK2(int lev)
{
    /// fill internal ghost values
    xvel_old[lev].FillBoundary();
    yvel_old[lev].FillBoundary();

    //! create temporary field data
    const amrex::Box &domain = geom[lev].Domain(); 
    amrex::BoxArray xba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(0));
    amrex::BoxArray yba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(1));
    amrex::MultiFab tmp_xvel_old(xba, dmap[lev], 1, Nghost);
    amrex::MultiFab tmp_yvel_old(yba, dmap[lev], 1, Nghost);
    if (lev > 0)
    {
        //! fill from old current level and coarse level
        //! need to pass the current data also
        amrex::Vector<amrex::Real> ctime(2);
        amrex::Vector<amrex::Real> ftime(2);

        amrex::Vector<amrex::MultiFab *> cxvel(1);
        amrex::Vector<amrex::MultiFab *> cyvel(1);

        //if (lev > 0)
        {
            ctime[0] = t_old;
            ctime[1] = t_new;

            cxvel[0] = &xvel_old[lev - 1];
            cyvel[0] = &yvel_old[lev - 1];
        }
        ftime[0] = t_old;
        ftime[1] = t_new;

        if(divergence_free_interpolation)
        {
            amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> tmp_vel_old = {&tmp_xvel_old, &tmp_yvel_old};
            amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> cvel;
            amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> fvel;

            cvel.resize(1);
            fvel.resize(1);
            cvel[0] = {&xvel_old[lev - 1], &yvel_old[lev - 1]};
            fvel[0] = {&xvel_old[lev], &yvel_old[lev]};

            FillPatch(lev, Time, tmp_vel_old, cvel, {Time}, fvel, {Time}, 0, 1);

            amrex::MultiFab::Copy(tmp_xvel_old, *tmp_vel_old[0], 0, 0, 1, Nghost);
            amrex::MultiFab::Copy(tmp_yvel_old, *tmp_vel_old[1], 0, 0, 1, Nghost);
        }
        else
        {
            FillPatch(lev, Time, tmp_xvel_old, {cxvel[0]}, {ctime[0]}, {&xvel_old[lev]}, {ftime[0]}, 0, 1);
            FillPatch(lev, Time, tmp_yvel_old, {cyvel[0]}, {ctime[0]}, {&yvel_old[lev]}, {ftime[0]}, 0, 1);
        }

        //! need to fill physical bcs
        XVelBoundaryConditions(lev, Time, tmp_xvel_old);
        YVelBoundaryConditions(lev, Time, tmp_yvel_old);
    }
    else
    {
        amrex::MultiFab::Copy(tmp_xvel_old, xvel_old[lev], 0, 0, 1, Nghost);
        amrex::MultiFab::Copy(tmp_yvel_old, yvel_old[lev], 0, 0, 1, Nghost);
    }
    
    const amrex::Real *dx = geom[lev].CellSize();
    const amrex::Real *prob_lo = geom[lev].ProbLo();

    amrex::iMultiFab *UMask;
    amrex::iMultiFab *VMask;
    amrex::iMultiFab *PMask;
    amrex::iMultiFab *Mask;

    amrex::iMultiFab UMask_NI;//NI stands for No Interface, need a better way of doing this
    amrex::iMultiFab VMask_NI;//NI stands for No Interface, need a better way of doing this
    amrex::iMultiFab PMask_NI;//NI stands for No Interface, need a better way of doing this
    amrex::iMultiFab Mask_NI;
        
    if(N_IF > 0 && lev == finest_level)
    {
        auto &&mask__= mask_[lev];
        UMask = &mask__->getUMask();
        VMask = &mask__->getVMask();
        PMask = &mask__->getPMask();
	Mask = &mask__->getMask();
    }
    else
    {
        UMask_NI.define(xba, dmap[lev], 1, Nghost);
        VMask_NI.define(yba, dmap[lev], 1, Nghost);
        PMask_NI.define(grids[lev], dmap[lev], 1, Nghost);
	Mask_NI.define(grids[lev], dmap[lev], 1, Nghost);
  
        UMask_NI.setVal(1);
        VMask_NI.setVal(1);
        PMask_NI.setVal(1);
	Mask_NI.setVal(1);
    }

    CFMask cfmask_rhs_;
    cfmask_rhs_.define(lev, geom, grids, dmap);

    for(int iscalar = 0;iscalar < nscalar;iscalar++)
    {
	//if(iscalar > 1 && iscalar < 8)
	//if(iscalar < 2 || iscalar > 5)
	if(iscalar == eps_max_num)
	    continue;
	if(advect_ref_cond && (iscalar >= F11_num && iscalar <= F33_num))
	    continue;
	if(iscalar == dmg_num)
	    continue;
	//amrex::Print()<<"scalar = "<<iscalar<<'\n';
        Scalars_old[lev][iscalar].FillBoundary();
        //! create temporary field data
        amrex::MultiFab tmp_Scalar_old(grids[lev], dmap[lev], 1, Nghost);
        if (lev > 0)
        {
            //! fill from old current level and coarse level
            //! need to pass the current data also
            amrex::Vector<amrex::MultiFab *> cScalar(1);
            
            amrex::Vector<amrex::Real> ctime(2);
            amrex::Vector<amrex::Real> ftime(2);
	    
            //if (lev > 0)
            {
                cScalar[0] = &Scalars_old[lev - 1][iscalar];
            }
            ftime[0] = t_old;
            ftime[1] = t_new;
	    
            FillPatch(lev, Time, tmp_Scalar_old, cScalar, {ctime[0]}, {&Scalars_old[lev][iscalar]}, {ftime[0]}, 0, 1);
	    
            ScalarBoundaryConditions(lev,iscalar, Time, tmp_Scalar_old);
        }
        else
        {
            amrex::MultiFab::Copy(tmp_Scalar_old, Scalars_old[lev][iscalar], 0, 0, 1, Nghost);
        }
    /// Advect Scalar field
        amrex::Real alpha_ = 1.0;
        for(amrex::MFIter mfi(Scalars[lev][iscalar]); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<amrex::Real const> const &xvel = tmp_xvel_old.const_array(mfi);
            amrex::Array4<amrex::Real const> const &yvel = tmp_yvel_old.const_array(mfi);
            amrex::Array4<amrex::Real const> const &scalarold = tmp_Scalar_old.const_array(mfi);
            amrex::Array4<amrex::Real> const &scalar = Scalars[lev][iscalar].array(mfi);
	    amrex::Array4<amrex::Real> const &frk_scalar = FRK_Scalars[lev][iscalar].array(mfi);
	    amrex::Array4<amrex::Real> const &F11 = Scalars_old[lev][2].array(mfi);
	    amrex::Array4<amrex::Real> const &F12 = Scalars_old[lev][3].array(mfi);
	    amrex::Array4<amrex::Real> const &F21 = Scalars_old[lev][4].array(mfi);
	    amrex::Array4<amrex::Real> const &F22 = Scalars_old[lev][5].array(mfi);
	    amrex::Array4<amrex::Real> const &Vel = U[lev].array(mfi);
	    amrex::Array4<amrex::Real> const &phi = Phi[lev].array(mfi);
        
            amrex::Array4<int> mask, obj_mask;
            if(N_IF == 0 || lev != finest_level)
                mask = PMask_NI.array(mfi);
            else
                mask = PMask->array(mfi);

            if(N_IF == 0 || lev != finest_level)
                obj_mask = Mask_NI.array(mfi);
            else
                obj_mask = Mask->array(mfi);

            amrex::Array4<amrex::Real const> dmg;
            if(DamageModel)
                dmg = Scalars_old[lev][dmg_num].array(mfi);

	    amrex::Array4<int const> const &cfmask = cfmask_rhs_.Mask().const_array(mfi);

            //Testing end
            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
		if(mask(i ,j ,k) == 1 && cfmask(i, j, k) != CFMask::covered)
                {

                    //amrex::Real Fij_x_L, Fij_x_R,Fij_y_L, Fij_y_R, Fij_x, Fij_y;
		    //amrex::Real uFij_x, vFij_y, H;

                    //amrex::Real uR = xvel(i + 1, j, k);
                    //amrex::Real uL = xvel(i, j, k);
                    amrex::Real u_Mid = 0.5 * (xvel(i, j, k) + xvel(i + 1, j, k));

                    //amrex::Real vR = yvel(i, j, k);
                    //amrex::Real vL = yvel(i, j + 1, k);
                    amrex::Real v_Mid = 0.5 * (yvel(i, j, k) + yvel(i, j + 1, k));

                    //WENO5_LS(Fij_x_L, Fij_x_R, Fij_y_L, Fij_y_R, i, j, dx, scalarold);
		    
                    //amrex::Real alpha = std::max(uL,uR);
                    //amrex::Real beta = std::max(vL,vR);

		    //if(obj_mask(i, j, k) == 3)
                    //{
                    //    Fij_x_R = (scalarold(i + 1, j, k) - scalarold(i, j, k)) / dx[0];
                    //    Fij_x_L = (scalarold(i, j, k) - scalarold(i - 1, j, k)) / dx[0];
                    //    Fij_y_R = (scalarold(i, j + 1, k) - scalarold(i, j, k)) / dx[1];
                    //    Fij_y_L = (scalarold(i, j, k) - scalarold(i, j - 1, k)) / dx[1];
                    //}
		    
                    //if(i <= domain.smallEnd(0) + 2 || i >= domain.bigEnd(0) - 2)
                    //{
                    //    Fij_x_R = (scalarold(i + 1, j, k) - scalarold(i, j, k)) / dx[0];
                    //    Fij_x_L = (scalarold(i, j, k) - scalarold(i - 1, j, k)) / dx[0];
                    //}
		    //if((!isAxisymmetric && j <= domain.smallEnd(1) + 2 )|| j >= domain.bigEnd(1) - 2)
                    //{
                    //    Fij_y_R = (scalarold(i, j + 1, k) - scalarold(i, j, k)) / dx[1];
                    //    Fij_y_L = (scalarold(i, j, k) - scalarold(i, j - 1, k)) / dx[1];
                    //}
		    
		    
		    //if(std::abs(u_Mid) > 10.0 || std::abs(v_Mid) > 10.0)
		    //{
                    //    Fij_x_R = (scalarold(i + 1, j, k) - scalarold(i, j, k)) / dx[0];
                    //    Fij_x_L = (scalarold(i, j, k) - scalarold(i - 1, j, k)) / dx[0];
                    //    Fij_y_R = (scalarold(i, j + 1, k) - scalarold(i, j, k)) / dx[1];
                    //    Fij_y_L = (scalarold(i, j, k) - scalarold(i, j - 1, k)) / dx[1];
		    //}
		    
                    //if(uR*uL >= 0.0 && vR*vL >= 0.0)
                    //{
                    //    if (u_Mid > 0.0)
                    //        Fij_x = Fij_x_L;
                    //    else
                    //        Fij_x = Fij_x_R;
                    //    if (v_Mid > 0.0)
                    //        Fij_y = Fij_y_L;
                    //    else
                    //        Fij_y = Fij_y_R;

                    //    H = u_Mid*Fij_x + v_Mid*Fij_y;
                    //}
                    //else if(uR*uL <= 0.0 && vR*vL <= 0.0)
                    //{
                    //    Fij_x = 0.5*(Fij_x_R + Fij_x_L);
                    //    Fij_y = 0.5*(Fij_y_R + Fij_x_L);

                    //    H = u_Mid*Fij_x + v_Mid*Fij_y
                    //        - alpha*0.5*(Fij_x_R - Fij_x_L)
                    //        - beta*0.5*(Fij_y_R - Fij_y_L);
                    //}
                    //else if(uR*uL <= 0.0)
                    //{
                    //    Fij_x = 0.5*(Fij_x_R + Fij_x_L);
                    //    if (v_Mid > 0.0)
                    //        Fij_y = Fij_y_L;
                    //    else
                    //        Fij_y = Fij_y_R;

                    //    H = u_Mid*Fij_x + v_Mid*Fij_y
                    //        - alpha*0.5*(Fij_x_R - Fij_x_L);
                    //}
                    //else if(vR*vL <= 0.0)
                    //{
                    //    if (u_Mid > 0.0)
                    //        Fij_x = Fij_x_L;
                    //    else
                    //        Fij_x = Fij_x_R;
                    //    Fij_y = 0.5*(Fij_y_R + Fij_x_L);

                    //    H = u_Mid*Fij_x + v_Mid*Fij_y
                    //        - beta*0.5*(Fij_y_R - Fij_y_L);
                    //}


                    amrex::Real Fij_e, Fij_w, Fij_n, Fij_s;
                    amrex::Real Fij_e_R,Fij_e_L,Fij_w_R,Fij_w_L;
                    amrex::Real Fij_n_R,Fij_n_L,Fij_s_R,Fij_s_L;
                    amrex::Real junk1,junk2;
		    if(Fij_order == 5)
	            {
			//amrex::Print()<<"Fifth order WENO for FIJ"<<"\n";
                        WENO5_Fij(Fij_e_L, Fij_e_R, Fij_n_L, Fij_n_R, i, j, dx, scalarold);
                        WENO5_Fij(Fij_w_L, Fij_w_R, junk1, junk2, i - 1, j, dx, scalarold);
                        WENO5_Fij(junk1, junk2, Fij_s_L, Fij_s_R, i, j - 1, dx, scalarold);
                    }
		    
                    amrex::Real scalar_e, scalar_w,scalar_n,scalar_s; 
		    amrex::Real fe, fw, fn, fs;
		    bool drop_order = false;
		    if(obj_mask(i, j, k) == 3 && phi(i , j , k ) <= -0.5)
		    {
		        drop_order = true;//drop order near interface if interface is in water/damaged region 
		    }
                    if(drop_order)
		    {
			//amrex::Print()<<"in 2nd order"<<'\n';
                        scalar_e = 0.5 * (scalarold(i, j, k) + scalarold(i + 1, j, k));
                        fe = scalar_e * xvel(i + 1, j, k) 
                                     + alpha_ * 0.5 * std::abs(xvel(i + 1, j, k)) * (scalarold(i, j, k) - scalarold(i + 1, j, k));
        
                        scalar_w = 0.5 * (scalarold(i, j, k) + scalarold(i - 1, j, k));
                        fw = scalar_w * xvel(i, j, k) 
                                     + alpha_ * 0.5 * std::abs(xvel(i, j, k)) * (scalarold(i - 1, j, k) - scalarold(i, j, k));
                    
                        scalar_n = 0.5 * (scalarold(i, j, k) + scalarold(i, j + 1, k));
                        fn = scalar_n * yvel(i, j + 1, k) 
                                     + alpha_ * 0.5 * std::abs(yvel(i, j + 1, k)) * (scalarold(i, j, k) - scalarold(i, j + 1, k));
        
                        scalar_s = 0.5 * (scalarold(i, j, k) + scalarold(i, j - 1, k));
                        fs = scalar_s * yvel(i, j, k) 
                                     + alpha_ * 0.5 * std::abs(yvel(i, j, k)) * (scalarold(i, j - 1, k) - scalarold(i, j, k));
		    }
		    else 
                    {
                        scalar_e = 0.5 * (Fij_e_L + Fij_e_R);
                        fe = scalar_e * xvel(i + 1, j, k)
                                     + alpha_ * 0.5 * std::abs(xvel(i + 1, j, k)) * (Fij_e_L - Fij_e_R);


                        scalar_w = 0.5 * (Fij_w_L + Fij_w_R);
                        fw = scalar_w * xvel(i, j, k)
                                     + alpha_ * 0.5 * std::abs(xvel(i, j, k)) * (Fij_w_L - Fij_w_R);


                        scalar_n = 0.5 * (Fij_n_L + Fij_n_R);
                        fn = scalar_n * yvel(i, j + 1, k)
                                     + alpha_ * 0.5 * std::abs(yvel(i, j + 1, k)) * (Fij_n_L - Fij_n_R);


                        scalar_s = 0.5 * (Fij_s_L + Fij_s_R);
                        fs = scalar_s * yvel(i, j, k)
                                     + alpha_ * 0.5 * std::abs(yvel(i, j, k)) * (Fij_s_L - Fij_s_R);
                    }


        
		    //Compute source for F_ij
                    amrex::Real source_Fij = 0.0;
                    amrex::Real ux = (xvel(i + 1, j, k) - xvel(i, j, k)) / dx[0];
                    amrex::Real uy = (xvel(i + 1, j + 1, k) - xvel(i + 1, j - 1, k) + xvel(i, j + 1, k) - xvel(i, j - 1, k)) / (4.0 * dx[1]);
                    amrex::Real vx = (yvel(i + 1, j + 1, k) - yvel(i - 1, j + 1, k) + yvel(i + 1, j, k) - yvel(i - 1, j, k)) / (4.0 * dx[0]);
                    amrex::Real vy = (yvel(i, j + 1, k) - yvel(i, j, k)) / dx[1];

		    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);

		    //if(phi(i, j, k) >=0)
		    {
	                if(iscalar == 2)
                            source_Fij = ux*F11(i, j, k) + uy*F21(i, j, k);  
	                else if(iscalar == 3)
                            source_Fij = ux*F12(i, j, k) + uy*F22(i, j, k);
                        else if(iscalar == 4)
                            source_Fij = vx*F11(i, j, k) + vy*F21(i, j, k);
                        else if(iscalar == 5)
                            source_Fij = vx*F12(i, j, k) + vy*F22(i, j, k);
	                else if(iscalar == 6)
                            source_Fij = v_Mid * scalarold(i, j, k)/y;
		    }
		    if(PhaseField)
		    {
		        source_Fij *= 0.5*(1.0 + phi(i, j, k));
		    }

                    if (isAxisymmetric)
                    {
                        //! box is cell centered
                        amrex::Real yn = prob_lo[1] + dx[1] * (j + 1);
                        amrex::Real ys = prob_lo[1] + dx[1] * j;
        
                        fn *= yn;
                        fs *= ys;
                        scalar(i, j, k) = frk_scalar(i, j, k) + dt * ((fw - fe) / dx[0] + (fs - fn) / (y * dx[1]) + source_Fij);
			//scalar(i, j, k) = frk_scalar(i, j, k) - dt * (u_Mid * Fij_x + v_Mid * Fij_y)  + dt * source_Fij;
			//if(j < 2)scalar(i, j, k) = frk_scalar(i, j, k) - dt * H  + dt * source_Fij;
                    }
                    else
                    {
                        scalar(i, j, k) = frk_scalar(i, j, k) +  dt * ((fw - fe) / dx[0] + (fs - fn) / dx[1] + source_Fij);
			//scalar(i, j, k) = frk_scalar(i, j, k) - dt * (u_Mid * Fij_x + v_Mid * Fij_y)  + dt * source_Fij;
                    }

                    /*if(DamageModel)
		    {
			amrex:: Real A = 0;
			if(iscalar == F11_num)
			   A = 1;
                        else if(iscalar == F12_num)
                           A = 0;
			else if(iscalar == F21_num)
                           A = 0;
			else if(iscalar == F22_num)
                           A = 1;
			else if(iscalar == F33_num)
                           A = 1;
		        scalar(i, j, k) += -1.0*dt*dmg(i, j, k)*(damage_coeff/dx[1])*(scalarold(i, j, k) - A);        
		    }*/
                }
            });
        }       
    }
}
} // namespace mycode
