#include "incFSI.H"
#include <WeightedENO.h>

namespace mycode
{

void incFSI::ComputeIntermediateVelocity(int lev)
{
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
            if(TempField)cTheta[0] = &Theta[lev - 1];
            if(PhaseField)cPhi[0] = &Phi[lev - 1];

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

    amrex::iMultiFab UMask_NI;//NI stands for No Interface, need a better way of doing this
    amrex::iMultiFab VMask_NI;//NI stands for No Interface, need a better way of doing this
    amrex::iMultiFab PMask_NI;//NI stands for No Interface, need a better way of doing this
        
    if(N_IF > 0)
    {
        auto &&mask__= mask_[lev];
        UMask = &mask__->getUMask();
        VMask = &mask__->getVMask();
        PMask = &mask__->getPMask();
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


    /// compute x-component
    for(amrex::MFIter mfi(xvel[lev]); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        amrex::Array4<amrex::Real const> const &u = tmp_xvel_old.const_array(mfi);
        amrex::Array4<amrex::Real const> const &v = tmp_yvel_old.const_array(mfi);
        amrex::Array4<amrex::Real const> const &P = tmp_Pressure.const_array(mfi);
        amrex::Array4<amrex::Real> const &ustar = xvel[lev].array(mfi);
        amrex::Box box = get_valid_face_box(lev, bx, 0);

        amrex::Array4<int> umask;

        if(N_IF == 0)
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
            if(umask(i, j, k) == 1)
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

                //Mu_e = Mu;

                amrex::Real fe = ue * ue - 2.0 * Mu_e * ux_e;
                /*amrex::Real fe = ue * ue 
                                 + alpha_ * 0.5 * abs(ue) * (u(i, j, k) - u(i + 1, j, k)) 
                                 - 2.0 * Mu * ux_e;*/   

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

                //Mu_w = Mu;

                amrex::Real fw = uw * uw - 2.0 * Mu_w * ux_w;
                //amrex::Real fw = uw * uw 
                //                 + alpha_ * 0.5 * abs(uw) * (u(i - 1, j, k) - u(i, j, k))
                //                 - 2.0 * Mu * ux_w;


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

                //Mu_n = Mu;
                amrex::Real fn = un * vn - Mu_n * (uy_n + vx_n);
                //amrex::Real fn = un * vn 
                //                 + alpha_ * 0.5 * abs(vn) * (u(i, j, k) - u(i, j + 1, k))
                //                 - Mu * (uy_n + vx_n);

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

                //Mu_s = Mu;
                amrex::Real fs = us * vs - Mu_s * (uy_s + vx_s);
                //amrex::Real fs =  us * vs 
                //                  + alpha_ * 0.5 * abs(vs) * (u(i, j - 1, k) - u(i, j , k))
                //                 - Mu * (uy_s + vx_s);
                
                //amrex::Print()<<"Mu_e = "<<Mu_e<<" , Mu_w = "<<Mu_w<<" , Mu_n = "<<Mu_n<<" , Mu_s = "<<Mu_s<<'\n';
                Mu_max = std::max(0.25*(Mu_e+Mu_w+Mu_n+Mu_s), Mu_max);
                
                if (isAxisymmetric)
                {
                    //TODO how to avoid singularity? i.e. what if y = 0
                    //! box is centered in y
                    amrex::Real yn = prob_lo[1] + dx[1] * (j + 1);
                    amrex::Real y  = prob_lo[1] + dx[1] * (j + 0.5);
                    amrex::Real ys = prob_lo[1] + dx[1] * j;
                
                    fn *= yn;
                    fs *= ys;
                    ustar(i, j, k) += dt * ( (fw - fe) / dx[0] + (fs - fn) / (y * dx[1]) );                    
                }
                else
                {
                    ustar(i, j, k) += dt * ( (fw - fe) / dx[0] + (fs - fn) / dx[1] );
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
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const &u = tmp_xvel_old.const_array(mfi);
        amrex::Array4<amrex::Real const> const &v = tmp_yvel_old.const_array(mfi);
        amrex::Array4<amrex::Real const> const &P = tmp_Pressure.const_array(mfi);
        amrex::Array4<amrex::Real> const &vstar = yvel[lev].array(mfi);
        amrex::Box box = get_valid_face_box(lev, bx, 1);

        amrex::Array4<int> vmask;
        if(N_IF == 0)
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
            if(vmask(i ,j ,k) == 1)
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

                amrex::Real Mu_e = visc_.GetViscosity(S11_e, S12_e, S22_e, phi_e, dmg_e);
                
                //Mu_e = Mu;
                amrex::Real fe = ue * ve - Mu_e * (uy_e + vx_e);
                //amrex::Real fe = ue * ve 
                //                 + alpha_ * 0.5 * abs(ue) * (v(i, j, k) - v(i + 1, j, k)) 
                //                 - Mu * (uy_e + vx_e);

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


                amrex::Real Mu_w = visc_.GetViscosity(S11_w, S12_w, S22_w, phi_w, dmg_w);

                //Mu_w = Mu;
                amrex::Real fw = uw * vw - Mu_w * (uy_w + vx_w);
                //amrex::Real fw = uw * vw 
                //                 + alpha_ * 0.5 * abs(uw) * (v(i - 1, j, k) - v(i, j, k))
                //                 - Mu * (uy_w + vx_w);

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

                amrex::Real Mu_n = visc_.GetViscosity(S11_n, S12_n, S22_n, phi_n, dmg_n);

                //Mu_n = Mu;
                amrex::Real fn = vn * vn - 2.0 * Mu_n * vy_n;
                //amrex::Real fn = vn * vn
                //                  + alpha_ * 0.5 * abs(vn) * (v(i, j, k) - v(i, j + 1, k))
                //                  - 2.0 * Mu * vy_n;

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

                amrex::Real Mu_s = visc_.GetViscosity(S11_s, S12_s, S22_s, phi_s, dmg_s);
                //Mu_s = Mu;
                amrex::Real fs = vs * vs - 2.0 * Mu_s * vy_s;
                
		Mu_max = std::max(0.25*(Mu_e+Mu_w+Mu_n+Mu_s), Mu_max);
                if (isAxisymmetric)
                {
                    //! box is nodal in y
                    amrex::Real yn = prob_lo[1] + dx[1] * (j + 0.5);
                    amrex::Real y = prob_lo[1] + dx[1] * j;
                    amrex::Real ys = prob_lo[1] + dx[1] * (j - 0.5);
                
                    fn *= yn;
                    fs *= ys;
                
                    //! avoid singularity
                    if (fabs(y) <  1e-10)
                    {
                        vstar(i, j, k) = 0.0;
                    }
                    else
                    {
                        vstar(i, j, k) += dt * ((fw - fe) / dx[0] + (fs - fn) / (y * dx[1]) - 2.0 * Mu * v(i, j, k) / (y * y));
                    }
                }
                else
                {
                    vstar(i, j, k) += dt * ((fw - fe) / dx[0] + (fs - fn) / dx[1]);
                }
                
                if (Projection)
                {
                    amrex::Real dpdy = (P(i, j, k) - P(i, j - 1, k)) / dx[1];
                    vstar(i, j, k) -= dt * dpdy;
                }
                /*if(i == 350 && (j == 512 || j == 511))
                {
                    amrex::Print()<<" i = "<<i<<" , j = "<<j<<'\n';
                    amrex::Print()<<"Predicct vel"<<'\n';
                    amrex::Print()<<"vstar("<<i<<", "<<j<<", "<<k<<") = "<<vstar(i, j, k)<<'\n';
                }*/
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
            amrex::Array4<amrex::Real> const &T = Theta[lev].array(mfi);
    
            amrex::Array4<int> mask;
            if(N_IF == 0)
                mask = PMask_NI.array(mfi);
            else
                mask = PMask->array(mfi);
    
    
            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                if(mask(i ,j ,k) == 1)
                {
                    amrex::Real theta_e = 0.5 * (Told(i, j, k) + Told(i + 1, j, k));
                    amrex::Real thetax_e = (Told(i + 1, j, k) - Told(i, j, k)) / dx[0];
                    k_ = visc_.GetConductivity(Mu, Prandtl_no);
                    if (mask(i,j,k) == 1 && mask(i+1,j,k) == 1)
                    {
                        k_max = std::max(k_, k_max);
                    }
                    amrex::Real fe = theta_e * xvel(i + 1, j, k) - k_ * thetax_e;
                    //amrex::Real fe = theta_e * xvel(i + 1, j, k) 
                    //                 + alpha_ * 0.5 * abs(xvel(i + 1, j, k)) * (Told(i, j, k) - Told(i + 1, j, k))
                    //                 - k_ * thetax_e;
    
                    amrex::Real theta_w = 0.5 * (Told(i, j, k) + Told(i - 1, j, k));
                    amrex::Real thetax_w = (Told(i, j, k) - Told(i - 1, j, k)) / dx[0];
                    k_ = visc_.GetConductivity(Mu, Prandtl_no);
                    if (mask(i, j, k) == 1 && mask(i - 1, j, k) == 1)
                    {
                        k_max = std::max(k_, k_max);
                    }
                    amrex::Real fw = theta_w * xvel(i, j, k) - k_ * thetax_w;
                    //amrex::Real fw = theta_w * xvel(i, j, k) 
                    //                 + alpha_ * 0.5 * abs(xvel(i, j, k)) * (Told(i - 1, j, k) - Told(i, j, k))
                    //                 - k_ * thetax_w;
                    
                    amrex::Real theta_n = 0.5 * (Told(i, j, k) + Told(i, j + 1, k));
                    amrex::Real thetay_n = (Told(i, j + 1, k) - Told(i, j, k)) / dx[1];
                    k_ = visc_.GetConductivity(Mu, Prandtl_no);
                    if (mask(i, j, k) == 1 && mask(i, j + 1, k) == 1)
                    {
                        k_max = std::max(k_, k_max);
                    }
                    amrex::Real fn = theta_n * yvel(i, j + 1, k) - k_ * thetay_n;
                    //amrex::Real fn = theta_n * yvel(i, j + 1, k) 
                    //                 + alpha_ * 0.5 * abs(yvel(i, j + 1, k)) * (Told(i, j, k) - Told(i, j + 1, k))
                    //                 - k_ * thetay_n;
    
                    amrex::Real theta_s = 0.5 * (Told(i, j, k) + Told(i, j - 1, k));
                    amrex::Real thetay_s = (Told(i, j, k) - Told(i, j - 1, k)) / dx[1];
                    k_ = visc_.GetConductivity(Mu, Prandtl_no);
                    if (mask(i, j, k) == 1 && mask(i, j - 1, k) == 1)
                    {
                        k_max = std::max(k_, k_max);
                    }
                    amrex::Real fs = theta_s * yvel(i, j, k) - k_ * thetay_s;
                    //amrex::Real fs = theta_s * yvel(i, j, k) 
                    //                 + alpha_ * 0.5 * abs(yvel(i, j, k)) * (Told(i, j - 1, k) - Told(i, j, k))
                    //                 - k_ * thetay_s;
    
                    amrex::Real ux = (xvel(i + 1, j, k) - xvel(i, j, k)) / dx[0];
                    amrex::Real uy = (xvel(i + 1, j + 1, k) - xvel(i + 1, j - 1, k) + xvel(i, j + 1, k) - xvel(i, j - 1, k)) / (4.0 * dx[1]);
                    amrex::Real vx = (yvel(i + 1, j + 1, k) - yvel(i - 1, j + 1, k) + yvel(i + 1, j, k) - yvel(i - 1, j, k)) / (4.0 * dx[0]);
                    amrex::Real vy = (yvel(i, j + 1, k) - yvel(i, j, k)) / dx[1];
                    amrex::Real dGamma_dt = visc_.ComputeGammaDot(ux, 0.5 * (uy + vx), vy);
                    amrex::Real Mu_ = visc_.GetViscosity(dGamma_dt);
                    //Mu_ = Mu;
                    if (isAxisymmetric)
                    {
                        //! box is cell centered
                        amrex::Real yn = prob_lo[1] + dx[1] * (j + 1);
                        amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                        amrex::Real ys = prob_lo[1] + dx[1] * j;
    
                        fn *= yn;
                        fs *= ys;
                        T(i, j, k) += dt * ((fw - fe) / dx[0] + (fs - fn) / (y * dx[1]) +  Mu_ * (u_ref*u_ref/cp_/T_ref) * dGamma_dt* dGamma_dt);
                    }
                    else
                    {
                        T(i, j, k) +=  dt * ((fw - fe) / dx[0] + (fs - fn) / dx[1] +  Mu_ * (u_ref*u_ref/cp_/T_ref) * dGamma_dt * dGamma_dt);
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
        for(amrex::MFIter mfi(Theta[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<amrex::Real const> const &xvel = tmp_xvel_old.const_array(mfi);
            amrex::Array4<amrex::Real const> const &yvel = tmp_yvel_old.const_array(mfi);
            amrex::Array4<amrex::Real const> const &phiold = tmp_Phi_old.const_array(mfi);
            amrex::Array4<amrex::Real> const &phi = Phi[lev].array(mfi);
        
            amrex::Array4<int> mask;
            if(N_IF == 0)
                mask = PMask_NI.array(mfi);
            else
                mask = PMask->array(mfi);
        
        
            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                if(mask(i ,j ,k) == 1)
                {
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
                    if (isAxisymmetric)
                    {
                        //! box is cell centered
                        amrex::Real yn = prob_lo[1] + dx[1] * (j + 1);
                        amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                        amrex::Real ys = prob_lo[1] + dx[1] * j;
        
                        fn *= yn;
                        fs *= ys;
                        phi(i, j, k) = phiold(i, j, k) + dt * ((fw - fe) / dx[0] + (fs - fn) / (y * dx[1]));
                    }
                    else
                    {
                        phi(i, j, k) = phiold(i, j, k) +  dt * ((fw - fe) / dx[0] + (fs - fn) / dx[1]);
                    }
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

void incFSI::AdvectScalars(int lev)
{
    /// fill internal ghost values
    xvel_old[lev].FillBoundary();
    yvel_old[lev].FillBoundary();

    //! create temporary field data
    
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

    amrex::iMultiFab UMask_NI;//NI stands for No Interface, need a better way of doing this
    amrex::iMultiFab VMask_NI;//NI stands for No Interface, need a better way of doing this
    amrex::iMultiFab PMask_NI;//NI stands for No Interface, need a better way of doing this
        
    if(N_IF > 0)
    {
        auto &&mask__= mask_[lev];
        UMask = &mask__->getUMask();
        VMask = &mask__->getVMask();
        PMask = &mask__->getPMask();
    }
    else
    {
        UMask_NI.define(xba, dmap[lev], 1, Nghost);
        VMask_NI.define(yba, dmap[lev], 1, Nghost);
        PMask_NI.define(grids[lev], dmap[lev], 1, Nghost);
  
        UMask_NI.setVal(1);
        VMask_NI.setVal(1);
        PMask_NI.setVal(1);
    }
    for(int iscalar = 0;iscalar < nscalar;iscalar++)
    {
	//if(iscalar > 1 && iscalar < 8)
	//if(iscalar < 2 || iscalar > 5)
	if(iscalar == eps_max_num)
	    continue;
	if(advect_ref_cond && (iscalar >= F11_num && iscalar <= F33_num))
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
	    amrex::Array4<amrex::Real> const &F11 = Scalars_old[lev][2].array(mfi);
	    amrex::Array4<amrex::Real> const &F12 = Scalars_old[lev][3].array(mfi);
	    amrex::Array4<amrex::Real> const &F21 = Scalars_old[lev][4].array(mfi);
	    amrex::Array4<amrex::Real> const &F22 = Scalars_old[lev][5].array(mfi);
        
            amrex::Array4<int> mask;
            if(N_IF == 0)
                mask = PMask_NI.array(mfi);
            else
                mask = PMask->array(mfi);
        
        
            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                //if((mask(i ,j ,k) > 0 && iscalar < 2) || (mask(i ,j ,k) == 1 && iscalar == 8))
		if(mask(i ,j ,k) == 1)
                {
                    amrex::Real scalar_e = 0.5 * (scalarold(i, j, k) + scalarold(i + 1, j, k));
                    //amrex::Real fe = scalarold(i, j, k) * xvel(i + 1, j, k); 
                    //amrex::Real fe = scalar_e * xvel(i + 1, j, k);
                    amrex::Real fe = scalar_e * xvel(i + 1, j, k) 
                                     + alpha_ * 0.5 * std::abs(xvel(i + 1, j, k)) * (scalarold(i, j, k) - scalarold(i + 1, j, k));
        
        
                    amrex::Real scalar_w = 0.5 * (scalarold(i, j, k) + scalarold(i - 1, j, k));
                    //amrex::Real fw = scalar_w * xvel(i, j, k);
                    amrex::Real fw = scalar_w * xvel(i, j, k) 
                                     + alpha_ * 0.5 * std::abs(xvel(i, j, k)) * (scalarold(i - 1, j, k) - scalarold(i, j, k));
        
                    
                    amrex::Real scalar_n = 0.5 * (scalarold(i, j, k) + scalarold(i, j + 1, k));
                    //amrex::Real fn = scalar_n * yvel(i, j + 1, k);
                    amrex::Real fn = scalar_n * yvel(i, j + 1, k) 
                                     + alpha_ * 0.5 * std::abs(yvel(i, j + 1, k)) * (scalarold(i, j, k) - scalarold(i, j + 1, k));
        
        
                    amrex::Real scalar_s = 0.5 * (scalarold(i, j, k) + scalarold(i, j - 1, k));
                    //amrex::Real fs = scalar_s * yvel(i, j, k);
                    amrex::Real fs = scalar_s * yvel(i, j, k) 
                                     + alpha_ * 0.5 * std::abs(yvel(i, j, k)) * (scalarold(i, j - 1, k) - scalarold(i, j, k));
        
		    //Compute source for F_ij
                    amrex::Real source_Fij = 0.0;
		    {
                        amrex::Real ux = (xvel(i + 1, j, k) - xvel(i, j, k)) / dx[0];
                        amrex::Real uy = (xvel(i + 1, j + 1, k) - xvel(i + 1, j - 1, k) + xvel(i, j + 1, k) - xvel(i, j - 1, k)) / (4.0 * dx[1]);
                        amrex::Real vx = (yvel(i + 1, j + 1, k) - yvel(i - 1, j + 1, k) + yvel(i + 1, j, k) - yvel(i - 1, j, k)) / (4.0 * dx[0]);
                        amrex::Real vy = (yvel(i, j + 1, k) - yvel(i, j, k)) / dx[1];

	                if(iscalar == 2)
			{
			    source_Fij = ux*F11(i, j, k) + uy*F21(i, j, k);
			}
			else if(iscalar == 3)
                        {
                            source_Fij = ux*F12(i, j, k) + uy*F22(i, j, k);
                        }
                        else if(iscalar == 4)
                        {
                            source_Fij = vx*F11(i, j, k) + vy*F21(i, j, k);
                        }
                        else if(iscalar == 5)
                        {
                            source_Fij = vx*F12(i, j, k) + vy*F22(i, j, k);
                        }
			else if(iscalar == 6)
                        {
                            source_Fij = 0.0;
                        }
		    }
                    if (isAxisymmetric && iscalar != F33_num)
                    {
                        //! box is cell centered
                        amrex::Real yn = prob_lo[1] + dx[1] * (j + 1);
                        amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                        amrex::Real ys = prob_lo[1] + dx[1] * j;
        
                        fn *= yn;
                        fs *= ys;
                        scalar(i, j, k) = scalarold(i, j, k) + dt * ((fw - fe) / dx[0] + (fs - fn) / (y * dx[1]) + source_Fij);
                    }
                    else
                    {
                        scalar(i, j, k) = scalarold(i, j, k) +  dt * ((fw - fe) / dx[0] + (fs - fn) / dx[1] + source_Fij);
                    }
                }
		if(i == 0 && j == 0)
		{
		    //amrex::Print()<<"i = "<<i<<" , j = "<<j<<'\n';
		    //amrex::Print()<<"scalar(i, j, k) = "<<scalar(i, j, k)<<'\n';
		    //amrex::Print()<<"scalarold(i, j, k) = "<<scalarold(i, j, k)<<'\n';
		}
		//else if( (mask(i ,j ,k) != 1 && iscalar == 8) )
		//else
		//{
		//    if(mask(i-1,j,k) == 1 )
	        //        scalar(i, j, k) = scalar(i - 1, j, k);
		//    else if(mask(i + 1,j,k) == 1 )
                //        scalar(i, j, k) = scalar(i + 1, j, k);
		//    else if(mask(i,j - 1,k) == 1 )
                //        scalar(i, j, k) = scalar(i, j - 1, k);
		//    else if(mask(i,j + 1,k) == 1 )
                //        scalar(i, j, k) = scalar(i, j + 1, k);
		//}
                //else if((mask(i ,j ,k) == 0 && iscalar < 2) )
                //{
		//    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
		//    amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                //    if(mask(i-1,j,k) > 0 )
                //        scalar(i, j, k) = scalar(i - 1, j, k);
                //    else if(mask(i + 1,j,k) > 0 )
                //        scalar(i, j, k) = scalar(i + 1, j, k);
                //    else if(mask(i,j - 1,k) > 0  )
                //        scalar(i, j, k) = scalar(i, j - 1, k);
                //    else if(mask(i,j + 1,k) > 0 )
                //        scalar(i, j, k) = scalar(i, j + 1, k);

		    /*if(iscalar == 0)
	                scalar(i, j, k) = x ;
		    else if(iscalar == 1)
			scalar(i, j, k) = y ;
		    else
			scalar(i, j, k) = 0.0 ;
		    */
                    //auto &solid = interfaces[0];
                    //phi(i, j, k) = 0.0;//P_interface[0] * solid->Volume() / (P_interface0[0] * solid->Volume_0());
                    //amrex::Print()<<"T(i, j, k) = "<<T(i, j, k)<<'\n';
                //}
            });
        }       
    }
}

void incFSI::ComputeCollocatedVelocityField(int lev)
{
    U[lev].setVal(0.0);
    const amrex::Box &domain = geom[lev].Domain();
    amrex::BoxArray xba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(0));
    amrex::BoxArray yba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(1));
    amrex::iMultiFab *UMask;
    amrex::iMultiFab *VMask;
    amrex::iMultiFab *PMask;

    amrex::iMultiFab UMask_NI;//NI stands for No Interface, need a better way of doing this
    amrex::iMultiFab VMask_NI;//NI stands for No Interface, need a better way of doing this
    amrex::iMultiFab PMask_NI;//NI stands for No Interface, need a better way of doing this

    if(N_IF > 0)
    {
        auto &&mask__= mask_[lev];
        UMask = &mask__->getUMask();
        VMask = &mask__->getVMask();
        PMask = &mask__->getPMask();
    }
    else
    {
        UMask_NI.define(xba, dmap[lev], 1, Nghost);
        VMask_NI.define(yba, dmap[lev], 1, Nghost);
        PMask_NI.define(grids[lev], dmap[lev], 1, Nghost);

        UMask_NI.setVal(1);
        VMask_NI.setVal(1);
        PMask_NI.setVal(1);
    }


    for(amrex::MFIter mfi(U[lev]); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const &u = xvel[lev].const_array(mfi);
        amrex::Array4<amrex::Real const> const &v = yvel[lev].const_array(mfi);
        amrex::Array4<amrex::Real> const &vel = U[lev].array(mfi);

        amrex::Array4<int const> umask,vmask,pmask;
        if(N_IF == 0)
        {
            umask = UMask_NI.const_array(mfi);
            vmask = VMask_NI.const_array(mfi);
            pmask = PMask_NI.const_array(mfi);
        }
        else        
        {
            umask = UMask->array(mfi);
            vmask = VMask->array(mfi);
            pmask = PMask->array(mfi);
        }

        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (pmask(i,j,k) == 1)
            {
                /// x component
                if (umask(i,j,k) == 1 && umask(i+1,j,k) == 1)
                {
                    vel(i, j, k, 0) = 0.5 * (u(i, j, k) + u(i + 1, j, k));
                }
                else if (umask(i,j,k) == 1 && umask(i+1,j,k) == 0)
                {
                    if (umask(i-1,j,k) == 1 && i > domain.smallEnd(0))
                    {
                        vel(i, j, k, 0) = 1.5 * (u(i, j, k) - 0.5 * u(i - 1, j, k));
                    }
                }
                else if (umask(i,j,k) == 0 && umask(i+1,j,k) == 1)
                {
                    if (umask(i+2,j,k) == 1 && i < domain.bigEnd(0))
                    {
                        vel(i, j, k, 0) = 1.5 * (u(i + 1, j, k) - 0.5 * u(i + 2, j, k));
                    }
                }
                else
                {
                    amrex::Abort("Problem computing collocated x component of velocity");
                }
        
                /// y compoonent
                if (vmask(i,j,k) == 1 && vmask(i,j+1,k) == 1)
                {
                    vel(i, j, k, 1) = 0.5 * (v(i, j, k) + v(i, j + 1, k));
                }
                else if (vmask(i,j,k) == 1 && vmask(i,j+1,k) == 0)
                {
                    if (vmask(i,j-1,k) == 1 && j > domain.smallEnd(1))
                    {
                        vel(i, j, k, 1) = 1.5 * (v(i, j, k) - 0.5 * v(i, j - 1, k));
                    }
                }
                else if (vmask(i,j,k) == 0 && vmask(i,j+1,k) == 1)
                {
                    if (vmask(i,j+2,k) == 1 && j < domain.bigEnd(1))
                    {
                        vel(i, j, k, 1) = 1.5 * (v(i, j + 1, k) - 0.5 * v(i, j + 2, k));
                    }
                }
                else
                {
                    amrex::Abort("Problem computing collocated y component of velocity");
                }
            }
        });
            
    }
    //if(lev == 0)
    U[lev].FillBoundary();
}

void incFSI::ComputeStaggeredVelocityField(int lev)
{
    const amrex::Real *dx = geom[lev].CellSize();
    const amrex::Real *prob_lo = geom[lev].ProbLo();

    for(amrex::MFIter mfi(U[lev]); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const &vel = U[lev].const_array(mfi);
        amrex::Array4<amrex::Real> const &u = xvel[lev].array(mfi);
        amrex::Array4<amrex::Real> const &v = yvel[lev].array(mfi);
        amrex::Array4<int const> const &umask = mask_[lev]->getUMask().const_array(mfi);
        amrex::Array4<int const> const &vmask = mask_[lev]->getVMask().const_array(mfi);        

        /// grown cc box by 1 cell at hi side
        amrex::Box gbx(bx);
        gbx.growHi(0, 1);
        gbx.growHi(1, 1);

        amrex::ParallelFor(gbx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            amrex::Real y = prob_lo[1] + dx[1] * (j);
            if (umask(i,j,k) != 1)
                u(i, j, k) = 0.5 * (vel(i, j, k, 0) + vel(i - 1, j, k, 0));
            if (vmask(i,j,k) != 1)
                v(i, j, k) = 0.5 * (vel(i, j, k, 1) + vel(i, j - 1, k, 1));        
            if(isAxisymmetric && y==0.0 )
                v(i, j, k) = 0.0;

        });
    }
}

} // namespace mycode
