#include "incFSI.H"
#include <WeightedENO.h>
#include <AdvectLS_helper.H>

namespace mycode
{

    void incFSI::ComputeIntermediateVelocity_RK(int lev, int RKStage)
    {
        /// fill internal ghost values
        //xvel_old[lev].FillBoundary();
        //yvel_old[lev].FillBoundary();
        //Pressure[lev].FillBoundary();
        //if(TempField)Theta_old[lev].FillBoundary();
        //if(PhaseField)Phi_old[lev].FillBoundary();
        
        const amrex::Real *dx = geom[lev].CellSize();
        const amrex::Real *prob_lo = geom[lev].ProbLo();
        amrex::BoxArray xba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(0));
        amrex::BoxArray yba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(1));
        /*
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
        
                amrex::MultiFab::Copy(xvel[lev], *tmp_vel_old[0], 0, 0, 1, Nghost);
                amrex::MultiFab::Copy(yvel[lev], *tmp_vel_old[1], 0, 0, 1, Nghost);
            }
            else
            {
                FillPatch(lev, Time, xvel[lev], {cxvel[0]}, {ctime[0]}, {&xvel_old[lev]}, {ftime[0]}, 0, 1);
                FillPatch(lev, Time, yvel[lev], {cyvel[0]}, {ctime[0]}, {&yvel_old[lev]}, {ftime[0]}, 0, 1);
            }
        
            //! no need of time interpolation of following
            FillPatch(lev, Time, Pressure[lev], cP, {ctime[0]}, {&Pressure[lev]}, {ftime[0]}, 0, 1);
            if(TempField) FillPatch(lev, Time, Theta[lev], cTheta, {ctime[0]}, {&Theta_old[lev]}, {ftime[0]}, 0, 1);
            if(PhaseField) FillPatch(lev, Time, Phi[lev], cPhi, {ctime[0]}, {&Phi_old[lev]}, {ftime[0]}, 0, 1);
        
        
            //! need to fill physical bcs
            XVelBoundaryConditions(lev, Time, xvel[lev]);
            YVelBoundaryConditions(lev, Time, yvel[lev]);
            PressureBoundaryConditions(lev, Time, Pressure[lev]);
            if(TempField) TemperatureBoundaryConditions(lev, Time, Theta[lev]);
            if(PhaseField) PhiBoundaryConditions(lev, Time, Phi[lev]);
        }
        else
        {
            amrex::MultiFab::Copy(xvel[lev], xvel_old[lev], 0, 0, 1, Nghost);
            amrex::MultiFab::Copy(yvel[lev], yvel_old[lev], 0, 0, 1, Nghost);
            if(TempField) amrex::MultiFab::Copy(Theta[lev], Theta_old[lev], 0, 0, 1, Nghost);
            if(PhaseField) amrex::MultiFab::Copy(Phi[lev], Phi_old[lev], 0, 0, 1, Nghost);            
        }        
        */
        amrex::MultiFab::Copy(xvel[lev], xvel_old[lev], 0, 0, 1, Nghost);
        amrex::MultiFab::Copy(yvel[lev], yvel_old[lev], 0, 0, 1, Nghost);
        if(TempField) amrex::MultiFab::Copy(Theta[lev], Theta_old[lev], 0, 0, 1, Nghost);
        if(PhaseField) amrex::MultiFab::Copy(Phi[lev], Phi_old[lev], 0, 0, 1, Nghost);
	if(DamageModel)
        {
	    for(int iscalar = 0;iscalar < nscalar;iscalar++)
	        amrex::MultiFab::Copy(Scalars[lev][iscalar], Scalars_old[lev][iscalar], 0, 0, 1, Nghost);
	}

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
        
        //amrex::Print()<<"RKStage = "<<RKStage<<'\n';

        for (int step = 1;step <= RKStage; step++)
        {
            {
                //amrex::Print()<<"RKStage = "<<RKStage<<" , step = "<<step<<'\n';
                //amrex::Print()<<"rk_a[RKStage - 1][step - 1] = "<<rk_a[RKStage - 1][step - 1]<<'\n';
                //amrex::Print()<<"rk_c[RKStage - 1] = "<<rk_c[RKStage - 1]<<'\n';
                ComputeRKSourceTerms(lev, RKStage, step);
		if(DamageModel) ComputeRKSourceTerms_DamageModel(lev, RKStage, step);
                /// compute x-component
                for(amrex::MFIter mfi(xvel[lev]); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx = mfi.validbox();
                    //amrex::Array4<amrex::Real> const &uold = xvel_old[lev].array(mfi);
                    amrex::Array4<amrex::Real> const &ustar = xvel[lev].array(mfi);
                    amrex::Array4<amrex::Real> const &src_u = Src_xvel[lev].array(mfi);
     
                    amrex::Box box = get_valid_face_box(lev, bx, 0);
     
                    amrex::Array4<int> umask;
     
                    if(N_IF == 0)
                        umask = UMask_NI.array(mfi);
                    else
                        umask = UMask->array(mfi);
     
                    amrex::ParallelFor(box,
                    [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        if(umask(i, j, k) == 1)
                        {
                            ustar(i, j, k) += dt * rk_a[RKStage - 1][step - 1] * src_u(i, j, k);
                        }
                    }); 
                }
     
                /// compute y-component
                for(amrex::MFIter mfi(yvel[lev]); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx = mfi.validbox(); 
                    //amrex::Array4<amrex::Real> const &vold = yvel_old[lev].array(mfi);
                    amrex::Array4<amrex::Real> const &vstar = yvel[lev].array(mfi);
                    amrex::Array4<amrex::Real> const &src_v = Src_yvel[lev].array(mfi);
     
                    amrex::Box box = get_valid_face_box(lev, bx, 1);
     
                    amrex::Array4<int> vmask;
                    if(N_IF == 0)
                        vmask = VMask_NI.array(mfi);
                    else
                        vmask = VMask->array(mfi);
     
                    amrex::ParallelFor(box,
                        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        if(vmask(i ,j ,k) == 1)
                        {
                            if (isAxisymmetric)
                            {
                                //! box is nodal in y
                                amrex::Real yn = prob_lo[1] + dx[1] * (j + 0.5);
                                amrex::Real y = prob_lo[1] + dx[1] * j;
                                amrex::Real ys = prob_lo[1] + dx[1] * (j - 0.5);
                                
                                //! avoid singularity
                                if (fabs(y) <  1e-10)
                                {
                                    vstar(i, j, k) = 0.0;
                                }
                                else
                                {
                                    vstar(i, j, k) += dt * rk_a[RKStage - 1][step - 1] * src_v(i, j, k);
                                }
                            }
                            else
                            {
                                vstar(i, j, k) += dt * rk_a[RKStage - 1][step - 1] * src_v(i, j, k);
                            }
                        }
                    });
                }
     
                if(TempField)
                {
                /// compute Temperature
                    for(amrex::MFIter mfi(Theta[lev]); mfi.isValid(); ++mfi)
                    {
                        const amrex::Box& bx = mfi.validbox();
                        //amrex::Array4<amrex::Real const> const &Told = Theta_old[lev].const_array(mfi);
                        amrex::Array4<amrex::Real> const &T = Theta[lev].array(mfi);
                        amrex::Array4<amrex::Real> const &src_T = Src_Theta[lev].array(mfi);
     
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
                                T(i, j, k) += dt * rk_a[RKStage - 1][step - 1] * src_T(i, j, k);
                            }
                        });
                    }           
                }
     
                if(PhaseField)
                {
                /// compute Temperature
                    for(amrex::MFIter mfi(Phi[lev]); mfi.isValid(); ++mfi)
                    {
                        const amrex::Box& bx = mfi.validbox();
                        //amrex::Array4<amrex::Real const> const &phiold = Phi_old[lev].const_array(mfi);
                        amrex::Array4<amrex::Real> const &phi = Phi[lev].array(mfi);
                        amrex::Array4<amrex::Real> const &src_phi = Src_Phi[lev].array(mfi);
     
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
                                phi(i, j, k) += dt * rk_a[RKStage - 1][step - 1] * src_phi(i, j, k);
                            }
                        });
                    }
                }

                if(DamageModel)
                {
                /// Advect F_ij
		    for(int iscalar = 0;iscalar < nscalar;iscalar++)
	            {
                        if(iscalar == eps_max_num)
                            continue;
                        if(advect_ref_cond && (iscalar >= F11_num && iscalar <= F33_num))
                            continue;

                        for(amrex::MFIter mfi(Scalars[lev][iscalar]); mfi.isValid(); ++mfi)
                        {
                            const amrex::Box& bx = mfi.validbox();
                            //amrex::Array4<amrex::Real const> const &phiold = Phi_old[lev].const_array(mfi);
                            amrex::Array4<amrex::Real> const &scalar = Scalars[lev][iscalar].array(mfi);
                            amrex::Array4<amrex::Real> const &src_scalar = Src_Scalars[lev][iscalar].array(mfi);
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
                                    scalar(i, j, k) += dt * rk_a[RKStage - 1][step - 1] * src_scalar(i, j, k);
                                }
                            });
                        }
		    }
                }
            }
        }
    }

    void incFSI::ComputeRKSourceTerms(int lev, int RKStage, int step)
    {
        //! create temporary field data
        Src_xvel[lev].setVal(0.0);
        Src_yvel[lev].setVal(0.0);
        if(TempField)Src_Theta[lev].setVal(0.0);
        if(PhaseField)Src_Phi[lev].setVal(0.0);
 
        //amrex::Print()<<"Copute source RK3:"<<"lev = "<<lev<<" , RKStage = "<<RKStage<<" , step = "<<step<<'\n';
 
        //if(RKStage == 3 && step == 1) return;
        
        amrex::BoxArray xba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(0));
        amrex::BoxArray yba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(1));
        amrex::MultiFab tmp_xvel(xba, dmap[lev], 1, Nghost);
        amrex::MultiFab tmp_yvel(yba, dmap[lev], 1, Nghost);
        amrex::MultiFab tmp_Pressure(grids[lev], dmap[lev], 1, Nghost);
        amrex::MultiFab tmp_Theta(grids[lev], dmap[lev], 1, Nghost);
        amrex::MultiFab tmp_Phi(grids[lev], dmap[lev], 1, Nghost);
        
        tmp_xvel.setVal(0.0);
        tmp_yvel.setVal(0.0);
        tmp_Pressure.setVal(0.0);
        tmp_Theta.setVal(0.0);
        tmp_Phi.setVal(0.0);
        amrex::Vector<amrex::MultiFab *> cP(1), cTheta(1), cPhi(1);
        amrex::Vector<amrex::MultiFab *> fP(1), fTheta(1), fPhi(1);

        amrex::Vector<amrex::Real> ctime(2);
        amrex::Vector<amrex::Real> ftime(2);

        amrex::Vector<amrex::MultiFab *> cxvel(1);
        amrex::Vector<amrex::MultiFab *> cyvel(1);

        amrex::Vector<amrex::MultiFab *> fxvel(1);
        amrex::Vector<amrex::MultiFab *> fyvel(1);

        amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> tmp_vel = {&tmp_xvel, &tmp_yvel};
        amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> cvel;
        amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> fvel;

        if(step == 1)
        {
            //amrex::Print()<<"source step: "<<step<<'\n';
            /// fill internal ghost values for old
            xvel_old[lev].FillBoundary();
            yvel_old[lev].FillBoundary();
            Pressure[lev].FillBoundary();
            if(TempField)Theta_old[lev].FillBoundary();
            if(PhaseField)Phi_old[lev].FillBoundary();

            if (lev > 0)
            {
                ctime[0] = t_old;
                ctime[1] = t_new;

                cP[0] = &Pressure[lev - 1];
                if(TempField)cTheta[0] = &Theta_old[lev - 1];
                if(PhaseField)cPhi[0] = &Phi_old[lev - 1];

                cxvel[0] = &xvel_old[lev - 1];
                cyvel[0] = &yvel_old[lev - 1];

                if(divergence_free_interpolation)
                {   
                    cvel.resize(1);
                    cvel[0] = {&xvel_old[lev - 1], &yvel_old[lev - 1]};
                }
            
                ftime[0] = t_old;
                ftime[1] = t_new;
                
                fP[0] = &Pressure[lev];
                if(TempField)fTheta[0] = &Theta_old[lev];
                if(PhaseField)fPhi[0] = &Phi_old[lev];
                
                fxvel[0] = &xvel_old[lev];
                fyvel[0] = &yvel_old[lev];
                
                if(divergence_free_interpolation)
                {
                    fvel.resize(1);
                    fvel[0] = {&xvel_old[lev], &yvel_old[lev]};
                }
            }
        }
        else if(step == 2)
        {
            //amrex::Print()<<"source step: "<<step<<'\n';
            /// fill internal ghost values for old
            RK1_xvel[lev].FillBoundary();
            RK1_yvel[lev].FillBoundary();
            RK1_p[lev].FillBoundary();
            if(TempField)RK1_Theta[lev].FillBoundary();
            if(PhaseField)RK1_Phi[lev].FillBoundary();

            if (lev > 0)
            {
                ctime[0] = t_old;
                ctime[1] = t_new;

                cP[0] = &RK1_p[lev - 1];
                if(TempField)cTheta[0] = &RK1_Theta[lev - 1];
                if(PhaseField)cPhi[0] = &RK1_Phi[lev - 1];

                cxvel[0] = &RK1_xvel[lev - 1];
                cyvel[0] = &RK1_yvel[lev - 1];

                if(divergence_free_interpolation)
                {
                    cvel.resize(1);
                    cvel[0] = {&RK1_xvel[lev - 1], &RK1_yvel[lev - 1]};
                }
                ftime[0] = t_old;
                ftime[1] = t_new;
                
                fP[0] = &RK1_p[lev];
                if(TempField)fTheta[0] = &RK1_Theta[lev];
                if(PhaseField)fPhi[0] = &RK1_Phi[lev];
                
                fxvel[0] = &RK1_xvel[lev];
                fyvel[0] = &RK1_yvel[lev];
                
                if(divergence_free_interpolation)
                {
                    fvel.resize(1);
                    fvel[0] = {&RK1_xvel[lev], &RK1_yvel[lev]};
                }
            }
        }
        else if(step == 3)
        {
            //amrex::Print()<<"source step: "<<step<<'\n';
            /// fill internal ghost values for old
            RK2_xvel[lev].FillBoundary();
            RK2_yvel[lev].FillBoundary();
            RK2_p[lev].FillBoundary();
            if(TempField)RK2_Theta[lev].FillBoundary();
            if(PhaseField)RK2_Phi[lev].FillBoundary();

            if (lev > 0)
            {
                ctime[0] = t_old;
                ctime[1] = t_new;

                cP[0] = &RK2_p[lev - 1];
                if(TempField)cTheta[0] = &RK2_Theta[lev - 1];
                if(PhaseField)cPhi[0] = &RK2_Phi[lev - 1];

                cxvel[0] = &RK2_xvel[lev - 1];
                cyvel[0] = &RK2_yvel[lev - 1];

                if(divergence_free_interpolation)
                {
                    cvel.resize(1);
                    cvel[0] = {&RK2_xvel[lev - 1], &RK2_yvel[lev - 1]};
                }
                ftime[0] = t_old;
                ftime[1] = t_new;
                
                fP[0] = &RK2_p[lev];
                if(TempField)fTheta[0] = &RK2_Theta[lev];
                if(PhaseField)fPhi[0] = &RK2_Phi[lev];
                
                fxvel[0] = &RK2_xvel[lev];
                fyvel[0] = &RK2_yvel[lev];
                
                if(divergence_free_interpolation)
                {
                    fvel.resize(1);
                    fvel[0] = {&RK2_xvel[lev], &RK2_yvel[lev]};
                }
            }
        }
        else if(step == 4)
        {
            //amrex::Print()<<"source step: "<<step<<'\n';
            /// fill internal ghost values for old
            RK3_xvel[lev].FillBoundary();
            RK3_yvel[lev].FillBoundary();
            RK3_p[lev].FillBoundary();
            if(TempField)RK3_Theta[lev].FillBoundary();
            if(PhaseField)RK3_Phi[lev].FillBoundary();

            if (lev > 0)
            {
                ctime[0] = t_old;
                ctime[1] = t_new;

                cP[0] = &RK3_p[lev - 1];
                if(TempField)cTheta[0] = &RK3_Theta[lev - 1];
                if(PhaseField)cPhi[0] = &RK3_Phi[lev - 1];

                cxvel[0] = &RK3_xvel[lev - 1];
                cyvel[0] = &RK3_yvel[lev - 1];

                if(divergence_free_interpolation)
                {
                    cvel.resize(1);
                    cvel[0] = {&RK3_xvel[lev - 1], &RK3_yvel[lev - 1]};
                }
                ftime[0] = t_old;
                ftime[1] = t_new;

                fP[0] = &RK3_p[lev];
                if(TempField)fTheta[0] = &RK3_Theta[lev];
                if(PhaseField)fPhi[0] = &RK3_Phi[lev];

                fxvel[0] = &RK3_xvel[lev];
                fyvel[0] = &RK3_yvel[lev];

                if(divergence_free_interpolation)
                {
                    fvel.resize(1);
                    fvel[0] = {&RK3_xvel[lev], &RK3_yvel[lev]};
                }
            }
        }

        if (lev > 0)
        {
            if(divergence_free_interpolation)
            {
                //amrex::Print()<<"In divergence_free_interpolation "<<'\n'; 

                FillPatch(lev, Time, tmp_vel, cvel, {Time}, fvel, {Time}, 0, 1);
        
                amrex::MultiFab::Copy(tmp_xvel, *tmp_vel[0], 0, 0, 1, Nghost);
                amrex::MultiFab::Copy(tmp_yvel, *tmp_vel[1], 0, 0, 1, Nghost);
            }
            else
            {
                FillPatch(lev, Time, tmp_xvel, {cxvel[0]}, {ctime[0]}, {fxvel[0]}, {ftime[0]}, 0, 1);
                FillPatch(lev, Time, tmp_yvel, {cyvel[0]}, {ctime[0]}, {fyvel[0]}, {ftime[0]}, 0, 1);
            }
        
            //! no need of time interpolation of following
            FillPatch(lev, Time, tmp_Pressure, cP, {ctime[0]}, fP, {ftime[0]}, 0, 1);
            if(TempField) FillPatch(lev, Time, tmp_Theta, cTheta, {ctime[0]}, fTheta, {ftime[0]}, 0, 1);
            if(PhaseField) FillPatch(lev, Time, tmp_Phi, cPhi, {ctime[0]}, fPhi, {ftime[0]}, 0, 1);
            //! need to fill physical bcs
            XVelBoundaryConditions(lev, Time, tmp_xvel);
            YVelBoundaryConditions(lev, Time, tmp_yvel);
            PressureBoundaryConditions(lev, Time, tmp_Pressure);
            if(TempField) TemperatureBoundaryConditions(lev, Time, tmp_Theta);
            if(PhaseField) PhiBoundaryConditions(lev, Time, tmp_Phi);
        }
        else
        {
            if(step == 1)
            {
                amrex::MultiFab::Copy(tmp_xvel, xvel_old[lev], 0, 0, 1, Nghost);
                amrex::MultiFab::Copy(tmp_yvel, yvel_old[lev], 0, 0, 1, Nghost);
                amrex::MultiFab::Copy(tmp_Pressure, Pressure[lev], 0, 0, 1, Nghost);
                if(TempField) amrex::MultiFab::Copy(tmp_Theta, Theta_old[lev], 0, 0, 1, Nghost);
                if(PhaseField) amrex::MultiFab::Copy(tmp_Phi, Phi_old[lev], 0, 0, 1, Nghost);
            }
            else if(step == 2)
            {
                amrex::MultiFab::Copy(tmp_xvel, RK1_xvel[lev], 0, 0, 1, Nghost);
                amrex::MultiFab::Copy(tmp_yvel, RK1_yvel[lev], 0, 0, 1, Nghost);
                amrex::MultiFab::Copy(tmp_Pressure, RK1_p[lev], 0, 0, 1, Nghost);
                if(TempField) amrex::MultiFab::Copy(tmp_Theta, RK1_Theta[lev], 0, 0, 1, Nghost);
                if(PhaseField) amrex::MultiFab::Copy(tmp_Phi, RK1_Phi[lev], 0, 0, 1, Nghost);
            }
            else if(step == 3)
            {
                amrex::MultiFab::Copy(tmp_xvel, RK2_xvel[lev], 0, 0, 1, Nghost);
                amrex::MultiFab::Copy(tmp_yvel, RK2_yvel[lev], 0, 0, 1, Nghost);
                amrex::MultiFab::Copy(tmp_Pressure, RK2_p[lev], 0, 0, 1, Nghost);
                if(TempField) amrex::MultiFab::Copy(tmp_Theta, RK2_Theta[lev], 0, 0, 1, Nghost);
                if(PhaseField) amrex::MultiFab::Copy(tmp_Phi, RK2_Phi[lev], 0, 0, 1, Nghost);
            }
            else if(step == 4)
            {
                amrex::MultiFab::Copy(tmp_xvel, RK3_xvel[lev], 0, 0, 1, Nghost);
                amrex::MultiFab::Copy(tmp_yvel, RK3_yvel[lev], 0, 0, 1, Nghost);
                amrex::MultiFab::Copy(tmp_Pressure, RK3_p[lev], 0, 0, 1, Nghost);
                if(TempField) amrex::MultiFab::Copy(tmp_Theta, RK3_Theta[lev], 0, 0, 1, Nghost);
                if(PhaseField) amrex::MultiFab::Copy(tmp_Phi, RK3_Phi[lev], 0, 0, 1, Nghost);
            }

        }

        
        const amrex::Real *dx = geom[lev].CellSize();
        const amrex::Real *prob_lo = geom[lev].ProbLo();
    
        //amrex::Print()<<"dx = "<<dx[0]<<" , "<<dx[1]<<'\n';
        //amrex::Print()<<"prob_lo = "<<prob_lo[0]<<" , "<<prob_lo[1]<<'\n'; 
    
        amrex::iMultiFab *UMask;
        amrex::iMultiFab *VMask;
        amrex::iMultiFab *PMask;
	amrex::iMultiFab *Mask;
    
        amrex::iMultiFab UMask_NI;//NI stands for No Interface, need a better way of doing this
        amrex::iMultiFab VMask_NI;//NI stands for No Interface, need a better way of doing this
        amrex::iMultiFab PMask_NI;//NI stands for No Interface, need a better way of doing this
	amrex::iMultiFab Mask_NI;
            
        if(N_IF > 0)
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
            PMask_NI.define(yba, dmap[lev], 1, Nghost);
	    Mask_NI.define(yba, dmap[lev], 1, Nghost);
      
            UMask_NI.setVal(1);
            VMask_NI.setVal(1);
            PMask_NI.setVal(1);
	    Mask_NI.setVal(1);
        }
    
        Mu_max = std::max(Mu_max,Mu); 
	//amrex::Print()<<"Mu = "<<Mu<<" , Mu_max = "<<Mu_max<<'\n';
        /// compute x-component
        for(amrex::MFIter mfi(Src_xvel[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
    
            amrex::Array4<amrex::Real const> const &u = tmp_xvel.const_array(mfi);
            amrex::Array4<amrex::Real const> const &v = tmp_yvel.const_array(mfi);
            amrex::Array4<amrex::Real const> const &P = tmp_Pressure.const_array(mfi);
            amrex::Array4<amrex::Real> const &src_u = Src_xvel[lev].array(mfi);
            amrex::Box box = get_valid_face_box(lev, bx, 0);
    
            amrex::Array4<int> umask;
    
            if(N_IF == 0)
                umask = UMask_NI.array(mfi);
            else
                umask = UMask->array(mfi);

            amrex::Array4<amrex::Real const> phi;
            if(PhaseField)
                phi = tmp_Phi.array(mfi);

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
                    amrex::Real phi_e;
                    if(PhaseField)
                        phi_e = phi(i, j, k);
                    else
                        phi_e = 0.0;
                    amrex::Real Mu_e = Mu;
		    if(PhaseField) Mu_e = visc_.GetViscosity(S11_e, S12_e, S22_e, phi_e);
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
                    amrex::Real phi_w;
                    if(PhaseField)
                        phi_w = phi(i - 1, j, k);
                    else
                        phi_w = 0.0;
                    amrex::Real Mu_w = Mu;
		    if(PhaseField)
		        Mu_w = visc_.GetViscosity(S11_w, S12_w, S22_w, phi_w);

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
                    amrex::Real phi_n;
                    if(PhaseField)
                        phi_n = 0.25*(phi(i, j, k) + phi(i, j+1, k) + phi(i+1, j, k) + phi(i+1, j+1, k));
                    else
                        phi_n = 0.0;
                    amrex::Real Mu_n = Mu;
                    if(PhaseField)
                        Mu_n = visc_.GetViscosity(S11_n, S12_n, S22_n, phi_n);
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
                    amrex::Real phi_s;
                    if(PhaseField)
                        phi_s = 0.25*(phi(i, j, k) + phi(i, j - 1, k) + phi(i - 1, j, k) + phi(i - 1, j - 1, k));
                    else
                        phi_s = 0.0;
                    amrex::Real Mu_s = Mu;
                    if(PhaseField)
                        Mu_s = visc_.GetViscosity(S11_s, S12_s, S22_s, phi_s);
                    amrex::Real fs = us * vs - Mu_s * (uy_s + vx_s);

                    Mu_max = std::max(0.25*(Mu_e + Mu_w + Mu_n + Mu_s), Mu_max);                    
                    if (isAxisymmetric)
                    {
                        //TODO how to avoid singularity? i.e. what if y = 0
                        //! box is centered in y
                        amrex::Real yn = prob_lo[1] + dx[1] * (j + 1);
                        amrex::Real y  = prob_lo[1] + dx[1] * (j + 0.5);
                        amrex::Real ys = prob_lo[1] + dx[1] * j;
                    
                        fn *= yn;
                        fs *= ys;
                        src_u(i, j, k) = ( (fw - fe) / dx[0] + (fs - fn) / (y * dx[1]) ); 
                    }
                    else
                    {
                        src_u(i, j, k) = ( (fw - fe) / dx[0] + (fs - fn) / dx[1] );
                    }
                    src_u(i, j, k) += gravity;
                    
                    if (Projection)
                    {
                        amrex::Real dpdx = (P(i, j, k) - P(i - 1, j, k)) / dx[0];
                        src_u(i, j, k) -= dpdx;
                    }
                }
            });
        }
    
        /// compute y-component
        for(amrex::MFIter mfi(yvel[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<amrex::Real const> const &u = tmp_xvel.const_array(mfi);
            amrex::Array4<amrex::Real const> const &v = tmp_yvel.const_array(mfi);
            amrex::Array4<amrex::Real const> const &P = tmp_Pressure.const_array(mfi);
            amrex::Array4<amrex::Real> const &src_v = Src_yvel[lev].array(mfi);
            amrex::Box box = get_valid_face_box(lev, bx, 1);
    
            amrex::Array4<int> vmask;
            if(N_IF == 0)
                vmask = VMask_NI.array(mfi);
            else
                vmask = VMask->array(mfi);
    
            amrex::Array4<amrex::Real const> phi;
            if(PhaseField)
                phi = tmp_Phi.array(mfi);

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
                    amrex::Real phi_e;
                    if(PhaseField)
                        phi_e = 0.25 * (phi(i, j, k) + phi(i + 1, j, k) + phi(i, j - 1, k) + phi(i + 1, j - 1, k));
                    else
                        phi_e = 0.0;
                    amrex::Real Mu_e = Mu;
                    if(PhaseField)
                        Mu_e = visc_.GetViscosity(S11_e, S12_e, S22_e, phi_e);
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
                    amrex::Real phi_w;
                    if(PhaseField)
                        phi_w = 0.25 * (phi(i - 1, j, k) + phi(i, j, k) + phi(i - 1, j - 1, k) + phi(i, j - 1, k));
                    else
                        phi_w = 0.0;
                    amrex::Real Mu_w = Mu;
                    if(PhaseField)
                        Mu_w = visc_.GetViscosity(S11_w, S12_w, S22_w, phi_w);
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
                    amrex::Real phi_n;
                    if(PhaseField)
                        phi_n = phi(i, j, k);
                    else
                        phi_n = 0.0;

                    amrex::Real Mu_n = Mu;
                    if(PhaseField)
                        Mu_n = visc_.GetViscosity(S11_n, S12_n, S22_n, phi_n);
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
                    amrex::Real phi_s;
                    if(PhaseField)
                        phi_s = phi(i, j - 1, k);
                    else
                        phi_s = 0.0;
                    amrex::Real Mu_s = Mu;
                    if(PhaseField)
                        Mu_s = visc_.GetViscosity(S11_s, S12_s, S22_s, phi_s);
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
                            src_v(i, j, k) = 0.0;
                        }
                        else
                        {
                            src_v(i, j, k) = ((fw - fe) / dx[0] + (fs - fn) / (y * dx[1]) - 2.0 * Mu_ * v(i, j, k) / (y * y));
                        }
                    }
                    else
                    {
                        src_v(i, j, k) = ((fw - fe) / dx[0] + (fs - fn) / dx[1]);
                    }
                    
                    if (Projection)
                    {
                        amrex::Real dpdy = (P(i, j, k) - P(i, j - 1, k)) / dx[1];
                        src_v(i, j, k) -= dpdy;
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
            amrex::Real alpha_ = 0.0;
            
            /// compute Temperature
            for(amrex::MFIter mfi(Theta[lev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<amrex::Real const> const &xvel = tmp_xvel.const_array(mfi);
                amrex::Array4<amrex::Real const> const &yvel = tmp_yvel.const_array(mfi);
                amrex::Array4<amrex::Real const> const &P = tmp_Pressure.const_array(mfi);
                amrex::Array4<amrex::Real const> const &Told = tmp_Theta.const_array(mfi);
                amrex::Array4<amrex::Real> const &src_T = Src_Theta[lev].array(mfi);
        
                amrex::Array4<int> mask;
                if(N_IF == 0)
                    mask = PMask_NI.array(mfi);
                else
                    mask = PMask->array(mfi);

                amrex::Array4<amrex::Real const> phi;
                if(PhaseField)
                    phi = tmp_Phi.array(mfi);
        
                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
                {
                    if(mask(i ,j ,k) == 1)
                    {
                        amrex::Real ux = (xvel(i + 1, j, k) - xvel(i, j, k)) / dx[0];
                        amrex::Real uy = (xvel(i + 1, j + 1, k) - xvel(i + 1, j - 1, k) + xvel(i, j + 1, k) - xvel(i, j - 1, k)) / (4.0 * dx[1]);
                        amrex::Real vx = (yvel(i + 1, j + 1, k) - yvel(i - 1, j + 1, k) + yvel(i + 1, j, k) - yvel(i - 1, j, k)) / (4.0 * dx[0]);
                        amrex::Real vy = (yvel(i, j + 1, k) - yvel(i, j, k)) / dx[1];
                        amrex::Real S11_ = ux;
                        amrex::Real S12_ = 0.5*(uy + vx);
                        amrex::Real S22_ = vy;

                        amrex::Real dGamma_dt = visc_.ComputeGammaDot(ux, 0.5 * (uy + vx), vy);
                        amrex::Real phi_;
                        if(PhaseField)
                            phi_ = phi(i, j, k);
                        else
                            phi_ = 0.0;
                        amrex::Real Mu_ = Mu;
                        if(PhaseField)
                            Mu_ = visc_.GetViscosity(S11_, S12_, S22_, phi_);

                        amrex::Real theta_e = 0.5 * (Told(i, j, k) + Told(i + 1, j, k));
                        amrex::Real thetax_e = (Told(i + 1, j, k) - Told(i, j, k)) / dx[0];
                        k_ = visc_.GetConductivity(Mu_, Prandtl_no);//incorrect implementation of heatflux with variable therma conductivity
                        if (mask(i,j,k) == 1 && mask(i+1,j,k) == 1)
                        {
                            k_max = std::max(k_, k_max);
                        }
                        //amrex::Real fe = theta_e * xvel(i + 1, j, k) - k_ * thetax_e;
                        amrex::Real fe = theta_e * xvel(i + 1, j, k) 
                                         + alpha_ * 0.5 * std::abs(xvel(i + 1, j, k)) * (Told(i, j, k) - Told(i + 1, j, k))
                                         - k_ * thetax_e;
        
                        amrex::Real theta_w = 0.5 * (Told(i, j, k) + Told(i - 1, j, k));
                        amrex::Real thetax_w = (Told(i, j, k) - Told(i - 1, j, k)) / dx[0];
                        k_ = visc_.GetConductivity(Mu_, Prandtl_no);
                        if (mask(i, j, k) == 1 && mask(i - 1, j, k) == 1)
                        {
                            k_max = std::max(k_, k_max);
                        }
                        //amrex::Real fw = theta_w * xvel(i, j, k) - k_ * thetax_w;
                        amrex::Real fw = theta_w * xvel(i, j, k) 
                                         + alpha_ * 0.5 *std:: abs(xvel(i, j, k)) * (Told(i - 1, j, k) - Told(i, j, k))
                                         - k_ * thetax_w;
                        
                        amrex::Real theta_n = 0.5 * (Told(i, j, k) + Told(i, j + 1, k));
                        amrex::Real thetay_n = (Told(i, j + 1, k) - Told(i, j, k)) / dx[1];
                        k_ = visc_.GetConductivity(Mu_, Prandtl_no);
                        if (mask(i, j, k) == 1 && mask(i, j + 1, k) == 1)
                        {
                            k_max = std::max(k_, k_max);
                        }
                        //amrex::Real fn = theta_n * yvel(i, j + 1, k) - k_ * thetay_n;
                        amrex::Real fn = theta_n * yvel(i, j + 1, k) 
                                         + alpha_ * 0.5 * std::abs(yvel(i, j + 1, k)) * (Told(i, j, k) - Told(i, j + 1, k))
                                         - k_ * thetay_n;
        
                        amrex::Real theta_s = 0.5 * (Told(i, j, k) + Told(i, j - 1, k));
                        amrex::Real thetay_s = (Told(i, j, k) - Told(i, j - 1, k)) / dx[1];
                        k_ = visc_.GetConductivity(Mu_, Prandtl_no);
                        if (mask(i, j, k) == 1 && mask(i, j - 1, k) == 1)
                        {
                            k_max = std::max(k_, k_max);
                        }
                        //amrex::Real fs = theta_s * yvel(i, j, k) - k_ * thetay_s;
                        amrex::Real fs = theta_s * yvel(i, j, k) 
                                         + alpha_ * 0.5 * std::abs(yvel(i, j, k)) * (Told(i, j - 1, k) - Told(i, j, k))
                                         - k_ * thetay_s;
        
                        if (isAxisymmetric)
                        {
                            //! box is cell centered
                            amrex::Real yn = prob_lo[1] + dx[1] * (j + 1);
                            amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                            amrex::Real ys = prob_lo[1] + dx[1] * j;
        
                            fn *= yn;
                            fs *= ys;
                            src_T(i, j, k) = ((fw - fe) / dx[0] + (fs - fn) / (y * dx[1]) +  Mu_ * (u_ref*u_ref/cp_/T_ref) * dGamma_dt* dGamma_dt);
                        }
                        else
                        {
                            src_T(i, j, k) = ((fw - fe) / dx[0] + (fs - fn) / dx[1] +  Mu_ * (u_ref*u_ref/cp_/T_ref) * dGamma_dt * dGamma_dt);
                        }
                    }
                    if(mask(i,j,k) != 1 && mask(i,j,k) != 2)
                    {
                        //auto &solid = interfaces[0];
                        src_T(i, j, k) =  0.0;//P_interface[0] * solid->Volume() / (P_interface0[0] * solid->Volume_0());
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
                amrex::Array4<amrex::Real const> const &xvel = tmp_xvel.const_array(mfi);
                amrex::Array4<amrex::Real const> const &yvel = tmp_yvel.const_array(mfi);
                amrex::Array4<amrex::Real const> const &phiold = tmp_Phi.const_array(mfi);
                amrex::Array4<amrex::Real> const &src_phi = Src_Phi[lev].array(mfi);
            
                amrex::Array4<int> pmask, mask;
                if(N_IF == 0)
                    pmask = PMask_NI.array(mfi);
                else
                    pmask = PMask->array(mfi);
            
                if(N_IF == 0)
                    mask = Mask_NI.array(mfi);
                else
                    mask = Mask->array(mfi);

                amrex::Array4<amrex::Real const> dmg;
                if(DamageModel)
                    dmg = Scalars_old[lev][dmg_num].array(mfi);
            
                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
                {
                    if(pmask(i ,j ,k) == 1)
                    {
                        amrex::Real Phix_L, Phix_R,Phiy_L, Phiy_R, Phix, Phiy;
                        amrex::Real u_Mid = 0.5 * (xvel(i, j, k) + xvel(i + 1, j, k));
                        amrex::Real v_Mid = 0.5 * (yvel(i, j, k) + yvel(i, j + 1, k));

                        WENO5_LS(Phix_L, Phix_R, Phiy_L, Phiy_R, i, j, dx, phiold);
			if(mask(i, j, k) == 3)
			{
			    //amrex::Print(-1)<<"near interface"<<'\n';
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
			    if(isAxisymmetric)
		            {
				amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
			        laplacian_phi += (1.0/y)*(phiold(i , j + 1 , k) - phiold(i , j - 1 , k))/(2.0*dx[1]);
				amrex::Real ny = 0.5*(phiold(i, j + 1, k) - phiold(i, j - 1, k))/(grad_phi*dx[1] + 1e-15 );
				curvat += (1.0/y)*ny;
			    }
                            RHS_phi = b*(laplacian_phi + (phiold(i,j,k)*(1.0 - phiold(i,j,k)*phiold(i,j,k))/(4.0 * dx[0]*dx[0])) - grad_phi*curvat);
                            //if(std::isnan(RHS_phi))
                            {
                                //amrex::Print()<<"RHS_phi = "<<RHS_phi<<'\n';
                                //amrex::Print()<<"curvat = "<<curvat<<'\n';
                                //exit(3);
                            }
                        }

                        //if (isAxisymmetric)
                        //{
                        //    //! box is cell centered
                        //    amrex::Real yn = prob_lo[1] + dx[1] * (j + 1);
                        //    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                        //    amrex::Real ys = prob_lo[1] + dx[1] * j;
            
                            //fn *= yn;
                            //fs *= ys;
                            //src_phi(i, j, k) = ((fw - fe) / dx[0] + (fs - fn) / (y * dx[1]));
			//    src_phi(i, j, k) = -1.0*(u_Mid * Phix + v_Mid * Phiy); 
                        //}
                        //else
                        //{
                        //    //src_phi(i, j, k) = ((fw - fe) / dx[0] + (fs - fn) / dx[1]);
			//    src_phi(i, j, k) = -1.0*(u_Mid * Phix + v_Mid * Phiy);
                        //
                        //}
			src_phi(i, j, k) = -1.0*(u_Mid * Phix + v_Mid * Phiy);
			if(sharpPhaseField)src_phi(i, j, k) += RHS_phi;
			if(DamageModel)src_phi(i, j, k) += -1.0*dmg(i, j, k)*(damage_coeff/dx[1])*(1.0 + phiold(i, j, k));
                    }
                });
            }       
        }
    }
} // namespace mycode
