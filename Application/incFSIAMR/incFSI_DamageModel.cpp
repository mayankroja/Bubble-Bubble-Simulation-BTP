#include "incFSI.H"
#include <cmath>
#include <complex>
#include <iostream>
#include <ostream>
#include <WeightedENO.h>
#include <AdvectLS_helper.H>
#include<SolveCubicEqn.h>

namespace mycode
{

    void incFSI::ComputeDeformationGradientTensor2D(int lev)
    {
        Scalars[lev][0].FillBoundary();
        Scalars[lev][1].FillBoundary();
        Scalars[lev][dmg_num].setVal(0.0);
        const amrex::Real *prob_lo = geom[lev].ProbLo();
        const amrex::Real *dx = geom[lev].CellSize();
    
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
    
    
    
        for(amrex::MFIter mfi(Scalars[lev][0]); mfi.isValid(); ++mfi)
        {
           amrex::Array4<int> mask;
           if(N_IF == 0)
               mask = PMask_NI.array(mfi);
           else
               mask = PMask->array(mfi);
    
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<amrex::Real const> const &X = Scalars[lev][0].const_array(mfi);
            amrex::Array4<amrex::Real const> const &Y = Scalars[lev][1].const_array(mfi);
            amrex::Array4<amrex::Real> const &F11 = Scalars[lev][2].array(mfi);
            amrex::Array4<amrex::Real> const &F12 = Scalars[lev][3].array(mfi);
            amrex::Array4<amrex::Real> const &F21 = Scalars[lev][4].array(mfi);
            amrex::Array4<amrex::Real> const &F22 = Scalars[lev][5].array(mfi);
            amrex::Array4<amrex::Real> const &F33 = Scalars[lev][6].array(mfi);
            amrex::Array4<amrex::Real> const &ep_max = Scalars[lev][7].array(mfi);
    	    amrex::Array4<amrex::Real> const &dmg = Scalars[lev][8].array(mfi);
	    amrex::Array4<amrex::Real> const &frk_F11 = FRK_Scalars[lev][2].array(mfi);
	    amrex::Array4<amrex::Real> const &frk_F12 = FRK_Scalars[lev][3].array(mfi);
	    amrex::Array4<amrex::Real> const &frk_F21 = FRK_Scalars[lev][4].array(mfi);
	    amrex::Array4<amrex::Real> const &frk_F22 = FRK_Scalars[lev][5].array(mfi);
	    amrex::Array4<amrex::Real> const &frk_F33 = FRK_Scalars[lev][6].array(mfi);
	    amrex::Array4<amrex::Real> const &phi = Phi[lev].array(mfi);
    
    
            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                /* 
    	        amrex::Real F_inv_11, F_inv_12, F_inv_21, F_inv_22;
    
                 
                amrex::Real Psix_L, Psix_R, Psiy_L, Psiy_R, SignPsi;
                {
                    //WENO5_LS(Psix_L, Psix_R, Psiy_L, Psiy_R, i, j, dx, X);
                }
                {
                    Psix_R = (X(i + 1, j, k) - X(i, j, k)) / dx[0];
                    Psix_L = (X(i, j, k) - X(i - 1, j, k)) / dx[0];
                    Psiy_R = (X(i, j + 1, k) - X(i, j, k)) / dx[1];
                    Psiy_L = (X(i, j, k) - X(i, j - 1, k)) / dx[1];
    
    		//if(X(i + 1, j, k)*X(i - 1, j, k) > 0.0)
    		//    F_inv_11 = 0.5*(Psix_L + Psix_R);
    		//else
    		{
                        if(std::abs(Psix_L)<std::abs(Psix_R))
                            F_inv_11 = Psix_L;
                        else
                            F_inv_11 = Psix_R;
                    }
                    //if(X(i, j + 1, k)*X(i, j - 1, k) > 0.0)
                    //    F_inv_21 = 0.5*(Psiy_L + Psiy_R);
                    //else
                    {
                        if(std::abs(Psiy_L)<std::abs(Psiy_R))
                            F_inv_21 = Psiy_L;
                        else
                            F_inv_21 = Psiy_R;
                    }
                }
                {
                    Psix_R = (Y(i + 1, j, k) - Y(i, j, k)) / dx[0];
                    Psix_L = (Y(i, j, k) - Y(i - 1, j, k)) / dx[0];
                    Psiy_R = (Y(i, j + 1, k) - Y(i, j, k)) / dx[1];
                    Psiy_L = (Y(i, j, k) - Y(i, j - 1, k)) / dx[1];
    
                    //if(Y(i + 1, j, k)*Y(i - 1, j, k) > 0.0)
                    //    F_inv_12 = 0.5*(Psix_L + Psix_R);
                    //else
                    {   
                        if(std::abs(Psix_L)<std::abs(Psix_R))
                            F_inv_12 = Psix_L;
                        else
                            F_inv_12 = Psix_R;
                    }
                    //if(Y(i, j + 1, k)*Y(i, j - 1, k) > 0.0)
                    //    F_inv_22 = 0.5*(Psiy_L + Psiy_R);
                    //else
                    {
                        if(std::abs(Psiy_L)<std::abs(Psiy_R))
                            F_inv_22 = Psiy_L;
                        else
                            F_inv_22 = Psiy_R;
                    }
                }
    
    	    
                //}
    	        amrex::Real Det_F_inv = F_inv_11*F_inv_22 - F_inv_12*F_inv_21 + 1e-20;
    	        if(advect_ref_cond)
    	        {
                    F11(i, j, k) = (1.0/Det_F_inv)*F_inv_22;
                    F12(i, j, k) = (-1.0/Det_F_inv)*F_inv_12;
                    F21(i, j, k) = (-1.0/Det_F_inv)*F_inv_21;
                    F22(i, j, k) = (1.0/Det_F_inv)*F_inv_11;
    	        }
		*/
                //Test Fij from exact solution
	        /* 
	        {
	            auto &solid = interfaces[finest_level][0];
                    amrex::Real x_gc = prob_lo[0] + dx[0] * (i + 0.5);
                    amrex::Real y_gc = prob_lo[1] + dx[1] * (j + 0.5);
                    amrex::Real r_gc = std::hypot(x_gc - solid->Xcp(), y_gc - solid->Ycp());
                    amrex::Real theta = std::acos((x_gc - solid->Xcp())/r_gc);
		    amrex::Real r_b_0 = std::pow(solid->Volume_0()/(4.0/3.0*M_PI),(1.0/3.0));
		    amrex::Real r_b = std::pow(solid->Volume()/(4.0/3.0*M_PI),(1.0/3.0));
		    amrex::Real R = std::pow((r_gc*r_gc*r_gc - (r_b*r_b*r_b - r_b_0*r_b_0*r_b_0)),1.0/3.0);
                    amrex::Real R_by_r = R/r_gc;///std::pow((solid->Volume_0() / solid->Volume()), 1.0/3.0); 

                    //amrex::Print()<<"r_gc = "<<r_gc<<" , theta = "<<theta<<'\n';

                    F11(i, j, k) = R_by_r*R_by_r*std::pow(std::cos(theta),2.0) + 1.0/R_by_r * std::pow(std::sin(theta),2.0);
                    F12(i, j, k) = (R_by_r*R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);
                    F21(i, j, k) = (R_by_r*R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);
                    F22(i, j, k) = R_by_r*R_by_r*std::pow(std::sin(theta),2.0) + 1.0/R_by_r * std::pow(std::cos(theta),2.0);
                    F33(i, j, k) = 1.0/R_by_r;
                    if(std::isnan(F11(i, j, k)) ||
                           std::isnan(F12(i, j, k)) ||
                           std::isnan(F21(i, j, k)) ||
                           std::isnan(F22(i, j, k)) ||
                           std::isnan(F33(i, j, k)))
		    {
		        amrex::Print()<<"i, j, k, lev = "<<i<<" , "<<j<<" , "<<k<<" , "<<lev<<'\n';
		    }
	        }*/	
                //Right Cauchy-Green Tensor
    	        amrex::Real C11 = F11(i, j, k)*F11(i, j, k) + F21(i, j, k)*F21(i, j, k);
    	        amrex::Real C12 = F11(i, j, k)*F12(i, j, k) + F21(i, j, k)*F22(i, j, k);
    	        amrex::Real C21 = F12(i, j, k)*F11(i, j, k) + F22(i, j, k)*F21(i, j, k);
    	        amrex::Real C22 = F12(i, j, k)*F12(i, j, k) + F22(i, j, k)*F22(i, j, k);
    	        amrex::Real C33 = F33(i, j, k)*F33(i, j, k);
    
    	        amrex::Real tr_C = C11+C22;
    	        amrex::Real det_C = C11*C22 - C21*C12;
    	        amrex::Real first_inv = C11+C22+C33;
    	        amrex::Real second_inv = C22*C33 + C11*C33 + C11*C22 - C12*C21;
    	        amrex::Real third_inv = C33*(C11*C22 - C12*C21);
    
    	        amrex::Real lambda_1;// = 0.5*(tr_C + std::sqrt(tr_C*tr_C - 4.0*det_C));
    	        amrex::Real lambda_2;// = 0.5*(tr_C - std::sqrt(tr_C*tr_C - 4.0*det_C));
    	        amrex::Real lambda_3;
    
    	        //amrex::Print()<<"first_inv = "<<first_inv<<" , second_inv = "<<second_inv<<" , third_inv = "<<third_inv<<'\n';
    	        //SolveCubicEqn(lambda_1, lambda_2, lambda_3, first_inv, second_inv, third_inv); 
    	        SolveQuadraticEqn(lambda_1,lambda_2, tr_C, det_C);
	        lambda_3 = C33;
                if(std::isnan(lambda_1)) lambda_1 = 1.0;
                if(std::isnan(lambda_2)) lambda_2 = 1.0;
                if(std::isnan(lambda_3)) lambda_3 = 1.0;
    
    	        //amrex::Print()<<"lambda_1 = "<<lambda_1<<", lambda_2 = "<<lambda_2<<" , lambda_3 = "<<lambda_3<<'\n';
    
    	        ep_max(i, j, k) = std::max(0.5*(lambda_1 - 1.0),
    			                   0.5*(lambda_2 - 1.0));

    	        ep_max(i, j, k) = std::max(ep_max(i, j, k),0.5*(lambda_3 - 1.0));
    
    	        //if(ep_max(i, j, k) > visc_.Lambda_f()) 
                //    dmg(i, j, k) = 1.0;
	        amrex::Real max_stretch  = std::max(std::sqrt(std::abs(lambda_1)),std::sqrt(std::abs(lambda_2)));
	        max_stretch = std::max(std::sqrt(std::abs(lambda_3)),max_stretch);
	        if(max_stretch > visc_.Lambda_f())
                    dmg(i, j, k) = 1.0;

    	        if(mask(i, j, k) != 1)
                {
                    ep_max(i, j, k) = 0.0;
                }
            });
        }
    }

    void incFSI::ComputeRKSourceTerms_DamageModel(int lev, int RKStage, int step)
    {
	//Fill patch for Xvel and Yvel
        amrex::BoxArray xba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(0));
        amrex::BoxArray yba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(1));
        amrex::MultiFab tmp_xvel(xba, dmap[lev], 1, Nghost);
        amrex::MultiFab tmp_yvel(yba, dmap[lev], 1, Nghost);
        
        tmp_xvel.setVal(0.0);
        tmp_yvel.setVal(0.0);

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

            if (lev > 0)
            {
                ctime[0] = t_old;
                ctime[1] = t_new;

                cxvel[0] = &xvel_old[lev - 1];
                cyvel[0] = &yvel_old[lev - 1];

                if(divergence_free_interpolation)
                {   
                    cvel.resize(1);
                    cvel[0] = {&xvel_old[lev - 1], &yvel_old[lev - 1]};
                }
            
                ftime[0] = t_old;
                ftime[1] = t_new;
                
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
            if (lev > 0)
            {
                ctime[0] = t_old;
                ctime[1] = t_new;

                cxvel[0] = &RK1_xvel[lev - 1];
                cyvel[0] = &RK1_yvel[lev - 1];

                if(divergence_free_interpolation)
                {
                    cvel.resize(1);
                    cvel[0] = {&RK1_xvel[lev - 1], &RK1_yvel[lev - 1]};
                }
                ftime[0] = t_old;
                ftime[1] = t_new;
                
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

            if (lev > 0)
            {
                ctime[0] = t_old;
                ctime[1] = t_new;

                cxvel[0] = &RK2_xvel[lev - 1];
                cyvel[0] = &RK2_yvel[lev - 1];

                if(divergence_free_interpolation)
                {
                    cvel.resize(1);
                    cvel[0] = {&RK2_xvel[lev - 1], &RK2_yvel[lev - 1]};
                }
                ftime[0] = t_old;
                ftime[1] = t_new;
                
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

            if (lev > 0)
            {
                ctime[0] = t_old;
                ctime[1] = t_new;

                cxvel[0] = &RK3_xvel[lev - 1];
                cyvel[0] = &RK3_yvel[lev - 1];

                if(divergence_free_interpolation)
                {
                    cvel.resize(1);
                    cvel[0] = {&RK3_xvel[lev - 1], &RK3_yvel[lev - 1]};
                }
                ftime[0] = t_old;
                ftime[1] = t_new;

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
            //! need to fill physical bcs
            XVelBoundaryConditions(lev, Time, tmp_xvel);
            YVelBoundaryConditions(lev, Time, tmp_yvel);
        }
        else
        {
            if(step == 1)
            {
                amrex::MultiFab::Copy(tmp_xvel, xvel_old[lev], 0, 0, 1, Nghost);
                amrex::MultiFab::Copy(tmp_yvel, yvel_old[lev], 0, 0, 1, Nghost);
            }
            else if(step == 2)
            {
                amrex::MultiFab::Copy(tmp_xvel, RK1_xvel[lev], 0, 0, 1, Nghost);
                amrex::MultiFab::Copy(tmp_yvel, RK1_yvel[lev], 0, 0, 1, Nghost);
            }
            else if(step == 3)
            {
                amrex::MultiFab::Copy(tmp_xvel, RK2_xvel[lev], 0, 0, 1, Nghost);
                amrex::MultiFab::Copy(tmp_yvel, RK2_yvel[lev], 0, 0, 1, Nghost);
            }
            else if(step == 4)
            {
                amrex::MultiFab::Copy(tmp_xvel, RK3_xvel[lev], 0, 0, 1, Nghost);
                amrex::MultiFab::Copy(tmp_yvel, RK3_yvel[lev], 0, 0, 1, Nghost);
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
    
	    //amrex::Print()<<"Mu = "<<Mu<<" , Mu_max = "<<Mu_max<<'\n';
        for(int iscalar = 0;iscalar < nscalar;iscalar++)
        {
	    //Set Source to zero
	    Src_Scalars[lev][iscalar].setVal(0.0);
            if(iscalar == eps_max_num)
                continue;
            if(advect_ref_cond && (iscalar >= F11_num && iscalar <= F33_num))
                continue;
            if(iscalar == dmg_num)
                continue;

		    
            amrex::MultiFab tmp_Scalar(grids[lev], dmap[lev], 1, Nghost);
            
            tmp_Scalar.setVal(0.0);
            amrex::Vector<amrex::MultiFab *> cScalar(1);
            amrex::Vector<amrex::MultiFab *> fScalar(1);
		    
            amrex::Vector<amrex::Real> ctime(2);
            amrex::Vector<amrex::Real> ftime(2);
            
            ctime[0] = t_old;
            ctime[1] = t_new;
            
            ftime[0] = t_old;
            ftime[1] = t_new;        
            if(step == 1)
            {
                Scalars_old[lev][iscalar].FillBoundary();
		    
                if (lev > 0)
                {
                    cScalar[0] = &Scalars_old[lev - 1][iscalar];
                    fScalar[0] = &Scalars_old[lev][iscalar];
                }
            }
            else if(step == 2)
            {
                RK1_Scalars[lev][iscalar].FillBoundary();
		    
                if (lev > 0)
                {
                    cScalar[0] = &RK1_Scalars[lev - 1][iscalar];
                    fScalar[0] = &RK1_Scalars[lev][iscalar];
                }
            }
            else if(step == 3)
            {
                RK2_Scalars[lev][iscalar].FillBoundary();
		    
                if (lev > 0)
                {
                    cScalar[0] = &RK2_Scalars[lev - 1][iscalar];
                    fScalar[0] = &RK2_Scalars[lev][iscalar];
                }
            }
            else if(step == 4)
            {
                RK3_Scalars[lev][iscalar].FillBoundary();
		    
                if (lev > 0)
                {
                    cScalar[0] = &RK3_Scalars[lev - 1][iscalar];
                    fScalar[0] = &RK3_Scalars[lev][iscalar];
                }
            }
		    
            if (lev > 0)
            {
		    
                FillPatch(lev, Time, tmp_Scalar, cScalar, {ctime[0]}, fScalar, {ftime[0]}, 0, 1);
                //! need to fill physical bcs
                ScalarBoundaryConditions(lev,iscalar, Time, tmp_Scalar);
            }
            else
            {
                if(step == 1)
                {
                    amrex::MultiFab::Copy(tmp_Scalar, Scalars_old[lev][iscalar], 0, 0, 1, Nghost);
                }
                else if(step == 2)
                {
                    amrex::MultiFab::Copy(tmp_Scalar, RK1_Scalars[lev][iscalar], 0, 0, 1, Nghost);
                }
                else if(step == 3)
                {
                    amrex::MultiFab::Copy(tmp_Scalar, RK2_Scalars[lev][iscalar], 0, 0, 1, Nghost);
                }
                else if(step == 4)
                {
                    amrex::MultiFab::Copy(tmp_Scalar, RK3_Scalars[lev][iscalar], 0, 0, 1, Nghost);
                }
            }


        /// Advect Scalar field
            amrex::Real alpha_ = 1.0;
            for(amrex::MFIter mfi(Scalars[lev][iscalar]); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<amrex::Real const> const &xvel = tmp_xvel.const_array(mfi);
                amrex::Array4<amrex::Real const> const &yvel = tmp_yvel.const_array(mfi);
                amrex::Array4<amrex::Real const> const &scalarold = tmp_Scalar.const_array(mfi);
		amrex::Array4<amrex::Real> const &src_scalar = Src_Scalars[lev][iscalar].array(mfi);
	        amrex::Array4<amrex::Real> F11;
	        amrex::Array4<amrex::Real> F12;
	        amrex::Array4<amrex::Real> F21;
	        amrex::Array4<amrex::Real> F22;

                if(step == 1)
                {
                    F11 =  Scalars_old[lev][2].array(mfi);
		    F12 =  Scalars_old[lev][3].array(mfi);
		    F21 =  Scalars_old[lev][4].array(mfi);
		    F22 =  Scalars_old[lev][5].array(mfi);
                }
                else if(step == 2)
                {
                    F11 =  RK1_Scalars[lev][2].array(mfi);
                    F12 =  RK1_Scalars[lev][3].array(mfi);
                    F21 =  RK1_Scalars[lev][4].array(mfi);
                    F22 =  RK1_Scalars[lev][5].array(mfi);
                }
                else if(step == 3)
                {
                    F11 =  RK2_Scalars[lev][2].array(mfi);
                    F12 =  RK2_Scalars[lev][3].array(mfi);
                    F21 =  RK2_Scalars[lev][4].array(mfi);
                    F22 =  RK2_Scalars[lev][5].array(mfi);
                }
                else if(step == 4)
                {
                    F11 =  RK3_Scalars[lev][2].array(mfi);
                    F12 =  RK3_Scalars[lev][3].array(mfi);
                    F21 =  RK3_Scalars[lev][4].array(mfi);
                    F22 =  RK3_Scalars[lev][5].array(mfi);
                }

                amrex::Array4<int> pmask, mask;
                if(N_IF == 0)
                    pmask = PMask_NI.array(mfi);
                else
                    pmask = PMask->array(mfi);

                if(N_IF == 0)
                    mask = Mask_NI.array(mfi);
                else
                    mask = Mask->array(mfi);
            
                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
                {
                    //if((mask(i ,j ,k) > 0 && iscalar < 2) || (mask(i ,j ,k) == 1 && iscalar == 8))
	    	    if(pmask(i ,j ,k) == 1)
                    {
                        /*
                        amrex::Real Sx_L, Sx_R, Sy_L, Sy_R, Sx, Sy;
                        amrex::Real u_Mid = 0.5 * (xvel(i, j, k) + xvel(i + 1, j, k));
                        amrex::Real v_Mid = 0.5 * (yvel(i, j, k) + yvel(i, j + 1, k));

                        WENO5_LS(Sx_L, Sx_R, Sy_L, Sy_R, i, j, dx, scalarold);
                        if(mask(i, j, k) == 3)
                        {
                            //amrex::Print(-1)<<"near interface"<<'\n';
                            Sx_R = (scalarold(i + 1, j, k) - scalarold(i, j, k)) / dx[0];
                            Sx_L = (scalarold(i, j, k) - scalarold(i - 1, j, k)) / dx[0];
                            Sy_R = (scalarold(i, j + 1, k) - scalarold(i, j, k)) / dx[1];
                            Sy_L = (scalarold(i, j, k) - scalarold(i, j - 1, k)) / dx[1];
                        }
                        if (u_Mid > 0.0)
                            Sx = Sx_L;
                        else
                            Sx = Sx_R;
                        if (v_Mid > 0.0)
                            Sy = Sy_L;
                        else
                            Sy = Sy_R;
                        */
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
                            src_scalar(i, j, k) = ((fw - fe) / dx[0] + (fs - fn) / (y * dx[1]) + source_Fij);
			    //src_scalar(i, j, k) = (-1.0*(u_Mid * Sx + v_Mid * Sy) + source_Fij);
                        }
                        else
                        {
                            src_scalar(i, j, k) = ((fw - fe) / dx[0] + (fs - fn) / dx[1] + source_Fij);
			    //src_scalar(i, j, k) = (-1.0*(u_Mid * Sx + v_Mid * Sy) + source_Fij);
                        }
                    }
                });
            }       
        }        
    }

} // namespace mycode
