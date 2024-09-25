#include "incFSI.H"
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <CFMask.H>

namespace mycode
{
    void incFSI::Scheme_RK2()
    {
        Regrid();
        CopyFRKVariables();
	 
        for (int RKStage = 1; RKStage <= RKOrder; RKStage++)
        {
            amrex::Print()<<"RKStage = "<<RKStage<<'\n';
            
            if(N_IF > 0)
            {
                //for (int lev = 0; lev <= finest_level; lev++)
                {
		            int lev = finest_level;
                    for (auto &&solid : interfaces[lev])
                    {
                        if(RKStage == 1)
                        {
			                solid->setVolume_prev();
                            solid->copyFRK_Psi();
                            solid->copyRK1_Psi();
                            solid->TubeIdentification();
                            mask_[lev]->GhostCellIdentfication();
                            solid->TubeAdvectLevelSet_RK2(xvel[lev], yvel[lev], dt); 
                            if (Iter % ReinitInt == 0)
                            {
                                if(use_algoim)
                                    solid->Reinit_algoim();
                                else
				                {
                                    if(Reinit_Method == 1)
                                        solid->Reinit();
                                    else if(Reinit_Method == 2)
				                        solid->Reinit2();
				                }
                            }
                        }
                        else if(RKStage == 2)
                        {
                            solid->copyRK1_Psi();
                            solid->TubeAdvectLevelSet_RK2(xvel[lev], yvel[lev], dt);
                            solid->RK2Avg_Psi();

			                if (LS_reg_iter != 0 && Iter % LS_reg_iter == 0)
		                    {
                                //solid->Regularization2(500);
                                amrex::Print()<<"Applying interface smoothing"<<'\n';
                                solid->Regularization();
                                if(Reinit_Method == 1)
                                    solid->Reinit();
                                else if(Reinit_Method == 2)
                                    solid->Reinit2();
                            }
                            if (Iter % ReinitInt == 0)
                            {
                                if(use_algoim)
                                    solid->Reinit_algoim();
                                else
                                {
                                     if(Reinit_Method == 1)
                                         solid->Reinit();
                                     else if(Reinit_Method == 2)
                                         solid->Reinit2();
                                }
                            }
                        }
                    }
                    //for (auto &&solid : interfaces[lev])
                    for (int i = 0; i < N_IF; i++)
                    {
            			auto &&solid = interfaces[lev][i];
                        solid->TubeIdentification();
			            solid->Compute_Normal_Curvature();
                    }

                    for (int i = 0; i < N_IF; i++)
		            {
                        auto &&solid = interfaces[lev][i];
                        auto &&solid_f = interfaces[finest_level][i];
                        const amrex::Real *dx = geom[finest_level].CellSize();
                         
                        while(solid_f->max_kappa_() > 1.0/(dx[0]))
                        {
                            amrex::Print(-1)<<"max_kappa = "<<solid_f->max_kappa_()<<" , 1.0/(1.0*dx[0]) = "<<1.0/(dx[0])<<'\n'<<"applying regularization"<<'\n';
                            solid->Regularization2(100);
                            solid->Reinit2();
                            
                            //solid->RemoveSpeckles();
                            //if(solid->has_spekle == 1)
                            //{
                            //    amrex::Print(-1)<<"Found speclke"<<'\n';
                            //    if(Reinit_Method == 1)
                            //        solid->Reinit();
                            //    else if(Reinit_Method == 2)
                            //        solid->Reinit2();
                            //}
                            
                            solid->TubeIdentification();
                        }
		        	
                    }
                    
                    mask_[lev]->GhostCellIdentfication();
                }

		        {
                    auto &&solid = interfaces[finest_level][0];
                    solid->DetectUnderresolvedJetAlongAxis();
                    int count = 0;
		            while((solid->has_jet  == 1)&& (count < 10))
                    {
                        amrex::Print()<<"Removing jet cell by cell"<<'\n';
                        amrex::Print()<<"solid->has_jet = "<<solid->has_jet<<'\n';
                        solid->RemoveJetCellByCell();
			            solid->Regularization();
                        solid->Reinit2();
                        solid->Regularization2(500);//remove high curvature underresolved films
                        solid->Reinit2();
                        solid->TubeIdentification();
                        mask_[finest_level]->GhostCellIdentfication();
			            solid->Compute_Normal_Curvature();
                        solid->DetectUnderresolvedJetAlongAxis();
                        count++;
                    }
                }
		
                
                AverageDownInterface();
                //amrex::Print()<<"After AverageDownInterface"<<'\n';
                
                ComputeBubbleVolume(false);
                //for (int lev = 0; lev <= finest_level; lev++)
                {
		             int lev = finest_level;
                    ComputeCutFaceVel(lev);
                    ComputeCutFacePressure(lev);
                }
		        AverageDown();
	        	
                //if(Iter == 23001 && RKStage == 2)
                //{
                    //amrex::Print()<<"Plotting"<<'\n';
                    
                    //auto &&solid = interfaces[finest_level][0];
                    //solid->DetectUnderresolvedJetAlongAxis();
                    //amrex::Print()<<"jet = "<<solid->has_jet<<'\n';
                    //int count = 0;
                    //do
                    //{
                    //	amrex::Print()<<"Removing jet cell by cell"<<'\n';
                    //	amrex::Print()<<"solid->has_jet = "<<solid->has_jet<<'\n';
                        //    solid->RemoveJetCellByCell();
                    //	solid->Reinit();
                    //	solid->Regularization2(500);
                    //	solid->Reinit();
                    //	solid->TubeIdentification();
                    //	mask_[finest_level]->GhostCellIdentfication();
                    //	solid->DetectUnderresolvedJetAlongAxis();
                    //	count++;
                    //    }while(solid->has_jet || (count > 10));
                    //solid->Regularization2();
                    
                    //auto &&solid = interfaces[finest_level][0];
                    //solid->Regularization();
                    //solid->Reinit2();
                    //WriteFile();
                    //WriteFileTecplot();
                    //ComputeCutFaceVel(finest_level);
                    //ComputeCutFacePressure(finest_level);
                    //WriteInterface();
                    //exit(1);      
                //}
            }
             
            /// copy new to old
            CopyNewToOld();
            //if(RKStage == 1) CopyFRKVariables();
            
            /// Prediction
            for (int lev = 0; lev <= finest_level; lev++)
                ComputeIntermediateVelocity_RK2(lev);
            if(DamageModel)
            {
                for (int lev = 0; lev <= finest_level; lev++)
                {
                    AdvectScalars_RK2(lev);
                }
            }
            
            ApplyBC();
            
            AverageDown();

            AMGPressurePoisson();

            FillPatchAroundBox();
            
            for (int lev = 0; lev <= finest_level; lev++)
            {
                if(N_IF > 0)
                { 
                    mask_[lev]->FillInGhost(Pressure[lev], Pintbcs);
                }
                ComputeFinalVelocityField(lev);
            }

            ApplyBC();
            
            AverageDown();

            if(RKStage == 2)
            {
                RK2Avg();
                ApplyBC();
            }
            
            for (int lev = 0; lev <= finest_level; lev++)
            {
                ComputeCollocatedVelocityField(lev);
            }
	    }

        Derive();
        if(DamageModel)
        {
            for (int lev = 0; lev <= finest_level; lev++)
            {
                ComputeDeformationGradientTensor2D(lev);
            }
        }
    }
    
    void incFSI::Scheme_TVDRK2()
    {
        //! Regrid
        amrex::Print()<<"TVDRK2"<<'\n';
        Regrid();
        //int RKStage = 1; 
        for (int RKStage = 1; RKStage <= RKOrder; RKStage++)
        {
            amrex::Print()<<"RKStage = "<<RKStage<<'\n';
            if(N_IF > 0)
            {
                for (int lev = 0; lev <= finest_level; lev++)
                {
                    for (auto &&solid : interfaces[lev])
                    {
                        if(RKStage == 1)
                        {
                            solid->copyFRK_Psi();
                            solid->TubeIdentification();
                            mask_[lev]->GhostCellIdentfication();
                            solid->TubeAdvectLevelSet(xvel[lev], yvel[lev], dt);                        
                            if (Iter % ReinitInt == 0)
                            {
                                if(use_algoim)
                                    solid->Reinit_algoim();
                                else
                                    solid->Reinit2();
                            }

                        }
                        else if(RKStage == 2)
                        {
                            solid->TubeAdvectLevelSet(xvel[lev], yvel[lev], dt);
                            solid->TVDRK2Avg_Psi();
                            if (Iter % ReinitInt == 0)
                            {
                                if(use_algoim)
                                    solid->Reinit_algoim();
                                else
                                    solid->Reinit2();
                            }
                            
                        }
                    }
                    for (auto &&solid : interfaces[lev])
                    {
                        solid->TubeIdentification();
                    }
                    mask_[lev]->GhostCellIdentfication();
                }
                ComputeBubbleVolume(false);
                for (int lev = 0; lev <= finest_level; lev++)
                {
                    ComputeCutFaceVel(lev);
                    ComputeCutFacePressure(lev);
                }
            }
            
            /// copy new to old
            CopyNewToOld();
            if(RKStage == 1) CopyFRKVariables();
            
            /// Prediction
            for (int lev = 0; lev <= finest_level; lev++)
                ComputeIntermediateVelocity(lev);
            
            ApplyBC();
            
            AverageDown();
            
            AMGPressurePoisson();
            
            for (int lev = 0; lev <= finest_level; lev++)
            {
                if(N_IF > 0)
                { 
                    mask_[lev]->FillInGhost(Pressure[lev], Pintbcs);
                }
                ComputeFinalVelocityField(lev);
            }
            
            //! Avg down face centered velocity
            ApplyBC();
            
            AverageDown();
    
            if(RKStage == 2)
            {
                TVD_RK2Avg();
                ApplyBC();
            }
            
            for (int lev = 0; lev <= finest_level; lev++)
            {
                if(N_IF == 0)ComputeCollocatedVelocityField(lev);
            }
        }
        Derive();
    }
    
    void incFSI::CopyFRKVariables()
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            amrex::MultiFab::Copy(FRK_xvel[lev], xvel[lev], 0, 0, 1, Nghost);
            amrex::MultiFab::Copy(FRK_yvel[lev], yvel[lev], 0, 0, 1, Nghost);
            amrex::MultiFab::Copy(FRK_p[lev], Pressure[lev], 0, 0, 1, Nghost);
            if(TempField)amrex::MultiFab::Copy(FRK_Theta[lev], Theta[lev], 0, 0, 1, Nghost);
            if(PhaseField)amrex::MultiFab::Copy(FRK_Phi[lev], Phi[lev], 0, 0, 1, Nghost);
            if(DamageModel)
            {
                for(int iscalar = 0; iscalar < nscalar; iscalar++)
                {
		    amrex::MultiFab::Copy(FRK_Scalars[lev][iscalar], Scalars[lev][iscalar], 0, 0, 1, Nghost);
                }
            }

            /*if(N_IF > 0)
            {
                for (auto &&solid : interfaces[lev])
                {
                    solid->copyFRK_Psi(); 
                }
            }*/
        }
    }
    
    void incFSI::TVD_RK2Avg()
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            for (amrex::MFIter mfi(Pressure[lev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();
    
                amrex::Array4<amrex::Real> const &p = Pressure[lev].array(mfi);
                amrex::Array4<amrex::Real> const &rk_p = FRK_p[lev].array(mfi);
    
                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    p(i, j, k) = 0.5 * (rk_p(i, j, k) + p(i, j, k));
                });
            }
    
            for (amrex::MFIter mfi(xvel[lev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();
    
                amrex::Array4<amrex::Real> const &u = xvel[lev].array(mfi);
                amrex::Array4<amrex::Real> const &rk_u = FRK_xvel[lev].array(mfi);
    
                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    u(i, j, k) = 0.5 * (rk_u(i, j, k) + u(i, j, k));
                });
            }
    
            for (amrex::MFIter mfi(yvel[lev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();
    
                amrex::Array4<amrex::Real> const &v = yvel[lev].array(mfi);
                amrex::Array4<amrex::Real> const &rk_v = FRK_yvel[lev].array(mfi);
    
                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    v(i, j, k) = 0.5 * (rk_v(i, j, k) + v(i, j, k));
                });
            }

            if(TempField)
            {
                for (amrex::MFIter mfi(Theta[lev]); mfi.isValid(); ++mfi)
                {
                    const amrex::Box &bx = mfi.validbox();

                    amrex::Array4<amrex::Real> const &T = Theta[lev].array(mfi);
                    amrex::Array4<amrex::Real> const &FRK_T = FRK_Theta[lev].array(mfi);

                    amrex::ParallelFor(bx,
                    [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        T(i, j, k) = 0.5 * (FRK_T(i, j, k) + T(i, j, k));
                    });
                }
            }

            if(PhaseField)
            {
                for (amrex::MFIter mfi(Phi[lev]); mfi.isValid(); ++mfi)
                {
                    const amrex::Box &bx = mfi.validbox();

                    amrex::Array4<amrex::Real> const &phi = Phi[lev].array(mfi);
                    amrex::Array4<amrex::Real> const &FRK_phi = FRK_Phi[lev].array(mfi);

                    amrex::ParallelFor(bx,
                    [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        phi(i, j, k) = 0.5 * (FRK_phi(i, j, k) + phi(i, j, k));
                    });
                }
            }

            if(DamageModel)
            {
                for(int iscalar = 0; iscalar < nscalar; iscalar++)
                {
                    if(iscalar == eps_max_num)
                        continue;
                    for (amrex::MFIter mfi(Scalars[lev][iscalar]); mfi.isValid(); ++mfi)
                    {
                        const amrex::Box &bx = mfi.validbox();
                        amrex::Array4<amrex::Real> const &scalar = Scalars[lev][iscalar].array(mfi);
                        amrex::Array4<amrex::Real> const &FRK_scalar = FRK_Scalars[lev][iscalar].array(mfi);
                        amrex::ParallelFor(bx,
                        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                        {
                            scalar(i, j, k) = 0.5 * (FRK_scalar(i, j, k) + scalar(i, j, k));
                        });
                    }
		}
            }
        }
    }
    
    void incFSI::RK2Avg()
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            for (amrex::MFIter mfi(Pressure[lev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();

                amrex::Array4<amrex::Real> const &p = Pressure[lev].array(mfi);
                amrex::Array4<amrex::Real> const &pold = Pressure_old[lev].array(mfi);

                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {   
                    p(i, j, k) = 0.5 * (pold(i, j, k) + p(i, j, k));
                });
            }

            for (amrex::MFIter mfi(xvel[lev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();

                amrex::Array4<amrex::Real> const &u = xvel[lev].array(mfi);
                amrex::Array4<amrex::Real> const &uold = xvel_old[lev].array(mfi);

                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {   
                    u(i, j, k) = 0.5 * (uold(i, j, k) + u(i, j, k));
                });
            }

            for (amrex::MFIter mfi(yvel[lev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();

                amrex::Array4<amrex::Real> const &v = yvel[lev].array(mfi);
                amrex::Array4<amrex::Real> const &vold = yvel_old[lev].array(mfi);

                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {   
                    v(i, j, k) = 0.5 * (vold(i, j, k) + v(i, j, k));
                });
            }

            if(TempField)
            {
                for (amrex::MFIter mfi(Theta[lev]); mfi.isValid(); ++mfi)
                {
                    const amrex::Box &bx = mfi.validbox();

                    amrex::Array4<amrex::Real> const &T = Theta[lev].array(mfi);
                    amrex::Array4<amrex::Real> const &T_old = Theta_old[lev].array(mfi);

                    amrex::ParallelFor(bx,
                    [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {   
                        T(i, j, k) = 0.5 * (T_old(i, j, k) + T(i, j, k));
                    });
                }
            }

            if(PhaseField)
            {
                for (amrex::MFIter mfi(Phi[lev]); mfi.isValid(); ++mfi)
                {
                    const amrex::Box &bx = mfi.validbox();

                    amrex::Array4<amrex::Real> const &phi = Phi[lev].array(mfi);
                    amrex::Array4<amrex::Real> const &phi_old = Phi_old[lev].array(mfi);

                    amrex::ParallelFor(bx,
                    [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        phi(i, j, k) = 0.5 * (phi_old(i, j, k) + phi(i, j, k));
                    });
                }
            }
            if(DamageModel)
            {
                for(int iscalar = 0; iscalar < nscalar; iscalar++)
                {
		    if(iscalar == eps_max_num)
			continue;
                    for (amrex::MFIter mfi(Scalars[lev][iscalar]); mfi.isValid(); ++mfi)
                    {
                        const amrex::Box &bx = mfi.validbox();

                        amrex::Array4<amrex::Real> const &scalar = Scalars[lev][iscalar].array(mfi);
                        amrex::Array4<amrex::Real> const &scalar_old = Scalars_old[lev][iscalar].array(mfi);
                        amrex::ParallelFor(bx,
                        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                        {
                            scalar(i, j, k) = 0.5 * (scalar_old(i, j, k) + scalar(i, j, k));
                        });
                    }
                }
            }
        }
    }


    void incFSI::RK2Avg_Psi(amrex::MultiFab& Phi_old, const amrex::MultiFab& Psi)
    {
    }
    
    void incFSI::CopyFRKPsi(amrex::MultiFab& Psi_old, const amrex::MultiFab& Psi)
    {
        if(N_IF > 0)
        {
            for (int lev = 0; lev <= finest_level; lev++)
            {
                for (auto &&solid : interfaces[lev])
                {
                    solid->copyFRK_Psi();
                }
            }
        }
    }



} /*End namespace mycode */

