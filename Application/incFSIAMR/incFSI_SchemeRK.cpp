#include "incFSI.H"
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <CFMask.H>

namespace mycode
{
    void incFSI::Scheme_RK()
    {
	/*
        int order = 4;
        int RKOrder_  = 1;
        //if (Iter % regrid_int == 0) order = 1;
        if(order == 1)
        {
            rk_a[0][0] = 1.0;
            rk_c[0] = 1.0;
            RKOrder_ = 1;
            if(N_IF > 0)
            {
                for (int lev = 0; lev <= finest_level; lev++)
                {
                    for (auto &&solid : interfaces[lev])
                    {
                        solid->rk_a[0][0] = 1.0;
                        solid->rk_c[0] = 1.0;
                    }
                }
            }
        }
        else if(order == 2)
        {
            rk_a[0][0] = 1.0;
            rk_a[1][0] = 0.5;
            rk_a[1][1] = 0.5;
            rk_c[0] = 1.0;
            rk_c[1] = 1.0;
            RKOrder_ = 2;
            if(N_IF > 0)
            {
                for (int lev = 0; lev <= finest_level; lev++)
                {
                    for (auto &&solid : interfaces[lev])
                    {   
                        solid->rk_a[0][0] = 1.0;
                        solid->rk_a[1][0] = 0.5;
                        solid->rk_a[1][1] = 0.5;
                        solid->rk_c[0] = 1.0;
                        solid->rk_c[1] = 1.0;
                    }
                }
            }
        }
        else if(order == 3)
        {
            rk_a[0][0] = 1.0/2.0;
            rk_a[1][0] = -1.0;
            rk_a[1][1] = 2.0;
            rk_a[2][0] = 1.0/6.0;
            rk_a[2][1] = 4.0/6.0;
            rk_a[2][2] = 1.0/6.0;

            rk_c[0] = 1.0/2.0;
            rk_c[1] = 1.0;
            rk_c[2] = 1.0;
            RKOrder_ = 3;
            if(N_IF > 0)
            {
                for (int lev = 0; lev <= finest_level; lev++)
                {       
                    for (auto &&solid : interfaces[lev])
                    {
                        solid->rk_a[0][0] = 1.0/2.0;
                        solid->rk_a[1][0] = -1.0;
                        solid->rk_a[1][1] = 2.0;
                        solid->rk_a[2][0] = 1.0/6.0;
                        solid->rk_a[2][1] = 4.0/6.0;
                        solid->rk_a[2][2] = 1.0/6.0;
                
                        solid->rk_c[0] = 1.0/2.0;
                        solid->rk_c[1] = 1.0;
                        solid->rk_c[2] = 1.0;
                    }
                }
            }
	    
            //rk_a[0][0] = 1.0/3.0;
            //rk_a[1][0] = -1.0;
            //rk_a[1][1] = 2.0;
            //rk_a[2][0] = 0.0;
            //rk_a[2][1] = 3.0/4.0;
            //rk_a[2][2] = 1.0/4.0;

            //rk_c[0] = 1.0/3.0;
            //rk_c[1] = 1.0;
            //rk_c[2] = 1.0;
            //RKOrder_ = 3;
            //if(N_IF > 0)
            //{
            //    for (int lev = 0; lev <= finest_level; lev++)
            //    {
            //        for (auto &&solid : interfaces[lev])
            //        {
            //            solid->rk_a[0][0] = 1.0/3.0;
            //            solid->rk_a[1][0] = -1.0;
            //            solid->rk_a[1][1] = 2.0;
            //            solid->rk_a[2][0] = 0.0;
            //            solid->rk_a[2][1] = 3.0/4.0;
            //            solid->rk_a[2][2] = 1.0/4.0;

            //            solid->rk_c[0] = 1.0/3.0;
            //            solid->rk_c[1] = 1.0;
            //            solid->rk_c[2] = 1.0;
            //        }
            //    }
            //}
            
        }
        else if(order == 4)
        {
            rk_a[0][0] = 1.0/2.0;
            rk_a[1][0] = 0.0;
            rk_a[1][1] = 1.0/2.0;
            rk_a[2][0] = 0.0;
            rk_a[2][1] = 0.0;
            rk_a[2][2] = 1.0;
            rk_a[3][0] = 1.0/6.0;
            rk_a[3][1] = 2.0/6.0;
            rk_a[3][2] = 2.0/6.0;
            rk_a[3][3] = 1.0/6.0;


            rk_c[0] = 1.0/2.0;
            rk_c[1] = 1.0/2.0;
            rk_c[2] = 1.0;
            rk_c[3] = 1.0; 
            RKOrder_ = 4;
            if(N_IF > 0)
            {
                for (int lev = 0; lev <= finest_level; lev++)
                {
                    for (auto &&solid : interfaces[lev])
                    {
                        solid->rk_a[0][0] = 1.0/2.0;
                        solid->rk_a[1][0] = 0.0;
                        solid->rk_a[1][1] = 1.0/2.0;
                        solid->rk_a[2][0] = 0.0;
                        solid->rk_a[2][1] = 0.0;
                        solid->rk_a[2][2] = 1.0;
                        solid->rk_a[3][0] = 1.0/6.0;
                        solid->rk_a[3][1] = 2.0/6.0;
                        solid->rk_a[3][2] = 2.0/6.0;
                        solid->rk_a[3][3] = 1.0/6.0;
                        
                        
                        solid->rk_c[0] = 1.0/2.0;
                        solid->rk_c[1] = 1.0/2.0;
                        solid->rk_c[2] = 1.0;
                        solid->rk_c[3] = 1.0;
                    }
                }
            }
        }
	*/
        //! Regrid
        //amrex::Print()<<"RK3"<<'\n';
        //amrex::Print()<<"Projection = "<<Projection<<"\n";
        //amrex::Print()<<"isAxisymmetric = "<<isAxisymmetric<<"\n";
        //if(Iter <= 658)
        Regrid();
        /*
        if(Iter == 658)
        {
            WriteFile();
            if(tecplot)
            {
                WriteFileTecplot();
                WriteFileTecplotW_Ghost();
            }
            //exit(1);
        }
        */
        //if(Iter >= 660) Mu = 0;
        /// copy new to old
        CopyNewToOld();
        for (int lev = 0; lev <= finest_level; lev++)
            for (auto &&solid : interfaces[lev])
                solid->copyFRK_Psi();

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
                            solid->TubeIdentification();
                            mask_[lev]->GhostCellIdentfication();
                            solid->TubeAdvectLevelSet_RK3(RKStage, xvel[lev], yvel[lev], dt); 
                            if (Iter % ReinitInt == 0 && RKOrder == 1)
                            {
                                if(use_algoim)
                                    solid->Reinit_algoim();
                                else
                                {
                                    if(lev == finest_level)
                                        solid->Reinit();
                                    else
                                        solid->Reinit2();
                                }
                            }
                            solid->copyRK1_Psi();

                        }
                        else if(RKStage == 2)
                        {
                            solid->TubeAdvectLevelSet_RK3(RKStage, xvel[lev], yvel[lev], dt);
                            
                            if (Iter % ReinitInt == 0 && RKOrder == 2)
                            {
                                if(use_algoim)
                                    solid->Reinit_algoim();
                                else
                                {
                                    if(lev == finest_level)
                                        solid->Reinit();
                                    else
                                        solid->Reinit2();
                                }
                            }
                            
                            solid->copyRK2_Psi();
                        }
                        else if(RKStage == 3)
                        {
                            solid->TubeAdvectLevelSet_RK3(RKStage, xvel[lev], yvel[lev], dt);
                            if (Iter % ReinitInt == 0 && RKOrder == 3)
                            {
                                if(use_algoim)
                                    solid->Reinit_algoim();
                                else
                                {
                                    if(lev == finest_level)
                                        solid->Reinit();
                                    else
                                        solid->Reinit2();
                                }
                            }
                            solid->copyRK3_Psi();
                        }
                        else if(RKStage == 4)
                        {
                            solid->TubeAdvectLevelSet_RK3(RKStage, xvel[lev], yvel[lev], dt);
                            if (Iter % ReinitInt == 0 && RKOrder == 4)
                            {
                                if(use_algoim)
                                    solid->Reinit_algoim();
                                else
                                {
                                    if(lev == finest_level)
                                        solid->Reinit();
                                    else
                                        solid->Reinit2();
                                }
                            }
                        }
                    }
                    for (int i = 0; i < N_IF; i++)
                    {
                        auto &&solid = interfaces[lev][i];
                        solid->TubeIdentification();
                        solid->Compute_Normal_Curvature();
                    }
                    //if(Iter == 8010 && RKStage == 4)
                    //{
                        //WriteFile();
                        //WriteFileTecplot();
                        //WriteInterface();
                        //exit(5);
                    //}

                    
                    for (int i = 0; i < N_IF; i++)
                    {
                        auto &&solid = interfaces[lev][i];
                        auto &&solid_f = interfaces[finest_level][i];
                        const amrex::Real *dx = geom[finest_level].CellSize();

                        solid->RemoveSpeckles();
                        if(solid->has_spekle == 1)
                        {
                            amrex::Print(-1)<<"Found speclke"<<'\n';
                            solid->Reinit();
                        }

                        if(solid_f->max_kappa_() > 1.0/(3.0*dx[0]))
                        {
                            amrex::Print(-1)<<"max_kappa = "<<solid_f->max_kappa_()<<" , 1.0/(3.0*dx[0]) = "<<1.0/(2.0*dx[0])<<'\n'<<"applying regularization"<<'\n';
                            solid->Regularization2(500);
                            if(lev == finest_level)
                               solid->Reinit();
                            else
                               solid->Reinit2();
                            //solid->RemoveSpeckles();
                            if(solid->has_spekle == 1)
                            {
                                amrex::Print(-1)<<"Found speclke"<<'\n';
                                //solid->Reinit();
                            }
                            solid->TubeIdentification();
                        }
                    }
                    //if(Iter == 8010 && RKStage == 4)
                    //{
                    //    WriteFile();
                    //    WriteFileTecplot();
                    //    //WriteInterface();
                    //    exit(6);
                    //}
                    
                    mask_[lev]->GhostCellIdentfication();
                }

                //if(Iter == 8010 && RKStage == 4)
                //{
                    //WriteFile();
                    //WriteFileTecplot();
                    //WriteInterface();
                    //exit(7);
                //}

                if(RKStage == RKOrder)
                {
	            for (auto &&solid : interfaces[finest_level])
                    {
                        solid->DetectUnderresolvedJetAlongAxis();
                        int count = 0;
                        if(solid->has_jet == 1)amrex::Print()<<"solid->has_jet = "<<solid->has_jet<<'\n';
                        while((solid->has_jet  == 1)&& (count < 10))
                        {
                            amrex::Print()<<"Removing jet cell by cell"<<'\n';
                            amrex::Print()<<"solid->has_jet = "<<solid->has_jet<<'\n';
                            solid->RemoveJetCellByCell();
                            solid->Reinit();
                            solid->Regularization2(500);
                            solid->Reinit();
                            solid->TubeIdentification();
                            mask_[finest_level]->GhostCellIdentfication();
                            solid->DetectUnderresolvedJetAlongAxis();
                            count++;
                        }
		    }
                }

		AverageDownInterface();
                ComputeBubbleVolume(false);
                //for (int lev = 0; lev <= finest_level; lev++)
                {
	            int lev = finest_level;
		    //amrex::Print()<<" lev = "<<lev<<'\n';
                    ComputeCutFaceVel(lev);
                    ComputeCutFacePressure(lev);
                }
		AverageDown();
            }
	    if(Iter == 7982 && RKStage == 4)
	    {
                //WriteFile();
                //WriteFileTecplot();
                //WriteInterface();
                //exit(1);
	    }
            
            /// Prediction
            for (int lev = 0; lev <= finest_level; lev++)
            {
                //ComputeIntermediateVelocity(lev);
                ComputeIntermediateVelocity_RK(lev, RKStage);
            }

            if(Iter == 7982 && RKStage == 3)
            {
                //WriteFile();
                //WriteFileTecplot();
                //WriteInterface();
                //exit(3);
            }
            
            ApplyBC();
            
            AverageDown();
            
            AMGPressurePoisson_RK(RKStage);

            FillPatchAroundBox();
            //ApplyBC();
            
            for (int lev = 0; lev <= finest_level; lev++)
            {
                if(N_IF > 0 & lev == finest_level)
                { 
                    mask_[lev]->FillInGhost(Pressure[lev], Pintbcs);
                }
                ComputeFinalVelocityField_RK(lev, RKStage);
            }

            //if(Iter == 8010 && RKStage == 3)
            //{
                //WriteFile();
                //WriteFileTecplot();
                //WriteInterface();
                //exit(4);
            //}
            
            //! Avg down face centered velocity
            ApplyBC();
            
            AverageDown();

            //FillPatchAroundBox();
            if(RKStage == 1) 
                CopyNewToRK1();
            if(RKStage == 2)
                CopyNewToRK2();
            if(RKStage == 3)
                CopyNewToRK3();
    
            for (int lev = 0; lev <= finest_level; lev++)
            {
                if(N_IF == 0)ComputeCollocatedVelocityField(lev);
            }
            //if(Iter == 8010 && RKStage == 3)
            //{
                //WriteFile();
                //WriteFileTecplot();
                //WriteInterface();
                //exit(4);
            //}
        }
        Derive();
        if(PhaseField)CheckConservationPhaseField();
        if(DamageModel)
        {
            for (int lev = 0; lev <= finest_level; lev++)
            {
                ComputeDeformationGradientTensor2D(lev);
            }
        }
    }

    void incFSI::CopyNewToRK1()
    {
        //amrex::Print()<<"CopyNewToRK1()"<<'\n';
        for (int lev = 0; lev <= finest_level; lev++)
        {
            amrex::MultiFab::Copy(RK1_xvel[lev], xvel[lev], 0, 0, 1, Nghost);
            amrex::MultiFab::Copy(RK1_yvel[lev], yvel[lev], 0, 0, 1, Nghost);
            amrex::MultiFab::Copy(RK1_p[lev], Pressure[lev], 0, 0, 1, Nghost);
            if(TempField) amrex::MultiFab::Copy(RK1_Theta[lev], Theta[lev], 0, 0, 1, Nghost);
            if(PhaseField) amrex::MultiFab::Copy(RK1_Phi[lev], Phi[lev], 0, 0, 1, Nghost);
            if(DamageModel)
            {
                for(int iscalar = 0; iscalar < nscalar; iscalar++)
                {
                    if(iscalar == eps_max_num)
                    continue;
                    amrex::MultiFab::Copy(RK1_Scalars[lev][iscalar], Scalars[lev][iscalar], 0, 0, 1, Nghost);
                }
            }
        }
    }

    void incFSI::CopyNewToRK2()
    {
        //amrex::Print()<<"CopyNewToRK2()"<<'\n';
        for (int lev = 0; lev <= finest_level; lev++)
        {
            amrex::MultiFab::Copy(RK2_xvel[lev], xvel[lev], 0, 0, 1, Nghost);
            amrex::MultiFab::Copy(RK2_yvel[lev], yvel[lev], 0, 0, 1, Nghost);
            amrex::MultiFab::Copy(RK2_p[lev], Pressure[lev], 0, 0, 1, Nghost);
            if(TempField) amrex::MultiFab::Copy(RK2_Theta[lev], Theta[lev], 0, 0, 1, Nghost);
            if(PhaseField) amrex::MultiFab::Copy(RK2_Phi[lev], Phi[lev], 0, 0, 1, Nghost);
            if(DamageModel)
            {   
                for(int iscalar = 0; iscalar < nscalar; iscalar++)
                {   
                    if(iscalar == eps_max_num)
                    continue;
                    amrex::MultiFab::Copy(RK2_Scalars[lev][iscalar], Scalars[lev][iscalar], 0, 0, 1, Nghost);
                }
            }
        }
    }

    void incFSI::CopyNewToRK3()
    {
        //amrex::Print()<<"CopyNewToRK2()"<<'\n';
        for (int lev = 0; lev <= finest_level; lev++)
        {
            amrex::MultiFab::Copy(RK3_xvel[lev], xvel[lev], 0, 0, 1, Nghost);
            amrex::MultiFab::Copy(RK3_yvel[lev], yvel[lev], 0, 0, 1, Nghost);
            amrex::MultiFab::Copy(RK3_p[lev], Pressure[lev], 0, 0, 1, Nghost);
            if(TempField) amrex::MultiFab::Copy(RK3_Theta[lev], Theta[lev], 0, 0, 1, Nghost);
            if(PhaseField) amrex::MultiFab::Copy(RK3_Phi[lev], Phi[lev], 0, 0, 1, Nghost);
            if(DamageModel)
            {   
                for(int iscalar = 0; iscalar < nscalar; iscalar++)
                {   
                    if(iscalar == eps_max_num)
                    continue;
                    amrex::MultiFab::Copy(RK3_Scalars[lev][iscalar], Scalars[lev][iscalar], 0, 0, 1, Nghost);
                }
            }
        }
     }

} /*End namespace mycode */

