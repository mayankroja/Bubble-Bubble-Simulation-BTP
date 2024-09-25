#include "incFSI.H"
#include <Mask.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MultiFabUtil.H>
#include <CFMask.H>

namespace mycode
{

    void incFSI::MakeInterface
    (
        int lev,  
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm
    )
    {
        mycode::amrexMesh mesh(ba, geom[lev], dm);
        if (N_IF > 0)
        {
            interfaces[lev].resize(N_IF);
            for (int i = 0; i < N_IF; i++)
            {
                if (IF_types[i] == "Cylinder" || IF_types[i] == "cylinder")
                {   
                    interfaces[lev][i] = mycode::MakeCylinder(mesh, IF_names[i]);
                }
                else if (IF_types[i] == "Ellipse" || IF_types[i] == "ellipse")
                {   
                    interfaces[lev][i] = mycode::MakeEllipse(mesh, IF_names[i]);
                }
                else if (IF_types[i] == "Foil" || IF_types[i] == "foil")
                {   
                    interfaces[lev][i] = mycode::MakeFoil(mesh, IF_names[i]);
                }
                else if (IF_types[i] == "Bubble" || IF_types[i] == "bubble")
                {   
                    interfaces[lev][i] = mycode::MakeBubble(mesh, IF_names[i]);
                }
                else
                {   
                    amrex::Abort("interface type not supported");
                }
            }      
          
            for (auto &&solid : interfaces[lev])
            {   
                if(!Restart)solid->PrescribeLevelSetMotion(0.0);//Creates levelset field at lev
                //solid->ComputeVolume(true);
                //amrex::Print()<<"Volume = "<<solid->Volume()<<'\n';
            }  
        }
    
        //Average down interface to lower lev
        if (lev > 0)
        {
            for (int i = 0; i < N_IF; i++)
            {		
                for(int ilev = lev ; ilev > 0;ilev--)
                {
                    auto &&solid_f = interfaces[ilev][i]; 
                    amrex::MultiFab &psi_f = solid_f->Psi();
    
                    auto &&solid_c = interfaces[ilev - 1][i];
                    amrex::MultiFab &psi_c = solid_c->Psi();
    
                    amrex::average_down(psi_f, psi_c, 0 , 1, ref_ratio[ilev - 1][0]);
                }
            }
        }
        
        /*//Compute bubble Volume from Multi level AMR mesh
        amrex::Vector<CFMask> cfmask_(finest_level + 1);
        for (int ilev = 0; ilev <= finest_level; ilev++)
        {
            cfmask_[ilev].define(ilev, geom, grids, dmap);
        }
        for (int i = 0; i < N_IF; i++)
        {
            for(int ilev = 0; ilev <= finest_level; ilev++)
            {
                auto &&solid = interfaces[ilev][i];
                solid->ComputeVolume(true, cfmask_[ilev].Mask());
                amrex::Print()<<"Volume of IF "<<i<<" at lev "<<ilev<<" = "<<solid->Volume()<<'\n';
            }
        }*/


      
        //! NOTE : if the interfaces vector empty this will simply give mask = 1 everywhere
    
	//if(lev == finest_level)
        //{
            mask_[lev].reset(new Mask);
    
            mask_[lev]->define(&mesh, &interfaces[lev]);
    
            //! compute mask and interface related qty
            if(!Restart)mask_[lev]->GhostCellIdentfication();

            ComputeIntProp(lev);
	//}

        if(RKOrder == 2)
        {
            for (auto &&solid : interfaces[lev])
            {
                solid->copyFRK_Psi();
            }
        }
            
    
    }  

    void incFSI::MakeNewLevelFromCoarse_interface
    (
        int lev,  
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm
    )
    {
        amrex::Print()<<"MakeNewLevelFromCoarse interface: lev = "<<lev<<'\n';
        mycode::amrexMesh mesh(ba, geom[lev], dm);
        amrex::MultiFab tmp_Psi(ba, dm, 1, Nghost);

        if (lev > 0)
        {
            //interfaces[lev].resize(N_IF);
            for (int i = 0; i < N_IF; i++)
            {
                if (IF_types[i] == "Cylinder" || IF_types[i] == "cylinder")
                {   
                    interfaces[lev][i] = mycode::MakeCylinder(mesh, IF_names[i]);
                }
                else if (IF_types[i] == "Ellipse" || IF_types[i] == "ellipse")
                {   
                    interfaces[lev][i] = mycode::MakeEllipse(mesh, IF_names[i]);
                }
                else if (IF_types[i] == "Foil" || IF_types[i] == "foil")
                {   
                    interfaces[lev][i] = mycode::MakeFoil(mesh, IF_names[i]);
                }
                else if (IF_types[i] == "Bubble" || IF_types[i] == "bubble")
                {   
                    interfaces[lev][i] = mycode::MakeBubble(mesh, IF_names[i]);
                }
                else
                {   
                    amrex::Abort("interface type not supported");
                }
            }

            for (int i = 0; i < N_IF; i++)
            {
                tmp_Psi.setVal(0.0);

                auto &&solid_c = interfaces[lev - 1][i];
                amrex::MultiFab &psi_c_old = solid_c->Psi();
   
                FillCoarsePatch(lev, time, tmp_Psi, psi_c_old, 0, 1);
                auto &&solid_f = interfaces[lev][i];
                solid_f->MakeInterfaceFromCoarse(mesh,tmp_Psi);
            }
            /*for (auto &&solid : interfaces[lev])
            {
                solid->ComputeVolume(true);
            }*/

        }


        //! NOTE : if the interfaces vector empty this will simply give mask = 1 everywhere
        //if(lev == finest_level)
        //{
            mask_[lev].reset(new Mask);

            mask_[lev]->define(&mesh, &interfaces[lev]);

            //! compute mask and interface related qty
            mask_[lev]->GhostCellIdentfication();
        //}
    }

    
    void incFSI::RemakeInterface
    (
        int lev,  
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm
    )
    {
        //amrex::Print()<<"Remake interface: lev = "<<lev<<'\n';
        mycode::amrexMesh mesh(ba, geom[lev], dm);
        amrex::MultiFab tmp_Psi(ba, dm, 1, Nghost);
        amrex::Vector<amrex::Real> ctime(2);
        amrex::Vector<amrex::Real> ftime(2);

        if (lev > 0)
        {
            ctime[0] = t_old;
            ctime[1] = t_new;

            ftime[0] = t_old;
            ftime[1] = t_new;

            for (int i = 0; i < N_IF; i++)
            {
                tmp_Psi.setVal(0.0);

                auto &&solid_f = interfaces[lev][i];
                amrex::MultiFab &psi_f_old = solid_f->Psi();

                auto &&solid_c = interfaces[lev - 1][i];
                amrex::MultiFab &psi_c_old = solid_c->Psi();

                /*for (amrex::MFIter mfi(psi_c_old); mfi.isValid(); ++mfi)
                {
                    const amrex::Box &bx = mfi.growntilebox();
                    const amrex::Box &bx_valid = mfi.validbox();
                    amrex::Array4<amrex::Real> const &f = psi_c_old.array(mfi);
                    //amrex::Print(-1)<<"bx = "<<bx<<'\n';
                    //amrex::Print(-1)<<"bx_valid = "<<bx_valid<<'\n';
                    //amrex::Print(-1)<<"dx = "<<dx[0]<<'\n';


                    /// copy Psi in phi
                    amrex::ParallelFor(bx,
                    [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        int ilocal = i - bx.smallEnd(0);
                        int jlocal = j - bx.smallEnd(1);
                        phi(ilocal, jlocal) = f(i, j, k);
                    });
                }*/

                FillPatch(lev, time, tmp_Psi, {&psi_c_old}, {ctime[0]}, {&psi_f_old}, {ftime[0]}, 0, 1);
                solid_f->ClearInterfaceData();
                solid_f->RemakeInterface(mesh,tmp_Psi);
                //solid_f->displayMesh();
            }
        }

        //if(lev == finest_level)
	//{
            //! NOTE : if the interfaces vector empty this will simply give mask = 1 everywhere
    
            mask_[lev].reset(new Mask);
    
            mask_[lev]->define(&mesh, &interfaces[lev]);
    
            //! compute mask and interface related qty
            mask_[lev]->GhostCellIdentfication();
        //}
    }  
    
    void incFSI::AverageDownInterface()
    {
        //amrex::Print()<<"Avg down interface: lev = "<<lev<<'\n';
        for (int i = 0; i < N_IF; i++)
        {
            for(int ilev = finest_level ; ilev > 0;ilev--)
            {
                auto &&solid_f = interfaces[ilev][i];
                amrex::MultiFab &psi_f = solid_f->Psi();

                auto &&solid_c = interfaces[ilev - 1][i];
                amrex::MultiFab &psi_c = solid_c->Psi();

                amrex::average_down(psi_f, psi_c, 0 , 1, ref_ratio[ilev - 1][0]);
		solid_c->Reinit2();
		solid_c->TubeIdentification();
		mask_[ilev - 1]->GhostCellIdentfication();
            }
        }
    }        
     
    void incFSI::ComputeCutFacePressure(int lev)
    {
   //Filling ghost values for PhaseField and Fijs need to be filled  before pressure because the phasefiled and damage at the interface is required for the calcualtion of pressure at the interface  
        if(PhaseField) 
	{
            mask_[lev]->FillInGhostPhaseField(Phi[lev], Phiintbcs);
	    mask_[lev]->ExtendPhaseFieldLSQ(Phi[lev]);
	}
	if(DamageModel) 
        {
	   mask_[lev]->FillInRefConfig(Scalars[lev][0], Scalars[lev][1], Scalars[lev][8]);
	   if(!advect_ref_cond)
           {
	      //for(int iscalar = 2 ; iscalar <= 6; iscalar++)
	      //     mask_[lev]->FillInGhostFij(iscalar,Scalars[lev][iscalar]);
	      if(PhaseField)
	      {
                  mask_[lev]->FillInGhostFij(Scalars[lev][2],
                                             Scalars[lev][3],
                                             Scalars[lev][4],
                                             Scalars[lev][5],
                                             Scalars[lev][6],
					     Phi[lev]);
              }
              else
              {
	          mask_[lev]->FillInGhostFij(Scalars[lev][2],
	      	                             Scalars[lev][3],
	      	         		     Scalars[lev][4],
	      			             Scalars[lev][5],
				             Scalars[lev][6]);
              }

	      mask_[lev]->ExtendFijLSQ(Scalars[lev][2],
                                       Scalars[lev][3],
                                       Scalars[lev][4],
                                       Scalars[lev][5],
                                       Scalars[lev][6]);
	      //mask_[lev]->FillInGhostFij(Scalars[lev][6]);
	   }
	}

        // Set interface pressure
        mask_[lev]->SetInterfacePressureVel();
        //! compute interpolation wts
        mask_[lev]->FillInGhost(Pressure[lev], Pintbcs);
        if(TempField) mask_[lev]->FillInGhostTheta(Theta[lev], Tintbcs);


        //ComputeAvgIntProp(lev);
    
    }
    
    void incFSI::ComputeCutFaceVel(int lev)
    {
        ComputeCollocatedVelocityField(lev);
        mask_[lev]->FillInVelocityComponents(U[lev]);
        mask_[lev]->ExtendVelocityLSQ(U[lev]);
        ComputeStaggeredVelocityField(lev);
    }
    
    void incFSI::ComputeIntProp(int lev)
    {
        const amrex::iMultiFab &PMask = mask_[lev]->getPMask();
        for (auto &&solid : interfaces[lev])
        {
            solid->ComputeAvgIntVel(xvel[lev],yvel[lev],PMask);
            IF_props[0]->Int_p = solid->AvgIntPres();
            if(Iter == 0) IF_props[0]->Int_p = solid->getP_interface0();
            IF_props[0]->Int_R = solid->AvgRadius();
            IF_props[0]->Int_R_dot = solid->AvgIntVel();
            IF_props[0]->Int_xcp = solid->Xcp();
            IF_props[0]->Int_xcp = solid->Ycp();
	    solid->SetTime(Time);
        }
    }

    void incFSI::ComputeBubbleVolume(bool init)
    {
        amrex::Vector<CFMask> cfmask_(finest_level + 1);
        for (int ilev = 0; ilev <= finest_level; ilev++)
        {
            cfmask_[ilev].define(ilev, geom, grids, dmap);
        }
        for (int i = 0; i < N_IF; i++)
        {
            amrex::Real vol_i = 0.0;
            amrex::Real vol_diff_i = 0.0;
            amrex::Real surf_area_i = 0.0;
            int resolution = 0;
            for(int ilev = 0; ilev <= finest_level; ilev++)
            {
                auto &&solid = interfaces[ilev][i];
                solid->ComputeVolume(true, cfmask_[ilev].Mask());
                vol_i += solid->Part_Volume();
                vol_diff_i += solid->Part_Vol_Diff();
                surf_area_i += solid->part_surface_area();
                resolution += solid->getPartResolution();
                // amrex::Print()<<"Volume of IF "<<i<<" at lev "<<ilev<<" = "<<solid->Part_Volume()<<'\n';
                // amrex::Print()<<"Surface area of IF "<<i<<" at lev "<<ilev<<" = "<<solid->part_surface_area()<<'\n';
                // amrex::Print()<<"Vol_Diff of IF "<<i<<" at lev "<<ilev<<" = "<<solid->Part_Vol_Diff()<<'\n';
            }
            //amrex::Print()<<"vol_i = "<<vol_i<<" , vol_diff = "<<vol_diff_i<<'\n';
            for(int ilev = 0; ilev <= finest_level; ilev++)
            {
                auto &&solid = interfaces[ilev][i];
                solid->SetVolume(vol_i);
                if(init)solid->SetVolume_0(vol_i);

                solid->setResolution(resolution);
                if(init)solid->setResolution0(resolution);   

                amrex::Real sphericity = 1.0;
                if(surf_area_i > 0)sphericity = std::pow(M_PI,1.0/3.0) * std::pow(6.0*vol_i,2.0/3.0)/surf_area_i;
                solid->setSphericity(sphericity);
                //solid->SetVolume(vol_diff_i);
                //if(init)solid->SetVolume_0(vol_diff_i);
                //amrex::Print()<<"Volume of IF "<<i<<" at lev "<<ilev<<" = "<<solid->Volume()<<'\n';
                //amrex::Print()<<"sphericity of IF "<<i<<" at lev "<<ilev<<" = "<<solid->getSphericity()<<'\n';
            }
        }

    }

    amrex::Real GetPhaseFieldRHS(amrex::Array4<amrex::Real const> const& Psi, int i, int j, int k)
    {
        return 0.0;
    }

    void incFSI::DetectBubbles(bool init)
    {
        for (int i = 0; i < N_IF; i++)
        {
            for(int ilev = 0; ilev <= finest_level; ilev++)
            {
                auto &&solid = interfaces[ilev][i];
                solid->DetectBubbles(true);
                amrex::Print()<<"Inside incFSI_interface : DetectBubbles()"<<"\n";
            }
        }
    }

} // namespace mycode
