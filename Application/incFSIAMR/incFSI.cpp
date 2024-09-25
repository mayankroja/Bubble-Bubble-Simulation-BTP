#include "incFSI.H"
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <CFMask.H>

namespace mycode
{

incFSI::incFSI(const std::vector<std::unique_ptr<interface_props>> &IF_ref):IF_props(IF_ref)
{
    Nghost = 1;
    TUBE_LAYERS = 10;
    Time = 0.0;
    dt = 5e-7;
    Iter = 0;
    Projection = false;
    int nlevels = max_level + 1;

    amrex::ParmParse pp;
    pp.query("Axisymmetric", isAxisymmetric);
    pp.query("RKOrder", RKOrder);
    pp.query("CFL",CFL);
    amrex::Print()<<"Axisymmetric calcualtion = "<<isAxisymmetric<<'\n';
    pp.get("FinalTime", FinalTime);
    pp.get("PlotInt", PlotInt);
    pp.get("ReinitInt", ReinitInt);
    pp.query("Restart", Restart);// default false
    pp.query("chk_int", chk_int);// default -1
    pp.query("chk_int_read", chk_int_read);
    pp.query("Mu", Mu);
    pp.query("gravity",gravity);
    if(gravity != 0.0)
        amrex::Print()<<"Gravity is on"<<'\n';
    pp.query("SIGMA",SIGMA);
    pp.get("Nghost", Nghost);
    pp.get("LAYERS",TUBE_LAYERS);
    pp.get("MaxIter", MaxIter);
    pp.query("dt", dt);
    pp.query("StopAtSteady", StopAtSteady);//default false
    pp.query("tag_vort", tag_vort);//default false
    pp.query("tag_region", tag_region);//default false
    pp.query("regrid_int", regrid_int);//default 2
    pp.query("SteadyTol", SteadyTol);//default 1e-5
    pp.query("divergence_free_interpolation",divergence_free_interpolation);//default false
    pp.query("Temperature_field",TempField);//default false
    pp.query("Phase_field",PhaseField);//default false
    pp.query("sharp_Phase_field",sharpPhaseField);
    pp.query("b_coeff",b_coeff);// A parameter for controlling the normal velocity of interface for phase-field interface sharpening
    pp.query("W_phasefield",W_phasefield);// A parameter for controlling the normal velocity of interface for phase-field interface sharpening
    pp.query("DamageModel",DamageModel);
    pp.query("damage_coeff",damage_coeff);
    pp.query("material_type",material_type);
    pp.query("Use_Algoim",use_algoim);//default false
    pp.query("Tecplot",tecplot);//default false
    pp.query("advect_ref_cond",advect_ref_cond);
    //Read plot options
    pp.query("plot_viscosity",plot_viscosity);
    pp.query("plot_strain_error",plot_strain_error);
    pp.query("LS_reg_iter",LS_reg_iter);
    pp.query("Reinit_Method",Reinit_Method);
    pp.query("Fij_order",Fij_order);
    pp.query("ComputeMomentumBalance",ComputeMomentumBalance_);

    grids_xvel.resize(nlevels);
    grids_yvel.resize(nlevels);

    xvel.resize(nlevels);
    xvel_old.resize(nlevels);

    yvel.resize(nlevels);
    yvel_old.resize(nlevels);

    U.resize(nlevels);
    Pressure.resize(nlevels);
    Pressure_old.resize(nlevels);
    RHS.resize(nlevels);

    if(TempField)
    {
        Theta.resize(nlevels);
        Theta_old.resize(nlevels);
    }

    if(PhaseField)
    {
        Phi.resize(nlevels); 
        Phi_old.resize(nlevels);
    }
    if(DamageModel)
    {
	amrex::Print()<<"Damage model is on."<<'\n';
    Scalars.resize(nlevels);
	for(int lev = 0;lev < nlevels;lev++)
    Scalars[lev].resize(nscalar);
    Scalars_old.resize(nlevels);
    for(int lev = 0;lev < nlevels;lev++)
    Scalars_old[lev].resize(nscalar);
	Scalars_bcf.resize(nscalar);
	lobc_Scalars.resize(nscalar);
        hibc_Scalars.resize(nscalar);	
    }

    Pstar.resize(nlevels);
    Vorticity.resize(nlevels);
    Umag.resize(nlevels);
    if (StopAtSteady)
        UmagOld.resize(nlevels);
    if (RKOrder == 2)
    {
        FRK_xvel.resize(nlevels);
        FRK_yvel.resize(nlevels);
        FRK_p.resize(nlevels);
        if(TempField) FRK_Theta.resize(nlevels);
        if(PhaseField) FRK_Phi.resize(nlevels);
        if(DamageModel) 
        {
	    FRK_Scalars.resize(nlevels);
            for(int lev = 0;lev < nlevels;lev++)
                FRK_Scalars[lev].resize(nscalar);
        }
    }
    if (RKOrder >= 3)
    {
        RK1_xvel.resize(nlevels);
        RK1_yvel.resize(nlevels);
        RK1_p.resize(nlevels);
        if(TempField) RK1_Theta.resize(nlevels);
        if(PhaseField) RK1_Phi.resize(nlevels);
        if(DamageModel)
        {
            RK1_Scalars.resize(nlevels);
            for(int lev = 0;lev < nlevels;lev++)
                RK1_Scalars[lev].resize(nscalar);
        }

        RK2_xvel.resize(nlevels);
        RK2_yvel.resize(nlevels);
        RK2_p.resize(nlevels);
        if(TempField) RK2_Theta.resize(nlevels);
        if(PhaseField) RK2_Phi.resize(nlevels);
        if(DamageModel)
        {
            RK2_Scalars.resize(nlevels);
            for(int lev = 0;lev < nlevels;lev++)
                RK2_Scalars[lev].resize(nscalar);
        }

        RK3_xvel.resize(nlevels);
        RK3_yvel.resize(nlevels);
        RK3_p.resize(nlevels);
        if(TempField) RK3_Theta.resize(nlevels);
        if(PhaseField) RK3_Phi.resize(nlevels);
        if(DamageModel)
        {
            RK3_Scalars.resize(nlevels);
            for(int lev = 0;lev < nlevels;lev++)
                RK3_Scalars[lev].resize(nscalar);
        }

        Src_xvel.resize(nlevels);
        Src_yvel.resize(nlevels);
        Src_p.resize(nlevels);
        if(TempField) Src_Theta.resize(nlevels);
        if(PhaseField) Src_Phi.resize(nlevels);
        if(DamageModel)
        {
            Src_Scalars.resize(nlevels);
            for(int lev = 0;lev < nlevels;lev++)
                Src_Scalars[lev].resize(nscalar);
        }
    }
    
    if (tag_region)
    {
        N_tag_regions = pp.countval("tag_region_name");
        tag_region_name.resize(N_tag_regions);
        pp.getarr("tag_region_name", tag_region_name);
        tag_region_lo.resize(AMREX_SPACEDIM);
        tag_region_hi.resize(AMREX_SPACEDIM);

        //! find the boxes corresponding to tag region on each level
        boxlist.resize(nlevels);
        for (int i = 0; i < N_tag_regions; i++)
        {
            amrex::ParmParse pp(tag_region_name[i]);
            pp.getarr("lo", tag_region_lo);
            pp.getarr("hi", tag_region_hi);

            amrex::RealBox rb(tag_region_lo.data(), tag_region_hi.data());

            for (int lev = 0; lev <= max_level; lev++)
            {
                BL_ASSERT(geom[lev].ProbDomain().contains(rb));
                const amrex::Real *dx = geom[lev].CellSize();
                const amrex::Real *prob_lo = geom[lev].ProbLo();
                int xlo = static_cast<int>((tag_region_lo[0] - prob_lo[0])/dx[0]);
                int ylo = static_cast<int>((tag_region_lo[1] - prob_lo[1])/dx[1]);
                int xhi = static_cast<int>((tag_region_hi[0] - prob_lo[0])/dx[0]);
                int yhi = static_cast<int>((tag_region_hi[1] - prob_lo[1])/dx[1]);
                amrex::IntVect lo(xlo, ylo);
                amrex::IntVect hi(xhi, yhi);
                amrex::Box bx_mark(lo, hi);

                boxlist[lev].push_back(bx_mark);
            }            
        }
    }

    setBCs(lobc_u, "lobc_u");
    setBCs(hibc_u, "hibc_u");
    setBCs(lobc_v, "lobc_v");
    setBCs(hibc_v, "hibc_v");
    setBCs(lobc_p, "lobc_p");
    setBCs(hibc_p, "hibc_p");
    if(TempField)
    {
        setBCs(lobc_T, "lobc_T");
        setBCs(hibc_T, "hibc_T");
    }
    if(PhaseField)
    {
        setBCs(lobc_phi, "lobc_Phi");
        setBCs(hibc_phi, "hibc_Phi");
    }
    //for(int iscalar; iscalar < nscalar; iscalar++)
    if(DamageModel)
    {
        setBCs(lobc_Scalars[0], "lobc_X");
        setBCs(hibc_Scalars[0], "hibc_X");

        setBCs(lobc_Scalars[1], "lobc_Y");
        setBCs(hibc_Scalars[1], "hibc_Y");

        setBCs(lobc_Scalars[2], "lobc_F11");
        setBCs(hibc_Scalars[2], "hibc_F11");

        setBCs(lobc_Scalars[3], "lobc_F12");
        setBCs(hibc_Scalars[3], "hibc_F12");

	setBCs(lobc_Scalars[4], "lobc_F21");
        setBCs(hibc_Scalars[4], "hibc_F21");

        setBCs(lobc_Scalars[5], "lobc_F22");
        setBCs(hibc_Scalars[5], "hibc_F22");

        setBCs(lobc_Scalars[6], "lobc_F33");
        setBCs(hibc_Scalars[6], "hibc_F33");

        setBCs(lobc_Scalars[7], "lobc_ep");
        setBCs(hibc_Scalars[7], "hibc_ep");	

	setBCs(lobc_Scalars[8], "lobc_DMG");
        setBCs(hibc_Scalars[8], "hibc_DMG");
    }
    //Read info for creating interface
    pp.get("N_interfaces", N_IF);
    Pintbcs.resize(N_IF);
    Uintbcs.resize(N_IF);
    Vintbcs.resize(N_IF);
    Tintbcs.resize(N_IF);
    Phiintbcs.resize(N_IF);
    intbcf.resize(N_IF);
    IF_names.resize(N_IF);
    IF_types.resize(N_IF);
    if (N_IF > 0)
    {
        pp.getarr("IF_names", IF_names);
        pp.getarr("IF_types", IF_types);

        setIntBCs(Pintbcs, "p_bc");
        setIntBCs(Uintbcs, "u_bc");
        setIntBCs(Vintbcs, "v_bc");
        if(TempField)setIntBCs(Tintbcs, "T_bc");
        if(PhaseField)setIntBCs(Phiintbcs, "Phi_bc");

        tag_interface = true;
    }
    interfaces.resize(nlevels);
    mask_.resize(nlevels);

    amrex::Print()<<"No. of levels = "<<max_level + 1<<'\n';
    amrex::Print()<<"No. of interfaces = "<<N_IF<<'\n';
}

incFSI::~incFSI() {}

void incFSI::Evolve()
{
    do
    {
        t_old = Time;
        Time += dt;
        t_new = Time;
        Iter++;

        amrex::Print() << "\n\nadvance time = " << Time << "\tIter = " << Iter << "\tdt = " << dt << "\n";
        amrex::PrintToFile("log") << "\n\nadvance time = " << Time << "\tIter = " << Iter << "\tdt = " << dt << "\n";

        FillPatchAroundBox();
 
        if(RKOrder == 1)
            Scheme();
        else if(RKOrder == 2)
            Scheme_RK2();
        else if(RKOrder >= 3)
            Scheme_RK();

        ComputeDt();

        //Write Output        
        if (Iter % PlotInt == 0 || Iter == 1)
        {
            WriteFile();
            if(tecplot) 
            {
                WriteFileTecplot();
                WriteFileTecplotW_Ghost();
            }
            if(N_IF > 0)WriteInterface();
	    if(PhaseField)CheckConservationPhaseField();
        }
        if(ComputeMomentumBalance_) ComputeMomentumBalance();

       //Write restart files
       if(Iter % chk_int == 0)
       {
           amrex::Print()<<"Writing chkpt"<<'\n';
           ChkFile();
       }

        amrex::PrintToFile("log") << "max vel in domain = " << MaxU << "\n";

        if (StopAtSteady)
        {
            amrex::PrintToFile("log") << "steady deviation = " << SteadyDeviation << "\n";

            if (IsSteady())
            {
                amrex::PrintToFile("log") << "steady state reached..\n";
                WriteFile();
                break;
            }
        }
        
        writeDataFile();
    } 
    while (Time < FinalTime && Iter < MaxIter);
}

void incFSI::AllocateMemory(int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    //amrex::Print()<<"Allocating for lev = "<<lev<<" , ba = "<<ba<<'\n';   

    amrex::BoxArray xba = amrex::convert(ba, amrex::IntVect::TheDimensionVector(0));
    amrex::BoxArray yba = amrex::convert(ba, amrex::IntVect::TheDimensionVector(1));

    grids_xvel[lev] = xba;
    grids_yvel[lev] = yba;

    /// face centered data
    xvel[lev].define(xba, dm, 1, Nghost);
    yvel[lev].define(yba, dm, 1, Nghost);

    xvel_old[lev].define(xba, dm, 1, Nghost);
    yvel_old[lev].define(yba, dm, 1, Nghost);

    /// cell centered data
    U[lev].define(ba, dm, AMREX_SPACEDIM, Nghost);
    Pressure[lev].define(ba, dm, 1, Nghost);
    Pressure_old[lev].define(ba, dm, 1, Nghost);
    Pstar[lev].define(ba, dm, 1, Nghost);
    RHS[lev].define(ba, dm, 1, 0);
    Vorticity[lev].define(ba, dm, 1, 0);
    if(TempField)
    {
        Theta[lev].define(ba, dm, 1, Nghost);
        Theta_old[lev].define(ba, dm, 1, Nghost);
    }
    if(PhaseField)
    {
        Phi[lev].define(ba, dm, 1, Nghost);
        Phi_old[lev].define(ba, dm, 1, Nghost);
    }
    if(DamageModel)
    {
        for(int iscalar = 0; iscalar < nscalar; iscalar++)
	{
	    Scalars[lev][iscalar].define(ba, dm, 1, Nghost);
	    Scalars_old[lev][iscalar].define(ba, dm, 1, Nghost);
	}
    }
    /// derived fields
    Umag[lev].define(ba, dm, 1, 0);
    if (StopAtSteady)
    {
        UmagOld[lev].define(ba, dm, 1, 0);
    }
    
    if (RKOrder == 2)
    {
        FRK_xvel[lev].define(xba, dm, 1, Nghost);
        FRK_yvel[lev].define(yba, dm, 1, Nghost);
        FRK_p[lev].define(ba, dm, 1, Nghost);
        if(TempField) FRK_Theta[lev].define(ba, dm, 1, Nghost);
        if(PhaseField) FRK_Phi[lev].define(ba, dm, 1, Nghost); 
	if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
            {
                FRK_Scalars[lev][iscalar].define(ba, dm, 1, Nghost);
            }
        }

    }
    if (RKOrder >= 3)
    {
        RK1_xvel[lev].define(xba, dm, 1, Nghost);
        RK1_yvel[lev].define(yba, dm, 1, Nghost);
        RK1_p[lev].define(ba, dm, 1, Nghost);
        if(TempField) RK1_Theta[lev].define(ba, dm, 1, Nghost);
        if(PhaseField) RK1_Phi[lev].define(ba, dm, 1, Nghost);
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
            {   
                RK1_Scalars[lev][iscalar].define(ba, dm, 1, Nghost);
            }
        }


        RK2_xvel[lev].define(xba, dm, 1, Nghost);
        RK2_yvel[lev].define(yba, dm, 1, Nghost);
        RK2_p[lev].define(ba, dm, 1, Nghost);
        if(TempField) RK2_Theta[lev].define(ba, dm, 1, Nghost);
        if(PhaseField) RK2_Phi[lev].define(ba, dm, 1, Nghost);
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
            {   
                RK2_Scalars[lev][iscalar].define(ba, dm, 1, Nghost);
            }
        }


        RK3_xvel[lev].define(xba, dm, 1, Nghost);
        RK3_yvel[lev].define(yba, dm, 1, Nghost);
        RK3_p[lev].define(ba, dm, 1, Nghost);
        if(TempField) RK3_Theta[lev].define(ba, dm, 1, Nghost);
        if(PhaseField) RK3_Phi[lev].define(ba, dm, 1, Nghost);
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
            {   
                RK3_Scalars[lev][iscalar].define(ba, dm, 1, Nghost);
            }
        }

        Src_xvel[lev].define(xba, dm, 1, Nghost);
        Src_yvel[lev].define(yba, dm, 1, Nghost);
        Src_p[lev].define(ba, dm, 1, Nghost);
        if(TempField) Src_Theta[lev].define(ba, dm, 1, Nghost);
        if(PhaseField) Src_Phi[lev].define(ba, dm, 1, Nghost);
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
            {   
                Src_Scalars[lev][iscalar].define(ba, dm, 1, Nghost);
            }
        }
    }
}

void incFSI::Initialize()
{
    amrex::Print()<<"Before Initialization"<<'\n';

    if(Restart)
    {
        ReadCheckpointFile();
	    ApplyBC();
        AMGPressurePoisson();
        FillPatchAroundBox();
        //WriteFile();
        //WriteFileTecplot();
    	//if(N_IF > 0)WriteInterface();
    }
    else
    {
        // start simulation from the beginning
        const amrex::Real time = 0.0;
        //! this will call MakeNewLevelFromScratch
        //! which does initialisation of data on each level
        InitFromScratch(time);
	    ApplyBC();
        if(N_IF > 0)
        {
            ComputeBubbleVolume(true);
            //for(int lev = 0; lev <= finest_level; lev++)
	    int lev = finest_level;
            {
                ComputeCutFaceVel(lev);
                ComputeCutFacePressure(lev);
                ComputeIntProp(lev);
            }
        }
        //! Multi-level Poisson solve
        AMGPressurePoisson();
        FillPatchAroundBox();
    }
    InitializeRKcoeffs();
    AverageDown();
    amrex::Print()<<"Initialization done..."<<'\n';
    //WriteFile();
    //amrex::Print()<<"after writing"<<'\n';
    //if(tecplot)WriteFileTecplot();
    //if(N_IF > 0)WriteInterface();
    //WriteFileTecplotW_Ghost();
    //amrex::ParallelDescriptor::Barrier();
    //for(int lev = 0; lev <finest_level; lev++ )
    //   FillPatchNearBoundary(lev); 
    //amrex::ParallelDescriptor::Barrier();
    //std::exit(5);
}

void incFSI::Scheme()
{
    //! Regrid
    Regrid();
    
    if(N_IF > 0)
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            for (auto &&solid : interfaces[lev])
            {
		solid->setVolume_prev();
                solid->TubeIdentification();
                mask_[lev]->GhostCellIdentfication();
                solid->TubeAdvectLevelSet(xvel[lev], yvel[lev], dt); 
                if (Iter % ReinitInt == 0)
                {
                    if(use_algoim) 
                         solid->Reinit_algoim();
                    else
                         solid->Reinit2();
                    
                    //solid->Reinit();
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

    /// Prediction
    for (int lev = 0; lev <= finest_level; lev++)
        ComputeIntermediateVelocity(lev);

    if(DamageModel)
    {
        for (int lev = 0; lev <= finest_level; lev++)
	{
	    AdvectScalars(lev);
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

    //! Avg down face centered velocity
    ApplyBC();

    AverageDown();

    if(DamageModel)
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            ComputeDeformationGradientTensor2D(lev);
        }
    }
    

    for (int lev = 0; lev <= finest_level; lev++)
    {
        if(N_IF == 0)ComputeCollocatedVelocityField(lev);
    }
    Derive();
}

void incFSI::SetInitialFlowField(int lev)
{
    xvel[lev].setVal(0.0);
    yvel[lev].setVal(0.0);
    xvel_old[lev].setVal(0.0);
    yvel_old[lev].setVal(0.0);
    U[lev].setVal(0.0);
    Pressure[lev].setVal(0.0);
    Pressure_old[lev].setVal(0.0);
    Pstar[lev].setVal(0.0);

    if(TempField)
    {
        Theta[lev].setVal(0.0);
        Theta_old[lev].setVal(0.0);
    }
    if(PhaseField)
    {
        Phi[lev].setVal(0.0);
        Phi_old[lev].setVal(0.0);
    }
    if(DamageModel)
    {
        for(int iscalar = 0;iscalar < nscalar; iscalar++)
        {
	    Scalars[lev][iscalar].setVal(0.0);
	    Scalars_old[lev][iscalar].setVal(0.0);
	}
    }

    const amrex::Real *prob_lo = geom[lev].ProbLo();
    const amrex::Real *dx = geom[lev].CellSize();

    amrex::iMultiFab *PMask;
    amrex::iMultiFab *UMask;
    amrex::iMultiFab *VMask;
    if(N_IF > 0)
    {
        auto &&mask__= mask_[lev];
        PMask = &mask__->getPMask();
        UMask = &mask__->getUMask();
        VMask = &mask__->getVMask();
    }

    //! initialise fields with given function
    if (u_init)
    {
        for (amrex::MFIter mfi(xvel[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox(); 
            amrex::Array4<amrex::Real> const &u = xvel[lev].array(mfi);
            amrex::Array4<amrex::Real> const &u_old = xvel_old[lev].array(mfi);
            amrex::Array4<int> umask;
            if(N_IF > 0) umask = UMask->array(mfi);
            //amrex::Array4<int const> const &umask = mask_[lev]->getUMask().const_array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real x = prob_lo[0] + dx[0] * i; 
                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                if(N_IF == 0)
                {
                    u(i, j, k) = u_init(x, y, Time);
                    u_old(i, j, k) = u_init(x, y, Time);
                }
                else
                {
                    if(umask(i,j,k) == 1)
                    {
                        u(i, j, k) = u_init(x, y, Time);
                        u_old(i, j, k) = u_init(x, y, Time);
                    }
                }
            });
        }
    }

    if (v_init)
    {
        for (amrex::MFIter mfi(yvel[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox(); 
            amrex::Array4<amrex::Real> const &v = yvel[lev].array(mfi);
            amrex::Array4<amrex::Real> const &v_old = yvel_old[lev].array(mfi);
            amrex::Array4<int> vmask;
            if(N_IF > 0) vmask = VMask->array(mfi);
            //amrex::Array4<int const> const &vmask = mask_[lev]->getVMask().const_array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                amrex::Real y = prob_lo[1] + dx[1] * j;
                if(N_IF == 0)
                {
                    v(i, j, k) = v_init(x, y, Time);
                    v_old(i, j, k) = v_init(x, y, Time);
                }
                else
                {   
                    if(vmask(i,j,k) == 1)
                    {   
                        v(i, j, k) = v_init(x, y, Time);
                        v_old(i, j, k) = v_init(x, y, Time);
                    }
                }
            });
        }
    }

    if (P_init)
    {
        for (amrex::MFIter mfi(Pressure[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real> const &P = Pressure[lev].array(mfi);
            amrex::Array4<amrex::Real> const &Pold = Pstar[lev].array(mfi);
            amrex::Array4<int> pmask;
            if(N_IF > 0) pmask = PMask->array(mfi);
            //amrex::Array4<int const> const &pmask = mask_[lev]->getPMask().const_array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                if(N_IF == 0)
                {   
                    P(i, j, k) = P_init(x, y, Time);
                    Pold(i, j, k) = P_init(x, y, Time);
                }
                else
                {
                    if(pmask(i,j,k) == 1)
                    {   
                        P(i, j, k) = P_init(x, y, Time);
                        Pold(i, j, k) = P_init(x, y, Time);
                    }
                }
            });
        }
    }

    if (Theta_init && TempField)
    {
        for (amrex::MFIter mfi(Theta[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real> const &T = Theta[lev].array(mfi);
            amrex::Array4<int> pmask;
            if(N_IF > 0) pmask = PMask->array(mfi);
            //amrex::Array4<int const> const &pmask = mask_[lev]->getPMask().const_array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                amrex::Real ds = 2.0 * std::min(dx[0], dx[1]);
                if(N_IF == 0)
                {
                    T(i, j, k) = Theta_init(x, y, ds);
                }
                else
                {
                    if(pmask(i,j,k) == 1)
                    {
                        T(i, j, k) = Theta_init(x, y, ds);
                    }
                }
            });
        }
    }

    if (Phi_init && PhaseField)
    {
        for (amrex::MFIter mfi(Phi[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real> const &phi = Phi[lev].array(mfi);
            amrex::Array4<int> pmask;
            if(N_IF > 0) pmask = PMask->array(mfi);
            //amrex::Array4<int const> const &pmask = mask_[lev]->getPMask().const_array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                amrex::Real ds = 2.0 * std::min(dx[0], dx[1]); 
                if(N_IF == 0)
                {
                    phi(i, j, k) = Phi_init(x, y, ds);
                }
                else
                {
                    if(pmask(i,j,k) == 1)
                    {
                        phi(i, j, k) = Phi_init(x, y, ds);
                    }
                }
            });
        }
    }
    if(DamageModel)
    {
        //Scalar 1 and 2 are X and Y of the reference 
        //Scalar 3 - 6 are the components of Deformation gradient tensor
	//Reference state is undeformed state, i.e., F = I
        for (amrex::MFIter mfi(Pressure[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real> const &sc1 = Scalars[lev][0].array(mfi);
	    amrex::Array4<amrex::Real> const &sc2 = Scalars[lev][1].array(mfi);
	    amrex::Array4<amrex::Real> const &sc3 = Scalars[lev][2].array(mfi);
	    amrex::Array4<amrex::Real> const &sc4 = Scalars[lev][3].array(mfi);
	    amrex::Array4<amrex::Real> const &sc5 = Scalars[lev][4].array(mfi);
	    amrex::Array4<amrex::Real> const &sc6 = Scalars[lev][5].array(mfi);
	    amrex::Array4<amrex::Real> const &sc7 = Scalars[lev][6].array(mfi);
            amrex::Array4<int> pmask;
            if(N_IF > 0) pmask = PMask->array(mfi);
            //amrex::Array4<int const> const &pmask = mask_[lev]->getPMask().const_array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                amrex::Real ds = 2.0 * std::min(dx[0], dx[1]);
                //if(N_IF == 0)
                {
                    sc1(i, j, k) = x;
                    sc2(i, j, k) = y;
		    sc3(i, j, k) = 1.0;
		    sc4(i, j, k) = 0.0;
		    sc5(i, j, k) = 0.0;
		    sc6(i, j, k) = 1.0;
		    sc7(i, j, k) = 1.0;
                }
                /*else
                {
                    if(pmask(i,j,k) == 1)
                    {
                        sc1(i, j, k) = x;
                        sc2(i, j, k) = y;
                        sc3(i, j, k) = 1.0;
                        sc4(i, j, k) = 0.0;
                        sc5(i, j, k) = 0.0;
                        sc6(i, j, k) = 1.0;
			sc7(i, j, k) = 1.0;
                    }
                }
		*/
            });
        }    

    }


    xvel[lev].FillBoundary();
    yvel[lev].FillBoundary();
    Pressure[lev].FillBoundary();
}

void incFSI::CopyNewToOld()
{
    for (int lev = 0; lev <= finest_level; lev++)
    {
        amrex::MultiFab::Copy(xvel_old[lev], xvel[lev], 0, 0, 1, Nghost);
        amrex::MultiFab::Copy(yvel_old[lev], yvel[lev], 0, 0, 1, Nghost);
        amrex::MultiFab::Copy(Pstar[lev], Pressure[lev], 0, 0, 1, Nghost);
        amrex::MultiFab::Copy(Pressure_old[lev], Pressure[lev], 0, 0, 1, Nghost);
        if(TempField)amrex::MultiFab::Copy(Theta_old[lev], Theta[lev], 0, 0, 1, Nghost);
        if(PhaseField)amrex::MultiFab::Copy(Phi_old[lev], Phi[lev], 0, 0, 1, Nghost);
	if(DamageModel)
	{
	    for(int iscalar = 0; iscalar < nscalar; iscalar++)
            {
                if(iscalar == eps_max_num)
                    continue;
	        amrex::MultiFab::Copy(Scalars_old[lev][iscalar], Scalars[lev][iscalar], 0, 0, 1, Nghost); 
	    }
	}

        if (StopAtSteady)
        {
            amrex::MultiFab::Copy(UmagOld[lev], Umag[lev], 0, 0, 1, 0);
        }
    }  
}

bool incFSI::IsSteady()
{
    SteadyDeviation = 0.0;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        for (amrex::MFIter mfi(Umag[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real> const &vel_mag = Umag[lev].array(mfi);
            amrex::Array4<amrex::Real> const &vel_mag_old = UmagOld[lev].array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                SteadyDeviation = std::max(SteadyDeviation, std::fabs(vel_mag(i, j, k) - vel_mag_old(i, j, k)));
            });
        }
    }

    amrex::ParallelDescriptor::ReduceRealMax(SteadyDeviation);
    return SteadyDeviation <= SteadyTol;
}

void incFSI::ComputeDt()
{
    amrex::Real dt_min = 1e100;
    amrex::Real umax = 0.0, vmax = 0.0;
    for (int lev = 0; lev <= finest_level; lev++)
    {
        //umax = U[lev].max(0);
        //vmax = U[lev].max(1);
        amrex::iMultiFab *PMask;
        amrex::iMultiFab PMask_NI;//NI stands for No Interface, need a better way of doing this
	CFMask cfmask_;
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

	cfmask_.define(lev, geom, grids, dmap);

        for (amrex::MFIter mfi(U[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<int const> mask;
	    amrex::Array4<int const> const &cfmask = cfmask_.Mask().const_array(mfi);

            if(N_IF > 0)
                mask = PMask->const_array(mfi);
            else
                mask = PMask_NI.const_array(mfi);

            amrex::Array4<amrex::Real const> const &vel = U[lev].const_array(mfi);
            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (mask(i, j, k) == 1 && cfmask(i,j,k) != CFMask::covered)
                {
                    umax = std::max(fabs(vel(i, j, k, 0)), umax);
                    vmax = std::max(fabs(vel(i, j, k, 1)), vmax);
                }
            });
        }
        /// get max over all processes
        amrex::ParallelDescriptor::ReduceRealMax(umax);
        amrex::ParallelDescriptor::ReduceRealMax(vmax);
        amrex::ParallelDescriptor::ReduceRealMax(Mu_max);
        if(TempField)amrex::ParallelDescriptor::ReduceRealMax(k_max);

	//amrex::Print(-1)<<"Mu_max = "<<Mu_max<<'\n';
	//amrex::Print(-1)<<"umax = "<<umax<<" , vmax = "<<vmax<<'\n';

        const amrex::Real *dx = geom[lev].CellSize();

	//amrex::Print(-1)<<"dx = "<<dx[0]<<" , "<<dx[1]<<'\n';

        //amrex::Real dt_diff = 0.2 * (dx[0] * dx[0] + dx[1] * dx[1]) / (2.0 * (Mu_max + 1e-15));
        //amrex::Real dt_diff = 0.2 * (1.0/Mu_max)/(1.0/(dx[0] * dx[0] + 1e-15) + 1.0/(dx[1] * dx[1] + 1e-15));
        //if(TempField) dt_diff = std::min(0.2 * (dx[0] * dx[0] + dx[1] * dx[1]) / (2.0 * (k_max + 1e-15)), dt_diff);
        //amrex::Real dt_adv = 0.4 / (umax / dx[0] + vmax / dx[1] + 1e-15);

        amrex::Real C_c = (umax / dx[0] + vmax / dx[1] + 1e-15);
        amrex::Real C_v = 2.0*((Mu_max + 1e-15)*(1.0/(dx[0]*dx[0]) + 1.0/(dx[1]*dx[1])));
	amrex::Real C_sft = 0.0;
        amrex::Real C_gravity = 0.0;
	if(lev == finest_level)
        {
            for (int i = 0; i < N_IF; i++)
            {
                auto &&solid = interfaces[lev][i];
		C_sft = std::max(C_sft,
                                 std::sqrt(SIGMA/std::pow(std::min(dx[0],dx[1]),3)));
	    }
            C_gravity = std::sqrt(std::abs(gravity)/dx[0]);
	}
        //if(TempField) C_v = std::max(C_v, 2.0*((k_max + 1e-15)*(1.0/(dx[0]*dx[0]) + 1.0/(dx[1]*dx[1]))));
	//amrex::Print(-1)<<"C_c = "<<C_c<<"; C_v = "<<C_v<<"; C_sft = "<<C_sft<<'\n';
  
        amrex::Real C_tot = 0.5*(C_c + C_v + 
			         std::sqrt( std::pow((C_c + C_v),2) + 4*std::pow(C_sft,2) + 4*std::pow(C_gravity,2)) );
        amrex::Real dt_lev = CFL/C_tot;
	//if(Iter == 1) dt_lev = 0.05/(C_c + C_v);

        //amrex::Real dt_lev = std::min(dt_adv, dt_diff);
        dt_min = std::min(dt_min, dt_lev);
    }
    dt = dt_min;
    if(PhaseField)
    {
        const amrex::Real *dx = geom[finest_level].CellSize();
	b_sharp = b_coeff * (3.0/10.0) * dx[0]*dx[0]/dt * (1 - 2.0 * 0.4); 
	//amrex::Print()<<"b_sharp = "<<b_sharp<<'\n';
    }
    //if (finest_level > 0 && Iter % regrid_int == (regrid_int - 1)) dt = dt/1000.0;
    //amrex::Print()<<"Iter % regrid_int = "<<Iter % regrid_int<<'\n';
}

void incFSI::ApplyBC()//At immersed interface and domain boundary
{
    FillPatchAroundBox();

    if(N_IF > 0) ComputeBubbleVolume(false);

    for (int lev = 0; lev <= finest_level; lev++)
    {
        if(N_IF > 0 && lev == finest_level)
        {
            ComputeCutFaceVel(lev);
            ComputeCutFacePressure(lev);
        }
    }

    if(N_IF > 0)ComputeIntProp(finest_level);

    for (int lev = 0; lev <= finest_level; lev++)
    {
        XVelBoundaryConditions(lev, Time, xvel[lev]);
        YVelBoundaryConditions(lev, Time, yvel[lev]);
        if(TempField)TemperatureBoundaryConditions(lev, Time, Theta[lev]);
        if(PhaseField)PhiBoundaryConditions(lev, Time, Phi[lev]);
	if(DamageModel)
        {
	    for(int iscalar = 0;iscalar < nscalar;iscalar++)
            {
		if(iscalar == eps_max_num)
		    continue;

	        ScalarBoundaryConditions(lev, iscalar, Time, Scalars[lev][iscalar]);
	    }
	}
        /*
        if(xvel[lev].contains_nan())
        {
            amrex::PrintToFile("log") << "xvel* on lev " << lev << " contains nan\n";
        }
        if(xvel[lev].contains_inf())
        {
            amrex::PrintToFile("log") << "xvel* on lev " << lev << " contains inf\n";
        }
        if(yvel[lev].contains_nan())
        {
            amrex::PrintToFile("log") << "yvel* on lev " << lev << " contains nan\n";
        }
        if(yvel[lev].contains_inf())
        {
            amrex::PrintToFile("log") << "yvel* on lev " << lev << " contains inf\n";
        }
        */
    }

}


/// set all sides to same function
void incFSI::SetXVelBoundaryFunction(BoundaryFunc f)
{
    for (size_t i = 0; i < 2*AMREX_SPACEDIM; i++)
    {
        u_bcf[i] = f;
    }
}

/// set side i to function f
void incFSI::SetXVelBoundaryFunction(int i, BoundaryFunc f)
{
    u_bcf[i] = f;
}

/// set all sides to same function
void incFSI::SetYVelBoundaryFunction(BoundaryFunc f)
{
    for (size_t i = 0; i < 2*AMREX_SPACEDIM; i++)
    {
        v_bcf[i] = f;
    }
}

/// set side i to function f
void incFSI::SetYVelBoundaryFunction(int i, BoundaryFunc f)
{
    v_bcf[i] = f;
}

/// set all sides to same function
void incFSI::SetPressureBoundaryFunction(BoundaryFunc f)
{
    for (size_t i = 0; i < 2*AMREX_SPACEDIM; i++)
    {
        p_bcf[i] = f;
    }
}

/// set side i to function f
void incFSI::SetPressureBoundaryFunction(int i, BoundaryFunc f)
{
    p_bcf[i] = f;
}

/// set all sides to same function
void incFSI::SetTemperatureBoundaryFunction(BoundaryFunc f)
{
    for (size_t i = 0; i < 2*AMREX_SPACEDIM; i++)
    {
        T_bcf[i] = f;
    }
}

/// set side i to function f
void incFSI::SetTemperatureBoundaryFunction(int i, BoundaryFunc f)
{
    T_bcf[i] = f;
}

/// set all sides to same function
void incFSI::SetPhiBoundaryFunction(BoundaryFunc f)
{
    for (size_t i = 0; i < 2*AMREX_SPACEDIM; i++)
    {
        Phi_bcf[i] = f;
    }
}

/// set side i to function f
void incFSI::SetPhiBoundaryFunction(int i, BoundaryFunc f)
{
    Phi_bcf[i] = f;
}

/// set all sides to same function
void incFSI::SetScalarsBoundaryFunction(int iscalar, BoundaryFunc f)
{
    for (size_t i = 0; i < 2*AMREX_SPACEDIM; i++)
    {
        Scalars_bcf[iscalar][i] = f;
    }
}

/// set side i to function f
void incFSI::SetScalarsBoundaryFunction(int iscalar, int i, BoundaryFunc f)
{
    Scalars_bcf[iscalar][i] = f;
}

/// return valid face box in direction 'dir' 
amrex::Box incFSI::get_valid_face_box(int lev, const amrex::Box& bx_i, int dir)
{
	amrex::Box bx_o(bx_i);
	std::array<amrex::LinOpBCType, AMREX_SPACEDIM> *lobc, *hibc;
	if(dir == 0)
	{
		lobc = &lobc_u;
		hibc = &hibc_u;
	}
	else if(dir==1)
	{
		lobc = &lobc_v;
		hibc = &hibc_v;	
	}
	
	const amrex::Box& domain = geom[lev].Domain();
	int d_lo = domain.smallEnd(dir);
	int d_hi = domain.bigEnd(dir);

	int b_lo = bx_o.smallEnd(dir);
	int b_hi = bx_o.bigEnd(dir);
	
	if( b_lo == d_lo && (*lobc)[dir] == amrex::LinOpBCType::Dirichlet )
	{
		bx_o.setSmall(dir, d_lo+1);
	}
	if( b_hi == (d_hi + 1) && (*hibc)[dir] == amrex::LinOpBCType::Dirichlet )
	{
		bx_o.setBig(dir, d_hi);
	}
	return bx_o;
}

void incFSI::setBCs(BCArray& bcarr, const char* s)
{
    amrex::ParmParse pp;

    amrex::Vector<std::string> bc(AMREX_SPACEDIM);
    pp.getarr(s, bc);

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        if (bc[i] == "dirichlet" || bc[i] == "Dirichlet")
        {
            bcarr[i] = amrex::LinOpBCType::Dirichlet;
        }
        else if (bc[i] == "neumann" || bc[i] == "Neumann")
        {
            bcarr[i] = amrex::LinOpBCType::Neumann;
        }
        else
        {
            amrex::Abort("setBcs: invalid bc specified");
        }
    }
}

void incFSI::setIntBCs(std::vector<amrex::LinOpBCType> &intbc, const char *s)
{
    for (size_t i = 0; i < N_IF; i++)
    {
        amrex::ParmParse pp(IF_names[i]);
        std::string bc;
        pp.get(s, bc);

        if (bc == "dirichlet" || bc == "Dirichlet")
        {
            intbc[i] = amrex::LinOpBCType::Dirichlet;
        }
        else if (bc == "neumann" || bc == "Neumann")
        {
            intbc[i] = amrex::LinOpBCType::Neumann;
        }
        else
        {
            amrex::Abort("setBcs: invalid bc specified");
        }
    }
}

void incFSI::InitializeRKcoeffs()
{
    if(RKOrder == 1)
    {
        rk_a[0][0] = 1.0;
        rk_c[0] = 1.0;
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
    else if(RKOrder == 2)
    {
        rk_a[0][0] = 1.0;
        rk_a[1][0] = 0.5;
        rk_a[1][1] = 0.5;
        rk_c[0] = 1.0;
        rk_c[1] = 1.0;
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
    else if(RKOrder == 3)
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
    }
    else if(RKOrder == 4)
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
}

} /*End namespace mycode */

