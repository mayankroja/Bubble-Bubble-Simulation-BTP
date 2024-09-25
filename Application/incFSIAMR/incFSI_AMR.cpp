#include "incFSI.H"
#include <AMReX_ParmParse.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_MultiFabUtil.H>

namespace mycode
{

void incFSI::Regrid()
{
    //! Regrid
    if (max_level > 0 && regrid_int > 0)  // We may need to regrid
    {
        if (Iter % regrid_int == 0)
        {
            regrid(0, Time);
            //amrex::Print()<<"Regrid"<<'\n';
            if(N_IF > 0)
            {
                ComputeBubbleVolume(false);
                DetectBubbles(true); 
                for (int lev = 0; lev <= finest_level; lev++)
                {
                    //FillGhost(lev); 
                    for (auto &&solid : interfaces[lev])
                    {   
                        solid->TubeIdentification();
                    }
                    mask_[lev]->GhostCellIdentfication();
		    if(lev == finest_level)
		    {
                        ComputeCutFaceVel(lev);
                        ComputeCutFacePressure(lev);
		    }
                }
            }
	    AverageDown();
        }
    }
}

void incFSI::AverageDown()
{
    for (int lev = finest_level; lev > 0; lev--)
    {
        amrex::Array<const amrex::MultiFab *, AMREX_SPACEDIM> vel = {&xvel[lev], &yvel[lev]};
        amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> velc = {&xvel[lev-1], &yvel[lev-1]};
        amrex::average_down_faces(vel, velc, ref_ratio[lev - 1], geom[lev - 1]);
	amrex::average_down(Pressure[lev], Pressure[lev - 1], 0, 1, ref_ratio[lev-1][0]);
        if(TempField)amrex::average_down(Theta[lev], Theta[lev - 1], 0, 1, ref_ratio[lev-1][0]);
        if(PhaseField) amrex::average_down(Phi[lev], Phi[lev - 1], 0, 1, ref_ratio[lev-1][0]);
	if(DamageModel)
	{
            for(int iscalar = 0;iscalar < nscalar; iscalar++)
            {
	        amrex::average_down(Scalars[lev][iscalar], Scalars[lev - 1][iscalar], 0, 1, ref_ratio[lev-1][0]);
	    }
        }
    }
}

void incFSI::FillPatchAroundBox()
{
    if(max_level > 0)
    {
        for(int lev = 0;lev <= finest_level;lev++)
        {
            //! create temp field data
            amrex::BoxArray xba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(0));
            amrex::BoxArray yba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(1));
            amrex::MultiFab tmp_xvel(xba, dmap[lev], 1, Nghost);
            amrex::MultiFab tmp_yvel(yba, dmap[lev], 1, Nghost);
            amrex::MultiFab tmp_Pressure(grids[lev], dmap[lev], 1, Nghost);

            tmp_xvel.setVal(0.0);
            tmp_yvel.setVal(0.0);
            tmp_Pressure.setVal(0.0);

            amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> tmp_vel = {&tmp_xvel, &tmp_yvel};
            amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> cvel;
            amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> fvel;

            cvel.resize(1);
            fvel.resize(1);
            cvel[0] = {&xvel[lev - 1], &yvel[lev - 1]};
            fvel[0] = {&xvel[lev], &yvel[lev]};
            {
                if(divergence_free_interpolation)
                {
                    FillPatch(lev, Time, tmp_vel, cvel, {Time}, fvel, {Time}, 0, 1);

                    std::swap(*tmp_vel[0], xvel[lev]);
                    std::swap(*tmp_vel[1], yvel[lev]);

                    XVelBoundaryConditions(lev, Time, xvel[lev]);
                    YVelBoundaryConditions(lev, Time, yvel[lev]);
                }
                else
                {
                    FillPatch(lev, Time, tmp_xvel, {&xvel[lev - 1]}, {Time}, {&xvel[lev]}, {Time}, 0, 1);
                    XVelBoundaryConditions(lev, Time, tmp_xvel);
                    std::swap(tmp_xvel, xvel[lev]);

                    FillPatch(lev, Time, tmp_yvel, {&yvel[lev - 1]}, {Time}, {&yvel[lev]}, {Time}, 0, 1);
                    YVelBoundaryConditions(lev, Time, tmp_yvel);
                    std::swap(tmp_yvel, yvel[lev]);
                }
                FillPatch(lev, Time, tmp_Pressure, {&Pressure[lev - 1]}, {Time}, {&Pressure[lev]}, {Time}, 0, 1);
                PressureBoundaryConditions(lev, Time, tmp_Pressure);
                std::swap(tmp_Pressure, Pressure[lev]);

                if(TempField)
                {
                    amrex::MultiFab tmp_Theta(grids[lev], dmap[lev], 1, Nghost);
                    tmp_Theta.setVal(0.0);
                    FillPatch(lev, Time, tmp_Theta, {&Theta[lev - 1]}, {Time}, {&Theta[lev]}, {Time}, 0, 1);
                    TemperatureBoundaryConditions(lev, Time, tmp_Theta);
                    std::swap(tmp_Theta, Theta[lev]);
                }
                if(PhaseField)
                {
                    amrex::MultiFab tmp_Phi(grids[lev], dmap[lev], 1, Nghost);
                    tmp_Phi.setVal(0.0);
                    FillPatch(lev, Time, tmp_Phi, {&Phi[lev - 1]}, {Time}, {&Phi[lev]}, {Time}, 0, 1);
                    PhiBoundaryConditions(lev, Time, tmp_Phi);
                    std::swap(tmp_Phi, Phi[lev]);
                }
                if(DamageModel)
                {
		    for(int iscalar = 0; iscalar < nscalar;iscalar++)
	            {
                        amrex::MultiFab tmp_Scalar(grids[lev], dmap[lev], 1, Nghost);
                        tmp_Scalar.setVal(0.0);
                        FillPatchScalars(lev, iscalar, Time, tmp_Scalar, {&Scalars[lev - 1][iscalar]}, {Time}, {&Scalars[lev][iscalar]}, {Time}, 0, 1);
                        ScalarBoundaryConditions(lev, iscalar, Time, tmp_Scalar);
                        std::swap(tmp_Scalar, Scalars[lev][iscalar]);
		    }
                }
            }
        }
    }     
}

void incFSI::MakeNewLevelFromScratch 
(
    int lev, 
    amrex::Real time, 
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm
)
{
    //! Allocate memory for this level and initialise data
    AllocateMemory(lev, ba, dm);
    if(N_IF > 0) MakeInterface(lev, ba, dm);

    //! initialise data on this level
    SetInitialFlowField(lev);
    if(lev > 0)
    {
        amrex::BoxArray xba = amrex::convert(ba, amrex::IntVect::TheDimensionVector(0));
        amrex::BoxArray yba = amrex::convert(ba, amrex::IntVect::TheDimensionVector(1));

        amrex::MultiFab tmp_Pressure(ba, dm, 1, Nghost);
        tmp_Pressure.setVal(0.0);
        FillPatch(lev, time, tmp_Pressure,{&Pressure[lev - 1]}, {0.0}, {&Pressure[lev]}, {0.0}, 0, 1);
        std::swap(tmp_Pressure, Pressure[lev]);

        if(TempField)
        {
            amrex::MultiFab tmp_Theta(ba, dm, 1, Nghost);
            tmp_Theta.setVal(0.0);
            FillPatch(lev, time, tmp_Theta,{&Theta[lev - 1]}, {0.0}, {&Theta[lev]}, {0.0}, 0, 1);
            std::swap(tmp_Theta, Theta[lev]);
        }
        if(PhaseField)
        {   
            amrex::MultiFab tmp_Phi(ba, dm, 1, Nghost);
            tmp_Phi.setVal(0.0);
            FillPatch(lev, time, tmp_Phi,{&Phi[lev - 1]}, {0.0}, {&Phi[lev]}, {0.0}, 0, 1);
            std::swap(tmp_Phi, Phi[lev]);
        }
	if(DamageModel)
        {
	    for(int iscalar = 0; iscalar < nscalar ; iscalar++)
            {
                amrex::MultiFab tmp_scalar(ba, dm, 1, Nghost);
                tmp_scalar.setVal(0.0);
                FillPatchScalars(lev, iscalar, time, tmp_scalar,{&Scalars[lev - 1][iscalar]}, {0.0}, {&Scalars[lev][iscalar]}, {0.0}, 0, 1);
                std::swap(tmp_scalar, Scalars[lev][iscalar]);
	    }
        }

        amrex::MultiFab tmp_xvel(xba, dm, 1, Nghost);
        tmp_xvel.setVal(0.0);

        amrex::MultiFab tmp_yvel(yba, dm, 1, Nghost);
        tmp_yvel.setVal(0.0);


        //amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> tmp_vel = {&tmp_xvel, &tmp_yvel};
        amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> tmp_vel;
  
        tmp_vel[0] = &tmp_xvel;
        tmp_vel[1] = &tmp_yvel;

        amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> cvel;
        amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> fvel;

        cvel.resize(1);
        fvel.resize(1);
        cvel[0] = {&xvel[lev - 1], &yvel[lev - 1]};
        fvel[0] = {&xvel[lev], &yvel[lev]};

        if(divergence_free_interpolation )
        {
            FillPatch(lev, time, tmp_vel, cvel, {0.0}, fvel, {0.0}, 0, 1);

            std::swap(*tmp_vel[0], xvel[lev]);
            std::swap(*tmp_vel[1], yvel[lev]);

            XVelBoundaryConditions(lev, time, xvel[lev]);
            YVelBoundaryConditions(lev, time, yvel[lev]);
        }
        else
        {
            FillPatch(lev, time, tmp_xvel, {&xvel[lev - 1]}, {0.0}, {&xvel[lev]}, {0.0}, 0, 1);
            std::swap(tmp_xvel, xvel[lev]);

            FillPatch(lev, time, tmp_yvel, {&yvel[lev - 1]}, {0.0}, {&yvel[lev]}, {0.0}, 0, 1);
            std::swap(tmp_yvel, yvel[lev]);
        }
    }
    XVelBoundaryConditions(lev, time, xvel[lev]);
    YVelBoundaryConditions(lev, time, yvel[lev]);
    PressureBoundaryConditions(lev, time, Pressure[lev]);
    if(TempField) TemperatureBoundaryConditions(lev, time, Theta[lev]);
    if(PhaseField) PhiBoundaryConditions(lev, time, Phi[lev]);
    /*if(N_IF > 0)
    {
        ComputeCutFaceVel(lev);
        ComputeCutFacePressure(lev);
        ComputeAvgIntProp(lev);
    }*/
}    

void incFSI::MakeNewLevelFromCoarse 
(
    int lev, 
    amrex::Real time, 
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm
)
{
    amrex::Print()<< "lev : " << lev <<'\t'<<"MakeNewLevelFromCoarse"<<'\n';
    if(N_IF > 0) MakeNewLevelFromCoarse_interface(lev, time, ba, dm);
    //! Allocate memory for this level
    AllocateMemory(lev, ba, dm);

    //! do a fill boundary on coarse data
    {
        xvel[lev - 1].FillBoundary();
        yvel[lev - 1].FillBoundary();
        xvel_old[lev - 1].FillBoundary();
        yvel_old[lev - 1].FillBoundary();
        if(TempField)
        {
            Theta[lev - 1].FillBoundary();
            Theta_old[lev - 1].FillBoundary();
        }
        if(PhaseField)
        {
            Phi[lev - 1].FillBoundary();
            Phi_old[lev - 1].FillBoundary();
        }

    }

    //! initialise data from coarse level
    //! the face_div_free mapper has some problem as it gives nans
    // amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> Vel,  Vel_c;
    // amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> Velold, Velold_c;
    // for (int i = 0; i < AMREX_SPACEDIM; i++)
    // {
    //     if (i == 0)
    //     {
    //         Vel[i] = &(xvel[lev]);
    //         Vel_c[i] = &(xvel[lev-1]);
    //         Velold[i] = &(xvel_old[lev]);
    //         Velold_c[i] = &(xvel_old[lev-1]);
    //     }
    //     else if (i == 1)
    //     {
    //         Vel[i] = &(yvel[lev]);
    //         Vel_c[i] = &(yvel[lev-1]);
    //         Velold[i] = &(yvel_old[lev]);
    //         Velold_c[i] = &(yvel_old[lev-1]);            
    //     }        
    // }
    // FillCoarsePatch(lev, time, Vel, Vel_c, 0, 1);
    // FillCoarsePatch(lev, time, Velold, Velold_c, 0, 1);

    //! try linear mapper
    FillCoarsePatch(lev, time, xvel[lev], xvel[lev - 1], 0, 1);
    FillCoarsePatch(lev, time, yvel[lev], yvel[lev - 1], 0, 1);
    FillCoarsePatch(lev, time, xvel_old[lev], xvel_old[lev - 1], 0, 1);
    FillCoarsePatch(lev, time, yvel_old[lev], yvel_old[lev - 1], 0, 1);

    ComputeCollocatedVelocityField(lev);

    FillCoarsePatch(lev, time, Pressure[lev], Pressure[lev-1], 0, 1);
    FillCoarsePatch(lev, time, Pressure_old[lev], Pressure_old[lev-1], 0, 1);
    FillCoarsePatch(lev, time, Pstar[lev], Pstar[lev-1], 0, 1);
    

    if(TempField) 
    {
        FillCoarsePatch(lev, time, Theta[lev], Theta[lev-1], 0, 1);
        FillCoarsePatch(lev, time, Theta_old[lev], Theta_old[lev-1], 0, 1);
    }
    if(PhaseField) 
    {  
        FillCoarsePatch(lev, time, Phi[lev], Phi[lev-1], 0, 1);
        FillCoarsePatch(lev, time, Phi_old[lev], Phi_old[lev-1], 0, 1);
    }

    // TODO : if physcbc of old data are required then?
    //! fill physbcs
    XVelBoundaryConditions(lev, time, xvel[lev]);
    YVelBoundaryConditions(lev, time, yvel[lev]);
    PressureBoundaryConditions(lev, time, Pressure[lev]);
    if(TempField) TemperatureBoundaryConditions(lev, time, Theta[lev]);
    if(PhaseField) PhiBoundaryConditions(lev, time, Pressure[lev]);
}

void incFSI::RemakeLevel 
(
    int lev, 
    amrex::Real time, 
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm
)
{
    //amrex::Print()<< "lev : " << lev <<'\t'<<"RemakeLevel"<<'\n';
    if(N_IF > 0) RemakeInterface(lev, time, ba, dm);
    //! create temp field data
    amrex::BoxArray xba = amrex::convert(ba, amrex::IntVect::TheDimensionVector(0));
    amrex::BoxArray yba = amrex::convert(ba, amrex::IntVect::TheDimensionVector(1));

    grids_xvel[lev] = xba;
    grids_yvel[lev] = yba;

    amrex::MultiFab tmp_xvel(xba, dm, 1, Nghost);
    amrex::MultiFab tmp_yvel(yba, dm, 1, Nghost);
    amrex::MultiFab tmp_xvel_old(xba, dm, 1, Nghost);
    amrex::MultiFab tmp_yvel_old(yba, dm, 1, Nghost);

    amrex::MultiFab tmp_U(ba, dm, AMREX_SPACEDIM, Nghost);
    amrex::MultiFab tmp_Pressure(ba, dm, 1, Nghost);
    amrex::MultiFab tmp_Pressure_old(ba, dm, 1, Nghost);
    amrex::MultiFab tmp_Pstar(ba, dm, 1, Nghost);
    amrex::MultiFab tmp_RHS(ba, dm, 1, 0);
    amrex::MultiFab tmp_Umag(ba, dm, 1, 0);
    amrex::MultiFab tmp_Vorticity(ba, dm, 1, 0);
    amrex::MultiFab tmp_Theta(ba, dm, 1, Nghost);
    amrex::MultiFab tmp_Theta_old(ba, dm, 1, Nghost);
    amrex::MultiFab tmp_Phi(ba, dm, 1, Nghost);
    amrex::MultiFab tmp_Phi_old(ba, dm, 1, Nghost);

    //! fill from old current level and coarse level
    //! need to pass the current data also
    amrex::Vector<amrex::MultiFab *> cP(1), cPstar(1), cTheta(1), cPhi(1);
    amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> tmp_vel = {&tmp_xvel, &tmp_yvel};
    amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> tmp_vel_old = {&tmp_xvel_old, &tmp_yvel_old};

    amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> cvel;
    amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> fvel;
    amrex::Vector<amrex::Real> ctime(2);
    amrex::Vector<amrex::Real> ftime(2);

    amrex::Vector<amrex::MultiFab *> cxvel(2);
    amrex::Vector<amrex::MultiFab *> cyvel(2);

    amrex::Vector<amrex::MultiFab *> fxvel(2);
    amrex::Vector<amrex::MultiFab *> fyvel(2);

    cvel.resize(2);
    if (lev > 0)
    {
        cvel[0] = {&xvel_old[lev - 1], &yvel_old[lev - 1]}; // old data
        cvel[1] = {&xvel[lev - 1], &yvel[lev - 1]};         // new data
        ctime[0] = t_old;
        ctime[1] = t_new;

        cP[0] = &Pressure[lev - 1];
        cPstar[0] = &Pstar[lev - 1];
        if(TempField) cTheta[0] = &Theta[lev - 1];
        if(PhaseField) cPhi[0] = &Phi[lev - 1];

        cxvel[0] = &xvel_old[lev - 1];
        cxvel[1] = &xvel[lev - 1];
        cyvel[0] = &yvel_old[lev - 1];
        cyvel[1] = &yvel[lev - 1];
    }
    fvel.resize(2);
    fvel[0] = {&xvel_old[lev], &yvel_old[lev]}; // old data
    fvel[1] = {&xvel[lev], &yvel[lev]};         // new data
    ftime[0] = t_old;
    ftime[1] = t_new;

    fxvel[0] = &xvel_old[lev];
    fxvel[1] = &xvel[lev];
    fyvel[0] = &yvel_old[lev];
    fyvel[1] = &yvel[lev];

    if(divergence_free_interpolation)
    {
        FillPatch(lev, time, tmp_vel, cvel, ctime, fvel, ftime, 0, 1);
        amrex::MultiFab::Copy(tmp_xvel, *tmp_vel[0], 0, 0, 1, Nghost);
        amrex::MultiFab::Copy(tmp_yvel, *tmp_vel[1], 0, 0, 1, Nghost);

        //! for old time data not need of interpolation in time, so only use old time data
        FillPatch(lev, time, tmp_vel_old, {cvel[0]}, {ctime[0]}, {fvel[0]}, {ftime[0]}, 0, 1);
        amrex::MultiFab::Copy(tmp_xvel_old, *tmp_vel_old[0], 0, 0, 1, Nghost);
        amrex::MultiFab::Copy(tmp_yvel_old, *tmp_vel_old[1], 0, 0, 1, Nghost);        
    }
    else
    {
        FillPatch(lev, time, tmp_xvel, cxvel, ctime, fxvel, ftime, 0, 1);
        FillPatch(lev, time, tmp_yvel, cyvel, ctime, fyvel, ftime, 0, 1);

        FillPatch(lev, time, tmp_xvel_old, {cxvel[0]}, {ctime[0]}, {&xvel_old[lev]}, {ftime[0]}, 0, 1);
        FillPatch(lev, Time, tmp_yvel_old, {cyvel[0]}, {ctime[0]}, {&yvel_old[lev]}, {ftime[0]}, 0, 1);
    }



    //! no need of time interpolation of following
    FillPatch(lev, time, tmp_Pressure, cP, {ctime[0]}, {&Pressure[lev]}, {ftime[0]}, 0, 1);
    FillPatch(lev, time, tmp_Pstar, cPstar, {ctime[0]}, {&Pstar[lev]}, {ftime[0]}, 0, 1);
    if(TempField) 
    {
        FillPatch(lev, time, tmp_Theta, cTheta, {ctime[0]}, {&Theta[lev]}, {ftime[0]}, 0, 1);
        FillPatch(lev, time, tmp_Theta_old, {&Theta_old[lev - 1]}, {ctime[0]}, {&Theta_old[lev]}, {ftime[0]}, 0, 1);
    }
    if(PhaseField) 
    {
        FillPatch(lev, time, tmp_Phi, cPhi, {ctime[0]}, {&Phi[lev]}, {ftime[0]}, 0, 1);
        FillPatch(lev, time, tmp_Phi_old, {&Phi_old[lev - 1]}, {ctime[0]}, {&Phi[lev]}, {ftime[0]}, 0, 1);
    }

    //! these are not required to interpolate from coarse
    tmp_U.setVal(0.0);
    tmp_RHS.setVal(0.0);
    tmp_Umag.setVal(0.0);
    tmp_Vorticity.setVal(0.0);

    //! now swap newly created level data with current level data
    std::swap(tmp_xvel, xvel[lev]);
    std::swap(tmp_yvel, yvel[lev]);
    std::swap(tmp_xvel_old, xvel_old[lev]);
    std::swap(tmp_yvel_old, yvel_old[lev]);
    std::swap(tmp_Pressure_old, Pressure_old[lev]);
    std::swap(tmp_U, U[lev]);
    std::swap(tmp_Pressure, Pressure[lev]);
    std::swap(tmp_Pstar, Pstar[lev]);
    std::swap(tmp_RHS, RHS[lev]);
    std::swap(tmp_Umag, Umag[lev]);
    std::swap(tmp_Vorticity, Vorticity[lev]);
    if(TempField) 
    {
        std::swap(tmp_Theta, Theta[lev]);
        std::swap(tmp_Theta_old, Theta_old[lev]);
    }
    if(PhaseField) 
    {
        std::swap(tmp_Phi, Phi[lev]);
        std::swap(tmp_Phi_old, Phi_old[lev]);
    }

    if(RKOrder == 2)
    {
        FRK_xvel[lev].clear();
        FRK_yvel[lev].clear();
        FRK_p[lev].clear();
        if(TempField) FRK_Theta[lev].clear();
        if(PhaseField) FRK_Phi[lev].clear();
	if(DamageModel)
	{
	    for(int iscalar = 0; iscalar < nscalar; iscalar++)
                FRK_Scalars[lev][iscalar].clear();
	}

        FRK_xvel[lev].define(xba, dm, 1, Nghost);
        FRK_yvel[lev].define(yba, dm, 1, Nghost);
        FRK_p[lev].define(ba, dm, 1, Nghost);
        if(TempField) FRK_Theta[lev].define(ba, dm, 1, Nghost);
        if(PhaseField) FRK_Phi[lev].define(ba, dm, 1, Nghost);
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
                FRK_Scalars[lev][iscalar].define(ba, dm, 1, Nghost);
        }

    }

    if(RKOrder >= 3)
    {
        RK1_xvel[lev].clear();
        RK1_yvel[lev].clear();
        RK1_p[lev].clear();
        if(TempField) RK1_Theta[lev].clear();
        if(PhaseField) RK1_Phi[lev].clear();
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
                RK1_Scalars[lev][iscalar].clear();
        }

        RK1_xvel[lev].define(xba, dm, 1, Nghost);
        RK1_yvel[lev].define(yba, dm, 1, Nghost);
        RK1_p[lev].define(ba, dm, 1, Nghost);
        if(TempField) RK1_Theta[lev].define(ba, dm, 1, Nghost);
        if(PhaseField) RK1_Phi[lev].define(ba, dm, 1, Nghost);
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
                RK1_Scalars[lev][iscalar].define(ba, dm, 1, Nghost);
        }

        RK2_xvel[lev].clear();
        RK2_yvel[lev].clear();
        RK2_p[lev].clear();
        if(TempField) RK2_Theta[lev].clear();
        if(PhaseField) RK2_Phi[lev].clear();
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
                RK2_Scalars[lev][iscalar].clear();
        }

        RK2_xvel[lev].define(xba, dm, 1, Nghost);
        RK2_yvel[lev].define(yba, dm, 1, Nghost);
        RK2_p[lev].define(ba, dm, 1, Nghost);
        if(TempField) RK2_Theta[lev].define(ba, dm, 1, Nghost);
        if(PhaseField) RK2_Phi[lev].define(ba, dm, 1, Nghost);
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
                RK2_Scalars[lev][iscalar].define(ba, dm, 1, Nghost);
        }

        RK3_xvel[lev].clear();
        RK3_yvel[lev].clear();
        RK3_p[lev].clear();
        if(TempField) RK3_Theta[lev].clear();
        if(PhaseField) RK3_Phi[lev].clear();
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
                RK3_Scalars[lev][iscalar].clear();
        }

        RK3_xvel[lev].define(xba, dm, 1, Nghost);
        RK3_yvel[lev].define(yba, dm, 1, Nghost);
        RK3_p[lev].define(ba, dm, 1, Nghost);
        if(TempField) RK3_Theta[lev].define(ba, dm, 1, Nghost);
        if(PhaseField) RK3_Phi[lev].define(ba, dm, 1, Nghost);
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
                RK3_Scalars[lev][iscalar].define(ba, dm, 1, Nghost);
        }

        Src_xvel[lev].clear();
        Src_yvel[lev].clear();
        Src_p[lev].clear();
        if(TempField) Src_Theta[lev].clear();
        if(PhaseField) Src_Phi[lev].clear();
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
                Src_Scalars[lev][iscalar].clear();
        }

        Src_xvel[lev].define(xba, dm, 1, Nghost);
        Src_yvel[lev].define(yba, dm, 1, Nghost);
        Src_p[lev].define(ba, dm, 1, Nghost);
        if(TempField) Src_Theta[lev].define(ba, dm, 1, Nghost);
        if(PhaseField) Src_Phi[lev].define(ba, dm, 1, Nghost);
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
                Src_Scalars[lev][iscalar].define(ba, dm, 1, Nghost);
        }
    }

    if(DamageModel)
    {
	for(int iscalar = 0; iscalar < nscalar; iscalar++)
	{
            amrex::MultiFab tmp_Scalar(ba, dm, 1, Nghost);
            amrex::MultiFab tmp_Scalar_old(ba, dm, 1, Nghost);
	    amrex::Vector<amrex::MultiFab *> cScalar(1);
	    if(lev > 0)
	        cScalar[0] = &Scalars[lev - 1][iscalar];

            FillPatchScalars(lev, iscalar, time, tmp_Scalar, cScalar, {ctime[0]}, {&Scalars[lev][iscalar]}, {ftime[0]}, 0, 1);
            FillPatchScalars(lev, iscalar, time, tmp_Scalar_old, {&Scalars_old[lev - 1][iscalar]}, {ctime[0]}, {&Scalars_old[lev][iscalar]}, {ftime[0]}, 0, 1);
            std::swap(tmp_Scalar, Scalars[lev][iscalar]);
            std::swap(tmp_Scalar_old, Scalars_old[lev][iscalar]);
	}
    }

    if (StopAtSteady)
    {
        amrex::MultiFab tmp_UmagOld(ba, dm, 1, 0);
        tmp_UmagOld.setVal(0.0);
        std::swap(tmp_UmagOld, UmagOld[lev]);
    }

    // TODO : if physcbc of old data are required then?
    //! fill physbcs
    XVelBoundaryConditions(lev, time, xvel[lev]);
    YVelBoundaryConditions(lev, time, yvel[lev]);
    PressureBoundaryConditions(lev, time, Pressure[lev]);
    if(TempField) TemperatureBoundaryConditions(lev, time, Theta[lev]);
    if(PhaseField) PhiBoundaryConditions(lev, time, Phi[lev]);
    if(DamageModel)
    {
        for(int iscalar = 0; iscalar < nscalar; iscalar++)
	    ScalarBoundaryConditions(lev,iscalar, time, Scalars[lev][iscalar]);
    }
}

void incFSI::ClearLevel(int lev)
{
    //! clear field data on this level
    amrex::Print() << "lev : " << lev << "\tClearLevel\n";
    xvel[lev].clear();
    yvel[lev].clear();
    U[lev].clear();
    Pressure[lev].clear();
    RHS[lev].clear();
    xvel_old[lev].clear();
    yvel_old[lev].clear();
    Pressure_old[lev].clear();
    Pstar[lev].clear();
    Umag[lev].clear();
    if(TempField)
    {
        Theta[lev].clear();
        Theta_old[lev].clear();
    }
    if(PhaseField)
    {
        Phi[lev].clear();
        Phi_old[lev].clear();
    }
    if(DamageModel)
    {
        for(int iscalar = 0; iscalar < nscalar; iscalar++)
	{
	    Scalars[lev][iscalar].clear();
	    Scalars_old[lev][iscalar].clear();
	}
    }
    if(RKOrder == 2)
    {
        FRK_xvel[lev].clear();
        FRK_yvel[lev].clear();
        FRK_p[lev].clear();
        if(TempField) FRK_Theta.clear();
        if(PhaseField) FRK_Phi.clear();
	if(DamageModel)
        {
	    for(int iscalar = 0; iscalar < nscalar; iscalar++)
                FRK_Scalars[lev][iscalar].clear();   
	}
    }
    if(RKOrder >= 3)
    {
        RK1_xvel[lev].clear();
        RK1_yvel[lev].clear();
        RK1_p[lev].clear();
        if(TempField) RK1_Theta.clear();
        if(PhaseField) RK1_Phi.clear();
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
                RK1_Scalars[lev][iscalar].clear();
        }
        
        RK2_xvel[lev].clear();
        RK2_yvel[lev].clear();
        RK2_p[lev].clear();
        if(TempField) RK2_Theta.clear();
        if(PhaseField) RK2_Phi.clear();    
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
                RK2_Scalars[lev][iscalar].clear();
        }

        RK3_xvel[lev].clear();
        RK3_yvel[lev].clear();
        RK3_p[lev].clear();
        if(TempField) RK3_Theta.clear();
        if(PhaseField) RK3_Phi.clear();
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
                RK3_Scalars[lev][iscalar].clear();
        }
        
        Src_xvel[lev].clear();
        Src_yvel[lev].clear();
        Src_p[lev].clear();
        if(TempField) Src_Theta.clear();
        if(PhaseField) Src_Phi.clear();              
        if(DamageModel)
        {
            for(int iscalar = 0; iscalar < nscalar; iscalar++)
                Src_Scalars[lev][iscalar].clear();
        }
    }    
    
    if (StopAtSteady)
        UmagOld[lev].clear();
}

void incFSI::ErrorEst 
(
    int lev, 
    amrex::TagBoxArray& tags, 
    amrex::Real time, 
    int ngrow
)
{
    //amrex::Print() << "lev : " << lev << "\tErrorEst\n";
    static bool first = true;
    static amrex::Vector<amrex::Real> vorterr;
    //! Tag tube around interface
    
    if (tag_interface)
    {
        const int   tagval = amrex::TagBox::SET;
        const amrex::Real *dx = geom[lev].CellSize();
        const amrex::Real *prob_lo = geom[lev].ProbLo();
        const amrex::Real *prob_hi = geom[lev].ProbHi();
        const amrex::Box &domain = geom[lev].Domain();

        const amrex::MultiFab& state = Pressure[lev];
        for (auto &&solid : interfaces[lev])
        {
            amrex::MultiFab &psi = solid->Psi();
            amrex::iMultiFab &objmask = solid->Mask();
            const amrex::MultiFab& state = Pressure[lev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
            {
                for (amrex::MFIter mfi(psi,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx  = mfi.tilebox();
                    amrex::Array4<amrex::Real> const &Psi = psi.array(mfi);
                    amrex::Array4<int> const &obj_mask = objmask.array(mfi);
                    const auto tagfab  = tags.array(mfi);
                    //! mark the intersection box for tagging
                    amrex::ParallelFor(bx,
                    [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
			//if(std::abs(Psi(i ,j ,k)) < 20*0.5*(dx[0] + dx[1]))
                        if(obj_mask(i, j, k) > 0 || Psi(i,j,k) > 0.0)
                            tagfab(i, j, k) = tagval;
                    });
                }
            }
	}
    } 
    if(PhaseField)
    {
        const int   tagval = amrex::TagBox::SET;
	const amrex::Real *dx = geom[lev].CellSize();
        const amrex::MultiFab& state = Phi[lev];
#ifdef AMREX_USE_OMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
        {
            for (amrex::MFIter mfi(state,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx  = mfi.tilebox();
                amrex::Array4<amrex::Real const> const &phi = state.array(mfi);
                const auto tagfab  = tags.array(mfi);
                //! mark the intersection box for tagging
                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    //if(obj_mask(i, j, k) > 0)
                    //if(phi(i,j,k) < 0.9 && phi(i, j, k) > 0.1)
                    //    tagfab(i, j, k) = tagval;
		    //grad_phi = std::hypot(std::abs(phi(i+1,j,k) - phi(i-1,j,k)), std::abs(phi(i,j+1,k) - phi(i,j-1,k)));
		    amrex::Real grad_phi = (1.0/dx[0])*std::hypot(0.5*(phi(i + 1, j, k) - phi(i - 1, j, k)), 0.5*(phi(i, j + 1, k) - phi(i, j - 1, k)));
		    //if(std::abs(phi(i+1,j,k) - phi(i-1,j,k)) > 0.1 || std::abs(phi(i,j+1,k) - phi(i,j-1,k)) > 0.1 )
		    if(grad_phi > 0.1)
		        tagfab(i, j, k) = tagval;
                    if(phi(i, j, k) < 0)
                        tagfab(i, j, k) = tagval;
                });
            }
        }
    }
    // only do this during the first call to ErrorEst
    if (first)
    {
        first = false;
        // read in an array of "phierr", which is the tagging threshold
        // tag values of "Vorticity" which are greater than vorterr
        // for that particular level
        amrex::ParmParse pp;
        if (tag_vort)
        {
            int n = pp.countval("vorterr");
            if (n > 0) {
                pp.getarr("vorterr", vorterr, 0, n);
            }
        }        
    }

    if (tag_vort && lev >= vorterr.size()) return;

    //! vorticity based tagging
    if(tag_vort)
    {
        //    const int clearval = TagBox::CLEAR;
        const int   tagval = amrex::TagBox::SET;

        const amrex::MultiFab& state = Vorticity[lev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
        {

            for (amrex::MFIter mfi(state,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx  = mfi.tilebox();
                const auto statefab = state.array(mfi);
                const auto tagfab  = tags.array(mfi);
                amrex::Real vorterror = vorterr[lev];

                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    //! mark all the cells in interface
                    if (fabs(statefab(i,j,k)) >= vorterror)
                    {
                        tagfab(i, j, k) = tagval;
                    }
                });
            }
        }
    }
    //! multiple regions
    if (tag_region)
    {
        const int   tagval = amrex::TagBox::SET;
        const amrex::MultiFab& state = Pressure[lev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
        {
            for (amrex::MFIter mfi(state,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx  = mfi.tilebox();
                const auto statefab = state.array(mfi);
                const auto tagfab  = tags.array(mfi);
                
                //! iterate over tag_region boxes
                for (auto &&box : boxlist[lev])
                {
                    //! check if the current box 'bx' has valid intersection with 
                    //! box in the boxlist of tag_region boxes
                    amrex::Box isect = bx & box;
                    if (isect.ok())
                    {
                        //! mark the intersection box for tagging
                        amrex::ParallelFor(isect,
                        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                        {
                            tagfab(i, j, k) = tagval;
                        });
                    }
                }
            }
        }
    } 
}    

// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void incFSI::FillPatch 
(
    int lev, 
    amrex::Real time, 
    amrex::MultiFab& mf,                        // dest mf
    const amrex::Vector<amrex::MultiFab*>& cmf, // coarse mf
    const amrex::Vector<amrex::Real>& ct,       // coarse time
    const amrex::Vector<amrex::MultiFab*>& fmf, // fine mf
    const amrex::Vector<amrex::Real>& ft,       // fine time
    int icomp, int ncomp
)
{
    amrex::Interpolater *mapper;
    amrex::Vector<amrex::BCRec> bcs(1);
    if (mf.is_cell_centered())
    {
        mapper = &amrex::cell_cons_interp;

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            bcs[0].setLo(idim, amrex::BCType::int_dir);
            bcs[0].setHi(idim, amrex::BCType::int_dir);
            //bcs[0].setLo(idim, amrex::BCType::reflect_even);
            //bcs[0].setHi(idim, amrex::BCType::reflect_even);
        }


        if (lev == 0)
        {
            amrex::PhysBCFunctNoOp physbc;
            //amrex::GpuBndryFuncFab<MyExtBCFill> gpu_bndry_func;
            //amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill>> physbc(geom[lev],bcs,gpu_bndry_func);
            amrex::FillPatchSingleLevel(mf, time, fmf, ft, 0, icomp, ncomp,
                                            geom[lev], physbc, 0);
        }
        else
        {
            //amrex::GpuBndryFuncFab<MyExtBCFill> gpu_bndry_func;
            //amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill>> cphysbc(geom[lev-1],bcs,gpu_bndry_func);
            //amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill>> fphysbc(geom[lev],bcs,gpu_bndry_func);
            amrex::PhysBCFunctNoOp cphysbc, fphysbc;
            amrex::FillPatchTwoLevels(mf, time, cmf, ct, fmf, ft,
                                          0, icomp, ncomp, geom[lev-1], geom[lev],
                                          cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                          mapper, bcs, 0);            
        }

    }
    else if 
    (
        mf.ixType() == amrex::IndexType(amrex::IntVect::TheDimensionVector(0)) || 
        mf.ixType() == amrex::IndexType(amrex::IntVect::TheDimensionVector(1))
    )
    {
        mapper = &amrex::face_linear_interp;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            bcs[0].setLo(idim, amrex::BCType::int_dir);
            bcs[0].setHi(idim, amrex::BCType::int_dir);
            //bcs[0].setLo(idim, amrex::BCType::foextrap);
            //bcs[0].setHi(idim, amrex::BCType::foextrap);
            //bcs[0].setLo(idim, amrex::BCType::ext_dir);
            //bcs[0].setHi(idim, amrex::BCType::ext_dir);
        }
        if (lev == 0)
        {
            amrex::PhysBCFunctNoOp physbc;
            //amrex::GpuBndryFuncFab<MyExtBCFill> gpu_bndry_func(MyExtBCFill{});
            //amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill>> physbc(geom[lev],bcs,gpu_bndry_func);
            amrex::FillPatchSingleLevel(mf, time, fmf, ft, 0, icomp, ncomp,
                                            geom[lev], physbc, 0);
        }
        else
        {
            //amrex::GpuBndryFuncFab<MyExtBCFill> gpu_bndry_func(MyExtBCFill{});
            //amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill>> cphysbc(geom[lev-1],bcs,gpu_bndry_func);
            //amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill>> fphysbc(geom[lev],bcs,gpu_bndry_func);
            amrex::PhysBCFunctNoOp cphysbc, fphysbc;
            amrex::FillPatchTwoLevels(mf, time, cmf, ct, fmf, ft,
                                          0, icomp, ncomp, geom[lev-1], geom[lev],
                                          cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                          mapper, bcs, 0);
        }

    }

    if (lev == 0)
    {
        //amrex::PhysBCFunctNoOp physbc;
        //amrex::GpuBndryFuncFab<MyExtBCFill> gpu_bndry_func(MyExtBCFill{});
        //amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill>> physbc(geom[lev-1],bcs,gpu_bndry_func);
        //amrex::FillPatchSingleLevel(mf, time, fmf, ft, 0, icomp, ncomp,
        //                               geom[lev], physbc, 0);
        /*
        if(amrex::Gpu::inLaunchRegion())
        {
            amrex::GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<AmrCoreFill> > physbc(geom[lev],bcs,gpu_bndry_func);
            amrex::FillPatchSingleLevel(mf, time, fmf, ft, 0, icomp, ncomp,
                                        geom[lev], physbc, 0);
        }
        else
        {
            amrex::CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
            amrex::PhysBCFunct<amrex::CpuBndryFuncFab> physbc(geom[lev],bcs,bndry_func);
            amrex::FillPatchSingleLevel(mf, time, fmf, ft, 0, icomp, ncomp,
                                        geom[lev], physbc, 0);
        }
        */
    }
    else
    {
        /*amrex::PhysBCFunctNoOp cphysbc, fphysbc;
        amrex::FillPatchTwoLevels(mf, time, cmf, ct, fmf, ft,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                      mapper, bcs, 0);
        */
        /*amrex::GpuBndryFuncFab<MyExtBCFill> gpu_bndry_func(MyExtBCFill{});
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill>> cphysbc(geom[lev-1],bcs,gpu_bndry_func);
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill>> fphysbc(geom[lev],bcs,gpu_bndry_func);
        amrex::FillPatchTwoLevels(mf, time, cmf, ct, fmf, ft,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                      mapper, bcs, 0);
        */
        /*
        if(amrex::Gpu::inLaunchRegion())
        {
            amrex::GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[lev],bcs,gpu_bndry_func);

            amrex::FillPatchTwoLevels(mf, time, cmf, ct, fmf, ft,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                      mapper, bcs, 0);
        }
        else
        {
            amrex::CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
            amrex::PhysBCFunct<amrex::CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
            amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);

            amrex::FillPatchTwoLevels(mf, time, cmf, ct, fmf, ft,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                      mapper, bcs, 0);
        }
        */
    }
}

void incFSI::FillPatchScalars
(
    int lev, 
    int iscalar,
    amrex::Real time, 
    amrex::MultiFab& mf,                        // dest mf
    const amrex::Vector<amrex::MultiFab*>& cmf, // coarse mf
    const amrex::Vector<amrex::Real>& ct,       // coarse time
    const amrex::Vector<amrex::MultiFab*>& fmf, // fine mf
    const amrex::Vector<amrex::Real>& ft,       // fine time
    int icomp, int ncomp
)
{
    amrex::Interpolater *mapper;
    //amrex::Vector<amrex::BCRec> bcs(1);
    amrex::Array<amrex::Vector<amrex::BCRec>,AMREX_SPACEDIM> bcs;
    for (int n = 0; n < AMREX_SPACEDIM; ++n)
    {
        bcs[n].resize(1);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            bcs[n][0].setLo(idim, amrex::BCType::int_dir);
            bcs[n][0].setHi(idim, amrex::BCType::int_dir);
        }
    }

    if (mf.is_cell_centered())
    {
	//amrex::Print()<<"Lev = "<<lev<<" , iscalar = "<<iscalar<<'\n';
        mapper = &amrex::cell_cons_interp;

        //for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        //{
        //    bcs[0].setLo(idim, amrex::BCType::ext_dir);
        //    bcs[0].setHi(idim, amrex::BCType::ext_dir);
            //bcs[0].setLo(idim, amrex::BCType::reflect_even);
            //bcs[0].setHi(idim, amrex::BCType::reflect_even);
        //}


        if (lev == 0)
        {
            amrex::PhysBCFunctNoOp physbc;
            //amrex::GpuBndryFuncFab<MyExtBCFill> gpu_bndry_func;
            //amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill>> physbc(geom[lev],bcs,gpu_bndry_func);
            amrex::FillPatchSingleLevel(mf, time, fmf, ft, 0, icomp, ncomp,
                                            geom[lev], physbc, 0);
        }
        else
        {
            //amrex::GpuBndryFuncFab<MyExtBCFill> gpu_bndry_func;
            //amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill>> cphysbc(geom[lev-1],bcs,gpu_bndry_func);
            //amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill>> fphysbc(geom[lev],bcs,gpu_bndry_func);
            //amrex::PhysBCFunctNoOp cphysbc, fphysbc;
            amrex::GpuBndryFuncFab<MyExtBCFillScalars> gpu_bndry_func(
                                                    MyExtBCFillScalars{
                                                        iscalar,grids[lev],
                                                        lobc_Scalars[iscalar],hibc_Scalars[iscalar],
                                                        Scalars_bcf[iscalar],
                                                        isAxisymmetric});//default GpuBndryFuncFab
	    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFillScalars>> cphysbc(geom[lev-1],bcs[0],gpu_bndry_func);
            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFillScalars>> fphysbc(geom[lev],bcs[0],gpu_bndry_func);

            amrex::FillPatchTwoLevels(mf, time, cmf, ct, fmf, ft,
                                          0, icomp, ncomp, geom[lev-1], geom[lev],
                                          cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                          mapper, bcs[0], 0);            
        }

    }
    else if 
    (
        mf.ixType() == amrex::IndexType(amrex::IntVect::TheDimensionVector(0)) || 
        mf.ixType() == amrex::IndexType(amrex::IntVect::TheDimensionVector(1))
    )
    {
        mapper = &amrex::face_linear_interp;
        //for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        //{
        //    bcs[0].setLo(idim, amrex::BCType::int_dir);
        //    bcs[0].setHi(idim, amrex::BCType::int_dir);
            //bcs[0].setLo(idim, amrex::BCType::foextrap);
            //bcs[0].setHi(idim, amrex::BCType::foextrap);
            //bcs[0].setLo(idim, amrex::BCType::ext_dir);
            //bcs[0].setHi(idim, amrex::BCType::ext_dir);
        //}
        if (lev == 0)
        {
            amrex::PhysBCFunctNoOp physbc;
            //amrex::GpuBndryFuncFab<MyExtBCFill> gpu_bndry_func(MyExtBCFill{});
            //amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill>> physbc(geom[lev],bcs,gpu_bndry_func);
            amrex::FillPatchSingleLevel(mf, time, fmf, ft, 0, icomp, ncomp,
                                            geom[lev], physbc, 0);
        }
        else
        {
            //amrex::GpuBndryFuncFab<MyExtBCFillScalars> gpu_bndry_func(
            //                                        MyExtBCFillScalars{
            //                                            iscalar,grids[lev],
            //                                            lobc_Scalars[iscalar],hibc_Scalars[iscalar],
            //                                            Scalars_bcf[iscalar],
            //                                            isAxisymmetric});//default GpuBndryFuncFab

            //amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFillScalars>> cphysbc(geom[lev-1],bcs,gpu_bndry_func);
            //amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFillScalars>> fphysbc(geom[lev],bcs,gpu_bndry_func);
            amrex::PhysBCFunctNoOp cphysbc, fphysbc;
            amrex::FillPatchTwoLevels(mf, time, cmf, ct, fmf, ft,
                                          0, icomp, ncomp, geom[lev-1], geom[lev],
                                          cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                          mapper, bcs[0], 0);
        }
    }
}

//! fill face centered velocity
void incFSI::FillPatch
(
    int lev, 
    amrex::Real time, 
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> const& mf,                  // destination mf, contains u and v for 2D
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>& cmf,  // coarse mf, contains u and v for 2D
    const amrex::Vector<amrex::Real>& ct,                                      // coarse time
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>& fmf,  // fine mf, contains u and v for 2D
    const amrex::Vector<amrex::Real>& ft,                                      // fine time
    int icomp, int ncomp
)
{
    amrex::BoxArray xba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(0));
    amrex::BoxArray yba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(1));

    //! No need to fill external boundary values
    amrex::Array<amrex::Vector<amrex::BCRec>,AMREX_SPACEDIM> bcs;
    for (int n = 0; n < AMREX_SPACEDIM; ++n)
    {
        bcs[n].resize(1);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            bcs[n][0].setLo(idim, amrex::BCType::int_dir);
            bcs[n][0].setHi(idim, amrex::BCType::int_dir);
        }
    }
    /*{
        bcs[0].resize(1);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            bcs[0][0].setLo(idim, amrex::BCType::hoextrap);
            bcs[0][0].setHi(idim, amrex::BCType::hoextrap);
        }

        bcs[1].resize(1);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            bcs[1][0].setLo(idim, amrex::BCType::foextrap);
            bcs[1][0].setHi(idim, amrex::BCType::foextrap);
        }
    }*/
    /*
    amrex::Array<amrex::Vector<amrex::BCRec>,AMREX_SPACEDIM> bcs;
    for (int n = 0; n < AMREX_SPACEDIM; ++n)
    {   
        bcs[n].resize(1);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {   
            if (idim == n) 
            { // we could drop the test and set everything to either hoextrap or foextrap.
                bcs[n][0].setLo(idim, amrex::BCType::hoextrap);
                bcs[n][0].setHi(idim, amrex::BCType::hoextrap);
                //bcs[n][0].setLo(idim, amrex::BCType::foextrap);
                //bcs[n][0].setHi(idim, amrex::BCType::foextrap);
            }
            else 
            {
                bcs[n][0].setLo(idim, amrex::BCType::foextrap);
                bcs[n][0].setHi(idim, amrex::BCType::foextrap);
            }
        }
    }
    */

    if (lev == 0)
    {
        amrex::PhysBCFunctNoOp physbc;
        //amrex::GpuBndryFuncFab<MyExtBCFill> gpu_bndry_func;//default GpuBndryFuncFab
        //amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill> > physbc(geom[lev],bcs[0],gpu_bndry_func);
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            amrex::Vector<amrex::MultiFab *> fmf_tmp(fmf.size());
            for (size_t i = 0; i < fmf_tmp.size(); i++)
            {
                fmf_tmp[i] = fmf[i][dir];
            }
            amrex::FillPatchSingleLevel(*mf[dir], time, fmf_tmp, ft, 0, icomp, ncomp,
                                        geom[lev], physbc, 0);
        }
    }
    else
    {
        /*amrex::Interpolater* mapper = &amrex::face_divfree_interp;
        amrex::Array<amrex::PhysBCFunctNoOp, AMREX_SPACEDIM> cphysbc_arr;
        amrex::Array<amrex::PhysBCFunctNoOp, AMREX_SPACEDIM> fphysbc_arr;

        amrex::FillPatchTwoLevels(mf, time, cmf, ct, fmf, ft,
                            0, icomp, ncomp, geom[lev-1], geom[lev],
                            cphysbc_arr, 0, fphysbc_arr, 0, refRatio(lev-1),
                            mapper, bcs, 0);
        */
        //Implementation with GpuBndryFuncFab
        amrex::Interpolater* mapper = &amrex::face_divfree_interp;
        //amrex::GpuBndryFuncFab<MyExtBCFill> gpu_bndry_func;//default GpuBndryFuncFab
        //amrex::GpuBndryFuncFab<MyExtBCFill> gpu_bndry_func;//default GpuBndryFuncFab
        amrex::GpuBndryFuncFab<MyExtBCFill> gpu_bndry_func_u(
                                            MyExtBCFill{
                                                        1,xba,
                                                        lobc_u,hibc_u,
                                                        u_bcf,
                                                        isAxisymmetric});//default GpuBndryFuncFab
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill> > cphysbc_u(geom[lev-1],bcs[0],gpu_bndry_func_u);
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill> > fphysbc_u(geom[lev],bcs[0],gpu_bndry_func_u);

        amrex::GpuBndryFuncFab<MyExtBCFill> gpu_bndry_func_v(
                                            MyExtBCFill{2,yba,
                                                        lobc_v,hibc_v,
                                                        v_bcf,
                                                        isAxisymmetric});//default GpuBndryFuncFab
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill> > cphysbc_v(geom[lev-1],bcs[0],gpu_bndry_func_v);
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill> > fphysbc_v(geom[lev],bcs[0],gpu_bndry_func_v);

        amrex::Array<amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill> >, AMREX_SPACEDIM> cphysbc_arr = {cphysbc_u, cphysbc_v};
        amrex::Array<amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill> >, AMREX_SPACEDIM> fphysbc_arr = {fphysbc_u, fphysbc_v};
        
        amrex::FillPatchTwoLevels(mf, time, cmf, ct, fmf, ft,
                            0, icomp, ncomp, geom[lev-1], geom[lev],
                            cphysbc_arr, 0, fphysbc_arr, 0, refRatio(lev-1),
                            mapper, bcs, 0);
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void incFSI::FillCoarsePatch 
(
    int lev, 
    amrex::Real time, 
    amrex::MultiFab& mf, 
    const amrex::MultiFab& cmf, 
    int icomp, 
    int ncomp
)
{
    BL_ASSERT(lev > 0);

    amrex::Interpolater *mapper;
    if (mf.is_cell_centered())
    {
        mapper = &amrex::cell_cons_interp;
    }
    else if 
    (
        mf.ixType() == amrex::IndexType(amrex::IntVect::TheDimensionVector(0)) || 
        mf.ixType() == amrex::IndexType(amrex::IntVect::TheDimensionVector(1))
    )
    {
        mapper = &amrex::face_linear_interp;
    }

    //! No need to fill external boundary values
    amrex::Vector<amrex::BCRec> bcs(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        bcs[0].setLo(idim, amrex::BCType::int_dir);
        bcs[0].setHi(idim, amrex::BCType::int_dir);
    }

    //amrex::PhysBCFunctNoOp cphysbc, fphysbc;
    amrex::GpuBndryFuncFab<MyExtBCFill> gpu_bndry_func;
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill>> cphysbc(geom[lev-1],bcs,gpu_bndry_func);
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<MyExtBCFill>> fphysbc(geom[lev],bcs,gpu_bndry_func);

    amrex::InterpFromCoarseLevel(mf, time, cmf, 0, icomp, ncomp, geom[lev-1], geom[lev],
                                cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                mapper, bcs, 0);

    /*
    if(amrex::Gpu::inLaunchRegion())
    {
        amrex::GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[lev],bcs,gpu_bndry_func);

        amrex::InterpFromCoarseLevel(mf, time, cmf, 0, icomp, ncomp, geom[lev-1], geom[lev],
                                     cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                     mapper, bcs, 0);
    }
    else
    {
        amrex::CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
        amrex::PhysBCFunct<amrex::CpuBndryFuncFab> cphysbc(geom[lev-1],{bcs},bndry_func);
        amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc(geom[lev],{bcs},bndry_func);

        amrex::InterpFromCoarseLevel(mf, time, cmf, 0, icomp, ncomp, geom[lev-1], geom[lev],
                                     cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                     mapper, bcs, 0);
    }
    */
}

void incFSI::FillCoarsePatch 
(
    int lev, 
    amrex::Real time,
    amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> const &mf,
    const amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> &cmf,
    int icomp, 
    int ncomp
)
{
    BL_ASSERT(lev > 0);

    amrex::Interpolater *mapper = &amrex::face_divfree_interp;

    //! No need to fill external boundary values
    amrex::Array<amrex::Vector<amrex::BCRec>,AMREX_SPACEDIM> bcs;
    for (int n = 0; n < AMREX_SPACEDIM; ++n)
    {
        bcs[n].resize(1);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            bcs[n][0].setLo(idim, amrex::BCType::int_dir);
            bcs[n][0].setHi(idim, amrex::BCType::int_dir);
        }
    }

    amrex::Array<amrex::PhysBCFunctNoOp, AMREX_SPACEDIM> cphysbc_arr;
    amrex::Array<amrex::PhysBCFunctNoOp, AMREX_SPACEDIM> fphysbc_arr;


    amrex::InterpFromCoarseLevel(mf, time, cmf, 0, icomp, ncomp, geom[lev-1], geom[lev],
                                cphysbc_arr, 0, fphysbc_arr, 0, refRatio(lev-1),
                                mapper, bcs, 0);

    // for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        amrex::MultiFab *mfptr = mf[0];
        for (amrex::MFIter mfi(*mfptr); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            auto& f = (*mfptr)[mfi];
            amrex::PrintToFile("fine_fab") << "bx : " << bx << "\n";
            amrex::PrintToFile("fine_fab") << f << "\n";
        }
    }

    {
        amrex::MultiFab *mfptr = cmf[0];
        for (amrex::MFIter mfi(*mfptr); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            auto& f = (*mfptr)[mfi];
            amrex::PrintToFile("crse_fab") << "bx : " << bx << "\n";
            amrex::PrintToFile("crse_fab") << f << "\n";
        }
    }

    /*
    if(amrex::Gpu::inLaunchRegion())
    {
        amrex::GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[lev-1],bcs[0],gpu_bndry_func);
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[lev],bcs[0],gpu_bndry_func);

        amrex::Array<amrex::PhysBCFunct<amrex::GpuBndryFuncFab<AmrCoreFill> >, AMREX_SPACEDIM> cphysbc_arr = {cphysbc, cphysbc};
        amrex::Array<amrex::PhysBCFunct<amrex::GpuBndryFuncFab<AmrCoreFill> >, AMREX_SPACEDIM> fphysbc_arr = {fphysbc, fphysbc};

        amrex::InterpFromCoarseLevel(mf, time, cmf, 0, icomp, ncomp, geom[lev-1], geom[lev],
                                     cphysbc_arr, 0, fphysbc_arr, 0, refRatio(lev-1),
                                     mapper, bcs, 0);
    }
    else
    {
        amrex::CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
        amrex::PhysBCFunct<amrex::CpuBndryFuncFab> cphysbc(geom[lev-1],bcs[0],bndry_func);
        amrex::PhysBCFunct<amrex::CpuBndryFuncFab> fphysbc(geom[lev],bcs[0],bndry_func);

        amrex::Array<amrex::PhysBCFunct<amrex::CpuBndryFuncFab>, AMREX_SPACEDIM> cphysbc_arr = {cphysbc, cphysbc};
        amrex::Array<amrex::PhysBCFunct<amrex::CpuBndryFuncFab>, AMREX_SPACEDIM> fphysbc_arr = {fphysbc, fphysbc};

        // amrex::Array<amrex::PhysBCFunctNoOp, AMREX_SPACEDIM> cphysbc_arr;
        // amrex::Array<amrex::PhysBCFunctNoOp, AMREX_SPACEDIM> fphysbc_arr;

        amrex::InterpFromCoarseLevel(mf, time, cmf, 0, icomp, ncomp, geom[lev-1], geom[lev],
                                     cphysbc_arr, 0, fphysbc_arr, 0, refRatio(lev-1),
                                     mapper, bcs, 0);
    }
    */

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        if (cmf[i]->contains_inf())
        {
            amrex::PrintToFile("log") << "cvel " << i << " contains inf\n";
        }
        if (cmf[i]->contains_nan())
        {
            amrex::PrintToFile("log") << "cvel " << i << " contains nan\n";
        }
        if (mf[i]->contains_inf())
        {
            amrex::PrintToFile("log") << "vel " << i << " contains inf\n";
        }
        if (mf[i]->contains_nan())
        {
            amrex::PrintToFile("log") << "vel " << i << " contains nan\n";
        }
    }
}

void incFSI::FillPatchNearBoundary(int lev)
{

    //! create temp field data
    amrex::BoxArray xba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(0));
    amrex::BoxArray yba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(1));
    //amrex::MultiFab tmp_xvel(xba, dmap[lev], 1, Nghost);
    //amrex::MultiFab tmp_yvel(yba, dmap[lev], 1, Nghost);
    //amrex::MultiFab tmp_Pressure(grids[lev], dmap[lev], 1, Nghost);
    
    //tmp_xvel.setVal(0.0);
    //tmp_yvel.setVal(0.0);
    //tmp_Pressure.setVal(0.0);
    
    int myproc = amrex::ParallelDescriptor::MyProc();  // Return the rank
    int nprocs = amrex::ParallelDescriptor::NProcs();  // Return the number of processes
    
    //amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> tmp_vel = {&tmp_xvel, &tmp_yvel};
    //amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> cvel;
    //amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> fvel;
    
    //cvel.resize(1);
    //fvel.resize(1);
    
    //xvel[lev - 1].setVal(0.0);
    //xvel[lev].setVal(0.0);
    //yvel[lev - 1].setVal(0.0);
    //yvel[lev].setVal(0.0);
    int proc = 5;
    
    const amrex::Box &domain_c = geom[lev].Domain();
    const amrex::Real *prob_lo_c = geom[lev].ProbLo();
    const amrex::Real *prob_hi_c = geom[lev].ProbHi();
    const amrex::Real *dx_c = geom[lev].CellSize();

    //amrex::Box grown_domain(domain);
    //grown_domain.grow(1);
    
    amrex::Print(proc)<<"lv = "<<lev<<", myproc = "<< myproc <<" , nprocs = "<<nprocs<<"****************************"<<'\n';
    //amrex::Print(proc)<<"grown domain = "<<grown_domain<<'\n'; 
    amrex::Print(proc)<<"Grids at lev "<<lev<<" : "<<grids[lev]<<'\n';
    amrex::Print(proc)<<"domain at lev "<<lev<<" : "<<geom[lev].Domain()<<'\n';
    amrex::Print(proc)<<"dmap at lev "<<lev<<" : "<<dmap[lev]<<'\n';
    amrex::Print(proc)<<"probLo at lev"<<lev<<" : "<<prob_lo_c[0]<<" , "<<prob_lo_c[1]<<'\n';
    amrex::Print(proc)<<"probHi at lev"<<lev<<" : "<<prob_hi_c[0]<<" , "<<prob_hi_c[1]<<'\n';
    amrex::Print(proc)<<"dx at lev"<<lev<<" : "<<dx_c[0]<<" , "<<dx_c[1]<<'\n';

    const amrex::Box &domain_f = geom[lev + 1].Domain();
    const amrex::Real *prob_lo_f = geom[lev + 1].ProbLo();
    const amrex::Real *prob_hi_f = geom[lev + 1].ProbHi();
    const amrex::Real *dx_f = geom[lev + 1].CellSize();

    amrex::Print(proc)<<"Grids at lev "<<lev+1<<" : "<<grids[lev + 1]<<'\n';
    amrex::Print(proc)<<"domain at lev "<<lev+1<<" : "<<geom[lev+1].Domain()<<'\n';
    amrex::Print(proc)<<"dmap at lev "<<lev+1<<" : "<<dmap[lev+1]<<'\n';
    amrex::Print(proc)<<"probLo at lev"<<lev + 1<<" : "<<prob_lo_f[0]<<" , "<<prob_lo_f[1]<<'\n';
    amrex::Print(proc)<<"probHi at lev"<<lev + 1<<" : "<<prob_hi_f[0]<<" , "<<prob_hi_f[1]<<'\n';
    amrex::Print(proc)<<"dx at lev"<<lev + 1<<" : "<<dx_f[0]<<" , "<<dx_f[1]<<'\n';
    //amrex::Print()<<"xba = "<<xba<<'\n';
    //amrex::Print()<<"yba = "<<yba<<'\n';
    for (amrex::MFIter mfi(xvel[lev]); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        amrex::Box fbx = get_valid_face_box(lev, bx, 0);
        amrex::Print(-1)<<"myproc = "<<myproc<<", lev = "<<lev<<" , xvel box = "<<bx<<'\n';
        //amrex::Print(proc)<<"face box = "<<fbx<<'\n'; 
    
        amrex::Box gbx_xvel(bx);
        //gbx_xvel.grow(1);
     
        //amrex::Print(proc)<<"grown box = "<<gbx_xvel<<'\n'; 
    
    
        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            //amrex::Print(proc)<<"i = "<<i<<" , j = "<<j<<" , k = "<<k<<"\n";
            //amrex::Print(proc)<<
        });

        for (amrex::BoxIterator bxi(bx); bxi.ok(); ++bxi)
        {
            //amrex::Print(-1)<<" bxi = "<<bxi.m_boxHi<<'\n';
        }
    }

    amrex::BoxList bxlst_crs = grids[lev].boxList(); 
    amrex::Vector<amrex::Box> bxvec_crs = bxlst_crs.data();

    for(amrex::Vector<amrex::Box>::iterator bx_ptr_crs = bxvec_crs.begin(); bx_ptr_crs < bxvec_crs.end(); bx_ptr_crs++)
    {
        amrex::Print(proc)<<"bxvec : coarse = "<<*bx_ptr_crs<<'\n';
    }

    amrex::BoxList bxlst_fine = grids[lev + 1].boxList();
    amrex::Vector<amrex::Box> bxvec_fine = bxlst_fine.data();

    for(amrex::Vector<amrex::Box>::iterator bx_ptr_fine = bxvec_fine.begin(); bx_ptr_fine < bxvec_fine.end(); bx_ptr_fine++)
    {
        amrex::Print(proc)<<"bxvec : fine = "<<*bx_ptr_fine<<'\n';
    }
    amrex::Vector<int> m_pmap_c = dmap[lev].ProcessorMap();
    amrex::Print(proc)<<"dmap[lev] = "<<m_pmap_c[proc]<<'\n';
    amrex::Vector<int> m_pmap_f = dmap[lev+1].ProcessorMap();
    amrex::Print(proc)<<"dmap[lev+1] = "<<m_pmap_f[proc]<<'\n';

    amrex::Vector<int> test_vec(nprocs);
   
    for(int i = 0; i < test_vec.size(); i++)
        test_vec[i] = 0;

    test_vec[myproc] = myproc;

    amrex::Print(-1)<<"test_vec[myproc] = "<<test_vec[myproc]<<'\n';
   
    amrex::Vector<int> r_test_vec(nprocs);
    if(myproc == 1)
    {
        amrex::ParallelDescriptor::Send( &test_vec[0], nprocs, 5, 0);
    }
    else if(myproc == 5)
    {
        amrex::ParallelDescriptor::Recv( r_test_vec, 1, 0);
    }
 
    for(int i = 0; i < r_test_vec.size(); i++)
        amrex::Print(5)<<"r_test_vec["<<i<<"] = "<<r_test_vec[i]<<'\n';


    //XVelBoundaryConditions(lev, Time, xvel[lev - 1]);
    //YVelBoundaryConditions(lev, Time, yvel[lev - 1]);
    //cvel[0] = {&xvel[lev - 1], &yvel[lev - 1]};
    //XVelBoundaryConditions(lev, Time, xvel[lev]);
    //YVelBoundaryConditions(lev, Time, yvel[lev]);
    //fvel[0] = {&xvel[lev], &yvel[lev]};
    //{
        //if(divergence_free_interpolation)
        //{
        //    FillPatch(lev, Time, tmp_vel, cvel, {Time}, fvel, {Time}, 0, 1);
        //    std::swap(tmp_vel[0][0], xvel[lev]);
        //    std::swap(tmp_vel[0][1], yvel[lev]);
            //XVelBoundaryConditions(lev, Time, xvel[lev]);
            //YVelBoundaryConditions(lev, Time, yvel[lev]);
        //}
        //else
        //{
        //    FillPatch(lev, Time, tmp_xvel, {&xvel[lev - 1]}, {Time}, {&xvel[lev]}, {Time}, 0, 1);
        //    XVelBoundaryConditions(lev, Time, tmp_xvel);
        //    std::swap(tmp_xvel, xvel[lev]);
    
        //    FillPatch(lev, Time, tmp_yvel, {&yvel[lev - 1]}, {Time}, {&yvel[lev]}, {Time}, 0, 1);
        //    YVelBoundaryConditions(lev, Time, tmp_yvel);
        //    std::swap(tmp_yvel, yvel[lev]);
        //}
    
        //FillPatch(lev, Time, tmp_Pressure, {&Pressure[lev - 1]}, {Time}, {&Pressure[lev]}, {Time}, 0, 1);
        //PressureBoundaryConditions(lev, Time, tmp_Pressure);
        //std::swap(tmp_Pressure, Pressure[lev]);
    
    //}    
}

} // namespace mycode
