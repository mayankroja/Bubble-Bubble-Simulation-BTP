#include <incFSI.H>
#include <MLPoisson.H>
#include <AMReX_MultiFabUtil.H>

namespace mycode
{

void incFSI::AMGPressurePoisson_RK(int RKStage)
{
    //! compute RHS on all levels
    for (int lev = 0; lev <= finest_level; lev++)
    {
        const amrex::Real *dxinv = geom[lev].InvCellSize();
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

        for(amrex::MFIter mfi(RHS[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<amrex::Real> const &rhs = RHS[lev].array(mfi);
            amrex::Array4<amrex::Real const> const &u = xvel[lev].const_array(mfi);
            amrex::Array4<amrex::Real const> const &v = yvel[lev].const_array(mfi);
            amrex::Array4<int const> pmask;
            if(N_IF == 0)
                pmask = PMask_NI.const_array(mfi);
            else
                pmask = PMask->const_array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                if (pmask(i, j ,k) == 1)
                {   
                    if (isAxisymmetric)
                    {   
                        amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                        amrex::Real yn = prob_lo[1] + dx[1] * (j + 1);
                        amrex::Real ys = prob_lo[1] + dx[1] * j;

                        rhs(i, j, k) = (u(i + 1, j, k) - u(i, j, k)) / dx[0]
                                 + (yn * v(i, j + 1, k) - ys * v(i, j, k)) / (y * dx[1]);
                    }
                    else
                    {
                        rhs(i, j, k) = dxinv[0] * (u(i + 1, j, k) - u(i, j, k))
                                     + dxinv[1] * (v(i, j + 1, k) - v(i, j, k));
                    }

                    if (!Projection)
                    {
                        rhs(i, j, k) /= (rk_c[RKStage - 1] * dt);
                    }
                }
            });
        }
    }

    //! use MLPoisson
    MPI_Comm global_comm = amrex::ParallelContext::CommunicatorAll();

    MLPoisson PEqn(global_comm, &geom, &grids, &dmap, &interfaces);
    //MLPoisson PEqn(global_comm, &geom, &grids, &dmap);

    PEqn.setBCs(lobc_p, hibc_p);

    for (size_t i = 0; i < 2 * AMREX_SPACEDIM; i++)
    {
        PEqn.setBoundaryFunc(i, p_bcf[i]);
    }

    amrex::Vector<amrex::iMultiFab> PMask_NI;
    amrex::Vector<amrex::iMultiFab *> PMask;//I don't know how to work with this
    PMask_NI.resize(finest_level + 1);
    PMask.resize(finest_level + 1); 

    for (int lev = 0; lev <= finest_level; lev++)
    {
        if(N_IF == 0)
        {
            PMask_NI[lev].define(Pressure[lev].boxArray(), Pressure[lev].DistributionMap(), 1, Nghost);
            PMask_NI[lev].setVal(1);
        }
        else
        {
            PMask[lev] = &mask_[lev]->getPMask();
        }
    }


    if(N_IF == 0)
        PEqn.setMask(amrex::GetVecOfPtrs(PMask_NI));
    else
        PEqn.setMask(PMask);
    PEqn.setRhs(amrex::GetVecOfPtrs(RHS));
    PEqn.setSoln(amrex::GetVecOfPtrs(Pressure));
    
    PEqn.solve(0.0);
    
    poisson_iter = PEqn.getNIterations();
    system_size = PEqn.getSystemSize();
    poisson_tol = PEqn.getResidue();
    
    //Fill patch around the boxes, probably not required
    /*for (int lev = 0;lev <= finest_level;lev++)
    {
        if(lev > 0)
        {
            amrex::MultiFab tmp_Pressure(grids[lev], dmap[lev], 1, Nghost);
            FillPatch(lev, Time, tmp_Pressure,{&Pressure[lev - 1]}, {0.0}, {&Pressure[lev]}, {0.0}, 0, 1);
            std::swap(tmp_Pressure, Pressure[lev]);
        }
    }*/
    //! fill the ghost cells on fine levels with quadratic interpolation
    if(N_IF > 0)
    {
        //for (int lev = 0; lev <= finest_level; lev++)
        {
            int lev = finest_level;
            bool IsBubble = false;
            std::vector<amrex::Array4<amrex::Real const>> psi_fab(interfaces[lev].size()) ;
            amrex::iMultiFab &PMask = mask_[lev]->getPMask();
            for (amrex::MFIter mfi(PMask); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();
                amrex::Array4<amrex::Real> const &P = Pressure[lev].array(mfi);
                amrex::Array4<int const> const &mask = PMask.const_array(mfi);
            
                for (size_t i = 0; i < interfaces[lev].size(); i++)
                {
                    psi_fab[i] = interfaces[lev][i]->Psi().const_array(mfi);
                }
            
                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (mask(i,j,k) == 0)
                    {
                        IsBubble = false;
                        for (size_t p = 0; p < interfaces[lev].size(); p++)
                        {
                            if (psi_fab[p](i, j, k) > 0.0)
                            {
                                if (IsBubble)
                                {
                                    amrex::Print() << " Encountered merger!\n";
                                    exit(1);
                                }
                                //P(i, j, k) = interfaces[lev][p]->getP_interface();
                                IsBubble = true;
                            }
                        }
                        if (!IsBubble)
                        {
                            amrex::Print() << " Didn't find any bubble!\n";
                            exit(1);
                        }
                    }
                });
            }
        }
    }

    //! Avg down Pressure
    for (int lev = finest_level; lev > 0; lev--)
    {
        amrex::average_down(Pressure[lev], Pressure[lev - 1], 0, 1, ref_ratio[lev-1][0]);
    }


    PEqn.fillBoundaryFine();

    //if(N_IF > 0)
    //    for (int lev = 0; lev <= finest_level; lev++) 
    //        mask_[lev]->FillInGhost(Pressure[lev], Pintbcs);//Fill the pressure and velocity of the cut cells and 
                           //fill the first layer of ghost cell with pressure
}

void incFSI::ComputeFinalVelocityField_RK(int lev, int RKStage)
{
    PressureBoundaryConditions(lev, Time, Pressure[lev]);

    const amrex::Real *dx = geom[lev].CellSize();
    const amrex::Real *prob_lo = geom[lev].ProbLo();
    amrex::BoxArray xba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(0));
    amrex::BoxArray yba = amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(1));


    amrex::iMultiFab *UMask;
    amrex::iMultiFab *VMask;

    amrex::iMultiFab UMask_NI;//NI stands for No Interface, need a better way of doing this
    amrex::iMultiFab VMask_NI;//NI stands for No Interface, need a better way of doing this

    if(N_IF > 0)
    {
        auto &&mask__= mask_[lev];
        UMask = &mask__->getUMask();
        VMask = &mask__->getVMask();
    }
    else
    {
        UMask_NI.define(xba, dmap[lev], 1, Nghost);
        VMask_NI.define(yba, dmap[lev], 1, Nghost);

        UMask_NI.setVal(1);
        VMask_NI.setVal(1);
    }


    //amrex::Print()<<"grids["<<lev<<"] = "<<grids[lev]<<'\n';
    /// correct x-component
    amrex::Real gradPx;
    for (amrex::MFIter mfi(xvel[lev]); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const &P = Pressure[lev].const_array(mfi);
        amrex::Array4<amrex::Real> const &u = xvel[lev].array(mfi);

        amrex::Box box = get_valid_face_box(lev, bx, 0);
        
        //amrex::Print()<<"x bx = "<<bx<<" , face box = "<<box<<"\n";

        amrex::Array4<int const> umask;
        if(N_IF == 0)
            umask = UMask_NI.const_array(mfi);
        else
            umask = UMask->const_array(mfi);
        

        amrex::ParallelFor(box,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
        {
            if(umask(i ,j ,k) == 1)
            {
                gradPx = (P(i, j, k) - P(i - 1, j, k)) / dx[0];
                if (Projection)
                {
                    u(i, j, k) -= gradPx;
                }
                else
                {
                    u(i, j, k) -= rk_c[RKStage - 1] * dt * gradPx;
                }
            }            
        });
    }

    /// compute y-component
    amrex::Real gradPy;
    for (amrex::MFIter mfi(yvel[lev]); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const &P = Pressure[lev].const_array(mfi);
        amrex::Array4<amrex::Real> const &v = yvel[lev].array(mfi);

        amrex::Box box = get_valid_face_box(lev, bx, 1);
  
        //amrex::Print()<<"y bx = "<<bx<<" , face box = "<<box<<"\n";

        amrex::Array4<int const> vmask;
        if(N_IF == 0)
            vmask = VMask_NI.const_array(mfi);
        else
            vmask = VMask->const_array(mfi);


        amrex::ParallelFor(box,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
        {
            if(vmask(i ,j ,k) == 1)
            {
                gradPy = (P(i, j, k) - P(i, j - 1, k)) / dx[1];
                if (Projection)
                {
                    v(i, j, k) -= gradPy;
                }
                else
                {
                    v(i, j, k) -= rk_c[RKStage - 1] * dt * gradPy;
                }            
            }
            /*if(i == 350 && (j == 511 || j == 512))
            {
                amrex::Print()<<"Correct y Vel"<<'\n';
                amrex::Print()<<"i = "<<i<<" , j = "<<j<<'\n';
                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                amrex::Real y = prob_lo[1] + dx[1] * j;
                amrex::Print()<<"x = "<<x<<" , y = "<<y<<'\n';
                amrex::Print()<<"v("<<i<<", "<<j<<", "<<k<<") = "<<v(i, j, k)<<'\n';
            }*/
        });
    }

    if (Projection)
    {
        for (amrex::MFIter mfi(Pressure[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real> const &p = Pressure[lev].array(mfi);
            amrex::Array4<amrex::Real const> const &pstar = Pstar[lev].const_array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                p(i, j, k) = pstar(i, j, k) + p(i, j, k) / dt;
            });
        }
        PressureBoundaryConditions(lev, Time, Pressure[lev]);
    }
}

} // namespace mycode
