#include <MLPoisson.H>
#include <HypreIJ.H>
//#include <PetscIJ.H>

#include <AMReX_IArrayBox.H>
#include <CFStencil.H>

namespace mycode
{

void MLPoisson::addGhostContro
(
    const amrex::Box &domain,
    const amrex::IntVect &gcell,            /// the cut cell
    Stencil &st,                            /// st to which wts need to be added
    amrex::Array4<int const> const &cellid, /// the traverse index data
    const amrex::Real &fact,                /// scale factor
    const InterceptData &icpt,              /// intercept data @ cut cell (to get computed wts)
    amrex::Array4<int const> const &mask    /// PMask
)
{
    amrex::Box bx(gcell, gcell);
    bx.grow(icpt.stencil_);
    amrex::Box bx_isect = bx & domain;
    //! set the rhs contro from interpolation wts corresponding to cut pts
    for (int i = 0; i < icpt.n_intercepts; i++)
    {
        st.addToSource(-icpt.Interpolation_weights[i] * icpt.P[i] * fact);
    }

    int i = icpt.n_intercepts;
    for (amrex::BoxIterator bit(bx_isect); bit.ok(); ++bit)
    {
        const amrex::IntVect &iv = bit();
        if (mask(iv) == 1)
        {
            st.addToColEntry(cellid(iv), icpt.Interpolation_weights[i] * fact);
            i++;
        }
    }
}

void MLPoisson::solve(amrex::Real Time)
{
    amrex::Real start_time = amrex::second();

    // check if singular
    // search if any dirichlet bcs present
    auto itplo = std::find(lobc.begin(), lobc.end(), amrex::LinOpBCType::Dirichlet);
    auto itphi = std::find(hibc.begin(), hibc.end(), amrex::LinOpBCType::Dirichlet);
    auto itp_Int = std::find(intbc.begin(), intbc.end(), amrex::LinOpBCType::Dirichlet);

    if (itplo == lobc.end() && itphi == hibc.end() && itp_Int == intbc.end())
        singular_ = true;

    MLTraverseIndex tridx(*geom_, *grid_, *dmap_, mask_, soln_[0]->nGrow()+1);
    global_sys_size = tridx.getWorldSize();

    int nlevels = geom_->size();
    amrex::Vector<CFMask> cfmask_(nlevels);
    for (int ilev = 0; ilev < nlevels; ilev++)
    {
        cfmask_[ilev].define(ilev, *geom_, *grid_, *dmap_);
        sys_size_lev[ilev] = tridx.getLevSize(ilev);
    }

    if (solver_lib == LinearSolver::linear_solver_lib::HYPRE)
    {
        solver = std::make_unique<HypreIJ>(comm_, geom_, grid_, dmap_, &cfmask_, mask_, &tridx);
    }
    else
    {
        //solver = std::make_unique<PetscIJ>(comm_, geom_, grid_, dmap_, &cfmask_, mask_, &tridx);
        //if (singular_)
        //{
        //    solver->singularSystem(true);
        //}
    }

    for (int ilev = 0; ilev < nlevels; ilev++)
    {
        /// setup system and solve
        const amrex::Real *dx      = (*geom_)[ilev].CellSize();
        const amrex::Real *dxinv   = (*geom_)[ilev].InvCellSize();
        const amrex::Box &domain   = (*geom_)[ilev].Domain();
        const amrex::Real *prob_lo = (*geom_)[ilev].ProbLo();
        const amrex::Real *prob_hi = (*geom_)[ilev].ProbHi();

        amrex::Real scale_fact = dx[0] * dx[0];

        //! set values of A, x, b
        {
            amrex::Real E_coeff = dxinv[0] * dxinv[0];
            amrex::Real W_coeff = dxinv[0] * dxinv[0];
            amrex::Real N_coeff = dxinv[1] * dxinv[1];
            amrex::Real S_coeff = dxinv[1] * dxinv[1];

            // scale both sides of equation by h^2
            E_coeff = 1.0;
            W_coeff = 1.0;
            N_coeff *= scale_fact;
            S_coeff *= scale_fact;
            amrex::Real P_coeff = -(E_coeff + W_coeff + N_coeff + S_coeff);

            amrex::Array4<int const> mask;
            amrex::IArrayBox mask_fab;

            //! get the index of unknowns on curr, crse and fine level
            std::unique_ptr<amrex::iMultiFab> tridx_crse;
            std::unique_ptr<amrex::iMultiFab> tridx_fine;
            const amrex::iMultiFab &tridx_curr = tridx.getTrIndex(ilev);

            if (ilev > 0)
            {
                tridx_crse = tridx.getTrIndex(ilev, ilev - 1);
            }
            if (ilev < nlevels-1)
            {
                tridx_fine = tridx.getTrIndex(ilev, ilev + 1);
            }

            amrex::Array4<int const> cell_id_crse;
            amrex::Array4<int const> cell_id_fine;

            for (amrex::MFIter mfi(*soln_[ilev]); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();

                amrex::Array4<amrex::Real> const &phi = soln_[ilev]->array(mfi);
                amrex::Array4<amrex::Real const> const &rhs = rhs_[ilev]->const_array(mfi);
                amrex::Array4<int const> const &cellid = tridx_curr.const_array(mfi);
                amrex::Array4<int const> const &cfmask = cfmask_[ilev].Mask().const_array(mfi);

                if (tridx_crse)
                {
                    cell_id_crse = tridx_crse->const_array(mfi);
                }
                if (tridx_fine)
                {
                    cell_id_fine = tridx_fine->const_array(mfi);
                }

                if (mask_[ilev])
                {
                    mask = mask_[ilev]->const_array(mfi);
                }
                else
                {
                    mask_fab.resize(bx);
                    mask_fab.setVal(1);
                    mask = mask_fab.const_array();
                }

                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (mask(i, j, k) == 1 && cfmask(i,j,k) != CFMask::covered)
                    {
                        //! Assuming y as radial dir and x as z for axisymmetric problems
                        //! modify the N and S coefficients
                        if (isAxisymmetric)
                        {
                            //Assuming uniform cartesian grid in amrex
                            amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
			    amrex::Real yn = prob_lo[1] + dx[1] * (j + 1.0);
                            S_coeff = dxinv[1] * dxinv[1] - 1.0 / (2.0 * y * dx[1]);
                            N_coeff = dxinv[1] * dxinv[1] + 1.0 / (2.0 * y * dx[1]);
                            P_coeff = -1.0*(2.0/(dx[0] * dx[0]) + 2.0/(dx[1] * dx[1]))*dx[0] * dx[0];
                            
                            S_coeff *= dx[0] * dx[0];
                            N_coeff *= dx[0] * dx[0];

			    if(j == 0)
			    {
			        //amrex::Print()<<"At axis: i = "<<i<<" , "<<j<<" , y = "<<y<<'\n';

                                S_coeff =  0.0;
                                N_coeff = dxinv[1] * dxinv[1] * (yn/y);
				P_coeff = -1.0*(2.0/(dx[0] * dx[0]) + (yn/y) * 1.0/(dx[1] * dx[1]))*dx[0] * dx[0];

                                S_coeff *= dx[0] * dx[0];
                                N_coeff *= dx[0] * dx[0];
			    }
			    if(y < 0.0)
		            {
			       amrex::Print(-1)<<" Poisson solver is not set up for full domain axisymmetric calcualtion"<<'\n';
			       amrex::Print(-1)<<" Check in Src/MLPoisson/MLPoisson_solve.cpp"<<'\n';
			       exit(1);
			    }
                            
                            /*if(fabs(y) < 0.9*dx[1])
                            {
                                if(y >= 0)
                                {
                                    S_coeff = dxinv[1] * dxinv[1];
                                    N_coeff = dxinv[1] * dxinv[1] + 1.0 / (y * dx[1]);
                            
                                    S_coeff *= dx[0] * dx[0];
                                    N_coeff *= dx[0] * dx[0];
                            
                                    P_coeff = (-2.0/(dx[0] * dx[0]) - 2.0/(dx[1] * dx[1]))*dx[0] * dx[0] - 1.0 / (  y * dx[1])* dx[0] * dx[0];
                                }
                                else if(y < 0)
                                {
                                    N_coeff = dxinv[1] * dxinv[1];
                                    S_coeff = dxinv[1] * dxinv[1] - 1.0 / (y * dx[1]);
                            
                                    S_coeff *= dx[0] * dx[0];
                                    N_coeff *= dx[0] * dx[0];
                            
                                    P_coeff = (-2.0/(dx[0] * dx[0]) - 2.0/(dx[1] * dx[1]))*dx[0] * dx[0] + 1.0 / (  y * dx[1])*dx[0] * dx[0];
                                }
                            }*/
                        }

                        amrex::IntVect diag_iv(AMREX_D_DECL(i, j, k));
                        amrex::IntVect bottom_iv(AMREX_D_DECL(i, j - 1, k));
                        amrex::IntVect left_iv(AMREX_D_DECL(i - 1, j, k));
                        amrex::IntVect right_iv(AMREX_D_DECL(i + 1, j, k));
                        amrex::IntVect top_iv(AMREX_D_DECL(i, j + 1, k));
                        
                        amrex::IntVect test_iv(AMREX_D_DECL(305, 511, 0));
                        //amrex::IntVect test_iv1(AMREX_D_DECL(234, 271, 0));

                        int P = cellid(i, j, k);
                        int W = cellid(i - 1, j, k);
                        int S = cellid(i, j - 1, k);
                        int E = cellid(i + 1, j, k);
                        int N = cellid(i, j + 1, k);

                        /// create Stencil object
                        /// this is done to simplify setup Matrix row entries when we have cut cell with too many wts
                        Stencil st;

                        /// P_coeff
                        st.push_back_colEntry(P, P_coeff);

                        /// scale rhs by h^2
                        st.setSource(rhs(i, j, k) * scale_fact);
                        st.setRowId(P);

                        //! W_coeff
                        {
                            /// W outside domain
                            if (i - 1 < domain.smallEnd(0))
                            {
                                amrex::Real x = prob_lo[0];
                                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                                amrex::Real P_b = boundary_funcs[0](x, y, Time);
                                if (lobc[0] == amrex::LinOpBCType::Neumann)
                                {
                                    /// P_w = P_p - dx * P_b
                                    st.addToColEntry(P, W_coeff);
                                    st.addToSource(P_b * dx[0] * W_coeff);
                                }
                                else if (lobc[0] == amrex::LinOpBCType::Dirichlet)
                                {
                                    /// P_w = 2*P_b - P_p
                                    st.addToColEntry(P, -W_coeff);
                                    st.addToSource(-2.0 * P_b * W_coeff);
                                }
                                else
                                {
                                    amrex::Abort("AMG MLPoisson : invalid bc left");
                                }
                            }

                            //! W is out of this level grid
                            //! as i-1 is not out of domain but left cell is outside grid
                            //! => coarse mesh @ left ... lev/lev-1 bndry
                            else if (!(*grid_)[ilev].contains(left_iv))
                            {
                                int dir = 0, cross_dir = 1, side = 0;
                                mycode::Stencil gst = getGhostStencil_fine(cross_dir, dir, side, left_iv, 
                                                                          (*grid_)[ilev], domain, cellid, cell_id_crse);

                                //! flux on west face = phi_p - phi_g
                                //! => 1/dx^2 * (phi_g)
                                st += gst * W_coeff;
                            }

                            //! W is covered by fine mesh
                            //! lev/lev+1 bndry
                            else if (cfmask(left_iv) == CFMask::covered )
                            {
                                int dir = 0, cross_dir = 1, side = 1;
                                const amrex::Real *dxf = (*geom_)[ilev + 1].CellSize();
                                amrex::IntVect Wcell_f_p(left_iv[0] * 2 + 1, left_iv[1] * 2 + 1);
                                amrex::IntVect Wcell_f_m(left_iv[0] * 2 + 1, left_iv[1] * 2);
                                mycode::Stencil phi_gp = getGhostStencil_coarse(cross_dir, dir, side, diag_iv,
                                                                                Wcell_f_p, cfmask,
                                                                                cellid, cell_id_fine);
                                mycode::Stencil phi_gm = getGhostStencil_coarse(cross_dir, dir, side, diag_iv,
                                                                                Wcell_f_m, cfmask,
                                                                                cellid, cell_id_fine);

                                //! flux on W face .. (phi_ghost - W_fine)/dx_fine
                                phi_gp.addToColEntry(cell_id_fine(Wcell_f_p), -1.0);
                                phi_gm.addToColEntry(cell_id_fine(Wcell_f_m), -1.0);

                                mycode::Stencil Wflux = (phi_gp + phi_gm) / (2.0 * dxf[0]);

                                //! as coeff corresponding to col P contains -2*W_coeff
                                //! remove W_coeff as W cell is not present
                                st.addToColEntry(P, W_coeff);

                                //! add flux contro of W
                                st += Wflux * (-scale_fact / dx[0]);
                            }

                            /// W is cutcell
                            else if (mask(i - 1, j, k) == 2)
                            {
                                /// find the interceptData corresponding to W
                                /// this could be made faster at the expense of memory by having
                                /// Indicator (which is mapping between cell index to index in vector of interceptData)

                                //! the intercept data is specific to each interface
                                //! iterate over interfaces and find the intercept data for corresponding cell
                                //! need a better way, but for now working

                                for (auto &&solid : (*interfaces)[ilev])
                                {
                                    const auto &icpt_data = solid->getInterceptData()[mfi];
                                    if (!icpt_data.empty())
                                    {
                                        auto itr = std::find_if(icpt_data.begin(), icpt_data.end(),
                                                                [&](const InterceptData &idt) {
                                                                    return idt.cellid_ == left_iv;
                                                                });

                                        if (itr != icpt_data.end())
                                        {
                                            addGhostContro(domain, left_iv, st, cellid, W_coeff, *itr, mask);
                                        }                                
                                    }
                                }
                            }

                            /// W is regular cell
                            else
                            {
                                st.addToColEntry(W, W_coeff);
                            }
                        }

                        //! S coeff
                        {
                            /// S out
                            if (j - 1 < domain.smallEnd(1))
                            {
                                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                                amrex::Real y = prob_lo[1];
                                amrex::Real P_b = boundary_funcs[1](x, y, Time);
                                if (lobc[1] == amrex::LinOpBCType::Neumann)
                                {
                                    /// P_s = P_p - P_b * dy
                                    st.addToColEntry(P, S_coeff);
                                    st.addToSource(P_b * dx[1] * S_coeff);
                                }
                                else if (lobc[1] == amrex::LinOpBCType::Dirichlet)
                                {
                                    /// P_s = 2.0*P_b - P_p
                                    st.addToColEntry(P, -S_coeff);
                                    st.addToSource(-2.0 * P_b * S_coeff);
                                }
                                else
                                {
                                    amrex::Abort("AMG MLPoisson : invalid bc bottom");
                                }
                            }

                            //! S is out of this level grid
                            //! as j-1 is not out of domain but bottom cell is outside grid
                            //! => coarse mesh @ bottom ... lev/lev-1 bndry
                            else if (!(*grid_)[ilev].contains(bottom_iv))
                            {
                                int dir = 1, cross_dir = 0, side = 0;
                                mycode::Stencil gst = getGhostStencil_fine(cross_dir, dir, side, bottom_iv, 
                                                                          (*grid_)[ilev], domain, cellid, cell_id_crse);

                                st += gst * S_coeff;
                            }

                            //! S is covered by fine mesh
                            //! lev/lev+1 bndry
                            else if (cfmask(bottom_iv) == CFMask::covered)
                            {
                                int dir = 1, cross_dir = 0, side = 1;
                                const amrex::Real *dxf = (*geom_)[ilev + 1].CellSize();
                                amrex::IntVect Scell_f_p(bottom_iv[0] * 2 + 1, bottom_iv[1] * 2 + 1);
                                amrex::IntVect Scell_f_m(bottom_iv[0] * 2, bottom_iv[1] * 2 + 1);
                                mycode::Stencil phi_gp = getGhostStencil_coarse(cross_dir, dir, side, diag_iv,
                                                                                Scell_f_p, cfmask,
                                                                                cellid, cell_id_fine);
                                mycode::Stencil phi_gm = getGhostStencil_coarse(cross_dir, dir, side, diag_iv,
                                                                                 Scell_f_m, cfmask,
                                                                                 cellid, cell_id_fine);

                                //! flux on S face .. (phi_ghost - S_fine)/dy_fine
                                phi_gp.addToColEntry(cell_id_fine(Scell_f_p), -1.0);
                                phi_gm.addToColEntry(cell_id_fine(Scell_f_m), -1.0);

                                mycode::Stencil Sflux = (phi_gp + phi_gm) / (2.0 * dxf[1]);

                                //! as coeff corresponding to col P contains -2*S_coeff
                                //! remove S_coeff as S cell is not present
                                st.addToColEntry(P, S_coeff);

                                //! add flux contro of S
                                st += Sflux * (-scale_fact / dx[1]);
                            }

                            /// S is cut cell
                            else if (mask(i, j - 1, k) == 2)
                            {
                                /// find the interceptData corresponding to S
                                for (auto &&solid : (*interfaces)[ilev])
                                {
                                    const auto &icpt_data = solid->getInterceptData()[mfi];

                                    if (!icpt_data.empty())
                                    {
                                        auto itr = std::find_if(icpt_data.begin(), icpt_data.end(),
                                                                [&](const InterceptData &idt) {
                                                                    return idt.cellid_ == bottom_iv;
                                                                });
                                        //if(diag_iv == test_iv)
                                        //{
                                        //    std::cout<<"test_iv = "<<test_iv[0]<<" , "<<test_iv[1]<<'\n';
                                        //    std::cout<<"bottom_iv = "<<bottom_iv[0]<<" , "<<bottom_iv[1]<<'\n';
                                            //std::cout<<"itr = "<<*itr<<'\n';
                                        //}
                                        if (itr != icpt_data.end())
                                        {
                                            addGhostContro(domain, bottom_iv, st, cellid, S_coeff, *itr, mask);
                                            //if(diag_iv == test_iv)
                                            //   std::cout<<"adding bottom_iv to ghost "<<'\n';
                                        }
                                    }
                                }
                            }

                            /// S is regular cell
                            else
                            {
                                st.addToColEntry(S, S_coeff);
                            }
                        }

                        //! E coeff
                        {
                            // E outside domain
                            if (i + 1 > domain.bigEnd(0))
                            {
                                amrex::Real x = prob_hi[0];
                                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                                amrex::Real P_b = boundary_funcs[2](x, y, Time);

                                if (hibc[0] == amrex::LinOpBCType::Neumann)
                                {
                                    /// P_E = P_b * dx + P_p
                                    st.addToColEntry(P, E_coeff);
                                    st.addToSource(-P_b * dx[0] * E_coeff);
                                }
                                else if (hibc[0] == amrex::LinOpBCType::Dirichlet)
                                {
                                    /// P_E = 2.0*P_b - P_p
                                    st.addToColEntry(P, -E_coeff);
                                    st.addToSource(-2.0 * P_b * E_coeff);
                                }
                                else
                                {
                                    amrex::Abort("AMG MLPoisson : invalid bc right");
                                }
                            }
                            
                            //! E is out of this level grid
                            //! as i+1 is not out of domain but right cell is outside grid
                            //! => coarse mesh @ right ... lev/lev-1 bndry
                            else if (!(*grid_)[ilev].contains(right_iv))
                            {
                                int dir = 0, cross_dir = 1, side = 1;
                                mycode::Stencil gst = getGhostStencil_fine(cross_dir, dir, side, right_iv, 
                                                                           (*grid_)[ilev], domain, cellid, cell_id_crse);

                                st += gst * E_coeff;
                            }

                            //! E is covered by fine mesh
                            //! lev/lev+1 bndry
                            else if (cfmask(right_iv) == CFMask::covered)
                            {
                                int dir = 0, cross_dir = 1, side = 0;
                                const amrex::Real *dxf = (*geom_)[ilev + 1].CellSize();
                                amrex::IntVect Ecell_f_p(right_iv[0] * 2, right_iv[1] * 2 + 1);
                                amrex::IntVect Ecell_f_m(right_iv[0] * 2, right_iv[1] * 2);
                                mycode::Stencil phi_gp = getGhostStencil_coarse(cross_dir, dir, side, diag_iv,
                                                                                Ecell_f_p, cfmask,
                                                                                cellid, cell_id_fine);
                                mycode::Stencil phi_gm = getGhostStencil_coarse(cross_dir, dir, side, diag_iv,
                                                                                Ecell_f_m, cfmask,
                                                                                cellid, cell_id_fine);

                                //! flux on E face = -(phi_ghost - E_fine)/dx_fine
                                phi_gp.addToColEntry(cell_id_fine(Ecell_f_p), -1.0);
                                phi_gm.addToColEntry(cell_id_fine(Ecell_f_m), -1.0);

                                mycode::Stencil Eflux = (phi_gp + phi_gm) / (2.0 * dxf[0]);

                                //! as coeff corresponding to col P contains -2*E_coeff
                                //! remove E_coeff as E cell is not present
                                st.addToColEntry(P, E_coeff);

                                //! add flux contro of E
                                //! here -ve sign is due to fine face flux 
                                st += Eflux * (-scale_fact / dx[0]);
                            }

                            /// E is cut cell
                            else if (mask(i + 1, j, k) == 2)
                            {
                                /// find the interceptData corresponding to E
                                for (auto &&solid : (*interfaces)[ilev])
                                {
                                    const auto &icpt_data = solid->getInterceptData()[mfi];

                                    if (!icpt_data.empty())
                                    {
                                        auto itr = std::find_if(icpt_data.begin(), icpt_data.end(),
                                                                [&](const InterceptData &idt) {
                                                                    return idt.cellid_ == right_iv;
                                                                });
                                        //if(diag_iv == test_iv1)
                                        //{
                                        //    std::cout<<"test_iv1 = "<<test_iv1[0]<<" , "<<test_iv1[1]<<'\n';
                                        //    std::cout<<"right_iv = "<<right_iv[0]<<" , "<<right_iv[1]<<'\n';
                                        //    //std::cout<<"itr = "<<*itr<<'\n';
                                        //}
                                        if (itr != icpt_data.end())
                                        {
                                            addGhostContro(domain, right_iv, st, cellid, E_coeff, *itr, mask);
                                            //if(diag_iv == test_iv1)
                                            //    std::cout<<"adding right_iv to ghost "<<'\n';
                                        }
                                    }
                                }
                            }

                            /// E is regular cell
                            else
                            {
                                st.addToColEntry(E, E_coeff);
                            }
                        }

                        /// N coeff
                        {
                            // N out
                            if (j + 1 > domain.bigEnd(1))
                            {
                                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                                amrex::Real y = prob_hi[1];
                                amrex::Real P_b = boundary_funcs[3](x, y, Time);
                                if (hibc[1] == amrex::LinOpBCType::Neumann)
                                {
                                    /// P_n = P_b * dy + P_p
                                    st.addToColEntry(P, N_coeff);
                                    st.addToSource(-P_b * dx[1] * N_coeff);
                                }
                                else if (hibc[1] == amrex::LinOpBCType::Dirichlet)
                                {
                                    /// P_n = 2.0 * P_b - P_p
                                    st.addToColEntry(P, -N_coeff);
                                    st.addToSource(-2.0 * P_b * N_coeff);
                                }
                                else
                                {
                                    amrex::Abort("AMG MLPoisson : invalid bc top");
                                }
                                //if(diag_iv == test_iv)
                                //{
                                //    amrex::Print()<<"cell = "<<diag_iv[0]<<" , "<<diag_iv[1]<<'\n';
                                //    amrex::Print()<<"P bc = "<<P_b<<'\n';
                                //    amrex::Print()<<"hibc[1] = "<<hibc[1]<<'\n';
                                //}
                            }

                            //! N is out of this level grid
                            //! as j+1 is not out of domain but top cell is outside grid
                            //! => coarse mesh @ top ... lev/lev-1 bndry
                            else if (!(*grid_)[ilev].contains(top_iv))
                            {
                                int dir = 1, cross_dir = 0, side = 1;
                                mycode::Stencil gst = getGhostStencil_fine(cross_dir, dir, side, top_iv, 
                                                                          (*grid_)[ilev], domain, cellid, cell_id_crse);

                                st += gst * N_coeff;  
                            }

                            //! W is covered by fine mesh
                            //! lev/lev+1 bndry
                            else if (cfmask(top_iv) == CFMask::covered)
                            {
                                int dir = 1, cross_dir = 0, side = 0;
                                const amrex::Real *dxf = (*geom_)[ilev + 1].CellSize();
                                amrex::IntVect Ncell_f_p(top_iv[0] * 2 + 1, top_iv[1] * 2);
                                amrex::IntVect Ncell_f_m(top_iv[0] * 2, top_iv[1] * 2);
                                mycode::Stencil phi_gp = getGhostStencil_coarse(cross_dir, dir, side, diag_iv,
                                                                                Ncell_f_p, cfmask,
                                                                                cellid, cell_id_fine);
                                mycode::Stencil phi_gm = getGhostStencil_coarse(cross_dir, dir, side, diag_iv,
                                                                                Ncell_f_m, cfmask,
                                                                                cellid, cell_id_fine);

                                //! flux on N face .. -(phi_ghost - N_fine)/dy_fine
                                phi_gp.addToColEntry(cell_id_fine(Ncell_f_p), -1.0);
                                phi_gm.addToColEntry(cell_id_fine(Ncell_f_m), -1.0);
                                mycode::Stencil Nflux = (phi_gp + phi_gm) / (2.0 * dxf[1]);

                                //! as coeff corresponding to col P contains -2*N_coeff
                                //! remove N_coeff as N cell is not present
                                st.addToColEntry(P, N_coeff);

                                //! add flux contro of N
                                st += Nflux * (-scale_fact / dx[1]);
                            }

                            /// N is cut cell
                            else if (mask(i, j + 1, k) == 2)
                            {
                                /// find the interceptData corresponding to N
                                for (auto &&solid : (*interfaces)[ilev])
                                {
                                    const auto &icpt_data = solid->getInterceptData()[mfi];

                                    if (!icpt_data.empty())
                                    {
                                        auto itr = std::find_if(icpt_data.begin(), icpt_data.end(),
                                                                [&](const InterceptData &idt) noexcept {
                                                                    return idt.cellid_ == top_iv;
                                                                });

                                        if (itr != icpt_data.end())
                                        {
                                            addGhostContro(domain, top_iv, st, cellid, N_coeff, *itr, mask);
                                        }
                                    }
                                }
                            }

                            /// N is regular cell
                            else
                            {
                                st.addToColEntry(N, N_coeff);
                            }
                        }

                        if (singular_ && P == 0)
                        {
                            int row = cellid(i, j, k);
                            amrex::Real b_val = 0.0;
                            amrex::Real x_val = phi(i, j, k);

                            int cols[] = {P};
                            double values[] = {1.0};
                            int nnz = 1;

                            solver->set_A_Val(&nnz, &row, cols, values);
                            solver->set_x_Val(&row, &x_val);
                            solver->set_b_Val(&row, &b_val);
                        }
                        else
                        {
                            int row = cellid(i, j, k);
                            amrex::Real b_val = st.getSource();
                            amrex::Real x_val = phi(i, j, k);

                            std::vector<int> cols = st.getColIdVec();
                            std::vector<amrex::Real> values = st.getCoeffVec();
                            int nnz = st.size();

                            solver->set_A_Val(&nnz, &row, cols.data(), values.data());
                            solver->set_x_Val(&row, &x_val);
                            solver->set_b_Val(&row, &b_val);
                        }
                    }
                });
            }
        }
    }

    solver->solve(soln_);
    residue_ = solver->getResidue();
    num_iterations = solver->getNIterations();
    

    amrex::Real stop_time = amrex::second() - start_time;
    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    amrex::PrintToFile("log") << "MLPoisson solve Time : " << stop_time << std::endl;
}

} // namespace mycode
