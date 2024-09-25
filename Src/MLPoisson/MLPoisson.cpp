#include "MLPoisson.H"
#include <AMReX_ParmParse.H>

namespace mycode
{

MLPoisson::MLPoisson
(
    MPI_Comm comm,
    amrex::Vector<amrex::Geometry> *geom,
    amrex::Vector<amrex::BoxArray> *grid,
    amrex::Vector<amrex::DistributionMapping> *dmap,
    amrex::Vector<std::vector<std::unique_ptr<Interface>>>* IF
)
:
comm_(comm),
geom_(geom),
grid_(grid),
dmap_(dmap),
interfaces(IF)
{
    /// read params
    amrex::ParmParse pp("LinearSolver");
    std::string lib = "hypre";
    pp.query("solver_lib", lib);
    if (lib == "HYPRE" || lib == "hypre" || lib == "Hypre")
    {
        solver_lib = LinearSolver::linear_solver_lib::HYPRE;
    }
    else if (lib == "PETSC" || lib == "petsc" || lib == "Petsc")
    {
        solver_lib = LinearSolver::linear_solver_lib::PETSC;
    }

    amrex::ParmParse pp1;
    pp1.query("Axisymmetric", isAxisymmetric);

    mask_.resize(geom_->size());
    sys_size_lev.resize(geom->size());
}

MLPoisson::~MLPoisson()
{}

void MLPoisson::getCoarseWts_fine
(
    int cross_dir,
    const amrex::IntVect &gcell, // fine ghost cell
    const amrex::BoxArray &fgrid,
    const amrex::Box& fdomain,
    amrex::Array<amrex::Real, 3> &wts,
    amrex::Array<amrex::IntVect, 3> &cells
)
{
    amrex::IntVect g_plus(gcell);
    g_plus.shift(cross_dir, 2);
    amrex::IntVect g_minus(gcell);
    g_minus.shift(cross_dir, -2);

    for (size_t i = 0; i < 3; i++)
    {
        cells[i] = gcell;
        cells[i].coarsen(2);
    }

    //! case 3 ... interpolation stencil shifted in -ve dir
    if (fgrid.contains(g_plus) || !fdomain.contains(g_plus))
    {
        cells[1].shift(cross_dir, -1);
        cells[2].shift(cross_dir, -2);
        //! at bottom
        if (gcell[cross_dir] % 2 == 0)
        {
            wts[0] = 21.0 / 32.0;
            wts[1] = 14.0 / 32.0;
            wts[2] = -3.0 / 32.0;
        }
        else
        {
            wts[0] = 45.0 / 32.0;
            wts[1] = -18.0 / 32.0;
            wts[2] = 5.0 / 32.0;
        }
    }
    //! case 2 ... interpolation stencil shifted in +ve dir
    else if (fgrid.contains(g_minus) || !fdomain.contains(g_minus))
    {
        cells[1].shift(cross_dir, 1);
        cells[2].shift(cross_dir, 2);
        //! at bottom
        if (gcell[cross_dir] % 2 == 0)
        {
            wts[0] = 45.0 / 32.0;
            wts[1] = -18.0 / 32.0;
            wts[2] = 5.0 / 32.0;            
        }
        else
        {
            wts[0] = 21.0 / 32.0;
            wts[1] = 14.0 / 32.0;
            wts[2] = -3.0 / 32.0;
        }
    }
    else
    {
        cells[0].shift(cross_dir, -1);
        cells[2].shift(cross_dir, 1);
        //! at bottom
        if (gcell[cross_dir] % 2 == 0)
        {
            wts[0] = 5.0 / 32.0;
            wts[1] = 30.0 / 32.0;
            wts[2] = -3.0 / 32.0;            
        }
        else
        {
            wts[0] = -3.0 / 32.0;
            wts[1] = 30.0 / 32.0;
            wts[2] = 5.0 / 32.0;
        }
    }
}

void MLPoisson::getCoarseWts_coarse
(
    int cross_dir,
    const amrex::IntVect &pcell,   // P cell
    const amrex::IntVect &pcell_f, // P cell fine
    amrex::Array4<int const> const &cfmask,
    amrex::Array<amrex::Real, 3> &wts,
    amrex::Array<amrex::IntVect, 3> &cells
)
{
    amrex::IntVect p_plus(pcell);
    p_plus.shift(cross_dir, 1);
    amrex::IntVect p_minus(pcell);
    p_minus.shift(cross_dir, -1);

    for (size_t i = 0; i < 3; i++)
    {
        cells[i] = pcell;
    }

    bool at_top = true;
    if (pcell_f[cross_dir] % 2 == 0)
    {
        at_top = false;
    }

    //! case 3 ... interpolation stencil shifted in -ve dir
    if (cfmask(p_plus) == CFMask::covered || cfmask(p_plus) == CFMask::out_of_grid)
    {
        cells[1].shift(cross_dir, -1);
        cells[2].shift(cross_dir, -2);
        
        if (at_top)
        {
            wts[0] = 45.0 / 32.0;
            wts[1] = -18.0 / 32.0;
            wts[2] = 5.0 / 32.0;
        }
        else
        {
            wts[0] = 21.0 / 32.0;
            wts[1] = 14.0 / 32.0;
            wts[2] = -3.0 / 32.0;
        }        
    }
    //! case 2 ... interpolation stencil shifted in +ve dir
    else if (cfmask(p_minus) == CFMask::covered || cfmask(p_minus) == CFMask::out_of_grid)
    {
        cells[1].shift(cross_dir, 1);
        cells[2].shift(cross_dir, 2);

        if (at_top)
        {
            wts[0] = 21.0 / 32.0;
            wts[1] = 14.0 / 32.0;
            wts[2] = -3.0 / 32.0;
        }
        else
        {
            wts[0] = 45.0 / 32.0;
            wts[1] = -18.0 / 32.0;
            wts[2] = 5.0 / 32.0;        
        }
    }
    else
    {
        cells[0].shift(cross_dir, -1);
        cells[2].shift(cross_dir, 1);

        if (at_top)
        {
            wts[0] = -3.0 / 32.0;
            wts[1] = 30.0 / 32.0;
            wts[2] = 5.0 / 32.0;
        }
        else
        {
            wts[0] = 5.0 / 32.0;
            wts[1] = 30.0 / 32.0;
            wts[2] = -3.0 / 32.0; 
        }
    }
}

mycode::Stencil MLPoisson::getGhostStencil_fine
(
    int cross_dir,
    int shift_dir,
    int side,                       // lo or hi (i.e. 0 => cell P is at lo face, 1 => at hi face) 
    const amrex::IntVect &gcell,    // fine ghost cell
    const amrex::BoxArray& fgrid,
    const amrex::Box& fdomain,
    amrex::Array4<int const> const& cell_id,
    amrex::Array4<int const> const& cell_id_crse
)
{
    amrex::Array<amrex::Real, 3> crse_wts;
    amrex::Array<amrex::IntVect, 3> crse_cells;
    getCoarseWts_fine(cross_dir, gcell, fgrid, fdomain, crse_wts, crse_cells);
    mycode::Stencil st;

    //! ghost cell stencil
    for (int n = 0; n < 3; n++)
    {
        st.addToColEntry(cell_id_crse(crse_cells[n]), (8.0/15.0) * crse_wts[n]);
    }

    amrex::IntVect Pcell(gcell), NbrCell(gcell);
    
    if (side == 0)
    {
        Pcell.shift(shift_dir, 1);
        NbrCell.shift(shift_dir, 2);
    }
    else
    {
        Pcell.shift(shift_dir, -1);
        NbrCell.shift(shift_dir, -2);
    }

    st.addToColEntry(cell_id(Pcell), 10.0 / 15.0);
    st.addToColEntry(cell_id(NbrCell), -3.0 / 15.0);

    st.setRowId(cell_id(Pcell));
    st.setSource(0.0);

    return st;
}

mycode::Stencil MLPoisson::getGhostStencil_coarse
(
    int cross_dir,
    int shift_dir,
    int side,                       // lo or hi (i.e. 0 => cell P is at lo face, 1 => at hi face) 
    const amrex::IntVect &pcell,    // P cell on coarse
    const amrex::IntVect &pcell_f,    // P cell fine
    amrex::Array4<int const> const &cfmask,
    amrex::Array4<int const> const& cell_id,
    amrex::Array4<int const> const& cell_id_fine
)
{
    amrex::Array<amrex::Real, 3> crse_wts;
    amrex::Array<amrex::IntVect, 3> curr_cells;
    getCoarseWts_coarse(cross_dir, pcell, pcell_f, cfmask, crse_wts, curr_cells);
    mycode::Stencil st;

    st.setRowId(cell_id(pcell));
    st.setSource(0.0);

    //! ghost cell stencil
    for (int n = 0; n < 3; n++)
    {
        st.addToColEntry(cell_id(curr_cells[n]), (8.0/15.0) * crse_wts[n]);
    }

    amrex::IntVect Pcell(pcell_f), NbrCell(pcell_f);
    
    if (side == 0)
    {
        NbrCell.shift(shift_dir, 1);
    }
    else
    {
        NbrCell.shift(shift_dir, -1);
    }
    st.addToColEntry(cell_id_fine(Pcell), 10.0 / 15.0);
    st.addToColEntry(cell_id_fine(NbrCell), -3.0 / 15.0);

    return st;
}

//! fill the fine level ghost cells using quadratic interpolation
//! this will be required for getting grad of solution
//! as in the projection of velocity
void MLPoisson::fillBoundaryFine()
{
    int nlevels = geom_->size();
    for (int ilev = 0; ilev < nlevels; ilev++)
    {
        soln_[ilev]->FillBoundary();
    }

    for (int ilev = 1; ilev < nlevels; ilev++)
    {
        //! need crse solution
        //! boxArray of fine level coarsened
        const amrex::BoxArray &fba = soln_[ilev]->boxArray();
        amrex::BoxArray cba(fba);
        cba.coarsen(2);

        //! dmap is same as fine level, thus desired data on each process
        //! as the coarse level covered cells are not the valid cells for 
        //! solution, only desired data is that of ghost cells
        int nghost = soln_[ilev]->nGrow();
        amrex::MultiFab crse_soln(cba, soln_[ilev]->DistributionMap(), 1, nghost);

        //! parallel copy
        crse_soln.copy(*soln_[ilev-1], 0, 0, 1, nghost, nghost);

        const amrex::Box &domain   = (*geom_)[ilev].Domain();
        auto dlo = amrex::lbound(domain);
        auto dhi = amrex::ubound(domain);
        amrex::Array<amrex::Real, 3> crse_wts;
        amrex::Array<amrex::IntVect, 3> crse_cells;

        for (amrex::MFIter mfi(*soln_[ilev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<amrex::Real> const &phif = soln_[ilev]->array(mfi);
            amrex::Array4<amrex::Real> const &phic = crse_soln.array(mfi);

            auto lo = amrex::lbound(bx);
            auto hi = amrex::ubound(bx);

            int k = 0;

            //! left side
            for (int j = lo.y; j <= hi.y; j++)
            {
                int i = lo.x;
                amrex::IntVect gcell(AMREX_D_DECL(i - 1, j, k));
                //! interpolate at c/f bndry only
                if (i > dlo.x && !(*grid_)[ilev].contains(gcell)) 
                {
                    int cross_dir = 1;
                    getCoarseWts_fine(cross_dir, gcell, (*grid_)[ilev], domain, crse_wts, crse_cells);

                    amrex::Real val = 0.0;
                    for (int n = 0; n < 3; n++)
                    {
                        val += (8.0/15.0) * crse_wts[n] * phic(crse_cells[n]);
                    }

                    val += (10.0 / 15.0) * phif(i, j, k) - (3.0 / 15.0) * phif(i + 1, j, k);
                    phif(gcell) = val;
                }
            }

            //! bottom
            for (int i = lo.x; i <= hi.x; i++)
            {
                int j = lo.y;
                amrex::IntVect gcell(AMREX_D_DECL(i, j - 1, k));
                if (j > dlo.y && !(*grid_)[ilev].contains(gcell))
                {
                    int cross_dir = 0;
                    
                    getCoarseWts_fine(cross_dir, gcell, (*grid_)[ilev], domain, crse_wts, crse_cells);

                    amrex::Real val = 0.0;
                    for (int n = 0; n < 3; n++)
                    {
                        val += (8.0/15.0) * crse_wts[n] * phic(crse_cells[n]);
                    }

                    val += (10.0 / 15.0) * phif(i, j, k) - (3.0 / 15.0) * phif(i, j + 1, k);
                    phif(gcell) = val;
                }
            }

            //! right side
            for (int j = lo.y; j <= hi.y; j++)
            {
                int i = hi.x;
                amrex::IntVect gcell(AMREX_D_DECL(i + 1, j, k));
                if (i < dhi.x && !(*grid_)[ilev].contains(gcell))
                {
                    int cross_dir = 1;
                    getCoarseWts_fine(cross_dir, gcell, (*grid_)[ilev], domain, crse_wts, crse_cells);

                    amrex::Real val = 0.0;
                    for (int n = 0; n < 3; n++)
                    {
                        val += (8.0/15.0) * crse_wts[n] * phic(crse_cells[n]);
                    }

                    val += (10.0 / 15.0) * phif(i, j, k) - (3.0 / 15.0) * phif(i - 1, j, k);
                    phif(gcell) = val;
                }
            }

            //! top
            for (int i = lo.x; i <= hi.x; i++)
            {
                int j = hi.y;
                amrex::IntVect gcell(AMREX_D_DECL(i, j + 1, k));
                if (j < dhi.y && !(*grid_)[ilev].contains(gcell))
                {
                    int cross_dir = 0;
                    getCoarseWts_fine(cross_dir, gcell, (*grid_)[ilev], domain, crse_wts, crse_cells);

                    amrex::Real val = 0.0;
                    for (int n = 0; n < 3; n++)
                    {
                        val += (8.0/15.0) * crse_wts[n] * phic(crse_cells[n]);
                    }

                    val += (10.0 / 15.0) * phif(i, j, k) - (3.0 / 15.0) * phif(i, j - 1, k);
                    phif(gcell) = val;
                }
            }
        }
    }
}

void MLPoisson::setBoundaryFunc(BoundaryFunc f)
{
    for (size_t i = 0; i < boundary_funcs.size(); i++)
    {
        boundary_funcs[i] = f;
    }
}

void MLPoisson::setBoundaryFunc(int i, BoundaryFunc f)
{
    boundary_funcs[i] = f;
}

} /*End namespace mycode */

