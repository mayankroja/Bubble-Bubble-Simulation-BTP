#include "LinearSolver.H"
#include <AMReX_ParmParse.H>

namespace mycode
{

LinearSolver::~LinearSolver()
{
    Mat_A = nullptr;
    Vec_x = nullptr;
    Vec_b = nullptr;
}

void LinearSolver::define
(
		MPI_Comm comm,
		amrex::Vector<amrex::Geometry> *geom,
		amrex::Vector<amrex::BoxArray> *grid,
		amrex::Vector<amrex::DistributionMapping> *dmap,
		amrex::Vector<CFMask> *cfmask,
		MLTraverseIndex *tri
)
{
    int n = geom->size();
    amrex::Vector<amrex::iMultiFab *> mask(n);
    for (size_t i = 0; i < n; i++)
    {
        mask[i] = nullptr;
    }
    define(comm, geom, grid, dmap, cfmask, mask, tri);
}

void LinearSolver::define
(
    MPI_Comm comm,
    amrex::Vector<amrex::Geometry> *geom,
    amrex::Vector<amrex::BoxArray> *grid,
    amrex::Vector<amrex::DistributionMapping> *dmap,
    amrex::Vector<CFMask> *cfmask,
    amrex::Vector<amrex::iMultiFab*>& mask,
    MLTraverseIndex *tri
)
{
    nlevels = geom->size();
    MaskPtr_.resize(nlevels);

    geom_ = geom;
    grids_ = grid;
    dmap_ = dmap;
    cfmask_ = cfmask;

    comm_ = comm;

    for (size_t ilev = 0; ilev < nlevels; ilev++)
    {
        MaskPtr_[ilev] = mask[ilev];

    }

    if (tri)
    {
        tridx = tri;    
    }
    else
    {
        amrex::ParmParse pp;
        int nghost = 1;
        pp.query("Nghost", nghost);
        tridx->define(*geom_, *grids_, *dmap_, nghost);
    }

    amrex::ParmParse pp("LinearSolver");

    pp.get("max_iter", max_iter);
    pp.get("rel_tol", rel_tolerance);
    pp.get("abs_tol", abs_tolerance);
    pp.query("print_system", print_system);
}

void LinearSolver::solve(amrex::Vector<amrex::MultiFab*> &soln)
{
    amrex::Real start_time = amrex::second();
    assembleSystem();
    solverSetupAndSolve();
    getSolution(soln);

    amrex::Real stop_time = amrex::second() - start_time;
    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::ReduceRealMax(stop_time, IOProc);

    amrex::PrintToFile("log") << "LinearSolver:: solve Time : " << stop_time << std::endl;
}

void LinearSolver::getSolution(amrex::Vector<amrex::MultiFab*> &soln)
{
    amrex::Array4<int const> mask;
    amrex::IArrayBox mask_fab;

    for (int ilev = 0; ilev < nlevels; ilev++)
    {
        for (amrex::MFIter mfi(*soln[ilev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            int npts = bx.numPts();
            std::vector<int> rows;
            rows.reserve(npts);

            amrex::Array4<int const> const &cellid = tridx->getTrIndex(ilev).const_array(mfi);
            amrex::Array4<amrex::Real> const &soln_fab = soln[ilev]->array(mfi);
            amrex::Array4<int const> const &cfmask = (*cfmask_)[ilev].Mask().const_array(mfi);

            if (MaskPtr_[ilev])
            {
                mask = MaskPtr_[ilev]->const_array(mfi);
            }
            else
            {
                mask_fab.resize(bx);
                mask_fab.setVal(1);
                mask = mask_fab.const_array();
            }

            int nvalues = 0;
            amrex::Loop(bx, [&](int i, int j, int k) 
            {
                if (mask(i, j, k) == 1 && cfmask(i,j,k) != CFMask::covered)
                {
                    rows.emplace_back(cellid(i, j, k));
                    nvalues++;
                }
            });

            rows.shrink_to_fit();
            std::vector<amrex::Real> xfab(nvalues);

            getVecVal(Vec_x, nvalues, rows.data(), xfab.data());

            int n = 0;
            amrex::Loop(bx, [&](int i, int j, int k) 
            {
                if (mask(i, j, k) == 1 && cfmask(i,j,k) != CFMask::covered)
                {
                    soln_fab(i, j, k) = xfab[n];
                    n++;
                }
            });
        }
    }
}

} /*End namespace mycode */

