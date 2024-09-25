#include <MLTraverseIndex.H>

namespace mycode
{

MLTraverseIndex::MLTraverseIndex()
{}

MLTraverseIndex::MLTraverseIndex
(
    const amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<amrex::BoxArray>& grid,
    const amrex::Vector<amrex::DistributionMapping>& dmap,
    int nghost
)
{
    int nlevels = geom.size();
    amrex::Vector<amrex::iMultiFab *> mask(nlevels);

    for (size_t i = 0; i < nlevels; i++)
    {
        mask[i] = nullptr;
    }

    define(geom, grid, dmap, mask, nghost);
}

void MLTraverseIndex::define
(
    const amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<amrex::BoxArray>& grid,
    const amrex::Vector<amrex::DistributionMapping>& dmap,
    int nghost
)
{
    int nlevels = geom.size();
    amrex::Vector<amrex::iMultiFab *> mask(nlevels);

    for (size_t i = 0; i < nlevels; i++)
    {
        mask[i] = nullptr;
    }

    define(geom, grid, dmap, mask, nghost);
}

MLTraverseIndex::MLTraverseIndex
(
    const amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<amrex::BoxArray>& grid,
    const amrex::Vector<amrex::DistributionMapping>& dmap,
	amrex::Vector<amrex::iMultiFab*> mask,
    int nghost
)
{
    define(geom, grid, dmap, mask, nghost);
}

MLTraverseIndex::~MLTraverseIndex()
{}

void MLTraverseIndex::define
(
    const amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<amrex::BoxArray>& grid,
    const amrex::Vector<amrex::DistributionMapping>& dmap,
    amrex::Vector<amrex::iMultiFab*> Pmask,
    int nghost
)
{
    int num_procs, myid;
    MPI_Comm comm_ = MPI_COMM_WORLD;

    MPI_Comm_size(comm_, &num_procs);
    MPI_Comm_rank(comm_, &myid);

    int nlevels = geom.size();
    int finest_lev = nlevels - 1;

    cell_id_.resize(nlevels);
    amrex::Vector<amrex::LayoutData<int>> nCells_grid(nlevels); // one box: total cells

    for (int lev = 0; lev < nlevels; lev++)
    {

        nCells_grid[lev].define(grid[lev], dmap[lev]);         // no of pts in a box
        cell_id_[lev].define(grid[lev], dmap[lev], 1, nghost); // index for cells over all boxes
    }

    int nCells_proc = 0; // no of cells per processor
    for (int lev = 0; lev < nlevels; lev++)
    {
#ifdef _OPENMP
#pragma omp parallel reduction(+:nCells_proc)
#endif
    {
        amrex::BaseFab<int> ifab; // index in a box
        for (amrex::MFIter mfi(cell_id_[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::BaseFab<int> &cellid_fab = cell_id_[lev][mfi];
            cellid_fab.setVal(std::numeric_limits<int>::lowest());
            {
                ifab.resize(bx);
                int *p = ifab.dataPtr();
                int npts = 0;

                if (Pmask[lev] != nullptr)
                {
                    amrex::Array4<int const> const &mask = Pmask[lev]->const_array(mfi);
                    if (lev < finest_lev)
                    {
                        amrex::Loop(bx, [&, p](int i, int j, int k) mutable
                        {
                            //! not solving on cells covered by fine cells
                            amrex::IntVect fcell(AMREX_D_DECL(2 * i, 2 * j, 2 * k));
                            if (mask(i, j, k) == 1 && !grid[lev + 1].contains(fcell))
                            {
                                *p = npts;
                                npts++;
                            }
                            p++;
                        });
                    }
                    else
                    {
                        amrex::Loop(bx, [&, p](int i, int j, int k) mutable 
                        {
                            if (mask(i, j, k) == 1)
                            {
                                *p = npts;
                                npts++;
                            }
                            p++;
                        });
                    }
                }
                else
                {
                    if (lev < finest_lev)
                    {
                        amrex::Loop(bx, [&, p](int i, int j, int k) mutable 
                        {
                            //! not solving on cells covered by fine cells
                            amrex::IntVect fcell(AMREX_D_DECL(2 * i, 2 * j, 2 * k));
                            if (!grid[lev + 1].contains(fcell))
                            {
                                *p = npts;
                                npts++;
                            }
                            p++;
                        });
                    }
                    else
                    {
                        amrex::Loop(bx, [&, p](int i, int j, int k) mutable 
                        {
                            *p = npts;
                            npts++;
                            p++;
                        });
                    }
                }
                cellid_fab.copy(ifab, bx);
                nCells_grid[lev][mfi] = npts;
                nCells_proc += npts;
            }
        }
    }

    }
    

    amrex::Vector<int> nNodes_allprocs(num_procs);
    MPI_Allgather(&nCells_proc, sizeof(int), MPI_CHAR,
                  nNodes_allprocs.data(), sizeof(int), MPI_CHAR,
                  comm_);
    proc_begin_ = 0;
    for (int i = 0; i < myid; ++i) 
    {
        proc_begin_ += nNodes_allprocs[i];
    }

    nNodes_world_ = 0;
    for (auto i : nNodes_allprocs) {
        nNodes_world_ += i;
    }
    proc_end_ = proc_begin_;

    offset_.resize(nlevels);
    nNodes_lev.resize(nlevels);
    for (int lev = 0; lev < nlevels; lev++)
    {
        nNodes_lev[lev] = 0;
        offset_[lev].define(grid[lev], dmap[lev]);

        for (amrex::MFIter mfi(nCells_grid[lev]); mfi.isValid(); ++mfi)
        {
            offset_[lev][mfi] = proc_end_;
            proc_end_ += nCells_grid[lev][mfi];
            nNodes_lev[lev] += nCells_grid[lev][mfi];
        }

        amrex::ParallelDescriptor::ReduceIntSum(nNodes_lev[lev]);
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(proc_end_ == proc_begin_+nCells_proc,
                                     "MLTraverseIndex::define: how did this happen?");

    for (int lev = 0; lev < nlevels; lev++)
    {

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (amrex::MFIter mfi(cell_id_[lev],true); mfi.isValid(); ++mfi)
    {
        cell_id_[lev][mfi].plus(offset_[lev][mfi], mfi.tilebox());
    }
    cell_id_[lev].FillBoundary(geom[lev].periodicity());
    
    }
}

std::unique_ptr<amrex::iMultiFab> MLTraverseIndex::getTrIndex(int to_lev, int from_lev)
{
    std::unique_ptr<amrex::iMultiFab> tridx;
    int nghost = cell_id_[to_lev].nGrow();
    const amrex::DistributionMapping &dmap = cell_id_[to_lev].DistributionMap();
    //! 'from_lev' is a coarse level
    if (to_lev > from_lev)
    {
        //! boxArray of fine level coarsened
        const amrex::BoxArray &fba = cell_id_[to_lev].boxArray();
        amrex::BoxArray cba(fba);
        cba.coarsen(2);

        //! dmap is same as fine level, thus desired data on each process
        //! as the coarse level covered cells are not the valid cells for 
        //! solution, only desired data is that of ghost cells
        tridx = std::make_unique<amrex::iMultiFab>(cba, dmap, 1, nghost);

        //! parallel copy
        tridx->copy(cell_id_[from_lev], 0, 0, 1, nghost, nghost);
    }
    //! 'from_lev' is fine
    else
    {
        // const amrex::BoxArray &fba = cell_id_[from_lev].boxArray();      
        //! refine the coarse level box array
        const amrex::BoxArray &cba = cell_id_[to_lev].boxArray();   
        amrex::BoxArray fba(cba);
        fba.refine(2);

        //! construct iMultiFab with refined boxArray of coarse
        //! and dmap with coarse... 
        //! this will ensure that the data at c/f interface is the desired data from fine
        tridx = std::make_unique<amrex::iMultiFab>(fba, dmap, 1, nghost);

        //! parallel copy
        //! this will copy the data from fine level onto coarse level
        //! as the boxArray is refined coarse level boxArray
        //! it will have more cells than that of fine level ( since grids are nested properly )
        //! the cells which have no corresponding cells on fine level will be filled zero
        tridx->copy(cell_id_[from_lev], 0, 0, 1, nghost, nghost);
    }
    return tridx;
}

} // End namespace mycode

