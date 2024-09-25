#include "CFMask.H"

namespace mycode
{

CFMask::CFMask()
{
}

CFMask::CFMask
(
    int lev,
    const amrex::Vector<amrex::Geometry> &geom,
    const amrex::Vector<amrex::BoxArray> &grid,
    const amrex::Vector<amrex::DistributionMapping> &dmap,
    int nghost
)
{
    define(lev, geom, grid, dmap, nghost);
}

void CFMask::define
(
    int lev,
    const amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<amrex::BoxArray> &grid,
    const amrex::Vector<amrex::DistributionMapping> &dmap,
    int nghost
)
{
    int finest_level = grid.size() - 1;
    mask_.define(grid[lev], dmap[lev], 1, nghost);

    const amrex::BoxArray &cba = grid[lev];
    for (amrex::MFIter mfi(mask_); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.growntilebox();
        amrex::Array4<int> const &mask = mask_.array(mfi);

        if (lev < finest_level)
        {
            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                amrex::IntVect cell(AMREX_D_DECL(i, j, k));
                amrex::IntVect fcell(AMREX_D_DECL(2*i, 2*j, 2*k)); //! assuming refinment ratio 2
                if (!cba.contains(cell))
                {
                    mask(i, j, k) = out_of_grid;
                }
                else if (grid[lev+1].contains(fcell))
                {
                    mask(i, j, k) = covered;
                }
                else
                {
                    mask(i, j, k) = valid;
                }
            });
        }
        else
        {
            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                amrex::IntVect cell(AMREX_D_DECL(i, j, k));
                if (!cba.contains(cell))
                {
                    mask(i, j, k) = out_of_grid;
                }
                else
                {
                    mask(i, j, k) = valid;
                }
            });
        }
    }
    //! getermine c/f bndry cell
    if (lev < finest_level)
    {
        // amrex::AllPrintToFile("log") << "lev " << lev << "\n";
        for (amrex::MFIter mfi(mask_); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<int> const &mask = mask_.array(mfi);

            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
            {
                if (mask(i, j, k) == valid)
                {
                    //! check 4 nbrs
                    //! if any of them is covered => at c/f
                    if (
                        mask(i + 1, j, k) == covered ||
                        mask(i, j + 1, k) == covered ||
                        mask(i - 1, j, k) == covered ||
                        mask(i, j - 1, k) == covered)
                    {
                        mask(i, j, k) = cfbndry;
                    }
                }
            });
        }
    }

    mask_.FillBoundary(geom[lev].periodicity());
}

CFMask::~CFMask()
{
}

} /*End namespace mycode */

