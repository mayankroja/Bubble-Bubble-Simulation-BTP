#include "amrexMesh.H"
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>

#include <AMReX_PlotFileUtil.H>

namespace mycode
{

/**
 * @brief Construct a new amrexMesh object using reference
 *        to BoxArray, Geometry and DistributionMapping
 * 
 */
amrexMesh::amrexMesh(const amrex::BoxArray& grid,
			  const amrex::Geometry& geom,
			  const amrex::DistributionMapping& dmap)
:
grid_(grid),
geom_(geom),
dmap_(dmap)
{
//	amrex::Print()<<"amrexMesh constructor 1"<<'\n';
}               


/**
 * @brief Construct a new amrexMesh object
 * 
 */
amrexMesh::amrexMesh()
{
    amrex::Print()<<"amrexMesh constructor 2"<<'\n';
    amrex::ParmParse pp("mesh");

    amrex::Vector<int> n_cell, max_grid_size;
    amrex::Vector<amrex::Real> real_lo, real_hi;

    n_cell.resize(AMREX_SPACEDIM);
    max_grid_size.resize(AMREX_SPACEDIM);
    real_lo.resize(AMREX_SPACEDIM);
    real_hi.resize(AMREX_SPACEDIM);

    pp.getarr("n_cell", n_cell, 0, AMREX_SPACEDIM);
    pp.getarr("max_grid_size", max_grid_size, 0, AMREX_SPACEDIM);
    pp.getarr("real_lo", real_lo, 0, AMREX_SPACEDIM);
    pp.getarr("real_hi", real_hi, 0, AMREX_SPACEDIM);

    amrex::IntVect dom_lo(AMREX_D_DECL(0, 0, 0));
    amrex::IntVect dom_hi(AMREX_D_DECL(n_cell[0] - 1, n_cell[1] - 1, n_cell[2] - 1));
    amrex::Box domain(dom_lo, dom_hi);
    grid_.define(domain);

    amrex::IntVect chop(AMREX_D_DECL(max_grid_size[0], max_grid_size[1], max_grid_size[2]));
    grid_.maxSize(chop);

    amrex::RealBox real_box({AMREX_D_DECL(real_lo[0], real_lo[1], real_lo[2])},
                            {AMREX_D_DECL(real_hi[0], real_hi[1], real_hi[2])});

    amrex::Vector<int> is_periodic(AMREX_SPACEDIM, 0);
    pp.queryarr("is_periodic", is_periodic, 0, AMREX_SPACEDIM);

    geom_.define(domain, &real_box, amrex::CoordSys::cartesian, is_periodic.data());
    dmap_.define(grid_);
}

amrexMesh::~amrexMesh()
{
}

namespace
{
    // utility to skip to next line in Header
    void GotoNextLine(std::istream &is)
    {
        constexpr std::streamsize bl_ignore_max{100000};
        is.ignore(bl_ignore_max, '\n');
    }

    void GetCheckFileName(std::string &check_file)
    {
        //! check status in InputR/status.dat and Input/status.dat
        int stat = 0, statR = 0;
        amrex::Real T = 0.0, TR = 0.0;
        // Input/status
        {
            std::string File("Input/status.dat");
            std::string word;
            std::ifstream dataFile(File);

            std::string lineContents;
            while (getline(dataFile, lineContents))
            {
                std::stringstream ss(lineContents);
                ss >> word; // time
                T = std::stod(word);
                ss >> word; // status
                stat = std::stoi(word);
            }
            dataFile.close();
        }
        // InputR/status
        {
            std::string File("InputR/status.dat");
            std::string word;
            std::ifstream dataFile(File);

            std::string lineContents;
            while (getline(dataFile, lineContents))
            {
                std::stringstream ss(lineContents);
                ss >> word; // time
                TR = std::stod(word);
                ss >> word; // status
                statR = std::stoi(word);
            }
            dataFile.close();
        }

        if (stat == 1 && statR == 1)
        {
            if (T > TR)
            {
                check_file = "Input/chk";
            }
            else
            {
                check_file = "InputR/chk";
            }
        }
        else if (stat == 1)
        {
            check_file = "Input/chk";
        }
        else if (statR == 1)
        {
            check_file = "InputR/chk";
        }
        else
        {
            amrex::Abort("chk files are currupted.. can not restart");
        }
    }

} // namespace

void amrexMesh::ReadFromCheckFile(std::string& chk_file)
{
    GetCheckFileName(chk_file);
    amrex::PrintToFile("log") << "Read boxArray from checkpoint " << chk_file << "\n";

    // Header
    std::string File(chk_file + "/Header");

    amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::GetIOBufferSize());

    amrex::Vector<char> fileCharPtr;
    amrex::ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    int finest_level;
    is >> finest_level;
    GotoNextLine(is);

    int Iter;
    amrex::Real dt, Time;
    // read in istep
    is >> Iter;

    // read in dt
    is >> dt;

    // read in time
    is >> Time;

    // read BoxArray
    amrex::BoxArray ba;
    ba.readFrom(is);
    grid_ = ba;

    amrex::DistributionMapping dm{ba, amrex::ParallelDescriptor::NProcs()};
    dmap_ = dm;
}

} // namespace amrexH
