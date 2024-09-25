#include "Bubble.H"
#include <AMReX_ParmParse.H>

namespace mycode
{

Bubble::Bubble
(
    const amrexMesh &mesh,
    const std::string& name
) 
: 
InterfaceAdvect(mesh, name)
{
    //amrex::Print()<<"Bubble constructor"<<'\n';
    amrex::ParmParse pp(name);
    amrex::Vector<amrex::Real> cp(AMREX_SPACEDIM);
    pp.getarr("cp", cp);
    xcp = cp[0];
    ycp = cp[1];

    pp.get("radius", radius);
    pp.get("is_fluid_in", is_fluid_in);
    pp.get("P_interface0",P_interface0);
    pp.get("is_advect_ls",is_advect_ls);
    pp.query("body_force",bdy_frc);
    pp.query("body_force_type",bdy_frc_type);
    P_interface = P_interface0;
}

Bubble::~Bubble() {}

//TODO : origin is not supposed to be center or pivot.. it is point on the interface
//TODO : for now not used in Bubble case
amrex::RealArray Bubble::getOrigin()
{
    amrex::RealArray origin;
    origin[0] = xcp;
    origin[1] = ycp;
    return origin;
}

void Bubble::PrescribeLevelSetMotion(const amrex::Real& t)
{
    ComputePsi();
}

void Bubble::ComputeLevelSetMotion(const amrex::Real &t, const amrex::Real &dt, int RKStage)
{
    amrex::Abort("Bubble: ComputeLevelSetMotion(...) is specific to prescribed level set.. use TubeAdvectLevelSet()");
}

void Bubble::ComputePsi()
{
    //amrex::Print()<<"Making new bubble"<<'\n';
    const amrex::Real *prob_lo = mesh_.geometry().ProbLo();
    const amrex::Real *dx = mesh_.geometry().CellSize();

    int m_sign(is_fluid_in ? 1.0 : -1.0);

    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const &f = psi.array(mfi);

        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
        {
            amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
            amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);

            amrex::Real x1 = x - xcp;
            amrex::Real y1 = y - ycp;

            amrex::Real r1 = std::hypot(x1, y1);
            r1 = r1 - radius;

            f(i, j, k) = m_sign * r1;
	    //f(i, j, k) = m_sign * x1;
        });
    }

    //Reinit();
    PhaseFieldBC();
}

std::unique_ptr<Bubble> MakeBubble
(
    const amrexMesh &mesh,
    const std::string &name
)
{
    //amrex::Print()<<"Going to make new bubble"<<'\n';
    auto f = std::make_unique<Bubble>(mesh, name);
    return f;
}

} /*End namespace mycode */

