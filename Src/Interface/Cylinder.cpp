#include "Cylinder.H"
#include <AMReX_ParmParse.H>

namespace mycode
{

Cylinder::Cylinder
(
    const amrexMesh &mesh,
    const std::string& name
) 
: 
Interface(mesh, name)
{
    amrex::ParmParse pp(name);
    amrex::Vector<amrex::Real> cp(AMREX_SPACEDIM);
    pp.getarr("cp", cp);
    xcp = cp[0];
    ycp = cp[1];
    xcp0 = xcp;
    ycp0 = ycp;

    pp.get("thetacp", thetacp);
    pp.get("dtheta", dtheta);
    pp.get("Amp", Amp);
    pp.get("Freq", Freq);
    pp.get("Lp", Lp);
    pp.get("L0", L0);
    pp.get("radius", radius);

    vel.resize(AMREX_SPACEDIM);
    pp.getarr("vel", vel);

    thetacp *= M_PI / 180;
    Amp *= M_PI / 180;

    is_advect_ls = false;
}

Cylinder::~Cylinder() {}

//TODO : origin is not supposed to be center or pivot.. it is point on the interface
//TODO : for now not used in cylinder case
amrex::RealArray Cylinder::getOrigin()
{
    amrex::RealArray origin;
    origin[0] = xcp;
    origin[1] = ycp;
    return origin;
}

void Cylinder::PrescribeLevelSetMotion(const amrex::Real& t)
{
    thetacp += 0.0 * Amp * sin(Freq * 2.0 * M_PI * t + 0.5 * M_PI);
    dxcpdt = 0.0;
    dycpdt = 0.0;
    dthetacpdt = 0.0 * Freq * 2.0 * M_PI * Amp * cos(Freq * 2.0 * M_PI * t + 0.5 * M_PI);
    d2thetacpdt2 = -0.0 * (Freq * 2.0 * M_PI) * (Freq * 2.0 * M_PI) * thetacp;

    ComputePsi();
}

void Cylinder::ComputeLevelSetMotion(const amrex::Real &t, const amrex::Real &dt, int RKStage)
{
    thetacp = 0.0 * Amp * sin(Freq * 2.0 * M_PI * t + 0.5 * M_PI);
    dxcpdt = vel[0];
    dycpdt = vel[1];
    dthetacpdt = 0.0 * Freq * 2.0 * M_PI * Amp * cos(Freq * 2.0 * M_PI * t + 0.5 * M_PI);
    d2thetacpdt2 = -0.0*(Freq * 2.0 * M_PI) * (Freq * 2.0 * M_PI) * thetacp;

    xcp = xcp0 + dxcpdt * t;
    ycp = ycp0 + dycpdt * t;

    ComputePsi();
}

void Cylinder::ComputePsi()
{
    const amrex::Real *prob_lo = mesh_.geometry().ProbLo();
    const amrex::Real *dx = mesh_.geometry().CellSize();

    for(amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const &f = psi.array(mfi);

        amrex::ParallelFor(bx,
        [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
        {
            amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
            amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);

            amrex::Real xcm = xcp + Lp * cos(thetacp);
            amrex::Real ycm = ycp + Lp * sin(thetacp);
            amrex::Real xc0 = xcm - L0 * cos(thetacp);
            amrex::Real yc0 = ycm - L0 * sin(thetacp);

            amrex::Real x2 = x - xc0;
            amrex::Real y2 = y - yc0;

            // Rotation
            amrex::Real x1 = x2 * cos(thetacp) + y2 * sin(thetacp);
            amrex::Real y1 = y2 * cos(thetacp) - x2 * sin(thetacp);

            amrex::Real r1 = std::hypot(x1, y1);
            r1 = r1 - radius;

            f(i, j, k) = -r1;
        });
    }

    //Reinit();
    PhaseFieldBC();
}

void Cylinder::TubeIdentification()
{
    amrex::Abort("Cylinder::TubeIdentification is specific to advected interfaces..");
}

/*void Cylinder::Regularization()
{
    amrex::Abort("No regularization for cylinder");
}
*/

void Cylinder::TubeAdvectLevelSet(const amrex::MultiFab &u, const amrex::MultiFab &v, const amrex::Real &dt)
{
    amrex::Abort("Cylinder::TubeAdvectLevelSet is specific to advected interfaces..");
}

std::unique_ptr<Cylinder> MakeCylinder
(
    const amrexMesh &mesh,
    const std::string &name
)
{
    auto f = std::make_unique<Cylinder>(mesh, name);
    return f;
}

} /*End namespace mycode */

