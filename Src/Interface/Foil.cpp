#include "Foil.H"
#include <AMReX_ParmParse.H>

namespace mycode
{

Foil::Foil
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

    pp.get("thetacp", thetacp);
    pp.get("dtheta", dtheta);
    pp.get("Amp", Amp);
    pp.get("Freq", Freq);
    pp.get("Phase", Phase);
    pp.get("rhoArea", rhoArea);
    pp.get("Lp", Lp);
    pp.get("L0", L0);
    pp.get("PreFactor", PreFactor);
    pp.get("N_Shape_1", N_Shape_1);
    pp.get("N_Shape_2", N_Shape_2);

    thetacp *= M_PI / 180;
    Amp *= M_PI / 180;
    Phase *= M_PI / 180;

    Len_foil = 1.0 / (1.0 + std::pow(0.5, 0.5 * (2.0 - N_Shape_2) / N_Shape_1));
    xcp0 = xcp;
    ycp0 = ycp;

    is_advect_ls = false;
}

Foil::~Foil() {}

amrex::RealArray Foil::getOrigin()
{
    amrex::RealArray origin;
    amrex::Real L_origin = pow(2.0, ((N_Shape_2 - 2.0) / (2.0 * N_Shape_1)));
    origin[0] = xc0 - L_origin * cos(thetacp);
    origin[1] = yc0 - L_origin * sin(thetacp);
    return origin;
}

void Foil::PrescribeLevelSetMotion(const amrex::Real& t)
{
    thetacp += Amp * sin(Freq * 2.0 * M_PI * t + Phase);
    dxcpdt = 0.0;
    dycpdt = 0.0;
    dthetacpdt = Freq * 2.0 * M_PI * Amp * cos(Freq * 2.0 * M_PI * t + Phase);
    d2thetacpdt2 = -(Freq * 2.0 * M_PI) * (Freq * 2.0 * M_PI) * thetacp;

    ComputePsi();
}

void Foil::ComputeLevelSetMotion(const amrex::Real &t, const amrex::Real &dt, int RKStage)
{
    thetacp = Amp * sin(Freq * 2.0 * M_PI * t + Phase);
    dthetacpdt = Freq * 2.0 * M_PI * Amp * cos(Freq * 2.0 * M_PI * t + Phase);
    d2thetacpdt2 = -(Freq * 2.0 * M_PI) * (Freq * 2.0 * M_PI) * thetacp;

    if (DOF > AMREX_SPACEDIM || DOF < 0 )
    {
        amrex::Abort("DOF should be <= space dim");
    }

    if (DOF > 0)
    {
        if (RKStage == 1)
        {
            xcp0 = xcp;
            dxcpdt0 = dxcpdt;

            if (DOF > 1)
            {
                ycp0 = ycp;
                dycpdt0 = dycpdt;
            }            
        }

        xcp += dt * dxcpdt;
        dxcpdt += dt * (FPx + FVx) / rhoArea;
        dxcpdt += dt * (Lp * sin(thetacp) * d2thetacpdt2 + Lp * cos(thetacp) * dthetacpdt * dthetacpdt);

        if (DOF > 1)
        {
            ycp += dt * dycpdt;
            dycpdt += dt * (FPy + FVy) / rhoArea;
            dycpdt += dt * (-Lp * cos(thetacp) * d2thetacpdt2 + Lp * sin(thetacp) * dthetacpdt * dthetacpdt);
        }

        if (RKStage == 2)
        {
            xcp = 0.5 * (xcp + xcp0);
            dxcpdt = 0.5 * (dxcpdt + dxcpdt0);

            if (DOF > 1)
            {
                ycp = 0.5 * (ycp + ycp0);
                dycpdt = 0.5 * (dycpdt + dycpdt0);
            }            
        }
    }

    ComputePsi();
}

void Foil::ComputePsi()
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
            xc0 = xcm - L0 * cos(thetacp);
            yc0 = ycm - L0 * sin(thetacp);

            amrex::Real x2 = x - xc0;
            amrex::Real y2 = y - yc0;
		
            // Rotation
            amrex::Real x1 = x2 * cos(thetacp) + y2 * sin(thetacp);
            amrex::Real y1 = y2 * cos(thetacp) - x2 * sin(thetacp);

            amrex::Real r1 = std::hypot(x1, y1);
            amrex::Real alpha1 = std::asin( std::fabs(y1)/(r1 + 1.0E-15));

            if(x1 < 0.0 && y1 < 0.0)
                alpha1 = alpha1 + M_PI;
            else if (x1 < 0.0 && y1 >= 0.0)
                alpha1 = M_PI - alpha1;
            else if(x1 >= 0.0 && y1 < 0.0)
                alpha1 = 2.0 * M_PI - alpha1;
                
            amrex::Real L1 = Len_foil / (std::pow(std::fabs(cos(alpha1 / 4.0)), N_Shape_2) + std::pow(std::fabs(sin(alpha1 / 4.0)), N_Shape_2));
            amrex::Real Den = 0.5 * (tanh(PreFactor * (alpha1 - dtheta)) - tanh(PreFactor * (alpha1 - (2.0 * M_PI - dtheta))));
            L1 = Den * (std::pow(L1, 0.5) - std::pow(Len_foil, 0.5)) + std::pow(Len_foil, 0.5);

            r1 = std::pow(r1, 0.5 * N_Shape_1) - L1;
            f(i, j, k) = -r1;
        });
    }

    //Reinit_algoim();
    Reinit();
    PhaseFieldBC();
}

void Foil::TubeIdentification()
{
    amrex::Abort("Foil::TubeIdentification is specific to advected interfaces..");
}

void Foil::TubeAdvectLevelSet(const amrex::MultiFab &u, const amrex::MultiFab &v, const amrex::Real &dt)
{
    amrex::Abort("Foil::TubeAdvectLevelSet is specific to advected interfaces..");
}

std::unique_ptr<Foil> MakeFoil(const amrexMesh &mesh,
							   const std::string &name)
{
    auto f = std::make_unique<Foil>(mesh, name);
    return f;
}

/*void Foil::Regularization()
{
    const amrex::Real *dx = mesh_.geometry().CellSize();
    for (amrex::MFIter mfi(Mask_); mfi.isValid(); ++mfi)
    {
        amrex::Array4<amrex::Real> const &Psi = psi.array(mfi);
        auto &index = Index_[mfi];
        for (auto &&cell : index)
        {
            int i = cell[0], j = cell[1], k = 0;
            if(cell[1] == 0)
            {
                //amrex::Print()<<"cell = "<<cell<<'\n';
                //amrex::Print()<<"Psi(cell) = "<<Psi(cell)<<'\n';

                if(Psi(cell) < 0.0 && std::abs(Psi(cell)) < .5*dx[1])
                {
                    if(Psi(i, 1 , 0) > 0.0 )
                    {
                        amrex::Print()<<"cell to be regularized = "<<cell<<" , psi = "<<Psi(cell)<<'\n';
                        Psi(cell) = -0.1*Psi(cell);
                    }
                }
            }
        }
    }

}
*/

} /*End namespace mycode */

