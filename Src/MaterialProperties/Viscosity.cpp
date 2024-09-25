#include "Viscosity.H"
#include <AMReX_ParmParse.H>

namespace mycode
{

Viscosity::Viscosity()
{
    amrex::ParmParse pp("viscosity");
    pp.get("lambda", lambda);
    pp.get("Mu_max", Mu_max);
    pp.get("Mu_min", Mu_min);
    pp.get("N", N);
    pp.query("shear mdulus",eta);

    amrex::ParmParse pp1("viscosity_polymer");
    pp1.query("lambda", lambda1);
    pp1.query("Mu_max", Mu1_max);
    pp1.query("Mu_min", Mu1_min);
    pp1.query("N", N1);
    pp1.query("shear modulus", eta1);
    pp1.query("lambda_f",lambda_f);

}

Viscosity::~Viscosity()
{
}

double Viscosity::GetViscosity(const double &S11, const double &S12, const double &S22)
{
    double Gamma_Dot = ComputeGammaDot(S11, S12, S22);
    return Mu_min + (Mu_max - Mu_min) / std::pow((1.0 + (lambda * Gamma_Dot) * (lambda * Gamma_Dot)), N / 2.0);
}

double Viscosity::GetViscosity(const double &Gamma_Dot)
{
    return Mu_min + (Mu_max - Mu_min) / std::pow((1.0 + (lambda * Gamma_Dot) * (lambda * Gamma_Dot)), N / 2.0);
}

double Viscosity::GetViscosity(const double &S11, const double &S12, const double &S22, const double &phi)
{
    double Gamma_Dot = ComputeGammaDot(S11, S12, S22);
    double Mu0 =  Mu_min + (Mu_max - Mu_min) / std::pow((1.0 + (lambda * Gamma_Dot) * (lambda * Gamma_Dot)), N / 2.0);
    double Mu1 =  Mu1_min + (Mu1_max - Mu1_min) / std::pow((1.0 + (lambda1 * Gamma_Dot) * (lambda1 * Gamma_Dot)), N1 / 2.0);
    //return (phi*Mu1 + (1 - phi)*Mu0);
    return (0.5*(1 + phi)*Mu1 + 0.5*(1 - phi)*Mu0);
}

double Viscosity::GetViscosity(const double &S11, const double &S12, const double &S22, const double &phi, const double &dmg)
{
    double Gamma_Dot = ComputeGammaDot(S11, S12, S22);
    double Mu0 =  Mu_min + (Mu_max - Mu_min) / std::pow((1.0 + (lambda * Gamma_Dot) * (lambda * Gamma_Dot)), N / 2.0);
    //double Mu1 =  (1.0 - dmg)*(Mu1_min + (Mu1_max - Mu1_min) / std::pow((1.0 + (lambda1 * Gamma_Dot) * (lambda1 * Gamma_Dot)), N1 / 2.0)) + dmg*Mu0;
    double Mu1 =  Mu1_min + (Mu1_max - Mu1_min) / std::pow((1.0 + (lambda1 * Gamma_Dot) * (lambda1 * Gamma_Dot)), N1 / 2.0);
    //return (phi*Mu1 + (1 - phi)*Mu0);
    return (0.5*(1 + phi)*Mu1 + 0.5*(1 - phi)*Mu0);
}

double Viscosity::GetViscosity(const double &Gamma_Dot, const double &phi)
{
    double Mu0 =  Mu_min + (Mu_max - Mu_min) / std::pow((1.0 + (lambda * Gamma_Dot) * (lambda * Gamma_Dot)), N / 2.0);
    double Mu1 =  Mu1_min + (Mu1_max - Mu1_min) / std::pow((1.0 + (lambda1 * Gamma_Dot) * (lambda1 * Gamma_Dot)), N1 / 2.0);
    //return (phi*Mu1 + (1 - phi)*Mu0);
    return (0.5*(1 + phi)*Mu1 + 0.5*(1 - phi)*Mu0);
}

double Viscosity::GetViscosityFromGammaDot(const double &Gamma_Dot, const double &phi, const double &dmg)
{
    double Mu0 =  Mu_min + (Mu_max - Mu_min) / std::pow((1.0 + (lambda * Gamma_Dot) * (lambda * Gamma_Dot)), N / 2.0);
    double Mu1 =  (1.0 - dmg)*(Mu1_min + (Mu1_max - Mu1_min) / std::pow((1.0 + (lambda1 * Gamma_Dot) * (lambda1 * Gamma_Dot)), N1 / 2.0)) + dmg*Mu0;
    //return (phi*Mu1 + (1 - phi)*Mu0);
    return (0.5*(1 + phi)*Mu1 + 0.5*(1 - phi)*Mu0);
}


double Viscosity::ComputeGammaDot(const double &S11, const double &S12, const double &S22)
{
    return std::sqrt(2.0 * (S11 * S11 + S22 * S22 + 2.0 * S12 * S12 + (S11 + S22) * (S11 + S22)));
}

double Viscosity::ShearModulus(const double &phi)
{
    //return (phi*Mu1 + (1 - phi)*Mu0);
    double Eta0 = eta;
    double Eta1 = eta1;
    return (0.5*(1 + phi)*Eta1 + 0.5*(1 - phi)*Eta0);
}


} /*End namespace mycode */

