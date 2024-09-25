#include <Viscosity.H>

namespace mycode
{

    double Viscosity::GetConductivity(double Mu, double Prandtl_no)
    {
        return Mu*(1.0/Prandtl_no);
    }  
} // namespace mycode
