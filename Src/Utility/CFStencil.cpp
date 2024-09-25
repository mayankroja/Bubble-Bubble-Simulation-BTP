#include "CFStencil.H"

namespace mycode
{

CFStencil::CFStencil(stencil_type type)
:
type_(type)
{
    amrex::Vector<amrex::Real> X(3);
    if (type_ == shifted_minus)
    {
        X[0] = 0.0;
        X[1] = -1.0;
        X[2] = -2.0;
    }
    else if (type_ == shifted_plus)
    {
        X[0] = 0.0;
        X[1] = 1.0;
        X[2] = 2.0;
    }
    else
    {
        X[0] = -1.0;
        X[1] = 0.0;
        X[2] = 1.0;
    }
    
    coarse_wts_p = LagrangeInterploation(X, 0.25);
    coarse_wts_m = LagrangeInterploation(X, -0.25);

    //! scale wts by 8/15 which is wt corresponding to lev pts

    for (size_t i = 0; i < 3; i++)
    {
        coarse_wts_p[i] *= 8.0 / 15.0;
        coarse_wts_m[i] *= 8.0 / 15.0;
    }
}

CFStencil::~CFStencil()
{
}

} /*End namespace mycode */

