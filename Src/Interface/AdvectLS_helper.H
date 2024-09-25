#ifndef _Adevct_LS_Helper_
#define _Adevct_LS_Helper_

#include <AMReX_Array4.H>

namespace mycode
{

void WENO5_LS
(
    double &Psix_L, double &Psix_R,
    double &Psiy_L, double &Psiy_R,
    int i, int j,
    const amrex::Real *dx,
    amrex::Array4<amrex::Real const> const &Psi
);

void WENO5_Fij
(
    double &Psix_L, double &Psix_R,
    double &Psiy_L, double &Psiy_R,
    int i, int j,
    const amrex::Real *dx,
    amrex::Array4<amrex::Real const> const &Psi
);

} // namespace mycode

#endif
