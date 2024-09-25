#ifndef _Lagrange_Interpolation_H_
#define _Lagrange_Interpolation_H_

#include<AMReX_Vector.H>

namespace mycode
{
    amrex::Vector<amrex::Real> LagrangeInterploation(const amrex::Vector<amrex::Real> &data, const amrex::Real &Xp);

} // namespace mycode

#endif