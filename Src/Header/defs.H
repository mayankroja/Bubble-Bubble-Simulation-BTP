#ifndef _DEFS_H_
#define _DEFS_H_

#include <AMReX_Array.H>
#include <AMReX_LO_BCTYPES.H>
#include <functional>

namespace mycode
{

    /// function to set boundary values at interface
    /// functions of space, normal and time
    typedef std::function<amrex::Real(AMREX_D_DECL(amrex::Real, amrex::Real, amrex::Real), amrex::RealArray, amrex::Real)> IntBoundaryFunc;
    
    /// function to set outer boundary values
    /// function of space and time
    typedef std::function<amrex::Real(AMREX_D_DECL(amrex::Real, amrex::Real, amrex::Real), amrex::Real)> BoundaryFunc;
    
    /// function to set internal field values
    /// function of space and time
    typedef std::function<amrex::Real(AMREX_D_DECL(amrex::Real, amrex::Real, amrex::Real), amrex::Real)> FieldFunc;
    
    /// boundary condtion array in dir
    typedef amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> BCArray;

} // namespace mycode



#endif