#ifndef _CFStencil_H_
#define _CFStencil_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Wed Apr 21 21:26:04 2021

	* compute wts of ghost cell @ lev/lev-1 interface
	* using quadratic interpolation

************************************************************************/

#include <LagrangeInterpolation.H>
#include <AMReX_IntVect.H>
#include <utility>

namespace mycode
{

class CFStencil
{
// public members
public:

	enum stencil_type
	{
		shifted_minus,
		shifted_plus,
		regular
	};

	// constructor:
	CFStencil(stencil_type type);

	// destructor:
	~CFStencil();

	/**************************** functions *********************************/

	amrex::Vector<amrex::Real> &coarseStencilPlus() { return coarse_wts_p; }
	const amrex::Vector<amrex::Real> &coarseStencilPlus() const { return coarse_wts_p; }
	
	amrex::Vector<amrex::Real> &coarseStencilMinus() { return coarse_wts_m; }
	const amrex::Vector<amrex::Real> &coarseStencilMinus() const { return coarse_wts_m; }

	/************************************************************************/

	// private members
private:
/**************************** functions *********************************/


/************************************************************************/


/****************************   data    *********************************/

	stencil_type type_;
	amrex::Vector<amrex::Real> coarse_wts_p; // evaluated at + side
	amrex::Vector<amrex::Real> coarse_wts_m; // evaluated at - side

	//! fine wts are fixed for ref_ratio 2 : 1/15 * (8, 10, -3)
	// amrex::Array<amrex::Real, 3> fine_wts;     // 

	/************************************************************************/

}; /*End class CFStencil */

} /*End namespace mycode */

#endif
