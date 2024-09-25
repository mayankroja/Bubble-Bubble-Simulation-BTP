#ifndef _CFMask_H_
#define _CFMask_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Tue Apr 20 11:59:13 2021

	* Determine if the cells are covered by fine cell

************************************************************************/

#include <AMReX_iMultiFab.H>

namespace mycode
{

class CFMask
{
// public members
public:

	enum MaskVal
	{
		covered,
		valid,
		cfbndry,
		out_of_grid
	};

	// default constructor:
	CFMask();

	// constructor:
	CFMask( int lev,
			const amrex::Vector<amrex::Geometry>& geom,
			const amrex::Vector<amrex::BoxArray>& grid,
			const amrex::Vector<amrex::DistributionMapping>& dmap,
			int nghost = 1);

// destructor:
	~CFMask();

/**************************** functions *********************************/
	// constructor:
	void define(int lev,
				const amrex::Vector<amrex::Geometry> &geom,
				const amrex::Vector<amrex::BoxArray> &grid,
				const amrex::Vector<amrex::DistributionMapping> &dmap,
				int nghost = 1);

	auto &Mask() { return mask_; }
	const auto &Mask() const { return mask_; }

	/************************************************************************/

// private members
private:

	/****************************   data    *********************************/
	amrex::iMultiFab mask_;

	/************************************************************************/

}; /*End class CFMask */

} /*End namespace mycode */

#endif
