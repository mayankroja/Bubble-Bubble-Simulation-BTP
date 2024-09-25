#ifndef _Bubble_H_
#define _Bubble_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Wed Dec 30 19:13:41 2020


************************************************************************/

#include <InterfaceAdvect.H>
#include <AMReX_iMultiFab.H>

namespace mycode
{

class Bubble : public InterfaceAdvect
{
// public members
public:

// constructor:
	Bubble(const amrexMesh &mesh, const std::string& name);

	// destructor:
	~Bubble();

/**************************** functions *********************************/

	void PrescribeLevelSetMotion(const amrex::Real& t);

	void ComputeLevelSetMotion(const amrex::Real &t, const amrex::Real &dt, int RKStage);

	amrex::RealArray getOrigin();

	/************************************************************************/

	// private members
private:
	/**************************** functions *********************************/
	void ComputePsi();

	/************************************************************************/

	/****************************   data    *********************************/

	amrex::Real radius;

	/************************************************************************/

}; /*End class Bubble */

std::unique_ptr<Bubble> MakeBubble(const amrexMesh &, const std::string &);

} /*End namespace mycode */

#endif
