#ifndef _InterfaceAdvect_H_
#define _InterfaceAdvect_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Wed Dec 30 19:13:41 2020


************************************************************************/

#include <Interface.H>

namespace mycode
{

class InterfaceAdvect : public Interface
{
// public members
public:

// constructor:
	InterfaceAdvect(const amrexMesh &mesh, const std::string& name);

	// destructor:
	~InterfaceAdvect();

/**************************** functions *********************************/

	void PrescribeLevelSetMotion(const amrex::Real& t) = 0;

	void ComputeLevelSetMotion(const amrex::Real &t, const amrex::Real &dt, int RKStage) = 0;

	amrex::RealArray getOrigin() = 0;

	//! additional functions for advection of levelset
	void TubeIdentification();
	void TubeAdvectLevelSet(const amrex::MultiFab &u, const amrex::MultiFab &v, const amrex::Real &dt);
        void TubeAdvectLevelSet_RK2(const amrex::MultiFab &u, const amrex::MultiFab &v, const amrex::Real &dt);
        void TubeAdvectLevelSet_RK3(int RKStage, const amrex::MultiFab &u, const amrex::MultiFab &v, const amrex::Real &dt);

	/************************************************************************/
	// private members
private:
	/**************************** functions *********************************/

	void ComputePsi();
	bool nbrTTube(amrex::Array4<amrex::Real> const &psi, const amrex::Box &bx);

	/************************************************************************/

}; /*End class InterfaceAdvect */


} /*End namespace mycode */

#endif
