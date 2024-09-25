#ifndef _Cylinder_H_
#define _Cylinder_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Wed Dec 30 19:13:41 2020


************************************************************************/

#include <Interface.H>

namespace mycode
{

class Cylinder : public Interface
{
// public members
public:

// constructor:
	Cylinder(const amrexMesh &mesh, const std::string& name);

	// destructor:
	~Cylinder();

/**************************** functions *********************************/

	void PrescribeLevelSetMotion(const amrex::Real& t);

	void ComputeLevelSetMotion(const amrex::Real &t, const amrex::Real &dt, int RKStage);

	amrex::RealArray getOrigin();

	//! functions for advected levelset
	void TubeIdentification();
	void TubeAdvectLevelSet(const amrex::MultiFab &u, const amrex::MultiFab &v, const amrex::Real &dt);
        void TubeAdvectLevelSet_RK2(const amrex::MultiFab &u, const amrex::MultiFab &v, const amrex::Real &dt){}
        void TubeAdvectLevelSet_RK3(int RKStage, const amrex::MultiFab &u, const amrex::MultiFab &v, const amrex::Real &dt){}
	//void Regularization();
	/************************************************************************/

	// private members
private:
	/**************************** functions *********************************/

	void ComputePsi();

	/************************************************************************/

	/****************************   data    *********************************/
	amrex::Real Lp, L0, radius;
	amrex::Vector<amrex::Real> vel;

/************************************************************************/

}; /*End class Cylinder */

std::unique_ptr<Cylinder> MakeCylinder(const amrexMesh &, const std::string &);

} /*End namespace mycode */

#endif
