#ifndef _Ellipse_H_
#define _Ellipse_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Wed Dec 30 19:13:41 2020


************************************************************************/

#include <Interface.H>

namespace mycode
{

class Ellipse : public Interface
{
// public members
public:

// constructor:
	Ellipse(const amrexMesh &mesh, const std::string& name);

	// destructor:
	~Ellipse();

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
	amrex::Real Lp, L0;
	amrex::Vector<amrex::Real> vel;
	amrex::Vector<amrex::Real> radii;

	/************************************************************************/

}; /*End class Ellipse */

std::unique_ptr<Ellipse> MakeEllipse(const amrexMesh &, const std::string &);

} /*End namespace mycode */

#endif
