#ifndef _Foil_H_
#define _Foil_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Wed Dec 30 19:13:41 2020


************************************************************************/

#include <Interface.H>

namespace mycode
{

class Foil : public Interface
{
// public members
public:

// constructor:
	Foil(const amrexMesh &, const std::string &);

	// destructor:
	~Foil();

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
	amrex::Real N_Shape_1, N_Shape_2;
	amrex::Real Lp, L0, Len_foil, PreFactor;
	amrex::Real Phase;
	amrex::Real rhoArea;

/************************************************************************/

}; /*End class Foil */

std::unique_ptr<Foil> MakeFoil(const amrexMesh &, const std::string &);

} /*End namespace mycode */

#endif
