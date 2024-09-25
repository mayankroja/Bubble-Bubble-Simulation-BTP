#ifndef _PetscIJ_H_
#define _PetscIJ_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Wed Mar 31 21:50:33 2021


************************************************************************/

#include <LinearSolver.H>
#include <petscksp.h>

namespace mycode
{

class PetscIJ : public LinearSolver
{
// public members
public:
	enum solver
	{
		GMRES,
		CG,
		BICGS
	};

	PetscIJ(
		MPI_Comm,
		amrex::Vector<amrex::Geometry> *geom,
		amrex::Vector<amrex::BoxArray> *grid,
		amrex::Vector<amrex::DistributionMapping> *dmap,
		amrex::Vector<CFMask> *cfmask,
		amrex::Vector<amrex::iMultiFab *> &mask,
		MLTraverseIndex *tri = nullptr);

	PetscIJ(
		MPI_Comm,
		amrex::Vector<amrex::Geometry> *geom,
		amrex::Vector<amrex::BoxArray> *grid,
		amrex::Vector<amrex::DistributionMapping> *dmap,
		amrex::Vector<CFMask> *cfmask,
		MLTraverseIndex *tri = nullptr);

	// destructor:
	~PetscIJ();

	/**************************** functions *********************************/

	//! set 1 row at a time
	virtual void set_A_Val(int *ncols, int *row_id, int *col_id, amrex::Real *values);
	virtual void set_x_Val(int *row_id, amrex::Real *values);
	virtual void set_b_Val(int *row_id, amrex::Real *values);
	virtual void assembleSystem();

	void setSolver(solver sid) { ls_solver_ = sid; }

	/************************************************************************/

	// protected members
protected:
/**************************** functions *********************************/
	virtual void getVecVal(void *x, int nrows, int *row_id, amrex::Real *values);
	virtual void preAssembleSystem();
	virtual void solverSetupAndSolve();

	void init();

	/************************************************************************/

	/****************************   data    *********************************/

	int regular_stencil_size = 2 * AMREX_SPACEDIM + 1;
	KSP PETSCsolver = nullptr;
	Mat PETSCMat_A = nullptr;
	Vec PETSCVec_b = nullptr;
	Vec PETSCVec_x = nullptr;
	solver ls_solver_;

	/************************************************************************/

}; /*End class PetscIJ */

} /*End namespace mycode */

#endif
