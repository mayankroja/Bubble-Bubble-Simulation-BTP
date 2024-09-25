#ifndef _HypreIJ_H_
#define _HypreIJ_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Wed Mar 31 21:50:33 2021


************************************************************************/

#include <LinearSolver.H>
#include "HYPRE.h"
#include "_hypre_utilities.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"
#include "HYPRE_krylov.h"

namespace mycode
{

class HypreIJ : public LinearSolver
{
// public members
public:
	enum solver
	{
		GMRES,
		AMG,
		PCG,
		PCG_AMG
	};

	HypreIJ(
		MPI_Comm,
		amrex::Vector<amrex::Geometry> *geom,
		amrex::Vector<amrex::BoxArray> *grid,
		amrex::Vector<amrex::DistributionMapping> *dmap,
		amrex::Vector<CFMask> *cfmask,
		amrex::Vector<amrex::iMultiFab *> &mask,
		MLTraverseIndex *tri = nullptr);

	HypreIJ(
		MPI_Comm,
		amrex::Vector<amrex::Geometry> *geom,
		amrex::Vector<amrex::BoxArray> *grid,
		amrex::Vector<amrex::DistributionMapping> *dmap,
		amrex::Vector<CFMask> *cfmask,
		MLTraverseIndex *tri = nullptr);

	// constructor:
	// HypreIJ(MPI_Comm, 
	// 		const amrex::Geometry *geom,
	// 		const amrex::BoxArray *grid,
	// 		const amrex::DistributionMapping *dmap, 
	// 		amrex::iMultiFab *mask = nullptr, 
	// 		MLTraverseIndex *tri = nullptr);

	// destructor:
	~HypreIJ();

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

	void AMGSolver();
	void GMRESSolver();
	void PCGSolver(bool preconditioner = false);

	void init();

	/************************************************************************/

	/****************************   data    *********************************/

	HYPRE_IJMatrix HYPREMat_A = NULL;
	HYPRE_IJVector HYPREVec_b = NULL;
	HYPRE_IJVector HYPREVec_x = NULL;
	HYPRE_Solver HYPREsolver = NULL;
	HYPRE_Solver HYPREprecond = NULL;

	HYPRE_ParCSRMatrix par_A = NULL;
	HYPRE_ParVector par_b = NULL;
	HYPRE_ParVector par_x = NULL;
	solver ls_solver_;
	
	/************************************************************************/

}; /*End class HypreIJ */

} /*End namespace mycode */

#endif
