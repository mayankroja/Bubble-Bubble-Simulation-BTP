#ifndef _LinearSolver_H_
#define _LinearSolver_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Wed Mar 31 19:27:22 2021

	Interface for hypre and petsc linear solvers.
	The matrix A , vectors x and b should be filled from application code,
	one row at a time.

************************************************************************/

// #include <amrexMesh.H>
#include <MLTraverseIndex.H>
#include <AMReX_MultiFab.H>
#include <CFMask.H>

namespace mycode
{

class LinearSolver
{
// public members
public:
	enum linear_solver_lib
	{
		HYPRE,
		PETSC
	};

	// destructor:
	virtual ~LinearSolver();

/**************************** functions *********************************/
	void define
	(
		MPI_Comm,
		amrex::Vector<amrex::Geometry> *geom,
		amrex::Vector<amrex::BoxArray> *grid,
		amrex::Vector<amrex::DistributionMapping> *dmap,
		amrex::Vector<CFMask> *cfmask,
		amrex::Vector<amrex::iMultiFab *> &mask,
		MLTraverseIndex *tri = nullptr
	);


	void define(
		MPI_Comm,
		amrex::Vector<amrex::Geometry> *geom,
		amrex::Vector<amrex::BoxArray> *grid,
		amrex::Vector<amrex::DistributionMapping> *dmap,
		amrex::Vector<CFMask> *cfmask,
		MLTraverseIndex *tri = nullptr);

	void solve(amrex::Vector<amrex::MultiFab*> &soln);

	void setVerbose(int v) { verbose = v; }
	void printSystem(bool b) { print_system = b; }
	void singularSystem(bool b) { is_singular_ = b; }

	//! set 1 row at a time
	virtual void set_A_Val(int *ncols, int *row_id, int *col_id, amrex::Real *values) = 0;
	virtual void set_x_Val(int *row_id, amrex::Real *values) = 0;
	virtual void set_b_Val(int *row_id, amrex::Real *values) = 0;
	virtual void assembleSystem() = 0;

	int getSystemSize()
	{
		return tridx->getWorldSize();
	}
	
	int getNIterations() { return num_iterations; }
	amrex::Real getResidue() { return res; }

	/************************************************************************/

	// protected members
protected:
	/**************************** functions *********************************/

	void getSolution(amrex::Vector<amrex::MultiFab*> &soln);

	virtual void getVecVal(void *x, int nrows, int *row_id, amrex::Real *values) = 0;
	virtual void preAssembleSystem() = 0;
	virtual void solverSetupAndSolve() = 0;

	/************************************************************************/

	/****************************   data    *********************************/
	void *Mat_A = nullptr;
	void *Vec_x = nullptr;
	void *Vec_b = nullptr;

	amrex::Vector<amrex::Geometry> *geom_;
	amrex::Vector<amrex::BoxArray> *grids_;
	amrex::Vector<amrex::DistributionMapping> *dmap_;

	amrex::Vector< amrex::iMultiFab*> MaskPtr_;
	MLTraverseIndex* tridx;
	amrex::Vector<CFMask> *cfmask_;

	int nlevels;

	int max_iter = 100;
	double rel_tolerance = 1e-10;
	double abs_tolerance = 0.0;
	amrex::Real res;
	int num_iterations;

	MPI_Comm comm_ = MPI_COMM_NULL;
	int verbose = 0;
	bool print_system = false;
	linear_solver_lib ls_lib_;
	bool is_singular_ = false;

	/************************************************************************/

}; /*End class LinearSolver */

} /*End namespace mycode */

#endif
