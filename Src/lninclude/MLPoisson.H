#ifndef _MLPoisson_H_
#define _MLPoisson_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Fri Jan 22 10:58:21 2021

	* Multi-Level poisson solver
	* If the level is not coasest level then the soln_
	* should be initialised by interpolating from coarse data
	* including ghost cells at c/f bndry

************************************************************************/

// #include <amrexMesh.H>
#include <AMReX_MultiFab.H>
#include <AMReX_LO_BCTYPES.H>

#include <Stencil.H>
#include <Interface.H>
#include <InterceptData.H>
#include <defs.H>
#include <LinearSolver.H>

namespace mycode
{

class MLPoisson
{
// public members
public:

	// constructor:
	MLPoisson(MPI_Comm comm,
			  amrex::Vector<amrex::Geometry> *geom,
			  amrex::Vector<amrex::BoxArray> *grid,
			  amrex::Vector<amrex::DistributionMapping> *dmap,
			  amrex::Vector<std::vector<std::unique_ptr<Interface>>>* IF = nullptr);

	// // constructor:
	// MLPoisson(amrexMesh* mesh,
	// 		const std::vector<std::unique_ptr<Interface>> &);

	// destructor:
	~MLPoisson();

	/**************************** functions *********************************/
	
	/// set boundary conditions on outer domain
	void setBCs(const BCArray &lbc, const BCArray &hbc)
	{
		lobc = lbc;
		hibc = hbc;
	}

	/// set boundary condition on EB
	void setIntBCs(const std::vector<amrex::LinOpBCType> &bc) { intbc = bc; }

	/// set function to define interface boundary values
	void setIntBoundaryFunc(IntBoundaryFunc f) { int_boundary_fun = f; }

	void setBoundaryFunc(BoundaryFunc f);
	void setBoundaryFunc(int i, BoundaryFunc f);

	void setRhs(amrex::Vector<amrex::MultiFab*> mf) { rhs_ = mf; }

	void setSoln(amrex::Vector<amrex::MultiFab*> mf) { soln_ = mf; }

	void setMask(amrex::Vector<amrex::iMultiFab*> mask) { mask_ = mask; }

	void solve(amrex::Real Time);

	void fillBoundaryFine();

	int getSystemSize() { return global_sys_size; }

	int getSystemSize(int lev) { return sys_size_lev[lev]; }
	
	amrex::Real getResidue() { return residue_; }
	
	int getNIterations() { return num_iterations; }

	/************************************************************************/

	// private members
private:
	/**************************** functions *********************************/

	void setupSystem();
	void addGhostContro
	(
                const amrex::Box &domain,
		const amrex::IntVect &gcell,			/// the cut cell
		Stencil &st,							/// st to which wts need to be added
		amrex::Array4<int const> const &cellid, /// the traverse index data
		const amrex::Real &fact,				/// scale factor
		const InterceptData &icpt,				/// intercept data @ cut cell (to get computed wts)
		amrex::Array4<int const> const &mask	/// PMask
	);

	/**
	 *  get the wts for coarse cells using quadratic interpolation
	 *  when on fine level and nbr is on coarse level
	 */
	void getCoarseWts_fine
	(
		int cross_dir,
		const amrex::IntVect &gcell, // fine ghost cell
		const amrex::BoxArray &fgrid,
		const amrex::Box& fdomain,
		amrex::Array<amrex::Real, 3> &wts,
		amrex::Array<amrex::IntVect, 3> &cells
	);

	/**
	 *  get the wts for coarse cells using quadratic interpolation
	 *  when on coarse level and nbr is on fine level
	 */
	void getCoarseWts_coarse
	(
		int cross_dir,
		const amrex::IntVect &pcell,   // P cell
		const amrex::IntVect &pcell_f, // P cell fine
		amrex::Array4<int const> const &cfmask,
		amrex::Array<amrex::Real, 3> &wts,
		amrex::Array<amrex::IntVect, 3> &cells
	);

	/**
	 *  get the coefficients when on fine level and nbr is on corser level  
	 * 
	 */
	mycode::Stencil getGhostStencil_fine
	(
		int cross_dir,
		int shift_dir,
		int side,                       // lo or hi (i.e. 0 => cell P is at lo face, 1 => at hi face) 
		const amrex::IntVect &gcell,    // fine ghost cell
		const amrex::BoxArray& fgrid,
		const amrex::Box& fdomain,      // to check if at domain bndry
		amrex::Array4<int const> const& cell_id,
		amrex::Array4<int const> const& cell_id_crse
	);

	/**
	 *  get the coefficients when nbr is on finer level  
	 * 
	 */
	mycode::Stencil getGhostStencil_coarse
	(
		int cross_dir,
		int shift_dir,
		int side,                       // lo or hi (i.e. 0 => cell P is at lo face, 1 => at hi face) 
		const amrex::IntVect &pcell,    // P cell on coarse
		const amrex::IntVect &pcell_f,    // P cell fine
		amrex::Array4<int const> const &cfmask,
		amrex::Array4<int const> const& cell_id,
		amrex::Array4<int const> const& cell_id_fine
	);

	/************************************************************************/

	/****************************   data    *********************************/
	MPI_Comm comm_;
	amrex::Vector<amrex::Geometry> *geom_;
	amrex::Vector<amrex::BoxArray> *grid_;
	amrex::Vector<amrex::DistributionMapping> *dmap_;

	/// rhs which is set by application code
	amrex::Vector<amrex::MultiFab*> rhs_;

	/// solution initialised by application code
	/// after solving Ax = b, x values are stored in soln_
	amrex::Vector<amrex::MultiFab*> soln_;

	/// mask to identify regular, ghost and covered cells
	amrex::Vector<amrex::iMultiFab*> mask_;

	/// Collection of interfaces in application
	amrex::Vector<std::vector<std::unique_ptr<Interface>>> *interfaces;

	BCArray lobc;
	BCArray hibc;
	std::vector<amrex::LinOpBCType> intbc;

	IntBoundaryFunc int_boundary_fun;
	amrex::Array<BoundaryFunc, 2*AMREX_SPACEDIM> boundary_funcs;

	/// soln parms
	amrex::Real abs_tol_;
	amrex::Real rel_tol_;
	int maxIter_;

	const int stencil_ = 4;
	bool singular_ = false;
	int global_sys_size;
	amrex::Vector<int> sys_size_lev;
	amrex::Real residue_;
	int num_iterations;

	bool isAxisymmetric = false;

	std::unique_ptr<LinearSolver> solver;
	LinearSolver::linear_solver_lib solver_lib;
	/************************************************************************/

}; /*End class MLPoisson */

} /*End namespace mycode */

#endif
