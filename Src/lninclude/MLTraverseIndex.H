#ifndef _MLTraverseIndex_H_
#define _MLTraverseIndex_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Wed Nov  4 12:15:51 2020

	Determine the global sequence of unknowns in x of the system Ax = b 

************************************************************************/

#include <AMReX_iMultiFab.H>
#include <AMReX_LO_BCTYPES.H>

namespace mycode
{

class MLTraverseIndex
{
// public members
public:

// default constructor:
	MLTraverseIndex();

// constructor:
	MLTraverseIndex
	(
		const amrex::Vector<amrex::Geometry>& geom,
		const amrex::Vector<amrex::BoxArray>& grid,
	    const amrex::Vector<amrex::DistributionMapping>& dmap,
	    amrex::Vector<amrex::iMultiFab*> mask,
		int nghost = 1
	);

	MLTraverseIndex(
		const amrex::Vector<amrex::Geometry> &geom,
		const amrex::Vector<amrex::BoxArray> &grid,
		const amrex::Vector<amrex::DistributionMapping> &dmap,
		int nghost = 1);

	// destructor:
	~MLTraverseIndex();

/**************************** functions *********************************/
	/**
	 * @brief compute cell_id and offset for cc data
	 * 
	 * @param geom 
	 * @param grid 
	 * @param dmap
	 * @param mask : cc indicator function
	 */
	void define
	(
		const amrex::Vector<amrex::Geometry>& geom,
		const amrex::Vector<amrex::BoxArray>& grid,
	    const amrex::Vector<amrex::DistributionMapping>& dmap,
	    amrex::Vector<amrex::iMultiFab*> mask,
		int nghost = 1
	);

	void define(
		const amrex::Vector<amrex::Geometry> &geom,
		const amrex::Vector<amrex::BoxArray> &grid,
		const amrex::Vector<amrex::DistributionMapping> &dmap,
		int nghost = 1);

	/**
	 * @brief Get the Offset object
	 * 
	 * @return amrex::LayoutData<int>& 
	 */
	auto& getOffset() { return offset_; }
	const auto& getOffset() const { return offset_; }

	/**
	 * @brief Get the Traverse Index object
	 * 
	 * @return amrex::iMultiFab& 
	 */
	auto& getTrIndex() { return cell_id_; }
	const auto& getTrIndex () const { return cell_id_; }

	auto &getTrIndex(int lev) { return cell_id_[lev]; }
	const auto &getTrIndex(int lev) const { return cell_id_[lev]; }

	/**
	 * @brief Get the traverse index on level 'to_lev' from level 'from_lev'
	 *        The TraverseIndex on 'from_lev' may have different DistributionMap
	 *        So to have the info abt nbr cell row number while forming System A
	 *        We need the copy on 'to_lev' 
	 * 
	 * 		  It is assumed that the levels are adjcent, i.e lev/lev-1 or lev/lev+1
	 *        also refinement ratio is assumed to be 2
	 * 
	 */

	std::unique_ptr<amrex::iMultiFab> getTrIndex(int to_lev, int from_lev);

	const int getProcBegin() { return proc_begin_; }
	const int getProcBegin() const { return proc_begin_; }

	const int getProcEnd() { return proc_end_; }
	const int getProcEnd() const { return proc_end_; }

	const int getWorldSize() { return nNodes_world_; } 
	const int getWorldSize() const { return nNodes_world_; }

	const int getLevSize(int lev) { return nNodes_lev[lev]; }
	const int getLevSize(int lev) const { return nNodes_lev[lev]; }

	const amrex::Vector<int>& getLevSize() { return nNodes_lev; }
	const amrex::Vector<int>& getLevSize() const { return nNodes_lev; }

	/************************************************************************/

	// private members
private:

/****************************   data    *********************************/
	amrex::Vector<amrex::LayoutData<int>> offset_;
	int proc_begin_;
	int proc_end_;
	
	// required by petsc
	int nNodes_world_;
	int global_begin_;

	amrex::Vector<int> nNodes_lev;

	amrex::Vector<amrex::iMultiFab> cell_id_;
/************************************************************************/

}; /*End class MLTraverseIndex */


} /*End namespace mycode */

#endif
