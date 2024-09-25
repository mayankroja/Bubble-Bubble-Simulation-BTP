#ifndef _amrexMesh_H_
#define _amrexMesh_H_

/************************************************************************

Author: Mahesh K Sawardekar
Date: Fri Dec 27 17:33:57 2019

    * wrapper for AMReX Geometry, BoxArray, distribution mapping..

************************************************************************/

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>

namespace mycode
{

/**
 * @brief wrapper class for AMReX Geometry, BoxArray, distribution mapping..
 * 
 */
class amrexMesh
{
// public members
public:

// default constructor:
	amrexMesh();

	amrexMesh(const amrex::BoxArray& grid,
			  const amrex::Geometry& geom,
			  const amrex::DistributionMapping& dmap);

	// destructor:
	~amrexMesh();

/**************************** functions *********************************/

	void SetBoxArray(const amrex::BoxArray &ba_in) { grid_ = ba_in; }
	void SetDistributionMap(const amrex::DistributionMapping &dmap_in) { dmap_ = dmap_in; }
	void ReadFromCheckFile(std::string &file);

	/**
	 * @brief get boxArray reference
	 * 
	 * @return const amrex::BoxArray& 
	 */
	const amrex::BoxArray& grid() { return grid_; }
	const amrex::BoxArray& grid() const { return grid_; }


	/**
	 * @brief get geometry reference
	 * 
	 * @return const amrex::Geometry& 
	 */
	const amrex::Geometry& geometry() { return geom_; }
	const amrex::Geometry& geometry() const { return geom_; }


	/**
	 * @brief get DistributionMapping reference
	 * 
	 * @return const amrex::DistributionMapping& 
	 */
	const amrex::DistributionMapping& dmap() { return dmap_; }
	const amrex::DistributionMapping& dmap() const { return dmap_; }

/************************************************************************/

// private members
private:

/****************************   data    *********************************/
	amrex::BoxArray grid_;
	amrex::Geometry geom_;
	amrex::DistributionMapping dmap_;

/************************************************************************/

}; /*End class amrexMesh */

} /*End namespace mycode */

#endif
