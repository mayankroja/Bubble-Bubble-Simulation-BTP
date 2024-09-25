#ifndef _InterceptData_H_
#define _InterceptData_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Fri Jan  1 10:32:22 2021


************************************************************************/

#include <AMReX_IntVect.H>
#include <AMReX_Array.H>
#include <utility>

namespace mycode
{

struct InterceptData
{
	/// the cell index of cut cell
	amrex::IntVect cellid_;
	std::vector<amrex::IntVect> nbr_cellid_; 

        /// Flag decides if the levleset in this cell need regularization
	int r_flag = 0;//default 0 = no, 1 = yes
	//int s_flag = 0;//default 0 = no, 1 = yes; for local smoothening

	/// how many faces of this cell are cut
	int n_intercepts, n_nbr;

        /// The size of stencil for constructing least-square field for this cell
        int stencil_ = 4;
        /// Order of least-sq approximation in that cell
        int PORDER = 4;

	/// dir of cut face 0 => x, 1 => y
	amrex::Array<int, 2 * AMREX_SPACEDIM> type_;

	/// data related to location of cut and for normal computation
	amrex::Array<amrex::Real, 2 * AMREX_SPACEDIM> frac_;
	amrex::Array<amrex::Real, 2 * AMREX_SPACEDIM> psix_;
	amrex::Array<amrex::Real, 2 * AMREX_SPACEDIM> psiy_;
	amrex::Array<amrex::Real, 2 * AMREX_SPACEDIM> kappa_;
	amrex::Array<amrex::Real, 2 * AMREX_SPACEDIM> norm_shear_;
	amrex::Array<amrex::Real, 2 * AMREX_SPACEDIM> tan_shear_;
        amrex::Array<amrex::Real, 2 * AMREX_SPACEDIM> Gamma_dot_Int;


	/// values of field variables on interface
	/// P, u, v @ interface
	amrex::Array<amrex::Real, 2 * AMREX_SPACEDIM> P, u, v, P_comp, P_Int, T_Int, U_comp, V_comp, X_Visc_Int, Y_Visc_Int,F11_Int,F12_Int,F21_Int,F22_Int,F33_Int, Dmg_Int, phi_Int,emax_Int;

	/// 
	std::vector<amrex::Real> Interpolation_weights;

	
	/// default constructor
        /// default constructor
        InterceptData()
        {
	    if(AMREX_SPACEDIM == 2)
            {
	        Dmg_Int = {0.0,0.0};
		phi_Int = {0.0,0.0};
		emax_Int = {0.0,0.0};
	    }
            if(AMREX_SPACEDIM == 3)
            {
                Dmg_Int = {0.0,0.0,0.0};
                phi_Int = {0.0,0.0,0.0};
		emax_Int = {0.0,0.0,0.0};
            }
	}
	/*
	InterceptData()
	:
	local_index(0),
	cellid_({0,0}),
	n_intercepts(0),
	type_({0}),
	frac_({0.0}),
	psix_({0.0}),
	psiy_({0.0}),
	kappa_({0.0}),
	P(0.0),
	u(0.0),
	v(0.0)
	{}
	*/
        /*
	/// copy constructor
	InterceptData(const InterceptData& id)
	:
	local_index(id.local_index),
	cellid_(id.cellid_),
	n_intercepts(id.n_intercepts),
	type_(id.type_),
	frac_(id.frac_),
	psix_(id.psix_),
	psiy_(id.psiy_),
	kappa_(id.kappa_),
	P(id.P),
	u(id.u),
	v(id.v)
	{}

	/// move constructor
	InterceptData(InterceptData&& id)
	:
	local_index(std::exchange(id.local_index, 0)),
	cellid_(std::move(id.cellid_)),
	n_intercepts(std::exchange(id.n_intercepts, 0)),
	type_(std::move(id.type_)),
	frac_(std::move(id.frac_)),
	psix_(std::move(id.psix_)),
	psiy_(std::move(id.psiy_)),
	kappa_(std::move(id.kappa_)),
	P(std::exchange(id.P, 0.0)),
	u(std::exchange(id.u, 0.0)),
	v(std::exchange(id.v, 0.0))
	{}
	*/
};

} /*End namespace mycode */

#endif
