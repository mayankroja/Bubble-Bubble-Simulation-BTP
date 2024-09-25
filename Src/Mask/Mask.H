#ifndef _Mask_H_
#define _Mask_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Tue Mar 16 10:29:07 2021

	* class to identify ghost cells and to compute the interfaceData qty.

************************************************************************/

#include <AMReX_iMultiFab.H>
#include <Interface.H>
#include <InterceptData.H>
#include <AMReX_LO_BCTYPES.H>
#include <Viscosity.H>

namespace mycode
{

class Mask
{
// public members
public:
	Mask() { }

	// constructor:
	Mask(amrexMesh*, std::vector<std::unique_ptr<Interface>>*);

	// destructor:
	~Mask();

	/**************************** functions *********************************/

	void define(amrexMesh*, std::vector<std::unique_ptr<Interface>>*);
	void GhostCellIdentfication();
        void SetInterfacePressureVel();
	void ComputeDrag();
	void FillInGhost(amrex::MultiFab&, const std::vector<amrex::LinOpBCType> &intbcs);
        void FillInGhostTheta(amrex::MultiFab&, const std::vector<amrex::LinOpBCType> &intbcs);
        void FillInGhostPhaseField(amrex::MultiFab&, const std::vector<amrex::LinOpBCType> &intbcs);
	void ExtendPhaseFieldLSQ(amrex::MultiFab&);
        void FillInVelocityComponents(amrex::MultiFab& U);
	void FillInRefConfig(amrex::MultiFab &X, amrex::MultiFab &Y, amrex::MultiFab &DMG);
	void FillInGhostFij(int &iscalar, amrex::MultiFab&);
	void FillInGhostFij(amrex::MultiFab&, amrex::MultiFab&, amrex::MultiFab&, amrex::MultiFab&, amrex::MultiFab&);
	void FillInGhostFij(amrex::MultiFab&, amrex::MultiFab&, amrex::MultiFab&, amrex::MultiFab&, amrex::MultiFab&, amrex::MultiFab&);
        void ExtendVelocityLSQ(amrex::MultiFab& U);
	void ExtendFijLSQ(amrex::MultiFab &mfF11, amrex::MultiFab &mfF12, amrex::MultiFab &mfF21, amrex::MultiFab &mfF22, amrex::MultiFab &mfF33);
	void getFaceMask(amrex::iMultiFab &fmask, int dir);

	amrex::iMultiFab &getMask() { return mask_; }
	const amrex::iMultiFab &getMask() const { return mask_; }

        amrex::iMultiFab &getPMask() { return PMask; }
        const amrex::iMultiFab &getPMask() const { return PMask; }


        amrex::iMultiFab &getUMask() { return UMask; }
        const amrex::iMultiFab &getUMask() const { return UMask; }

        amrex::iMultiFab &getVMask() { return VMask; }
        const amrex::iMultiFab &getVMask() const { return VMask; }

	amrex::LayoutData< std::array< std::vector<amrex::IntVect>, 12 > >& getIndexAdditionalLayers()
	{
		return Index_Additional_Layers_LD;
	}

	/************************************************************************/
        Viscosity visc_;
        amrex::Real Mu;
// private members
private :
	/**************************** functions *********************************/
        bool isAxisymmetric = false;
        int nghost;

	void MaskBC();

	int Cubic_Solve(double &x_, double a_, double b_, double c_, double d_);
	void Compute_Normal_Curvature
	(
		amrex::Real *dPsi,
		amrex::Real &Curv,
		int i, int j,
		amrex::Array4<amrex::Real const> const &Psi,
		const amrex::Real &deltax
	);
        void Compute_Normal_Curvature
        (
                amrex::Real *dPsi,
                amrex::Real &Curv,
                int i, int j,
                amrex::Array4<amrex::Real const> const &Psi,
                const amrex::Real &deltax,
                const amrex::Real* dx,
                const amrex::Real* prob_lo,
                bool axisymmetric
        );


	void Compute_LSQ_Weights(double &weights_, double dist_);

	void LSQ_Parameters
	(
		double &P_,
		double &x_,
		double &y_,
		double &w_,
		int i1, int j1,
		int i, int j,
		amrex::Array4<amrex::Real const> const &P,
		const amrex::Real *prob_lo,
		const amrex::Real *dx
	);

        void LSQ_Parameters_Fij
        (
                double &Fij_,
		int &iscalar,
                double &x_,
                double &y_,
                double &w_,
                int i1, int j1,
                int i, int j,
                amrex::Array4<amrex::Real const> const &P,
                const amrex::Real *prob_lo,
                const amrex::Real *dx,
		const amrex::Box &domain
        );


	void QR_LS_Pressure_Second_Order
	(
		bool &sing,
		int Num_,
		int Num_i_,
		double *x_, double *y_,
		double *P_,
		double *w_,
		amrex::Array4<amrex::Real> const &P,
		InterceptData &idata
	);

	void QR_LS_Pressure_Third_Order
	(
		bool &sing,
		int Num_,
		int Num_i_,
		double *x_,
		double *y_,
		double *P_,
		double *w_,
		amrex::Array4<amrex::Real> const &P,
		InterceptData &idata
	);

	void AllMCQR_LS_Pressure_Fourth_Order
	(
		bool &sing,
		int Num_,
		int Num_i_,
		double *x_,
		double *y_,
		double *P_,
		double *w_,
		amrex::Array4<amrex::Real> const &P,
		InterceptData &idata,
		const amrex::LinOpBCType &bc
	);

        void QR_LS_Temperature_Second_Order
        (
                bool &sing,
                int Num_,
                int Num_i_,
                double *x_, double *y_,
                double *P_,
                double *w_,
                amrex::Array4<amrex::Real> const &P,
                InterceptData &idata
        );

        void QR_LS_Temperature_Third_Order
        (
                bool &sing,
                int Num_,
                int Num_i_,
                double *x_,
                double *y_,
                double *P_,
                double *w_,
                amrex::Array4<amrex::Real> const &P,
                InterceptData &idata
        );


	/************************************************************************/

	/****************************   data    *********************************/
        amrex::Vector<std::string> IF_names;
        amrex::Vector<std::string> IF_types;
	amrex::iMultiFab mask_;
	amrex::iMultiFab maskTemp_;
        amrex::iMultiFab PMask, UMask, VMask, PMaskTemp;
	amrex::iMultiFab nbr_flag_;// for tagging cells outside interface which has been assigned as a nbr to the cell in intercept
	amrexMesh *mesh_;
	amrex::Geometry geom_;
	amrex::MultiFab dummyF;
	

	amrex::LayoutData<std::array<std::vector<amrex::IntVect>, 12>> Index_Additional_Layers_LD;
	int LAYERS = 12;

	std::vector<std::unique_ptr<Interface>> *interfaces;

	bool CUBIC_SOLVE = true;
	bool DamageModel = false;
	bool PhaseField = false;

	int stencil_ = 4;
	int PORDER = 4;

	const double exp_ = 2.0;

        //Material Props required for Interface treatment
        amrex::Real SIGMA = 0.0;

	int LSQ_order = 4;//order of least square resconstruction  
	int Fij_order = 3;//order of least square resconstruction of the deformation gradient tensor

	/************************************************************************/

}; /*End class Mask */

} /*End namespace mycode */

#endif
