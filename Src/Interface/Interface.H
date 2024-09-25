#ifndef _Interface_H_
#define _Interface_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Wed Dec 30 11:31:30 2020


************************************************************************/

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <amrexMesh.H>
#include <InterceptData.H>
#include <CFMask.H>

namespace mycode
{

class Interface
{
public:

	Interface(const amrexMesh &, const std::string& name);

        void ClearInterfaceData();
        void RemakeInterface(const amrexMesh &newmesh, const amrex::MultiFab &Psi);
        void MakeInterfaceFromCoarse(const amrexMesh &newmesh, const amrex::MultiFab &Psi);

	virtual ~Interface();

	/**************************** functions *********************************/

	//! member functions common to all
	/////////////////////////////////////////////////////////////////////////

	void PhaseFieldBC();
	void Reinit_algoim();
        void Reinit2();
	amrex::MultiFab &Psi() { return psi; }
	amrex::MultiFab &Kappa_() { return kappa_; }
	amrex::MultiFab &Normal_() { return normal; }
	amrex::Real &max_kappa_() { return max_kappa_interface; }
	const amrex::MultiFab &Psi() const { return psi; }
        void setPsi(const amrex::MultiFab &newpsi){ amrex::MultiFab::Copy(psi, newpsi, 0, 0, 1, nghost); }
        amrex::MultiFab &FRK_Psi() { return FRK_psi; }
        const amrex::MultiFab &FRK_Psi() const { return FRK_psi; }
        void copyFRK_Psi();
        void copyRK1_Psi();
        void copyRK2_Psi();
        void copyRK3_Psi();
        void setFRK_Psi(const amrex::MultiFab &newFRKpsi){ amrex::MultiFab::Copy(FRK_psi, newFRKpsi, 0, 0, 1, nghost); }
        void TVDRK2Avg_Psi();
        void RK2Avg_Psi();

	amrex::LayoutData<std::vector<InterceptData>> &getInterceptData() 
	{
		return intercept_data; 
	}
	
	const amrex::LayoutData<std::vector<InterceptData>> &getInterceptData() const 
	{ 
		return intercept_data; 
	}

	bool isAdvectLS() { return is_advect_ls; }
        amrex::Real Volume_prev(){ return Volume_prev_;}
        void setVolume_prev(){ Volume_prev_ = Volume_;}

	//////////////////////////////////////////////////////////////////////////

	//! virtual functions to be implemented in derived classes
	//////////////////////////////////////////////////////////////////////////

	virtual void PrescribeLevelSetMotion(const amrex::Real& t) = 0;

	//! functions for prescribed level set
	virtual void ComputeLevelSetMotion(const amrex::Real& t, const amrex::Real& dt, int RKStage) = 0;
	virtual amrex::RealArray getOrigin() = 0;

	//! functions for advected levelset
	virtual void TubeIdentification() = 0;
	virtual void TubeAdvectLevelSet(const amrex::MultiFab &u, const amrex::MultiFab &v, const amrex::Real &dt) = 0;
        virtual void TubeAdvectLevelSet_RK2(const amrex::MultiFab &u, const amrex::MultiFab &v, const amrex::Real &dt) = 0;
        virtual void TubeAdvectLevelSet_RK3(int RKStage, const amrex::MultiFab &u, const amrex::MultiFab &v, const amrex::Real &dt) = 0;
	//virtual void Regularization() = 0;


	//////////////////////////////////////////////////////////////////////////

	//! member functions specific to prescribed level set
	//////////////////////////////////////////////////////////////////////////

	void SetXcp(const amrex::Real &x) { xcp = x; }
	void SetYcp(const amrex::Real &y) { ycp = y; }
	void SetThetacp(const amrex::Real &t) { thetacp = t; }
	void SetDdt_xcp(const amrex::Real &u) { dxcpdt = u; }
	void SetDdt_ycp(const amrex::Real &v) { dycpdt = v; }
	void SetDdt_thetacp(const amrex::Real &w) { dthetacpdt = w; }
	void SetD2dt2_thetacp(const amrex::Real &w) { d2thetacpdt2 = w; }

	int &N_Intercept() { return N_intercepts; }
	const int& N_Intercept() const { return N_intercepts; }
	const amrex::Real &Xcp() { return xcp; }
	const amrex::Real &Ycp() { return ycp; }
	const amrex::Real &Thetacp() { return thetacp; }
	const amrex::Real &DTheta() { return dtheta; }
	const amrex::Real &ddt_xcp() { return dxcpdt; }
	const amrex::Real &ddt_ycp() { return dycpdt; }
	const amrex::Real &ddt_thetacp() { return dthetacpdt; }
	const amrex::Real &d2dt2_thetacp() { return d2thetacpdt2; }
	const amrex::Real &frequency() { return Freq; }

	//////////////////////////////////////////////////////////////////////////

	//! member functions specific to advected level set
	//////////////////////////////////////////////////////////////////////////

	void setLBeta(const amrex::Real &lb) { L_BETA = lb; }
	void setLGamma(const amrex::Real &lg) { L_GAMMA = lg; }

	amrex::iMultiFab &Mask() { return Mask_; }
	auto &Index() { return Index_; }
	int &TTube() { return TTube_; }
	int &NTube() { return NTube_; }

	void ComputeVolume(bool intit = false);
	void ComputeVolume(bool intit , const amrex::iMultiFab &cfmask_);
	void SetVolume(amrex::Real Vol){ Volume_ = Vol;}
	void SetVolume_0(amrex::Real Vol){ Volume0_ = Vol;}
	amrex::Real Volume() { return Volume_; }
	amrex::Real Volume_0() { return Volume0_; }
	amrex::Real Vol_Diff() { return Vol_Diff_; }
	amrex::Real Vol_Diff_0() { return Vol_Diff0_; }
	int getResolution0() { return resolution0_;}
	int getResolution() { return resolution_;}
	void SetTime(amrex::Real Time){ time = Time;}
	amrex::Real Time() { return time; }

	amrex::Real Part_Volume() { return Part_Volume_; }
	amrex::Real Part_Volume_0() { return Part_Volume0_; }
	amrex::Real Part_Vol_Diff() { return Part_Vol_Diff_; }
	amrex::Real Part_Vol_Diff_0() { return Part_Vol_Diff0_; }
	int getPartResolution0() { return Part_resolution0;}
	int getPartResolution() { return Part_resolution;}

	amrex::Real surface_area() { return surf_area; }
	amrex::Real part_surface_area() { return part_surf_area; }

	void setSphericity(amrex::Real val){ sphericity = val;}
	amrex::Real getSphericity(){ return sphericity;}

	void setVolume(const amrex::Real &v) { Volume_ = v; }
	void setVolume0(const amrex::Real &v) { Volume0_ = v; }

	void setResolution(const int &v) { resolution_ = v; }
	void setResolution0(const int &v) { resolution0_ = v; }

	//////////////////////////////////////////////////////////////////////////

	const std::string &name() { return Name; }

	/****************************   data    *********************************/

	//! data specific to prescribed level set
	/// forces
	amrex::Real FPx, FPy, FVx, FVy, MomP_z, MomV_z;

	//! Fourier qty
	struct FourierQty
	{
		amrex::Real FPx, FPy, FVx, FVy, MomP_z, MomV_z;
	}Fourier;

	int has_jet = 0;
	int has_spekle = 0;
	

	/************************************************************************/
        //Pressure at interface for Bubble simualtions
        amrex::Real P_interface0,P_interface;
        int bdy_frc_type = 0;
        int const_pressure_gradient = 1;
        int const_force = 2;
	amrex::Real bdy_frc = 0.0;//body force on bubble
	amrex::Real F_visc_nn = 0.0;
	amrex::Real F_visc_nt = 0.0;
	amrex::Real F_p = 0.0;
	amrex::Real F_tot = 0.0;

        void setP_interface(const amrex::Real &newP_interface){ P_interface = newP_interface;}
        amrex::Real getP_interface(){ return P_interface;}
        amrex::Real getP_interface0(){ return P_interface0;}


        void ComputeAvgIntVel(const amrex::MultiFab &u, const amrex::MultiFab &v, const amrex::iMultiFab &PMask);
        amrex::Real AvgRadius() { return AvgRadius_;}
        amrex::Real AvgIntVel() { return AvgIntVel_;}
        amrex::Real AvgIntPres() { return AvgIntPres_;}

	void SetAvgIntStrain(const amrex::Real &ep_avg) { AvgIntStrain_ = ep_avg; }
	amrex::Real AvgIntStrain() { return AvgIntStrain_;}
 
        amrexMesh getMesh() {return mesh_;}

	void Compute_Normal_Curvature();

        void displayMesh();

	void Regularization();

	void DetectUnderresolvedJetAlongAxis();

	void RemoveJetCellByCell();

	void Regularization2(int MAX_ITER);

	void RemoveSpeckles();

        /**************************** functions *********************************/

        void Reinit();

        int getNumReintIter() { return LS_Iter; }

        amrex::Real getLSReinitTolerance() { return LS_REINIT_TOL; }

        /************************************************************************/
        double rk_a[4][4] = {{1.0/2.0,0.0, 0.0, 0.0}, {-1.0, 2.0, 0.0, 0.0}, {1.0/6.0, 4.0/6.0, 1.0/6.0, 0.0}, {0.0,0.0, 0.0, 0.0}};
        double rk_c[4] = {1.0/2.0, 1.0, 1.0, 0.0};

protected:
	/**************************** functions *********************************/

        void computeVolumeAxisymmetric(bool);
        void computeVolumeAxisymmetric(bool, const amrex::iMultiFab &cfmask_);

	amrex::Real ComputePhiFromPsi(const amrex::Real&);
        //amrex::Real ComputeVolumeFraction(const amrex::Real&, const amrex::Real&, const amrex::Real&,
        //                                  const amrex::Real&, const amrex::Real&, const amrex::Real&,
        //                                  const amrex::Real&,
        //                                  const amrex::Real&, const amrex::Real&, const amrex::Real&);
        amrex::Real Function2(amrex::Real);

	/****************************   data    *********************************/

	//! boxArray, geometry, dmap
	amrexMesh mesh_;

	//! name of interface to read from inputs
	std::string Name;

	//! is the fluid inside or outside of interface
	bool is_fluid_in;

	//! interface to be advected?
	bool is_advect_ls;

        //! axisymmetric
        bool isAxisymmetric = false;

	//! data specific to prescribed level set
	/////////////////////////////////////////////////////////////////////////
	
	amrex::Real xcp, ycp, thetacp, dtheta;
	amrex::Real dxcpdt, dycpdt, dthetacpdt, d2thetacpdt2;
	amrex::Real xcp0, ycp0, dxcpdt0, dycpdt0;
	amrex::Real xc0, yc0;
	amrex::Real Amp, Freq;

	/// the total number of cut cells
	int N_intercepts;
	int LAYERS = 12;
        int nghost = 5;
        int RKOrder = 1;
	amrex::MultiFab psi, psi_old, psi_prev, source_reinit,FRK_psi, RK1_psi, RK2_psi, RK3_psi, kappa_,normal;
	amrex::LayoutData<std::vector<InterceptData> > intercept_data;

	//! flags related to DOF
	int DOF = 0;

	///////////////////////////////////////////////////////////////////////////

	//! data specific to advected level set
	///////////////////////////////////////////////////////////////////////////
	amrex::iMultiFab Mask_;
	amrex::LayoutData<std::vector<amrex::IntVect>> Index_;
	amrex::LayoutData<std::vector<amrex::Real>> Source;
	//! these are global no of cells to be advected and reinit
	//! for parallel ... need to identify from mask
	int TTube_, NTube_;
	int NTubes_global;

	//! Volume of interface
	amrex::Real Volume_, Volume0_, Vol_Diff_, Vol_Diff0_,Volume_prev_;
	amrex::Real Part_Volume_, Part_Volume0_, Part_Vol_Diff_, Part_Vol_Diff0_;
	amrex::Real part_surf_area = 0.0;
	amrex::Real surf_area = 0.0;
	amrex::Real sphericity = 1.0;

	//! average bubble radius(R) and speed of interface(R_DOT)
	amrex::Real AvgRadius_,AvgIntVel_,AvgIntPres_,AvgIntStrain_;
	amrex::Real time;
	int Part_resolution0, Part_resolution;
	int resolution0_, resolution_;

	//!
	amrex::Real L_BETA, L_GAMMA;

private:
        amrex::Real max_reinit_tol = 1.0e-6;
	amrex::Real max_kappa = 0.0;
	amrex::Real max_kappa_interface = 0.0;
        void compute_source_reinit();

        void compute_source_regularization();

        amrex::Real minmod(amrex::Real a, amrex::Real b)
        {
            if(a*b <= 0) return 0.0;
            else return std::copysign(std::min(std::fabs(a),std::fabs(b)),a);
        };

	//////////////////////////////////////////////////////////////////////////
        /**************************** LSReinit functions *********************************/

        void LSGammaIdentification();
        void LSSIdentification();
        void TubeComputeLevelSetRHS();
        void LSCIdentification();
        void LSSum();

        void AlgoimReinit();

        /************************************************************************/

        /****************************   data    *********************************/

        //const amrexMesh &mesh_;
        //Interface *interface_;

        int LS_Iter, LS_MAXITER;
        amrex::Real LS_REINIT_TOL;

        bool use_FMM;
        //int LAYERS;

        //! fields
        amrex::iMultiFab LS_Gamma;
        amrex::iMultiFab LS_S; // 4 comp
        // amrex::iMultiFab LS_S, LS_C; // 4 comp
        amrex::MultiFab LS_R, LS_Dealta, LS_Sum;
        amrex::MultiFab LSReinit_Source; // 4 comp

        /************************************************************************/


	/************************************************************************/

}; /*End class Interface */

} /*End namespace mycode */

#endif
