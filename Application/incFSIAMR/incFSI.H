#ifndef _incFSI_H_
#define _incFSI_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Mon Dec 21 16:40:52 2020


************************************************************************/

#include <amrexMesh.H>
#include <AMReX_MultiFab.H>
#include <AMReX_LO_BCTYPES.H>
#include <functional>
#include <defs.H>

#include <Interface.H>
#include <Cylinder.H>
#include <Ellipse.H>
#include <Foil.H>
#include <Bubble.H>
#include <InterceptData.H>

#include <Mask.H>

#include <AMReX_AmrCore.H>
#include <AMReX_PhysBCFunct.H>

namespace mycode
{

    struct interface_props
    {
    public:
        amrex::Real Int_p;
        amrex::Real Int_R;
        amrex::Real Int_R_dot;
        amrex::Real Int_xcp;
        amrex::Real Int_ycp;

    };

    class incFSI : public amrex::AmrCore
    {
    
    public:
    
        incFSI(const std::vector< std::unique_ptr<interface_props> >& IF_ref);
    
    // destructor:
        virtual ~incFSI();
    
    /**************************** functions *********************************/
    
        //! AMR related functions
        //! the virtual functions from AmrCore
        // Make a new level using provided BoxArray and DistributionMapping and
        // fill with interpolated coarse level data.
        // overrides the pure virtual function in AmrCore
        virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                             const amrex::DistributionMapping& dm) override;
    
        // Remake an existing level using provided BoxArray and DistributionMapping and
        // fill with existing fine and coarse data.
        // overrides the pure virtual function in AmrCore
        virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                  const amrex::DistributionMapping& dm) override;
    
        // Delete level data
        // overrides the pure virtual function in AmrCore
        virtual void ClearLevel (int lev) override;
    
        // Make a new level from scratch using provided BoxArray and DistributionMapping.
        // Only used during initialization.
        // overrides the pure virtual function in AmrCore
        virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                              const amrex::DistributionMapping& dm) override;
    
        // tag all cells for refinement
        // overrides the pure virtual function in AmrCore
        virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;
    
    	// compute a new multifab by coping in phi from valid region and filling ghost cells
        // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    	//! fill cell-centered mfs
        void FillPatch 
    	(
    		int lev, 
    		amrex::Real time, 
    		amrex::MultiFab& mf,                        // dest mf
    		const amrex::Vector<amrex::MultiFab*>& cmf, // coarse mf
    		const amrex::Vector<amrex::Real>& ct,       // coarse time
    		const amrex::Vector<amrex::MultiFab*>& fmf, // fine mf
    		const amrex::Vector<amrex::Real>& ft,       // fine time
    		int icomp, int ncomp
    	);
	//! fill cell-centered mfs
        void FillPatchScalars
        (
                int lev,
		int iscalar,
                amrex::Real time,
                amrex::MultiFab& mf,                        // dest mf
                const amrex::Vector<amrex::MultiFab*>& cmf, // coarse mf
                const amrex::Vector<amrex::Real>& ct,       // coarse time
                const amrex::Vector<amrex::MultiFab*>& fmf, // fine mf
                const amrex::Vector<amrex::Real>& ft,       // fine time
                int icomp, int ncomp
        );

    
    	//! fill face centered velocity
        void FillPatch
    	(
    		int lev, 
    		amrex::Real time, 
    		amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> const& mf,                  // dest mf
    		const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>& cmf,  // coarse mf
    		const amrex::Vector<amrex::Real>& ct,                                      // coarse time
    		const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>& fmf,  // fine mf
    		const amrex::Vector<amrex::Real>& ft,                                      // fine time
    		int icomp, int ncomp
    	);
    
        // fill an entire multifab by interpolating from the coarser level
        // this comes into play when a new level of refinement appears
        void FillCoarsePatch 
    	(
    		int lev, amrex::Real time, 
    		amrex::MultiFab& mf, 
    		const amrex::MultiFab& cmf, 
    		int icomp, int ncomp
    	);
    	
    	//! fill vel using div-free mapper
    	void FillCoarsePatch
    	(
    		int lev, amrex::Real time,
    		amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> const &mf,
    		const amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> &Cmf,
    		int icomp, int ncomp
    	);
    	//!
        void FillPatchAroundBox();
    
        void FillPatchNearBoundary(int lev);
    
        void Initialize();
  
        void InitializeRKcoeffs();
    
        void Evolve();

        void ApplyBC();
  
        void Regrid();

        void AverageDown();
    
        /// set Boundary functions
        void SetXVelBoundaryFunction(BoundaryFunc f);
    
        void SetXVelBoundaryFunction(int side, BoundaryFunc f);
    
        void SetYVelBoundaryFunction(BoundaryFunc f);
    
        void SetYVelBoundaryFunction(int side, BoundaryFunc f);
    
        void SetPressureBoundaryFunction(BoundaryFunc f);
    
        void SetPressureBoundaryFunction(int side, BoundaryFunc f);
    
        void SetTemperatureBoundaryFunction(BoundaryFunc f);
    
        void SetTemperatureBoundaryFunction(int side, BoundaryFunc f);
    
        void SetPhiBoundaryFunction(BoundaryFunc f);
    
        void SetPhiBoundaryFunction(int side, BoundaryFunc f);

        void SetScalarsBoundaryFunction(int iscalar, BoundaryFunc f);

        void SetScalarsBoundaryFunction(int iscalar, int side, BoundaryFunc f);
    
        void setInitXVel(FieldFunc f) { u_init = f; }
    
        void setInitYVel(FieldFunc f) { v_init = f; }
    
        void setInitPressure(FieldFunc f) { P_init = f; }
    
        void setInitTemperature(FieldFunc f) { Theta_init = f; }
    
        void setInitPhi(FieldFunc f) { Phi_init = f; }
    
        //void setIntBCs(std::vector<amrex::LinOpBCType> &intbc, const char *s);
    
        //void SetInterfaceValues(int lev);
        void ComputeIntProp(int lev);
    
        void ComputeBubbleVolume(bool init);

	    void CheckConservationPhaseField();

        int getReinitPt(){return chk_int_read;}

        void setReinitPt(int new_chk_int_read){chk_int_read=new_chk_int_read;}

        int getMaxIter(){return MaxIter;}

        void setMaxIter(int new_MaxIter){MaxIter=new_MaxIter;}
    
        amrex::Real Int_p,Int_R, Int_ur, Int_xcp,Int_ycp;
    
    
    /************************************************************************/
    
    private:
    /**************************** functions *********************************/
    
        void AllocateMemory(int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);
    
        void SetInitialFlowField(int lev);
    
        void XVelBoundaryConditions(int lev, amrex::Real time, amrex::MultiFab& mf);
    
        void YVelBoundaryConditions(int lev, amrex::Real time, amrex::MultiFab& mf);
    
        void PressureBoundaryConditions(int lev, amrex::Real time, amrex::MultiFab& mf);
    
        void TemperatureBoundaryConditions(int lev, amrex::Real time, amrex::MultiFab& mf);
     
        void PhiBoundaryConditions(int lev, amrex::Real time, amrex::MultiFab& mf);

        void ScalarBoundaryConditions(int lev,int iscalar, amrex::Real time, amrex::MultiFab& mf);
        
        void CollocatedVelocityBoundaryConditions(int lev, amrex::Real time, amrex::MultiFab& mf);
    
        void AMGPressurePoisson();

        void AMGPressurePoisson_RK(int RKStage);
    
        void Scheme();

        void Scheme_RK2();
   
        void Scheme_TVDRK2();

        void Scheme_RK();
    
        void CopyNewToOld();

        void CopyNewToRK1();

        void CopyNewToRK2();

        void CopyNewToRK3();
    
        void ComputeIntermediateVelocity(int lev);

        void AdvectScalars(int lev);

        void AdvectScalars_RK2(int lev);
	 
        void ComputeDeformationGradientTensor2D(int lev);

        void ComputeIntermediateVelocity_RK2(int lev);

        void ComputeIntermediateVelocity_RK(int lev, int RKStage);
 
        void ComputeRKSourceTerms(int lev, int RKStage, int step);

        void ComputeRKSourceTerms_DamageModel(int lev, int RKStage, int step);
    
        void ComputeCollocatedVelocityField(int lev);
    
        void ComputeStaggeredVelocityField(int lev);
    
        void ComputeFinalVelocityField(int lev);

        void ComputeFinalVelocityField_RK(int lev, int RKStage);
    
        void Derive();
    
        void ComputeDt();
    
        bool IsSteady();
    
        void WriteFile();
    	
    	void WriteFileParaview();
    	
        void writeDataFile();
    
        void WriteFileTecplot();
    
        void WriteFileTecplotW_Ghost();

        void ComputeMomentumBalance();
        
        void WriteInterface();
        
        amrex::Box get_valid_face_box(int lev, const amrex::Box& bx_i, int dir);
    
        void setBCs(BCArray &bcarr, const char *s);
    
        void CopyFRKVariables();

        void TVD_RK2Avg();
 
        void RK2Avg();

        void RK2Avg_Psi(amrex::MultiFab&, const amrex::MultiFab&);

        void CopyFRKPsi(amrex::MultiFab&, const amrex::MultiFab&);

        /***********Interface sharpening****************************************/
        amrex::Real GetPhaseFieldRHS(amrex::Array4<amrex::Real const> const& Psi, int i, int j, int k);

        /***********Function for embedded bc************************************/
        
        void MakeInterface( int lev,  const amrex::BoxArray& ba, const amrex::DistributionMapping& dm );
    
        void RemakeInterface( int lev, amrex::Real time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm );

        void MakeNewLevelFromCoarse_interface(int lev, amrex::Real time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

        void AverageDownInterface();
    
        void ComputeCutFaceVel(int lev);
        
        void ComputeCutFacePressure(int lev);
    
/************************************************************************/
    
        //! Restart
        void WriteCheckpointFile(const std::string &file);
    
        void ReadCheckpointFile();
    
        void ChkFile();

        std::pair<double,double> convertToPhysicalCoordinate(int i, int j, int lev); 

        std::pair<int,int> convertToGridCoordinate(double x, double y, int lev);
    
/************************************************************************/
/************************staggered grid info****************************/

        amrex::Vector<amrex::BoxArray> grids_xvel,grids_yvel,grids_zvel;
    
/****************************   data    *********************************/

        /// run parameters
        bool isAxisymmetric = false;
        amrex::Real dt, Time, FinalTime;
        int Iter, PlotInt, MaxIter, ReinitInt;
        bool Projection = false;    
        //!Restart
        bool Restart = false;
        std::string chk_file;
        int chk_int = -1;
        int chk_int_read = 0;
    
        ///solver options
        int RKOrder = 1;
        int Fij_order = 5;//WENO order for evolving deformation gradient tensor
        bool divergence_free_interpolation = true;
        bool TempField = false;
        bool PhaseField = false; 
        bool sharpPhaseField = true;
        bool use_algoim = false;
        bool tecplot = false;
        bool ComputeMomentumBalance_ = false;
        amrex::Real CFL = 0.5;
        amrex::Real b_sharp = 0.0;//b in eq.(12) Sun and beckermann, jcp, 2007 
        amrex::Real W_phasefield = 2.0;
        amrex::Real b_coeff = 1.0;
        bool DamageModel = false;
        amrex::Real damage_coeff = 0.2;//A parameter controlling the time-sclae of damage
        bool advect_ref_cond = false;
        int nscalar = 9;//X,Y,F11,F12,F21,F22,F33,max principal strain,Damage indicator
        int X_ref_num = 0;
        int Y_ref_num = 1;
        int F11_num = 2;
        int F12_num = 3;
        int F21_num = 4;
        int F22_num = 5;
        int F33_num = 6;
        int eps_max_num = 7;
        int dmg_num = 8;
        int Reinit_Method = 1;
    
//Material type
        int material_type = 1;
        int viscous_fluid = 1;
        int KelvinVoigt = 2;
//variables to plot
        bool plot_strain_rate = true;
        bool plot_strain_error = false;
        bool plot_viscosity = false;
        int LS_reg_iter = 0;
        
/// after how many steps to do regrid
        int regrid_int = 2;
        bool tag_vort = false;
        bool tag_region = false;
        bool tag_interface = false;
        int N_tag_regions;
        amrex::Vector<std::string> tag_region_name;
        amrex::Vector<amrex::Real> tag_region_lo, tag_region_hi;
        amrex::Vector<amrex::BoxList> boxlist;
    
        bool StopAtSteady = false;
        bool Steady = false;
        amrex::Real SteadyDeviation;
        amrex::Real SteadyTol = 1e-5;
        amrex::Real MaxU;
    
//! time for each level
        amrex::Real t_new, t_old;
    
        ///
        int Nghost,TUBE_LAYERS;
    
/// field variables
        amrex::Vector<amrex::MultiFab> xvel, yvel, U, Pressure, RHS, Theta, Phi;
        amrex::Vector<amrex::MultiFab> xvel_old, yvel_old, Pressure_old, Pstar, Theta_old, Phi_old;
        amrex::Vector<amrex::MultiFab> Umag, UmagOld, Vorticity, Dissipation;

        amrex::Vector<amrex::MultiFab> FRK_xvel, FRK_yvel, FRK_p, FRK_Theta, FRK_Phi;//Previous soln
        amrex::Vector<amrex::MultiFab> RK1_xvel, RK1_yvel, RK1_p, RK1_Theta, RK1_Phi;//First RK stage
        amrex::Vector<amrex::MultiFab> RK2_xvel, RK2_yvel, RK2_p, RK2_Theta, RK2_Phi;//Second RK stage
        amrex::Vector<amrex::MultiFab> RK3_xvel, RK3_yvel, RK3_p, RK3_Theta, RK3_Phi;//Second RK stage
        amrex::Vector<amrex::MultiFab> Src_xvel, Src_yvel, Src_p, Src_Theta, Src_Phi;//Second RK stage
        
        amrex::Vector<amrex::Vector<amrex::MultiFab>> FRK_Psi; // only if Advect level set
        amrex::Vector<amrex::Vector<amrex::MultiFab>> Scalars, Scalars_old, FRK_Scalars, RK1_Scalars, RK2_Scalars, RK3_Scalars, Src_Scalars;
        //amrex::Array<Vector<amrex::MultiFab, 2> FRK_Psi;

    
        //amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> flux;
    
/// properties
        amrex::Real Mu, Mu_max;
        amrex::Real gravity = 0.0;
        Viscosity visc_; // variable viscosity
        amrex::Real Polytropic_Index = 1.4;
        amrex::Real SIGMA = 1.450000000000e-01;
        amrex::Real lambda = 1.0;
        amrex::Real k_max;
        amrex::Real Prandtl_no = 6.13;
        amrex::Real u_ref = 10.0;
        amrex::Real T_ref = 367.0;
        amrex::Real cp_ = 4200.0;
        amrex::Real rho_ref = 1000.0;
    
    
/// BCs
/// on low faces (left and bottom)
        BCArray lobc_u;
        BCArray lobc_v;
        BCArray lobc_p;
        BCArray lobc_T;
        BCArray lobc_phi;
//amrex::Vector<std::unique_ptr<BCArray>> lobc_Scalars;
        amrex::Vector<BCArray> lobc_Scalars;
    
/// on hi faces (right and top)
        BCArray hibc_u;
        BCArray hibc_v;
        BCArray hibc_p;
        BCArray hibc_T;
        BCArray hibc_phi;
//amrex::Vector<std::unique_ptr<BCArray>> hibc_Scalars;
        amrex::Vector<BCArray> hibc_Scalars;
    
        /// Boundary condition function
        amrex::Array<BoundaryFunc, 2*AMREX_SPACEDIM> u_bcf;
        amrex::Array<BoundaryFunc, 2*AMREX_SPACEDIM> v_bcf;
        amrex::Array<BoundaryFunc, 2*AMREX_SPACEDIM> p_bcf;
        amrex::Array<BoundaryFunc, 2*AMREX_SPACEDIM> T_bcf;
        amrex::Array<BoundaryFunc, 2*AMREX_SPACEDIM> Phi_bcf;
	    amrex::Vector<amrex::Array<BoundaryFunc, 2*AMREX_SPACEDIM>> Scalars_bcf;
    
        /// initial conditions
        FieldFunc u_init, v_init, P_init, Theta_init, Phi_init;
    
    	
        int poisson_iter, system_size;
        amrex::Real poisson_tol;
        /********************Interface retaled data and function*******************/
        //! collection of interfaces
        int N_IF = 0;
        amrex::Vector<std::string> IF_names;
        amrex::Vector<std::string> IF_types;
        amrex::Vector<std::vector<std::unique_ptr<Interface>>> interfaces;
    
        //! mask
        amrex::Vector<std::unique_ptr<Mask>> mask_;
    
        //! interface boundary functions
        std::vector<IntBoundaryFunc> intbcf;
        //! interface boundary conditions
        std::vector<amrex::LinOpBCType> Pintbcs,Uintbcs,Vintbcs,Tintbcs,Phiintbcs;
    
        void setIntBCs(std::vector<amrex::LinOpBCType> &intbc, const char *s);
    
        const std::vector<std::unique_ptr<interface_props>> &IF_props ;
    
        /************************************************************************/
        /********************RK3 Coefficients**********************************************/
        double rk_a[4][4] = {{1.0/2.0,0.0, 0.0, 0.0}, {-1.0, 2.0, 0.0, 0.0}, {1.0/6.0, 4.0/6.0, 1.0/6.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
        double rk_c[4] = {1.0/2.0, 1.0, 1.0, 0.0};


        /**********Control volume calculations********************************************/
        amrex::Real mass_cv_old = 0.0;
        amrex::Real mass_cv = 0.0;
        amrex::Real x_mom_cv_old = 0.0;
        amrex::Real y_mom_cv_old = 0.0;
        amrex::Real x_mom_cv = 0.0;
        amrex::Real y_mom_cv = 0.0;

    
    }; /*End class incFSI */
    
    // Appply boundary condition for velocity while populaing AMR patch 
    struct MyExtBCFill 
    {
        int velflag;// 1 = u, 2 = v, 3 = w
        amrex::BoxArray ba;
        BCArray lobc;
        BCArray hibc;
        amrex::Array<BoundaryFunc, 2*AMREX_SPACEDIM> bcf;
        bool isAxisymmetric;
        AMREX_GPU_DEVICE
        void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& dest,
                         const int dcomp, const int numcomp,
                         amrex::GeometryData const& geom, const amrex::Real time,
                         const amrex::BCRec* bcr, const int bcomp,
                         const int orig_comp) const
        {
            using namespace amrex;
    
            // do something for external Dirichlet (BCType::ext_dir)
            const int i = iv[0];
            const int j = iv[1];
            const int k = 0;
            const amrex::Box &domain = geom.Domain();
            const amrex::Real *prob_lo = geom.ProbLo();
            const amrex::Real* prob_hi = geom.ProbHi();
            const amrex::Real* dx = geom.CellSize();
    
            if(velflag == 1)//Boundaru condition for u velocity
            {
                amrex::Real x = prob_lo[0] + dx[0] * i;
                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                if(i == domain.smallEnd(0))//Left boundary
                {
                    //int shift = i - domain.smallEnd(0);
                    //if(shift == 0)
                    {
                        if (lobc[0] == amrex::LinOpBCType::Dirichlet)
                            dest(i , j, k) = bcf[0](x, y, time);
                        else if (lobc[0] == amrex::LinOpBCType::Neumann)
                        {
                            //dest(i, j, k) = dest(i + 1, j, k) - dx[0] * bcf[0](x, y, time);
                            dest(i - 1, j, k) = dest(i + 1, j, k) - 2.0*dx[0] * bcf[0](x, y, time);
                        }
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else if(j == domain.smallEnd(1) - 1)//bottom boundary
                {
                    int shift = j - domain.smallEnd(1);
                    if(shift == -1)
                    {
                        if (lobc[1] == amrex::LinOpBCType::Dirichlet)
                            dest(i , j, k) = 2.0 * bcf[1](x, y, time) - dest(i , j + 1, k);
                        else if (lobc[1] == amrex::LinOpBCType::Neumann)
                            dest(i , j, k) = -1.0* dx[1] * bcf[1](x, y, time) + dest(i , j + 1, k);
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else if(i == domain.bigEnd(0) + 1)//right boundary
                {
                    //int shift = i - (domain.bigEnd(0) + 1);
                    //if(shift == 0)
                    {
                        if (hibc[0] == amrex::LinOpBCType::Dirichlet)
                            dest(i , j, k) = bcf[2](x, y, time) ;
                        else if (hibc[0] == amrex::LinOpBCType::Neumann)
                        { 
                            //dest(i, j, k) = dx[0] * bcf[2](x, y, time) + dest(i - 1, j, k); 
                            dest(i + 1, j, k) = 2.0 * dx[0] * bcf[2](x, y, time) + dest(i - 1, j, k);
                        }
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else if(j == domain.bigEnd(1) + 1)//top boundary
                {
                    //int shift = j - domain.bigEnd(1);
                    //if(shift == 1)
                    {
                        if (hibc[1] == amrex::LinOpBCType::Dirichlet)
                            dest(i , j, k) = 2.0 * bcf[3](x, y, time) + dest(i , j - 1 ,k);
                        else if (hibc[1] == amrex::LinOpBCType::Neumann)
                            dest(i , j, k) = dx[1] * bcf[3](x, y, time) + dest(i , j - 1 ,k);
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else
                    dest(i, j, k) = 0.0;
               //dest(i, j, k) = 0.0;
            }
            if(velflag == 2)//Boundary condition for V velocity
            {   
                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5); 
                amrex::Real y = prob_lo[1] + dx[1] * j;
                if(i == domain.smallEnd(0) - 1)//Left boundary
                {
                    //int shift = i - domain.smallEnd(0);
                    //if(shift == -1)
                    {
                        if (lobc[0] == amrex::LinOpBCType::Dirichlet)
                            dest(i , j, k) = 2.0 * bcf[0](x, y, time) - dest(i + 1, j, k);
                        else if (lobc[0] == amrex::LinOpBCType::Neumann)
                            dest(i , j, k) = -1.0* dx[0] * bcf[0](x, y, time) + dest(i + 1, j, k);
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else if(j == domain.smallEnd(1))//bottom boundary
                {
                    //int shift = j - domain.smallEnd(1);
                    //if(shift == 0)
                    {
                        if (lobc[1] == amrex::LinOpBCType::Dirichlet) 
                            dest(i , j, k) = bcf[1](x, y, time);
                        else if (lobc[1] == amrex::LinOpBCType::Neumann)
                        {   
                            if(isAxisymmetric)
                            {
                                dest(i, j, k) = 0.0;
                            }
                            //else
                            //{
                            //    dest(i , j, k) =  -1.0 * dx[1] * bcf[1](x, y, time) + dest(i, j + 1, k);
                            //}
                            dest(i , j - 1, k) = -2.0 * dx[1] * bcf[1](x, y, time) + dest(i, j + 1, k);
                        }
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else if(i == domain.bigEnd(0) + 1)//right boundary
                {
                    //int shift = i - domain.bigEnd(0) ;
                    //if(shift == 1)
                    {
                        if (hibc[0] == amrex::LinOpBCType::Dirichlet)
                            dest(i , j, k) = 2.0 * bcf[2](x, y, time) - dest(i - 1, j, k);
                        else if (hibc[0] == amrex::LinOpBCType::Neumann)
                            dest(i , j, k) = dx[0] * bcf[2](x, y, time) + dest(i - 1, j, k);
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else if(j == domain.bigEnd(1) + 1)//top boundary
                {
                    //int shift = j - (domain.bigEnd(1) + 1);
                    //if(shift == 0)
                    {
                        if (hibc[1] == amrex::LinOpBCType::Dirichlet) 
                            dest(i , j, k) = bcf[3](x, y, time);
                        else if (hibc[1] == amrex::LinOpBCType::Neumann)
                        {
                            //dest(i , j, k) = dx[1] * bcf[3](x, y, time) + dest(i, j - 1, k);
                            dest(i , j + 1, k) = 2.0 * dx[1] * bcf[3](x, y, time) + dest(i, j - 1, k);
                        }
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else
                    dest(i, j, k) = 0.0;
            }
            //dest(i, j, k) = 0.0;
        }
    };
    
    // Appply boundary condition for Scalars while populaing AMR patch 
    struct MyExtBCFillScalars 
    {
        int scalarflag;// 0 = X, 1 = Y,.....
        amrex::BoxArray ba;
        BCArray lobc;
        BCArray hibc;
        amrex::Array<BoundaryFunc, 2*AMREX_SPACEDIM> Scalars_bcf;
        bool isAxisymmetric;
        AMREX_GPU_DEVICE
        void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& scalar,
                         const int dcomp, const int numcomp,
                         amrex::GeometryData const& geom, const amrex::Real time,
                         const amrex::BCRec* bcr, const int bcomp,
                         const int orig_comp) const
        {
            using namespace amrex;
    
            // do something for external Dirichlet (BCType::ext_dir)
            const int i = iv[0];
            const int j = iv[1];
            const int k = 0;
            const amrex::Box &domain = geom.Domain();
            const amrex::Real *prob_lo = geom.ProbLo();
            const amrex::Real* prob_hi = geom.ProbHi();
            const amrex::Real* dx = geom.CellSize();
    
            {
		//amrex::Print()<<"demain = "<<domain<<'\n';
		//amrex::Print()<<" in Scalar bc amr:"<<i<<" , "<<j<<" , "<<k<<" , dx = "<<dx[0]<<'\n';
                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                if(i == domain.smallEnd(0) - 1)//Left boundary
                {
                    //int shift = i - domain.smallEnd(0);
                    //if(shift == 0)
		    x  = prob_lo[0];
                    {
                        if (lobc[0] == amrex::LinOpBCType::Dirichlet)
		        {
                            scalar(i, j, k) = 2.0 * Scalars_bcf[0](x, y, time) - scalar(i + 1, j, k);
                        }
                        else if (lobc[0] == amrex::LinOpBCType::Neumann)
                        {
			    scalar(i, j, k) = scalar(i + 1, j, k) - dx[0] * Scalars_bcf[0](x, y, time);
                        }
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else if(j == domain.smallEnd(1) - 1)//bottom boundary
                {
                    //int shift = j - domain.smallEnd(1);
		    y = prob_lo[1];
                    //if(shift == -1)
                    {
                        if (lobc[1] == amrex::LinOpBCType::Dirichlet)
			    scalar(i, j , k) = 2.0 * Scalars_bcf[1](x, y, time) - scalar(i, j + 1, k);
                        else if (lobc[1] == amrex::LinOpBCType::Neumann)
			    scalar(i, j , k) = scalar(i, j + 1, k) - dx[1] * Scalars_bcf[1](x, y, time);
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else if(i == domain.bigEnd(0) + 1)//right boundary
                {
                    //int shift = i - (domain.bigEnd(0) + 1);
		    x = prob_hi[0];
                    //if(shift == 0)
                    {
                        if (hibc[0] == amrex::LinOpBCType::Dirichlet)
			    scalar(i , j, k) = 2.0 * Scalars_bcf[2](x, y, time) - scalar(i - 1, j, k);
                        else if (hibc[0] == amrex::LinOpBCType::Neumann)
                        { 
			    scalar(i , j, k) = scalar(i - 1, j, k) + dx[0] * Scalars_bcf[2](x, y, time);
                        }
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else if(j == domain.bigEnd(1) + 1)//top boundary
                {
                    //int shift = j - domain.bigEnd(1);
		    y = prob_hi[1];
                    //if(shift == 1)
                    {
                        if (hibc[1] == amrex::LinOpBCType::Dirichlet)
			    scalar(i, j , k) = 2.0 * Scalars_bcf[3](x, y, time) - scalar(i, j - 1, k);
                        else if (hibc[1] == amrex::LinOpBCType::Neumann)
			    scalar(i, j , k) = scalar(i, j - 1, k) + dx[1] * Scalars_bcf[3](x, y, time);
                    }
                    //else
                    //    dest(i , j, k) = 0.0;
                }
                else
                    scalar(i, j, k) = 0.0;
               //dest(i, j, k) = 0.0;
            }
        }
    };

} /*End namespace mycode */

#endif
