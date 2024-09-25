#include "Mask.H"
#include <AMReX_ParmParse.H>
#include <SolveCubicEqn.h>

namespace mycode
{

    void Mask::define
    (
        amrexMesh *mesh, 
        std::vector<std::unique_ptr<Interface>> *IF
    )    
    {
        mesh_ = mesh;
        geom_ = mesh_->geometry();
        interfaces = IF;
        
        //int nghost;
        amrex::ParmParse pp;
        pp.get("Nghost", nghost);
        pp.get("SIGMA", SIGMA);
        pp.query("Mu",Mu);
        pp.get("Axisymmetric",isAxisymmetric);
        pp.get("LAYERS",LAYERS);
	pp.query("DamageModel",DamageModel);
	pp.query("Phase_field",PhaseField);
	pp.query("LSQ_order",LSQ_order);
	pp.query("Fij_order",Fij_order);

	int N_IF;
	pp.get("N_interfaces", N_IF);
        IF_names.resize(N_IF);
        IF_types.resize(N_IF);
        if (N_IF > 0)
        {
            pp.getarr("IF_names", IF_names);
            pp.getarr("IF_types", IF_types);
        }
    
        mask_.define(mesh_->grid(), mesh_->dmap(), 1, nghost);
        maskTemp_.define(mesh_->grid(), mesh_->dmap(), 1, nghost);
        PMask.define(mesh_->grid(), mesh_->dmap(), 1, nghost);
        UMask.define(mesh_->grid(), mesh_->dmap(), 1, nghost);
        VMask.define(mesh_->grid(), mesh_->dmap(), 1, nghost);
        Index_Additional_Layers_LD.define(mesh_->grid(), mesh_->dmap());
	dummyF.define(mesh_->grid(), mesh_->dmap(), 3, nghost);
	nbr_flag_.define(mesh_->grid(), mesh_->dmap(), 1, nghost);

        //LAYERS = 12;
    }
    
    Mask::Mask
    (
        amrexMesh *mesh, 
        std::vector<std::unique_ptr<Interface>> *IF
    )
    {
        define(mesh, IF);
    }
    
    Mask::~Mask()
    {
    }
    
    void Mask::GhostCellIdentfication()
    {
        amrex::ParallelDescriptor::Barrier();
        mask_.setVal(1);
	nbr_flag_.setVal(0);
        
        if (interfaces->empty())
        {
            return;
        }
        
        const amrex::Real *dx = geom_.CellSize();
        const amrex::Real *prob_lo = geom_.ProbLo();
        const amrex::Box &domain = geom_.Domain();
        
        amrex::Real dPsi_L[2], dPsi_R[2];
        amrex::Real frac, Curv_L, Curv_R;
        amrex::Real cu_;
        amrex::Real distance;
        const amrex::Real dist_TOL = dx[0];
        bool Is_Intersect, Retain;
        
        for (auto &&solid : *interfaces)
        {
            solid->TubeIdentification();
            amrex::iMultiFab &objMask = solid->Mask();
            int &N_Intercept = solid->N_Intercept();
            N_Intercept = 0;
            
            int l = 0;
            for (amrex::MFIter mfi(mask_); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();
                amrex::Array4<int> const &pmask = mask_.array(mfi);
                amrex::Array4<int> const &pmask_temp = maskTemp_.array(mfi);
                //amrex::Array4<int> const &indicator = Indicator.array(mfi);
                amrex::Array4<amrex::Real const> const &psi = solid->Psi().const_array(mfi);
                amrex::Array4<int> const &obj_mask = objMask.array(mfi);
		amrex::Array4<int> const &nbr_flag = nbr_flag_.array(mfi);
                
                auto &icpt_data = solid->getInterceptData()[mfi];
                
                //! NOTE : before computing intercept data and adding it to vector
                //!        we should clear already present data if any
                //!        otherwise in each iteration new interceptData will be added to vector
                //!        and it will keep on growing
                
                if (!icpt_data.empty())
                {
                    icpt_data.clear();
                    // icpt_data.shrink_to_fit();
                }
                else
                {
                    icpt_data.reserve(100);
                }
                
                //! we need intercept data at 1 grown box
                //! rather than communicating later.. just compute on grown box
                //! this will add only few extra cells to be computed (1 or 2 in a direction)
                amrex::Box gbx(bx);
                gbx.grow(1);
                    
                /// take intersection with domain box as we do not want to compute for cells outside domain
                amrex::Box isect = gbx & domain;
                
                amrex::ParallelFor(isect,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept 
                {
                    if (psi(i, j, k) >= 0.0 )
                    {
                        pmask(i, j, k) = 0;
                        Is_Intersect = false;
                        int i_ = 0;
                        InterceptData idata;
			            idata.n_nbr = 0;
                        amrex::IntVect iv(AMREX_D_DECL(i,j,k));				    
                        amrex::IntVect iv_e(AMREX_D_DECL(i + 1,j,k));
                        amrex::IntVect iv_w(AMREX_D_DECL(i - 1,j,k));
                        amrex::IntVect iv_n(AMREX_D_DECL(i,j + 1,k));
                        amrex::IntVect iv_s(AMREX_D_DECL(i,j - 1,k));
                        //if(obj_mask(i ,j ,k) > 0)
                        {
                            /// check 4 nbrs
                            if (psi(i - 1, j, k) < 0.0 && isect.contains(iv_w))
                            {
                                Is_Intersect = true;
                                frac = psi(i, j, k) / (psi(i, j, k) - psi(i - 1, j, k));
                                
                                if (CUBIC_SOLVE)
                                {
                                    frac = frac - 1.0;
                                    cu_ = Cubic_Solve(frac, psi(i - 2, j, k), psi(i - 1, j, k), psi(i, j, k), psi(i + 1, j, k));
                                    frac = frac + 1.0;
                                    if (cu_ == -1)
                                    {
                                        amrex::Print() << "\nX Intercept: " << frac << "\n";
                                        exit(1);
                                    }
                                }
                                    // sanity check....
                                if (frac < 0.0 || frac > 1.0)
                                {
                                amrex::Print() << "\nIntercept out of bounds!" << frac << "\n";
                                amrex::Print() << "frac linear: i,i-1" << psi(i, j, k) / (psi(i, j, k) - psi(i - 1, j, k)) << "\n";
                                    amrex::Print() << psi(i - 2, j, k) << "\t" << psi(i - 1, j, k) << "\t" << psi(i, j, k) << "\t" << psi(i + 1, j, k) << "\n";
                                            exit(1);
                                }
                                idata.type_[i_] = 0;
                                idata.frac_[i_] = frac;
                                // Note we are not storing curvature or normals
                                // This can lead to redundant calculations at
                                // overlapping nodes but memory requirement will be low
                                //Compute_Normal_Curvature(dPsi_L, Curv_L, i - 1, j, psi, dx[0]);
                                //Compute_Normal_Curvature(dPsi_R, Curv_R, i, j, psi, dx[0]);
                                Compute_Normal_Curvature(dPsi_L, Curv_L, i - 1, j, psi, dx[0], dx, prob_lo, isAxisymmetric);
                                Compute_Normal_Curvature(dPsi_R, Curv_R, i, j, psi, dx[0], dx, prob_lo, isAxisymmetric);
                                
                                idata.psix_[i_] = frac * dPsi_L[0] + (1.0 - frac) * dPsi_R[0];
                                idata.psiy_[i_] = frac * dPsi_L[1] + (1.0 - frac) * dPsi_R[1];
                                idata.kappa_[i_] = frac * Curv_L + (1.0 - frac) * Curv_R;
				if(nbr_flag(iv_w) == 0)
			        {
                                    idata.nbr_cellid_.push_back(iv_w);
				    nbr_flag(iv_w) = 1;
				    idata.n_nbr++;
				}
                                i_++;
                            }
                            
                            if (psi(i + 1, j, k) < 0.0 && isect.contains(iv_e))
                            {
                                Is_Intersect = true;
                                frac = psi(i, j, k) / (psi(i + 1, j, k) - psi(i, j, k));
                                if (CUBIC_SOLVE)
                                {
                                    cu_ = Cubic_Solve(frac, psi(i - 1, j, k), psi(i, j, k), psi(i + 1, j, k), psi(i + 2, j, k));
                                    if (cu_ == -1)
                                    {
                                        amrex::Print() << "\nX Intercept: " << frac << "\n";
                                        exit(1);
                                    }
                                }
                                // sanity check....
                                if (frac > 0.0 || frac < -1.0)
                                {
                                    amrex::Print() << "\nIntercept out of bounds!" << frac << "\n";
                                    amrex::Print() << "frac linear: i,i+1" << psi(i, j, k) / (psi(i + 1, j, k) - psi(i, j, k)) << "\n";
                                    amrex::Print() << psi(i - 1, j, k) << "\t" << psi(i, j, k) << "\t" << psi(i + 1, j, k) << "\t" << psi(i + 2, j, k) << "\n";
                                    exit(1);
                                }
                                // check for the distance....
                                if (i_ == 1)
                                {
                                    Retain = false;
                                    distance = fabs(frac - idata.frac_[0]);
                                    if (distance > dist_TOL)
                                        Retain = true;
                                }
                                else
                                    Retain = true;
                                
                                if (Retain)
                                {
                                    if ((i_ == 0) || (fabs(frac) > fabs(idata.frac_[0])))
                                    {
                                        idata.type_[i_] = 0;
                                        idata.frac_[i_] = frac;
                                        Compute_Normal_Curvature(dPsi_L, Curv_L, i, j, psi, dx[0],dx, prob_lo, isAxisymmetric);
                                        Compute_Normal_Curvature(dPsi_R, Curv_R, i + 1, j, psi, dx[0],dx, prob_lo, isAxisymmetric);

                                
                                        idata.psix_[i_] = -frac * dPsi_R[0] + (1.0 + frac) * dPsi_L[0];
                                        idata.psiy_[i_] = -frac * dPsi_R[1] + (1.0 + frac) * dPsi_L[1];
                                        idata.kappa_[i_] = -frac * Curv_R + (1.0 + frac) * Curv_L;
                                    }
                                    else
                                    {
                                        /// latest data will be in index 0
                                        idata.type_[i_] = idata.type_[i_ - 1];
                                        idata.frac_[i_] = idata.frac_[i_ - 1];
                                        idata.psix_[i_] = idata.psix_[i_ - 1];
                                        idata.psiy_[i_] = idata.psiy_[i_ - 1];
                                        idata.kappa_[i_] = idata.kappa_[i_ - 1];
                                
                                        idata.type_[i_ - 1] = 0;
                                        idata.frac_[i_ - 1] = frac;
                                        Compute_Normal_Curvature(dPsi_L, Curv_L, i, j, psi, dx[0],dx, prob_lo, isAxisymmetric);
                                        Compute_Normal_Curvature(dPsi_R, Curv_R, i + 1, j, psi, dx[0],dx, prob_lo, isAxisymmetric);
                                        idata.psix_[i_ - 1] = -frac * dPsi_R[0] + (1.0 + frac) * dPsi_L[0];
                                        idata.psiy_[i_ - 1] = -frac * dPsi_R[1] + (1.0 + frac) * dPsi_L[1];
                                        idata.kappa_[i_ - 1] = -frac * Curv_R + (1.0 + frac) * Curv_L;
                                    }
                                    i_++;
                                }
                                if(nbr_flag(iv_e) == 0)
                                {
                                    idata.nbr_cellid_.push_back(iv_e);
                                    nbr_flag(iv_e) = 1;
                                    idata.n_nbr++;
                                }
                            }
                            
                            if (psi(i, j - 1, k) < 0.0 && isect.contains(iv_s))
                            {
                                Is_Intersect = true;
                                frac = psi(i, j, k) / (psi(i, j, k) - psi(i, j - 1, k));
                                if (CUBIC_SOLVE)
                                {
                                    frac = frac - 1.0;
                                    cu_ = Cubic_Solve(frac, psi(i, j - 2, k), psi(i, j - 1, k), psi(i, j, k), psi(i, j + 1, k));
                                    frac = frac + 1.0;
                                    if (cu_ == -1)
                                    {
                                        amrex::Print() << "\nY Intercept: " << frac << "\n";
                                        exit(1);
                                    }
                                }
                                // sanity check....
                                if (frac < 0.0 || frac > 1.0)
                                {
                                    amrex::Print() << "\nIntercept out of bounds!" << frac << "\n";
                                    amrex::Print() << "frac linear: j,j-1" << psi(i, j, k) / (psi(i, j, k) - psi(i, j - 1, k)) << "\n";
                                    amrex::Print() << psi(i, j - 2, k) << "\t" << psi(i, j - 1, k) << "\t" << psi(i, j, k) << "\t" << psi(i, j + 1, k) << "\n";
                                    exit(1);
                                }
                                // check for the distance....
                                if (i_ > 0)
                                {
                                    Retain = false;
                                    distance = sqrt(frac * frac + idata.frac_[0] * idata.frac_[0]); // even if i_ = 2, 1 and 0 will be sorted. therefore compare only with index 0.
                                    if (distance > dist_TOL)
                                        Retain = true;
                                }
                                else
                                    Retain = true;
                                
                                if (Retain)
                                {
                                    if ((i_ == 0) || (fabs(frac) > fabs(idata.frac_[0])))
                                    {
                                        idata.type_[i_] = 1; // 1 for y..
                                        idata.frac_[i_] = frac;
                                        Compute_Normal_Curvature(dPsi_L, Curv_L, i, j - 1, psi, dx[1],dx, prob_lo, isAxisymmetric);
                                        Compute_Normal_Curvature(dPsi_R, Curv_R, i, j, psi, dx[1],dx, prob_lo, isAxisymmetric);
                                        idata.psix_[i_] = frac * dPsi_L[0] + (1.0 - frac) * dPsi_R[0];
                                        idata.psiy_[i_] = frac * dPsi_L[1] + (1.0 - frac) * dPsi_R[1];
                                        idata.kappa_[i_] = frac * Curv_L + (1.0 - frac) * Curv_R;
                                    }
                                    else
                                    {
                                        if (i_ == 1)
                                        {
                                            idata.type_[i_] = idata.type_[i_ - 1];
                                            idata.frac_[i_] = idata.frac_[i_ - 1];
                                            idata.psix_[i_] = idata.psix_[i_ - 1];
                                            idata.psiy_[i_] = idata.psiy_[i_ - 1];
                                            idata.kappa_[i_] = idata.kappa_[i_ - 1];
                                        }
                                        else if (i_ == 2)
                                        {
                                            idata.type_[i_] = idata.type_[i_ - 1];
                                            idata.frac_[i_] = idata.frac_[i_ - 1];
                                            idata.psix_[i_] = idata.psix_[i_ - 1];
                                            idata.psiy_[i_] = idata.psiy_[i_ - 1];
                                            idata.kappa_[i_] = idata.kappa_[i_ - 1];
                                
                                            idata.type_[i_ - 1] = idata.type_[i_ - 2];
                                            idata.frac_[i_ - 1] = idata.frac_[i_ - 2];
                                            idata.psix_[i_ - 1] = idata.psix_[i_ - 2];
                                            idata.psiy_[i_ - 1] = idata.psiy_[i_ - 2];
                                            idata.kappa_[i_ - 1] = idata.kappa_[i_ - 2];
                                        }
                                        idata.type_[0] = 1; // 1 for y..
                                        idata.frac_[0] = frac;
                                        Compute_Normal_Curvature(dPsi_L, Curv_L, i, j - 1, psi, dx[1], dx, prob_lo, isAxisymmetric);
                                        Compute_Normal_Curvature(dPsi_R, Curv_R, i, j, psi, dx[1], dx, prob_lo, isAxisymmetric);
                                        idata.psix_[0] = frac * dPsi_L[0] + (1.0 - frac) * dPsi_R[0];
                                        idata.psiy_[0] = frac * dPsi_L[1] + (1.0 - frac) * dPsi_R[1];
                                        idata.kappa_[0] = frac * Curv_L + (1.0 - frac) * Curv_R;
                                    }
                                    i_++;
                                }
                                if(nbr_flag(iv_s) == 0)
                                {
                                    idata.nbr_cellid_.push_back(iv_s);
                                    nbr_flag(iv_s) = 1;
                                    idata.n_nbr++;
                                }
                            }
                            
                            if (psi(i, j + 1, k) < 0.0 && isect.contains(iv_n))
                            {
                                Is_Intersect = true;
                                frac = psi(i, j, k) / (psi(i, j + 1, k) - psi(i, j, k));
                                if (CUBIC_SOLVE)
                                {
                                    cu_ = Cubic_Solve(frac, psi(i, j - 1, k), psi(i, j, k), psi(i, j + 1, k), psi(i, j + 2, k));
                                    if (cu_ == -1)
                                    {
                                        amrex::Print() << "\nY Intercept: " << frac << "\n";
                                        exit(1);
                                    }
                                }
                                // sanity check....
                                if (frac > 0.0 || frac < -1.0)
                                {
                                    amrex::Print() << "\nIntercept out of bounds!" << frac << "\n";
                                    std::cout<<" i = "<<i<<" , j = "<<j<<" , k = "<<k<<" , process = "<<amrex::ParallelDescriptor::IOProcessorNumber()<<'\n';
                                    amrex::Print() << "frac linear: j,j+1" << psi(i, j, k) / (psi(i, j + 1, k) - psi(i, j, k)) << "\n";
                                    amrex::Print() << psi(i, j - 1, k) << "\t" << psi(i, j, k) << "\t" << psi(i, j + 1, k) << "\t" << psi(i, j + 2, k) << "\n";
                                    exit(1);
                                }
                                // check for the distance....
                                if (i_ > 0)
                                {
                                    Retain = false;
                                    if (idata.type_[0] == 1)
                                    {
                                        distance = fabs(frac - idata.frac_[0]);
                                    }
                                    else
                                    {
                                        distance = sqrt(frac * frac + idata.frac_[0] * idata.frac_[0]); // even if i_ = 2, 1 and 0 will be sorted. therefore compare only with index 0.
                                    }
                                    if (distance > dist_TOL)
                                        Retain = true;
                                }
                                else
                                    Retain = true;
                                if (Retain)
                                {
                                    if ((i_ == 0) || (i_ == 3) || (fabs(frac) > fabs(idata.frac_[0])))
                                    {
                                        idata.type_[i_] = 1; // 1 for y..
                                        idata.frac_[i_] = frac;
                                        Compute_Normal_Curvature(dPsi_L, Curv_L, i, j, psi, dx[1], dx, prob_lo, isAxisymmetric);
                                        Compute_Normal_Curvature(dPsi_R, Curv_R, i, j + 1, psi, dx[1], dx, prob_lo, isAxisymmetric);
                                        idata.psix_[i_] = -frac * dPsi_R[0] + (1.0 + frac) * dPsi_L[0];
                                        idata.psiy_[i_] = -frac * dPsi_R[1] + (1.0 + frac) * dPsi_L[1];
                                        idata.kappa_[i_] = -frac * Curv_R + (1.0 + frac) * Curv_L;
                                    }
                                    else
                                    {
                                        if (i_ == 1)
                                        {
                                            idata.type_[i_] = idata.type_[i_ - 1];
                                            idata.frac_[i_] = idata.frac_[i_ - 1];
                                            idata.psix_[i_] = idata.psix_[i_ - 1];
                                            idata.psiy_[i_] = idata.psiy_[i_ - 1];
                                            idata.kappa_[i_] = idata.kappa_[i_ - 1];
                                        }
                                        else if (i_ == 2)
                                        {
                                            idata.type_[i_] = idata.type_[i_ - 1];
                                            idata.frac_[i_] = idata.frac_[i_ - 1];
                                            idata.psix_[i_] = idata.psix_[i_ - 1];
                                            idata.psiy_[i_] = idata.psiy_[i_ - 1];
                                            idata.kappa_[i_] = idata.kappa_[i_ - 1];
                                
                                            idata.type_[i_ - 1] = idata.type_[i_ - 2];
                                            idata.frac_[i_ - 1] = idata.frac_[i_ - 2];
                                            idata.psix_[i_ - 1] = idata.psix_[i_ - 2];
                                            idata.psiy_[i_ - 1] = idata.psiy_[i_ - 2];
                                            idata.kappa_[i_ - 1] = idata.kappa_[i_ - 2];
                                        }
                                        idata.type_[0] = 1; // 1 for y..
                                        idata.frac_[0] = frac;
                                        Compute_Normal_Curvature(dPsi_L, Curv_L, i, j, psi, dx[1], dx, prob_lo, isAxisymmetric);
                                        Compute_Normal_Curvature(dPsi_R, Curv_R, i, j + 1, psi, dx[1], dx, prob_lo, isAxisymmetric);
                                        idata.psix_[0] = -frac * dPsi_R[0] + (1.0 + frac) * dPsi_L[0];
                                        idata.psiy_[0] = -frac * dPsi_R[1] + (1.0 + frac) * dPsi_L[1];
                                        idata.kappa_[0] = -frac * Curv_R + (1.0 + frac) * Curv_L;
                                    }
                                    i_++;
                                }
                                if(nbr_flag(iv_n) == 0)
                                {
                                    idata.nbr_cellid_.push_back(iv_n);
                                    nbr_flag(iv_n) = 1;
                                    idata.n_nbr++;
                                }
                            }
                        }
                        if (Is_Intersect)
                        {
                            pmask(i, j, k) = 2;
                            idata.cellid_ = amrex::IntVect(AMREX_D_DECL(i, j, k));
                            idata.n_intercepts = i_;
                            //indicator(i, j, k) = l;
                            N_Intercept++;
                            l++;
                            
                            /// add idata to icpt_data vect
                            icpt_data.push_back(std::move(idata));
                        }				    
                    }
                    pmask_temp(i, j, k) = pmask(i, j, k);
		    /*if(i == 95 && j == 15)
	            {
		        amrex::Print()<<" i = "<<i<<" , j = "<<j<<" , pmask_temp = "<<pmask_temp(i, j, k)<<'\n';
			amrex::Print()<<"psi(i, j, k) = "<<psi(i, j, k)<<'\n';
		    }*/
                });
                
                icpt_data.shrink_to_fit();
            }
            amrex::ParallelDescriptor::ReduceIntSum(N_Intercept);
        }
        
        
        MaskBC();
        for (int outer = 0; outer < LAYERS; outer++)
        {
            for (amrex::MFIter mfi(mask_); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();
                amrex::Array4<int> const &pmask = mask_.array(mfi);
                amrex::Array4<int> const &pmask_temp = maskTemp_.array(mfi);
                auto &Index_Additional_Layers = Index_Additional_Layers_LD[mfi];
        
                if (!Index_Additional_Layers[outer].empty())
                {
                    Index_Additional_Layers[outer].clear();
                    // Index_Additional_Layers[outer].shrink_to_fit();
                }
                else
                {
                    Index_Additional_Layers[outer].reserve(100);
                }
        
                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        
                    if (pmask_temp(i, j, k) == 0)
                    {
                        if (outer == 0)
                        {
                            if (pmask_temp(i + 1, j, k) == 1 || pmask_temp(i - 1, j, k) == 1 || pmask_temp(i, j + 1, k) == 1 || pmask_temp(i, j - 1, k) == 1)
                            {
                                amrex::Print() << "\n Detected a problem with mask_ \n"
                                                << "PMasks assigned one and zero are neighbors\n"
						<<"cell = "<<i<<" , "<<j<<" , "<<k<<"\n"
						<<"Outer = "<<outer<<"\n"
						<<pmask_temp(i + 1, j, k)<<" , "<<pmask_temp(i - 1, j, k)<<" , "<<pmask_temp(i, j + 1, k)<<" , "<<pmask_temp(i, j - 1, k)<<"\n";
                                exit(1);
                            }
                            else if (pmask_temp(i + 1, j, k) > 1 || pmask_temp(i - 1, j, k) > 1 || pmask_temp(i, j + 1, k) > 1 || pmask_temp(i, j - 1, k) > 1)
                            {
                                pmask(i, j, k) = 100;
                                amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                                Index_Additional_Layers[outer].push_back(std::move(iv));
                            }
                        }
                        else
                        {
                            if (pmask_temp(i + 1, j, k) == 99 + outer || pmask_temp(i - 1, j, k) == 99 + outer || pmask_temp(i, j + 1, k) == 99 + outer || pmask_temp(i, j - 1, k) == 99 + outer)
                            {
                                pmask(i, j, k) = 100 + outer;
                                amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                                Index_Additional_Layers[outer].push_back(std::move(iv));
                            }
                        }
                    }
                });
        
                Index_Additional_Layers[outer].shrink_to_fit();
            }
            amrex::iMultiFab::Copy(maskTemp_, mask_, 0, 0, 1, 0);
            MaskBC();
        }
        amrex::iMultiFab::Copy(PMask, mask_, 0, 0, 1, nghost);
	PMask.FillBoundary();

        /// UMask
        UMask.setVal(1);
        for(amrex::MFIter mfi(UMask); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<int> const &umask = UMask.array(mfi);
            amrex::Array4<int const> const &pmask = PMask.const_array(mfi);
        
            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (pmask(i, j, k) != 1 && pmask(i - 1, j, k) != 1)
                {
                    umask(i, j, k) = 0;
                }
            });
        }
        UMask.FillBoundary();
        
        /// VMask
        VMask.setVal(1);
        for(amrex::MFIter mfi(VMask); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<int> const &vmask = VMask.array(mfi);
            amrex::Array4<int const> const &pmask = PMask.const_array(mfi);
        
            amrex::ParallelFor(bx,
            [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (pmask(i, j, k) != 1 && pmask(i, j - 1, k) != 1)
                {
                    vmask(i, j, k) = 0;
                }
            });
        }
        VMask.FillBoundary();
    }
    
    void Mask::SetInterfacePressureVel()
    {
        const amrex::Real *prob_lo = geom_.ProbLo();
        const amrex::Real *dx = geom_.CellSize();
        amrex::Real Polytropic_Index = 1.4;//hard coded!!
        amrex::Real mu_ = Mu;
        int p = 0;
        for (auto &&solid : *interfaces)
        {
            if (solid->isAdvectLS())
            {
                //solid->ComputeVolume(false);
                //solid->ComputeAvgIntVel(u,v,PMask);
                //amrex::Print()<<"solid->Volume() = "<<solid->Volume()<<'\n';
                amrex::Real resolution_threshold = 100.0;
                amrex::Real bubble_resolution = solid->getResolution();
                amrex::Real d_gamma = 0.0;//modified Polytropic index in under resolved bubble is increased by 2.0
                amrex::Real Polytropic_Index_modified = Polytropic_Index;
                if(bubble_resolution < 400.0) 
                {
		    Polytropic_Index_modified = Polytropic_Index + 0.5*d_gamma*
                                            (1.0 - std::tanh((bubble_resolution - resolution_threshold)/resolution_threshold));
                    //amrex::Print()<<"bubble is under resolved! bubble resolution = "<<bubble_resolution<<'\n';
                    //amrex::Print()<<"Polytropic_Index_modified = "<<Polytropic_Index_modified<<'\n'; 
                }
                if (IF_types[p] == "Bubble" || IF_types[p] == "bubble")
                {
                    amrex::Real P_int = solid->getP_interface0() * pow((solid->Volume_0() / solid->Volume()), Polytropic_Index_modified);
                    solid->setP_interface(P_int);
                }
                else
                {
    		    amrex::Real P_int = solid->getP_interface0();
                    solid->setP_interface(P_int); 
                }
		
            }
            //amrex::Print()<<" P_interface  = "<<solid->P_interface<<'\n';
            //amrex::Print()<<" P_interface0  = "<<solid->P_interface0<<'\n';
            //amrex::Print()<<"Volume  = "<<solid->Volume()<<'\n';
            //Switch for aoustic radiation force(ARF). ARF if turned off 
            //when bubble is trying to penetrate into the tissue
            amrex::Real apply_ARF = 1.0;
            //for (amrex::MFIter mfi(PMask); mfi.isValid(); ++mfi)
            //{
            //    auto &icpt_data = solid->getInterceptData()[mfi];
            //    for (auto &&idt : icpt_data)
            //    {
            //        for (int ni = 0; ni < idt.n_intercepts; ni++)
            //        {
            //            amrex::Real x = prob_lo[0] + dx[0] * (idt.cellid_[0] + 0.5);
            //            amrex::Real y = prob_lo[1] + dx[1] * (idt.cellid_[1] + 0.5);
            //            if (idt.type_[ni] == 0)
            //            {
            //                x += -idt.frac_[ni] * dx[0];
            //            }
            //            else if (idt.type_[ni] == 1)
            //            {
            //                y += -idt.frac_[ni] * dx[1];
            //            }
            //            amrex::Real r = std::hypot(x - solid->Xcp(), y - solid->Ycp());
            //            amrex::Real theta = std::acos(x/r);
                     
            //            //amrex::Print()<<"x = "<<x<<", y = "<<y<<" , theta = "<<theta<<"idt.phi_Int[ni] = "<<idt.phi_Int[ni]<<'\n';
                        
            //            //if(PhaseField && idt.phi_Int[ni] > 0.0 && theta >= 0.5*M_PI)
            //            //    apply_ARF = 0.0;// turn off ARF
            //        }
            //     }
            //}
            apply_ARF = 1.0;
            amrex::ParallelDescriptor::ReduceRealMin(apply_ARF);
            //amrex::Print()<<"Apply ARF = "<<apply_ARF<<'\n'; 

            for (amrex::MFIter mfi(PMask); mfi.isValid(); ++mfi)
            {
                auto &icpt_data = solid->getInterceptData()[mfi];
                
    
                for (auto &&idt : icpt_data)
                {
                    for (int ni = 0; ni < idt.n_intercepts; ni++)
                    {
                        amrex::Real x = prob_lo[0] + dx[0] * (idt.cellid_[0] + 0.5);
                        amrex::Real y = prob_lo[1] + dx[1] * (idt.cellid_[1] + 0.5);
                        if (idt.type_[ni] == 0)
                        {
                            x += -idt.frac_[ni] * dx[0];
                        }
                        else if (idt.type_[ni] == 1)
                        {
                            y += -idt.frac_[ni] * dx[1];
                        }
		   
			idt.Dmg_Int[ni] = 0.0;
                        /*
                        if(DamageModel)
                        {
                            amrex::Real ep_max;
                            amrex::Real C11 = idt.F11_Int[ni]*idt.F11_Int[ni] + idt.F21_Int[ni]*idt.F21_Int[ni];
                            amrex::Real C12 = idt.F11_Int[ni]*idt.F12_Int[ni] + idt.F21_Int[ni]*idt.F22_Int[ni];
                            amrex::Real C21 = idt.F12_Int[ni]*idt.F11_Int[ni] + idt.F22_Int[ni]*idt.F21_Ini[ni];                           
			    amrex::Real C22 = idt.F12_Int[ni]*idt.F12_Int[ni] + idt.F22_Int[ni]*idt.F22_Int[ni];
                            amrex::Real C33 = idt.F33_Int[ni]*idt.F33_Int[ni];

                            amrex::Real first_inv = C11+C22+C33;
                            amrex::Real second_inv = C22*C33 + C11*C33 + C11*C22 - C12*C21;
                            amrex::Real third_inv = C33*(C11*C22 - C12*C21);

                            amrex::Real lambda_1;// = 0.5*(tr_C + std::sqrt(tr_C*tr_C - 4.0*det_C));
                            amrex::Real lambda_2;// = 0.5*(tr_C - std::sqrt(tr_C*tr_C - 4.0*det_C));
                            amrex::Real lambda_3;
                            SolveCubicEqn(lambda_1, lambda_2, lambda_3, first_inv, second_inv, third_inv);
                            if(std::isnan(lambda_1)) lambda_1 = 1.0;
                            if(std::isnan(lambda_2)) lambda_2 = 1.0;
                            if(std::isnan(lambda_3)) lambda_3 = 1.0;
                            ep_max = std::max(std::abs(0.5*(lambda_1 - 1.0)),
                                               std::abs(0.5*(lambda_2 - 1.0)));
                            ep_max = std::max(ep_max, std::abs(0.5*(lambda_3 - 1.0)));

                            idt.emax_Int[ni] = ep_max;

                            if(ep_max > 2.0)
                                idt.Dmg_Int[ni] = 1.0;

                        }
                        */
			if(PhaseField) mu_ = visc_.GetViscosity(idt.Gamma_dot_Int[ni],idt.phi_Int[ni]);
			//amrex::Print()<<"mu_ = "<<mu_<<'\n';
                        if (solid->isAdvectLS())
                        {
                            if (IF_types[p] == "Bubble" || IF_types[p] == "bubble")
                            {
                                idt.P[ni] = solid->P_interface + SIGMA * idt.kappa_[ni] +   mu_ * idt.norm_shear_[ni];
                                amrex::Real dpdz = 0.0;// pressure at the interface due to acoustic radiation force
                                if(solid->bdy_frc_type == solid->const_pressure_gradient)
                                    dpdz = solid->bdy_frc;
                                else if(solid->bdy_frc_type == solid->const_force)
                                    dpdz = solid->bdy_frc/solid->Volume();
                                idt.P[ni] += apply_ARF * dpdz * (solid->Xcp() - x);
                                idt.T_Int[ni] = solid->P_interface * solid->Volume() / (solid->P_interface0 * solid->Volume_0());
				//amrex::Print()<<"dpdz = "<<dpdz<<'\n';
                                //if(idt.cellid_[0] == 202 && idt.cellid_[1] == 45)
                                //{
                                //    std::cout<<" Int_P  "<<ni<<"  = "<<idt.P[ni]<<" , norm_sh = "<<idt.norm_shear_[ni]<<"idt.kappa_[ni] = "<<idt.kappa_[ni]<<'\n';
                                //}
				//amrex::Print()<<"Location = "<<
				//amrex::Print()<<"idt.psix and psiy["<<ni<<"] = "<<idt.psix_[i_]<<" , "<<idt.psix_[i_]<<'\n';
				//amrex::Print()<<"idt.kappa_["<<ni<<"] = "<<idt.kappa_[ni]<<'\n';
                                //amrex::Print()<<"idt.norm_shear_["<<ni<<"] = "<<idt.norm_shear_[ni]<<'\n';


                            }
                            else
                            {
                                idt.P[ni] = solid->P_interface;
                                idt.T_Int[ni] = solid->P_interface * solid->Volume() / (solid->P_interface0 * solid->Volume_0());
                            }
                            //idt.T_Int[ni] = 1.0;
                            //idt.P[ni] = 20.0;
                        }
                        else
                        {
                            idt.P[ni] = 0.0; /// neumann
                            idt.u[ni] = 0.0; ///solid->ddt_xcp() - (y - solid->Ycp()) * solid->ddt_thetacp();
                            idt.v[ni] = 0.0; ///solid->ddt_ycp() + (x - solid->Xcp()) * solid->ddt_thetacp();
                        }   

                        //amrex::Real npi = 2.0*M_PI;
                        //amrex::Real npi2 = npi * npi;
                        //idt.P[ni] = sin(npi * x) * sin(npi * y);
                        //amrex::Print()<<"idt.P[ni]  = "<<idt.P[ni] <<'\n';                 
                    }
                }
            }
            //if (solid->isAdvectLS()) solid->ComputeAvgIntVel(u,v,PMask);
            p++;
        }
    }
    
    void Mask::FillInGhost(amrex::MultiFab &mf, const std::vector<amrex::LinOpBCType> &intbcs)
    {
        if (interfaces->empty())
        {
            return;
        }
    
        /// fill internal ghost values for pressure
        mf.FillBoundary();
    
        const amrex::Real *prob_lo = geom_.ProbLo();
        const amrex::Real *dx = geom_.CellSize();
        const amrex::Box &domain = geom_.Domain();
    
        amrex::Real sum_w, sum_sol;
        bool Sing_;
    
        int nsolid = 0;
        for (auto &&solid : *interfaces)
        {
            for (amrex::MFIter mfi(mask_); mfi.isValid(); ++mfi)
            {
                amrex::Array4<int const> const &pmask = PMask.const_array(mfi);
                amrex::Array4<amrex::Real> const &P = mf.array(mfi);
                auto &icpt_data = solid->getInterceptData()[mfi];
              
                for (auto &&idt : icpt_data)
                {
                    int dist_x1 = abs(idt.cellid_[0] - domain.smallEnd(0));
                    int dist_x2 = abs(idt.cellid_[0] - domain.bigEnd(0));
                    int dist_y1 = abs(idt.cellid_[1] - domain.smallEnd(1));
                    int dist_y2 = abs(idt.cellid_[1] - domain.bigEnd(1));
                    
                    //decide size of the lease-square stencil based on available points
                    if( dist_x1 >= 4 && dist_x2 >= 4 && dist_y1 >= 4 && dist_y2 >= 4)
                    {
                            idt.stencil_ = 4;
                    }
                    else if( dist_x1 >= 3 && dist_x2 >= 3 && dist_y1 >= 3 && dist_y2 >= 3)
                    {
                            idt.stencil_ = 3;
                            idt.PORDER = 3;
                            //amrex::Print()<<"setting stencil size to 3 "<<idt.cellid_<<'\n';
                    }
                    else if( dist_x1 >= 2 && dist_x2 >= 2 && dist_y1 >= 2 && dist_y2 >= 2)
                    {
                            idt.stencil_ = 2;
                            idt.PORDER = 2;
                            //amrex::Print()<<"setting stencil size to 2 "<<idt.cellid_<<'\n';
                    }
                    else
                    {
                            idt.stencil_ = 2;
                            idt.PORDER = 2;
                            //amrex::Print()<<"setting stencil size to 1 "<<idt.cellid_<<'\n';
                    }
                    //idt.stencil_ = 4;
                    //idt.PORDER = 3;

                    int N_Max = (2 * idt.stencil_ + 1) * (2 * idt.stencil_ + 1) + 4;
                    double x_[N_Max], y_[N_Max], w_[N_Max], P_[N_Max];
                    Sing_ = false; // Need this. Sing_ is only changed to true.
                    int ni;
                    for (ni = 0; ni < idt.n_intercepts; ni++)
                    {
                        x_[ni] = y_[ni] = 0.0;
                        if (idt.type_[ni] == 0)
                        {
                            x_[ni] = -idt.frac_[ni];
                        }
                        else if (idt.type_[ni] == 1)
                        {
                            y_[ni] = -idt.frac_[ni];
                        }
                        else
                        {
                            amrex::Print() << "Error in identifying interface type.\n";
                            exit(1);
                        }
                        P_[ni] = idt.P[ni];
                        Compute_LSQ_Weights(w_[ni], std::fabs(idt.frac_[ni]));
                    }
                    int Num_i_ = ni;
    
                    /// create grown box around cell
                    amrex::Box gbx(idt.cellid_, idt.cellid_);
                    gbx.grow(idt.stencil_);
                    
                    amrex::Box gbx_isect = gbx & domain; 
    
                    for (amrex::BoxIterator bit(gbx_isect); bit.ok(); ++bit)
                    {
                        const amrex::IntVect &iv = bit();
                        if (pmask(iv) == 1)
                        {
                            LSQ_Parameters(P_[ni], x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], P, prob_lo, dx);
                            ni++;
                        }
                    }
    
                    int Num_ = ni;
                    idt.Interpolation_weights.resize(Num_);
    
                    if (idt.PORDER == 4)
                    {
                        if (Num_ < 3 + Num_i_)
                        {
                            /// too few points have been found...can only do a simple average
                            sum_w = sum_sol = 0.0;
                            for (int i = 0; i < Num_; i++)
                            {
                                sum_w += w_[i];
                                sum_sol += w_[i] * P_[i];
                            }
                            P(idt.cellid_) = sum_sol / sum_w;
    
                            for (int i = 0; i < Num_; i++)
                            {
                                idt.Interpolation_weights[i] = w_[i] / sum_w;
                            }
                            //amrex::Print() << " Inside first ";
                        }
                        else if (Num_ < 6 + Num_i_)
                        {
                            QR_LS_Pressure_Second_Order(Sing_, Num_, Num_i_, x_, y_, P_, w_, P, idt);
			    //if(P(idt.cellid_) < 0.0) Sing_ = true;
                            if (Sing_)
                            {
                                sum_w = sum_sol = 0.0;
                                for (int i = 0; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol += w_[i] * P_[i];
                                }
                                P(idt.cellid_) = sum_sol / sum_w;
    
                                for (int i = 0; i < Num_; i++)
                                {
                                    idt.Interpolation_weights[i] = w_[i] / sum_w;
                                }
                            }
                            //amrex::Print() << " Second first ";
                        }
                        else if (Num_ < 10 + Num_i_)
                        {
                            QR_LS_Pressure_Third_Order(Sing_, Num_, Num_i_, x_, y_, P_, w_, P, idt);
                            if (Sing_)
                            {
                                QR_LS_Pressure_Second_Order(Sing_, Num_, Num_i_, x_, y_, P_, w_, P, idt);
                            }
                            if (Sing_)
                            {
                                sum_w = sum_sol = 0.0;
                                for (int i = 0; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol += w_[i] * P_[i];
                                }
                                P(idt.cellid_) = sum_sol / sum_w;
    
                                for (int i = 0; i < Num_; i++)
                                {
                                    idt.Interpolation_weights[i] = w_[i] / sum_w;
                                }
                            }
                            //amrex::Print() << " Third first ";
                        }
                        else
                        {
                            AllMCQR_LS_Pressure_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, P_, w_, P, idt, intbcs[nsolid]);
			    //if(P(idt.cellid_) < 0.0) Sing_ = true;
                            if (Sing_)
                            {
                                QR_LS_Pressure_Third_Order(Sing_, Num_, Num_i_, x_, y_, P_, w_, P, idt);
                                //amrex::Print() << " Inside second ";
				//if(P(idt.cellid_) < 0.0) Sing_ = true;
                            }
                            if (Sing_)
                            {
                                QR_LS_Pressure_Second_Order(Sing_, Num_, Num_i_, x_, y_, P_, w_, P, idt);
                                //amrex::Print() << " Inside first ";
				//if(P(idt.cellid_) < 0.0) Sing_ = true;
                            }
                            if (Sing_)
                            {
                                sum_w = sum_sol = 0.0;
                                for (int i = 0; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol += w_[i] * P_[i];
                                }
                                P(idt.cellid_) = sum_sol / sum_w;
    
                                for (int i = 0; i < Num_; i++)
                                {
                                    idt.Interpolation_weights[i] = w_[i] / sum_w;
                                }
                                //amrex::Print() << " Inside zero ";
                            }
                        }
                    }
                    else if (idt.PORDER == 3)
                    {
                        if (Num_ < 3 + Num_i_)
                        {
                            /// too few points have been found...can only do a simple average
                            sum_w = sum_sol = 0.0;
                            for (int i = 0; i < Num_; i++)
                            {
                                sum_w += w_[i];
                                sum_sol += w_[i] * P_[i];
                            }
                            P(idt.cellid_) = sum_sol / sum_w;
    
                            for (int i = 0; i < Num_; i++)
                            {
                                idt.Interpolation_weights[i] = w_[i] / sum_w;
                            }
                            //amrex::Print() << " Inside first ";
                        }
                        else if (Num_ < 6 + Num_i_)
                        {
                            QR_LS_Pressure_Second_Order(Sing_, Num_, Num_i_, x_, y_, P_, w_, P, idt);
                            if (Sing_)
                            {
                                sum_w = sum_sol = 0.0;
                                for (int i = 0; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol += w_[i] * P_[i];
                                }
                                P(idt.cellid_) = sum_sol / sum_w;
    
                                for (int i = 0; i < Num_; i++)
                                {
                                    idt.Interpolation_weights[i] = w_[i] / sum_w;
                                }
                            }
                            //amrex::Print() << " Second first ";
                        }
                        else
                        {
                            QR_LS_Pressure_Third_Order(Sing_, Num_, Num_i_, x_, y_, P_, w_, P, idt);
                            if (Sing_)
                            {
                                QR_LS_Pressure_Second_Order(Sing_, Num_, Num_i_, x_, y_, P_, w_, P, idt);
                            }
                            if (Sing_)
                            {
                                sum_w = sum_sol = 0.0;
                                for (int i = 0; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol += w_[i] * P_[i];
                                }
                                P(idt.cellid_) = sum_sol / sum_w;
    
                                for (int i = 0; i < Num_; i++)
                                {
                                    idt.Interpolation_weights[i] = w_[i] / sum_w;
                                }
                            }
                            //amrex::Print() << " Third first ";
                        }
                    }
                    else if (idt.PORDER == 2)
                    {
                        if (Num_ < 3 + Num_i_)//ideally should be 3
                        {
                            /// too few points have been found...can only do a simple average
                            sum_w = sum_sol = 0.0;
                            for (int i = 0; i < Num_; i++)
                            {
                                sum_w += w_[i];
                                sum_sol += w_[i] * P_[i];
                            }
                            P(idt.cellid_) = sum_sol / sum_w;
    
                            for (int i = 0; i < Num_; i++)
                            {
                                idt.Interpolation_weights[i] = w_[i] / sum_w;
                            }
                            //amrex::Print() << " Inside first at "<<idt.cellid_[0]<<" , "<<idt.cellid_[1]<<'\n';
                        }
                        else
                        {
                            QR_LS_Pressure_Second_Order(Sing_, Num_, Num_i_, x_, y_, P_, w_, P, idt);
                            if (Sing_)
                            {
                                sum_w = sum_sol = 0.0;
                                for (int i = 0; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol += w_[i] * P_[i];
                                }
                                P(idt.cellid_) = sum_sol / sum_w;

    
                                for (int i = 0; i < Num_; i++)
                                {
                                    idt.Interpolation_weights[i] = w_[i] / sum_w;
                                }
                            }
                            //amrex::Print() << " Second first at "<<idt.cellid_[0]<<" , "<<idt.cellid_[1]<<'\n';
                        }
                    }
		    
                    //if(idt.cellid_[0] == 1014 && idt.cellid_[1] == 7)
                    //{
                    //    amrex::Print(-1)<<"i = "<<idt.cellid_[0]<<" , j = "<<idt.cellid_[1]<<'\n';
                    //    amrex::Print(-1)<<"idt.stencil_ = "<<idt.stencil_<<'\n';
                    //    amrex::Print(-1)<<"idt.PORDER = "<<idt.PORDER<<"; Num_ = "<<Num_<<" , Num_i_ = "<<Num_i_<<'\n';
                    //    amrex::Print(-1)<<"P(idt.cellid_) = "<<P(idt.cellid_)<<'\n';
                        //if(P(idt.cellid_) < 0) exit(7);
                    //}
                }
            }
            nsolid++;
        }
    
        /// as values are updated
        /// fill the internal ghost values
        mf.FillBoundary();
    }

    void Mask::FillInGhostTheta(amrex::MultiFab &mf,const std::vector<amrex::LinOpBCType> &intbcs)
    {
        /// fill internal ghost values for pressure
        mf.FillBoundary();
    
        const amrex::Real *prob_lo = geom_.ProbLo();
        const amrex::Real *dx = geom_.CellSize();
        const amrex::Box &domain = geom_.Domain();
    
        amrex::Real sum_w, sum_sol;
        bool Sing_;
    
        int nsolid = 0;
        for (auto &&solid : *interfaces)
        {
            for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                amrex::Array4<int const> const &pmask = PMask.const_array(mfi);
                amrex::Array4<amrex::Real> const &T = mf.array(mfi);
                auto &icpt_data = solid->getInterceptData()[mfi];
    
                for (auto &&idt : icpt_data)
                {
                    int N_Max = (2 * idt.stencil_ + 1) * (2 * idt.stencil_ + 1) + 4;
                    double x_[N_Max], y_[N_Max], w_[N_Max], T_[N_Max];
                    Sing_ = false; // Need this. Sing_ is only changed to true.
                    int ni;
                    for (ni = 0; ni < idt.n_intercepts; ni++)
                    {
                        x_[ni] = y_[ni] = 0.0;
                        if (idt.type_[ni] == 0)
                        {
                            x_[ni] = -idt.frac_[ni];
                        }
                        else if (idt.type_[ni] == 1)
                        {
                            y_[ni] = -idt.frac_[ni];
                        }
                        else
                        {
                            amrex::PrintToFile("log") << "Error in identifying interface type.\n";
                            exit(1);
                        }
                        T_[ni] = idt.T_Int[ni];
                        Compute_LSQ_Weights(w_[ni], std::fabs(idt.frac_[ni]));
                    }
                    int Num_i_ = ni;
    
                    /// create grown box around cell
                    amrex::Box gbx(idt.cellid_, idt.cellid_);
                    gbx.grow(idt.stencil_);
                    amrex::Box gbx_isect = gbx & domain;
    
                    for (amrex::BoxIterator bit(gbx_isect); bit.ok(); ++bit)
                    {
                        const amrex::IntVect &iv = bit();
                        if (pmask(iv) == 1)
                        {
                            LSQ_Parameters(T_[ni], x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], T, prob_lo, dx);
                            ni++;
                        }
                    }
    
                    int Num_ = ni;
    
                    if (idt.PORDER >= 3)
                    {
                        if (Num_ < 3 + Num_i_)
                        {
                            /// too few points have been found...can only do a simple average
                            sum_w = sum_sol = 0.0;
                            for (int i = 0; i < Num_; i++)
                            {
                                sum_w += w_[i];
                                sum_sol += w_[i] * T_[i];
                            }
                            T(idt.cellid_) = sum_sol / sum_w;
    
                        }
                        else if (Num_ < 10 + Num_i_)
                        {
                            QR_LS_Temperature_Second_Order(Sing_, Num_, Num_i_, x_, y_, T_, w_, T, idt);
                            if (Sing_)
                            {
                                sum_w = sum_sol = 0.0;
                                for (int i = 0; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol += w_[i] * T_[i];
                                }
                                T(idt.cellid_) = sum_sol / sum_w;
    
                            }
                            //amrex::PrintToFile("log") << " Second first ";
                        }
                        else
                        {
                            QR_LS_Temperature_Third_Order(Sing_, Num_, Num_i_, x_, y_, T_, w_, T, idt);
                            if (Sing_)
                            {
                                QR_LS_Temperature_Second_Order(Sing_, Num_, Num_i_, x_, y_, T_, w_, T, idt);
                            }
                            if (Sing_)
                            {
                                sum_w = sum_sol = 0.0;
                                for (int i = 0; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol += w_[i] * T_[i];
                                }
                                T(idt.cellid_) = sum_sol / sum_w;
                            }
                            //amrex::PrintToFile("log") << " Third first ";
                        }
                    }
    
                    else if (idt.PORDER == 2)
                    {
                        if (Num_ < 6 + Num_i_)
                        {
                            /// too few points have been found...can only do a simple average
                            sum_w = sum_sol = 0.0;
                            for (int i = 0; i < Num_; i++)
                            {
                                sum_w += w_[i];
                                sum_sol += w_[i] * T_[i];
                            }
                            T(idt.cellid_) = sum_sol / sum_w;
                        }
                        else
                        {
                            QR_LS_Temperature_Second_Order(Sing_, Num_, Num_i_, x_, y_, T_, w_, T, idt);
                            if (Sing_)
                            {
                                sum_w = sum_sol = 0.0;
                                for (int i = 0; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol += w_[i] * T_[i];
                                }
                                T(idt.cellid_) = sum_sol / sum_w;
                            }
                            //amrex::PrintToFile("log") << " Second first ";
                        }
                    }
    		/*if(idt.cellid_[0] == 2 && idt.cellid_[1] == 138 && Iter == 1730 && debug)
    		{
    			amrex::Print()<<"cell : "<<idt.cellid_<<'\n';
                            amrex::Print()<<"temperature : "<< T(idt.cellid_)<<'\n';
                            amrex::Print()<<"Num_ = "<<Num_<<" , Num_i_ = "<<Num_i_<<'\n';
                            for(int ii = 0;ii<= Num_;ii++)
                            {
                                    amrex::Print()<<"ii = "<<ii<<"\t";
                            	amrex::Print()<<"x_ = "<<x_[ii]<<"\t";
                            	amrex::Print()<<"y_ = "<<y_[ii]<<"\t";
                            	amrex::Print()<<"T_ = "<<T_[ii]<<"\t";
                                    amrex::Print()<<"idt.Interpolation_weights = "<<idt.Interpolation_weights[ii]<<'\t';
                            	amrex::Print()<<"w_ = "<<w_[ii]<<"\n";
                            }
                            std::exit(9);
    		}*/
                }
            }
            nsolid++;
        }
    
        /// as pressure values are updated
        /// fill the internal ghost values
        mf.FillBoundary();
    }
    
    /*
    void Mask::FillInGhostPhaseField(amrex::MultiFab &mf, const std::vector<amrex::LinOpBCType> &intbcs)
    {
        /// fill internal ghost values for pressure
        mf.FillBoundary();
    
        const amrex::Real *prob_lo = geom_.ProbLo();
        const amrex::Real *dx = geom_.CellSize();
        const amrex::Box &domain = geom_.Domain();
    
        amrex::Real sum_w, sum_sol;
        bool Sing_;
    
        int nsolid = 0;
        for (auto &&solid : *interfaces)
        {
            for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                amrex::Array4<int const> const &pmask = PMask.const_array(mfi);
                amrex::Array4<amrex::Real> const &Phi = mf.array(mfi);
                auto &icpt_data = solid->getInterceptData()[mfi];
    
           
    
                for (auto &&idt : icpt_data)
                {
                    int N_Max = (2 * idt.stencil_ + 1) * (2 * idt.stencil_ + 1) + 4;
                    double x_[N_Max], y_[N_Max], w_[N_Max], Phi_[N_Max];
                    Sing_ = false; // Need this. Sing_ is only changed to true.
                    int ni;
                    for (ni = 0; ni < idt.n_intercepts; ni++)
                    {
                        x_[ni] = y_[ni] = 0.0;
                        if (idt.type_[ni] == 0)
                        {
                            x_[ni] = -idt.frac_[ni];
                        }
                        else if (idt.type_[ni] == 1)
                        {
                            y_[ni] = -idt.frac_[ni];
                        }
                        else
                        {
                            amrex::PrintToFile("log") << "Error in identifying interface type.\n";
                            exit(1);
                        }
                        Phi_[ni] = 0.0;//idt.T_Int[ni];
                        Compute_LSQ_Weights(w_[ni], std::fabs(idt.frac_[ni]));
                    }
                    int Num_i_ = ni;
    
                    /// create grown box around cell
                    amrex::Box gbx(idt.cellid_, idt.cellid_);
                    gbx.grow(idt.stencil_);
                    amrex::Box gbx_isect = gbx & domain;
    
                    for (amrex::BoxIterator bit(gbx_isect); bit.ok(); ++bit)
                    {
                        const amrex::IntVect &iv = bit();
                        if (pmask(iv) == 1)
                        {
                            LSQ_Parameters(Phi_[ni], x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], Phi, prob_lo, dx);
                            ni++;
                        }
                    }
    
                    int Num_ = ni;
    
		    int PORDER_ = 2;
                    if (PORDER_ >= 3)
                    {
                        if (Num_ < 3 + Num_i_)
                        {
                            /// too few points have been found...can only do a simple average
                            sum_w = sum_sol = 0.0;
                            for (int i = 0; i < Num_; i++)
                            {
                                sum_w += w_[i];
                                sum_sol += w_[i] * Phi_[i];
                            }
                            Phi(idt.cellid_) = sum_sol / sum_w;
    
                        }
                        else if (Num_ < 10 + Num_i_)
                        {
                            QR_LS_Temperature_Second_Order(Sing_, Num_, Num_i_, x_, y_, Phi_, w_, Phi, idt);
                            if (Sing_)
                            {
                                sum_w = sum_sol = 0.0;
                                for (int i = 0; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol += w_[i] * Phi_[i];
                                }
                                Phi(idt.cellid_) = sum_sol / sum_w;
    
                            }
                            //amrex::PrintToFile("log") << " Second first ";
                        }
                        else
                        {
                            QR_LS_Temperature_Third_Order(Sing_, Num_, Num_i_, x_, y_, Phi_, w_, Phi, idt);
                            if (Sing_)
                            {
                                QR_LS_Temperature_Second_Order(Sing_, Num_, Num_i_, x_, y_, Phi_, w_, Phi, idt);
                            }
                            if (Sing_)
                            {
                                sum_w = sum_sol = 0.0;
                                for (int i = 0; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol += w_[i] * Phi_[i];
                                }
                                Phi(idt.cellid_) = sum_sol / sum_w;
                            }
                            //amrex::PrintToFile("log") << " Third first ";
                        }
                    }
    
                    else if (PORDER_ == 2)
                    {
                        if (Num_ < 3 + Num_i_)
                        {
                            /// too few points have been found...can only do a simple average
                            sum_w = sum_sol = 0.0;
                            for (int i = 0; i < Num_; i++)
                            {
                                sum_w += w_[i];
                                sum_sol += w_[i] * Phi_[i];
                            }
                            Phi(idt.cellid_) = sum_sol / sum_w;
                            for (ni = 0; ni < idt.n_intercepts; ni++)
                            {
                                idt.phi_Int[ni] = sum_sol / sum_w;
                            }
                        }
                        else
                        {
                            QR_LS_Temperature_Second_Order(Sing_, Num_, Num_i_, x_, y_, T_, w_, T, idt);
                            if (Sing_)
                            {
                                sum_w = sum_sol = 0.0;
                                for (int i = Num_i_; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol += w_[i] * Phi_[i];
                                }
                                Phi(idt.cellid_) = sum_sol / sum_w;
				for (ni = 0; ni < idt.n_intercepts; ni++)
				{
				    idt.phi_Int[ni] = sum_sol / sum_w;
				}
                            }
                            //amrex::PrintToFile("log") << " Second first ";
                        }
                    }
		     
    		    if(idt.cellid_[0] == 115 && idt.cellid_[1] == 128)
    		    {
    			amrex::Print()<<"cell : "<<idt.cellid_<<'\n';
                            amrex::Print()<<"temperature : "<< T(idt.cellid_)<<'\n';
                            amrex::Print()<<"Num_ = "<<Num_<<" , Num_i_ = "<<Num_i_<<'\n';
                            for(int ii = 0;ii<= Num_;ii++)
                            {
                                    //amrex::Print()<<"ii = "<<ii<<"\t";
                            	    //amrex::Print()<<"x_ = "<<x_[ii]<<"\t";
                            	    //amrex::Print()<<"y_ = "<<y_[ii]<<"\t";
                            	    //amrex::Print()<<"T_ = "<<T_[ii]<<"\n";
                            }
                            //std::exit(9);
    		    }
                }
            }
            nsolid++;
        }
    
        /// as pressure values are updated
        /// fill the internal ghost values
        mf.FillBoundary();
    }            
    */
    
    void Mask::getFaceMask(amrex::iMultiFab &fmask, int dir)
    {
        fmask.setVal(1);
        if (interfaces->empty())
        {
            return;
        }
    
        for (amrex::MFIter mfi(fmask); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            amrex::Array4<int> const &fm = fmask.array(mfi);
            amrex::Array4<int const> const &pmask = mask_.const_array(mfi);
    
            if (dir == 0)
            {
                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (pmask(i, j, k) != 1 && pmask(i - 1, j, k) != 1)
                    {
                        fm(i, j, k) = 0;
                    }
                });
            }
            else if (dir == 1)
            {
                amrex::ParallelFor(bx,
                [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (pmask(i, j, k) != 1 && pmask(i, j - 1, k) != 1)
                    {
                        fm(i, j, k) = 0;
                    }
                });
            }        
        }
        fmask.FillBoundary();
    }

} /*End namespace mycode */

