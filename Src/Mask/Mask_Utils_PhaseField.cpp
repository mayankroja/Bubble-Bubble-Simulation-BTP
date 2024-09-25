#include "Mask.H"
#include <AMReX_ParmParse.H>
#include <qr.h>
#include <WeightedENO.h>
#include <Viscosity.H>

namespace mycode
{

    namespace AdvectLSIntPhaseField
    {    
        void QR_LS_PhaseField_Second_Order
        (
            bool& sing,
            int Num_, 
            double* x_, double* y_,
            double* Phi_,
            double* w_,
            amrex::Array4<amrex::Real> const& Phi,
	    const amrex::IntVect& cellid_
        )
        {
            int i_, j_ ;  
            double **Mat, *Sol, *Vec_  ;
        
            Allocate_2D_R(Mat,Num_,3) ; Sol = new double[3] ;
	    Vec_ = new double[Num_] ;
        
            for(i_ = 0; i_ < Num_ ; i_++)
            {
		int j_ = i_ ;
                Mat[j_][0] = w_[i_] ; // 1 
                Mat[j_][1] = w_[i_]*x_[i_] ; // x
                Mat[j_][2] = w_[i_]*y_[i_] ; // y
                Vec_[j_] = w_[i_]*Phi_[i_] ; // P
            } 
            QRdcmp<double> QR(Mat, Num_ , 3);
            QR.solve(Vec_, Sol);
        
            /// P(x,y) = a + bx + cy
            /// @ cut cell the (x,y) = (0,0) as all x_, y_ are computed wrt cut cell
            /// so P = a @ cut cell
            Phi(cellid_) = Sol[0];
            Phi(cellid_) = std::min(1.0, Phi(cellid_));
            Phi(cellid_) = std::max(-1.0, Phi(cellid_));

            for(i_ = 0 ; i_ < Num_ ;i_++) { delete [] Mat[i_] ; }
                delete [] Mat ; 
                delete [] Sol ; delete [] Vec_ ;
        }

        void QR_LS_PhaseField_Second_Order_Neumann
        (
            bool& sing,
            int Num_,
            int Num_i_,
            double* x_, double* y_,
            double* Phi_,
            double* w_,
            amrex::Array4<amrex::Real> const& Phi,
            InterceptData& idata
        )
        {
            const amrex::LinOpBCType bc = amrex::LinOpBCType::Dirichlet;
            //int i_, j_ ;
            //double **Mat, *Sol, *Vec_  ;

            //Allocate_2D_R(Mat,Num_ - Num_i_,3) ; Sol = new double[3] ; 
            //Vec_ = new double[Num_ - Num_i_] ;

            int i_, j_, k_, N_, Size_Loc = 3;
            double **Mat;    // Least squares and constraint matrices
            double **Q, **R; // Q R decomosition of Mat;
            double **Q1, **Q2T, **Q2T_Q, **Q2T_R;
            double *Sol, *Vec_, *Constr_, *Sol_u_, *Sol_w_, *Q1Tb; // *Q2T_QTQ1Tb ;
            double dx, dy, du, dv, nx, ny, dP;
            auto &PsiX = idata.psix_;
            auto &PsiY = idata.psiy_;

            N_ = Num_ - Num_i_;

            Allocate_2D_R(Mat, Num_, Size_Loc);
            Sol = new double[Size_Loc];
            Vec_ = new double[ N_];
            Allocate_2D_R(Q,  Num_, Size_Loc);
            Allocate_2D_R(R, Size_Loc, Size_Loc);
            Allocate_2D_R(Q1,  N_, Size_Loc);
            Allocate_2D_R(Q2T, Size_Loc,  Num_i_);
            Allocate_2D_R(Q2T_Q, Size_Loc,  Num_i_);
            Allocate_2D_R(Q2T_R, Num_i_, Num_i_);
            Constr_ = new double[Num_i_];
            Sol_u_ = new double[Num_i_];
            Sol_w_ = new double[Num_i_];
            Q1Tb = new double[Size_Loc]; // Q2T_QTQ1Tb = new double[Num_i_] ;

            for(i_ = Num_i_ ; i_ < Num_ ; i_++)
            {
                int j_ = i_ - Num_i_;
                Mat[j_][0] = w_[i_] ; // 1 
                Mat[j_][1] = w_[i_]*x_[i_] ; // x
                Mat[j_][2] = w_[i_]*y_[i_] ; // y
                Vec_[j_] = w_[i_]*Phi_[i_] ; // P
                /*
                if(idata.cellid_[0] == 16402 && idata.cellid_[1] == 7)
                        {
                    amrex::Print()<<"Mat["<<j_<<"_][:] = "<<Mat[j_][0]/w_[i_]<<" , "<<Mat[j_][1]/w_[i_]<<" , "<<Mat[j_][2]/w_[i_]<<'\n';
                    amrex::Print()<<"Vec_["<<j_<<"] = "<<Vec_[j_]/w_[i_]<<"\n";
                }
                */
            }
            if (bc == amrex::LinOpBCType::Neumann)
            {
                for (i_ = 0; i_ < Num_i_; i_++)
                {
                    dx = x_[i_];
                    dy = y_[i_];
                    j_ = N_ + i_;

                    nx = PsiX[i_] / sqrt(PsiX[i_] * PsiX[i_] + PsiY[i_] * PsiY[i_]);
                    ny = sqrt(1.0 - ny * ny);
                    nx = PsiX[i_];
                    ny = PsiY[i_];
                    Mat[j_][0] = 0.0;                               // 1
                    Mat[j_][1] = nx;                                // x
                    Mat[j_][2] = ny;                                // y
                    Constr_[i_] = 0.0;
                    /*
                    if(idata.cellid_[0] == 16402 && idata.cellid_[1] == 7)
                            {
                                amrex::Print()<<"Mat["<<j_<<"_][:] = "<<Mat[j_][0]<<" , "<<Mat[j_][1]<<" , "<<Mat[j_][2]<<'\n';
                        amrex::Print()<<"Constr_["<<i_<<"] = "<<Constr_[i_]<<"\n";
                    }
                    */
                }
            }
            else if (bc == amrex::LinOpBCType::Dirichlet)
            {
                for (i_ = 0; i_ < Num_i_; i_++)
                {
                    dx = x_[i_];
                    dy = y_[i_];
                    dP = Phi_[i_];
                    j_ = Num_ - Num_i_ + i_;
                    Mat[j_][0] = 1.0;             // 1
                    Mat[j_][1] = dx;              // x
                    Mat[j_][2] = dy;              // y
                    Constr_[i_] = dP;
                }
            }
            //QRdcmp<double> QR(Mat, Num_ - Num_i_, 3);
            //QR.solve(Vec_, Sol);

            QRdcmp<double> QR(Mat, Num_, Size_Loc);
            QR.get_Q(Q);
            QR.get_R(R);

	    //amrex::Print()<<"point : "<<idata.cellid_[0]<<" , "<<idata.cellid_[1]<<'\n';
	    //amrex::Print()<<"Num_ = "<<Num_<<"Num_i_ = "<<Num_i_<<'\n';

            for (j_ = 0; j_ < Size_Loc; j_++)
            {
                for (i_ = 0; i_ < N_; i_++)
                    Q1[i_][j_] = Q[i_][j_];
                for (i_ = 0; i_ < Num_i_; i_++)
                {
                    Q2T[j_][i_] = Q[N_ + i_][j_];
		    //amrex::Print()<<"Q2T["<<j_<<"]["<<i_<<"] = "<<Q2T[j_][i_]<<'\n';
                }
            }
            QRdcmp<double> QR1(Q2T, Size_Loc, Num_i_);
            QR1.get_Q(Q2T_Q);
            QR1.get_R(Q2T_R);

            // solve for (Q2T_R)T u = Constr ;
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                Sol_u_[i_] = Constr_[i_] / Q2T_R[i_][i_];
                for (j_ = 0; j_ < i_; j_++)
                    Sol_u_[i_] -= Sol_u_[j_] * Q2T_R[j_][i_] / Q2T_R[i_][i_];
            }

            // Find the vector Q1T b.
            for (i_ = 0; i_ < Size_Loc; i_++)
            {
                Q1Tb[i_] = 0.0;
                for (j_ = 0; j_ <N_; j_++)
                    Q1Tb[i_] += Q1[j_][i_] * Vec_[j_];
            }

            // Find the vector 2.0*(Q2T_Q)T Q1Tb -2.0*Sol_u_;
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                Vec_[i_] = -2.0 * Sol_u_[i_];
                for (j_ = 0; j_ < Size_Loc; j_++)
                    Vec_[i_] += 2.0 * Q2T_Q[j_][i_] * Q1Tb[j_];
            }
            // solve (Q2T_R) w = RHS from above
            for (i_ = Num_i_ - 1; i_ >= 0; i_--)
            {
                Sol_w_[i_] = Vec_[i_] / Q2T_R[i_][i_];
                for (j_ = i_ + 1; j_ < Num_i_; j_++)
                    Sol_w_[i_] -= Sol_w_[j_] * Q2T_R[i_][j_] / Q2T_R[i_][i_];
            }

            // Find the vector (Q1T_b) - Q2Tw/2;
            // Find the vector (Q1T_b) - Q2Tw/2;
            for (i_ = 0; i_ < Size_Loc; i_++)
            {
                Vec_[i_] = Q1Tb[i_];
                for (j_ = 0; j_ < Num_i_; j_++)
                    Vec_[i_] -= 0.5 * Q2T[i_][j_] * Sol_w_[j_];
            }
            // solve (Mat_R) w = RHS from above
            for (i_ = 2; i_ >= 0; i_--)
            {
                Sol[i_] = Vec_[i_] / R[i_][i_];
                for (j_ = i_ + 1; j_ < Size_Loc; j_++)
                    Sol[i_] -= Sol[j_] * R[i_][j_] / R[i_][i_];
            }


            /// P(x,y) = a + bx + cy
            /// @ cut cell the (x,y) = (0,0) as all x_, y_ are computed wrt cut cell
            /// so P = a @ cut cell
            Phi(idata.cellid_) = Sol[0];
            /*
	    if(idata.cellid_[0] == 16402 && idata.cellid_[1] == 7)
	    {
		for (i_ = 2; i_ >= 0; i_--)
		{
	            amrex::Print()<<"Vec_["<<i_<<"] = "<<Vec_[i_]<<"\n";
	            for (j_ = 2; j_ >= 0; j_--)
	            {
		        amrex::Print()<<"R["<<i_<<"]["<<j_<<"] = "<<R[i_][j_]<<"\n";
		    }
		}
		amrex::Print()<<"cell = "<<idata.cellid_[0]<<" , "<<idata.cellid_[1]<<'\n';
	        amrex::Print()<<"Phi = "<<Phi(idata.cellid_)<<'\n';
		amrex::Print()<<Sol[0]<<" , "<<Sol[1]<<" , "<<Sol[2]<<"\n";
	    }
	    */

            Phi(idata.cellid_) = std::min(1.0, Phi(idata.cellid_));
            Phi(idata.cellid_) = std::max(-1.0, Phi(idata.cellid_));

            if (bc == amrex::LinOpBCType::Neumann)
            {
                for (i_ = 0; i_ < Num_i_; i_++)
                {
                    dx = x_[i_];
                    dy = y_[i_];
                    idata.phi_Int[i_] = Sol[0] + Sol[1] * dx + Sol[2] * dy;
                    idata.phi_Int[i_] = std::min(1.0,idata.phi_Int[i_]);
                    idata.phi_Int[i_] = std::max(-1.0,idata.phi_Int[i_]);
                }
            }
            else if (bc == amrex::LinOpBCType::Dirichlet)
            {
                for (i_ = 0; i_ < Num_i_; i_++)
                    idata.phi_Int[i_] = -1.0;
            }



	    //if(Phi(idata.cellid_) < 0) Phi(idata.cellid_) = -1.0;

            //for(i_ = 0 ; i_ < Num_ - Num_i_ ;i_++) { delete [] Mat[i_] ; }
            //    delete [] Mat ;
            //    delete [] Sol ; delete [] Vec_ ;
            for (i_ = 0; i_ < Num_; i_++)
            {
                delete[] Mat[i_];
                delete[] Q[i_];
            }
            for (i_ = 0; i_ < Size_Loc; i_++)
            {
                delete[] R[i_];
                delete[] Q2T[i_];
                delete[] Q2T_Q[i_];
            }
            for (i_ = 0; i_ < N_; i_++)
            {
                delete[] Q1[i_];
            }
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                delete[] Q2T_R[i_];
            }
            delete[] Mat;
            delete[] Q;
            delete[] R;
            delete[] Q1;
            delete[] Q2T;
            delete[] Q2T_Q;
            delete[] Q2T_R;
            delete[] Sol;
            delete[] Vec_;
            delete[] Constr_;
            delete[] Sol_u_;
            delete[] Sol_w_;
            delete[] Q1Tb;
            /*
            int i_, j_, k_, N_;
            double **Mat;    // Least squares and constraint matrices
            double **Q, **R; // Q R decomosition of Mat;
            double **Q1, **Q2T, **Q2T_Q, **Q2T_R;
            double *Sol, *Vec_, *Constr_, *Sol_u_, *Sol_w_, *Q1Tb; // *Q2T_QTQ1Tb ;
            double dx, dy, du, dv, nx, ny;
            Allocate_2D_R(Q,  Num_, 3);
            Allocate_2D_R(R, 3, 3);
            Allocate_2D_R(Q1,  N_, 3);
            Allocate_2D_R(Q2T, 3,  Num_i_);
            Allocate_2D_R(Q2T_Q, 3,  Num_i_);
            Allocate_2D_R(Q2T_R, Num_i_, Num_i_);
	    */
        }

    } // namespace AdvectLSIntPhaseField
    
    
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
                        Phi_[ni] = -1.0;//idt.T_Int[ni];
                        Compute_LSQ_Weights(w_[ni], std::fabs(idt.frac_[ni]));
                    }
                    int Num_i_ = ni;
		    //if(Num_i_ > 2)
		    //{
		    //    amrex::Print()<<"weird point:"<<idt.cellid_[0]<<" , "<<idt.cellid_[1]<<'\n';
		    //}
    
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
                    if (PORDER_ == 2)
                    {
                        if (Num_ < 11 + Num_i_ || Num_i_ > 2)
                        {
                            /// too few points have been found...can only do a simple average
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
                        else
                        {
                            //amrex::Print()<<"Second order for Phi"<<"\n";
                            AdvectLSIntPhaseField::QR_LS_PhaseField_Second_Order_Neumann(Sing_, Num_, Num_i_, x_, y_, Phi_, w_, Phi, idt);
                            if(Sing_)
                            {
                                sum_w = sum_sol = 0.0;
                                for(int i = Num_i_; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol += w_[i] * Phi_[i];
                                }
                                Phi(idt.cellid_) = sum_sol / sum_w;
                                for(ni = 0; ni < idt.n_intercepts; ni++)
                                {
                                    idt.phi_Int[ni] = sum_sol / sum_w;
                                }
                            }
                            //amrex::PrintToFile("log") << " Second first ";
                        }
                    }
		     
    		    //if(idt.cellid_[0] == 16454 && idt.cellid_[1] == 0)
    		    //{
    			    //amrex::Print()<<"cell : "<<idt.cellid_<<'\n';
                            //amrex::Print()<<"Phi : "<< Phi(idt.cellid_)<<'\n';
                            //amrex::Print()<<"Num_ = "<<Num_<<" , Num_i_ = "<<Num_i_<<'\n';
                            //for(int ii = 0;ii< Num_;ii++)
                            //{
                                    //amrex::Print()<<"ii = "<<ii<<" , w_[ii] = "<<w_[ii]<<"\t";
                            	    //amrex::Print()<<", x_ = "<<x_[ii]<<"\t";
                            	    //amrex::Print()<<", y_ = "<<y_[ii]<<"\t";
                            	    //amrex::Print()<<", Phi_ = "<<Phi_[ii]<<"\n";
                            //}
                            //std::exit(9);
    		    //}

                }
            }
            nsolid++;
        }
    
        /// as pressure values are updated
        /// fill the internal ghost values
        mf.FillBoundary();
    }         


    void Mask::ExtendPhaseFieldLSQ(amrex::MultiFab &mfPhi)
    {
        if (interfaces->empty())
        {
            return;
        }
        const amrex::Real *prob_lo = geom_.ProbLo();
        const amrex::Real *dx = geom_.CellSize();
        const amrex::Box &domain = geom_.Domain();
    
        bool Sing_;
        double sum_w, sum_sol_u, sum_sol_v;
    
        int N_Max = (2 * stencil_ + 1) * (2 * stencil_ + 1) + 4;
        double x_[N_Max], y_[N_Max], w_[N_Max], Phi_[N_Max];
    
        for (int outer = 0; outer < LAYERS; outer++)
        {
            mfPhi.FillBoundary();
            for (amrex::MFIter mfi(PMask); mfi.isValid(); ++mfi)
            {
                amrex::Array4<int> const &pmask = PMask.array(mfi);
                amrex::Array4<amrex::Real> const &Phi = mfPhi.array(mfi);

                auto &Index_Additional_Layers = Index_Additional_Layers_LD[mfi];
    
                for (auto &&cell : Index_Additional_Layers[outer])
                {
                    Sing_ = false;
                    amrex::Box bx(cell, cell);
                    bx.grow(stencil_);
                    amrex::Box bx_isect = bx & domain;
		    int iscalar = 0;

                    int i_ = 0;
                    for (amrex::BoxIterator bit(bx_isect); bit.ok(); ++bit)
                    {
                        const amrex::IntVect &iv = bit();
                        if (outer == 0)
                        {
                            if (pmask(iv) >= 1 && pmask(iv) <= 20)
                            {
                                if (iv == cell)
                                {
                                }
                                else
                                {
                                    LSQ_Parameters(Phi_[i_], x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], Phi, prob_lo, dx);
                                    i_++;
                                }
                            }
                        }
                        else
                        {
                            if (pmask(iv) >= 1 && pmask(iv) <= 99+outer)
                            {
                                if (iv == cell)
                                {
                                }
                                else
                                {
                                    LSQ_Parameters(Phi_[i_], x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], Phi, prob_lo, dx);
                                    i_++;
                                }
                            }
                        }
                    }
                    int Num_ = i_;

	            {
			//amrex::Print()<<"Second order extension"<<'\n';
		        {
			    AdvectLSIntPhaseField::QR_LS_PhaseField_Second_Order(Sing_, Num_, x_, y_, Phi_, w_, Phi, cell);
			}
                    }
		    /*
                    if(cell[0] == 16402 && cell[1] == 7)
                    {
                        amrex::Print(-1)<<"i = "<<cell[0]<<" , j = "<<cell[1]<<'\n';
                        amrex::Print(-1)<<"pmask = "<<pmask(cell)<<'\n';
                        amrex::Print(-1)<<"outer = "<<outer<<'\n';
                        amrex::Print(-1)<<"Num_ = "<<Num_<<'\n';
                        amrex::Print(-1)<<"Phi = "<<Phi(cell)<<'\n';
                        
                        for(int ii = 0; ii<Num_;ii++)
                        {
                            amrex::Print(-1)<<"ii = "<<ii<<'\n';
                            amrex::Print(-1)<<"x_[ii] = "<<x_[ii]<<" , y_[ii] = "<<y_[ii]<<'\n';
                            amrex::Print(-1)<<"Phi_[ii] = "<<Phi_[ii]<<'\n';
                        }
                    }*/

                }
            }
            //CollocatedVelBoundaryConditions();
        }
    }    
} /*End namespace mycode */

