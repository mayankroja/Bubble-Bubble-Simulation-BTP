#include "Mask.H"
#include <AMReX_ParmParse.H>
#include <qr.h>
#include <WeightedENO.h>
#include <Viscosity.H>

namespace mycode
{
    namespace
    {
        // Sign function
        int I_DSign(double x)
        {
            if (x >= 0.0)
                return 1;
            else
                return -1;
        }
        
        double Sign(double x)
        {
            if (x >= 0.0)
                return 1.0;
            else
                return -1.0;
        }

        const double exp_ = 2.0;
                
        int Cubic_Solve(double& x_, double a_, double b_, double c_, double d_) 
        {
        	int iter_cu = 0 ;
        	double a, b , c , d, Res, den, x_o ;
        
        	x_o = x_ ;
        	if( fabs(x_o) > 1.0 ) {
        		//amrex::PrintToFile("log") << "\n" << "Second-order estimate out of bounds! " << "\n" ;
        		exit(1) ;
        	}
        	// Note the negative signs compared to Nourgialev, .. h * theta = x_i - xI
        	a = -( 3.0*(b_ - c_) - a_ + d_ )/6.0 ;
        	b = (c_ + a_ - 2.0*b_)/2.0 ;
        	c = (3.0*b_ + 2.0*a_ - 6.0*c_ + d_)/6.0 ;
        	d = b_ ;
        
        	double x_L, x_U, F_L, F_U, x_m, F_m ;
                x_L = 0.0 ; x_U = Sign(x_) ; x_m = x_ ; bool STOP = false ;
                F_L = a*x_L*x_L*x_L + b*x_L*x_L + c*x_L + d ;
                F_U = a*x_U*x_U*x_U + b*x_U*x_U + c*x_U + d ;
                F_m = a*x_m*x_m*x_m + b*x_m*x_m + c*x_m + d ;
        
                if(F_L*F_m < 0.0) { x_U = x_m ; F_U = F_m ; }
                else if(F_U*F_m < 0.0) {x_L = x_m ; F_L = F_m ; }
                else {
                        //amrex::PrintToFile("log") << "\n Found the root without iterations!" << "\n" ;
                        if( fabs(F_L) < 1.0E-15 ) { x_ = x_L ; STOP = true ; }
                        else if( fabs(F_U) < 1.0E-15 ) { x_ = x_U ; STOP = true ; }
                        else if( fabs(F_m) < 1.0E-15 ) { x_ = x_m ; STOP = true ; }
                }
        
        	do {
        
                        if( fabs(F_L) < 1.0E-15 ) { x_ = x_L ; STOP = true ; }
                        else if( fabs(F_U) < 1.0E-15 ) { x_ = x_U ; STOP = true ; }
                        else {
                                x_m = 0.5*(x_L + x_U) ;
                                F_m = a*x_m*x_m*x_m + b*x_m*x_m + c*x_m + d ;
                                if(fabs(F_m) < 1.0E-15) { x_ = x_m ; STOP = true; }
                                if( I_DSign(F_m) == I_DSign(F_L) ) { x_L = x_m ; F_L = F_m ; }
                                if( I_DSign(F_m) == I_DSign(F_U) ) { x_U = x_m ; F_U = F_m ; }
                        }
                        iter_cu++ ;
                //        amrex::PrintToFile("log") << iter_cu << "\t" << F_m << "\t" << x_m <<  "\n" ;
                } while(iter_cu < 10 && !STOP) ;
                x_ = 0.5*(x_L + x_U) ;
        
        	do {
        		den = 3.0*a*x_*x_ + 2.0*b*x_ + c ;
        		if( fabs(den) < 1.0E-12) { return -1 ; }
        		Res = (a*x_*x_*x_ + b*x_*x_ + c*x_ + d)/den ;
        		x_ -= Res ;
        		iter_cu++ ;	
        	} while(iter_cu < 21 && fabs(Res) > 1.0E-12) ;
        
        	if(fabs(x_) > 1.0) {
        		//amrex::PrintToFile("log") << "\n" << "Cubic solve out of bounds!" << "\n" ;
        		//amrex::PrintToFile("log") << "\n" << "Number of Newton's iterations: " << iter_cu << "\n" ;
        		//amrex::PrintToFile("log") << x_o << "\t" << x_ << "\n" ;
        		//amrex::PrintToFile("log") << a_ << "\t" << b_ << "\t" << c_ << "\t" << d_ << "\n" ; 
        		//amrex::PrintToFile("log") << "Cubic: " << a << "\t" << b << "\t" << c << "\t" << d << "\n" ; 
        		exit(1) ; 
        	}
        	
        	return iter_cu ;
        }
        
        void Compute_Normal_Curvature
        (
            amrex::Real* dPsi, 
            amrex::Real& Curv, 
            int i, int j,
            amrex::Array4<amrex::Real const> const& Psi,
            const amrex::Real& deltax,
            const amrex::Real* dx,
            const amrex::Real* prob_lo,
            bool axisymmetric
        ) 
        {
        	amrex::Real Psi_xx, Psi_xy, Psi_yy ;
            int k = 0;
            dPsi[0] = 0.5 * (Psi(i + 1, j, k) - Psi(i - 1, j, k));
            dPsi[1] = 0.5 * (Psi(i, j + 1, k) - Psi(i, j - 1, k));
            Psi_xx = Psi(i + 1, j, k) - 2.0 * Psi(i, j, k) + Psi(i - 1, j, k);
            Psi_yy = Psi(i, j + 1, k) - 2.0 * Psi(i, j, k) + Psi(i, j - 1, k);
            Psi_xy = 0.25 * (Psi(i + 1, j + 1, k) - Psi(i + 1, j - 1, k) - Psi(i - 1, j + 1, k) + Psi(i - 1, j - 1, k));
        
            Curv = (Psi_xx * (dPsi[1] * dPsi[1]) + Psi_yy * (dPsi[0] * dPsi[0]) - 2.0 * dPsi[0] * dPsi[1] * Psi_xy) / (deltax * pow(dPsi[0] * dPsi[0] + dPsi[1] * dPsi[1] + 1.0E-14, 1.5));
        
            if (axisymmetric)
            {
                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                Curv += (1.0 / y) * dPsi[1] / sqrt(dPsi[0] * dPsi[0] + dPsi[1] * dPsi[1] + 1.0E-14);
            }
        }
        
        void Compute_LSQ_Weights(double& weights_, double dist_) 
        {
            weights_ = 1.0 / pow((dist_ * dist_ + 1.0), exp_);
        }
        
        
        void Fij_LSQ_Parameters(
            double &F11_,double &F12_,
            double &F21_,double &F22_,
	    double &F33_,
            double &x_,
            double &y_,
            double &w_,
            int i1, int j1,
            int i, int j,
            amrex::Array4<amrex::Real const> const &F11, /// cc vel with 2 comp
	    amrex::Array4<amrex::Real const> const &F12,
	    amrex::Array4<amrex::Real const> const &F21,
	    amrex::Array4<amrex::Real const> const &F22,
	    amrex::Array4<amrex::Real const> const &F33,
            const amrex::Real *prob_lo,
            const amrex::Real *dx)
        {
            amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
            amrex::Real x1 = prob_lo[0] + dx[0] * (i1 + 0.5);
            amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
            amrex::Real y1 = prob_lo[1] + dx[1] * (j1 + 0.5);
        
            x_ = (x1 - x) / dx[0];
            y_ = (y1 - y) / dx[1];
            F11_ = F11(i1, j1, 0);
            F12_ = F12(i1, j1, 0);
	    F21_ = F21(i1, j1, 0);
            F22_ = F22(i1, j1, 0);
	    F33_ = F33(i1, j1, 0);
        
            Compute_LSQ_Weights(w_, std::sqrt(x_ * x_ + y_ * y_));
        }
        
        void RefConfg_LSQ_Parameters(
            double &X_,
            double &Y_,
            double &x_,
            double &y_,
            double &w_,
            int i1, int j1,
            int i, int j,
            amrex::Array4<amrex::Real const> const &X_ref, /// cc vel with 2 comp
            amrex::Array4<amrex::Real const> const &Y_ref, /// cc vel with 2 comp
            const amrex::Real *prob_lo,
            const amrex::Real *dx)
        {
            amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
            amrex::Real x1 = prob_lo[0] + dx[0] * (i1 + 0.5);
            amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
            amrex::Real y1 = prob_lo[1] + dx[1] * (j1 + 0.5);
        
            x_ = (x1 - x) / dx[0];
            y_ = (y1 - y) / dx[1];
            X_ = X_ref(i1, j1, 0, 0);
            Y_ = Y_ref(i1, j1, 0, 1);
        
            Compute_LSQ_Weights(w_, std::sqrt(x_ * x_ + y_ * y_));
        }
    }
    namespace FillFij
    {
        void QR_Second_Order
        (
            bool& sing,
            int Num_,
            int Num_i_,
            double* x_, double* y_,
            double* T_,
            double* w_,
            amrex::Array4<amrex::Real> const& T,
            InterceptData& idata,
	    int iscalar
        )
        {
            sing = false;
            int i_, j_ ;
            double **Mat, *Sol, *Vec_  ;
	    double **R;

            Allocate_2D_R(Mat,Num_ - Num_i_,3) ; 
            Allocate_2D_R(R, 3, 3);

	    Sol = new double[3] ; 
	    Vec_ = new double[Num_ - Num_i_] ;

            for(i_ = Num_i_ ; i_ < Num_ ; i_++)
            {
	        j_ = i_ - Num_i_;	
                Mat[j_][0] = w_[i_] ; // 1 
                Mat[j_][1] = w_[i_]*x_[i_] ; // x
                Mat[j_][2] = w_[i_]*y_[i_] ; // y
                Vec_[j_] = w_[i_]*T_[i_] ; // P
            }
            QRdcmp<double> QR(Mat, Num_ - Num_i_, 3);
            QR.solve(Vec_, Sol);

	    QR.get_R(R);

            /// P(x,y) = a + bx + cy
            /// @ cut cell the (x,y) = (0,0) as all x_, y_ are computed wrt cut cell
            /// so P = a @ cut cell
            T(idata.cellid_) = Sol[0];

            for (i_ = 0; i_ < Num_i_; i_++)
            {
	        if(iscalar == 2)
	            idata.F11_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_];
                else if(iscalar == 3)
                    idata.F12_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_];
	        else if(iscalar == 4)
                    idata.F21_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_];
                else if(iscalar == 5)
                    idata.F22_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_];
	        else if(iscalar == 6)
                    idata.F33_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_];
	    }

            if(std::isnan(Sol[0])) sing = true;
            for(i_ = 0 ; i_ < 3; i_++)
            {
                if(std::abs(R[i_][i_]) < 1.0e-6) sing = true;
            }

            for(i_ = 0 ; i_ < Num_ - Num_i_ ;i_++) { delete [] Mat[i_] ; }
            delete [] Mat ;
            delete [] Sol ; delete [] Vec_ ;

            for (i_ = 0; i_ < 3; i_++)
            {
                delete[] R[i_];
	    }
	    delete[] R;
        }
        
        
        void QR_Second_Order_F11F21_2D
        (
            bool& sing,
            int Num_,
            int Num_i_,
            double* x_, double* y_,
            double* F11_,
	    double* F21_,
            double* w_,
            amrex::Array4<amrex::Real> const& F11,
	    amrex::Array4<amrex::Real> const& F21,
            InterceptData& idata
        )
        {
            int i_, j_, N_ ;
            double **Mat, *Sol, *Vec_  ;
	    double dx, dy, dF11, dF21, nx, ny;

	    N_ = Num_ - Num_i_;
            Allocate_2D_R(Mat, 2*N_,5) ; 
	    Sol = new double[5] ; 
	    Vec_ = new double[ 2*N_] ;

	    //amrex::Print()<<"cell "<<idata.cellid_[0]<<" , "<<idata.cellid_[1]<<" , "<<Num_<<" , "<<Num_i_<<'\n';

            for(i_ = 0 ; i_ < N_ ; i_++)
            {
		j_ = i_ + Num_i_;

                Mat[i_][0] = w_[j_] ; // 1 
                Mat[i_][1] = w_[j_]*x_[j_] ; // x
                Mat[i_][2] = w_[j_]*y_[j_] ; // y
                Mat[i_][3] = 0.0 ; // y
		Mat[i_][4] = 0.0 ; // y

                Mat[i_ + N_][0] = 0.0 ; // 1 
                Mat[i_ + N_][1] = -1.0*w_[j_] * y_[j_] ; // x
                Mat[i_ + N_][2] = 0.0 ; // y
                Mat[i_ + N_][3] = w_[j_] ; // y
                Mat[i_ + N_][4] = w_[j_] * x_[j_] ; // y


                Vec_[i_] = w_[j_]*F11_[j_] ; // F11
		Vec_[i_ + N_] = w_[j_]*F21_[j_] ; // F11
            }
            QRdcmp<double> QR(Mat, 2*N_ , 5);
            QR.solve(Vec_, Sol);

            /// P(x,y) = a + bx + cy
            /// @ cut cell the (x,y) = (0,0) as all x_, y_ are computed wrt cut cell
            /// so P = a @ cut cell
            F11(idata.cellid_) = Sol[0];
	    F21(idata.cellid_) = Sol[3];

            for (i_ = 0; i_ < Num_i_; i_++)
            {
	            idata.F11_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_];
                    idata.F21_Int[i_] = Sol[3] - Sol[1] * y_[i_] + Sol[4] * x_[i_];
	    }

            for(i_ = 0 ; i_ <  2*N_  ;i_++) { delete [] Mat[i_] ; }
                delete [] Mat ;
                delete [] Sol ; delete [] Vec_ ;

            sing = false;
            if(std::isnan(Sol[0])) sing = true;
        }        
        

        void QR_Third_Order
        (
            bool& sing,
            int Num_,
            int Num_i_,
            double* x_, double* y_,
            double* T_,
            double* w_,
            amrex::Array4<amrex::Real> const& T,
            InterceptData& idata,
	    int iscalar
        )
        {
	    sing = false;
            int i_, j_ ;
            double **Mat, *Sol, *Vec_  ;
	    double **Q, **R;

	    //amrex::Print()<<"Num_ = "<<Num_<<", num_i_ = "<<Num_i_<<'\n';

            Allocate_2D_R(Mat,Num_ - Num_i_,6) ; 
	    //Allocate_2D_R(Q, Num_ - Num_i_, 6) ;
	    Allocate_2D_R(R, 6, 6) ;
	    Sol = new double[6] ; 

            //Allocate_2D_R(Mat,Num_ - Num_i_,3) ;
	    //Allocate_2D_R(Q, Num_ - Num_i_, 3) ;
            //Allocate_2D_R(R, 3, 3) ;
            //Sol = new double[3] ;

	    Vec_ = new double[Num_ - Num_i_] ;

            for(i_ = Num_i_ ; i_ < Num_ ; i_++)
            {
		//amrex::Print()<<"i_ = "<<i_<<'\n';
	        j_ = i_ - Num_i_;	
                Mat[j_][0] = w_[i_] ; // 1
                Mat[j_][1] = w_[i_]*x_[i_] ; // x
                Mat[j_][2] = w_[i_]*y_[i_] ; // y
                Mat[j_][3] = w_[i_]*x_[i_]*x_[i_] ; // x^2
                Mat[j_][4] = w_[i_]*x_[i_]*y_[i_] ; // xy
                Mat[j_][5] = w_[i_]*y_[i_]*y_[i_] ; // y^2
                Vec_[j_] = w_[i_]*T_[i_] ; // x
            }
            QRdcmp<double> QR(Mat, Num_ - Num_i_, 6);

            QR.solve(Vec_, Sol);
            QR.get_R(R);

            /// P(x,y) = a + bx + cy
            /// @ cut cell the (x,y) = (0,0) as all x_, y_ are computed wrt cut cell
            /// so P = a @ cut cell
            T(idata.cellid_) = Sol[0];
            if(std::isnan(Sol[0])) sing = true;
            for(i_ = 0 ; i_ < 6; i_++)
            {
                if(std::abs(R[i_][i_]) < 1.0e-6) sing = true;
            }

            for (i_ = 0; i_ < Num_i_; i_++)
            {
		
                if(iscalar == 2)
	                idata.F11_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_]
			                + Sol[3] * x_[i_] * x_[i_] + Sol[4] * x_[i_] * y_[i_]
					        + Sol[5] * y_[i_] * y_[i_];
                else if(iscalar == 3)
                    idata.F12_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_]
			                + Sol[3] * x_[i_] * x_[i_] + Sol[4] * x_[i_] * y_[i_]
                                        + Sol[5] * y_[i_] * y_[i_];
                else if(iscalar == 4)
                    idata.F21_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_]
                                        + Sol[3] * x_[i_] * x_[i_] + Sol[4] * x_[i_] * y_[i_]
                                        + Sol[5] * y_[i_] * y_[i_];
                else if(iscalar == 5)
                    idata.F22_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_]
                                        + Sol[3] * x_[i_] * x_[i_] + Sol[4] * x_[i_] * y_[i_]
                                        + Sol[5] * y_[i_] * y_[i_];
                else if(iscalar == 6)
                    idata.F33_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_]
                                        + Sol[3] * x_[i_] * x_[i_] + Sol[4] * x_[i_] * y_[i_]
                                        + Sol[5] * y_[i_] * y_[i_];
            
            }

            for(i_ = 0 ; i_ < Num_ - Num_i_ ;i_++) { delete [] Mat[i_] ; }
            delete [] Mat ;
            delete [] Sol ; delete [] Vec_ ;

            for (i_ = 0; i_ < 6; i_++)
            {
                delete[] R[i_];
            }
            delete[] R;
        }   

        void QR_Third_Order_DivFr
        (
            bool& sing,
            int Num_,
            int Num_i_,
            double* x_, double* y_,
            double* F1_,
            double* F2_,
            double* w_,
            amrex::Array4<amrex::Real> const& F1,
            amrex::Array4<amrex::Real> const& F2,
            InterceptData& idata,
            int iscalar1, int iscalar2
        )
        {
            int i_, j_, N_;
            double **Mat, *Sol, *Vec_  ;
	    //double **Q, **R;

	    //amrex::Print()<<"Num_ = "<<Num_<<", num_i_ = "<<Num_i_<<'\n';

            Allocate_2D_R(Mat,2 * (Num_ - Num_i_),9) ; 
	    //Allocate_2D_R(Q, Num_ - Num_i_, 6) ;
	    //Allocate_2D_R(R, 6, 6) ;
	    Sol = new double[9] ; 

            //Allocate_2D_R(Mat,Num_ - Num_i_,3) ;
	    //Allocate_2D_R(Q, Num_ - Num_i_, 3) ;
            //Allocate_2D_R(R, 3, 3) ;
            //Sol = new double[3] ;

	    Vec_ = new double[2 * (Num_ - Num_i_)] ;
	    N_ = Num_ - Num_i_;

            for(i_ = Num_i_ ; i_ < Num_ ; i_++)
            {
		//amrex::Print()<<"i_ = "<<i_<<'\n';
	        j_ = i_ - Num_i_;	
                amrex::Real dx = x_[i_];
	        amrex::Real dy = y_[i_];	

                //Mat[j_][0] = w_[i_] ; // 1
                //Mat[j_][1] = w_[i_]*x_[i_] ; // x
                //Mat[j_][2] = w_[i_]*y_[i_] ; // y
                //Mat[j_][3] = w_[i_]*x_[i_]*x_[i_] ; // x^2
                //Mat[j_][4] = w_[i_]*x_[i_]*y_[i_] ; // xy
                //Mat[j_][5] = w_[i_]*y_[i_]*y_[i_] ; // y^2
		/*
                Mat[j_][0] = w_[i_]*1.0; // w
                Mat[j_][1] = w_[i_]*dx;  // x
                Mat[j_][2] = w_[i_]*dy;  // y
                Mat[j_][3] = 0.0;
                Mat[j_][4] = 0.0;

                Mat[j_ + N_][0] = 0.0; // w
                Mat[j_ + N_][1] = -w_[i_]*dy; // x
                Mat[j_ + N_][2] = 0.0; // y
                Mat[j_ + N_][4] = 1.0;
                Mat[j_ + N_][5] = w_[i_]*dx;
                */
	        Mat[j_][0] = w_[i_] ;
                Mat[j_][1] = w_[i_]*dx;  // x
                Mat[j_][2] = w_[i_]*dy;  // y
                Mat[j_][3] = w_[i_]*dx * dx;
                Mat[j_][4] = w_[i_]*dx * dy;
                Mat[j_][5] = w_[i_]*dy * dy;
                Mat[j_][6] = 0.0;
                Mat[j_][7] = 0.0;
                Mat[j_][8] = 0.0;

                Mat[j_ + N_][0] = 0.0; // w
                Mat[j_ + N_][1] = -w_[i_]*dy; // x
                Mat[j_ + N_][2] = 0.0; // y
                Mat[j_ + N_][3] = -w_[i_]*2.0 * dx * dy;
                Mat[j_ + N_][4] = -w_[i_]*0.5 * dy * dy;
                Mat[j_ + N_][5] = 0.0;
                Mat[j_ + N_][6] = w_[i_]*1.0;
                Mat[j_ + N_][7] = w_[i_]*dx;
                Mat[j_ + N_][8] = w_[i_]*dx * dx;

                Vec_[j_] = w_[i_]*F1_[i_] ; // x
		Vec_[j_ + N_] = w_[i_]*F2_[i_] ; // x
            }
            QRdcmp<double> QR(Mat, 2 * (Num_ - Num_i_), 9);

	    //QRdcmp<double> QR(Mat, Num_ - Num_i_, 3);
            QR.solve(Vec_, Sol);
	    //QR.get_Q(Q);
            //QR.get_R(R);

            /// P(x,y) = a + bx + cy
            /// @ cut cell the (x,y) = (0,0) as all x_, y_ are computed wrt cut cell
            /// so P = a @ cut cell
            F1(idata.cellid_) = Sol[0];
	    F2(idata.cellid_) = Sol[6];

            for (i_ = 0; i_ < Num_i_; i_++)
            {
		amrex::Real Int_F1 = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_]
                                     + Sol[3] * x_[i_] * x_[i_] + Sol[4] * x_[i_] * y_[i_]
                                     + Sol[5] * y_[i_] * y_[i_];	
                if(iscalar1 == 2)
	            idata.F11_Int[i_] = Int_F1; 
                else if(iscalar1 == 3)
                    idata.F12_Int[i_] = Int_F1; 
                else if(iscalar1 == 4)
                    idata.F21_Int[i_] = Int_F1; 
                else if(iscalar1 == 5)
                    idata.F22_Int[i_] = Int_F1;
                else if(iscalar1 == 6)
                    idata.F33_Int[i_] = Int_F1; 

                amrex::Real Int_F2 = -1.0*Sol[1]*y_[i_] - 2.0*Sol[3]*x_[i_]*y_[i_] - 0.5*Sol[4]*y_[i_]*y_[i_]
                                        + Sol[6] + Sol[7]*x_[i_] + Sol[8]*x_[i_]*x_[i_];
                if(iscalar2 == 2)
                    idata.F11_Int[i_] = Int_F2; 
                else if(iscalar2 == 3)
                    idata.F12_Int[i_] = Int_F2;
                else if(iscalar2 == 4)
                    idata.F21_Int[i_] = Int_F2;
                else if(iscalar2 == 5)
                    idata.F22_Int[i_] = Int_F2;
                else if(iscalar2 == 6)
                    idata.F33_Int[i_] = Int_F2; 
            }

            for(i_ = 0 ; i_ < 2 * (Num_ - Num_i_) ;i_++) { delete [] Mat[i_] ; }
                delete [] Mat ;
                delete [] Sol ; delete [] Vec_ ;

	    sing = false;
            if(std::isnan(Sol[0])) sing = true;
        }   

        void QR_Fourth_Order
        (
            bool& sing,
            int Num_,
            int Num_i_,
            double* x_, double* y_,
            double* T_,
            double* w_,
            amrex::Array4<amrex::Real> const& T,
            InterceptData& idata,
            int iscalar
        )
        {
	    sing = false;
            int i_, j_ ;
            double **Mat, *Sol, *Vec_  ;
            double **Q, **R;

            //amrex::Print()<<"Num_ = "<<Num_<<", num_i_ = "<<Num_i_<<'\n';

            Allocate_2D_R(Mat,Num_ - Num_i_,10) ;
            //Allocate_2D_R(Q, Num_ - Num_i_, 6) ;
            Allocate_2D_R(R, 10, 10) ;
            Sol = new double[10] ;

            //Allocate_2D_R(Mat,Num_ - Num_i_,3) ;
            //Allocate_2D_R(Q, Num_ - Num_i_, 3) ;
            //Allocate_2D_R(R, 3, 3) ;
            //Sol = new double[3] ;

            Vec_ = new double[Num_ - Num_i_] ;

            for(i_ = Num_i_ ; i_ < Num_ ; i_++)
            {
                //amrex::Print()<<"i_ = "<<i_<<'\n';
                j_ = i_ - Num_i_;
                Mat[j_][0] = w_[i_] ; // 1
                Mat[j_][1] = w_[i_]*x_[i_] ; // x
                Mat[j_][2] = w_[i_]*y_[i_] ; // y
                Mat[j_][3] = w_[i_]*x_[i_]*x_[i_] ; // x^2
                Mat[j_][4] = w_[i_]*x_[i_]*y_[i_] ; // xy
                Mat[j_][5] = w_[i_]*y_[i_]*y_[i_] ; // y^2
                Mat[j_][6] = w_[i_]*x_[i_]*x_[i_]*x_[i_] ;
                Mat[j_][7] = w_[i_]*x_[i_]*x_[i_]*y_[i_] ;
                Mat[j_][8] = w_[i_]*x_[i_]*y_[i_]*y_[i_] ;
                Mat[j_][9] = w_[i_]*y_[i_]*y_[i_]*y_[i_] ;
                Vec_[j_] = w_[i_]*T_[i_] ; // x
            }
            QRdcmp<double> QR(Mat, Num_ - Num_i_, 10);

            //QRdcmp<double> QR(Mat, Num_ - Num_i_, 3);
            QR.solve(Vec_, Sol);
	    QR.get_R(R);
	    T(idata.cellid_) = Sol[0];
            for (i_ = 0; i_ < Num_i_; i_++)
            {

                if(iscalar == 2)
                    idata.F11_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_]
                                        + Sol[3] * x_[i_] * x_[i_] + Sol[4] * x_[i_] * y_[i_]
                                        + Sol[5] * y_[i_] * y_[i_]
					+ Sol[6] * x_[i_] * x_[i_] * x_[i_]
					+ Sol[7] * x_[i_] * x_[i_] * y_[i_]
					+ Sol[8] * x_[i_] * y_[i_] * y_[i_]
					+ Sol[9] * y_[i_] * y_[i_] * y_[i_];
                else if(iscalar == 3)
                    idata.F12_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_]
                                        + Sol[3] * x_[i_] * x_[i_] + Sol[4] * x_[i_] * y_[i_]
                                        + Sol[5] * y_[i_] * y_[i_]
                                        + Sol[6] * x_[i_] * x_[i_] * x_[i_]
                                        + Sol[7] * x_[i_] * x_[i_] * y_[i_]
                                        + Sol[8] * x_[i_] * y_[i_] * y_[i_]
                                        + Sol[9] * y_[i_] * y_[i_] * y_[i_];
                else if(iscalar == 4)
                    idata.F21_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_]
                                        + Sol[3] * x_[i_] * x_[i_] + Sol[4] * x_[i_] * y_[i_]
                                        + Sol[5] * y_[i_] * y_[i_]
                                        + Sol[6] * x_[i_] * x_[i_] * x_[i_]
                                        + Sol[7] * x_[i_] * x_[i_] * y_[i_]
                                        + Sol[8] * x_[i_] * y_[i_] * y_[i_]
                                        + Sol[9] * y_[i_] * y_[i_] * y_[i_];
                else if(iscalar == 5)
                    idata.F22_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_]
                                        + Sol[3] * x_[i_] * x_[i_] + Sol[4] * x_[i_] * y_[i_]
                                        + Sol[5] * y_[i_] * y_[i_]
                                        + Sol[6] * x_[i_] * x_[i_] * x_[i_]
                                        + Sol[7] * x_[i_] * x_[i_] * y_[i_]
                                        + Sol[8] * x_[i_] * y_[i_] * y_[i_]
                                        + Sol[9] * y_[i_] * y_[i_] * y_[i_];
                else if(iscalar == 6)
                    idata.F33_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_]
                                        + Sol[3] * x_[i_] * x_[i_] + Sol[4] * x_[i_] * y_[i_]
                                        + Sol[5] * y_[i_] * y_[i_]
                                        + Sol[6] * x_[i_] * x_[i_] * x_[i_]
                                        + Sol[7] * x_[i_] * x_[i_] * y_[i_]
                                        + Sol[8] * x_[i_] * y_[i_] * y_[i_]
                                        + Sol[9] * y_[i_] * y_[i_] * y_[i_];

            }

            if(std::isnan(Sol[0])) sing = true;
            for(i_ = 0 ; i_ < 10; i_++)
            {
                if(std::abs(R[i_][i_]) < 1.0e-5) sing = true;
            }

            for(i_ = 0 ; i_ < Num_ - Num_i_ ;i_++) { delete [] Mat[i_] ; }
                delete [] Mat ;
                delete [] Sol ; delete [] Vec_ ;

            for (i_ = 0; i_ < 10; i_++)
            {
                delete[] R[i_];
            }
            delete[] R;
        }

        void LS_QR_Second_Order_F11F21_2D
        (
            bool& sing,
            int Num_,
            int Num_i_,
            double* x_, double* y_,
            double* F11_,
            double* F21_,
            double* w_,
            amrex::Array4<amrex::Real> const& F11,
            amrex::Array4<amrex::Real> const& F21,
	    const amrex::IntVect& cellid_
        )
        {
            int i_, j_, N_ ;
            double **Mat, *Sol, *Vec_  ;
            double dx, dy, dF11, dF21, nx, ny;

            N_ = Num_ - Num_i_;
            Allocate_2D_R(Mat, 2*N_,5) ;
            Sol = new double[5] ;
            Vec_ = new double[ 2*N_] ;

            //amrex::Print()<<"cell "<<idata.cellid_[0]<<" , "<<idata.cellid_[1]<<" , "<<Num_<<" , "<<Num_i_<<'\n';

            for(i_ = 0 ; i_ < N_ ; i_++)
            {
                j_ = i_ ;

                Mat[i_][0] = w_[j_] ; // 1 
                Mat[i_][1] = w_[j_]*x_[j_] ; // x
                Mat[i_][2] = w_[j_]*y_[j_] ; // y
                Mat[i_][3] = 0.0 ; // y
                Mat[i_][4] = 0.0 ; // y

                Mat[i_ + N_][0] = 0.0 ; // 1 
                Mat[i_ + N_][1] = -1.0*w_[j_] * y_[j_] ; // x
                Mat[i_ + N_][2] = 0.0 ; // y
                Mat[i_ + N_][3] = w_[j_] ; // y
                Mat[i_ + N_][4] = w_[j_] * x_[j_] ; // y


                Vec_[i_] = w_[j_]*F11_[j_] ; // F11
                Vec_[i_ + N_] = w_[j_]*F21_[j_] ; // F11
            }
            QRdcmp<double> QR(Mat, 2*N_ , 5);
            QR.solve(Vec_, Sol);

            /// P(x,y) = a + bx + cy
            /// @ cut cell the (x,y) = (0,0) as all x_, y_ are computed wrt cut cell
            /// so P = a @ cut cell
            F11(cellid_) = Sol[0];
            F21(cellid_) = Sol[3];


            for(i_ = 0 ; i_ <  2*N_  ;i_++) { delete [] Mat[i_] ; }
                delete [] Mat ;
                delete [] Sol ; delete [] Vec_ ;
        }

        void LS_QR_Second_Order
        (
            bool& sing,
            int Num_,
            double* x_, double* y_,
            double* T_,
            double* w_,
            amrex::Array4<amrex::Real> const& T,
            const amrex::IntVect& cellid_, 
            int iscalar
        )
        {
	    sing = false;
            int i_, j_ ;
            double **Mat, *Sol, *Vec_  ;
	    double **R;

            Allocate_2D_R(Mat,Num_ ,3) ;
	    Allocate_2D_R(R, 3, 3);
            Sol = new double[3] ;
            Vec_ = new double[Num_ ] ;

            for(i_ = 0 ; i_ < Num_ ; i_++)
            {
                Mat[i_][0] = w_[i_] ; // 1 
                Mat[i_][1] = w_[i_]*x_[i_] ; // x
                Mat[i_][2] = w_[i_]*y_[i_] ; // y
                Vec_[i_] = w_[i_]*T_[i_] ; // P
            }
            QRdcmp<double> QR(Mat, Num_, 3);
            QR.solve(Vec_, Sol);
	    QR.get_R(R);

            /// P(x,y) = a + bx + cy
            /// @ cut cell the (x,y) = (0,0) as all x_, y_ are computed wrt cut cell
            /// so P = a @ cut cell
            T(cellid_) = Sol[0];

            if(std::isnan(Sol[0])) sing = true;
            for(i_ = 0 ; i_ < 3; i_++)
            {
                if(std::abs(R[i_][i_]) < 1.0e-6) sing = true;
            }

            for(i_ = 0 ; i_ < Num_;i_++) { delete [] Mat[i_] ; }
                delete [] Mat ;
                delete [] Sol ; delete [] Vec_ ;

            for (i_ = 0; i_ < 3; i_++)
            {
                delete[] R[i_];
            }
            delete[] R;
        }

        void LS_QR_Third_Order
        (
            bool& sing,
            int Num_,
            double* x_, double* y_,
            double* T_,
            double* w_,
            amrex::Array4<amrex::Real> const& T,
            const amrex::IntVect& cellid_, 
            int iscalar
        )
        {
	    sing = false;
            int i_, j_ ;
            double **Mat, *Sol, *Vec_  ;
	    double **R;

            //amrex::Print()<<"Num_ = "<<Num_<<", num_i_ = "<<Num_i_<<'\n';

            Allocate_2D_R(Mat,Num_,6);
	    Allocate_2D_R(R, 6, 6);
            Sol = new double[6];

            Vec_ = new double[Num_] ;

            for(i_ = 0 ; i_ < Num_ ; i_++)
            {
                Mat[i_][0] = w_[i_] ; // 1
                Mat[i_][1] = w_[i_]*x_[i_] ; // x
                Mat[i_][2] = w_[i_]*y_[i_] ; // y
                Mat[i_][3] = w_[i_]*x_[i_]*x_[i_] ; // x^2
                Mat[i_][4] = w_[i_]*x_[i_]*y_[i_] ; // xy
                Mat[i_][5] = w_[i_]*y_[i_]*y_[i_] ; // y^2
                Vec_[i_] = w_[i_]*T_[i_] ; // x
            }
            QRdcmp<double> QR(Mat, Num_, 6);

            QR.solve(Vec_, Sol);
            QR.get_R(R);

            T(cellid_) = Sol[0];
            if(std::isnan(Sol[0])) sing = true;
            for(i_ = 0 ; i_ < 6; i_++)
            {
                if(std::abs(R[i_][i_]) < 1.0e-6) sing = true;
            }

            for(i_ = 0 ; i_ < Num_ ;i_++) { delete [] Mat[i_] ; }
                delete [] Mat ;
                delete [] Sol ; delete [] Vec_ ;

            for (i_ = 0; i_ < 6; i_++)
            {
                delete[] R[i_];
            }
            delete[] R;
        }
     
        void LS_QR_Third_Order_DivFr
        (
            bool& sing,
            int Num_,
            double* x_, double* y_,
            double* F1_,
            double* F2_,
            double* w_,
            amrex::Array4<amrex::Real> const& F1,
            amrex::Array4<amrex::Real> const& F2,
	    const amrex::IntVect& cellid_,
            int iscalar1, int iscalar2
        )
        {
            int i_, j_, N_;
            double **Mat, *Sol, *Vec_  ;
	    //double **Q, **R;

	    //amrex::Print()<<"Num_ = "<<Num_<<", num_i_ = "<<Num_i_<<'\n';

            Allocate_2D_R(Mat,2 * Num_ ,9) ; 
	    //Allocate_2D_R(Q, Num_ - Num_i_, 6) ;
	    //Allocate_2D_R(R, 6, 6) ;
	    Sol = new double[9] ; 

            //Allocate_2D_R(Mat,Num_ - Num_i_,3) ;
	    //Allocate_2D_R(Q, Num_ - Num_i_, 3) ;
            //Allocate_2D_R(R, 3, 3) ;
            //Sol = new double[3] ;

	    Vec_ = new double[2 * Num_] ;

            for(i_ = 0; i_ < Num_ ; i_++)
            {
		//amrex::Print()<<"i_ = "<<i_<<'\n';
                amrex::Real dx = x_[i_];
                amrex::Real dy = y_[i_];	

	            Mat[i_][0] = w_[i_] ;
                Mat[i_][1] = w_[i_]*dx;  // x
                Mat[i_][2] = w_[i_]*dy;  // y
                Mat[i_][3] = w_[i_]*dx * dx;
                Mat[i_][4] = w_[i_]*dx * dy;
                Mat[i_][5] = w_[i_]*dy * dy;
                Mat[i_][6] = 0.0;
                Mat[i_][7] = 0.0;
                Mat[i_][8] = 0.0;

                Mat[i_ + Num_][0] = 0.0; // w
                Mat[i_ + Num_][1] = -w_[i_]*dy; // x
                Mat[i_ + Num_][2] = 0.0; // y
                Mat[i_ + Num_][3] = -w_[i_]*2.0 * dx * dy;
                Mat[i_ + Num_][4] = -w_[i_]*0.5 * dy * dy;
                Mat[i_ + Num_][5] = 0.0;
                Mat[i_ + Num_][6] = w_[i_]*1.0;
                Mat[i_ + Num_][7] = w_[i_]*dx;
                Mat[i_ + Num_][8] = w_[i_]*dx * dx;

                Vec_[i_] = w_[i_]*F1_[i_] ; // x
		        Vec_[i_ + Num_] = w_[i_]*F2_[i_] ; // x
            }
            QRdcmp<double> QR(Mat, 2 * Num_ , 9);

	    //QRdcmp<double> QR(Mat, Num_ - Num_i_, 3);
            QR.solve(Vec_, Sol);
	    //QR.get_Q(Q);
            //QR.get_R(R);

            /// P(x,y) = a + bx + cy
            /// @ cut cell the (x,y) = (0,0) as all x_, y_ are computed wrt cut cell
            /// so P = a @ cut cell
            F1(cellid_) = Sol[0];
    	    F2(cellid_) = Sol[6];

            for(i_ = 0 ; i_ < 2 * Num_  ;i_++) { delete [] Mat[i_] ; }
                delete [] Mat ;
                delete [] Sol ; delete [] Vec_ ;
        }  
        
        void LS_QR_Fourth_Order
        (
            bool& sing,
            int Num_,
            double* x_, double* y_,
            double* T_,
            double* w_,
            amrex::Array4<amrex::Real> const& T,
            const amrex::IntVect& cellid_, 
            int iscalar
        )
        {
            int i_, j_ ;
            double **Mat, *Sol, *Vec_  ;

            //amrex::Print()<<"Num_ = "<<Num_<<", num_i_ = "<<Num_i_<<'\n';

            Allocate_2D_R(Mat,Num_,10) ;
            Sol = new double[10] ;

            Vec_ = new double[Num_] ;

            for(i_ = 0 ; i_ < Num_ ; i_++)
            {
                Mat[i_][0] = w_[i_] ; // 1
                Mat[i_][1] = w_[i_]*x_[i_] ; // x
                Mat[i_][2] = w_[i_]*y_[i_] ; // y
                Mat[i_][3] = w_[i_]*x_[i_]*x_[i_] ; // x^2
                Mat[i_][4] = w_[i_]*x_[i_]*y_[i_] ; // xy
                Mat[i_][5] = w_[i_]*y_[i_]*y_[i_] ; // y^2
                Mat[i_][6] = w_[i_]*x_[i_]*x_[i_]*x_[i_] ;
                Mat[i_][7] = w_[i_]*x_[i_]*x_[i_]*y_[i_] ;
                Mat[i_][8] = w_[i_]*x_[i_]*y_[i_]*y_[i_] ;
                Mat[i_][9] = w_[i_]*y_[i_]*y_[i_]*y_[i_] ;
                Vec_[i_] = w_[i_]*T_[i_] ; // x
            }
            QRdcmp<double> QR(Mat, Num_ , 10);

            QR.solve(Vec_, Sol);

            T(cellid_) = Sol[0];

            for(i_ = 0 ; i_ < Num_ ;i_++) { delete [] Mat[i_] ; }
                delete [] Mat ;
                delete [] Sol ; delete [] Vec_ ;
        }


        namespace Axisymmetric
        {
            void Fij_LSQ_Parameters(
                double &F11_,double &F12_,
                double &F21_,double &F22_,
                double &F33_,
                double &x_,
                double &y_,
                double &w_,
                int i1, int j1,
                int i, int j,
                amrex::Array4<amrex::Real const> const &F11, /// cc vel with 2 comp
                amrex::Array4<amrex::Real const> const &F12,
                amrex::Array4<amrex::Real const> const &F21,
                amrex::Array4<amrex::Real const> const &F22,
                amrex::Array4<amrex::Real const> const &F33,
                const amrex::Real *prob_lo,
                const amrex::Real *dx,
                const amrex::Box &domain)
            {
                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                amrex::Real x1 = prob_lo[0] + dx[0] * (i1 + 0.5);
                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                amrex::Real y1 = prob_lo[1] + dx[1] * (j1 + 0.5);

                x_ = (x1 - x) / dx[0];
                y_ = (y1 - y) / dx[1];
                F11_ = F11(i1, j1, 0);
                F12_ = F12(i1, j1, 0);
                F21_ = F21(i1, j1, 0);
                F22_ = F22(i1, j1, 0);
                F33_ = F33(i1, j1, 0);
                if(j1 < domain.smallEnd(1))
                {
                    //amrex::Print()<<"i = "<<i<<" , j = "<<j<<" , i1 = "<<i1<<" , j1 = "<<j1<<'\n';
                    j1 = domain.smallEnd(1) - j1 - 1;
                    //amrex::Print()<<"new j1 = "<<j1<<'\n';
		    F11_ = F11(i1, j1, 0);
                    F12_ = -1.0*F12(i1, j1, 0);
                    F21_ = -1.0*F21(i1, j1, 0);
                    F22_ = F22(i1, j1, 0);
                    F33_ = F33(i1, j1, 0);
                }

                Compute_LSQ_Weights(w_, std::sqrt(x_ * x_ + y_ * y_));
            }
            
            void RefConfg_LSQ_Parameters(
                double &X_,
                double &Y_,
                double &x_,
                double &y_,
                double &w_,
                int i1, int j1,
                int i, int j,
                amrex::Array4<amrex::Real const> const &X_ref, 
                amrex::Array4<amrex::Real const> const &Y_ref,
                const amrex::Real *prob_lo,
                const amrex::Real *dx,
                const amrex::Box &domain)
            {
                amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
                amrex::Real x1 = prob_lo[0] + dx[0] * (i1 + 0.5);
                amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
                amrex::Real y1 = prob_lo[1] + dx[1] * (j1 + 0.5);

                x_ = (x1 - x) / dx[0];
                y_ = (y1 - y) / dx[1];
                X_ = X_ref(i1, j1, 0, 0);
                Y_ = Y_ref(i1, j1, 0, 1);
                if(j1 < domain.smallEnd(1))
                {
                    //amrex::Print()<<"i = "<<i<<" , j = "<<j<<" , i1 = "<<i1<<" , j1 = "<<j1<<'\n';
                    j1 = domain.smallEnd(1) - j1 - 1;
                    //amrex::Print()<<"new j1 = "<<j1<<'\n';
                    X_ = X_ref(i1, j1, 0, 0);
                    Y_ = -1.0*Y_ref(i1, j1, 0, 1);
                }

                Compute_LSQ_Weights(w_, std::sqrt(x_ * x_ + y_ * y_));
            }

            void QR_Second_Order_DivF_z
            (
                bool& sing,
                int Num_,
                int Num_i_,
		const amrex::IntVect &cell,
                double* x_, double* y_,
                double* F11_,
		double* F21_,
                double* w_,
                amrex::Array4<amrex::Real> const& F11,
		amrex::Array4<amrex::Real> const& F21,
                InterceptData& idata,
	        const amrex::Real* deltax,
                const amrex::Real* prob_lo
            )
            {
                int i_, j_;
                double **Mat; // Least squares and constraint matrices
                double *Sol, *Vec_;
                double dx, dy, dF11, dF21, nx, ny, y0;

                Allocate_2D_R(Mat, 2 * (Num_ - Num_i_), 3);
                Sol = new double[3];
                Vec_ = new double[2 * (Num_ - Num_i_)];
	
                //Allocate_2D_R(Mat,  (Num_ - Num_i_), 4);
                //Sol = new double[4];
                //Vec_ = new double[ (Num_ - Num_i_)];

                amrex::Real y = prob_lo[1] + deltax[1] * (cell[1] + 0.5);
                y0 = y / deltax[1];

                for (i_ = Num_i_; i_ < Num_; i_++)
                {
                    dx = x_[i_];
                    dy = y_[i_];
                    dF11 = F11_[i_];
                    dF21 = F21_[i_];

		    j_ = i_ - Num_i_;
                    Mat[j_][0] = w_[i_];      // w
                    Mat[j_][1] = w_[i_] * dx; // x
                    Mat[j_][2] = w_[i_] * dy; // y
		    //Mat[j_][3] = -1.0 * w_[i_] * dx * dx; // y

                    Mat[j_ + Num_ - Num_i_][0] = 0.0;
                    Mat[j_ + Num_ - Num_i_][1] = -w_[i_] * (y0 + dy) / 2.0;
                    Mat[j_ + Num_ - Num_i_][2] = 0.0;
		    //Mat[j_ + Num_ - Num_i_][3] = w_[i_] * dx * (y0 + dy) * (y0 + dy);


                    Vec_[j_] = w_[i_] * dF11;        // x
                    Vec_[j_ + Num_ - Num_i_] = w_[i_] * dF21; // x
                 }
                 QRdcmp<double> QR(Mat, 2 * (Num_ - Num_i_), 3);
		 //QRdcmp<double> QR(Mat, (Num_ - Num_i_), 4);
                 QR.solve(Vec_, Sol);

                 F11(cell) = Sol[0];
                 F21(cell) = -0.5 * y0 * Sol[1];
	
                 for (i_ = 0; i_ < Num_i_; i_++)
                 {       
                        idata.F11_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_];
                        idata.F21_Int[i_] = -0.5 * Sol[1] * (y0 + y_[i_]);
                 }

		 for (i_ = 0; i_ < 2 * (Num_ - Num_i_); i_++)
                 {
                      delete[] Mat[i_];
                 }
                 delete[] Mat;
                 delete[] Sol;
                 delete[] Vec_;
            }
           
            void QR_Second_OrderF12F22_axisymmetric
            (
                bool& sing,
                int Num_,
                int Num_i_,
	        const amrex::IntVect &cell,
                double* x_, double* y_,
                double* F12_,
	        double* F22_,
                double* w_,
                amrex::Array4<amrex::Real> const& F12,
	        amrex::Array4<amrex::Real> const& F22,
                InterceptData& idata,
                const amrex::Real* deltax,
                const amrex::Real* prob_lo
            )
            {
                int i_, j_;
                double **Mat; // Least squares and constraint matrices
                double *Sol, *Vec_;
                double dx, dy, dF12, dF22, nx, ny, y0;

                Allocate_2D_R(Mat, 2 * (Num_ - Num_i_), 3);
                Sol = new double[3];
                Vec_ = new double[2 * (Num_ - Num_i_)];

                amrex::Real y = prob_lo[1] + deltax[1] * (cell[1] + 0.5);
                y0 = y / deltax[1];

                for (i_ = Num_i_; i_ < Num_; i_++)
                {
                    dx = x_[i_];
                    dy = y_[i_];
                    dF12 = F12_[i_];
                    dF22 = F22_[i_];

		            j_ = i_ - Num_i_;
                    Mat[j_][0] = w_[i_];      // w
                    Mat[j_][1] = w_[i_] * dx; // x
                    Mat[j_][2] = w_[i_] * dy; // y

                    Mat[j_ + Num_ - Num_i_][0] = 0.0;
                    Mat[j_ + Num_ - Num_i_][1] = -w_[i_] * (y0 + dy) / 2.0;
                    Mat[j_ + Num_ - Num_i_][2] = 0.0;

                    Vec_[j_] = w_[i_] * dF12;        // x
                    Vec_[j_ + Num_ - Num_i_] = w_[i_] * dF22; // x
                 }
                 QRdcmp<double> QR(Mat, 2 * (Num_ - Num_i_), 3);
                 QR.solve(Vec_, Sol);

                 F12(cell) = Sol[0];
                 F22(cell) = -0.5 * y0 * Sol[1];
	
                 for (i_ = 0; i_ < Num_i_; i_++)
                 {       
                        idata.F12_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_];
                        idata.F22_Int[i_] = -0.5 * Sol[1] * (y0 + y_[i_]);
                 }

                 for (i_ = 0; i_ < 2 * (Num_ - Num_i_); i_++)
                 {
                      delete[] Mat[i_];
                 }
                 delete[] Mat;
                 delete[] Sol;
                 delete[] Vec_;
            }
         
            void QR_Second_OrderF12F22F33
            (
                bool& sing,
                int Num_,
                int Num_i_,
                const amrex::IntVect &cell,
                double* x_, double* y_,
                double* F12_,
                double* F22_,
                double* F33_,
                double* w_,
                amrex::Array4<amrex::Real> const& F12,
                amrex::Array4<amrex::Real> const& F22,
                amrex::Array4<amrex::Real> const& F33,
                InterceptData& idata,
                const amrex::Real* deltax,
                const amrex::Real* prob_lo
            )
            {
		
                int i_, j_;
                double **Mat; // Least squares and constraint matrices
                double *Sol, *Vec_;
                double dx, dy, dF12, dF22, dF33, nx, ny, y0;
                int N_ = Num_ - Num_i_;
                Allocate_2D_R(Mat, 3 * N_, 3);
                Sol = new double[3];
                Vec_ = new double[3 * N_];

                amrex::Real y = prob_lo[1] + deltax[1] * (cell[1] + 0.5);
                y0 = y / deltax[1];

                for (i_ = Num_i_; i_ < Num_; i_++)
                {
                    dx = x_[i_];
                    dy = y_[i_];
                    dF12 = F12_[i_];
                    dF22 = F22_[i_];
                    dF33 = F33_[i_];
 
                    j_ = i_ - Num_i_;
                    Mat[j_][0] = w_[i_];      // w
                    Mat[j_][1] = w_[i_] * dx; // x
                    Mat[j_][2] = w_[i_] * dy; // y

		            Mat[j_ + N_][0] = 1.0 * w_[i_];      // w
                    Mat[j_ + N_][1] = -1.0 * w_[i_] * (y0 + dy); // x
                    Mat[j_ + N_][2] = w_[i_] * dx; // y

                    Mat[j_ + 2*N_][0] = w_[i_];      // w
                    Mat[j_ + 2*N_][1] = -1.0 * w_[i_] * (y0 + dy); // x
                    Mat[j_ + 2*N_][2] = w_[i_] * dx; // y

                    Vec_[j_] = w_[i_] * dF12;        // x
                    Vec_[j_ + N_] = w_[i_] * dF22; // x
		            Vec_[j_ + 2*N_] = w_[i_] * dF33;
                 }
                 QRdcmp<double> QR(Mat, 3 * N_ , 3);
                 QR.solve(Vec_, Sol);

                 F12(cell) = Sol[0];
                 F22(cell) = Sol[0] - y0 * Sol[1];
		         F33(cell) = Sol[0] - y0 * Sol[1];
	
                 for (i_ = 0; i_ < Num_i_; i_++)
                 {       
                    idata.F12_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * (y_[i_]);
			        idata.F22_Int[i_] = Sol[0] - Sol[1] * (y0 + y_[i_]) + Sol[2] * x_[i_];
		            idata.F33_Int[i_] = Sol[0] - Sol[1] * (y0 + y_[i_]) + Sol[2] * x_[i_];
                    //idata.F22_Int[i_] = 3.0 * x_[i_] * Sol[0] + Sol[1] * (y0 + y_[i_]);
			        //idata.F33_Int[i_] = x_[i_] * Sol[0] + Sol[1] * (y0 + y_[i_]);
                 }

                 for (i_ = 0; i_ < 3 * N_ ; i_++)
                 {
                      delete[] Mat[i_];
                 }
                 delete[] Mat;
                 delete[] Sol;
                 delete[] Vec_;
                /*	
                int i_, j_;
                double **Mat; // Least squares and constraint matrices
                double *Sol, *Vec_;
                double dx, dy, dF12, dF22, dF33, nx, ny, y0;
                int N_ = Num_ - Num_i_;
                Allocate_2D_R(Mat, 3 * N_, 3);
                Sol = new double[3];
                Vec_ = new double[3 * N_];

                amrex::Real y = prob_lo[1] + deltax[1] * (cell[1] + 0.5);
                y0 = y / deltax[1];

                for (i_ = Num_i_; i_ < Num_; i_++)
                {
                    dx = x_[i_];
                    dy = y_[i_];
                    dF12 = F12_[i_];
                    dF22 = F22_[i_];
                    dF33 = F33_[i_];

                    j_ = i_ - Num_i_;
                    Mat[j_][0] = w_[i_];      // w
                    Mat[j_][1] = w_[i_] * dx; // x
                    Mat[j_][2] = w_[i_] * dx * (dy + y0); // y

                    Mat[j_ + N_][0] = w_[i_];      // w
                    Mat[j_ + N_][1] = w_[i_] * (y0 + dy); // x
                    Mat[j_ + N_][2] = 0.0; // y

                    Mat[j_ + 2*N_][0] = w_[i_] * 3.0 * (y + y0) * (y + y0);      // w
                    Mat[j_ + 2*N_][1] = w_[i_] * (y0 + dy); // x
                    Mat[j_ + 2*N_][2] = 3.0*w_[i_]; // y

                    Vec_[j_] = w_[i_] * dF12;        // x
                    Vec_[j_ + N_] = w_[i_] * dF22; // x
                    Vec_[j_ + 2*N_] = w_[i_] * dF33;
                 }
                 QRdcmp<double> QR(Mat, 3 * N_ , 3);
                 QR.solve(Vec_, Sol);

                 F12(cell) = Sol[0];
                 F22(cell) = Sol[0] + Sol[1] * y0;
                 F33(cell) = Sol[0] * 3.0 * y0 * y0 + Sol[1] * 3.0 * y0 + + 3.0*Sol[2];

                 for (i_ = 0; i_ < Num_i_; i_++)
                 {
                        idata.F12_Int[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * x_[i_] * (y0 + y_[i_]);
                        idata.F22_Int[i_] = Sol[0] + Sol[1] *(y0 + y_[i_]);
                        idata.F33_Int[i_] = Sol[0] * 3.0 * (y0 + y_[i_]) * (y0 + y_[i_]) + Sol[1] * (y0 + y_[i_]) + 3.0*Sol[2];
                 }

                 for (i_ = 0; i_ < 3 * N_ ; i_++)
                 {
                      delete[] Mat[i_];
                 }
                 delete[] Mat;
                 delete[] Sol;
                 delete[] Vec_;
		 */
            }
        
 
        
        } // namespace Axisymmetric
        
        
    } // namespace FillFij	

    void Mask::LSQ_Parameters_Fij
    (
        double &Fij_,
	    int &iscalar,
        double &x_,
        double &y_,
        double &w_,
        int i1, int j1,
        int i, int j,
        amrex::Array4<amrex::Real const> const &Fij,
        const amrex::Real *prob_lo,
        const amrex::Real *dx,
        const amrex::Box &domain
    )
    {
        amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
        amrex::Real x1 = prob_lo[0] + dx[0] * (i1 + 0.5);
        amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
        amrex::Real y1 = prob_lo[1] + dx[1] * (j1 + 0.5);

        x_ = x1 - x;
        x_ /= dx[0];
        y_ = y1 - y;
        y_ /= dx[1];
        Fij_ = Fij(i1, j1, 0);
        if(j1 < domain.smallEnd(1))
        {
            //amrex::Print()<<"i = "<<i<<" , j = "<<j<<" , i1 = "<<i1<<" , j1 = "<<j1<<'\n';
            j1 = 2*domain.smallEnd(1) - j1 - 1;
            //amrex::Print()<<"new j1 = "<<j1<<'\n';
	        if(iscalar == 3 || iscalar == 4)
                Fij_ = -1.0*Fij(i1, j1, 0);
	        else
                Fij_ = 1.0*Fij(i1, j1, 0);
        }
        Compute_LSQ_Weights(w_, std::sqrt(x_ * x_ + y_ * y_));
    }

    void Mask::FillInGhostFij(int &iscalar, amrex::MultiFab &mf)
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
                    int Num_i_ = 0;
    
                    /// create grown box around cell
                    amrex::Box gbx(idt.cellid_, idt.cellid_);
                    gbx.grow(idt.stencil_);
                    //amrex::Box gbx_isect = gbx & domain;
                    int ni = 0; 
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
                        T_[ni] = 0.0;
                        Compute_LSQ_Weights(w_[ni], std::fabs(idt.frac_[ni]));
                    }
		    //ni = 0;

		            Num_i_ = ni;
                    for (amrex::BoxIterator bit(gbx); bit.ok(); ++bit)
                    {
                        const amrex::IntVect &iv = bit();
                        if (pmask(iv) == 1)
                        {
		                    if(isAxisymmetric)
                                LSQ_Parameters_Fij(T_[ni], iscalar, x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], T, prob_lo, dx, domain);
			                else
		                        LSQ_Parameters(T_[ni], x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], T, prob_lo, dx);
                            ni++;
                        }
                    }
    
                    int Num_ = ni;
    
		    int PORDER_ = Fij_order;
		    //PORDER_ = idt.PORDER;
                    if (PORDER_ >= 3)
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
			    FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, T_, w_, T, idt, iscalar);
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
                            FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, T_, w_, T, idt, iscalar);
                            if (Sing_)
                            {
		                FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, T_, w_, T, idt, iscalar);
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
                    else if (PORDER_ == 2)
                    {
			//amrex::Print()<<"In second order"<<'\n';
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
                        else
                        {
			    FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, T_, w_, T, idt, iscalar);
			    amrex::Real y = prob_lo[1] + dx[1] * (idt.cellid_[1] + 0.5);

			    if(Sing_)
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
		            bool debug = false;
                    if(idt.cellid_[0] == 115 && idt.cellid_[1] == 0 && debug)
                    {
                        amrex::Print()<<"cell : "<<idt.cellid_<<'\n';
                        amrex::Print()<<"iscalar = "<<iscalar<<" , Fij : "<< T(idt.cellid_)<<'\n';
                        amrex::Print()<<"Num_ = "<<Num_<<" , Num_i_ = "<<Num_i_<<'\n';
                        for(int ii = 0;ii< Num_;ii++)
                        {
                            amrex::Print()<<"ii = "<<ii<<"\t";
                            amrex::Print()<<" x_ = "<<x_[ii]<<"\t";
                            amrex::Print()<<" y_ = "<<y_[ii]<<"\t";
                            amrex::Print()<<" T_ = "<<T_[ii]<<"\t";
                            amrex::Print()<<" idt.Interpolation_weights = "<<idt.Interpolation_weights[ii]<<'\t';
                            amrex::Print()<<" w_ = "<<w_[ii]<<"\n";
                        }
                        std::exit(9);
                    }
                }
            }
            nsolid++;
        }
    
        /// as pressure values are updated
        /// fill the internal ghost values
        mf.FillBoundary();
    }      

    void Mask::FillInGhostFij(amrex::MultiFab &mfF11, amrex::MultiFab &mfF12, amrex::MultiFab &mfF21, amrex::MultiFab &mfF22, amrex::MultiFab &mfF33)
    {
        /// fill internal ghost values for pressure
        mfF11.FillBoundary();
        mfF12.FillBoundary();
        mfF21.FillBoundary();
        mfF22.FillBoundary();
	    mfF33.FillBoundary();
        
        int PORDER_ = Fij_order;//Higher order does not work for Fij
        int PORDER_MIN = 2;
        int iscalar = 0;
        bool use_div_free = false;
        {
            if (interfaces->empty())
            {
                return;
            }
            const amrex::Real *prob_lo = geom_.ProbLo();
            const amrex::Real *dx = geom_.CellSize();
            const amrex::Box &domain = geom_.Domain();
		    
            bool Sing_;
            double sum_w, sum_sol_F11, sum_sol_F21,sum_sol_F12,sum_sol_F22,sum_sol_F33;
		    
            int N_Max = (2 * stencil_ + 1) * (2 * stencil_ + 1) + 4;
            double x_[N_Max], y_[N_Max], w_[N_Max], F11_[N_Max], F12_[N_Max], F21_[N_Max], F22_[N_Max], F33_[N_Max];
		    
            for (auto &&solid : *interfaces)
            {
                int i_beg = 0;
                if (!solid->isAdvectLS())
                    return;
                for (amrex::MFIter mfi(mfF11); mfi.isValid(); ++mfi)
                {
                    amrex::Array4<int const> const &pmask = PMask.const_array(mfi);
                    amrex::Array4<amrex::Real> const &F11 = mfF11.array(mfi);
                    amrex::Array4<amrex::Real> const &F12 = mfF12.array(mfi);
                    amrex::Array4<amrex::Real> const &F21 = mfF21.array(mfi);
                    amrex::Array4<amrex::Real> const &F22 = mfF22.array(mfi);
                    amrex::Array4<amrex::Real> const &F33 = mfF33.array(mfi);

                    auto &icpt_data = solid->getInterceptData()[mfi];
                    amrex::IntVect test_iv(AMREX_D_DECL(115, 128, 0));
		    
                    for (auto &&idt : icpt_data)
                    {
                        Sing_ = false; // Need this. Sing_ is only changed to true.
                        int ni = 0;
		        	
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
                                //amrex::PrintToFile("log") << "Error in identifying interface type.\n";
                                exit(1);
                            }
                            //U_[ni] = idt.u[ni];
                            //V_[ni] = idt.v[ni];
                            Compute_LSQ_Weights(w_[ni], std::fabs(idt.frac_[ni]));
                        }
			
                        int Num_i_ = ni;

		                /*	
                        int dist_x1 = abs(idt.cellid_[0] - domain.smallEnd(0));
                        int dist_x2 = abs(idt.cellid_[0] - domain.bigEnd(0));
                        int dist_y1 = abs(idt.cellid_[1] - domain.smallEnd(1));
                        int dist_y2 = abs(idt.cellid_[1] - domain.bigEnd(1));

                        //decide size of the lease-square stencil based on available points
                        if( dist_x1 >= 4 && dist_x2 >= 4 && dist_y1 >= 4 && dist_y2 >= 4)
                        {
                            //idt.stencil_ = 4;
                        }
                        else if( dist_x1 >= 3 && dist_x2 >= 3 && dist_y1 >= 3 && dist_y2 >= 3)
                        {
                            //idt.stencil_ = 3;
                            PORDER_ = 3;
                            //amrex::Print()<<"setting stencil size to 3 "<<idt.cellid_<<'\n';
                        }
                        else if( dist_x1 >= 2 && dist_x2 >= 2 && dist_y1 >= 2 && dist_y2 >= 2)
                        {
                            //idt.stencil_ = 2;
                            PORDER_ = 3;
                            //amrex::Print()<<"setting stencil size to 2 "<<idt.cellid_<<'\n';
                        }
                        else
                        {
                            //idt.stencil_ = 2;
                            PORDER_ = 3;
                            //amrex::Print()<<"setting stencil size to 1 "<<idt.cellid_<<'\n';
                        }
			            */

		    
                        /// create grown box around cell
                        amrex::Box gbx(idt.cellid_, idt.cellid_);
                        gbx.grow(stencil_);
		    
                        amrex::Box gbx_isect = gbx & domain;
		    
                        for (amrex::BoxIterator bit(gbx_isect); bit.ok(); ++bit)
                        {
                            const amrex::IntVect &iv = bit();
                            if (pmask(iv) == 1)
                            {
                                if(isAxisymmetric)
				                {
                                    //FillFij::Axisymmetric::Fij_LSQ_Parameters(F11_[ni], F12_[ni], F21_[ni], F22_[ni], F33_[ni], x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], F11, F12, F21, F22, F33, prob_lo, dx,domain);
			                        {
                                        iscalar = 2;
                                        LSQ_Parameters_Fij(F11_[ni], iscalar, x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], F11, prob_lo, dx, domain);
                                    }
                                    {
                                        iscalar = 3;
                                        LSQ_Parameters_Fij(F12_[ni], iscalar, x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], F12, prob_lo, dx, domain);
                                    }
                                    {
                                        iscalar = 4;
                                        LSQ_Parameters_Fij(F21_[ni], iscalar, x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], F21, prob_lo, dx, domain);
                                    }
                                    {
                                        iscalar = 5;
                                        LSQ_Parameters_Fij(F22_[ni], iscalar, x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], F22, prob_lo, dx, domain);
                                    }
                                    {
                                        iscalar = 6;
                                        LSQ_Parameters_Fij(F33_[ni], iscalar, x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], F33, prob_lo, dx, domain);
                                    }
				                }
                                else
                                    Fij_LSQ_Parameters(F11_[ni], F12_[ni], F21_[ni], F22_[ni], F33_[ni], x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], F11, F12, F21, F22, F33, prob_lo, dx);
                                ni++;
                            }
                            /*if( idt.cellid_ == test_iv)
                            {
                                amrex::Print()<<" ni = "<<ni<<'\n';
                                amrex::Print()<<" iv = "<<iv[0]<<" , "<<iv[1]<<'\n';
                                amrex::Print()<<"in domain = "<<domain.contains(iv)<<'\n';
                            }*/
                        }
		    
                        int Num_ = ni;
		    
                        if (solid->isAdvectLS())
                        {
                            i_beg = Num_i_;
                        }

			            if (PORDER_ == 4)
                        {
                            if (Num_ < 3 + Num_i_)
                            {
                                // too few points have been found...can only do a simple average
                                sum_w = sum_sol_F11 = sum_sol_F12 = sum_sol_F21 = sum_sol_F22 = sum_sol_F33 = 0.0;
                                for (int i = i_beg; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol_F11 += w_[i] * F11_[i];
                                    sum_sol_F12 += w_[i] * F12_[i];
                                    sum_sol_F21 += w_[i] * F21_[i];
                                    sum_sol_F22 += w_[i] * F22_[i];
                                    sum_sol_F33 += w_[i] * F33_[i];
                                }
                                F11(idt.cellid_) = sum_sol_F11 / sum_w;
                                F12(idt.cellid_) = sum_sol_F12 / sum_w;
				                F21(idt.cellid_) = sum_sol_F21 / sum_w;
                                F22(idt.cellid_) = sum_sol_F22 / sum_w;
				                F33(idt.cellid_) = sum_sol_F33 / sum_w;
		    
                                //amrex::PrintToFile("log") << " Inside first \n";
                            }
                            else if(Num_ < 6 + Num_i_)
                            {
                                if (isAxisymmetric)
                                {
                                    FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                    FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                }
                                else
                                {
                                    FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                        
                                }
                                if (Sing_)
                                {
                                    // too few points have been found...can only do a simple average
                                    sum_w = sum_sol_F11 = sum_sol_F12 = sum_sol_F21 = sum_sol_F22 = sum_sol_F33 = 0.0;
                                    for (int i = i_beg; i < Num_; i++)
                                    {
                                        sum_w += w_[i];
                                        sum_sol_F11 += w_[i] * F11_[i];
                                        sum_sol_F12 += w_[i] * F12_[i];
                                        sum_sol_F21 += w_[i] * F21_[i];
                                        sum_sol_F22 += w_[i] * F22_[i];
                                        sum_sol_F33 += w_[i] * F33_[i];
                                    }
                                    F11(idt.cellid_) = sum_sol_F11 / sum_w;
                                    F12(idt.cellid_) = sum_sol_F12 / sum_w;
                                    F21(idt.cellid_) = sum_sol_F21 / sum_w;
                                    F22(idt.cellid_) = sum_sol_F22 / sum_w;
                                    F33(idt.cellid_) = sum_sol_F33 / sum_w;
                                }
                            }
                            else if(Num_ < 10 + Num_i_)
                            {
                                if (isAxisymmetric)
                                {
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                    if(Sing_)
                                    {
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                    }
                                }
                                else
                                {
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);

                                    //FillFij::QR_Third_Order_DivFr(Sing_, Num_, Num_i_, x_, y_, F11_, F21_, w_, F11, F21, idt, 2,4);
                                    //FillFij::QR_Third_Order_DivFr(Sing_, Num_, Num_i_, x_, y_, F12_, F22_, w_, F12, F22, idt, 3,5);
                                    if(Sing_)
                                    {
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                    }
                                }
                            }
                            else
                            {
		                		//amrex::Print()<<"4th order gfm"<<'\n';
                                if (isAxisymmetric)
                                {
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                    //
                                    FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                    FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                    if(Sing_)
                                    {
					                    Sing_ = false;
                                        FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                        FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                    }
                                    if(Sing_)
                                    {
					                    Sing_ = false;
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                    }
                                }
                                else
                                {
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                    //FillFij::QR_Third_Order_DivFr(Sing_, Num_, Num_i_, x_, y_, F11_, F21_, w_, F11, F21, idt, 2,4);
                                    //FillFij::QR_Third_Order_DivFr(Sing_, Num_, Num_i_, x_, y_, F12_, F22_, w_, F12, F22, idt, 3,5);

                                    FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                    //FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                    if(Sing_)
                                    {
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);

                                    }
                                }
                            }
                        }
                        else if (PORDER_ == 3)
                        {
                            if (Num_ < 3 + Num_i_)
                            {
                                // too few points have been found...can only do a simple average
                                sum_w = sum_sol_F11 = sum_sol_F12 = sum_sol_F21 = sum_sol_F22 = sum_sol_F33 = 0.0;
                                for (int i = i_beg; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol_F11 += w_[i] * F11_[i];
                                    sum_sol_F12 += w_[i] * F12_[i];
                                    sum_sol_F21 += w_[i] * F21_[i];
                                    sum_sol_F22 += w_[i] * F22_[i];
                                    sum_sol_F33 += w_[i] * F33_[i];
                                }
                                F11(idt.cellid_) = sum_sol_F11 / sum_w;
                                F12(idt.cellid_) = sum_sol_F12 / sum_w;
                				F21(idt.cellid_) = sum_sol_F21 / sum_w;
                                F22(idt.cellid_) = sum_sol_F22 / sum_w;
				                F33(idt.cellid_) = sum_sol_F33 / sum_w;
		    
                                //amrex::PrintToFile("log") << " Inside first \n";
                            }
			                else if(Num_ < 6 + Num_i_)
		                    {
                                if (isAxisymmetric)
			                    {
                                     FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                     FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                     FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                     FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                     FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                     if(std::isnan(F11(idt.cellid_))||
                                        std::isnan(F12(idt.cellid_))||
                                        std::isnan(F21(idt.cellid_))||
                                        std::isnan(F22(idt.cellid_))||
                                        std::isnan(F33(idt.cellid_)))
                                         Sing_ = true;
                                     //else
                                     //    Sing_ = false;
				                }
                                else
                                {
                                    FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                    if(std::isnan(F11(idt.cellid_))||
                                    std::isnan(F12(idt.cellid_))||
                                    std::isnan(F21(idt.cellid_))||
                                    std::isnan(F22(idt.cellid_)))
                                        Sing_ = true;
                                    else
                                        Sing_ = false;
                                }
                                if (Sing_)
                                {
                                    // too few points have been found...can only do a simple average
                                    sum_w = sum_sol_F11 = sum_sol_F12 = sum_sol_F21 = sum_sol_F22 = sum_sol_F33 = 0.0;
                                    for (int i = i_beg; i < Num_; i++)
                                    {
                                        sum_w += w_[i];
                                        sum_sol_F11 += w_[i] * F11_[i];
                                        sum_sol_F12 += w_[i] * F12_[i];
                                        sum_sol_F21 += w_[i] * F21_[i];
                                        sum_sol_F22 += w_[i] * F22_[i];
                                        sum_sol_F33 += w_[i] * F33_[i];
                                    }
                                    F11(idt.cellid_) = sum_sol_F11 / sum_w;
                                    F12(idt.cellid_) = sum_sol_F12 / sum_w;
                                    F21(idt.cellid_) = sum_sol_F21 / sum_w;
                                    F22(idt.cellid_) = sum_sol_F22 / sum_w;
                                    F33(idt.cellid_) = sum_sol_F33 / sum_w;
                                }
                            }   
                            else
                            {
                            if (isAxisymmetric)
                            {
                                FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                if(std::isnan(F11(idt.cellid_))||
                                    std::isnan(F12(idt.cellid_))||
                                    std::isnan(F21(idt.cellid_))||
                                    std::isnan(F22(idt.cellid_))||
                                    std::isnan(F33(idt.cellid_)))
                                    Sing_ = true;
                                    if(Sing_)
                                    {
					                    amrex::Print(-1)<<"Dorpping F_ij extension order to 2nd"<<'\n';
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                        if(std::isnan(F11(idt.cellid_))||
                                           std::isnan(F12(idt.cellid_))||
                                           std::isnan(F21(idt.cellid_))||
                                           std::isnan(F22(idt.cellid_))||
                                           std::isnan(F33(idt.cellid_)))
                                           Sing_ = true;
                                        
                                        if(Sing_)
                                        {
					                        amrex::Print(-1)<<"Dorpping F_ij extension order to 1st"<<'\n';
                                            sum_w = sum_sol_F11 = sum_sol_F12 = sum_sol_F21 = sum_sol_F22 = sum_sol_F33 = 0.0;
                                            for (int i = i_beg; i < Num_; i++)
                                            {
                                                sum_w += w_[i];
                                                sum_sol_F11 += w_[i] * F11_[i];
                                                sum_sol_F12 += w_[i] * F12_[i];
                                                sum_sol_F21 += w_[i] * F21_[i];
                                                sum_sol_F22 += w_[i] * F22_[i];
                                                sum_sol_F33 += w_[i] * F33_[i];
                                            }
                                            F11(idt.cellid_) = sum_sol_F11 / sum_w;
                                            F12(idt.cellid_) = sum_sol_F12 / sum_w;
                                            F21(idt.cellid_) = sum_sol_F21 / sum_w;
                                            F22(idt.cellid_) = sum_sol_F22 / sum_w;
                                            F33(idt.cellid_) = sum_sol_F33 / sum_w;
                                        }
                                    }
                                }
                                else
                                {
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);

				    if(Sing_)
				    {
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
				    }
				}
                            }
                        }
			else if(PORDER_ == 2)
			{
                //amrex::Print()<<"Second order gfm"<<'\n';
                if (Num_ < 3 + Num_i_)
                {
                    // too few points have been found...can only do a simple average
                    sum_w = sum_sol_F11 = sum_sol_F12 = sum_sol_F21 = sum_sol_F22 = sum_sol_F33 = 0.0;
                    for (int i = i_beg; i < Num_; i++)
                    {
                        sum_w += w_[i];
                        sum_sol_F11 += w_[i] * F11_[i];
                        sum_sol_F12 += w_[i] * F12_[i];
                        sum_sol_F21 += w_[i] * F21_[i];
                        sum_sol_F22 += w_[i] * F22_[i];
                        sum_sol_F33 += w_[i] * F33_[i];
                    }
                    F11(idt.cellid_) = sum_sol_F11 / sum_w;
                    F12(idt.cellid_) = sum_sol_F12 / sum_w;
                    F21(idt.cellid_) = sum_sol_F21 / sum_w;
                    F22(idt.cellid_) = sum_sol_F22 / sum_w;
                    F33(idt.cellid_) = sum_sol_F33 / sum_w;

                    //amrex::PrintToFile("log") << " Inside first \n";
                }
                else
                {
                    if (isAxisymmetric)
                    {
        //FillFij::Axisymmetric::QR_Second_Order_DivF_z(Sing_, Num_, Num_i_, idt.cellid_, x_, y_, F11_, F21_, w_, F11, F21, idt, dx, prob_lo);
                            FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                            FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                            FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                            FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                            FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                    }
                    else
                    {
                            FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                            FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                            FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                            FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);

                    }
                }
			}
                if(std::isnan(F11(idt.cellid_))||
                    std::isnan(F12(idt.cellid_))||
                    std::isnan(F21(idt.cellid_))||
                    std::isnan(F22(idt.cellid_))||
                    std::isnan(F33(idt.cellid_)))
                {
                    amrex::Print(-1)<<"Found nan in extension of F_ij"<<'\n';
                    amrex::Print()<<"cell : "<<idt.cellid_<<'\n';
                    amrex::Print()<<"Num_ = "<<Num_<<" , Num_i_ = "<<Num_i_<<'\n';
                    amrex::Print()<<"F11(idt.cellid_) : "<< F11(idt.cellid_)<<'\n';
                    amrex::Print()<<"F12(idt.cellid_) : "<< F12(idt.cellid_)<<'\n';
                    amrex::Print()<<"F21(idt.cellid_) : "<< F21(idt.cellid_)<<'\n';
                    amrex::Print()<<"F22(idt.cellid_) : "<< F22(idt.cellid_)<<'\n';
                    amrex::Print()<<"F33(idt.cellid_) : "<< F33(idt.cellid_)<<'\n';
                    for(int ii = 0; ii<Num_;ii++)
                    {
                        amrex::Print(-1)<<"ii = "<<ii<<'\n';
                        amrex::Print(-1)<<"x_[ii] = "<<x_[ii]<<" , y_[ii] = "<<y_[ii]<<'\n';
                        amrex::Print(-1)<<"F11_[ii] = "<<F11_[ii]<<'\n';
                    }
                    exit(9);
                }
                        /* 
                        amrex::Real x_gc = prob_lo[0] + dx[0] * (idt.cellid_[0] + 0.5);
                        amrex::Real y_gc = prob_lo[1] + dx[1] * (idt.cellid_[1] + 0.5);
                        amrex::Real r_gc = std::hypot(x_gc - solid->Xcp(), y_gc - solid->Ycp());
                        amrex::Real theta = std::acos((x_gc - solid->Xcp())/r_gc);
			if(y_gc < 0)
		            theta = 2*M_PI - theta; 
			amrex::Real R_by_r,R;
			if(isAxisymmetric)
			{
		            //amrex::Real R = std::pow(solid->Volume_0()/(4.0/3.0*M_PI),(1.0/3.0));
                    //R_by_r = std::pow(solid->Volume_0()/(4.0/3.0*M_PI),(1.0/3.0))/r_gc;///std::pow((solid->Volume_0() / solid->Volume()), 1.0/3.0); 
                    amrex::Real R = std::pow(r_gc*r_gc*r_gc + solid->Time(),1.0/3.0);
                    amrex::Real R_by_r = R/r_gc;
                    F11(idt.cellid_) = R_by_r*R_by_r*std::pow(std::cos(theta),2.0) + 1.0/R_by_r * std::pow(std::sin(theta),2.0);
                    F12(idt.cellid_) = (R_by_r*R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);
                    F21(idt.cellid_) = (R_by_r*R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);
                    F22(idt.cellid_) = R_by_r*R_by_r*std::pow(std::sin(theta),2.0) + 1.0/R_by_r * std::pow(std::cos(theta),2.0);
                    F33(idt.cellid_) = 1.0/R_by_r;
			}
			else
			{
			    amrex::Real R = std::pow(r_gc*r_gc + solid->Time(),1.0/2.0);
                            amrex::Real R_by_r = R/r_gc;

                            F11(idt.cellid_) = R_by_r*std::pow(std::cos(theta),2.0) + 1.0/R_by_r * std::pow(std::sin(theta),2.0);
                            F12(idt.cellid_) = (R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);;
                            F21(idt.cellid_) = (R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);;
                            F22(idt.cellid_) = R_by_r*std::pow(std::sin(theta),2.0) + 1.0/R_by_r * std::pow(std::cos(theta),2.0);
			}
                        for (int i_ = 0; i_ < idt.n_intercepts; i_++)
                        {
                            idt.F11_Int[i_] = F11(idt.cellid_);
                            idt.F12_Int[i_] = F12(idt.cellid_);
                            idt.F21_Int[i_] = F21(idt.cellid_);
                            idt.F22_Int[i_] = F22(idt.cellid_);
                            idt.F33_Int[i_] = F33(idt.cellid_);
                        }
                        */		
                        //Vel(idt.cellid_, 0) = 0.0;//sum_sol_u / sum_w;
                        //Vel(idt.cellid_, 1) = 0.0;//sum_sol_v / sum_w;
			/*
                        if(idt.cellid_ == test_iv)
                        {
                                amrex::Print()<<"cell : "<<idt.cellid_<<'\n';
                                amrex::Print()<<"Vel(idt.cellid_) : "<< Vel(idt.cellid_, 0)<<" , "<<Vel(idt.cellid_, 1)<<'\n';
                                amrex::Print()<<"Num_ = "<<Num_<<" , Num_i_ = "<<Num_i_<<'\n';
                                //for (ni = 0; ni < idt.n_intercepts; ni++)
                                //{
				    
                                    //amrex::Print()<<"idt.norm_shear_["<< ni <<"] = "<<idt.norm_shear_[ni]<<'\n';
                                //}
                                for(int ii = 0;ii<= Num_;ii++)
                                {
                                        amrex::Print()<<"ii = "<<ii<<"\t";
                                        amrex::Print()<<"x_ = "<<x_[ii]<<"\t";
                                        amrex::Print()<<"y_ = "<<y_[ii]<<"\n";
                                        amrex::Print()<<"U_(x_, 0) = "<<U_[ii]<<"\n";
				        amrex::Print()<<"V_(y_, 1) = "<<V_[ii]<<"\n";
                                }
                                //std::exit(9);
                        }
		        */
                    }    
                }
            }
	}
    }

    void Mask::FillInGhostFij(amrex::MultiFab &mfF11, amrex::MultiFab &mfF12, amrex::MultiFab &mfF21, amrex::MultiFab &mfF22, amrex::MultiFab &mfF33, amrex::MultiFab &mfPhi)
    {
        /// fill internal ghost values for pressure
        mfF11.FillBoundary();
        mfF12.FillBoundary();
        mfF21.FillBoundary();
        mfF22.FillBoundary();
	    mfF33.FillBoundary();
        
        int PORDER_ = 1;//Set first order for all poitns, later decide if higher 
                        // order can be used based on the least-square cloud 
        int PORDER_MIN = 2;
        int iscalar = 0;
        bool use_div_free = false;
        {
            if (interfaces->empty())
            {
                return;
            }
            const amrex::Real *prob_lo = geom_.ProbLo();
            const amrex::Real *dx = geom_.CellSize();
            const amrex::Box &domain = geom_.Domain();
		    
            bool Sing_;
            double sum_w, sum_sol_F11, sum_sol_F21,sum_sol_F12,sum_sol_F22,sum_sol_F33;
		    
            int N_Max = (2 * stencil_ + 1) * (2 * stencil_ + 1) + 4;
            double x_[N_Max], y_[N_Max], w_[N_Max], F11_[N_Max], F12_[N_Max], F21_[N_Max], F22_[N_Max], F33_[N_Max];
		    
            for (auto &&solid : *interfaces)
            {
                int i_beg = 0;
                if (!solid->isAdvectLS())
                    return;
                for (amrex::MFIter mfi(mfF11); mfi.isValid(); ++mfi)
                {
                    amrex::Array4<int const> const &pmask = PMask.const_array(mfi);
                    amrex::Array4<amrex::Real> const &F11 = mfF11.array(mfi);
                    amrex::Array4<amrex::Real> const &F12 = mfF12.array(mfi);
                    amrex::Array4<amrex::Real> const &F21 = mfF21.array(mfi);
                    amrex::Array4<amrex::Real> const &F22 = mfF22.array(mfi);
                    amrex::Array4<amrex::Real> const &F33 = mfF33.array(mfi);
                    amrex::Array4<amrex::Real> const &phi = mfPhi.array(mfi);

                    auto &icpt_data = solid->getInterceptData()[mfi];
                    amrex::IntVect test_iv(AMREX_D_DECL(115, 128, 0));
		    
                    for (auto &&idt : icpt_data)
                    {
			            //if(phi(idt.cellid_) < -0.9)
                        //    PORDER_ = 1;//drop to first order if the the interface is in damaged region
                        Sing_ = false; // Need this. Sing_ is only changed to true.
                        int ni = 0;
		        	
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
                                //amrex::PrintToFile("log") << "Error in identifying interface type.\n";
                                exit(1);
                            }
                            //U_[ni] = idt.u[ni];
                            //V_[ni] = idt.v[ni];
                            Compute_LSQ_Weights(w_[ni], std::fabs(idt.frac_[ni]));
                        }
			
                        int Num_i_ = ni;

		        /*	
                        int dist_x1 = abs(idt.cellid_[0] - domain.smallEnd(0));
                        int dist_x2 = abs(idt.cellid_[0] - domain.bigEnd(0));
                        int dist_y1 = abs(idt.cellid_[1] - domain.smallEnd(1));
                        int dist_y2 = abs(idt.cellid_[1] - domain.bigEnd(1));

                        //decide size of the lease-square stencil based on available points
                        if( dist_x1 >= 4 && dist_x2 >= 4 && dist_y1 >= 4 && dist_y2 >= 4)
                        {
                            //idt.stencil_ = 4;
                        }
                        else if( dist_x1 >= 3 && dist_x2 >= 3 && dist_y1 >= 3 && dist_y2 >= 3)
                        {
                            //idt.stencil_ = 3;
                            PORDER_ = 3;
                            //amrex::Print()<<"setting stencil size to 3 "<<idt.cellid_<<'\n';
                        }
                        else if( dist_x1 >= 2 && dist_x2 >= 2 && dist_y1 >= 2 && dist_y2 >= 2)
                        {
                            //idt.stencil_ = 2;
                            PORDER_ = 3;
                            //amrex::Print()<<"setting stencil size to 2 "<<idt.cellid_<<'\n';
                        }
                        else
                        {
                            //idt.stencil_ = 2;
                            PORDER_ = 3;
                            //amrex::Print()<<"setting stencil size to 1 "<<idt.cellid_<<'\n';
                        }
			*/

		    
                        /// create grown box around cell
                        amrex::Box gbx(idt.cellid_, idt.cellid_);
                        gbx.grow(stencil_);
		    
                        amrex::Box gbx_isect = gbx & domain;
		    
                        for (amrex::BoxIterator bit(gbx_isect); bit.ok(); ++bit)
                        {
                            const amrex::IntVect &iv = bit();
                            if (pmask(iv) == 1)
                            {
				                if(phi(iv[0],iv[1],0) >= -0.5 )
				                    PORDER_ = Fij_order;//Use higher order if any of the cloud points are
                                                        //in near the gel-water interface	
                                if(isAxisymmetric)
				                {
                                    //FillFij::Axisymmetric::Fij_LSQ_Parameters(F11_[ni], F12_[ni], F21_[ni], F22_[ni], F33_[ni], x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], F11, F12, F21, F22, F33, prob_lo, dx,domain);
			                        {
                                        iscalar = 2;
                                        LSQ_Parameters_Fij(F11_[ni], iscalar, x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], F11, prob_lo, dx, domain);
                                    }
                                    {
                                        iscalar = 3;
                                        LSQ_Parameters_Fij(F12_[ni], iscalar, x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], F12, prob_lo, dx, domain);
                                    }
                                    {
                                        iscalar = 4;
                                        LSQ_Parameters_Fij(F21_[ni], iscalar, x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], F21, prob_lo, dx, domain);
                                    }
                                    {
                                        iscalar = 5;
                                        LSQ_Parameters_Fij(F22_[ni], iscalar, x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], F22, prob_lo, dx, domain);
                                    }
                                    {
                                        iscalar = 6;
                                        LSQ_Parameters_Fij(F33_[ni], iscalar, x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], F33, prob_lo, dx, domain);
                                    }
				                }
                                else
                                    Fij_LSQ_Parameters(F11_[ni], F12_[ni], F21_[ni], F22_[ni], F33_[ni], x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], F11, F12, F21, F22, F33, prob_lo, dx);
                                ni++;
                            }
                            /*if( idt.cellid_ == test_iv)
                            {
                                amrex::Print()<<" ni = "<<ni<<'\n';
                                amrex::Print()<<" iv = "<<iv[0]<<" , "<<iv[1]<<'\n';
                                amrex::Print()<<"in domain = "<<domain.contains(iv)<<'\n';
                            }*/
                        }
		    
                        int Num_ = ni;
		    
                        if (solid->isAdvectLS())
                        {
                            i_beg = Num_i_;
                        }

			if (PORDER_ == 4)
                        {
                            if (Num_ < 3 + Num_i_)
                            {
                                // too few points have been found...can only do a simple average
                                sum_w = sum_sol_F11 = sum_sol_F12 = sum_sol_F21 = sum_sol_F22 = sum_sol_F33 = 0.0;
                                for (int i = i_beg; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol_F11 += w_[i] * F11_[i];
				    sum_sol_F12 += w_[i] * F12_[i];
				    sum_sol_F21 += w_[i] * F21_[i];
				    sum_sol_F22 += w_[i] * F22_[i];
				    sum_sol_F33 += w_[i] * F33_[i];
                                }
                                F11(idt.cellid_) = sum_sol_F11 / sum_w;
                                F12(idt.cellid_) = sum_sol_F12 / sum_w;
				F21(idt.cellid_) = sum_sol_F21 / sum_w;
                                F22(idt.cellid_) = sum_sol_F22 / sum_w;
				F33(idt.cellid_) = sum_sol_F33 / sum_w;
		    
                                //amrex::PrintToFile("log") << " Inside first \n";
                            }
			    else if(Num_ < 6 + Num_i_)
		            {
                                if (isAxisymmetric)
			        {
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
				}
			        else
			        {
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
				
				}
                                if (Sing_)
                                {
                                    // too few points have been found...can only do a simple average
                                    sum_w = sum_sol_F11 = sum_sol_F12 = sum_sol_F21 = sum_sol_F22 = sum_sol_F33 = 0.0;
                                    for (int i = i_beg; i < Num_; i++)
                                    {
                                        sum_w += w_[i];
                                        sum_sol_F11 += w_[i] * F11_[i];
                                        sum_sol_F12 += w_[i] * F12_[i];
                                        sum_sol_F21 += w_[i] * F21_[i];
                                        sum_sol_F22 += w_[i] * F22_[i];
                                        sum_sol_F33 += w_[i] * F33_[i];
                                    }
                                    F11(idt.cellid_) = sum_sol_F11 / sum_w;
                                    F12(idt.cellid_) = sum_sol_F12 / sum_w;
                                    F21(idt.cellid_) = sum_sol_F21 / sum_w;
                                    F22(idt.cellid_) = sum_sol_F22 / sum_w;
                                    F33(idt.cellid_) = sum_sol_F33 / sum_w;
                                }
			    }
                            else if(Num_ < 10 + Num_i_)
                            {
				if (isAxisymmetric)
				{
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                    if(Sing_)
                                    {
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                    }
				}
				else
				{
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);

                                    //FillFij::QR_Third_Order_DivFr(Sing_, Num_, Num_i_, x_, y_, F11_, F21_, w_, F11, F21, idt, 2,4);
				    //FillFij::QR_Third_Order_DivFr(Sing_, Num_, Num_i_, x_, y_, F12_, F22_, w_, F12, F22, idt, 3,5);
				    if(Sing_)
				    {
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);

				    }
				}
                            }
                            else
                            {
				//amrex::Print()<<"4th order gfm"<<'\n';
                                if (isAxisymmetric)
                                {
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                    //
                                    FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                    FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                    if(Sing_)
                                    {
					Sing_ = false;
                                        FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                        FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                    }
                                    if(Sing_)
                                    {
					Sing_ = false;
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                    }
                                }
                                else
                                {
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    //FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                    //FillFij::QR_Third_Order_DivFr(Sing_, Num_, Num_i_, x_, y_, F11_, F21_, w_, F11, F21, idt, 2,4);
                                    //FillFij::QR_Third_Order_DivFr(Sing_, Num_, Num_i_, x_, y_, F12_, F22_, w_, F12, F22, idt, 3,5);

                                    FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                    //FillFij::QR_Fourth_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                    if(Sing_)
                                    {
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);

                                    }
                                }
                            }
                        }
                        else if (PORDER_ == 3)
                        {
                            if (Num_ < 3 + Num_i_)
                            {
                                // too few points have been found...can only do a simple average
                                sum_w = sum_sol_F11 = sum_sol_F12 = sum_sol_F21 = sum_sol_F22 = sum_sol_F33 = 0.0;
                                for (int i = i_beg; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol_F11 += w_[i] * F11_[i];
				    sum_sol_F12 += w_[i] * F12_[i];
				    sum_sol_F21 += w_[i] * F21_[i];
				    sum_sol_F22 += w_[i] * F22_[i];
				    sum_sol_F33 += w_[i] * F33_[i];
                                }
                                F11(idt.cellid_) = sum_sol_F11 / sum_w;
                                F12(idt.cellid_) = sum_sol_F12 / sum_w;
				F21(idt.cellid_) = sum_sol_F21 / sum_w;
                                F22(idt.cellid_) = sum_sol_F22 / sum_w;
				F33(idt.cellid_) = sum_sol_F33 / sum_w;
		    
                                //amrex::PrintToFile("log") << " Inside first \n";
                            }
			    else if(Num_ < 6 + Num_i_)
		            {
                                if (isAxisymmetric)
			        {
                                     FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                     FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                     FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                     FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                     FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                     if(std::isnan(F11(idt.cellid_))||
                                        std::isnan(F12(idt.cellid_))||
                                        std::isnan(F21(idt.cellid_))||
                                        std::isnan(F22(idt.cellid_))||
                                        std::isnan(F33(idt.cellid_)))
                                         Sing_ = true;
                                     //else
                                     //    Sing_ = false;
				}
			        else
			        {
                                     FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                     FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                     FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                     FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                     if(std::isnan(F11(idt.cellid_))||
                                        std::isnan(F12(idt.cellid_))||
                                        std::isnan(F21(idt.cellid_))||
                                        std::isnan(F22(idt.cellid_)))
                                         Sing_ = true;
                                     else
                                         Sing_ = false;

				
				}

                                if (Sing_)
                                {
                                    // too few points have been found...can only do a simple average
                                    sum_w = sum_sol_F11 = sum_sol_F12 = sum_sol_F21 = sum_sol_F22 = sum_sol_F33 = 0.0;
                                    for (int i = i_beg; i < Num_; i++)
                                    {
                                        sum_w += w_[i];
                                        sum_sol_F11 += w_[i] * F11_[i];
                                        sum_sol_F12 += w_[i] * F12_[i];
                                        sum_sol_F21 += w_[i] * F21_[i];
                                        sum_sol_F22 += w_[i] * F22_[i];
                                        sum_sol_F33 += w_[i] * F33_[i];
                                    }
                                    F11(idt.cellid_) = sum_sol_F11 / sum_w;
                                    F12(idt.cellid_) = sum_sol_F12 / sum_w;
                                    F21(idt.cellid_) = sum_sol_F21 / sum_w;
                                    F22(idt.cellid_) = sum_sol_F22 / sum_w;
                                    F33(idt.cellid_) = sum_sol_F33 / sum_w;
                                }
			    }
                            else
                            {
				if (isAxisymmetric)
				{
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                    if(std::isnan(F11(idt.cellid_))||
                                       std::isnan(F12(idt.cellid_))||
                                       std::isnan(F21(idt.cellid_))||
                                       std::isnan(F22(idt.cellid_))||
                                       std::isnan(F33(idt.cellid_)))
                                       Sing_ = true;
                                    if(Sing_)
                                    {
					amrex::Print(-1)<<"Dorpping F_ij extension order to 2nd"<<'\n';
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                        if(std::isnan(F11(idt.cellid_))||
                                           std::isnan(F12(idt.cellid_))||
                                           std::isnan(F21(idt.cellid_))||
                                           std::isnan(F22(idt.cellid_))||
                                           std::isnan(F33(idt.cellid_)))
                                           Sing_ = true;
                                        
                                        if(Sing_)
                                        {
					    amrex::Print(-1)<<"Dorpping F_ij extension order to 1st"<<'\n';
                                            sum_w = sum_sol_F11 = sum_sol_F12 = sum_sol_F21 = sum_sol_F22 = sum_sol_F33 = 0.0;
                                            for (int i = i_beg; i < Num_; i++)
                                            {
                                                sum_w += w_[i];
                                                sum_sol_F11 += w_[i] * F11_[i];
                                                sum_sol_F12 += w_[i] * F12_[i];
                                                sum_sol_F21 += w_[i] * F21_[i];
                                                sum_sol_F22 += w_[i] * F22_[i];
                                                sum_sol_F33 += w_[i] * F33_[i];
                                            }
                                            F11(idt.cellid_) = sum_sol_F11 / sum_w;
                                            F12(idt.cellid_) = sum_sol_F12 / sum_w;
                                            F21(idt.cellid_) = sum_sol_F21 / sum_w;
                                            F22(idt.cellid_) = sum_sol_F22 / sum_w;
                                            F33(idt.cellid_) = sum_sol_F33 / sum_w;
                                        }

                                    }
				}
				else
				{
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                    FillFij::QR_Third_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);

				    if(Sing_)
				    {
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
				    }
				}
                            }
                        }
			else if(PORDER_ == 2)
			{
		            //amrex::Print()<<"Second order gfm"<<'\n';
                            if (Num_ < 3 + Num_i_)
                            {
                                // too few points have been found...can only do a simple average
                                sum_w = sum_sol_F11 = sum_sol_F12 = sum_sol_F21 = sum_sol_F22 = sum_sol_F33 = 0.0;
                                for (int i = i_beg; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol_F11 += w_[i] * F11_[i];
                                    sum_sol_F12 += w_[i] * F12_[i];
                                    sum_sol_F21 += w_[i] * F21_[i];
                                    sum_sol_F22 += w_[i] * F22_[i];
                                    sum_sol_F33 += w_[i] * F33_[i];
                                }
                                F11(idt.cellid_) = sum_sol_F11 / sum_w;
                                F12(idt.cellid_) = sum_sol_F12 / sum_w;
                                F21(idt.cellid_) = sum_sol_F21 / sum_w;
                                F22(idt.cellid_) = sum_sol_F22 / sum_w;
                                F33(idt.cellid_) = sum_sol_F33 / sum_w;

                                //amrex::PrintToFile("log") << " Inside first \n";
                            }
                            else
                            {
                                if (isAxisymmetric)
                                {
					//FillFij::Axisymmetric::QR_Second_Order_DivF_z(Sing_, Num_, Num_i_, idt.cellid_, x_, y_, F11_, F21_, w_, F11, F21, idt, dx, prob_lo);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F33_, w_, F33, idt, 6);
                                }
                                else
                                {
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F11_, w_, F11, idt, 2);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F12_, w_, F12, idt, 3);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F21_, w_, F21, idt, 4);
                                        FillFij::QR_Second_Order(Sing_, Num_, Num_i_, x_, y_, F22_, w_, F22, idt, 5);

                                }
                            }
			}
			else if(PORDER_ == 1)
			{
			    //amrex::Print()<<"First order GFM"<<'\n';
                            sum_w = sum_sol_F11 = sum_sol_F12 = sum_sol_F21 = sum_sol_F22 = sum_sol_F33 = 0.0;
                            for (int i = i_beg; i < Num_; i++)
                            {
                                sum_w += w_[i];
                                sum_sol_F11 += w_[i] * F11_[i];
                                sum_sol_F12 += w_[i] * F12_[i];
                                sum_sol_F21 += w_[i] * F21_[i];
                                sum_sol_F22 += w_[i] * F22_[i];
                                sum_sol_F33 += w_[i] * F33_[i];
                            }
                            F11(idt.cellid_) = sum_sol_F11 / sum_w;
                            F12(idt.cellid_) = sum_sol_F12 / sum_w;
                            F21(idt.cellid_) = sum_sol_F21 / sum_w;
                            F22(idt.cellid_) = sum_sol_F22 / sum_w;
                            F33(idt.cellid_) = sum_sol_F33 / sum_w;
			}

                        if(std::isnan(F11(idt.cellid_))||
                           std::isnan(F12(idt.cellid_))||
                           std::isnan(F21(idt.cellid_))||
                           std::isnan(F22(idt.cellid_))||
                           std::isnan(F33(idt.cellid_)))
                        {
                            amrex::Print(-1)<<"Found nan in extension of F_ij"<<'\n';
                            amrex::Print()<<"cell : "<<idt.cellid_<<'\n';
                            amrex::Print()<<"Num_ = "<<Num_<<" , Num_i_ = "<<Num_i_<<'\n';
                            amrex::Print()<<"F11(idt.cellid_) : "<< F11(idt.cellid_)<<'\n';
                            amrex::Print()<<"F12(idt.cellid_) : "<< F12(idt.cellid_)<<'\n';
                            amrex::Print()<<"F21(idt.cellid_) : "<< F21(idt.cellid_)<<'\n';
                            amrex::Print()<<"F22(idt.cellid_) : "<< F22(idt.cellid_)<<'\n';
                            amrex::Print()<<"F33(idt.cellid_) : "<< F33(idt.cellid_)<<'\n';
                            for(int ii = 0; ii<Num_;ii++)
                            {
                                amrex::Print(-1)<<"ii = "<<ii<<'\n';
                                amrex::Print(-1)<<"x_[ii] = "<<x_[ii]<<" , y_[ii] = "<<y_[ii]<<'\n';
                                amrex::Print(-1)<<"F11_[ii] = "<<F11_[ii]<<'\n';
                            }
			    exit(9);
			}
                        /* 
                        amrex::Real x_gc = prob_lo[0] + dx[0] * (idt.cellid_[0] + 0.5);
                        amrex::Real y_gc = prob_lo[1] + dx[1] * (idt.cellid_[1] + 0.5);
                        amrex::Real r_gc = std::hypot(x_gc - solid->Xcp(), y_gc - solid->Ycp());
                        amrex::Real theta = std::acos((x_gc - solid->Xcp())/r_gc);
			if(y_gc < 0)
		            theta = 2*M_PI - theta; 
			amrex::Real R_by_r,R;
			if(isAxisymmetric)
			{
		            //amrex::Real R = std::pow(solid->Volume_0()/(4.0/3.0*M_PI),(1.0/3.0));
                            //R_by_r = std::pow(solid->Volume_0()/(4.0/3.0*M_PI),(1.0/3.0))/r_gc;///std::pow((solid->Volume_0() / solid->Volume()), 1.0/3.0); 
			    amrex::Real R = std::pow(r_gc*r_gc*r_gc + solid->Time(),1.0/3.0);
			    amrex::Real R_by_r = R/r_gc;
                            F11(idt.cellid_) = R_by_r*R_by_r*std::pow(std::cos(theta),2.0) + 1.0/R_by_r * std::pow(std::sin(theta),2.0);
                            F12(idt.cellid_) = (R_by_r*R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);
                            F21(idt.cellid_) = (R_by_r*R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);
                            F22(idt.cellid_) = R_by_r*R_by_r*std::pow(std::sin(theta),2.0) + 1.0/R_by_r * std::pow(std::cos(theta),2.0);
                            F33(idt.cellid_) = 1.0/R_by_r;
			}
			else
			{
			    amrex::Real R = std::pow(r_gc*r_gc + solid->Time(),1.0/2.0);
                            amrex::Real R_by_r = R/r_gc;

                            F11(idt.cellid_) = R_by_r*std::pow(std::cos(theta),2.0) + 1.0/R_by_r * std::pow(std::sin(theta),2.0);
                            F12(idt.cellid_) = (R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);;
                            F21(idt.cellid_) = (R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);;
                            F22(idt.cellid_) = R_by_r*std::pow(std::sin(theta),2.0) + 1.0/R_by_r * std::pow(std::cos(theta),2.0);
			}
                        for (int i_ = 0; i_ < idt.n_intercepts; i_++)
                        {
                            idt.F11_Int[i_] = F11(idt.cellid_);
                            idt.F12_Int[i_] = F12(idt.cellid_);
                            idt.F21_Int[i_] = F21(idt.cellid_);
                            idt.F22_Int[i_] = F22(idt.cellid_);
                            idt.F33_Int[i_] = F33(idt.cellid_);
                        }
                        */		
                        //Vel(idt.cellid_, 0) = 0.0;//sum_sol_u / sum_w;
                        //Vel(idt.cellid_, 1) = 0.0;//sum_sol_v / sum_w;
			/*
                        if(idt.cellid_ == test_iv)
                        {
                                amrex::Print()<<"cell : "<<idt.cellid_<<'\n';
                                amrex::Print()<<"Vel(idt.cellid_) : "<< Vel(idt.cellid_, 0)<<" , "<<Vel(idt.cellid_, 1)<<'\n';
                                amrex::Print()<<"Num_ = "<<Num_<<" , Num_i_ = "<<Num_i_<<'\n';
                                //for (ni = 0; ni < idt.n_intercepts; ni++)
                                //{
				    
                                    //amrex::Print()<<"idt.norm_shear_["<< ni <<"] = "<<idt.norm_shear_[ni]<<'\n';
                                //}
                                for(int ii = 0;ii<= Num_;ii++)
                                {
                                        amrex::Print()<<"ii = "<<ii<<"\t";
                                        amrex::Print()<<"x_ = "<<x_[ii]<<"\t";
                                        amrex::Print()<<"y_ = "<<y_[ii]<<"\n";
                                        amrex::Print()<<"U_(x_, 0) = "<<U_[ii]<<"\n";
				        amrex::Print()<<"V_(y_, 1) = "<<V_[ii]<<"\n";
                                }
                                //std::exit(9);
                        }
		        */
                    }    
                }
            }
	}
    }


    void Mask::ExtendFijLSQ(amrex::MultiFab &mfF11, amrex::MultiFab &mfF12, amrex::MultiFab &mfF21, amrex::MultiFab &mfF22, amrex::MultiFab &mfF33)
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
        double x_[N_Max], y_[N_Max], w_[N_Max], F11_[N_Max], F12_[N_Max], F21_[N_Max], F22_[N_Max];
	double F33_[N_Max];
    
        for (int outer = 0; outer < LAYERS; outer++)
        {
            mfF11.FillBoundary();
	    mfF21.FillBoundary();
	    mfF12.FillBoundary();
	    mfF22.FillBoundary();
	    mfF33.FillBoundary();
            for (amrex::MFIter mfi(PMask); mfi.isValid(); ++mfi)
            {
                amrex::Array4<int> const &pmask = PMask.array(mfi);
                amrex::Array4<amrex::Real> const &F11 = mfF11.array(mfi);
                amrex::Array4<amrex::Real> const &F12 = mfF12.array(mfi);
                amrex::Array4<amrex::Real> const &F21 = mfF21.array(mfi);
                amrex::Array4<amrex::Real> const &F22 = mfF22.array(mfi);
                amrex::Array4<amrex::Real> const &F33 = mfF33.array(mfi);

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
				    if(isAxisymmetric)
                                    {
                                        {
                                            iscalar = 2;
                                            LSQ_Parameters_Fij(F11_[i_], iscalar, x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], F11, prob_lo, dx, domain);
                                        }
                                        {
                                            iscalar = 3;
                                            LSQ_Parameters_Fij(F12_[i_], iscalar, x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], F12, prob_lo, dx, domain);
                                        }
                                        {
                                            iscalar = 4;
                                            LSQ_Parameters_Fij(F21_[i_], iscalar, x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], F21, prob_lo, dx, domain);
                                        }
                                        {
                                            iscalar = 5;
                                            LSQ_Parameters_Fij(F22_[i_], iscalar, x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], F22, prob_lo, dx, domain);
                                        }
                                        {
                                            iscalar = 6;
                                            LSQ_Parameters_Fij(F33_[i_], iscalar, x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], F33, prob_lo, dx, domain);
                                        }
                                    }
                                    else
                                        Fij_LSQ_Parameters(F11_[i_], F12_[i_], F21_[i_], F22_[i_], F33_[i_], x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], F11, F12, F21, F22, F33, prob_lo, dx);

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
                                    if(isAxisymmetric)
                                    {
                                        {
                                            iscalar = 2;
                                            LSQ_Parameters_Fij(F11_[i_], iscalar, x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], F11, prob_lo, dx, domain);
                                        }
                                        {
                                            iscalar = 3;
                                            LSQ_Parameters_Fij(F12_[i_], iscalar, x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], F12, prob_lo, dx, domain);
                                        }
                                        {
                                            iscalar = 4;
                                            LSQ_Parameters_Fij(F21_[i_], iscalar, x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], F21, prob_lo, dx, domain);
                                        }
                                        {
                                            iscalar = 5;
                                            LSQ_Parameters_Fij(F22_[i_], iscalar, x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], F22, prob_lo, dx, domain);
                                        }
                                        {
                                            iscalar = 6;
                                            LSQ_Parameters_Fij(F33_[i_], iscalar, x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], F33, prob_lo, dx, domain);
                                        }
                                    }
                                    else
                                        Fij_LSQ_Parameters(F11_[i_], F12_[i_], F21_[i_], F22_[i_], F33_[i_], x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], F11, F12, F21, F22, F33, prob_lo, dx);

                                    i_++;
                                }
                            }
                        }
                    }
                    int Num_ = i_;

                    int PORDER_ = Fij_order; 
                    if (outer <= 3 && PORDER_ >= 3)
                    {
			/*
			if(Num_ > 1000)
                        {
                            if(isAxisymmetric)
                            {
		                FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F11_, w_, F11, cell, 2);
				FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F12_, w_, F12, cell, 3);
				FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F21_, w_, F21, cell, 4);
				FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F22_, w_, F22, cell, 5);
				FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F33_, w_, F33, cell, 6);
                            }
			    else
                            {
                                FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F11_, w_, F11, cell, 2);                        
                                FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F12_, w_, F12, cell, 3);
                                FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F21_, w_, F21, cell, 4);
                                FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F22_, w_, F22, cell, 5);
                            }
			}
			else if(Num_ <=1000 && Num_ > 6)
		        */
			if(Num_ > 6)
                        {
                            if(isAxisymmetric)
                            {
                                FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F11_, w_, F11, cell, 2);
                                FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F12_, w_, F12, cell, 3);
                                FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F21_, w_, F21, cell, 4);
                                FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F22_, w_, F22, cell, 5);                        
                                FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F33_, w_, F33, cell, 6);
                                if(std::isnan(F11(cell))||
                                   std::isnan(F12(cell))||
                                   std::isnan(F21(cell))||
                                   std::isnan(F22(cell))||
                                   std::isnan(F33(cell)))
                                   Sing_ = true;

                            }
                            else
                            {
                                FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F11_, w_, F11, cell, 2);
                                FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F12_, w_, F12, cell, 3);
                                FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F21_, w_, F21, cell, 4);
                                FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F22_, w_, F22, cell, 5);
                                if(std::isnan(F11(cell))||
                                   std::isnan(F12(cell))||
                                   std::isnan(F21(cell))||
                                   std::isnan(F22(cell)))
                                   Sing_ = true;

                            }
                            if(Sing_)
                            {
				amrex::Print()<<"Dropping order of extesion"<<'\n';
                                if(isAxisymmetric)
                                {
                                    FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F11_, w_, F11, cell, 2);
                                    FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F12_, w_, F12, cell, 3);
                                    FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F21_, w_, F21, cell, 4);
                                    FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F22_, w_, F22, cell, 5);
                                    FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F33_, w_, F33, cell, 6);
                                }
                                else
                                {
                                    FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F11_, w_, F11, cell, 2);
                                    FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F12_, w_, F12, cell, 3);
                                    FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F21_, w_, F21, cell, 4);
                                    FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F22_, w_, F22, cell, 5);
                                }
                            }
                        }
                        else if(Num_ <=6)
                        {       
                            if(isAxisymmetric)
                            {
                                FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F11_, w_, F11, cell, 2);
                                FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F12_, w_, F12, cell, 3);
                                FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F21_, w_, F21, cell, 4);
                                FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F22_, w_, F22, cell, 5); 
                                FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F33_, w_, F33, cell, 6);
                            }
                            else
                            {
                                FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F11_, w_, F11, cell, 2);
                                FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F12_, w_, F12, cell, 3);
                                FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F21_, w_, F21, cell, 4);
                                FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F22_, w_, F22, cell, 5);
                                //FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F11_, w_, F11, cell, 2);
                                //FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F12_, w_, F12, cell, 3);
                                //FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F21_, w_, F21, cell, 4);
                                //FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F22_, w_, F22, cell, 5);

                            }
                        }
		    }
		    else
	            {
			//amrex::Print()<<"Second order extension"<<'\n';
	                if(isAxisymmetric) 
		        {
                            FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F11_, w_, F11, cell, 2);
                            FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F12_, w_, F12, cell, 3);
                            FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F21_, w_, F21, cell, 4);
                            FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F22_, w_, F22, cell, 5);
                            FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F33_, w_, F33, cell, 6);
                            //FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F11_, w_, F11, cell, 2);
                            //FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F12_, w_, F12, cell, 3);
                            //FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F21_, w_, F21, cell, 4);
                            //FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F22_, w_, F22, cell, 5);
                            //FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F33_, w_, F33, cell, 6);
			    //FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F11_, w_, F11, cell, 2);
                            //FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F12_, w_, F12, cell, 3);
                            //FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F21_, w_, F21, cell, 4);
                            //FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F22_, w_, F22, cell, 5);
                            //FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F33_, w_, F33, cell, 6);
			}
                        else
		        {
                            FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F11_, w_, F11, cell, 2);
                            FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F12_, w_, F12, cell, 3);
                            FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F21_, w_, F21, cell, 4);
                            FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F22_, w_, F22, cell, 5);
                            //FillFij::LS_QR_Second_Order(Sing_, Num_, x_, y_, F33_, w_, F33, cell, 6);

                            //FillFij::LS_QR_Second_Order_F11F21_2D(Sing_, Num_, 0, x_, y_, F11_, F21_, w_, F11, F21, cell);
                            //FillFij::LS_QR_Second_Order_F11F21_2D(Sing_, Num_, 0, x_, y_, F12_, F22_, w_, F12, F22, cell);
			    //
                            //FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F11_, w_, F11, cell, 2);
                            //FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F12_, w_, F12, cell, 3);
                            //FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F21_, w_, F21, cell, 4);
                            //FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F22_, w_, F22, cell, 5);
                            //FillFij::LS_QR_Third_Order(Sing_, Num_, x_, y_, F33_, w_, F33, cell, 6);
                            //FillFij::LS_QR_Third_Order_DivFr(Sing_, Num_, x_, y_, F11_, F21_, w_, F11, F21, cell, 2, 4);
                            //FillFij::LS_QR_Third_Order_DivFr(Sing_, Num_, x_, y_, F12_, F22_, w_, F12, F22, cell, 3, 5);
			    //
                            //FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F11_, w_, F11, cell, 2);
                            //FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F12_, w_, F12, cell, 3);
                            //FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F21_, w_, F21, cell, 4);
                            //FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F22_, w_, F22, cell, 5);
                            //FillFij::LS_QR_Fourth_Order(Sing_, Num_, x_, y_, F33_, w_, F33, cell, 6);
			}
                    }
	            //Analytical	
                    /*	
                    auto &&solid = (*interfaces)[0];
                    {
                        amrex::Real x_gc = prob_lo[0] + dx[0] * (cell[0] + 0.5);
                        amrex::Real y_gc = prob_lo[1] + dx[1] * (cell[1] + 0.5);
                        amrex::Real r_gc = std::hypot(x_gc - solid->Xcp(), y_gc - solid->Ycp());
                        amrex::Real theta = std::acos((x_gc - solid->Xcp())/r_gc);
                        amrex::Real R, R_by_r;
                        if(y_gc<0)
                            theta = 2*M_PI - theta;
			if(isAxisymmetric)
                        {
                            //amrex::Real R_b_0 = std::pow(solid->Volume_0()/(4.0/3.0*M_PI),(1.0/3.0));
                            //amrex::Real R_b = std::pow(solid->Volume()/(4.0/3.0*M_PI),(1.0/3.0));

		            //R_by_r = std::pow(solid->Volume_0()/(4.0/3.0*M_PI),(1.0/3.0))/r_gc;///std::pow((solid->Volume_0() / solid->Volume()), 1.0/3.0);
                            //amrex::Real R = std::pow(r_gc*r_gc*r_gc + (std::pow(R_b_0,3.0) - std::pow(R_b,3.0)),1.0/3.0);
                            R = std::pow(r_gc*r_gc*r_gc + solid->Time(),1.0/3.0);
                            amrex::Real R_by_r = R/r_gc;

                            F11(cell) = R_by_r*R_by_r*std::pow(std::cos(theta),2.0) + 1.0/R_by_r * std::pow(std::sin(theta),2.0);
                            F12(cell) = (R_by_r*R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);
                            F21(cell) = (R_by_r*R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);
                            F22(cell) = R_by_r*R_by_r*std::pow(std::sin(theta),2.0) + 1.0/R_by_r * std::pow(std::cos(theta),2.0);
                            F33(cell) = 1.0/R_by_r;
			}
			else
		        {
                            R = std::pow(r_gc*r_gc + solid->Time(),1.0/2.0);
                            R_by_r = R/r_gc;
                            //amrex::Print()<<"r_gc = "<<r_gc<<" , theta = "<<theta<<'\n';
                            F11(cell) = R_by_r*std::pow(std::cos(theta),2.0) + 1.0/R_by_r * std::pow(std::sin(theta),2.0);
                            F12(cell) = (R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);
                            F21(cell) = (R_by_r - 1.0/R_by_r)*std::sin(theta)*std::cos(theta);
                            F22(cell) = R_by_r*std::pow(std::sin(theta),2.0) + 1.0/R_by_r * std::pow(std::cos(theta),2.0);
                            //F33(cell) = 1.0/R_by_r;
			}
                    }
		    */
                    /* 
                    if(cell[0] == 32867 && cell[1] == 1)
                    {
                        amrex::Print(-1)<<"i = "<<cell[0]<<" , j = "<<cell[1]<<'\n';
			amrex::Print(-1)<<"pmask = "<<pmask(cell)<<'\n';
                        amrex::Print(-1)<<"outer = "<<outer<<'\n';
                        amrex::Print(-1)<<"Num_ = "<<Num_<<'\n';
			amrex::Print(-1)<<"F12 = "<<F12(cell)<<'\n';
			
			for(int ii = 0; ii<Num_;ii++)
		        {
			    amrex::Print(-1)<<"ii = "<<ii<<'\n';
			    amrex::Print(-1)<<"x_[ii] = "<<x_[ii]<<" , y_[ii] = "<<y_[ii]<<'\n';
			    amrex::Print(-1)<<"F12_[ii] = "<<F12_[ii]<<'\n';
			}
                    }
		    */
		    
		    
		    
		    
                }
            }
            //CollocatedVelBoundaryConditions();
        }
    }
}
