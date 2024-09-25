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
        
        
        void Vel_LSQ_Parameters(
            double &U_,
            double &V_,
            double &x_,
            double &y_,
            double &w_,
            int i1, int j1,
            int i, int j,
            amrex::Array4<amrex::Real const> const &U, /// cc vel with 2 comp
            const amrex::Real *prob_lo,
            const amrex::Real *dx)
        {
            amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
            amrex::Real x1 = prob_lo[0] + dx[0] * (i1 + 0.5);
            amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
            amrex::Real y1 = prob_lo[1] + dx[1] * (j1 + 0.5);
        
            x_ = (x1 - x) / dx[0];
            y_ = (y1 - y) / dx[1];
            U_ = U(i1, j1, 0, 0);
            V_ = U(i1, j1, 0, 1);
        
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
        
        
        void Int_LS_Vel_Second_Order(
            bool &sing,
            int Num_,
            int Num_i_,
            InterceptData &idata,
            amrex::Array4<amrex::Real> const &U,
            double *x_,
            double *y_,
            double *U_,
            double *V_,
            double *w_)
        {
            int i_, j_, k_, N_;
            double **Mat;    // Least squares and constraint matrices
            double **Q, **R; // Q R decomosition of Mat;
            double **Q1, **Q2T, **Q2T_Q, **Q2T_R;
            double *Sol, *Vec_, *Constr_, *Sol_u_, *Sol_w_, *Q1Tb; // *Q2T_QTQ1Tb ;
            double dx, dy, du, dv, nx, ny;
        
            N_ = Num_ - Num_i_;
        
            Allocate_2D_R(Mat, 2 * Num_, 5);
            Sol = new double[5];
            Vec_ = new double[2 * N_];
            Allocate_2D_R(Q, 2 * Num_, 5);
            Allocate_2D_R(R, 5, 5);
            Allocate_2D_R(Q1, 2 * N_, 5);
            Allocate_2D_R(Q2T, 5, 2 * Num_i_);
            Allocate_2D_R(Q2T_Q, 5, 2 * Num_i_);
            Allocate_2D_R(Q2T_R, 2 * Num_i_, 2 * Num_i_);
            Constr_ = new double[2 * Num_i_];
            Sol_u_ = new double[2 * Num_i_];
            Sol_w_ = new double[2 * Num_i_];
            Q1Tb = new double[5]; // Q2T_QTQ1Tb = new double[Num_i_] ;
        
            for (i_ = Num_i_; i_ < Num_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                du = U_[i_];
                dv = V_[i_];
        
                j_ = i_ - Num_i_;
                Mat[j_][0] = w_[i_];      // w
                Mat[j_][1] = w_[i_] * dx; // x
                Mat[j_][2] = w_[i_] * dy; // y
                Mat[j_][3] = 0.0;
                Mat[j_][4] = 0.0;
        
                Mat[j_ + N_][0] = 0.0;          // w
                Mat[j_ + N_][1] = -w_[i_] * dy; // x
                Mat[j_ + N_][2] = 0.0;          // y
                Mat[j_ + N_][3] = w_[i_];
                Mat[j_ + N_][4] = w_[i_] * dx;
        
                Vec_[j_] = w_[i_] * du;      // x
                Vec_[j_ + N_] = w_[i_] * dv; // x
            }
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                du = U_[i_];
                dv = V_[i_];
                j_ = 2 * N_ + i_;
        
                Mat[j_][0] = 1.0; // w
                Mat[j_][1] = dx;  // x
                Mat[j_][2] = dy;  // y
                Mat[j_][3] = 0.0;
                Mat[j_][4] = 0.0;
        
                Mat[j_ + Num_i_][0] = 0.0; // w
                Mat[j_ + Num_i_][1] = -dy; // x
                Mat[j_ + Num_i_][2] = 0.0; // y
                Mat[j_ + Num_i_][3] = 1.0;
                Mat[j_ + Num_i_][4] = dx;
        
                Constr_[i_] = du;
                Constr_[i_ + Num_i_] = dv;
            }
            QRdcmp<double> QR(Mat, 2 * Num_, 5);
            QR.get_Q(Q);
            QR.get_R(R);
        
            for (j_ = 0; j_ < 5; j_++)
            {
                for (i_ = 0; i_ < 2 * N_; i_++)
                    Q1[i_][j_] = Q[i_][j_];
                for (i_ = 0; i_ < 2 * Num_i_; i_++)
                    Q2T[j_][i_] = Q[2 * N_ + i_][j_];
            }
            QRdcmp<double> QR1(Q2T, 5, 2 * Num_i_);
            QR1.get_Q(Q2T_Q);
            QR1.get_R(Q2T_R);
        
            // solve for (Q2T_R)T u = Constr ;
            for (i_ = 0; i_ < 2 * Num_i_; i_++)
            {
                Sol_u_[i_] = Constr_[i_] / Q2T_R[i_][i_];
                for (j_ = 0; j_ < i_; j_++)
                    Sol_u_[i_] -= Sol_u_[j_] * Q2T_R[j_][i_] / Q2T_R[i_][i_];
            }
        
            // Find the vector Q1T b.
            for (i_ = 0; i_ < 5; i_++)
            {
                Q1Tb[i_] = 0.0;
                for (j_ = 0; j_ < 2 * N_; j_++)
                    Q1Tb[i_] += Q1[j_][i_] * Vec_[j_];
            }
        
            // Find the vector 2.0*(Q2T_Q)T Q1Tb -2.0*Sol_u_;
            for (i_ = 0; i_ < 2 * Num_i_; i_++)
            {
                Vec_[i_] = -2.0 * Sol_u_[i_];
                for (j_ = 0; j_ < 5; j_++)
                    Vec_[i_] += 2.0 * Q2T_Q[j_][i_] * Q1Tb[j_];
            }
            // solve (Q2T_R) w = RHS from above
            for (i_ = 2 * Num_i_ - 1; i_ >= 0; i_--)
            {
                Sol_w_[i_] = Vec_[i_] / Q2T_R[i_][i_];
                for (j_ = i_ + 1; j_ < 2 * Num_i_; j_++)
                    Sol_w_[i_] -= Sol_w_[j_] * Q2T_R[i_][j_] / Q2T_R[i_][i_];
            }
        
            // Find the vector (Q1T_b) - Q2Tw/2;
            for (i_ = 0; i_ < 5; i_++)
            {
                Vec_[i_] = Q1Tb[i_];
                for (j_ = 0; j_ < 2 * Num_i_; j_++)
                    Vec_[i_] -= 0.5 * Q2T[i_][j_] * Sol_w_[j_];
            }
            // solve (Mat_R) w = RHS from above
            for (i_ = 4; i_ >= 0; i_--)
            {
                Sol[i_] = Vec_[i_] / R[i_][i_];
                for (j_ = i_ + 1; j_ < 5; j_++)
                    Sol[i_] -= Sol[j_] * R[i_][j_] / R[i_][i_];
            }
        
            U(idata.cellid_, 0) = Sol[0];
            U(idata.cellid_, 1) = Sol[3];
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                idata.U_comp[i_] = Sol[0] + Sol[1] * dx + Sol[2] * dy;
                idata.V_comp[i_] = Sol[3] + Sol[4] * dx - Sol[1] * dy;
            }
        
            for (i_ = 0; i_ < 2 * Num_; i_++)
            {
                delete[] Mat[i_];
                delete[] Q[i_];
            }
            for (i_ = 0; i_ < 5; i_++)
            {
                delete[] R[i_];
                delete[] Q2T[i_];
                delete[] Q2T_Q[i_];
            }
            for (i_ = 0; i_ < 2 * N_; i_++)
            {
                delete[] Q1[i_];
            }
            for (i_ = 0; i_ < 2 * Num_i_; i_++)
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
        }
        
        void Int_LS_Vel_Third_Order(
            bool &sing,
            int Num_,
            int Num_i_,
            InterceptData &idata,
            amrex::Array4<amrex::Real> const &U,
            double *x_,
            double *y_,
            double *U_,
            double *V_,
            double *w_)
        {
        
            int i_, j_, k_, N_;
            double **Mat;    // Least squares and constraint matrices
            double **Q, **R; // Q R decomosition of Mat;
            double **Q1, **Q2T, **Q2T_Q, **Q2T_R;
            double *Sol, *Vec_, *Constr_, *Sol_u_, *Sol_w_, *Q1Tb; // *Q2T_QTQ1Tb ;
            double dx, dy, du, dv, nx, ny;
        
            N_ = Num_ - Num_i_;
        
            Allocate_2D_R(Mat, 2 * Num_, 9);
            Sol = new double[9];
            Vec_ = new double[2 * N_];
            Allocate_2D_R(Q, 2 * Num_, 9);
            Allocate_2D_R(R, 9, 9);
            Allocate_2D_R(Q1, 2 * N_, 9);
            Allocate_2D_R(Q2T, 9, 2 * Num_i_);
            Allocate_2D_R(Q2T_Q, 9, 2 * Num_i_);
            Allocate_2D_R(Q2T_R, 2 * Num_i_, 2 * Num_i_);
            Constr_ = new double[2 * Num_i_];
            Sol_u_ = new double[2 * Num_i_];
            Sol_w_ = new double[2 * Num_i_];
            Q1Tb = new double[9]; // Q2T_QTQ1Tb = new double[Num_i_] ;
        
            for (i_ = Num_i_; i_ < Num_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                du = U_[i_];
                dv = V_[i_];
        
                j_ = i_ - Num_i_;
                Mat[j_][0] = w_[i_];           // w
                Mat[j_][1] = w_[i_] * dx;      // x
                Mat[j_][2] = w_[i_] * dy;      // y
                Mat[j_][3] = w_[i_] * dx * dx; // x^2
                Mat[j_][4] = w_[i_] * dx * dy; // xy
                Mat[j_][5] = w_[i_] * dy * dy; // y^2
                Mat[j_][6] = 0.0;              // v_0
                Mat[j_][7] = 0.0;              // v_1
                Mat[j_][8] = 0.0;              // v_11
        
                Mat[j_ + N_][0] = 0.0;                     // w
                Mat[j_ + N_][1] = -w_[i_] * dy;            // -v_2*dy
                Mat[j_ + N_][2] = 0.0;                     // y
                Mat[j_ + N_][3] = -2.0 * w_[i_] * dx * dy; // v_12 = -2.0*u_11
                Mat[j_ + N_][4] = -0.5 * w_[i_] * dy * dy;
                Mat[j_ + N_][5] = 0.0;
                Mat[j_ + N_][6] = w_[i_];           // v_0
                Mat[j_ + N_][7] = w_[i_] * dx;      // v_1
                Mat[j_ + N_][8] = w_[i_] * dx * dx; // v_11
        
                Vec_[j_] = w_[i_] * du;      // x
                Vec_[j_ + N_] = w_[i_] * dv; // x
            }
        
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                du = U_[i_];
                dv = V_[i_];
                j_ = 2 * N_ + i_;
        
                Mat[j_][0] = 1.0; // w
                Mat[j_][1] = dx;  // x
                Mat[j_][2] = dy;  // y
                Mat[j_][3] = dx * dx;
                Mat[j_][4] = dx * dy;
                Mat[j_][5] = dy * dy;
                Mat[j_][6] = 0.0;
                Mat[j_][7] = 0.0;
                Mat[j_][8] = 0.0;
        
                Mat[j_ + Num_i_][0] = 0.0; // w
                Mat[j_ + Num_i_][1] = -dy; // x
                Mat[j_ + Num_i_][2] = 0.0; // y
                Mat[j_ + Num_i_][3] = -2.0 * dx * dy;
                Mat[j_ + Num_i_][4] = -0.5 * dy * dy;
                Mat[j_ + Num_i_][5] = 0.0;
                Mat[j_ + Num_i_][6] = 1.0;
                Mat[j_ + Num_i_][7] = dx;
                Mat[j_ + Num_i_][8] = dx * dx;
        
                Constr_[i_] = du;
                Constr_[i_ + Num_i_] = dv;
            }
            QRdcmp<double> QR(Mat, 2 * Num_, 9);
            QR.get_Q(Q);
            QR.get_R(R);
        
            for (j_ = 0; j_ < 9; j_++)
            {
                for (i_ = 0; i_ < 2 * N_; i_++)
                    Q1[i_][j_] = Q[i_][j_];
                for (i_ = 0; i_ < 2 * Num_i_; i_++)
                    Q2T[j_][i_] = Q[2 * N_ + i_][j_];
            }
            QRdcmp<double> QR1(Q2T, 9, 2 * Num_i_);
            QR1.get_Q(Q2T_Q);
            QR1.get_R(Q2T_R);
        
            // solve for (Q2T_R)T u = Constr ;
            for (i_ = 0; i_ < 2 * Num_i_; i_++)
            {
                Sol_u_[i_] = Constr_[i_] / Q2T_R[i_][i_];
                for (j_ = 0; j_ < i_; j_++)
                    Sol_u_[i_] -= Sol_u_[j_] * Q2T_R[j_][i_] / Q2T_R[i_][i_];
            }
        
            // Find the vector Q1T b.
            for (i_ = 0; i_ < 9; i_++)
            {
                Q1Tb[i_] = 0.0;
                for (j_ = 0; j_ < 2 * N_; j_++)
                    Q1Tb[i_] += Q1[j_][i_] * Vec_[j_];
            }
        
            // Find the vector 2.0*(Q2T_Q)T Q1Tb -2.0*Sol_u_;
            for (i_ = 0; i_ < 2 * Num_i_; i_++)
            {
                Vec_[i_] = -2.0 * Sol_u_[i_];
                for (j_ = 0; j_ < 9; j_++)
                    Vec_[i_] += 2.0 * Q2T_Q[j_][i_] * Q1Tb[j_];
            }
            // solve (Q2T_R) w = RHS from above
            for (i_ = 2 * Num_i_ - 1; i_ >= 0; i_--)
            {
                Sol_w_[i_] = Vec_[i_] / Q2T_R[i_][i_];
                for (j_ = i_ + 1; j_ < 2 * Num_i_; j_++)
                    Sol_w_[i_] -= Sol_w_[j_] * Q2T_R[i_][j_] / Q2T_R[i_][i_];
            }
        
            // Find the vector (Q1T_b) - Q2Tw/2;
            for (i_ = 0; i_ < 9; i_++)
            {
                Vec_[i_] = Q1Tb[i_];
                for (j_ = 0; j_ < 2 * Num_i_; j_++)
                    Vec_[i_] -= 0.5 * Q2T[i_][j_] * Sol_w_[j_];
            }
            // solve (Mat_R) w = RHS from above
            for (i_ = 8; i_ >= 0; i_--)
            {
                Sol[i_] = Vec_[i_] / R[i_][i_];
                for (j_ = i_ + 1; j_ < 9; j_++)
                    Sol[i_] -= Sol[j_] * R[i_][j_] / R[i_][i_];
            }
        
            U(idata.cellid_, 0) = Sol[0];
            U(idata.cellid_, 1) = Sol[6];
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                idata.U_comp[i_] = Sol[0] + Sol[1] * dx + Sol[2] * dy + Sol[3] * dx * dx + Sol[4] * dx * dy + Sol[5] * dy * dy;
                idata.V_comp[i_] = Sol[6] + Sol[7] * dx - Sol[1] * dy + Sol[8] * dx * dx - 2.0 * Sol[3] * dx * dy - 0.5 * Sol[4] * dy * dy;
            }
        
            for (i_ = 0; i_ < 2 * Num_; i_++)
            {
                delete[] Mat[i_];
                delete[] Q[i_];
            }
            for (i_ = 0; i_ < 9; i_++)
            {
                delete[] R[i_];
                delete[] Q2T[i_];
                delete[] Q2T_Q[i_];
            }
            for (i_ = 0; i_ < 2 * N_; i_++)
            {
                delete[] Q1[i_];
            }
            for (i_ = 0; i_ < 2 * Num_i_; i_++)
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
        }
        
        void Int_LS_Vel_Fourth_Order(
            bool &sing,
            int Num_,
            int Num_i_,
            InterceptData &idata,
            amrex::Array4<amrex::Real> const &U,
            double *x_,
            double *y_,
            double *U_,
            double *V_,
            double *w_,
            const amrex::Real &MU_0,
            const amrex::Real *deltax)
        {
            int i_, j_, k_, N_;
            double **Mat;    // Least squares and constraint matrices
            double **Q, **R; // Q R decomosition of Mat;
            double **Q1, **Q2T, **Q2T_Q, **Q2T_R;
            double *Sol, *Vec_, *Constr_, *Sol_u_, *Sol_w_, *Q1Tb; // *Q2T_QTQ1Tb ;
            double dx, dy, du, dv, nx, ny;
        
            N_ = Num_ - Num_i_;
        
            Allocate_2D_R(Mat, 2 * Num_, 14);
            Sol = new double[14];
            Vec_ = new double[2 * N_];
            Allocate_2D_R(Q, 2 * Num_, 14);
            Allocate_2D_R(R, 14, 14);
            Allocate_2D_R(Q1, 2 * N_, 14);
            Allocate_2D_R(Q2T, 14, 2 * Num_i_);
            Allocate_2D_R(Q2T_Q, 14, 2 * Num_i_);
            Allocate_2D_R(Q2T_R, 2 * Num_i_, 2 * Num_i_);
            Constr_ = new double[2 * Num_i_];
            Sol_u_ = new double[2 * Num_i_];
            Sol_w_ = new double[2 * Num_i_];
            Q1Tb = new double[14]; // Q2T_QTQ1Tb = new double[Num_i_] ;
        
            for (i_ = Num_i_; i_ < Num_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                du = U_[i_];
                dv = V_[i_];
        
                j_ = i_ - Num_i_;
        
                Mat[j_][0] = w_[i_];           // w
                Mat[j_][1] = w_[i_] * dx;      // x
                Mat[j_][2] = w_[i_] * dy;      // y
                Mat[j_][3] = w_[i_] * dx * dx; // x^2
                Mat[j_][4] = w_[i_] * dx * dy; // xy
                Mat[j_][5] = w_[i_] * dy * dy; // y^2
        
                Mat[j_][6] = w_[i_] * dx * dx * dx; // v_0
                Mat[j_][7] = w_[i_] * dx * dx * dy; // v_1
                Mat[j_][8] = w_[i_] * dx * dy * dy; // v_11
                Mat[j_][9] = w_[i_] * dy * dy * dy; // v_11
        
                Mat[j_][10] = 0.0; // v_11
                Mat[j_][11] = 0.0; // v_11
                Mat[j_][12] = 0.0; // v_11
                Mat[j_][13] = 0.0; // v_11
        
                Mat[j_ + N_][0] = 0.0;                     // w
                Mat[j_ + N_][1] = -w_[i_] * dy;            // -v_2*dy
                Mat[j_ + N_][2] = 0.0;                     // y
                Mat[j_ + N_][3] = -2.0 * w_[i_] * dx * dy; // v_12 = -2.0*u_11
                Mat[j_ + N_][4] = -0.5 * w_[i_] * dy * dy;
                Mat[j_ + N_][5] = 0.0;
        
                Mat[j_ + N_][6] = -3.0 * w_[i_] * dx * dx * dy;
                Mat[j_ + N_][7] = -w_[i_] * dx * dy * dy;
                Mat[j_ + N_][8] = -w_[i_] * dy * dy * dy / 3.0;
                Mat[j_ + N_][9] = 0.0;
        
                Mat[j_ + N_][10] = w_[i_];           // v_0
                Mat[j_ + N_][11] = w_[i_] * dx;      // v_1
                Mat[j_ + N_][12] = w_[i_] * dx * dx; // v_11
                Mat[j_ + N_][13] = w_[i_] * dx * dx * dx;
        
                Vec_[j_] = w_[i_] * du;      // x
                Vec_[j_ + N_] = w_[i_] * dv; // x
            }
        
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                du = U_[i_];
                dv = V_[i_];
                j_ = 2 * N_ + i_;
        
                Mat[j_][0] = 1.0;     // w
                Mat[j_][1] = dx;      // x
                Mat[j_][2] = dy;      // y
                Mat[j_][3] = dx * dx; // x^2
                Mat[j_][4] = dx * dy; // xy
                Mat[j_][5] = dy * dy; // y^2
        
                Mat[j_][6] = dx * dx * dx; // v_0
                Mat[j_][7] = dx * dx * dy; // v_1
                Mat[j_][8] = dx * dy * dy; // v_11
                Mat[j_][9] = dy * dy * dy; // v_11
        
                Mat[j_][10] = 0.0; // v_11
                Mat[j_][11] = 0.0; // v_11
                Mat[j_][12] = 0.0; // v_11
                Mat[j_][13] = 0.0; // v_11
        
                Mat[j_ + Num_i_][0] = 0.0;            // w
                Mat[j_ + Num_i_][1] = -dy;            // -v_2*dy
                Mat[j_ + Num_i_][2] = 0.0;            // y
                Mat[j_ + Num_i_][3] = -2.0 * dx * dy; // v_12 = -2.0*u_11
                Mat[j_ + Num_i_][4] = -0.5 * dy * dy;
                Mat[j_ + Num_i_][5] = 0.0;
        
                Mat[j_ + Num_i_][6] = -3.0 * dx * dx * dy;
                Mat[j_ + Num_i_][7] = -dx * dy * dy;
                Mat[j_ + Num_i_][8] = -dy * dy * dy / 3.0;
                Mat[j_ + Num_i_][9] = 0.0;
        
                Mat[j_ + Num_i_][10] = 1.0;     // v_0
                Mat[j_ + Num_i_][11] = dx;      // v_1
                Mat[j_ + Num_i_][12] = dx * dx; // v_11
                Mat[j_ + Num_i_][13] = dx * dx * dx;
        
                Constr_[i_] = du;
                Constr_[i_ + Num_i_] = dv;
            }
            QRdcmp<double> QR(Mat, 2 * Num_, 14);
            QR.get_Q(Q);
            QR.get_R(R);
        
            for (j_ = 0; j_ < 14; j_++)
            {
                for (i_ = 0; i_ < 2 * N_; i_++)
                    Q1[i_][j_] = Q[i_][j_];
                for (i_ = 0; i_ < 2 * Num_i_; i_++)
                    Q2T[j_][i_] = Q[2 * N_ + i_][j_];
            }
            QRdcmp<double> QR1(Q2T, 14, 2 * Num_i_);
            QR1.get_Q(Q2T_Q);
            QR1.get_R(Q2T_R);
        
            // solve for (Q2T_R)T u = Constr ;
            for (i_ = 0; i_ < 2 * Num_i_; i_++)
            {
                Sol_u_[i_] = Constr_[i_] / Q2T_R[i_][i_];
                for (j_ = 0; j_ < i_; j_++)
                    Sol_u_[i_] -= Sol_u_[j_] * Q2T_R[j_][i_] / Q2T_R[i_][i_];
            }
        
            // Find the vector Q1T b.
            for (i_ = 0; i_ < 14; i_++)
            {
                Q1Tb[i_] = 0.0;
                for (j_ = 0; j_ < 2 * N_; j_++)
                    Q1Tb[i_] += Q1[j_][i_] * Vec_[j_];
            }
        
            // Find the vector 2.0*(Q2T_Q)T Q1Tb -2.0*Sol_u_;
            for (i_ = 0; i_ < 2 * Num_i_; i_++)
            {
                Vec_[i_] = -2.0 * Sol_u_[i_];
                for (j_ = 0; j_ < 14; j_++)
                    Vec_[i_] += 2.0 * Q2T_Q[j_][i_] * Q1Tb[j_];
            }
            // solve (Q2T_R) w = RHS from above
            for (i_ = 2 * Num_i_ - 1; i_ >= 0; i_--)
            {
                Sol_w_[i_] = Vec_[i_] / Q2T_R[i_][i_];
                for (j_ = i_ + 1; j_ < 2 * Num_i_; j_++)
                    Sol_w_[i_] -= Sol_w_[j_] * Q2T_R[i_][j_] / Q2T_R[i_][i_];
            }
        
            // Find the vector (Q1T_b) - Q2Tw/2;
            for (i_ = 0; i_ < 14; i_++)
            {
                Vec_[i_] = Q1Tb[i_];
                for (j_ = 0; j_ < 2 * Num_i_; j_++)
                    Vec_[i_] -= 0.5 * Q2T[i_][j_] * Sol_w_[j_];
            }
            // solve (Mat_R) w = RHS from above
            for (i_ = 13; i_ >= 0; i_--)
            {
                Sol[i_] = Vec_[i_] / R[i_][i_];
                for (j_ = i_ + 1; j_ < 14; j_++)
                    Sol[i_] -= Sol[j_] * R[i_][j_] / R[i_][i_];
            }
        
            U(idata.cellid_, 0) = Sol[0];
            U(idata.cellid_, 1) = Sol[10];
            double ux, uy, vx, vy;
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
        
                idata.U_comp[i_] = Sol[0] + Sol[1] * dx + Sol[2] * dy + Sol[3] * dx * dx + Sol[4] * dx * dy + Sol[5] * dy * dy;
                idata.U_comp[i_] += Sol[6] * dx * dx * dx + Sol[7] * dx * dx * dy + Sol[8] * dx * dy * dy + Sol[9] * dy * dy * dy;
                idata.V_comp[i_] = Sol[10] + Sol[11] * dx - Sol[1] * dy + Sol[12] * dx * dx - 2.0 * Sol[3] * dx * dy - 0.5 * Sol[4] * dy * dy;
                idata.V_comp[i_] += Sol[13] * dx * dx * dx - 3.0 * Sol[6] * dx * dx * dy - Sol[7] * dx * dy * dy - Sol[8] * dy * dy * dy / 3.0;
        
                ux = Sol[1] + 2.0 * Sol[3] * dx + Sol[4] * dy + 3.0 * Sol[6] * dx * dx + 2.0 * Sol[7] * dx * dy + Sol[8] * dy * dy;
                uy = Sol[2] + Sol[4] * dx + 2.0 * Sol[5] * dy + Sol[7] * dx * dx + 2.0 * Sol[8] * dx * dy + 3.0 * Sol[9] * dy * dy;
                vx = Sol[11] + 2.0 * Sol[12] * dx - 2.0 * Sol[3] * dy + 3.0 * Sol[13] * dx * dx - 6.0 * Sol[6] * dx * dy - Sol[7] * dy * dy;
                vy = -ux;
        
                idata.X_Visc_Int[i_] = 2.0 * (ux)*idata.psix_[i_] + (uy + vx) * idata.psiy_[i_];
                idata.Y_Visc_Int[i_] = (uy + vx) * idata.psix_[i_] + 2.0 * (vy)*idata.psiy_[i_];
        
                idata.X_Visc_Int[i_] /= sqrt(idata.psix_[i_] * idata.psix_[i_] + idata.psiy_[i_] * idata.psiy_[i_] + 1.0E-12);
                idata.Y_Visc_Int[i_] /= sqrt(idata.psix_[i_] * idata.psix_[i_] + idata.psiy_[i_] * idata.psiy_[i_] + 1.0E-12);
                idata.X_Visc_Int[i_] *= MU_0 / deltax[0];
                idata.Y_Visc_Int[i_] *= MU_0 / deltax[1];
            }
        
            for (i_ = 0; i_ < 2 * Num_; i_++)
            {
                delete[] Mat[i_];
                delete[] Q[i_];
            }
            for (i_ = 0; i_ < 14; i_++)
            {
                delete[] R[i_];
                delete[] Q2T[i_];
                delete[] Q2T_Q[i_];
            }
            for (i_ = 0; i_ < 2 * N_; i_++)
            {
                delete[] Q1[i_];
            }
            for (i_ = 0; i_ < 2 * Num_i_; i_++)
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
        }
        
        void LS_Vel_Second_Order
        (
            bool &sing, 
            int Num_, 
            const amrex::IntVect& cell, 
            amrex::Array4<amrex::Real> const &U,
            double *x_, 
            double *y_, 
            double *U_, 
            double *V_, 
            double *w_
        )
        {
            int i_, j_;
            double **Mat; // Least squares and constraint matrices
            double *Sol, *Vec_;
            double dx, dy, du, dv, nx, ny;
        
            Allocate_2D_R(Mat, 2 * Num_, 5);
            Sol = new double[5];
            Vec_ = new double[2 * Num_];
        
            for (i_ = 0; i_ < Num_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                du = U_[i_];
                dv = V_[i_];
        
                j_ = i_;
                Mat[j_][0] = w_[i_];      // w
                Mat[j_][1] = w_[i_] * dx; // x
                Mat[j_][2] = w_[i_] * dy; // y
                Mat[j_][3] = 0.0;
                Mat[j_][4] = 0.0;
        
                Mat[j_ + Num_][0] = 0.0;          // w
                Mat[j_ + Num_][1] = -w_[i_] * dy; // x
                Mat[j_ + Num_][2] = 0.0;          // y
                Mat[j_ + Num_][3] = w_[i_];
                Mat[j_ + Num_][4] = w_[i_] * dx;
        
                Vec_[j_] = w_[i_] * du;        // x
                Vec_[j_ + Num_] = w_[i_] * dv; // x
            }
            QRdcmp<double> QR(Mat, 2 * Num_, 5);
            QR.solve(Vec_, Sol);
        
            U(cell, 0) = Sol[0];
            U(cell, 1) = Sol[3];
        
            for (i_ = 0; i_ < 2 * Num_; i_++)
            {
                delete[] Mat[i_];
            }
            delete[] Mat;
            delete[] Sol;
            delete[] Vec_;
        }
    } // namespace 
/**********************************************************************/

    namespace AdvectLSIntVel
    {
    
        void Int_LS_Vel_Fourth_Order
        (
            bool &sing, 
            int Num_, int Num_i_, 
            InterceptData &idata,
            amrex::Array4<amrex::Real> const &U,
            double *x_, 
            double *y_, 
            double *U_, 
            double *V_, 
            double *w_,
            const amrex::Real &MU_0,
            const amrex::Real *deltax
        )
        {
            int i_, j_, k_, N_, Size_Loc = 14;
            double **Mat;    // Least squares and constraint matrices
            double **Q, **R; // Q R decomosition of Mat;
            double **Q1, **Q2T, **Q2T_Q, **Q2T_R;
            double *Sol, *Vec_, *Constr_, *Sol_u_, *Sol_w_, *Q1Tb; // *Q2T_QTQ1Tb ;
            double dx, dy, du, dv, nx, ny, dGamma_dt, MU;
        
            N_ = Num_ - Num_i_;
        
            Allocate_2D_R(Mat, 2 * Num_ - Num_i_, Size_Loc);
            Sol = new double[Size_Loc];
            Vec_ = new double[2 * N_];
            Allocate_2D_R(Q, 2 * Num_ - Num_i_, Size_Loc);
            Allocate_2D_R(R, Size_Loc, Size_Loc);
            Allocate_2D_R(Q1, 2 * N_, Size_Loc);
            Allocate_2D_R(Q2T, Size_Loc, Num_i_);
            Allocate_2D_R(Q2T_Q, Size_Loc, Num_i_);
            Allocate_2D_R(Q2T_R, Num_i_, Num_i_);
            Constr_ = new double[Num_i_];
            Sol_u_ = new double[Num_i_];
            Sol_w_ = new double[Num_i_];
            Q1Tb = new double[Size_Loc]; // Q2T_QTQ1Tb = new double[Num_i_] ;
        
            for (i_ = Num_i_; i_ < Num_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                du = U_[i_];
                dv = V_[i_];
        
                j_ = i_ - Num_i_;
        
                Mat[j_][0] = w_[i_];           // w
                Mat[j_][1] = w_[i_] * dx;      // x
                Mat[j_][2] = w_[i_] * dy;      // y
                Mat[j_][3] = w_[i_] * dx * dx; // x^2
                Mat[j_][4] = w_[i_] * dx * dy; // xy
                Mat[j_][5] = w_[i_] * dy * dy; // y^2
        
                Mat[j_][6] = w_[i_] * dx * dx * dx; // v_0
                Mat[j_][7] = w_[i_] * dx * dx * dy; // v_1
                Mat[j_][8] = w_[i_] * dx * dy * dy; // v_11
                Mat[j_][9] = w_[i_] * dy * dy * dy; // v_11
        
                Mat[j_][10] = 0.0; // v_11
                Mat[j_][11] = 0.0; // v_11
                Mat[j_][12] = 0.0; // v_11
                Mat[j_][13] = 0.0; // v_11
        
                Mat[j_ + N_][0] = 0.0;                     // w
                Mat[j_ + N_][1] = -w_[i_] * dy;            // -v_2*dy
                Mat[j_ + N_][2] = 0.0;                     // y
                Mat[j_ + N_][3] = -2.0 * w_[i_] * dx * dy; // v_12 = -2.0*u_11
                Mat[j_ + N_][4] = -0.5 * w_[i_] * dy * dy;
                Mat[j_ + N_][5] = 0.0;
        
                Mat[j_ + N_][6] = -3.0 * w_[i_] * dx * dx * dy;
                Mat[j_ + N_][7] = -w_[i_] * dx * dy * dy;
                Mat[j_ + N_][8] = -w_[i_] * dy * dy * dy / 3.0;
                Mat[j_ + N_][9] = 0.0;
        
                Mat[j_ + N_][10] = w_[i_];           // v_0
                Mat[j_ + N_][11] = w_[i_] * dx;      // v_1
                Mat[j_ + N_][12] = w_[i_] * dx * dx; // v_11
                Mat[j_ + N_][13] = w_[i_] * dx * dx * dx;
        
                Vec_[j_] = w_[i_] * du;      // x
                Vec_[j_ + N_] = w_[i_] * dv; // x
            }
            double TwoPsixPsiy, Psixsq_minus_Psiysq;
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                j_ = 2 * N_ + i_;
        
                TwoPsixPsiy = 2.0 * idata.psix_[i_] * idata.psiy_[i_];
                Psixsq_minus_Psiysq = (idata.psix_[i_] * idata.psix_[i_] - idata.psiy_[i_] * idata.psiy_[i_]);
        
                Mat[j_][0] = 0.0;                                                      // w
                Mat[j_][1] = -2.0 * TwoPsixPsiy;                                       // x
                Mat[j_][2] = Psixsq_minus_Psiysq;                                      // y
                Mat[j_][3] = -4.0 * dx * TwoPsixPsiy - 2.0 * dy * Psixsq_minus_Psiysq; // x^2
                Mat[j_][4] = -2.0 * dy * TwoPsixPsiy + dx * Psixsq_minus_Psiysq;       // xy
                Mat[j_][5] = 2.0 * dy * Psixsq_minus_Psiysq;                           // y^2
        
                Mat[j_][6] = -6.0 * dx * dx * TwoPsixPsiy - 6.0 * dx * dy * Psixsq_minus_Psiysq;       // v_0
                Mat[j_][7] = -4.0 * dx * dy * TwoPsixPsiy + (dx * dx - dy * dy) * Psixsq_minus_Psiysq; // v_1
                Mat[j_][8] = -2.0 * dy * dy * TwoPsixPsiy + 2.0 * dx * dy * Psixsq_minus_Psiysq;       // v_11
                Mat[j_][9] = 3.0 * dy * dy * Psixsq_minus_Psiysq;                                      // v_11
        
                Mat[j_][10] = 0.0;                                 // v_11
                Mat[j_][11] = Psixsq_minus_Psiysq;                 // v_11
                Mat[j_][12] = 2.0 * dx * Psixsq_minus_Psiysq;      // v_11
                Mat[j_][13] = 3.0 * dx * dx * Psixsq_minus_Psiysq; // v_11
        
                Constr_[i_] = 0.0;
            }
            QRdcmp<double> QR(Mat, 2 * Num_ - Num_i_, Size_Loc);
            QR.get_Q(Q);
            QR.get_R(R);
        
            for (j_ = 0; j_ < Size_Loc; j_++)
            {
                for (i_ = 0; i_ < 2 * N_; i_++)
                    Q1[i_][j_] = Q[i_][j_];
                for (i_ = 0; i_ < Num_i_; i_++)
                    Q2T[j_][i_] = Q[2 * N_ + i_][j_];
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
                for (j_ = 0; j_ < 2 * N_; j_++)
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
            for (i_ = 0; i_ < Size_Loc; i_++)
            {
                Vec_[i_] = Q1Tb[i_];
                for (j_ = 0; j_ < Num_i_; j_++)
                    Vec_[i_] -= 0.5 * Q2T[i_][j_] * Sol_w_[j_];
            }
            // solve (Mat_R) w = RHS from above
            for (i_ = Size_Loc - 1; i_ >= 0; i_--)
            {
                Sol[i_] = Vec_[i_] / R[i_][i_];
                for (j_ = i_ + 1; j_ < Size_Loc; j_++)
                    Sol[i_] -= Sol[j_] * R[i_][j_] / R[i_][i_];
            }
        
            U(idata.cellid_, 0) = Sol[0];
            U(idata.cellid_, 1) = Sol[10];
            double ux, uy, vx, vy;
            for (i_ = 0; i_ < Num_i_; i_++)
            {
        
                dx = x_[i_];
                dy = y_[i_];
        
                TwoPsixPsiy = 2.0 * idata.psix_[i_] * idata.psiy_[i_];
                Psixsq_minus_Psiysq = (idata.psix_[i_] * idata.psix_[i_] - idata.psiy_[i_] * idata.psiy_[i_]);
        
                ux = Sol[1] + 2.0 * Sol[3] * dx + Sol[4] * dy + 3.0 * Sol[6] * dx * dx + 2.0 * Sol[7] * dx * dy + Sol[8] * dy * dy;
                uy = Sol[2] + Sol[4] * dx + 2.0 * Sol[5] * dy + Sol[7] * dx * dx + 2.0 * Sol[8] * dx * dy + 3.0 * Sol[9] * dy * dy;
                vx = Sol[11] + 2.0 * Sol[12] * dx - 2.0 * Sol[3] * dy + 3.0 * Sol[13] * dx * dx - 6.0 * Sol[6] * dx * dy - Sol[7] * dy * dy;
                ux /= deltax[0];
                uy /= deltax[1];
                vx /= deltax[0];
                vy = -ux;
        
                //	dGamma_dt = Compute_Gamma_Dot(ux, 0.5*(uy + vx), vy) ;
                //	MU = Get_Viscosity(dGamma_dt) ;
                MU = MU_0;
        
		double S11 = ux;
		double S12 = 0.5*(uy + vx);
		double S22 = vy;
		idata.Gamma_dot_Int[i_] = std::sqrt(2.0 * (S11 * S11 + S22 * S22 + 2.0 * S12 * S12 + (S11 + S22) * (S11 + S22)));
                idata.tan_shear_[i_] = Psixsq_minus_Psiysq * (uy + vx) + TwoPsixPsiy * (vy - ux);
                idata.norm_shear_[i_] = 2.0 * ux * idata.psix_[i_] * idata.psix_[i_] + 2.0 * vy * idata.psiy_[i_] * idata.psiy_[i_] + 2.0 * (vx + uy) * idata.psix_[i_] * idata.psiy_[i_];
                idata.norm_shear_[i_] /= (idata.psix_[i_] * idata.psix_[i_] + idata.psiy_[i_] * idata.psiy_[i_] + 1.0E-12);
        
                idata.norm_shear_[i_] *= MU;
            }
        
            for (i_ = 0; i_ < 2 * Num_ - Num_i_; i_++)
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
            for (i_ = 0; i_ < 2 * N_; i_++)
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
        }
        
        void Int_LS_Vel_Third_Order
        (
            bool &sing, 
            int Num_, int Num_i_, 
            InterceptData &idata,
            amrex::Array4<amrex::Real> const &U, 
            double *x_, 
            double *y_, 
            double *U_, 
            double *V_, 
            double *w_,
            const amrex::Real &MU_0,
            const amrex::Real *deltax
        )
        {
            int i_, j_, k_, N_, Size_Loc = 9;
            double **Mat;    // Least squares and constraint matrices
            double **Q, **R; // Q R decomosition of Mat;
            double **Q1, **Q2T, **Q2T_Q, **Q2T_R;
            double *Sol, *Vec_, *Constr_, *Sol_u_, *Sol_w_, *Q1Tb; // *Q2T_QTQ1Tb ;
            double dx, dy, du, dv, nx, ny, dGamma_dt, MU;
        
            N_ = Num_ - Num_i_;
        
            Allocate_2D_R(Mat, 2 * Num_ - Num_i_, Size_Loc);
            Sol = new double[Size_Loc];
            Vec_ = new double[2 * N_];
            Allocate_2D_R(Q, 2 * Num_ - Num_i_, Size_Loc);
            Allocate_2D_R(R, Size_Loc, Size_Loc);
            Allocate_2D_R(Q1, 2 * N_, Size_Loc);
            Allocate_2D_R(Q2T, Size_Loc, Num_i_);
            Allocate_2D_R(Q2T_Q, Size_Loc, Num_i_);
            Allocate_2D_R(Q2T_R, Num_i_, Num_i_);
            Constr_ = new double[Num_i_];
            Sol_u_ = new double[Num_i_];
            Sol_w_ = new double[Num_i_];
            Q1Tb = new double[Size_Loc]; // Q2T_QTQ1Tb = new double[Num_i_] ;
        
            for (i_ = Num_i_; i_ < Num_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                du = U_[i_];
                dv = V_[i_];
        
                j_ = i_ - Num_i_;
        
                Mat[j_][0] = w_[i_];           // w
                Mat[j_][1] = w_[i_] * dx;      // x
                Mat[j_][2] = w_[i_] * dy;      // y
                Mat[j_][3] = w_[i_] * dx * dx; // x^2
                Mat[j_][4] = w_[i_] * dx * dy; // xy
                Mat[j_][5] = w_[i_] * dy * dy; // y^2
        
                Mat[j_][6] = 0.0; // v_0
                Mat[j_][7] = 0.0; // v_1
                Mat[j_][8] = 0.0; // v_11
        
                Mat[j_ + N_][0] = 0.0;                     // w
                Mat[j_ + N_][1] = -w_[i_] * dy;            // -v_2*dy
                Mat[j_ + N_][2] = 0.0;                     // y
                Mat[j_ + N_][3] = -2.0 * w_[i_] * dx * dy; // v_12 = -2.0*u_11
                Mat[j_ + N_][4] = -0.5 * w_[i_] * dy * dy;
                Mat[j_ + N_][5] = 0.0;
        
                Mat[j_ + N_][6] = w_[i_];           // v_0
                Mat[j_ + N_][7] = w_[i_] * dx;      // v_1
                Mat[j_ + N_][8] = w_[i_] * dx * dx; // v_11
        
                Vec_[j_] = w_[i_] * du;      // x
                Vec_[j_ + N_] = w_[i_] * dv; // x
            }
            double TwoPsixPsiy, Psixsq_minus_Psiysq;
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                j_ = 2 * N_ + i_;
        
                TwoPsixPsiy = 2.0 * idata.psix_[i_] * idata.psiy_[i_];
                Psixsq_minus_Psiysq = (idata.psix_[i_] * idata.psix_[i_] - idata.psiy_[i_] * idata.psiy_[i_]);
        
                Mat[j_][0] = 0.0;                                                      // w
                Mat[j_][1] = -2.0 * TwoPsixPsiy;                                       // x
                Mat[j_][2] = Psixsq_minus_Psiysq;                                      // y
                Mat[j_][3] = -4.0 * dx * TwoPsixPsiy - 2.0 * dy * Psixsq_minus_Psiysq; // x^2
                Mat[j_][4] = -2.0 * dy * TwoPsixPsiy + dx * Psixsq_minus_Psiysq;       // xy
                Mat[j_][5] = 2.0 * dy * Psixsq_minus_Psiysq;                           // y^2
        
                Mat[j_][6] = 0.0;                            // v_11
                Mat[j_][7] = Psixsq_minus_Psiysq;            // v_11
                Mat[j_][8] = 2.0 * dx * Psixsq_minus_Psiysq; // v_11
        
                Constr_[i_] = 0.0;
            }
            QRdcmp<double> QR(Mat, 2 * Num_ - Num_i_, Size_Loc);
            QR.get_Q(Q);
            QR.get_R(R);
        
            for (j_ = 0; j_ < Size_Loc; j_++)
            {
                for (i_ = 0; i_ < 2 * N_; i_++)
                    Q1[i_][j_] = Q[i_][j_];
                for (i_ = 0; i_ < Num_i_; i_++)
                    Q2T[j_][i_] = Q[2 * N_ + i_][j_];
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
                for (j_ = 0; j_ < 2 * N_; j_++)
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
            for (i_ = 0; i_ < Size_Loc; i_++)
            {
                Vec_[i_] = Q1Tb[i_];
                for (j_ = 0; j_ < Num_i_; j_++)
                    Vec_[i_] -= 0.5 * Q2T[i_][j_] * Sol_w_[j_];
            }
            // solve (Mat_R) w = RHS from above
            for (i_ = Size_Loc - 1; i_ >= 0; i_--)
            {
                Sol[i_] = Vec_[i_] / R[i_][i_];
                for (j_ = i_ + 1; j_ < Size_Loc; j_++)
                    Sol[i_] -= Sol[j_] * R[i_][j_] / R[i_][i_];
            }
        
            U(idata.cellid_, 0) = Sol[0];
            U(idata.cellid_, 1) = Sol[6];
            double ux, uy, vx, vy;
            for (i_ = 0; i_ < Num_i_; i_++)
            {
        
                dx = x_[i_];
                dy = y_[i_];
        
                TwoPsixPsiy = 2.0 * idata.psix_[i_] * idata.psiy_[i_];
                Psixsq_minus_Psiysq = (idata.psix_[i_] * idata.psix_[i_] - idata.psiy_[i_] * idata.psiy_[i_]);
        
                ux = Sol[1] + 2.0 * Sol[3] * dx + Sol[4] * dy;
                uy = Sol[2] + Sol[4] * dx + 2.0 * Sol[5] * dy;
                vx = Sol[7] + 2.0 * Sol[8] * dx - 2.0 * Sol[3] * dy;
        
                ux /= deltax[0];
                uy /= deltax[1];
                vx /= deltax[0];
                vy = -ux;
        
                //	dGamma_dt = Compute_Gamma_Dot(ux, 0.5*(uy + vx), vy) ;
                //	MU = Get_Viscosity(dGamma_dt) ;
                MU = MU_0;
        
                double S11 = ux;
                double S12 = 0.5*(uy + vx);
                double S22 = vy;
                idata.Gamma_dot_Int[i_] = std::sqrt(2.0 * (S11 * S11 + S22 * S22 + 2.0 * S12 * S12 + (S11 + S22) * (S11 + S22)));
                idata.tan_shear_[i_] = Psixsq_minus_Psiysq * (uy + vx) + TwoPsixPsiy * (vy - ux);
                idata.norm_shear_[i_] = 2.0 * ux * idata.psix_[i_] * idata.psix_[i_] + 2.0 * vy * idata.psiy_[i_] * idata.psiy_[i_] + 2.0 * (vx + uy) * idata.psix_[i_] * idata.psiy_[i_];
                idata.norm_shear_[i_] /= (idata.psix_[i_] * idata.psix_[i_] + idata.psiy_[i_] * idata.psiy_[i_] + 1.0E-12);
        
                idata.norm_shear_[i_] *= MU;
            }
        
            for (i_ = 0; i_ < 2 * Num_ - Num_i_; i_++)
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
            for (i_ = 0; i_ < 2 * N_; i_++)
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
        }
        
        void Int_LS_Vel_Second_Order
        (
            bool &sing, 
            int Num_, int Num_i_,
            InterceptData &idata,
            amrex::Array4<amrex::Real> const &U, 
            double *x_, 
            double *y_, 
            double *U_, 
            double *V_, 
            double *w_,
            const amrex::Real &MU_0,
            const amrex::Real *deltax
        )
        {
            int i_, j_, k_, N_, Size_Loc = 5;
            double **Mat;    // Least squares and constraint matrices
            double **Q, **R; // Q R decomosition of Mat;
            double **Q1, **Q2T, **Q2T_Q, **Q2T_R;
            double *Sol, *Vec_, *Constr_, *Sol_u_, *Sol_w_, *Q1Tb; // *Q2T_QTQ1Tb ;
            double dx, dy, du, dv, nx, ny, dGamma_dt, MU;
        
            N_ = Num_ - Num_i_;
        
            Allocate_2D_R(Mat, 2 * Num_ - Num_i_, Size_Loc);
            Sol = new double[Size_Loc];
            Vec_ = new double[2 * N_];
            Allocate_2D_R(Q, 2 * Num_ - Num_i_, Size_Loc);
            Allocate_2D_R(R, Size_Loc, Size_Loc);
            Allocate_2D_R(Q1, 2 * N_, Size_Loc);
            Allocate_2D_R(Q2T, Size_Loc, Num_i_);
            Allocate_2D_R(Q2T_Q, Size_Loc, Num_i_);
            Allocate_2D_R(Q2T_R, Num_i_, Num_i_);
            Constr_ = new double[Num_i_];
            Sol_u_ = new double[Num_i_];
            Sol_w_ = new double[Num_i_];
            Q1Tb = new double[Size_Loc]; // Q2T_QTQ1Tb = new double[Num_i_] ;
        
            for (i_ = Num_i_; i_ < Num_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                du = U_[i_];
                dv = V_[i_];
        
                j_ = i_ - Num_i_;
        
                Mat[j_][0] = w_[i_];      // w
                Mat[j_][1] = w_[i_] * dx; // x
                Mat[j_][2] = w_[i_] * dy; // y
        
                Mat[j_][3] = 0.0; // v_0
                Mat[j_][4] = 0.0; // v_1
        
                Mat[j_ + N_][0] = 0.0;          // w
                Mat[j_ + N_][1] = -w_[i_] * dy; // -v_2*dy
                Mat[j_ + N_][2] = 0.0;          // y
        
                Mat[j_ + N_][3] = w_[i_];      // v_0
                Mat[j_ + N_][4] = w_[i_] * dx; // v_1
        
                Vec_[j_] = w_[i_] * du;      // x
                Vec_[j_ + N_] = w_[i_] * dv; // x
            }
            double TwoPsixPsiy, Psixsq_minus_Psiysq;
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                j_ = 2 * N_ + i_;
        
                TwoPsixPsiy = 2.0 * idata.psix_[i_] * idata.psiy_[i_];
                Psixsq_minus_Psiysq = (idata.psix_[i_] * idata.psix_[i_] - idata.psiy_[i_] * idata.psiy_[i_]);
        
                Mat[j_][0] = 0.0;                 // w
                Mat[j_][1] = -2.0 * TwoPsixPsiy;  // x
                Mat[j_][2] = Psixsq_minus_Psiysq; // y
        
                Mat[j_][3] = 0.0;                 // v_11
                Mat[j_][4] = Psixsq_minus_Psiysq; // v_11
        
                Constr_[i_] = 0.0;
            }
            QRdcmp<double> QR(Mat, 2 * Num_ - Num_i_, Size_Loc);
            QR.get_Q(Q);
            QR.get_R(R);
        
            for (j_ = 0; j_ < Size_Loc; j_++)
            {
                for (i_ = 0; i_ < 2 * N_; i_++)
                    Q1[i_][j_] = Q[i_][j_];
                for (i_ = 0; i_ < Num_i_; i_++)
                    Q2T[j_][i_] = Q[2 * N_ + i_][j_];
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
                for (j_ = 0; j_ < 2 * N_; j_++)
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
            for (i_ = 0; i_ < Size_Loc; i_++)
            {
                Vec_[i_] = Q1Tb[i_];
                for (j_ = 0; j_ < Num_i_; j_++)
                    Vec_[i_] -= 0.5 * Q2T[i_][j_] * Sol_w_[j_];
            }
            // solve (Mat_R) w = RHS from above
            for (i_ = Size_Loc - 1; i_ >= 0; i_--)
            {
                Sol[i_] = Vec_[i_] / R[i_][i_];
                for (j_ = i_ + 1; j_ < Size_Loc; j_++)
                    Sol[i_] -= Sol[j_] * R[i_][j_] / R[i_][i_];
            }
        
            U(idata.cellid_, 0) = Sol[0];
            U(idata.cellid_, 1) = Sol[3];
            double ux, uy, vx, vy;
            for (i_ = 0; i_ < Num_i_; i_++)
            {
        
                dx = x_[i_];
                dy = y_[i_];
        
                TwoPsixPsiy = 2.0 * idata.psix_[i_] * idata.psiy_[i_];
                Psixsq_minus_Psiysq = (idata.psix_[i_] * idata.psix_[i_] - idata.psiy_[i_] * idata.psiy_[i_]);
        
                ux = Sol[1];
                uy = Sol[2];
                vx = Sol[4];
        
                ux /= deltax[0];
                uy /= deltax[1];
                vx /= deltax[0];
                vy = -ux;
        
                //	dGamma_dt = Compute_Gamma_Dot(ux, 0.5*(uy + vx), vy) ;
                //	MU = Get_Viscosity(dGamma_dt) ;
                MU = MU_0;
        
                double S11 = ux;
                double S12 = 0.5*(uy + vx);
                double S22 = vy;
                idata.Gamma_dot_Int[i_] = std::sqrt(2.0 * (S11 * S11 + S22 * S22 + 2.0 * S12 * S12 + (S11 + S22) * (S11 + S22)));
                idata.tan_shear_[i_] = Psixsq_minus_Psiysq * (uy + vx) + TwoPsixPsiy * (vy - ux);
                idata.norm_shear_[i_] = 2.0 * ux * idata.psix_[i_] * idata.psix_[i_] + 2.0 * vy * idata.psiy_[i_] * idata.psiy_[i_] + 2.0 * (vx + uy) * idata.psix_[i_] * idata.psiy_[i_];
                idata.norm_shear_[i_] /= (idata.psix_[i_] * idata.psix_[i_] + idata.psiy_[i_] * idata.psiy_[i_] + 1.0E-12);
        
                idata.norm_shear_[i_] *= MU;
            }
        
            for (i_ = 0; i_ < 2 * Num_ - Num_i_; i_++)
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
            for (i_ = 0; i_ < 2 * N_; i_++)
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
        }
        
        namespace Axisymmetric
        {
            void Vel_LSQ_Parameters(
                double &U_,
                double &V_,
                double &x_,
                double &y_,
                double &w_,
                int i1, int j1,
                int i, int j,
                amrex::Array4<amrex::Real const> const &U, /// cc vel with 2 comp
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
                U_ = U(i1, j1, 0, 0);
                V_ = U(i1, j1, 0, 1);
                if(j1 < domain.smallEnd(1))
                {
                    //amrex::Print()<<"i = "<<i<<" , j = "<<j<<" , i1 = "<<i1<<" , j1 = "<<j1<<'\n';
                    j1 = domain.smallEnd(1) - j1 - 1;
                    //amrex::Print()<<"new j1 = "<<j1<<'\n';
                    U_ = U(i1, j1, 0, 0);
                    V_ = -1.0*U(i1, j1, 0, 1);
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

        
            void Int_LS_Vel_Fourth_Order
            (
                bool &sing,
                int Num_, int Num_i_,
                InterceptData &idata,
                amrex::Array4<amrex::Real> const &U,
                double *x_,
                double *y_,
                double *U_,
                double *V_,
                double *w_,
                // const amrex::Real &MU_0,
                Viscosity& visc,
                const amrex::Real *deltax,
                const amrex::Real *prob_lo
            )
            {
                int i_, j_, k_, N_, Size_Loc = 10;
                double **Mat;    // Least squares and constraint matrices
                double **Q, **R; // Q R decomosition of Mat;
                double **Q1, **Q2T, **Q2T_Q, **Q2T_R;
                double *Sol, *Vec_, *Constr_, *Sol_u_, *Sol_w_, *Q1Tb; // *Q2T_QTQ1Tb ;
                double dx, dy, du, dv, nx, ny, dGamma_dt, MU, y0;
            
                N_ = Num_ - Num_i_;
            
                Allocate_2D_R(Mat, 2 * Num_ - Num_i_, Size_Loc);
                Sol = new double[Size_Loc];
                Vec_ = new double[2 * N_];
                Allocate_2D_R(Q, 2 * Num_ - Num_i_, Size_Loc);
                Allocate_2D_R(R, Size_Loc, Size_Loc);
                Allocate_2D_R(Q1, 2 * N_, Size_Loc);
                Allocate_2D_R(Q2T, Size_Loc, Num_i_);
                Allocate_2D_R(Q2T_Q, Size_Loc, Num_i_);
                Allocate_2D_R(Q2T_R, Num_i_, Num_i_);
                Constr_ = new double[Num_i_];
                Sol_u_ = new double[Num_i_];
                Sol_w_ = new double[Num_i_];
                Q1Tb = new double[Size_Loc]; // Q2T_QTQ1Tb = new double[Num_i_] ;
            
                amrex::Real y = prob_lo[1] + deltax[1] * (idata.cellid_[1] + 0.5);
                y0 = y / deltax[1];
            
                for (i_ = Num_i_; i_ < Num_; i_++)
                {
                    dx = x_[i_];
                    dy = y_[i_];
                    du = U_[i_];
                    dv = V_[i_];
            
                    j_ = i_ - Num_i_;
            
                    Mat[j_][0] = w_[i_];           // w
                    Mat[j_][1] = w_[i_] * dx;      // x
                    Mat[j_][2] = w_[i_] * dy;      // y
                    Mat[j_][3] = w_[i_] * dx * dx; // x^2
                    Mat[j_][4] = w_[i_] * dx * dy; // xy
                    Mat[j_][5] = w_[i_] * dy * dy; // y^2
            
                    Mat[j_][6] = w_[i_] * dx * dx * dx; // v_0
                    Mat[j_][7] = w_[i_] * dx * dx * dy; // v_1
                    Mat[j_][8] = w_[i_] * dx * dy * dy; // v_11
                    Mat[j_][9] = w_[i_] * dy * dy * dy; // v_11
            
                    Mat[j_ + N_][0] = 0.0; // w
                    Mat[j_ + N_][1] = -w_[i_] * (y0 + dy) / 2.0;
                    Mat[j_ + N_][2] = 0.0; // y
            
                    Mat[j_ + N_][3] = -w_[i_] * dx * (y0 + dy);
                    Mat[j_ + N_][4] = w_[i_] * (y0 + dy) * (y0 - 2.0 * dy) / 6.0;
                    Mat[j_ + N_][5] = 0.0;
            
                    Mat[j_ + N_][6] = -3.0 * w_[i_] * dx * dx * (y0 + dy) / 2.0;
                    Mat[j_ + N_][7] = w_[i_] * dx * (y0 + dy) * (y0 - 2.0 * dy) / 3.0;
                    Mat[j_ + N_][8] = -w_[i_] * (y0 + dy) * (y0 * y0 - 2.0 * y0 * dy + 3.0 * dy * dy) / 12.0;
                    Mat[j_ + N_][9] = 0.0;
            
                    Vec_[j_] = w_[i_] * du;      // x
                    Vec_[j_ + N_] = w_[i_] * dv; // x
                }
                double TwoPsixPsiy, Psixsq_minus_Psiysq;
                for (i_ = 0; i_ < Num_i_; i_++)
                {
                    dx = x_[i_];
                    dy = y_[i_];
                    j_ = 2 * N_ + i_;
            
                    TwoPsixPsiy = 2.0 * idata.psix_[i_] * idata.psiy_[i_];
                    Psixsq_minus_Psiysq = (idata.psix_[i_] * idata.psix_[i_] - idata.psiy_[i_] * idata.psiy_[i_]);
            
                    Mat[j_][0] = 0.0;                                                            // w
                    Mat[j_][1] = -3.0 * TwoPsixPsiy / 2.0;                                       // x
                    Mat[j_][2] = Psixsq_minus_Psiysq;                                            // y
                    Mat[j_][3] = -3.0 * dx * TwoPsixPsiy - (y0 + dy) * Psixsq_minus_Psiysq;      // x^2
                    Mat[j_][4] = -(y0 + 10 * dy) * TwoPsixPsiy / 6.0 + dx * Psixsq_minus_Psiysq; // xy
                    Mat[j_][5] = 2.0 * dy * Psixsq_minus_Psiysq;                                 // y^2
            
                    Mat[j_][6] = -9.0 * dx * dx * TwoPsixPsiy / 2.0 - 3.0 * dx * (y0 + dy) * Psixsq_minus_Psiysq;                                          // v_0
                    Mat[j_][7] = -dx * (y0 + 10.0 * dy) * TwoPsixPsiy / 3.0 + (dx * dx + (y0 * y0 - dy * y0 - 2.0 * dy * dy) / 3.0) * Psixsq_minus_Psiysq; // v_1
                    Mat[j_][8] = ((y0 * y0 - 2.0 * dy * y0 - 21.0 * dy * dy) / 12.0) * TwoPsixPsiy + 2.0 * dx * dy * Psixsq_minus_Psiysq;                  // v_11
                    Mat[j_][9] = 3.0 * dy * dy * Psixsq_minus_Psiysq;                                                                                      // v_11
            
                    Constr_[i_] = 0.0;
                }
                QRdcmp<double> QR(Mat, 2 * Num_ - Num_i_, Size_Loc);
                QR.get_Q(Q);
                QR.get_R(R);
            
                for (j_ = 0; j_ < Size_Loc; j_++)
                {
                    for (i_ = 0; i_ < 2 * N_; i_++)
                        Q1[i_][j_] = Q[i_][j_];
                    for (i_ = 0; i_ < Num_i_; i_++)
                        Q2T[j_][i_] = Q[2 * N_ + i_][j_];
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
                    for (j_ = 0; j_ < 2 * N_; j_++)
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
                for (i_ = 0; i_ < Size_Loc; i_++)
                {
                    Vec_[i_] = Q1Tb[i_];
                    for (j_ = 0; j_ < Num_i_; j_++)
                        Vec_[i_] -= 0.5 * Q2T[i_][j_] * Sol_w_[j_];
                }
                // solve (Mat_R) w = RHS from above
                for (i_ = Size_Loc - 1; i_ >= 0; i_--)
                {
                    Sol[i_] = Vec_[i_] / R[i_][i_];
                    for (j_ = i_ + 1; j_ < Size_Loc; j_++)
                        Sol[i_] -= Sol[j_] * R[i_][j_] / R[i_][i_];
                }
            
                U(idata.cellid_, 0) = Sol[0];
                U(idata.cellid_, 1) = -Sol[1] * y0 / 2.0 + Sol[4] * y0 * y0 / 6.0 - Sol[8] * y0 * y0 * y0 / 12.0;
                double ux, uy, vx, vy;
                for (i_ = 0; i_ < Num_i_; i_++)
                {
            
                    dx = x_[i_];
                    dy = y_[i_];
            
                    TwoPsixPsiy = 2.0 * idata.psix_[i_] * idata.psiy_[i_];
                    Psixsq_minus_Psiysq = (idata.psix_[i_] * idata.psix_[i_] - idata.psiy_[i_] * idata.psiy_[i_]);
            
                    ux = Sol[1] + 2.0 * Sol[3] * dx + Sol[4] * dy + 3.0 * Sol[6] * dx * dx + 2.0 * Sol[7] * dx * dy + Sol[8] * dy * dy;
                    uy = Sol[2] + Sol[4] * dx + 2.0 * Sol[5] * dy + Sol[7] * dx * dx + 2.0 * Sol[8] * dx * dy + 3.0 * Sol[9] * dy * dy;
                    vx = -Sol[3] * (y0 + dy) - 3.0 * Sol[6] * dx * (y0 + dy) + Sol[7] * (y0 + dy) * (y0 - 2.0 * dy) / 3.0;
                    vy = -Sol[1] / 2.0 - Sol[3] * dx + Sol[4] * (-y0 - 4.0 * dy) / 6.0 - 3.0 * dx * dx * Sol[6] / 2.0 - Sol[7] * dx * (y0 + 4.0 * dy) / 3.0 + Sol[8] * (y0 * y0 / 12.0 - y0 * dy / 6.0 - 3.0 * dy * dy / 4.0);
            
                    ux /= deltax[0];
                    uy /= deltax[1];
                    vx /= deltax[0];
                    vy /= deltax[1];
            
                    dGamma_dt = visc.ComputeGammaDot(ux, 0.5 * (uy + vx), vy);
                    MU = visc.GetViscosity(dGamma_dt);
                    //amrex::Print()<<"MU = "<<MU<<'\n';
                    // MU = MU_0;
            
                    idata.Gamma_dot_Int[i_] = dGamma_dt;
                    idata.tan_shear_[i_] = Psixsq_minus_Psiysq * (uy + vx) + TwoPsixPsiy * (vy - ux);
                    idata.norm_shear_[i_] = 2.0 * ux * idata.psix_[i_] * idata.psix_[i_] + 2.0 * vy * idata.psiy_[i_] * idata.psiy_[i_] + 2.0 * (vx + uy) * idata.psix_[i_] * idata.psiy_[i_];
                    idata.norm_shear_[i_] /= (idata.psix_[i_] * idata.psix_[i_] + idata.psiy_[i_] * idata.psiy_[i_] + 1.0E-12);
   
                    //idata.norm_shear_[i_] *= MU;
            
                    //amrex::IntVect test_iv(AMREX_D_DECL(124, 139, 0));
                    //if(idata.cellid_ == test_iv)
                    //{
                    //    amrex::Print()<<"ux = "<<ux<<" , uy = "<<uy<<'\n';
                    //    amrex::Print()<<"vx = "<<vx<<" , vy = "<<vy<<'\n';
                    //    amrex::Print()<<"MU = "<<MU<<'\n';
                    //    amrex::Print()<<"idata.psix_[i_] = "<<idata.psix_[i_]<<" , idata.psiy_[i_] = "<<idata.psiy_[i_] <<'\n';
                    //    amrex::Print()<<"idata.norm_shear_[i_] = "<<idata.norm_shear_[i_]<<'\n';
                    //}

                }
            
                for (i_ = 0; i_ < 2 * Num_ - Num_i_; i_++)
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
                for (i_ = 0; i_ < 2 * N_; i_++)
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
            }
            
            void Int_LS_Vel_Third_Order
            (
                bool &sing,
                int Num_, int Num_i_,
                InterceptData &idata,
                amrex::Array4<amrex::Real> const &U,
                double *x_,
                double *y_,
                double *U_,
                double *V_,
                double *w_,
                // const amrex::Real &MU_0,
                Viscosity& visc,
                const amrex::Real *deltax,
                const amrex::Real *prob_lo
            )
            {
                int i_, j_, k_, N_, Size_Loc = 6;
                double **Mat;    // Least squares and constraint matrices
                double **Q, **R; // Q R decomosition of Mat;
                double **Q1, **Q2T, **Q2T_Q, **Q2T_R;
                double *Sol, *Vec_, *Constr_, *Sol_u_, *Sol_w_, *Q1Tb; // *Q2T_QTQ1Tb ;
                double dx, dy, du, dv, nx, ny, dGamma_dt, MU, y0;
            
                N_ = Num_ - Num_i_;
            
                Allocate_2D_R(Mat, 2 * Num_ - Num_i_, Size_Loc);
                Sol = new double[Size_Loc];
                Vec_ = new double[2 * N_];
                Allocate_2D_R(Q, 2 * Num_ - Num_i_, Size_Loc);
                Allocate_2D_R(R, Size_Loc, Size_Loc);
                Allocate_2D_R(Q1, 2 * N_, Size_Loc);
                Allocate_2D_R(Q2T, Size_Loc, Num_i_);
                Allocate_2D_R(Q2T_Q, Size_Loc, Num_i_);
                Allocate_2D_R(Q2T_R, Num_i_, Num_i_);
                Constr_ = new double[Num_i_];
                Sol_u_ = new double[Num_i_];
                Sol_w_ = new double[Num_i_];
                Q1Tb = new double[Size_Loc]; // Q2T_QTQ1Tb = new double[Num_i_] ;
            
                amrex::Real y = prob_lo[1] + deltax[1] * (idata.cellid_[1] + 0.5);
                y0 = y / deltax[1];
            
                for (i_ = Num_i_; i_ < Num_; i_++)
                {
                    dx = x_[i_];
                    dy = y_[i_];
                    du = U_[i_];
                    dv = V_[i_];
            
                    j_ = i_ - Num_i_;
            
                    Mat[j_][0] = w_[i_];           // w
                    Mat[j_][1] = w_[i_] * dx;      // x
                    Mat[j_][2] = w_[i_] * dy;      // y
                    Mat[j_][3] = w_[i_] * dx * dx; // x^2
                    Mat[j_][4] = w_[i_] * dx * dy; // xy
                    Mat[j_][5] = w_[i_] * dy * dy; // y^2
            
                    Mat[j_ + N_][0] = 0.0; // w
                    Mat[j_ + N_][1] = -w_[i_] * (y0 + dy) / 2.0;
                    Mat[j_ + N_][2] = 0.0; // y
            
                    Mat[j_ + N_][3] = -w_[i_] * dx * (y0 + dy);
                    Mat[j_ + N_][4] = w_[i_] * (y0 + dy) * (y0 - 2.0 * dy) / 6.0;
                    Mat[j_ + N_][5] = 0.0;
            
                    Vec_[j_] = w_[i_] * du;      // x
                    Vec_[j_ + N_] = w_[i_] * dv; // x
                }
                double TwoPsixPsiy, Psixsq_minus_Psiysq;
                for (i_ = 0; i_ < Num_i_; i_++)
                {
                    dx = x_[i_];
                    dy = y_[i_];
                    j_ = 2 * N_ + i_;
            
                    TwoPsixPsiy = 2.0 * idata.psix_[i_] * idata.psiy_[i_];
                    Psixsq_minus_Psiysq = (idata.psix_[i_] * idata.psix_[i_] - idata.psiy_[i_] * idata.psiy_[i_]);
            
                    Mat[j_][0] = 0.0;                                                            // w
                    Mat[j_][1] = -3.0 * TwoPsixPsiy / 2.0;                                       // x
                    Mat[j_][2] = Psixsq_minus_Psiysq;                                            // y
                    Mat[j_][3] = -3.0 * dx * TwoPsixPsiy - (y0 + dy) * Psixsq_minus_Psiysq;      // x^2
                    Mat[j_][4] = -(y0 + 10 * dy) * TwoPsixPsiy / 6.0 + dx * Psixsq_minus_Psiysq; // xy
                    Mat[j_][5] = 2.0 * dy * Psixsq_minus_Psiysq;                                 // y^2
            
                    Constr_[i_] = 0.0;
                }
                QRdcmp<double> QR(Mat, 2 * Num_ - Num_i_, Size_Loc);
                QR.get_Q(Q);
                QR.get_R(R);
            
                for (j_ = 0; j_ < Size_Loc; j_++)
                {
                    for (i_ = 0; i_ < 2 * N_; i_++)
                        Q1[i_][j_] = Q[i_][j_];
                    for (i_ = 0; i_ < Num_i_; i_++)
                        Q2T[j_][i_] = Q[2 * N_ + i_][j_];
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
                    for (j_ = 0; j_ < 2 * N_; j_++)
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
                for (i_ = 0; i_ < Size_Loc; i_++)
                {
                    Vec_[i_] = Q1Tb[i_];
                    for (j_ = 0; j_ < Num_i_; j_++)
                        Vec_[i_] -= 0.5 * Q2T[i_][j_] * Sol_w_[j_];
                }
                // solve (Mat_R) w = RHS from above
                for (i_ = Size_Loc - 1; i_ >= 0; i_--)
                {
                    Sol[i_] = Vec_[i_] / R[i_][i_];
                    for (j_ = i_ + 1; j_ < Size_Loc; j_++)
                        Sol[i_] -= Sol[j_] * R[i_][j_] / R[i_][i_];
                }
            
                U(idata.cellid_, 0) = Sol[0];
                U(idata.cellid_, 1) = -Sol[1] * y0 / 2.0 + Sol[4] * y0 * y0 / 6.0;
                double ux, uy, vx, vy;
                for (i_ = 0; i_ < Num_i_; i_++)
                {
            
                    dx = x_[i_];
                    dy = y_[i_];
            
                    TwoPsixPsiy = 2.0 * idata.psix_[i_] * idata.psiy_[i_];
                    Psixsq_minus_Psiysq = (idata.psix_[i_] * idata.psix_[i_] - idata.psiy_[i_] * idata.psiy_[i_]);
            
                    ux = Sol[1] + 2.0 * Sol[3] * dx + Sol[4] * dy;
                    uy = Sol[2] + Sol[4] * dx + 2.0 * Sol[5] * dy;
                    vx = -Sol[3] * (y0 + dy);
                    vy = -Sol[1] / 2.0 - Sol[3] * dx + Sol[4] * (-y0 - 4.0 * dy) / 6.0;
            
                    ux /= deltax[0];
                    uy /= deltax[1];
                    vx /= deltax[0];
                    vy /= deltax[1];
            
                    //TODO varoable viscosity
                    dGamma_dt = visc.ComputeGammaDot(ux, 0.5 * (uy + vx), vy);
                    MU = visc.GetViscosity(dGamma_dt);
                    // MU = MU_0;
            
                    idata.Gamma_dot_Int[i_] = dGamma_dt;
                    idata.tan_shear_[i_] = Psixsq_minus_Psiysq * (uy + vx) + TwoPsixPsiy * (vy - ux);
                    idata.norm_shear_[i_] = 2.0 * ux * idata.psix_[i_] * idata.psix_[i_] + 2.0 * vy * idata.psiy_[i_] * idata.psiy_[i_] + 2.0 * (vx + uy) * idata.psix_[i_] * idata.psiy_[i_];
                    idata.norm_shear_[i_] /= (idata.psix_[i_] * idata.psix_[i_] + idata.psiy_[i_] * idata.psiy_[i_] + 1.0E-12);
            
                    //idata.norm_shear_[i_] *= MU;
                }
            
                for (i_ = 0; i_ < 2 * Num_ - Num_i_; i_++)
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
                for (i_ = 0; i_ < 2 * N_; i_++)
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
            }
            
            void Int_LS_Vel_Second_Order
            (
                bool &sing,
                int Num_, int Num_i_,
                InterceptData &idata,
                amrex::Array4<amrex::Real> const &U,
                double *x_,
                double *y_,
                double *U_,
                double *V_,
                double *w_,
                // const amrex::Real &MU_0,
                Viscosity& visc,
                const amrex::Real *deltax,
                const amrex::Real *prob_lo
            )
            {
                int i_, j_, k_, N_, Size_Loc = 3;
                double **Mat;    // Least squares and constraint matrices
                double **Q, **R; // Q R decomosition of Mat;
                double **Q1, **Q2T, **Q2T_Q, **Q2T_R;
                double *Sol, *Vec_, *Constr_, *Sol_u_, *Sol_w_, *Q1Tb; // *Q2T_QTQ1Tb ;
                double dx, dy, du, dv, nx, ny, dGamma_dt, MU, y0;
            
                N_ = Num_ - Num_i_;
            
                Allocate_2D_R(Mat, 2 * Num_ - Num_i_, Size_Loc);
                Sol = new double[Size_Loc];
                Vec_ = new double[2 * N_];
                Allocate_2D_R(Q, 2 * Num_ - Num_i_, Size_Loc);
                Allocate_2D_R(R, Size_Loc, Size_Loc);
                Allocate_2D_R(Q1, 2 * N_, Size_Loc);
                Allocate_2D_R(Q2T, Size_Loc, Num_i_);
                Allocate_2D_R(Q2T_Q, Size_Loc, Num_i_);
                Allocate_2D_R(Q2T_R, Num_i_, Num_i_);
                Constr_ = new double[Num_i_];
                Sol_u_ = new double[Num_i_];
                Sol_w_ = new double[Num_i_];
                Q1Tb = new double[Size_Loc]; // Q2T_QTQ1Tb = new double[Num_i_] ;
            
                amrex::Real y = prob_lo[1] + deltax[1] * (idata.cellid_[1] + 0.5);
                y0 = y / deltax[1];
            
                for (i_ = Num_i_; i_ < Num_; i_++)
                {
                    dx = x_[i_];
                    dy = y_[i_];
                    du = U_[i_];
                    dv = V_[i_];
            
                    j_ = i_ - Num_i_;
            
                    Mat[j_][0] = w_[i_];      // w
                    Mat[j_][1] = w_[i_] * dx; // x
                    Mat[j_][2] = w_[i_] * dy; // y
            
                    Mat[j_ + N_][0] = 0.0; // w
                    Mat[j_ + N_][1] = -w_[i_] * (y0 + dy) / 2.0;
                    Mat[j_ + N_][2] = 0.0; // y
            
                    Vec_[j_] = w_[i_] * du;      // x
                    Vec_[j_ + N_] = w_[i_] * dv; // x
                }
                double TwoPsixPsiy, Psixsq_minus_Psiysq;
                for (i_ = 0; i_ < Num_i_; i_++)
                {
                    dx = x_[i_];
                    dy = y_[i_];
                    j_ = 2 * N_ + i_;
            
                    TwoPsixPsiy = 2.0 * idata.psix_[i_] * idata.psiy_[i_];
                    Psixsq_minus_Psiysq = (idata.psix_[i_] * idata.psix_[i_] - idata.psiy_[i_] * idata.psiy_[i_]);
            
                    Mat[j_][0] = 0.0;                      // w
                    Mat[j_][1] = -3.0 * TwoPsixPsiy / 2.0; // x
                    Mat[j_][2] = Psixsq_minus_Psiysq;      // y
            
                    Constr_[i_] = 0.0;
                }
                QRdcmp<double> QR(Mat, 2 * Num_ - Num_i_, Size_Loc);
                QR.get_Q(Q);
                QR.get_R(R);
            
                for (j_ = 0; j_ < Size_Loc; j_++)
                {
                    for (i_ = 0; i_ < 2 * N_; i_++)
                        Q1[i_][j_] = Q[i_][j_];
                    for (i_ = 0; i_ < Num_i_; i_++)
                        Q2T[j_][i_] = Q[2 * N_ + i_][j_];
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
                    for (j_ = 0; j_ < 2 * N_; j_++)
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
                for (i_ = 0; i_ < Size_Loc; i_++)
                {
                    Vec_[i_] = Q1Tb[i_];
                    for (j_ = 0; j_ < Num_i_; j_++)
                        Vec_[i_] -= 0.5 * Q2T[i_][j_] * Sol_w_[j_];
                }
                // solve (Mat_R) w = RHS from above
                for (i_ = Size_Loc - 1; i_ >= 0; i_--)
                {
                    Sol[i_] = Vec_[i_] / R[i_][i_];
                    for (j_ = i_ + 1; j_ < Size_Loc; j_++)
                        Sol[i_] -= Sol[j_] * R[i_][j_] / R[i_][i_];
                }
            
                U(idata.cellid_, 0) = Sol[0];
                U(idata.cellid_, 1) = -Sol[1] * y0 / 2.0;
                double ux, uy, vx, vy;
                for (i_ = 0; i_ < Num_i_; i_++)
                {
            
                    dx = x_[i_];
                    dy = y_[i_];
            
                    TwoPsixPsiy = 2.0 * idata.psix_[i_] * idata.psiy_[i_];
                    Psixsq_minus_Psiysq = (idata.psix_[i_] * idata.psix_[i_] - idata.psiy_[i_] * idata.psiy_[i_]);
            
                    ux = Sol[1];
                    uy = Sol[2];
                    vx = 0.0;
                    vy = -Sol[1] / 2.0;
            
                    ux /= deltax[0];
                    uy /= deltax[1];
                    vx /= deltax[0];
                    vy /= deltax[1];
            
                    //TODO  variable viscosity
                    dGamma_dt = visc.ComputeGammaDot(ux, 0.5 * (uy + vx), vy);
                    MU = visc.GetViscosity(dGamma_dt);
                    // MU = MU_0;
            
                    idata.Gamma_dot_Int[i_] = dGamma_dt;
                    idata.tan_shear_[i_] = Psixsq_minus_Psiysq * (uy + vx) + TwoPsixPsiy * (vy - ux);
                    idata.norm_shear_[i_] = 2.0 * ux * idata.psix_[i_] * idata.psix_[i_] + 2.0 * vy * idata.psiy_[i_] * idata.psiy_[i_] + 2.0 * (vx + uy) * idata.psix_[i_] * idata.psiy_[i_];
                    idata.norm_shear_[i_] /= (idata.psix_[i_] * idata.psix_[i_] + idata.psiy_[i_] * idata.psiy_[i_] + 1.0E-12);
            
                    //idata.norm_shear_[i_] *= MU;
                }
            
                for (i_ = 0; i_ < 2 * Num_ - Num_i_; i_++)
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
                for (i_ = 0; i_ < 2 * N_; i_++)
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
            }
        
        } // namespace Axisymmetric
        
        
    } // namespace AdvectLSIntVel
/**********************************************************************/
    namespace Axisymmetric
    {
    
    void LS_Vel_Fourth_Order
    (
        bool &sing,
        int Num_,
        const amrex::IntVect &cell,
        amrex::Array4<amrex::Real> const &U,
        double *x_,
        double *y_,
        double *U_,
        double *V_,
        double *w_,
        const amrex::Real* deltax,
        const amrex::Real* prob_lo
    )
    {
        int i_, j_;
        double **Mat, **R; // Least squares and constraint matrices
        double *Sol, *Vec_;
        double dx, dy, du, dv, nx, ny, y0;
    
        Allocate_2D_R(Mat, 2 * Num_, 10);
	Allocate_2D_R(R, 10, 10);
        Sol = new double[10];
        Vec_ = new double[2 * Num_];
    
        amrex::Real y = prob_lo[1] + deltax[1] * (cell[1] + 0.5);
        y0 = y / deltax[1];
    
        for (i_ = 0; i_ < Num_; i_++)
        {
            dx = x_[i_];
            dy = y_[i_];
            du = U_[i_];
            dv = V_[i_];
    
            j_ = i_;
            Mat[j_][0] = w_[i_];           // w
            Mat[j_][1] = w_[i_] * dx;      // x
            Mat[j_][2] = w_[i_] * dy;      // y
            Mat[j_][3] = w_[i_] * dx * dx; // x^2
            Mat[j_][4] = w_[i_] * dx * dy; // xy
            Mat[j_][5] = w_[i_] * dy * dy; // y^2
    
            Mat[j_][6] = w_[i_] * dx * dx * dx; // v_0
            Mat[j_][7] = w_[i_] * dx * dx * dy; // v_1
            Mat[j_][8] = w_[i_] * dx * dy * dy; // v_11
            Mat[j_][9] = w_[i_] * dy * dy * dy; // v_11
    
            Mat[j_ + Num_][0] = 0.0;
            Mat[j_ + Num_][1] = -w_[i_] * (y0 + dy) / 2.0;
            Mat[j_ + Num_][2] = 0.0;
            Mat[j_ + Num_][3] = -w_[i_] * dx * (y0 + dy);
            Mat[j_ + Num_][4] = w_[i_] * (y0 + dy) * (y0 - 2.0 * dy) / 6.0;
            Mat[j_ + Num_][5] = 0.0;
    
            Mat[j_ + Num_][6] = -3.0 * w_[i_] * dx * dx * (y0 + dy) / 2.0;
            Mat[j_ + Num_][7] = w_[i_] * dx * (y0 + dy) * (y0 - 2.0 * dy) / 3.0;
            Mat[j_ + Num_][8] = -w_[i_] * (y0 + dy) * (y0 * y0 - 2.0 * y0 * dy + 3.0 * dy * dy) / 12.0;
            Mat[j_ + Num_][9] = 0.0;
    
            Vec_[j_] = w_[i_] * du;        // x
            Vec_[j_ + Num_] = w_[i_] * dv; // x
        }
        QRdcmp<double> QR(Mat, 2 * Num_, 10);
        QR.solve(Vec_, Sol);
	QR.get_R(R);
    
        U(cell, 0) = Sol[0];
        U(cell, 1) = -0.5 * y0 * Sol[1] + y0 * y0 * Sol[4] / 6.0 - y0 * y0 * y0 * Sol[8] / 12.0;

        sing = false;
        for(i_ = 0 ; i_ < 10; i_++)
        {
            if(std::abs(R[i_][i_]) < 1.0e-5) sing = true;
        }
    
        for (i_ = 0; i_ < 2 * Num_; i_++)
        {
            delete[] Mat[i_];
        }
        for (i_ = 0; i_ < 10; i_++)
        {
            delete[] R[i_];
        }

	delete[] R;
        delete[] Mat;
        delete[] Sol;
        delete[] Vec_;
    }
    
    void LS_Vel_Third_Order
    (
        bool &sing,
        int Num_,
        const amrex::IntVect &cell,
        amrex::Array4<amrex::Real> const &U,
        double *x_,
        double *y_,
        double *U_,
        double *V_,
        double *w_,
        const amrex::Real* deltax,
        const amrex::Real* prob_lo
    )
    {
        int i_, j_;
        double **Mat, **R; // Least squares and constraint matrices
        double *Sol, *Vec_;
        double dx, dy, du, dv, nx, ny, y0;
    
        Allocate_2D_R(Mat, 2 * Num_, 6);
	Allocate_2D_R(R, 6, 6);
        Sol = new double[6];
        Vec_ = new double[2 * Num_];
    
        amrex::Real y = prob_lo[1] + deltax[1] * (cell[1] + 0.5);
        y0 = y / deltax[1];
    
        for (i_ = 0; i_ < Num_; i_++)
        {
            dx = x_[i_];
            dy = y_[i_];
            du = U_[i_];
            dv = V_[i_];
    
            j_ = i_;
            Mat[j_][0] = w_[i_];           // w
            Mat[j_][1] = w_[i_] * dx;      // x
            Mat[j_][2] = w_[i_] * dy;      // y
            Mat[j_][3] = w_[i_] * dx * dx; // x^2
            Mat[j_][4] = w_[i_] * dx * dy; // xy
            Mat[j_][5] = w_[i_] * dy * dy; // y^2
    
            Mat[j_ + Num_][0] = 0.0;
            Mat[j_ + Num_][1] = -w_[i_] * (y0 + dy) / 2.0;
            Mat[j_ + Num_][2] = 0.0;
            Mat[j_ + Num_][3] = -w_[i_] * dx * (y0 + dy);
            Mat[j_ + Num_][4] = w_[i_] * (y0 + dy) * (y0 - 2.0 * dy) / 6.0;
            Mat[j_ + Num_][5] = 0.0;
    
            Vec_[j_] = w_[i_] * du;        // x
            Vec_[j_ + Num_] = w_[i_] * dv; // x
        }
        QRdcmp<double> QR(Mat, 2 * Num_, 6);
        QR.solve(Vec_, Sol);
	QR.get_R(R);
    
        U(cell, 0) = Sol[0];
        U(cell, 1) = -0.5 * y0 * Sol[1] + y0 * y0 * Sol[4] / 6.0;
    
        sing = false;
        for(i_ = 0 ; i_ < 6; i_++)
        {
            if(std::abs(R[i_][i_]) < 1.0e-5) sing = true;
        }

        for (i_ = 0; i_ < 2 * Num_; i_++)
        {
            delete[] Mat[i_];
        }
        for (i_ = 0; i_ < 6; i_++)
        {
            delete[] R[i_];
        }
        delete[] R;
        delete[] Mat;
        delete[] Sol;
        delete[] Vec_;
    }
    
    void LS_Vel_Second_Order
    (
        bool &sing,
        int Num_,
        const amrex::IntVect &cell,
        amrex::Array4<amrex::Real> const &U,
        double *x_,
        double *y_,
        double *U_,
        double *V_,
        double *w_,
        const amrex::Real* deltax,
        const amrex::Real* prob_lo
    )
    {
        int i_, j_;
        double **Mat; // Least squares and constraint matrices
        double *Sol, *Vec_;
        double dx, dy, du, dv, nx, ny, y0;
    
        Allocate_2D_R(Mat, 2 * Num_, 3);
        Sol = new double[3];
        Vec_ = new double[2 * Num_];
    
        amrex::Real y = prob_lo[1] + deltax[1] * (cell[1] + 0.5);
        y0 = y / deltax[1];
    
        for (i_ = 0; i_ < Num_; i_++)
        {
            dx = x_[i_];
            dy = y_[i_];
            du = U_[i_];
            dv = V_[i_];
    
            j_ = i_;
            Mat[j_][0] = w_[i_];      // w
            Mat[j_][1] = w_[i_] * dx; // x
            Mat[j_][2] = w_[i_] * dy; // y
    
            Mat[j_ + Num_][0] = 0.0;
            Mat[j_ + Num_][1] = -w_[i_] * (y0 + dy) / 2.0;
            Mat[j_ + Num_][2] = 0.0;
    
            Vec_[j_] = w_[i_] * du;        // x
            Vec_[j_ + Num_] = w_[i_] * dv; // x
        }
        QRdcmp<double> QR(Mat, 2 * Num_, 3);
        QR.solve(Vec_, Sol);
    
        U(cell, 0) = Sol[0];
        U(cell, 1) = -0.5 * y0 * Sol[1];
    
        for (i_ = 0; i_ < 2 * Num_; i_++)
        {
            delete[] Mat[i_];
        }
        delete[] Mat;
        delete[] Sol;
        delete[] Vec_;
    }
    
    } // namespace Axisymmetric

/**********************************************************************/




    
    int Mask::Cubic_Solve(double &x_, double a_, double b_, double c_, double d_)
    {
        int iter_cu = 0;
        double a, b, c, d, Res, den, x_o;
    
        x_o = x_;
        if (fabs(x_o) > 1.0)
        {
            amrex::Print() << "\n"
                           << "Second-order estimate out of bounds! "
                           << "\n";
            exit(1);
        }
        // Note the negative signs compared to Nourgialev, .. h * theta = x_i - xI
        a = -(3.0 * (b_ - c_) - a_ + d_) / 6.0;
        b = (c_ + a_ - 2.0 * b_) / 2.0;
        c = (3.0 * b_ + 2.0 * a_ - 6.0 * c_ + d_) / 6.0;
        d = b_;
    
        double x_L, x_U, F_L, F_U, x_m, F_m;
        x_L = 0.0;
        x_U = Sign(x_);
        x_m = x_;
        bool STOP = false;
        F_L = a * x_L * x_L * x_L + b * x_L * x_L + c * x_L + d;
        F_U = a * x_U * x_U * x_U + b * x_U * x_U + c * x_U + d;
        F_m = a * x_m * x_m * x_m + b * x_m * x_m + c * x_m + d;
    
        if (F_L * F_m < 0.0)
        {
            x_U = x_m;
            F_U = F_m;
        }
        else if (F_U * F_m < 0.0)
        {
            x_L = x_m;
            F_L = F_m;
        }
        else
        {
            //amrex::Print() << "\n Found the root without iterations!"
            //               << "\n";
            if (fabs(F_L) < 1.0E-15)
            {
                x_ = x_L;
                STOP = true;
            }
            else if (fabs(F_U) < 1.0E-15)
            {
                x_ = x_U;
                STOP = true;
            }
            else if (fabs(F_m) < 1.0E-15)
            {
                x_ = x_m;
                STOP = true;
            }
        }
    
        do
        {
    
            if (fabs(F_L) < 1.0E-15)
            {
                x_ = x_L;
                STOP = true;
            }
            else if (fabs(F_U) < 1.0E-15)
            {
                x_ = x_U;
                STOP = true;
            }
            else
            {
                x_m = 0.5 * (x_L + x_U);
                F_m = a * x_m * x_m * x_m + b * x_m * x_m + c * x_m + d;
                if (fabs(F_m) < 1.0E-15)
                {
                    x_ = x_m;
                    STOP = true;
                }
                if (I_DSign(F_m) == I_DSign(F_L))
                {
                    x_L = x_m;
                    F_L = F_m;
                }
                if (I_DSign(F_m) == I_DSign(F_U))
                {
                    x_U = x_m;
                    F_U = F_m;
                }
            }
            iter_cu++;
            //        amrex::Print() << iter_cu << "\t" << F_m << "\t" << x_m <<  "\n" ;
        } while (iter_cu < 10 && !STOP);
        x_ = 0.5 * (x_L + x_U);
    
        do
        {
            den = 3.0 * a * x_ * x_ + 2.0 * b * x_ + c;
            if (fabs(den) < 1.0E-12)
            {
                return -1;
            }
            Res = (a * x_ * x_ * x_ + b * x_ * x_ + c * x_ + d) / den;
            x_ -= Res;
            iter_cu++;
        } while (iter_cu < 21 && fabs(Res) > 1.0E-12);
    
        if (fabs(x_) > 1.0)
        {
            amrex::Print() << "\n"
                           << "Cubic solve out of bounds!"
                           << "\n";
            amrex::Print() << "\n"
                           << "Number of Newton's iterations: " << iter_cu << "\n";
            amrex::Print() << x_o << "\t" << x_ << "\n";
            amrex::Print() << a_ << "\t" << b_ << "\t" << c_ << "\t" << d_ << "\n";
            amrex::Print() << "Cubic: " << a << "\t" << b << "\t" << c << "\t" << d << "\n";
            exit(1);
        }
    
        return iter_cu;
    }
    
    void Mask::Compute_Normal_Curvature
    (
        amrex::Real *dPsi,
        amrex::Real &Curv,
        int i, int j,
        amrex::Array4<amrex::Real const> const &Psi,
        const amrex::Real &deltax
    )
    {
        amrex::Real Psi_xx, Psi_xy, Psi_yy;
        int k = 0;
        dPsi[0] = 0.5 * (Psi(i + 1, j, k) - Psi(i - 1, j, k));
        dPsi[1] = 0.5 * (Psi(i, j + 1, k) - Psi(i, j - 1, k));
        Psi_xx = Psi(i + 1, j, k) - 2.0 * Psi(i, j, k) + Psi(i - 1, j, k);
        Psi_yy = Psi(i, j + 1, k) - 2.0 * Psi(i, j, k) + Psi(i, j - 1, k);
        Psi_xy = 0.25 * (Psi(i + 1, j + 1, k) - Psi(i + 1, j - 1, k) - Psi(i - 1, j + 1, k) + Psi(i - 1, j - 1, k));
    
        Curv = (Psi_xx * (dPsi[1] * dPsi[1]) + Psi_yy * (dPsi[0] * dPsi[0]) - 2.0 * dPsi[0] * dPsi[1] * Psi_xy) 
             / (deltax * pow(dPsi[0] * dPsi[0] + dPsi[1] * dPsi[1] + 1.0E-14, 1.5));
    }
    void Mask::Compute_Normal_Curvature
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
    
    void Mask::Compute_LSQ_Weights(double &weights_, double dist_)
    {
        weights_ = 1.0 / pow((dist_ * dist_ + 1.0), exp_);
    }
    
    void Mask::LSQ_Parameters
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
    
        P_ = P(i1, j1, 0);
        Compute_LSQ_Weights(w_, std::sqrt(x_ * x_ + y_ * y_));
    }

    void Mask::QR_LS_Pressure_Second_Order
    (
        bool& sing, 
        int Num_, 
        int Num_i_, 
        double* x_, double* y_, 
        double* P_, 
        double* w_,
        amrex::Array4<amrex::Real> const& P,
        InterceptData& idata
    ) 
    {
        int i_, j_ ; 
        double **Mat, *Sol, *Vec_  ;
    
        Allocate_2D_R(Mat,Num_,3) ; Sol = new double[3] ; Vec_ = new double[Num_] ;
    
        for(i_ = 0 ; i_ < Num_ ; i_++) 
        {
    	    Mat[i_][0] = w_[i_] ; // 1
    	    Mat[i_][1] = w_[i_]*x_[i_] ; // x
    	    Mat[i_][2] = w_[i_]*y_[i_] ; // y
    	    Vec_[i_] = w_[i_]*P_[i_] ; // P
        } 
        QRdcmp<double> QR(Mat, Num_, 3); 
        QR.solve(Vec_, Sol);
    
        /// P(x,y) = a + bx + cy
        /// @ cut cell the (x,y) = (0,0) as all x_, y_ are computed wrt cut cell
        /// so P = a @ cut cell
        P(idata.cellid_) = Sol[0];
    
        /// compute P @ intercepts
    	for(i_ = 0 ; i_ < Num_i_ ; i_++)
            idata.P_comp[i_] = Sol[0] + Sol[1] * x_[i_] + Sol[2] * y_[i_];
        
        /// 
        for (i_ = 0; i_ < Num_; i_++)
        {
            for(j_ = 0 ; j_ < Num_ ; j_++) Vec_[j_] = 0.0 ;
    	    Vec_[i_] = w_[i_] ;
    	    QR.solve(Vec_, Sol);
            idata.Interpolation_weights[i_] = Sol[0];
        }
    
        for(i_ = 0 ; i_ < Num_ ;i_++) { delete [] Mat[i_] ; } 
        delete [] Mat ; 
        delete [] Sol ; delete [] Vec_ ;
        sing = false;
    }
    
    void Mask::QR_LS_Pressure_Third_Order
    (
        bool& sing, 
        int Num_, 
        int Num_i_, 
        double* x_, 
        double* y_, 
        double* P_, 
        double* w_,
        amrex::Array4<amrex::Real> const& P,
        InterceptData& idata
    ) 
    {
        int i_, j_ ;
        double **Mat, *Sol, *Vec_  ;
    
        Allocate_2D_R(Mat,Num_,6) ; 
        Sol = new double[6] ; 
        Vec_ = new double[Num_] ;
    
        for(i_ = 0 ; i_ < Num_ ; i_++) 
        {
        	Mat[i_][0] = w_[i_] ; // 1
            Mat[i_][1] = w_[i_]*x_[i_] ; // x
    	    Mat[i_][2] = w_[i_]*y_[i_] ; // y
    	    Mat[i_][3] = Mat[i_][1]*x_[i_] ; // x^2
    	    Mat[i_][4] = Mat[i_][1]*y_[i_] ; // xy
    	    Mat[i_][5] = Mat[i_][2]*y_[i_] ; // y^2
    	    Vec_[i_] = w_[i_]*P_[i_] ; // x
        } 
        QRdcmp<double> QR(Mat, Num_, 6); 
        QR.solve(Vec_, Sol);
        P(idata.cellid_) = Sol[0];
        for(i_ = 0 ; i_ < Num_i_ ; i_++) idata.P_comp[i_] = Sol[0] + Sol[1]*x_[i_] + Sol[2]*y_[i_] + Sol[3]*x_[i_]*x_[i_] + Sol[4]*x_[i_]*y_[i_] + Sol[5]*y_[i_]*y_[i_] ;
        for(i_ = 0 ; i_ < Num_ ; i_++) 
        {
    	    for(j_ = 0 ; j_ < Num_ ; j_++) Vec_[j_] = 0.0 ;
    	    Vec_[i_] = w_[i_] ;
    	    QR.solve(Vec_, Sol);
            idata.Interpolation_weights[i_] = Sol[0];
        }
    
        for(i_ = 0 ; i_ < Num_ ;i_++) { delete [] Mat[i_] ; } 
        delete [] Mat ; 
        delete [] Sol ; 
        delete [] Vec_ ;
        sing = false;
    }
    
    void Mask::AllMCQR_LS_Pressure_Fourth_Order
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
    )
    {
        int i_, j_, k_;
        double **Mat;    // Least squares and constraint matrices
        double **Q, **R; // Q R decomosition of Mat;
        double **Q1, **Q2T, **Q2T_Q, **Q2T_R;
        double *Sol, *Vec_, *Constr_, *Sol_u_, *Sol_w_, *Q1Tb; // *Q2T_QTQ1Tb ;
        double dx, dy, dP, nx, ny;
    
        Allocate_2D_R(Mat, Num_, 10);
        Sol = new double[10];
        Vec_ = new double[Num_ - Num_i_];
        Allocate_2D_R(Q, Num_, 10);
        Allocate_2D_R(R, 10, 10);
        Allocate_2D_R(Q1, Num_ - Num_i_, 10);
        Allocate_2D_R(Q2T, 10, Num_i_);
        Allocate_2D_R(Q2T_Q, 10, Num_i_);
        Allocate_2D_R(Q2T_R, Num_i_, Num_i_);
        Constr_ = new double[Num_i_];
        Sol_u_ = new double[Num_i_];
        Sol_w_ = new double[Num_i_];
        Q1Tb = new double[10]; // Q2T_QTQ1Tb = new double[Num_i_] ;
    
        auto &PsiX = idata.psix_;
        auto &PsiY = idata.psiy_;
    
        for (i_ = Num_i_; i_ < Num_; i_++)
        {
            dx = x_[i_];
            dy = y_[i_];
            dP = P_[i_];
    
            j_ = i_ - Num_i_;
            Mat[j_][0] = w_[i_];          // w
            Mat[j_][1] = w_[i_] * dx;     // x
            Mat[j_][2] = w_[i_] * dy;     // y
            Mat[j_][3] = Mat[j_][1] * dx; // x^2
            Mat[j_][4] = Mat[j_][1] * dy; // xy
            Mat[j_][5] = Mat[j_][2] * dy; // y^2
            Mat[j_][6] = Mat[j_][3] * dx; // x^3
            Mat[j_][7] = Mat[j_][4] * dx; // x^2y
            Mat[j_][8] = Mat[j_][5] * dx; // xy^2
            Mat[j_][9] = Mat[j_][5] * dy; // y^3
            Vec_[j_] = w_[i_] * dP;       // x
        }
    
        if (bc == amrex::LinOpBCType::Dirichlet)
        {
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                dP = P_[i_];
                j_ = Num_ - Num_i_ + i_;
                Mat[j_][0] = 1.0;             // 1
                Mat[j_][1] = dx;              // x
                Mat[j_][2] = dy;              // y
                Mat[j_][3] = Mat[j_][1] * dx; // x^2
                Mat[j_][4] = Mat[j_][1] * dy; // xy
                Mat[j_][5] = Mat[j_][2] * dy; // y^2
                Mat[j_][6] = Mat[j_][3] * dx; // x^3
                Mat[j_][7] = Mat[j_][4] * dx; // x^2y
                Mat[j_][8] = Mat[j_][5] * dx; // xy^2
                Mat[j_][9] = Mat[j_][5] * dy; // y^3
                Constr_[i_] = dP;
            }
        }
        else if (bc == amrex::LinOpBCType::Neumann)
        {
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                dP = 0.0;
                j_ = Num_ - Num_i_ + i_;
    
                nx = PsiX[i_] / sqrt(PsiX[i_] * PsiX[i_] + PsiY[i_] * PsiY[i_]);
                ny = sqrt(1.0 - ny * ny);
                nx = PsiX[i_];
                ny = PsiY[i_];
                Mat[j_][0] = 0.0;                               // 1
                Mat[j_][1] = nx;                                // x
                Mat[j_][2] = ny;                                // y
                Mat[j_][3] = 2.0 * nx * dx;                     // x^2
                Mat[j_][4] = (nx * dy + ny * dx);               // xy
                Mat[j_][5] = 2.0 * ny * dy;                     // y^2
                Mat[j_][6] = 3.0 * nx * dx * dx;                // x^3
                Mat[j_][7] = 2.0 * nx * dx * dy + ny * dx * dx; // x^2y
                Mat[j_][8] = nx * dy * dy + 2.0 * ny * dx * dy; // xy^2
                Mat[j_][9] = 3.0 * ny * dy * dy;                // y^3
                Constr_[i_] = dP;
            }
        }
    
        QRdcmp<double> QR(Mat, Num_, 10);
        QR.get_Q(Q);
        QR.get_R(R);
    
        for (j_ = 0; j_ < 10; j_++)
        {
            for (i_ = 0; i_ < Num_ - Num_i_; i_++)
                Q1[i_][j_] = Q[i_][j_];
            for (i_ = 0; i_ < Num_i_; i_++)
                Q2T[j_][i_] = Q[Num_ - Num_i_ + i_][j_];
        }
        QRdcmp<double> QR1(Q2T, 10, Num_i_);
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
        for (i_ = 0; i_ < 10; i_++)
        {
            Q1Tb[i_] = 0.0;
            for (j_ = 0; j_ < Num_ - Num_i_; j_++)
                Q1Tb[i_] += Q1[j_][i_] * Vec_[j_];
        }
    
        // Find the vector 2.0*(Q2T_Q)T Q1Tb -2.0*Sol_u_;
        for (i_ = 0; i_ < Num_i_; i_++)
        {
            Vec_[i_] = -2.0 * Sol_u_[i_];
            for (j_ = 0; j_ < 10; j_++)
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
        for (i_ = 0; i_ < 10; i_++)
        {
            Vec_[i_] = Q1Tb[i_];
            for (j_ = 0; j_ < Num_i_; j_++)
                Vec_[i_] -= 0.5 * Q2T[i_][j_] * Sol_w_[j_];
        }
        // solve (Mat_R) w = RHS from above
        for (i_ = 9; i_ >= 0; i_--)
        {
            Sol[i_] = Vec_[i_] / R[i_][i_];
            for (j_ = i_ + 1; j_ < 10; j_++)
                Sol[i_] -= Sol[j_] * R[i_][j_] / R[i_][i_];
        }
    
        P(idata.cellid_) = Sol[0];
    
        if (bc == amrex::LinOpBCType::Dirichlet)
        {
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                idata.P_comp[i_] = Sol[0] + Sol[1] * dx + Sol[2] * dy + Sol[3] * dx * dx + Sol[4] * dx * dy + Sol[5] * dy * dy + Sol[6] * dx * dx * dx + Sol[7] * dx * dx * dy + Sol[8] * dx * dy * dy + Sol[9] * dy * dy * dy;
            }
        }
    
        else if (bc == amrex::LinOpBCType::Neumann)
        {
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                dx = x_[i_];
                dy = y_[i_];
                nx = PsiX[i_] / sqrt(PsiX[i_] * PsiX[i_] + PsiY[i_] * PsiY[i_]);
                ny = sqrt(1.0 - ny * ny);
                nx = PsiX[i_];
                ny = PsiY[i_];
                idata.P_comp[i_] = Sol[1] * nx + Sol[2] * ny + 2.0 * Sol[3] * nx * dx + Sol[4] * (nx * dy + ny * dx) + 2.0 * Sol[5] * ny * dy + 3.0 * Sol[6] * nx * dx * dx + Sol[7] * (2.0 * nx * dx * dy + ny * dx * dx) + Sol[8] * (2.0 * dx * dy * ny + dy * dy * nx) + 3.0 * Sol[9] * ny * dy * dy;
                idata.P_Int[i_] = Sol[0] + Sol[1] * dx + Sol[2] * dy + Sol[3] * dx * dx + Sol[4] * dx * dy + Sol[5] * dy * dy + Sol[6] * dx * dx * dx + Sol[7] * dx * dx * dy + Sol[8] * dx * dy * dy + Sol[9] * dy * dy * dy;
            }
        }
        // For the interface... Here the vector b = 0.
        for (k_ = 0; k_ < Num_i_; k_++)
        {
            for (j_ = 0; j_ < Num_i_; j_++)
                Constr_[j_] = 0.0;
            Constr_[k_] = 1.0;
    
            // solve for (Q2T_R)T u = Constr ;
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                Sol_u_[i_] = Constr_[i_] / Q2T_R[i_][i_];
                for (j_ = 0; j_ < i_; j_++)
                    Sol_u_[i_] -= Sol_u_[j_] * Q2T_R[j_][i_] / Q2T_R[i_][i_];
            }
            // Find the vector -2.0*Sol_u_;
            for (i_ = 0; i_ < Num_i_; i_++)
                Vec_[i_] = -2.0 * Sol_u_[i_];
            // solve (Q2T_R) w = RHS from above
            for (i_ = Num_i_ - 1; i_ >= 0; i_--)
            {
                Sol_w_[i_] = Vec_[i_] / Q2T_R[i_][i_];
                for (j_ = i_ + 1; j_ < Num_i_; j_++)
                    Sol_w_[i_] -= Sol_w_[j_] * Q2T_R[i_][j_] / Q2T_R[i_][i_];
            }
            // Find the vector - Q2Tw/2;
            for (i_ = 0; i_ < 10; i_++)
            {
                Vec_[i_] = 0.0;
                for (j_ = 0; j_ < Num_i_; j_++)
                    Vec_[i_] -= 0.5 * Q2T[i_][j_] * Sol_w_[j_];
            }
            // solve (Mat_R) w = RHS from above
            for (i_ = 9; i_ >= 0; i_--)
            {
                Sol[i_] = Vec_[i_] / R[i_][i_];
                for (j_ = i_ + 1; j_ < 10; j_++)
                    Sol[i_] -= Sol[j_] * R[i_][j_] / R[i_][i_];
            }
            idata.Interpolation_weights[k_] = Sol[0];
        }
        for (k_ = Num_i_; k_ < Num_; k_++)
        {
            for (j_ = Num_i_; j_ < Num_; j_++)
                Vec_[j_ - Num_i_] = 0.0;
            Vec_[k_ - Num_i_] = w_[k_];
    
            // Find the vector Q1T b.
            for (i_ = 0; i_ < 10; i_++)
            {
                Q1Tb[i_] = 0.0;
                for (j_ = 0; j_ < Num_ - Num_i_; j_++)
                    Q1Tb[i_] += Q1[j_][i_] * Vec_[j_];
            }
    
            // Find the vector 2.0*(Q2T_Q)T Q1Tb -2.0*Sol_u_;
            for (i_ = 0; i_ < Num_i_; i_++)
            {
                Vec_[i_] = 0.0;
                for (j_ = 0; j_ < 10; j_++)
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
            for (i_ = 0; i_ < 10; i_++)
            {
                Vec_[i_] = Q1Tb[i_];
                for (j_ = 0; j_ < Num_i_; j_++)
                    Vec_[i_] -= 0.5 * Q2T[i_][j_] * Sol_w_[j_];
            }
            // solve (Mat_R) w = RHS from above
            for (i_ = 9; i_ >= 0; i_--)
            {
                Sol[i_] = Vec_[i_] / R[i_][i_];
                for (j_ = i_ + 1; j_ < 10; j_++)
                    Sol[i_] -= Sol[j_] * R[i_][j_] / R[i_][i_];
            }
            idata.Interpolation_weights[k_] = Sol[0];
        }

        sing = false;
        for(i_ = 0 ; i_ < 10; i_++)
        {
            if(std::abs(R[i_][i_]) < 1.0e-6) sing = true;
        }
    
        for (i_ = 0; i_ < Num_; i_++)
        {
            delete[] Mat[i_];
            delete[] Q[i_];
        }
        for (i_ = 0; i_ < 10; i_++)
        {
            delete[] R[i_];
            delete[] Q2T[i_];
            delete[] Q2T_Q[i_];
        }
        for (i_ = 0; i_ < Num_ - Num_i_; i_++)
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
    }
   
    void Mask::QR_LS_Temperature_Second_Order
    (
        bool& sing,
        int Num_, 
        int Num_i_, 
        double* x_, double* y_,
        double* T_,
        double* w_,
        amrex::Array4<amrex::Real> const& T,
        InterceptData& idata
    )
    {
        int i_, j_ ;  
        double **Mat, *Sol, *Vec_  ;
    
        Allocate_2D_R(Mat,Num_,3) ; Sol = new double[3] ; Vec_ = new double[Num_] ;
    
        for(i_ = 0 ; i_ < Num_ ; i_++)
        {
            Mat[i_][0] = w_[i_] ; // 1 
            Mat[i_][1] = w_[i_]*x_[i_] ; // x
            Mat[i_][2] = w_[i_]*y_[i_] ; // y
            Vec_[i_] = w_[i_]*T_[i_] ; // P
        } 
        QRdcmp<double> QR(Mat, Num_, 3);
        QR.solve(Vec_, Sol);
    
        /// P(x,y) = a + bx + cy
        /// @ cut cell the (x,y) = (0,0) as all x_, y_ are computed wrt cut cell
        /// so P = a @ cut cell
        T(idata.cellid_) = Sol[0];
    
        for(i_ = 0 ; i_ < Num_ ;i_++) { delete [] Mat[i_] ; }
            delete [] Mat ; 
            delete [] Sol ; delete [] Vec_ ;
    }
    
    void Mask::QR_LS_Temperature_Third_Order
    (
        bool& sing,
        int Num_,
        int Num_i_,
        double* x_,
        double* y_,
        double* T_,
        double* w_,
        amrex::Array4<amrex::Real> const& T,
        InterceptData& idata
    )
    {
            int i_, j_ ;
            double **Mat, *Sol, *Vec_  ;
    
            Allocate_2D_R(Mat,Num_,6) ; Sol = new double[6] ; Vec_ = new double[Num_] ;
    
            for(i_ = 0 ; i_ < Num_ ; i_++) {
                    Mat[i_][0] = w_[i_] ; // 1
                    Mat[i_][1] = w_[i_]*x_[i_] ; // x
                    Mat[i_][2] = w_[i_]*y_[i_] ; // y
                    Mat[i_][3] = Mat[i_][1]*x_[i_] ; // x^2
                    Mat[i_][4] = Mat[i_][1]*y_[i_] ; // xy
                    Mat[i_][5] = Mat[i_][2]*y_[i_] ; // y^2
                    Vec_[i_] = w_[i_]*T_[i_] ; // x
            }
            QRdcmp<double> QR(Mat, Num_, 6);
            QR.solve(Vec_, Sol);
            T(idata.cellid_) = Sol[0];
    
            for(i_ = 0 ; i_ < Num_ ;i_++) { delete [] Mat[i_] ; }
            delete [] Mat ;
            delete [] Sol ; delete [] Vec_ ;
    }
    
        
    
    
    void Mask::MaskBC()
    {
        mask_.FillBoundary();
        maskTemp_.FillBoundary();
        const amrex::Box &domain = geom_.Domain();
        int ngrow = mask_.nGrow();
        for (amrex::MFIter mfi(mask_); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
    
            amrex::Array4<int> const &pmask = mask_.array(mfi);
            amrex::Array4<int> const &pmask_temp = maskTemp_.array(mfi);
    
            int k = 0; // as 2d
    
            amrex::Box gbx(bx);
            gbx.grow(ngrow);
    
            /// left
            if (bx.smallEnd(0) == domain.smallEnd(0))
            {
                //int ig = bx.smallEnd(0) - ngrow;
                //int ip = bx.smallEnd(0) + ngrow - 1;
    
                int ig = bx.smallEnd(0) - 1;
                for (int i = bx.smallEnd(0); i < bx.smallEnd(0) + ngrow; i++)
                {
                    for (int j = bx.smallEnd(1) - ngrow; j <= bx.bigEnd(1) + ngrow; j++)
                    {
                        //pmask(i-ngrow, j, k) = pmask(2 * ngrow - 1 - i, j, k);
                        pmask(ig, j, k) = pmask(i, j, k);
                        pmask_temp(ig, j, k) = pmask_temp(i, j, k);
                    }
                    ig--;
                }
            }
            /// bottom
            if (bx.smallEnd(1) == domain.smallEnd(1))
            {
                //int jm = bx.smallEnd(1) - ngrow;
                //int jp = bx.smallEnd(1) + ngrow - 1;
    
                int jg = bx.smallEnd(1) - 1;
                for (int j = bx.smallEnd(1); j < bx.smallEnd(1) + ngrow; j++)
                {
                    for (int i = bx.smallEnd(0) - ngrow; i <= bx.bigEnd(0) + ngrow; i++)
                    {
                        pmask(i, jg, k) = pmask(i, j, k);
                        pmask_temp(i, jg, k) = pmask_temp(i, j, k);
                    }
                    jg--;
                }
            }
            /// right
            if (bx.bigEnd(0) == domain.bigEnd(0))
            {
                int ig = bx.bigEnd(0) + 1;
                for (int i = bx.bigEnd(0); i > bx.bigEnd(0) - ngrow; i--)
                {
                    for (int j = bx.smallEnd(1) - ngrow; j <= bx.bigEnd(1) + ngrow; j++)
                    {
                        pmask(ig, j, k) = pmask(i, j, k);
                        pmask_temp(ig, j, k) = pmask_temp(i, j, k);
                    }
                    ig++;
                }
            }
            /// top
            if (bx.bigEnd(1) == domain.bigEnd(1))
            {
                int jg = bx.bigEnd(1) + 1;
                for (int j = bx.bigEnd(1); j > bx.bigEnd(1) - ngrow; j--)
                {
                    for (int i = bx.smallEnd(0) - ngrow; i <= bx.bigEnd(0) + ngrow; i++)
                    {
                        pmask(i, jg, k) = pmask(i, j, k);
                        pmask_temp(i, jg, k) = pmask_temp(i, j, k);
                    }
                    jg++;
                }
            }
        }
    }    

    void Mask::FillInVelocityComponents(amrex::MultiFab &U)
    {
        U.FillBoundary();

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
        double x_[N_Max], y_[N_Max], w_[N_Max], U_[N_Max], V_[N_Max];
    
        for (auto &&solid : *interfaces)
        {
            int i_beg = 0;
            if (!solid->isAdvectLS())
                return;
            for (amrex::MFIter mfi(U); mfi.isValid(); ++mfi)
            {
                amrex::Array4<int const> const &pmask = PMask.const_array(mfi);
                amrex::Array4<amrex::Real> const &Vel = U.array(mfi);
                auto &icpt_data = solid->getInterceptData()[mfi];
                amrex::IntVect test_iv(AMREX_D_DECL(115, 0, 0));
    
                for (auto &&idt : icpt_data)
                {
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
                            //amrex::PrintToFile("log") << "Error in identifying interface type.\n";
                            exit(1);
                        }
                        U_[ni] = idt.u[ni];
                        V_[ni] = idt.v[ni];
                        Compute_LSQ_Weights(w_[ni], std::fabs(idt.frac_[ni]));
                    }
                    int Num_i_ = ni;
    
                    /// create grown box around cell
                    amrex::Box gbx(idt.cellid_, idt.cellid_);
                    gbx.grow(stencil_);

                    amrex::Box gbx_isect = gbx & domain;
    
                    for (amrex::BoxIterator bit(gbx); bit.ok(); ++bit)
                    {
                        const amrex::IntVect &iv = bit();
                        if (pmask(iv) == 1)
                        {
                            if(isAxisymmetric)
                                AdvectLSIntVel::Axisymmetric::Vel_LSQ_Parameters(U_[ni], V_[ni], x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], Vel, prob_lo, dx,domain);
                            else
                                Vel_LSQ_Parameters(U_[ni], V_[ni], x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], Vel, prob_lo, dx);
                            ni++;
                        }
			
                        //if( idt.cellid_[0] == 16464 && idt.cellid_[1] == 0)
			//if(iv[0] == 16460 && iv[1] == -4)
                        //{
                        //    amrex::Print(-1)<<" ni = "<<ni<<'\n';
                        //    amrex::Print(-1)<<" iv = "<<iv[0]<<" , "<<iv[1]<<'\n';
			//    amrex::Print(-1)<<"pmask(iv) = "<<pmask(iv)<<'\n';
                        //}
			
                    }
    
                    int Num_ = ni;
    
                    if (solid->isAdvectLS())
                    {
                        i_beg = Num_i_;
                    }

                /*amrex::IntVect test_iv(AMREX_D_DECL(124, 139, 0));
                if(idt.cellid_ == test_iv)
                {
                    amrex::Print()<<"Cell = "<<idt.cellid_<<'\n';
                    amrex::Print()<<"Vel(idt.cellid_, 0) : "<< Vel(idt.cellid_, 0)<<'\n';
                    amrex::Print()<<"Num_ = "<<Num_<<" , Num_i_ = "<<Num_i_<<'\n';
                    for(int ii = 0;ii<= Num_;ii++)
                    {
                            amrex::Print()<<"ii = "<<ii<<"\t";
                            amrex::Print()<<"x_ = "<<x_[ii]<<"\t";
                            amrex::Print()<<"y_ = "<<y_[ii]<<"\t";
                            //amrex::Print()<<"Vel(idt.cellid_, 0) = "<<Vel(idt.cellid_, 0)<<"\t";
                            amrex::Print()<<"U_ = "<<U_[ii]<<'\t';
                            amrex::Print()<<"V_ = "<<V_[ii]<<'\n';
                    }
                    for (ni = 0; ni < idt.n_intercepts; ni++)
                    {
                        amrex::Print()<<"idt.norm_shear_["<< ni <<"] = "<<idt.norm_shear_[ni]<<'\n';
                        amrex::Print()<<"idt.psix_[ni] = "<<idt.psix_[ni]<<" , idt.psiy_[ni] = "<<idt.psiy_[ni]<<'\n';
                    }
                        //std::exit(9);
                }*/

    
                    if (PORDER == 4)
                    {
                        if (Num_ < 3 + Num_i_)
                        {
                            // too few points have been found...can only do a simple average
                            sum_w = sum_sol_u = sum_sol_v = 0.0;
                            for (int i = i_beg; i < Num_; i++)
                            {
                                sum_w += w_[i];
                                sum_sol_u += w_[i] * U_[i];
                                sum_sol_v += w_[i] * V_[i];
                            }
                            Vel(idt.cellid_, 0) = sum_sol_u / sum_w;
                            Vel(idt.cellid_, 1) = sum_sol_v / sum_w;
    
                            //amrex::PrintToFile("log") << " Inside first \n";
                        }
                        else if (Num_ < 6 + Num_i_)
                        {
                            if (solid->isAdvectLS())
                            {
                                if(isAxisymmetric)
                                    AdvectLSIntVel::Axisymmetric::Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, visc_, dx, prob_lo);
                                else 
                                    AdvectLSIntVel::Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, Mu, dx);
                            }
                            else
                                Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_);
                            if (Sing_)
                            {
                                sum_w = sum_sol_u = sum_sol_v = 0.0;
                                for (int i = i_beg; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol_u += w_[i] * U_[i];
                                    sum_sol_v += w_[i] * V_[i];
                                }
                                Vel(idt.cellid_, 0) = sum_sol_u / sum_w;
                                Vel(idt.cellid_, 1) = sum_sol_v / sum_w;
                            }
                            //amrex::PrintToFile("log") << " Second first \n";
                        }
                        else if (Num_ < 20 + Num_i_)
                        {
                            if (solid->isAdvectLS())
                            {
                                if (isAxisymmetric)
                                    AdvectLSIntVel::Axisymmetric::Int_LS_Vel_Third_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, visc_, dx, prob_lo);
                                else
                                    AdvectLSIntVel::Int_LS_Vel_Third_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, Mu, dx);
                            }
                            else
                                Int_LS_Vel_Third_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_);
                            if (Sing_)
                            {
                                if (solid->isAdvectLS())
                                {
                                    if (isAxisymmetric)
                                        AdvectLSIntVel::Axisymmetric::Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, visc_, dx, prob_lo);
                                    else
                                        AdvectLSIntVel::Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, Mu, dx);
                                }
                                else
                                    Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_);
                            }
                            if (Sing_)
                            {
                                sum_w = sum_sol_u = sum_sol_v = 0.0;
                                for (int i = i_beg; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol_u += w_[i] * U_[i];
                                    sum_sol_v += w_[i] * V_[i];
                                }
                                Vel(idt.cellid_, 0) = sum_sol_u / sum_w;
                                Vel(idt.cellid_, 1) = sum_sol_v / sum_w;
                            }
                            //amrex::PrintToFile("log") << " Third first\n";
                        }
                        else
                        {
                            if (solid->isAdvectLS())
                            {
                                if (isAxisymmetric)
                                    AdvectLSIntVel::Axisymmetric::Int_LS_Vel_Fourth_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, visc_, dx, prob_lo);
                                else
                                    AdvectLSIntVel::Int_LS_Vel_Fourth_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, Mu, dx);
                            }
                            else
                                Int_LS_Vel_Fourth_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, Mu, dx);
                            if (Sing_)
                            {
                                if (solid->isAdvectLS())
                                {
                                    if (isAxisymmetric)
                                        AdvectLSIntVel::Axisymmetric::Int_LS_Vel_Third_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, visc_, dx, prob_lo);
                                    else
                                        AdvectLSIntVel::Int_LS_Vel_Third_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, Mu, dx);
                                }
                                else
                                    Int_LS_Vel_Third_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_);
                                //amrex::PrintToFile("log") << " Inside second\n";
                            }
                            if (Sing_)
                            {
                                if (solid->isAdvectLS())
                                {
                                    if (isAxisymmetric)
                                        AdvectLSIntVel::Axisymmetric::Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, visc_, dx, prob_lo);
                                    else
                                        AdvectLSIntVel::Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, Mu, dx);
                                }
                                else
                                    Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_);
                                //amrex::PrintToFile("log") << " Inside first\n";
                            }
                            if (Sing_)
                            {
                                sum_w = sum_sol_u = sum_sol_v = 0.0;
                                for (int i_ = i_beg; i_ < Num_; i_++)
                                {
                                    sum_w += w_[i_];
                                    sum_sol_u += w_[i_] * U_[i_];
                                    sum_sol_v += w_[i_] * V_[i_];
                                }
                                Vel(idt.cellid_, 0) = sum_sol_u / sum_w;
                                Vel(idt.cellid_, 1) = sum_sol_v / sum_w;
                            }
                        }
                    }
                    else if (PORDER == 3)
                    {
                        if (Num_ < 3 + Num_i_)
                        {
                            // too few points have been found...can only do a simple average
                            sum_w = sum_sol_u = sum_sol_v = 0.0;
                            for (int i = i_beg; i < Num_; i++)
                            {
                                sum_w += w_[i];
                                sum_sol_u += w_[i] * U_[i];
                                sum_sol_v += w_[i] * V_[i];
                            }
                            Vel(idt.cellid_, 0) = sum_sol_u / sum_w;
                            Vel(idt.cellid_, 1) = sum_sol_v / sum_w;
    
                            //amrex::PrintToFile("log") << " Inside first \n";
                        }
                        else if (Num_ < 6 + Num_i_)
                        {
                            if (solid->isAdvectLS())
                            {
                                if (isAxisymmetric)
                                    AdvectLSIntVel::Axisymmetric::Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, visc_, dx, prob_lo);
                                else
                                    AdvectLSIntVel::Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, Mu, dx);
                            }
                            else
                                Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_);
                            if (Sing_)
                            {
                                sum_w = sum_sol_u = sum_sol_v = 0.0;
                                for (int i = i_beg; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol_u += w_[i] * U_[i];
                                    sum_sol_v += w_[i] * V_[i];
                                }
                                Vel(idt.cellid_, 0) = sum_sol_u / sum_w;
                                Vel(idt.cellid_, 1) = sum_sol_v / sum_w;
                            }
                            //amrex::PrintToFile("log") << " Second first \n";
                        }
                        else
                        {
                            if (solid->isAdvectLS())
                            {
                                if (isAxisymmetric)
                                    AdvectLSIntVel::Axisymmetric::Int_LS_Vel_Third_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, visc_, dx, prob_lo);
                                else
                                    AdvectLSIntVel::Int_LS_Vel_Third_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, Mu, dx);
                            }
                            else
                                Int_LS_Vel_Third_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_);
                            if (Sing_)
                            {
                                if (solid->isAdvectLS())
                                {
                                    if (isAxisymmetric)
                                        AdvectLSIntVel::Axisymmetric::Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, visc_, dx, prob_lo);
                                    else
                                        AdvectLSIntVel::Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, Mu, dx);
                                }
                                else
                                    Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_);
                            }
                            if (Sing_)
                            {
                                sum_w = sum_sol_u = sum_sol_v = 0.0;
                                for (int i = i_beg; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol_u += w_[i] * U_[i];
                                    sum_sol_v += w_[i] * V_[i];
                                }
                                Vel(idt.cellid_, 0) = sum_sol_u / sum_w;
                                Vel(idt.cellid_, 1) = sum_sol_v / sum_w;
                            }
                        }
                    }
                    else if (PORDER == 2)
                    {
                        if (Num_ < 3 + Num_i_)
                        {
                            // too few points have been found...can only do a simple average
                            sum_w = sum_sol_u = sum_sol_v = 0.0;
                            for (int i = i_beg; i < Num_; i++)
                            {
                                sum_w += w_[i];
                                sum_sol_u += w_[i] * U_[i];
                                sum_sol_v += w_[i] * V_[i];
                            }
                            Vel(idt.cellid_, 0) = sum_sol_u / sum_w;
                            Vel(idt.cellid_, 1) = sum_sol_v / sum_w;
    
                            //amrex::PrintToFile("log") << " Inside first \n";
                        }
                        else
                        {
                            if (solid->isAdvectLS())
                            {
                                if (isAxisymmetric)
                                    AdvectLSIntVel::Axisymmetric::Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, visc_, dx, prob_lo);
                                else
                                    AdvectLSIntVel::Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_, Mu, dx);
                            }
                            else
                                Int_LS_Vel_Second_Order(Sing_, Num_, Num_i_, idt, Vel, x_, y_, U_, V_, w_);
                            if (Sing_)
                            {
                                sum_w = sum_sol_u = sum_sol_v = 0.0;
                                for (int i = i_beg; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol_u += w_[i] * U_[i];
                                    sum_sol_v += w_[i] * V_[i];
                                }
                                Vel(idt.cellid_, 0) = sum_sol_u / sum_w;
                                Vel(idt.cellid_, 1) = sum_sol_v / sum_w;
                            }
                        }
                    }
                    //Vel(idt.cellid_, 0) = 0.0;//sum_sol_u / sum_w;
                    //Vel(idt.cellid_, 1) = 0.0;//sum_sol_v / sum_w;
		    //For Testing source and sink
		    /*
		    {
                        amrex::Real x_gc = prob_lo[0] + dx[0] * (idt.cellid_[0] + 0.5);
                        amrex::Real y_gc = prob_lo[1] + dx[1] * (idt.cellid_[1] + 0.5);
                        amrex::Real r = std::hypot(x_gc , y_gc );
                        amrex::Real theta = std::acos(x_gc/r);
			amrex::Real v_r = -1.0/(3.0*r*r);
                        amrex::Real v_theta = 0.0;

			Vel(idt.cellid_, 0) = v_r*cos(theta) - v_theta*sin(theta); 
			Vel(idt.cellid_, 0) = v_r*sin(theta) + v_theta*cos(theta);
			//amrex::Print()<<"cell : "<<idt.cellid_<<'\n';
                        //amrex::Print()<<"Vel(idt.cellid_) : "<< Vel(idt.cellid_, 0)<<" , "<<Vel(idt.cellid_, 1)<<'\n';

		    }
		    */
		    
                    /* 
                    if(idt.cellid_[0] == 16464 && idt.cellid_[1] == 0 )
                    {
                        amrex::Print(-1)<<"cell : "<<idt.cellid_<<'\n';
                        amrex::Print(-1)<<"Vel(idt.cellid_) : "<< Vel(idt.cellid_, 0)<<" , "<<Vel(idt.cellid_, 1)<<'\n';
                        amrex::Print(-1)<<"Num_ = "<<Num_<<" , Num_i_ = "<<Num_i_<<'\n';
			
                        for(int ii = 0;ii<= Num_;ii++)
                        {
                                amrex::Print(-1)<<"ii = "<<ii<<"\t";
                                amrex::Print(-1)<<"x_ = "<<x_[ii]<<"\t";
                                amrex::Print(-1)<<"y_ = "<<y_[ii]<<"\t";
                                amrex::Print(-1)<<"U_(ii) = "<<U_[ii]<<" , V_(ii) = "<<V_[ii]<<"\n";
                        }
			
                        //std::exit(9);
                    }
		    */ 
                }    
            }
        }      
    }
 
    void Mask::ExtendVelocityLSQ(amrex::MultiFab &U)
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
        double x_[N_Max], y_[N_Max], w_[N_Max], U_[N_Max], V_[N_Max];
    
        for (int outer = 0; outer < LAYERS; outer++)
        {
            U.FillBoundary();
            for (amrex::MFIter mfi(PMask); mfi.isValid(); ++mfi)
            {
                amrex::Array4<int> const &pmask = PMask.array(mfi);
                amrex::Array4<amrex::Real> const &Vel = U.array(mfi);
                auto &Index_Additional_Layers = Index_Additional_Layers_LD[mfi];
    
                for (auto &&cell : Index_Additional_Layers[outer])
                {
                    Sing_ = false;
                    amrex::Box bx(cell, cell);
                    bx.grow(stencil_);
                    amrex::Box bx_isect = bx & domain;

                    int i_ = 0;
                    for (amrex::BoxIterator bit(bx); bit.ok(); ++bit)
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
                                    //Vel_LSQ_Parameters(U_[i_], V_[i_], x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], Vel, prob_lo, dx);
                                    if(isAxisymmetric)
                                        AdvectLSIntVel::Axisymmetric::Vel_LSQ_Parameters(U_[i_], V_[i_], x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], Vel, prob_lo, dx,domain);
                                    else
                                        Vel_LSQ_Parameters(U_[i_], V_[i_], x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], Vel, prob_lo, dx);

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
                                    //Vel_LSQ_Parameters(U_[i_], V_[i_], x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], Vel, prob_lo, dx);
                                    if(isAxisymmetric)
                                        AdvectLSIntVel::Axisymmetric::Vel_LSQ_Parameters(U_[i_], V_[i_], x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], Vel, prob_lo, dx,domain);
                                    else
                                        Vel_LSQ_Parameters(U_[i_], V_[i_], x_[i_], y_[i_], w_[i_], iv[0], iv[1], cell[0], cell[1], Vel, prob_lo, dx);


                                    i_++;
                                }
                            }
                        }
                    }
                    int Num_ = i_;
    
                    if (outer <= LAYERS / 2)
                    {
                        if (isAxisymmetric)
			{
                            Axisymmetric::LS_Vel_Fourth_Order(Sing_, Num_, cell, Vel, x_, y_, U_, V_, w_, dx, prob_lo);
			    //amrex::Print()<<"here"<<'\n';
                        }
                        else
                            LS_Vel_Second_Order(Sing_, Num_, cell, Vel, x_, y_, U_, V_, w_);
                        if (Sing_)
                        {
                            if (isAxisymmetric)
                                Axisymmetric::LS_Vel_Third_Order(Sing_, Num_, cell, Vel, x_, y_, U_, V_, w_, dx, prob_lo);
                            else
                                LS_Vel_Second_Order(Sing_, Num_, cell, Vel, x_, y_, U_, V_, w_);
                        }
                        if (Sing_)
                        {
                            if (isAxisymmetric)
                                Axisymmetric::LS_Vel_Second_Order(Sing_, Num_, cell, Vel, x_, y_, U_, V_, w_, dx, prob_lo);
                            else
                                LS_Vel_Second_Order(Sing_, Num_, cell, Vel, x_, y_, U_, V_, w_);
                        }
                    }
		    else
                    {
                        if (isAxisymmetric)
                            Axisymmetric::LS_Vel_Third_Order(Sing_, Num_, cell, Vel, x_, y_, U_, V_, w_, dx, prob_lo);

                        else LS_Vel_Second_Order(Sing_, Num_, cell, Vel, x_, y_, U_, V_, w_);
                        if (Sing_)
                        {
                            if (isAxisymmetric)
                                Axisymmetric::LS_Vel_Second_Order(Sing_, Num_, cell, Vel, x_, y_, U_, V_, w_, dx, prob_lo);
                            else
                                LS_Vel_Second_Order(Sing_, Num_, cell, Vel, x_, y_, U_, V_, w_);
                        }
                    }
		    /*
                    if(cell[0] == 16462 && cell[1] == 0)
                    {
                        amrex::Print(-1)<<"i = "<<cell[0]<<" , j = "<<cell[1]<<'\n';
                        amrex::Print(-1)<<"outer = "<<outer<<'\n';
                        amrex::Print(-1)<<"Num_ = "<<Num_<<'\n';

                        for(int ii = 0; ii<Num_;ii++)
                        {
                            amrex::Print(-1)<<"ii = "<<ii<<'\n';
                            amrex::Print(-1)<<"x_[ii] = "<<x_[ii]<<" , y_[ii] = "<<y_[ii]<<'\n';
                            amrex::Print(-1)<<"U_[ii] = "<<U_[ii]<<'\n';
                        }
		    }
		    */
                }
            }
            //CollocatedVelBoundaryConditions();
        }
    }    
    
    void Mask::FillInRefConfig(amrex::MultiFab &X, amrex::MultiFab &Y, amrex::MultiFab &DMG)
    {
        X.FillBoundary();
        Y.FillBoundary();
	    DMG.FillBoundary();
        
        if (interfaces->empty())
        {
            return;
        }
        const amrex::Real *prob_lo = geom_.ProbLo();
        const amrex::Real *dx = geom_.CellSize();
        const amrex::Box &domain = geom_.Domain();
    
        bool Sing_;
        double sum_w, sum_sol_X, sum_sol_Y, sum_sol_DMG;
    
        for (auto &&solid : *interfaces)
        {
            int i_beg = 0;
            if (!solid->isAdvectLS())
                return;
            for (amrex::MFIter mfi(X); mfi.isValid(); ++mfi)
            {
                amrex::Array4<int const> const &pmask = PMask.const_array(mfi);
                amrex::Array4<amrex::Real> const &X_ref = X.array(mfi);
                amrex::Array4<amrex::Real> const &Y_ref = Y.array(mfi);
		amrex::Array4<amrex::Real> const &dmg = DMG.array(mfi);
                auto &icpt_data = solid->getInterceptData()[mfi];
                amrex::IntVect test_iv(AMREX_D_DECL(115, 0, 0));
    
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
	            int N_Max = (2 * stencil_ + 1) * (2 * stencil_ + 1) + 4;
                    double x_[N_Max], y_[N_Max], w_[N_Max], X_[N_Max], Y_[N_Max], DMG_[N_Max];

                    Sing_ = false; // Need this. Sing_ is only changed to true.
                    int ni;
                    //for (ni = 0; ni < idt.n_intercepts; ni++)
                    //{
                     //   x_[ni] = y_[ni] = 0.0;
                     //   if (idt.type_[ni] == 0)
                     //   {
                     //       x_[ni] = -idt.frac_[ni];
                     //   }
                     //   else if (idt.type_[ni] == 1)
                     //   {
                     //       y_[ni] = -idt.frac_[ni];
                     //   }
                     //   else
                     //   {
                     //       //amrex::PrintToFile("log") << "Error in identifying interface type.\n";
                     //       exit(1);
                     //   }
                     //   U_[ni] = idt.u[ni];
                     //   V_[ni] = idt.v[ni];
                     //   Compute_LSQ_Weights(w_[ni], std::fabs(idt.frac_[ni]));
                    //}
                    ni = 0; //Interfacial points are not included in the least square cloud
                    int Num_i_ = ni;
    
                    /// create grown box around cell
                    amrex::Box gbx(idt.cellid_, idt.cellid_);
                    gbx.grow(idt.stencil_);

                    amrex::Box gbx_isect = gbx & domain;
    
                    for (amrex::BoxIterator bit(gbx); bit.ok(); ++bit)
                    {
                        const amrex::IntVect &iv = bit();
                        if (pmask(iv) == 1)
                        {
                            //if(isAxisymmetric)
                            //    AdvectLSIntVel::Axisymmetric::RefConfg_LSQ_Parameters(X_[ni], Y_[ni], x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], X_ref, Y_ref, prob_lo, dx,domain);
                            //else
                            //    RefConfg_LSQ_Parameters(X_[ni], Y_[ni], x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], X_ref, Y_ref, prob_lo, dx);
			    LSQ_Parameters(X_[ni], x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], X_ref, prob_lo, dx);
			    LSQ_Parameters(Y_[ni], x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], Y_ref, prob_lo, dx);
			    LSQ_Parameters(DMG_[ni], x_[ni], y_[ni], w_[ni], iv[0], iv[1], idt.cellid_[0], idt.cellid_[1], dmg, prob_lo, dx);
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

                    /*amrex::IntVect test_iv(AMREX_D_DECL(124, 139, 0));
                    if(idt.cellid_ == test_iv)
                    {
                        amrex::Print()<<"Cell = "<<idt.cellid_<<'\n';
                        amrex::Print()<<"Vel(idt.cellid_, 0) : "<< Vel(idt.cellid_, 0)<<'\n';
                        amrex::Print()<<"Num_ = "<<Num_<<" , Num_i_ = "<<Num_i_<<'\n';
                        for(int ii = 0;ii<= Num_;ii++)
                        {
                                amrex::Print()<<"ii = "<<ii<<"\t";
                                amrex::Print()<<"x_ = "<<x_[ii]<<"\t";
                                amrex::Print()<<"y_ = "<<y_[ii]<<"\t";
                                //amrex::Print()<<"Vel(idt.cellid_, 0) = "<<Vel(idt.cellid_, 0)<<"\t";
                                amrex::Print()<<"U_ = "<<U_[ii]<<'\t';
                                amrex::Print()<<"V_ = "<<V_[ii]<<'\n';
                        }
                        for (ni = 0; ni < idt.n_intercepts; ni++)
                        {
                            amrex::Print()<<"idt.norm_shear_["<< ni <<"] = "<<idt.norm_shear_[ni]<<'\n';
                            amrex::Print()<<"idt.psix_[ni] = "<<idt.psix_[ni]<<" , idt.psiy_[ni] = "<<idt.psiy_[ni]<<'\n';
                        }
                            //std::exit(9);
                    }*/

                    int PORDER_ = 2;
		    //idt.PORDER = 2;
		    PORDER_ = idt.PORDER;
                    if (PORDER_ >= 3)
                    {
                        if (Num_ < 3 + Num_i_)
                        {
                            /// too few points have been found...can only do a simple average
                            sum_w = sum_sol_X = sum_sol_Y = sum_sol_DMG = 0.0;
                            for (int i = i_beg; i < Num_; i++)
                            {
                                sum_w += w_[i];
                                sum_sol_X += w_[i] * X_[i];
                                sum_sol_Y += w_[i] * Y_[i];
                                sum_sol_DMG += w_[i] * DMG_[i];
                            }
                            X_ref(idt.cellid_) = sum_sol_X / sum_w;
                            Y_ref(idt.cellid_) = sum_sol_Y / sum_w;
                            dmg(idt.cellid_) = sum_sol_DMG / sum_w;
                        }
                        else if (Num_ < 10 + Num_i_)
                        {
                            QR_LS_Temperature_Second_Order(Sing_, Num_, Num_i_, x_, y_, X_, w_, X_ref, idt);
                            QR_LS_Temperature_Second_Order(Sing_, Num_, Num_i_, x_, y_, Y_, w_, Y_ref, idt);
			    QR_LS_Temperature_Second_Order(Sing_, Num_, Num_i_, x_, y_, DMG_, w_, dmg, idt);
                            if (Sing_)
                            {
                                sum_w = sum_sol_X = sum_sol_Y = sum_sol_DMG = 0.0;
                                for (int i = i_beg; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol_X += w_[i] * X_[i];
                                    sum_sol_Y += w_[i] * Y_[i];
                                    sum_sol_DMG += w_[i] * DMG_[i];
                                }
                                X_ref(idt.cellid_) = sum_sol_X / sum_w;
                                Y_ref(idt.cellid_) = sum_sol_Y / sum_w;
                                dmg(idt.cellid_) = sum_sol_DMG / sum_w;
                            }
                            //amrex::PrintToFile("log") << " Second first ";
                        }
                        else
                        {
                            QR_LS_Temperature_Third_Order(Sing_, Num_, Num_i_, x_, y_, X_, w_, X_ref, idt);
                            QR_LS_Temperature_Third_Order(Sing_, Num_, Num_i_, x_, y_, Y_, w_, Y_ref, idt);
			    QR_LS_Temperature_Third_Order(Sing_, Num_, Num_i_, x_, y_, DMG_, w_, dmg, idt);
			    //amrex::Print(-1)<<"using third order at "<<idt.cellid_<<'\n';
                            if (Sing_)
                            {
	                        QR_LS_Temperature_Second_Order(Sing_, Num_, Num_i_, x_, y_, X_, w_, X_ref, idt);
                                QR_LS_Temperature_Second_Order(Sing_, Num_, Num_i_, x_, y_, Y_, w_, Y_ref, idt);
				QR_LS_Temperature_Second_Order(Sing_, Num_, Num_i_, x_, y_, DMG_, w_, dmg, idt);
                            }
                            if (Sing_)
                            {
                                sum_w = sum_sol_X = sum_sol_Y = sum_sol_DMG = 0.0;
                                for (int i = i_beg; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol_X += w_[i] * X_[i];
                                    sum_sol_Y += w_[i] * Y_[i];
                                    sum_sol_DMG += w_[i] * DMG_[i];
                                }
                                X_ref(idt.cellid_) = sum_sol_X / sum_w;
                                Y_ref(idt.cellid_) = sum_sol_Y / sum_w;
                                dmg(idt.cellid_) = sum_sol_DMG / sum_w;
                            }
                            //amrex::PrintToFile("log") << " Third first ";
                        }
                    }		    
		    else if (PORDER_ == 2)
                    {
                        if (Num_ < 3 + Num_i_)
                        {
                            // too few points have been found...can only do a simple average
                            sum_w = sum_sol_X = sum_sol_Y = sum_sol_DMG = 0.0;
                            for (int i = i_beg; i < Num_; i++)
                            {
                                sum_w += w_[i];
                                sum_sol_X += w_[i] * X_[i];
                                sum_sol_Y += w_[i] * Y_[i];
				sum_sol_DMG += w_[i] * DMG_[i];
                            }
                            X_ref(idt.cellid_) = sum_sol_X / sum_w;
                            Y_ref(idt.cellid_) = sum_sol_Y / sum_w;
			    dmg(idt.cellid_) = sum_sol_DMG / sum_w;
                            //amrex::PrintToFile("log") << " Inside first \n";
                        }
                        else
                        {
		            QR_LS_Temperature_Second_Order(Sing_, Num_, Num_i_, x_, y_, X_, w_, X_ref, idt);
			    QR_LS_Temperature_Second_Order(Sing_, Num_, Num_i_, x_, y_, Y_, w_, Y_ref, idt);
			    QR_LS_Temperature_Second_Order(Sing_, Num_, Num_i_, x_, y_, DMG_, w_, dmg, idt);
                            if (Sing_)
                            {
                                sum_w = sum_sol_X = sum_sol_Y = sum_sol_DMG = 0.0;
                                for (int i = i_beg; i < Num_; i++)
                                {
                                    sum_w += w_[i];
                                    sum_sol_X += w_[i] * X_[i];
                                    sum_sol_Y += w_[i] * Y_[i];
				    sum_sol_DMG += w_[i] * DMG_[i];
                                }
                                X_ref(idt.cellid_) = sum_sol_X / sum_w;
                                Y_ref(idt.cellid_) = sum_sol_Y / sum_w;
                                dmg(idt.cellid_) = sum_sol_DMG / sum_w;
                            }
                        }
                    }
                    //Vel(idt.cellid_, 0) = 0.0;//sum_sol_u / sum_w;
                    //Vel(idt.cellid_, 1) = 0.0;//sum_sol_v / sum_w;
                /*if(idt.cellid_ == test_iv)
                {
                        amrex::Print()<<"cell : "<<idt.cellid_<<'\n';
                        amrex::Print()<<"Vel(idt.cellid_) : "<< Vel(idt.cellid_, 0)<<" , "<<Vel(idt.cellid_, 1)<<'\n';
                        amrex::Print()<<"Num_ = "<<Num_<<" , Num_i_ = "<<Num_i_<<'\n';
                        for (ni = 0; ni < idt.n_intercepts; ni++)
                        {

                            amrex::Print()<<"idt.norm_shear_["<< ni <<"] = "<<idt.norm_shear_[ni]<<'\n';
                        }
                        //for(int ii = 0;ii<= Num_;ii++)
                        //{
                        //        amrex::Print()<<"ii = "<<ii<<"\t";
                        //        amrex::Print()<<"x_ = "<<x_[ii]<<"\t";
                        //        amrex::Print()<<"y_ = "<<y_[ii]<<"\t";
                        //        amrex::Print()<<"Vel(idt.cellid_, 0) = "<<Vel(idt.cellid_, 0)<<"\t";
                        //}
                        //std::exit(9);
                }*/

                }    
            }
        }      
    } 
} /*End namespace mycode */

