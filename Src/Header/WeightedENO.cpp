#include <WeightedENO.h>
#include <iostream>
#include <cmath>

MatrixPtr *AllocMatR(int size1, int size2) {
        MatrixPtr* v = new MatrixPtr[size1];
        for(int i=0;i<size1;i++)
                v[i] = new double[size2];
        return v;
}

void Allocate_2D_R(double**& m, int d1, int d2) {
        m=new double* [d1];
        for (int i=0; i<d1; ++i) {
                m[i]=new double [d2];
                for (int j=0; j<d2; ++j)
                        m[i][j]=0.0;
        }
}

void Allocate_3D_R(double***& m, int d1, int d2, int d3) {
        m=new double** [d1];
        for (int i=0; i<d1; ++i) {
                m[i]=new double* [d2];
                for (int j=0; j<d2; ++j) {
                        m[i][j]=new double [d3];
                        for (int k=0; k<d3; ++k)
                                m[i][j][k]=0.0;
                }
        }
}

void Allocate_4D_R(double****& m, int d1, int d2, int d3, int d4) {
        m=new double*** [d1];
        for (int i=0; i<d1; ++i) {
                m[i]=new double** [d2];
                for (int j=0; j<d2; ++j) {
                        m[i][j]=new double* [d3];
                        for (int k=0; k<d3; ++k) {
                        	m[i][j][k]=new double [d4];
                        	for (int l=0; l<d4; ++l) m[i][j][k][l]=0.0;
			}
                }
        }
}

void Allocate_2D_I(int**& m, int d1, int d2) {
        m=new int* [d1];
        for (int i=0; i<d1; ++i) {
                m[i]=new int [d2];
                for (int j=0; j<d2; ++j) {
                        m[i][j]= 0;
                }
        }
}

void Allocate_3D_I(int***& m, int d1, int d2, int d3) {
        m=new int** [d1];
        for (int i=0; i<d1; ++i) {
                m[i]=new int* [d2];
                for (int j=0; j<d2; ++j) {
                        m[i][j]=new int [d3];
                        for (int k=0; k<d3; ++k)
                                m[i][j][k]=0;
                }
        }
}

void Allocate_4D_I(int****& m, int d1, int d2, int d3, int d4) {
        m=new int*** [d1];
        for (int i=0; i<d1; ++i) {
                m[i]=new int** [d2];
                for (int j=0; j<d2; ++j) {
                        m[i][j]=new int* [d3];
                        for (int k=0; k<d3; ++k) {
                        	m[i][j][k]=new int [d4];
                        	for (int l=0; l<d4; ++l) m[i][j][k][l]=0;
			}
                }
        }
}

int ISign(int x) {
	if(x == 0) {
		std::cout << "Error computing Sign of zero !" << "\n" ;
		exit(1) ;
	}
	return x/abs(x) ;
}

double Max(double a, double b) {
        return (a > b ? a : b);
}

double SIGN(double a, double b) {
        return (b >= 0.0 ? fabs(a) : -fabs(a));
}

// Sign function
int I_DSign(double x) {
        if(x >= 0.0) return 1 ;
        else return -1 ;
}
double Sign(double x) {
        if(x >= 0.0) return 1.0 ;
        else return -1.0 ;
}

double Minimum2(double x, double y) {
        if(x > y) return y ;
        else return x;
}

double Maximum2(double x, double y) {
        if(x > y) return x ;
        else return y;
}

double Minimum3(double x, double y, double z) {
        if(x < y) return Minimum2(x,z) ; 
        else return Minimum2(y,z) ;
}

double Maximum3(double x, double y, double z) {
        if(x > y) return Maximum2(x,z) ;
        else return Maximum2(y,z) ;
}

double Maximum4(double a0, double a1, double a2, double a3) {
        if(a0 < a1) return Maximum3(a1,a2,a3) ;
        else return Maximum3(a0,a2,a3) ; 
}

double Maximum5(double a0, double a1, double a2, double a3, double a4) {
        if(a0 < a1) return Maximum4(a1,a2,a3,a4) ;
        else return Maximum4(a0,a2,a3,a4) ; 
}

double Maximum6(double a0, double a1, double a2, double a3, double a4, double a5) {
        if(a0 < a1) return Maximum5(a1,a2,a3,a4,a5) ;
        else return Maximum5(a0,a2,a3,a4,a5) ; 
}

double Minimum4(double a0, double a1, double a2, double a3) {
        if(a0 < a1) return Minimum3(a0,a2,a3) ;
        else return Minimum3(a1,a2,a3) ; 
}

double Minimum5(double a0, double a1, double a2, double a3, double a4) {
        if(a0 < a1) return Minimum4(a0,a2,a3,a4) ;
        else return Minimum4(a1,a2,a3,a4) ; 
}

double Minimum6(double a0, double a1, double a2, double a3, double a4, double a5) {
        if(a0 < a1) return Minimum5(a0,a2,a3,a4,a5) ;
        else return Minimum5(a1,a2,a3,a4,a5) ; 
}

// Minmod function with two variables
double MinMod2(double x, double y) {
        if( (x > 0) && (y > 0) ) {
                if( x > y ) return y ;
                else return x ;
        } else if( (x < 0) && (y < 0) ) {
                if( x < y ) return y ;
                else return x ;
        } else { return 0.0 ;}
}

//Minmod function with three variables
double MinMod3(double x, double y, double z) {
        if( (x > 0) && (y > 0) ) {
                if( x > y ) return MinMod2(y,z) ;
                else return MinMod2(x,z) ;
        } else if( (x < 0) && (y < 0) ) {
                if( x < y ) return MinMod2(y,z) ;
                else return MinMod2(x,z) ;
        } else { return 0.0 ;}
}

// Minmod function with four variables
double MinMod4(double a0, double a1, double a2, double a3) {
        if( (a0 > 0) && (a1 > 0) && (a2 > 0) && (a3 > 0) ) {
                return Minimum4(a0,a1,a2,a3) ;
        } else if( (a0 < 0) && (a1 < 0) && (a2 < 0) && (a3 < 0) ) {
                return -Minimum4(-a0, -a1, -a2, -a3) ;
        } else { return 0.0 ; }
}

double MinMod6(double a0, double a1, double a2, double a3, double a4, double a5) {
	if( (a0 > 0.0) && (a1 > 0.0) && (a2 > 0.0) && (a3 > 0.0) && (a4 > 0.0) && (a5 > 0.0) ) {
		return Minimum6(a0,a1,a2,a3,a4,a5) ;
	} else if( (a0 < 0.0) && (a1 < 0.0) && (a2 < 0.0) && (a3 < 0.0) && (a4 < 0.0) && (a5 < 0.0) ) {
		return -Minimum6(-a0,-a1,-a2,-a3,-a4,-a5) ;
	} else { return 0.0 ; }
}

// Median function
double Median(double x, double y, double z) {
        if( (y > x) && (z > x) ) {
                if(y > z) return z;     // y > z > x
                else return y ;         // z > y > x
        } else if( (y < x) && (z < x) ) {
                if(y < z) return z;     // x > z > y 
                else return y ;         // x > y > z
        } else { return x ;}
}

// Henricks mapping function for mapped WENO...
double MappedWENOFunction(double Omega, double C) {
	return Omega*(C + C*C - 3.0*C*Omega + Omega*Omega)/(C*C + Omega - 2.0*C*Omega) ;
}

void MappedWeights(double*& Omega_, double* d_, int m, int flag) {
	int i ;
	double *OmegaM, *alphaM, aM ;
	OmegaM = new double[m] ; alphaM = new double[m] ;
	aM = 0.0 ;
	for(i = 0 ; i < m ; i++) {
		if(flag == -1) alphaM[i] = MappedWENOFunction(Omega_[i],d_[m-1-i]) ;
		else if(flag == 1) alphaM[i] = MappedWENOFunction(Omega_[i],d_[i]) ;
		else {
			std::cout << "Error in flag! exiting" << "\n" ; exit(1) ;
		}
		aM += alphaM[i] ;
	}
	for(i = 0 ; i < m ; i++) {
		Omega_[i] = alphaM[i]/aM ;
	}
	delete [] OmegaM ; delete [] alphaM ;
}


// The notation etc. is followed from the NASA/CR-97-206253 Report by Chi-Wang Shu
// k denotes the order of accuracy. C_{rj} : r denotes the stencil and j is the summation index.
// The coefficients have been tabulated here for convenience, they can also be computed.
// Smoothness indicators should be computed in the code itself. Coefficients from Balsara and Shu.

void GetTabulatedWENOCoefficientsUniform(MatrixPtr*& Crj, double*& d, int k) {
	if(k == 1) {
		Crj[0][0] = 1.0 ;
		Crj[1][0] = 1.0 ; 
		d[0] = 1.0 ; // No smoothness indicator required here.
	} else if(k == 2) {
		Crj[0][0] = 3.0/2.0  ; Crj[0][1] = -1.0/2.0 ;
		Crj[1][0] = 1.0/2.0  ; Crj[1][1] = 1.0/2.0 ;
		Crj[2][0] = -1.0/2.0 ; Crj[2][1] = 3.0/2.0 ;
		d[0] = 2.0/3.0 ; d[1] = 1.0/3.0 ;
	} else if(k == 3) {
		Crj[0][0] = 11.0/6.0 ; Crj[0][1] = -7.0/6.0 ; Crj[0][2] = 1.0/3.0 ;
		Crj[1][0] = 1.0/3.0  ; Crj[1][1] = 5.0/6.0  ; Crj[1][2] = -1.0/6.0 ;
		Crj[2][0] = -1.0/6.0 ; Crj[2][1] = 5.0/6.0  ; Crj[2][2] = 1.0/3.0 ;
		Crj[3][0] = 1.0/3.0  ; Crj[3][1] = -7.0/6.0 ; Crj[3][2] = 11.0/6.0 ;
		d[0] = 3.0/10.0 ; d[1] = 3.0/5.0 ; d[2] = 1.0/10.0 ;
	} else if(k == 4) {
		Crj[0][0] = 25.0/12.0 ; Crj[0][1] = -23.0/12.0 ; Crj[0][2] = 13.0/12.0  ; Crj[0][3] = -1.0/4.0 ;
		Crj[1][0] = 1.0/4.0   ; Crj[1][1] = 13.0/12.0  ; Crj[1][2] = -5.0/12.0  ; Crj[1][3] = 1/12.0 ;
		Crj[2][0] = -1.0/12.0 ; Crj[2][1] = 7.0/12.0   ; Crj[2][2] = 7.0/12.0   ; Crj[2][3] = -1.0/12.0 ;
		Crj[3][0] = 1.0/12.0  ; Crj[3][1] = -5.0/12.0  ; Crj[3][2] = 13.0/12.0  ; Crj[3][3] = 1.0/4.0 ;
		Crj[4][0] = -1.0/4.0  ; Crj[4][1] = 13.0/12.0  ; Crj[4][2] = -23.0/12.0 ; Crj[4][3] = 25.0/12.0 ;
		d[0] = 4.0/35.0 ; d[1] = 18.0/35.0 ; d[2] = 12.0/35.0 ; d[3] = 1.0/35.0 ;
	} else if(k == 5) {
		Crj[0][0] = 137.0/60.0 ; Crj[0][1] = -163.0/60.0 ; Crj[0][2] = 137.0/60.0 ; Crj[0][3] = -21.0/20.0  ; Crj[0][4] = 1.0/5.0 ;
		Crj[1][0] = 1.0/5.0    ; Crj[1][1] = 77.0/60.0   ; Crj[1][2] = -43.0/60.0 ; Crj[1][3] = 17.0/60.0   ; Crj[1][4] = -1.0/20.0 ;
		Crj[2][0] = -1.0/20.0  ; Crj[2][1] = 9.0/20.0    ; Crj[2][2] = 47.0/60.0  ; Crj[2][3] = -13.0/60.0  ; Crj[2][4] = 1.0/30.0 ;
		Crj[3][0] = 1.0/30.0   ; Crj[3][1] = -13.0/60.0  ; Crj[3][2] = 47.0/60.0  ; Crj[3][3] = 9.0/20.0    ; Crj[3][4] = -1.0/20.0 ;
		Crj[4][0] = -1.0/20.0  ; Crj[4][1] = 17.0/60.0   ; Crj[4][2] = -43.0/60.0 ; Crj[4][3] = 77.0/60.0   ; Crj[4][4] = 1.0/5.0 ;
		Crj[5][0] = 1.0/5.0    ; Crj[5][1] = -21.0/20.0  ; Crj[5][2] = 137.0/60.0 ; Crj[5][3] = -163.0/60.0 ; Crj[5][4] = 137.0/60.0 ;
		d[0] = 5.0/126.0 ; d[1] = 20.0/63.0 ; d[2] = 10.0/21.0 ; d[3] = 10.0/63.0 ; d[4] = 1.0/126.0 ;
	} else if(k == 6) {
		Crj[0][0] = 49.0/20.0	; Crj[0][1] = -71.0/20.0	; Crj[0][2] = 79.0/20.0		; Crj[0][3] = -163.0/60.0	; Crj[0][4] = 31.0/30.0		; Crj[0][5] = -1.0/6.0 ;
		Crj[1][0] = 1.0/6.0	; Crj[1][1] = 29.0/20.0		; Crj[1][2] = -21.0/20.0	; Crj[1][3] = 37.0/60.0		; Crj[1][4] = -13.0/60.0	; Crj[1][5] = 1.0/30.0 ;
		Crj[2][0] = -1.0/30.0	; Crj[2][1] = 11.0/30.0		; Crj[2][2] = 19.0/20.0		; Crj[2][3] = -23.0/60.0	; Crj[2][4] = 7.0/60.0		; Crj[2][5] = -1.0/60.0 ;
		Crj[3][0] = 1.0/60.0	; Crj[3][1] = -2.0/15.0		; Crj[3][2] = 37.0/60.0		; Crj[3][3] = 37.0/60.0		; Crj[3][4] = -2.0/15.0		; Crj[3][5] = 1.0/60.0 ;
		Crj[4][0] = -1.0/60.0	; Crj[4][1] = 7.0/60.0		; Crj[4][2] = -23.0/60.0	; Crj[4][3] = 19.0/20.0		; Crj[4][4] = 11.0/30.0		; Crj[4][5] = -1.0/30.0 ;
		Crj[5][0] = 1.0/30.0	; Crj[5][1] = -13.0/60.0	; Crj[5][2] = 37.0/60.0		; Crj[5][3] = -21.0/20.0	; Crj[5][4] = 29.0/20.0		; Crj[5][5] = 1.0/6.0 ;
		Crj[6][0] = -1.0/6.0	; Crj[6][1] = 31.0/30.0		; Crj[6][2] = -163.0/60.0	; Crj[6][3] = 79.0/20.0		; Crj[6][4] = -71.0/20.0	; Crj[6][5] = 49.0/20.0 ;
		d[0] = 1.0/77.0 ; d[1] = 25.0/154.0 ; d[2] = 100.0/231.0 ; d[3] = 25.0/77.0 ; d[4] = 5.0/77.0 ; d[5] = 1.0/462.0 ;
	} else if(k == 7) { // Incomplete.
		Crj[0][0] = 363.0/140.0 ; Crj[0][1] = -617.0/140.0 ; Crj[0][2] = 853.0/140.0 ; Crj[0][3] = -2341.0/420.0 ; Crj[0][4] = 667.0/210.0 ; Crj[0][5] = -43.0/42.0 ; Crj[0][6] = 1.0/7.0 ;
		Crj[1][0] = 1.0/7.0 ; Crj[1][1] = 223.0/140.0 ; Crj[1][2] = -197.0/140.0 ; Crj[1][3] = 153.0/140.0 ; Crj[1][4] = -241.0/420.0 ; Crj[1][5] = 37.0/210.0 ; Crj[1][6] = -1.0/42.0 ;
		Crj[2][0] = -1.0/42.0 ; Crj[2][1] = 13.0/42.0 ; Crj[2][2] = 153.0/140.0 ; Crj[2][3] = -241.0/420.0 ; Crj[2][4] = 109.0/420.0 ; Crj[2][5] = -31.0/420.0 ; Crj[2][6] = 1.0/105.0 ;
		Crj[3][0] = 1.0/105.0 ; Crj[3][1] = -19.0/420.0 ; Crj[3][2] = 107.0/210.0 ; Crj[3][3] = 319.0/420.0 ; Crj[3][4] = -101.0/420.0 ; Crj[3][5] = 5.0/84.0 ; Crj[3][6] = -1.0/140.0 ;
		Crj[4][0] = -1.0/140.0 ; Crj[4][1] = 5.0/84.0 ; Crj[4][2] = -101.0/420.0 ; Crj[4][3] = 319.0/420.0 ; Crj[4][4] = 107.0/210.0 ; Crj[4][5] = -19.0/420.0 ; Crj[4][6] = 1.0/105.0 ;
		Crj[5][0] = 1.0/105.0 ; Crj[5][1] = -31.0/420.0 ; Crj[5][2] = 109.0/420.0 ; Crj[5][3] = -241.0/420.0 ; Crj[5][4] = 153.0/140.0 ; Crj[5][5] = 13.0/42.0 ; Crj[5][6] = -1.0/42.0 ;
		Crj[6][0] = -1.0/42.0 ; Crj[6][1] = 37.0/210.0 ; Crj[6][2] = -241.0/420.0 ; Crj[6][3] = 153.0/140.0 ; Crj[6][4] = -197.0/140.0 ; Crj[6][5] = 223.0/140.0 ; Crj[6][6] = 1.0/7.0 ;
		Crj[7][0] = 1.0/7.0 ; Crj[7][1] = -43.0/42.0 ; Crj[7][2] = 667.0/210.0 ; Crj[7][3] = -2341.0/420.0 ; Crj[7][4] = 853.0/140.0 ; Crj[7][5] = -667.0/210.0 ; Crj[7][6] = 363.0/140.0 ;
	}
}

// size of Beta is k, and that of V is 2k-1.
void GetTabulatedWENOSmoothnessIndicatorUniform(double*& Beta, double* V, int k) {
	if(k == 1) { // No indicator required  due to a single stencil.
		Beta[0] = 1.0 ; // Any value will do since we have only one choice.
	} else if(k == 2) {
		// (i-1) -> 0, i -> 1, (i+1) -> 2.
		Beta[0] = (V[2] - V[1])*(V[2] - V[1]) ;
		Beta[1] = (V[1] - V[0])*(V[1] - V[0]) ; 
	} else if(k == 3) {
		// (i-2) -> 0, (i-1) -> 1, i -> 2, (i+1) -> 3, (i+2) -> 4. 
		Beta[0] = (13.0/12.0)*(V[2] - 2.0*V[3] + V[4])*(V[2] - 2.0*V[3] + V[4]) + (1.0/4.0)*(3.0*V[2]-4.0*V[3]+V[4])*(3.0*V[2]-4.0*V[3]+V[4]) ;
		Beta[1] = (13.0/12.0)*(V[1] - 2.0*V[2] + V[3])*(V[1] - 2.0*V[2] + V[3]) + (1.0/4.0)*(V[1]-V[3])*(V[1]-V[3]) ;
		Beta[2] = (13.0/12.0)*(V[0] - 2.0*V[1] + V[2])*(V[0] - 2.0*V[1] + V[2]) + (1.0/4.0)*(V[0]-4.0*V[1]+3.0*V[2])*(V[0]-4.0*V[1]+3.0*V[2]) ;
	} else if(k == 4) {
		// (i-3) -> 0, (i-2) -> 1, (i-1) -> 2, i -> 3, (i+1) -> 4, (i+2) -> 5, (i+3) -> 7. 
		Beta[0] = 547*V[6]*V[6]  + V[5]*(7043*V[5] - 3882*V[6]) + V[4]*( 11003*V[4] - 17246*V[5] + 4642*V[6]) + V[3]*(2107*V[3] - 9402*V[4] + 7042*V[5] -1854*V[6] )  ;	
		Beta[1] = 267*V[5]*V[5] + V[4]*(2843*V[4] - 1642*V[5]) + V[3]*(3443*V[3] - 5966*V[4] + 1602*V[5]) + V[2]*(547*V[2] - 2522*V[3] + 1922*V[4] - 494*V[5]) ;	
		Beta[2] = 547*V[4]*V[4] + V[3]*(3443*V[3] - 2522*V[4]) + V[2]*(2843*V[2] - 5966*V[3] + 1922*V[4]) + V[1]*(267*V[1] - 1642*V[2] + 1602*V[3] - 494*V[4]) ;
		Beta[3] = 2107*V[3]*V[3] + V[2]*(11003*V[2] - 9402*V[3]) + V[1]*(7043*V[1] - 17246*V[2] + 7042*V[3]) + V[0]*(547*V[0] - 3882*V[1] + 4642*V[2] - 1854*V[3]) ;	
	} else if(k == 5) {
		Beta[0] = 22658*V[8]*V[8] + V[7]*(482963*V[7] - 208501*V[8]) + V[6]*(1521393*V[6] - 1704396*V[7] + 364863*V[8]) + V[5]*(1020563*V[5] - 2462076*V[6] + 1358458*V[7] - 288007*V[8]) + V[4]*(107918*V[4] - 649501*V[5] + 758823*V[6] - 411487*V[7] + 86329*V[8]) ;
		Beta[1] = 6908*V[7]*V[7] + V[6]*(138563*V[6] - 60871*V[7]) + V[5]*(406293*V[5] - 464976*V[6] + 99213*V[7]) + V[4]*(242723*V[4] - 611976*V[5] + 337018*V[6] - 70237*V[7]) + V[3]*(22658*V[3] - 140251*V[4] + 165153*V[5] - 88297*V[6] + 18079*V[7]) ;
		Beta[2] = 6908*V[6]*V[6] + V[5]*(104963*V[5] - 51001*V[6]) + V[4]*(231153*V[4] - 299076*V[5] + 67923*V[6]) + V[3]*(104963*V[3] - 299076*V[4] + 179098*V[5] - 38947*V[6]) + V[2]*(6908*V[2] - 51001*V[3] + 67923*V[4] - 38947*V[5] + 8209*V[6]) ;
		Beta[3] = 22658*V[5]*V[5] + V[4]*(242723*V[4] - 140251*V[5]) + V[3]*(406293*V[3] - 611976*V[4] + 165153*V[5]) + V[2]*(138563*V[2] - 464976*V[3] + 337018*V[4] - 88297*V[5]) + V[1]*(6908*V[1] - 60871*V[2] + 99213*V[3] - 70237*V[4] + 18079*V[5]) ;
		Beta[4] = 107918*V[4]*V[4] + V[3]*(1020563*V[3] - 649501*V[4]) + V[2]*(1521393*V[2] - 2462076*V[3] + 758823*V[4]) + V[1]*(482963*V[1] - 1704396*V[2] + 1358458*V[3] - 411487*V[4]) + V[0]*(22658*V[0] - 208501*V[1] + 364863*V[2] - 288007*V[3] + 86329*V[4]) ;
	} else if(k == 6) {
		Beta[0] = 1152561*V[10]*V[10] + V[9]*(36480687*V[9] - 12950184*V[10]) + V[8]*(190757572*V[8] - 166461044*V[9] + 29442256*V[10]) + V[7]*(260445372*V[7] - 444003904*V[8] + 192596472*V[9] - 33918804*V[10]) + V[6]*(94851237*V[6] - 311771244*V[7] + 262901672*V[8] - 113206788*V[9] + 19834350*V[10]) + V[5]*(6150211*V[5] - 47460464*V[6] + 76206736*V[7] - 63394124*V[8] + 27060170*V[9] - 4712740*V[10]) ;
		Beta[1] = 271779*V[9]*V[9] + V[8]*(8449957*V[8] - 3015728*V[9]) + V[7]*(43093692*V[7] - 37913324*V[8] + 6694608*V[9]) + V[6]*(56662212*V[6] - 97838784*V[7] + 42405032*V[8] - 7408908*V[9]) + V[5]*(19365967*V[5] - 65224244*V[6] + 55053752*V[7] - 23510468*V[8] + 4067018*V[9]) + V[4]*(1152561*V[4] - 9117992*V[5] + 14742480*V[6] - 12183636*V[7] + 5134574*V[8] - 880548*V[9]) ;
		Beta[2] = 139633*V[8]*V[8] + V[7]*(3824847*V[7] - 1429976*V[8]) + V[6]*(17195652*V[6] - 15880404*V[7] + 2863984*V[8]) + V[5]*(19510972*V[5] - 35817664*V[6] + 15929912*V[7] - 2792660*V[8]) + V[4]*(5653317*V[4] - 20427884*V[5] + 17905032*V[6] - 7727988*V[7] + 1325006*V[8]) + V[3]*(271779*V[3] - 2380800*V[4] + 4086352*V[5] - 3462252*V[6] + 1458762*V[7] - 245620*V[8]) ;
		Beta[3] = 271779*V[7]*V[7] + V[6]*(5653317*V[6] - 2380800*V[7]) + V[5]*(19510972*V[5] - 20427884*V[6] + 4086352*V[7]) + V[4]*(17195652*V[4] - 35817664*V[5] + 17905032*V[6] - 3462252*V[7]) + V[3]*(3824847*V[3] - 15880404*V[4] + 15929912*V[5] - 7727988*V[6] + 1458762*V[7]) + V[2]*(139633*V[2] - 1429976*V[3] + 2863984*V[4] - 2792660*V[5] + 1325006*V[6] - 245620*V[7]) ;
		Beta[4] = 1152561*V[6]*V[6] + V[5]*(19365967*V[5] - 9117992*V[6]) + V[4]*(56662212*V[4] - 65224244*V[5] + 14742480*V[6]) + V[3]*(43093692*V[3] -97838784*V[4] + 55053752*V[5] - 12183636*V[6]) + V[2]*(8449957*V[2] - 37913324*V[3] + 42405032*V[4] - 23510468*V[5] + 5134574*V[6]) + V[1]*(271779*V[1] - 3015728*V[2] + 6694608*V[3] - 7408908*V[4] + 4067018*V[5] - 880548*V[6]) ;
		Beta[5] = 6150211*V[5]*V[5] + V[4]*(94851237*V[4] - 47460464*V[5]) + V[3]*(260445372*V[3] - 311771244*V[4] + 76206736*V[5]) + V[2]*(190757572*V[2] - 444003904*V[3] + 262901672*V[4] - 63394124*V[5]) + V[1]*(36480687*V[1] - 166461044*V[2] + 192596472*V[3] - 113206788*V[4] + 27060170*V[5]) + V[0]*(1152561*V[0] - 12950184*V[1] + 29442256*V[2] - 33918804*V[3] + 19834350*V[4] - 4712740*V[5]) ;
	} else if(k == 7) {
	}
}

// Gaussian quadrature weights for the WENO reconstruction in two dimensions.
void GetWENOGaussianQuadratureWeights(double**& GWMj, double**& GWPj, double*& dM, double*& dP, int k) {
	if(k == 1) {
	} else if(k == 2) {
	} else if(k == 3) {

		// weights for j - \Delta j/(2\sqrt(3))
		GWMj[0][0] = 1.0 + sqrt(3.0)/4.0	; GWMj[0][1] = -1.0/sqrt(3.0)  	; GWMj[0][2] = sqrt(3.0)/12.0 ;		//  j,   j+1 , j+2
		GWMj[1][0] = sqrt(3.0)/12.0 		; GWMj[1][1] = 1.0  			; GWMj[1][2] = -sqrt(3.0)/12.0 ;	// j-1,   j  , j+1
		GWMj[2][0] = -sqrt(3.0)/12.0  		; GWMj[2][1] = 1.0/sqrt(3.0)	; GWMj[2][2] = 1.0 - sqrt(3.0)/4.0 ;	// j-2 , j-1 , j
		dM[0] = ( 210.0 - sqrt(3.0) )/1080.0 	; dM[1] = 11.0/18.0 		; dM[2] = ( 210.0 + sqrt(3.0) )/1080.0 ;
//		dM[0] = ( 210.0 + sqrt(3.0) )/1080.0 	; dM[1] = 11.0/18.0 		; dM[2] = ( 210.0 - sqrt(3.0) )/1080.0 ;

		// weights for j + \Delta j/(2\sqrt(3))
		GWPj[0][0] = 1.0 - sqrt(3.0)/4.0 	; GWPj[0][1] = 1.0/sqrt(3.0)  	; GWPj[0][2] = -sqrt(3.0)/12.0 ;	//  j,   j+1 , j+2
		GWPj[1][0] = -sqrt(3.0)/12.0		; GWPj[1][1] = 1.0			 	; GWPj[1][2] = sqrt(3.0)/12.0 ;		// j-1,   j  , j+1
		GWPj[2][0] =  sqrt(3.0)/12.0		; GWPj[2][1] = -1.0/sqrt(3.0) 	; GWPj[2][2] = 1.0 + sqrt(3.0)/4.0 ;	// j-2 , j-1 , j
		dP[0] = ( 210.0 + sqrt(3.0) )/1080.0 	; dP[1] = 11.0/18.0 		; dP[2] = ( 210.0 - sqrt(3.0) )/1080.0 ;
//		dP[0] = ( 210.0 - sqrt(3.0) )/1080.0 	; dP[1] = 11.0/18.0 		; dP[2] = ( 210.0 + sqrt(3.0) )/1080.0 ;

	} else if(k == 4) {

		// weights for j - \Delta j/(2\sqrt(3))
		GWMj[0][0] = 1.0 + 43.0*sqrt(3.0)/144.0 + sqrt(3.0)/432.0 	; GWMj[0][1] = -69.0*sqrt(3.0)/144.0 - 3.0*sqrt(3.0)/432.0	; GWMj[0][2] = 33.0*sqrt(3.0)/144.0 + 3.0*sqrt(3.0)/432.0 	; GWMj[0][3] = -7.0*sqrt(3.0)/144.0 - sqrt(3.0)/432.0 ;		//  j, j+1, j+2, j+3
		GWMj[1][0] = 7.0*sqrt(3.0)/144.0 + sqrt(3.0)/432.0 		; GWMj[1][1] = 1.0 + 15.0*sqrt(3.0)/144.0 - 3.0*sqrt(3.0)/432.0	; GWMj[1][2] = -27.0*sqrt(3.0)/144.0 + 3.0*sqrt(3.0)/432.0 	; GWMj[1][3] = 5.0*sqrt(3.0)/144.0 - sqrt(3.0)/432.0 ; 		// j-1,	j, j+1, j+2
		GWMj[2][0] = -5.0*sqrt(3.0)/144.0 + sqrt(3.0)/432.0  		; GWMj[2][1] = 27.0*sqrt(3.0)/144.0 - 3.0*sqrt(3.0)/432.0 	; GWMj[2][2] = 1.0 - 15.0*sqrt(3.0)/144.0 + 3.0*sqrt(3.0)/432.0	; GWMj[2][3] = -7.0*sqrt(3.0)/144.0 - 1.0*sqrt(3.0)/432.0 ;	// j-2, j-1, j, j+1
		GWMj[3][0] = 7.0*sqrt(3.0)/144.0 + sqrt(3.0)/432.0  		; GWMj[3][1] = -33.0*sqrt(3.0)/144.0 - 3.0*sqrt(3.0)/432.0 	; GWMj[3][2] = 69.0*sqrt(3.0)/144.0 + 3.0*sqrt(3.0)/432.0	; GWMj[3][3] = 1.0 - 43.0*sqrt(3.0)/144.0 - sqrt(3.0)/432.0 ;	// j-3, j-2, j-1, j
		dM[0] = ( 3717.0 - 50.0*sqrt(3.0) )/166320.0 	; dM[1] = (72.0*sqrt(3.0)/7.0)*( (889.0*sqrt(3.0)/63360.0) - (587.0/1995840.0) ) ; dM[2] = (72.0*sqrt(3.0)/7.0)*( (889.0*sqrt(3.0)/63360.0) + (587.0/1995840.0) ) ; dM[3] = ( 3717.0 + 50.0*sqrt(3.0) )/166320.0 ;
	//	dM[0] = ( 3717.0 + 50.0*sqrt(3.0) )/166320.0 	; dM[1] = (72.0*sqrt(3.0)/7.0)*( (889.0*sqrt(3.0)/63360.0) + (587.0/1995840.0) ) ; dM[2] = (72.0*sqrt(3.0)/7.0)*( (889.0*sqrt(3.0)/63360.0) - (587.0/1995840.0) ) ; dM[3] = ( 3717.0 - 50.0*sqrt(3.0) )/166320.0 ;

		// weights for j + \Delta j/(2\sqrt(3))
		GWPj[0][0] = 1.0 - 43.0*sqrt(3.0)/144.0 - sqrt(3.0)/432.0 	; GWPj[0][1] = 69.0*sqrt(3.0)/144.0 + 3.0*sqrt(3.0)/432.0  	; GWPj[0][2] = -33.0*sqrt(3.0)/144.0 - 3.0*sqrt(3.0)/432.0 	; GWPj[0][3] = 7.0*sqrt(3.0)/144.0 + sqrt(3.0)/432.0 ;		//  j, j+1, j+2, j+3
		GWPj[1][0] = -7.0*sqrt(3.0)/144.0 - sqrt(3.0)/432.0		; GWPj[1][1] = 1.0 - 15.0*sqrt(3.0)/144.0 + 3.0*sqrt(3.0)/432.0 ; GWPj[1][2] = 27.0*sqrt(3.0)/144.0 - 3.0*sqrt(3.0)/432.0 	; GWPj[1][3] = -5.0*sqrt(3.0)/144.0 + sqrt(3.0)/432.0 ;		// j-1, j,  j+1, j+2
		GWPj[2][0] = 5.0*sqrt(3.0)/144.0 - sqrt(3.0)/432.0 		; GWPj[2][1] = -27.0*sqrt(3.0)/144.0 + 3.0*sqrt(3.0)/432.0 	; GWPj[2][2] = 1.0 + 15.0*sqrt(3.0)/144.0 - 3.0*sqrt(3.0)/432.0 ; GWPj[2][3] = 7.0*sqrt(3.0)/144.0 + 1.0*sqrt(3.0)/432.0 ;	// j-2 ,j-1, j,  j+1
		GWPj[3][0] = -7.0*sqrt(3.0)/144.0 - sqrt(3.0)/432.0 		; GWPj[3][1] = 33.0*sqrt(3.0)/144.0 + 3.0*sqrt(3.0)/432.0 	; GWPj[3][2] = -69.0*sqrt(3.0)/144.0 - 3.0*sqrt(3.0)/432.0 	; GWPj[3][3] =	1.0 + 43.0*sqrt(3.0)/144.0 + sqrt(3.0)/432.0 ;  // j-3 ,j-2, j-1, j
		dP[0] = ( 3717.0 + 50.0*sqrt(3.0) )/166320.0 	; dP[1] = (72.0*sqrt(3.0)/7.0)*( (889.0*sqrt(3.0)/63360.0) + (587.0/1995840.0) ) ; dP[2] = (72.0*sqrt(3.0)/7.0)*( (889.0*sqrt(3.0)/63360.0) - (587.0/1995840.0) ) ; dP[3] = ( 3717.0 - 50.0*sqrt(3.0) )/166320.0 ;
	//	dP[0] = ( 3717.0 - 50.0*sqrt(3.0) )/166320.0 	; dP[1] = (72.0*sqrt(3.0)/7.0)*( (889.0*sqrt(3.0)/63360.0) - (587.0/1995840.0) ) ; dP[2] = (72.0*sqrt(3.0)/7.0)*( (889.0*sqrt(3.0)/63360.0) + (587.0/1995840.0) ) ; dP[3] = ( 3717.0 + 50.0*sqrt(3.0) )/166320.0 ;

	} 
}

double DenominatorProd(int m, int k){
	double product = 1.0 ;
	for(int l = 0 ; l <= k ; l++) {
		if(l != m) product *= (m - l) ;
	}
	return product ;
}

double NumeratorProd(int m, int l, int r, int k){
	double product = 1.0 ;
	for(int q = 0 ; q <= k ; q++) {
		if((q != m) && (q != l)) product *= (r - q) ;
	}
	return product ;
}

void ComputeWENOCoefficientsUniform(MatrixPtr*& Crj, double*& d, int k) {
	int r,j,m,l;
	double *CH ; // coefficients of higher order optimal scheme.
	double numerator, denominator ;
	for(r = 0 ; r < k + 1 ; r++ ) {
		for(j = 0 ; j < k ; j++) {
			Crj[r][j] = 0.0 ;
			for(m = j + 1 ; m <= k ; m++) {
				denominator = DenominatorProd(m,k);
				numerator = 0.0 ;
				for(l = 0 ; l <= k ; l++) {
					if(l != m) numerator += NumeratorProd(m,l,r,k);
				}
				Crj[r][j] += numerator/denominator ;
			}
		}
	}
	// now computation of optimal weights. 
	// compute coefficients of the 2k-1 optimal scheme first
	CH = new double[2*k-1] ;
	r = k ;
	for(j = 0 ; j < 2*k-1 ; j++) {
		CH[j] = 0.0 ;
		for(m = j + 1 ; m <= 2*k-1 ; m++) {
			denominator = DenominatorProd(m,2*k-1);
			numerator = 0.0 ;
			for(l = 0 ; l <= 2*k-1 ; l++) {
				if(l != m) numerator += NumeratorProd(m,l,r,2*k-1);
			}
			CH[j] += numerator/denominator ;
		}
	}
	d[k-1] = CH[0]/Crj[k][0] ;
	for(j = k-2 ; j >= 0 ; j--) {
		d[j] = CH[k-1-j]/Crj[j+1][0] ;
		for(m = 1 ; m <= (k-1) - j ; m++) d[j] -= d[j+m]*Crj[j+m+1][m]/Crj[j+1][0] ;
	}
	delete [] CH ;
}

double IrrationalPlusNumeratorProd(int m, int l, int r, int k){
	double product = 1.0 ;
	for(int q = 0 ; q <= k ; q++) {
		if((q != m) && (q != l)) product *= (r - q + 0.5/sqrt(3.0) + 0.5 ) ;
	}
	return product ;
}

double IrrationalMinusNumeratorProd(int m, int l, int r, int k){
	double product = 1.0 ;
	for(int q = 0 ; q <= k ; q++) {
		if((q != m) && (q != l)) product *= (r - q - 0.5/sqrt(3.0) + 0.5 ) ;
	}
	return product ;
}

void ComputeWENOGaussianQuadratureWeights(MatrixPtr*& GWMj, MatrixPtr*& GWPj, double*& dM, double*& dP, int k) {
	int r,j,m,l;
	double *CHP, *CHM ; // coefficients of higher order optimal scheme.
	
	double Plus_numerator, Minus_numerator, denominator ;
	for(r = 0 ; r < k  ; r++ ) {
		for(j = 0 ; j < k ; j++) {
			GWMj[r][j] = GWPj[r][j] = 0.0 ;
			for(m = j + 1 ; m <= k ; m++) {
				denominator = DenominatorProd(m,k);
				Plus_numerator = Minus_numerator = 0.0 ;
				for(l = 0 ; l <= k ; l++) {
					if(l != m) {
						Plus_numerator += IrrationalPlusNumeratorProd(m,l,r,k);
						Minus_numerator += IrrationalMinusNumeratorProd(m,l,r,k);
					}
				}
				GWMj[r][j] += Minus_numerator/denominator ;
				GWPj[r][j] += Plus_numerator/denominator ;
			}
		}
	}
	// now computation of optimal weights. 
	// compute coefficients of the 2k-1 optimal scheme first
	CHP = new double[2*k-1] ; CHM = new double[2*k-1] ;
	r = k-1 ;
	for(j = 0 ; j < 2*k-1 ; j++) {
		CHM[j] = 0.0 ; CHP[j] = 0.0 ;
		for(m = j + 1 ; m <= 2*k-1 ; m++) {
			denominator = DenominatorProd(m,2*k-1);
			Plus_numerator = Minus_numerator = 0.0 ;
			for(l = 0 ; l <= 2*k-1 ; l++) {
				if(l != m) {
					Plus_numerator += IrrationalPlusNumeratorProd(m,l,r,2*k-1);
					Minus_numerator += IrrationalMinusNumeratorProd(m,l,r,2*k-1);
				}
			}
			CHM[j] += Minus_numerator/denominator ;
			CHP[j] += Plus_numerator/denominator ;
		} 
		std::cout << j << "\t" << CHM[j] << "\t" << CHP[j] << "\n" ;
	} 
	dM[k-1] = CHM[0]/GWMj[k-1][0] ; 
	dP[k-1] = CHP[0]/GWPj[k-1][0] ;
	for(j = k-2 ; j >= 0 ; j--) {
		dM[j] = CHM[k-1-j]/GWMj[j][0] ; dP[j] = CHP[k-1-j]/GWPj[j][0] ;
		for(m = 1 ; m <= (k-1) - j ; m++) {
			dM[j] -= dM[j+m]*GWMj[j+m][m]/GWMj[j][0] ;
			dP[j] -= dP[j+m]*GWPj[j+m][m]/GWPj[j][0] ;
		}
	} 
	delete [] CHM ; delete [] CHP ;
}

// Construct QR decomposition of matrix a which is n X n.
// Upper Triangular matrix R is returned in upper triangle of a, except for diagonals of R which are returned in d.
// sing returns 1 if singularity is encountered , else 0.
void QRdecompose(double** a, int n, double*& c, double*& d, bool& sing) {
        int i,j,k;
        double scale, sigma, sum, tau;

        c = new  double[n];
        d = new  double[n];
        sing = false;
        for(k = 0 ; k < n-1 ; k++) {
                scale = 0.0 ;
                for(i = k ; i < n ; i++) scale = Max(scale,fabs(a[i][k]));
                if (scale == 0.0) {
                        sing = true;
                        c[k] = d[k] = 0.0 ;
                }
                else {
                        for(i = k ; i < n ; i++) a[i][k] = a[i][k]/scale;
                        sum = 0.0 ;
                        for(i = k ; i < n ; i++) sum += a[i][k]*a[i][k];
                        sigma = SIGN(sqrt(sum), a[k][k]);
                        a[k][k] += sigma;
                        c[k] = sigma*a[k][k] ;
                        d[k] = -scale*sigma;
                        for(j = k+1 ; j < n ; j++) {
                                sum = 0.0 ;
                                for(i = k ; i < n ; i++) sum += a[i][k]*a[i][j];
                                tau = sum/c[k];
                                for(i = k ; i < n ; i++) a[i][j] -= tau*a[i][k];
                        }
                }
        }
        d[n-1] = a[n-1][n-1] ;
        if(fabs(d[n-1]) < 1.0E-12) sing = true;
}

// Solves Rx = B, R is upper triangular stored in a and d.
void rsolv(  double** a, int n,   double* d,   double*& b) {
        int i, j;
          double sum;

        b[n-1] = b[n-1]/d[n-1];
        for(i = n-2 ; i >= 0 ; i--) {
                sum = 0.0;
                for(j = i+1 ; j < n ; j++) sum += a[i][j]*b[j];
                b[i] = (b[i]-sum)/d[i];
        }
}

// Solves AX = B. a[][], c[], d[] are input from output of QRdecompose and are not modified.
// b is right hand side known vector and is overwritten by solution.
void qrsolv(  double** a, int n,   double* c,   double* d,   double*& b) {
        int i, j;
        double sum, tau;
        for(j = 0 ; j < n-1 ; j++) {
                sum = 0.0 ;
                for(i = j ; i < n ; i++) sum += a[i][j]*b[i];
                tau = sum/c[j];
                for(i = j ; i < n ; i++) b[i] -= tau*a[i][j];
        }
        rsolv(a,n,d,b);
}

// Find inverse of a matrix using QR factorization
void QRInverse(  double**& a_inv,   double** a, int n) {
        bool sing;
        int i, j;
        double *c, *d, *Solution, **Temp_Mat ;

        Solution = new   double[n] ;
        Allocate_2D_R(Temp_Mat,n,n);
        for(i = 0 ; i < n ; i++) {
                for(j = 0 ; j < n ; j++) Temp_Mat[i][j] = a[i][j] ;
        }

        QRdecompose(Temp_Mat,n,c,d,sing) ;
        if(sing) {
                std::cout << "Exiting due to singular matrix in QR Inverse function";
                exit(1);
        } else {
                for(i = 0 ; i < n ; i++) {
                        for(j = 0 ; j < n ; j++) {
                                if(j == i) Solution[j] = 1.0;
                                else Solution[j] = 0.0 ;
                        }
                        qrsolv(Temp_Mat,n,c,d,Solution) ;
                        for(j = 0 ; j < n ; j++) a_inv[j][i] = Solution[j] ;
                }
        }

        delete [] c ; delete [] d ;     delete [] Solution ;
        for(i = 0 ; i < n ; i++ ) delete [] Temp_Mat[i] ;
        delete [] Temp_Mat ;
}

// Solution of a linear system AX = B using Gaussian elimination
void Factor(  double** a, int n, int*& npivot,  double& det, bool& sing) {
        int i0;
        double* s;
        double* c;
        s = new   double[n];
        c = new   double[n];

        det = 1.0; sing = false ;
        for(int p=0;p<n;p++) {
                s[p] = fabs(a[p][0]);
                for(int q=1;q<n;q++) {
                        if(fabs(a[p][q])>s[p])
                                s[p] = fabs(a[p][q]);
                }
        }

        for(int k=0;k<(n-1);k++) {

                c[k] = fabs(a[k][k])/s[k];
                i0 = k;
                for(int i=k+1;i<n;i++) {
                        if((fabs(a[i][k])/s[i])>c[k]) {
                                i0 = i;
                                c[k] = (fabs(a[i][k])/s[i]);
                        }
                }
                npivot[k] = i0;
                if(fabs(c[k]) <= 1E-16) {
                        sing = true;
                /*        std::cout << "\n ERROR: Singular Matrix";
			std::cout << "\n" ;
			for(int i_ = 0 ; i_ < n ; i_++) {
				for(int j_ = 0 ; j_ < n ; j_++) std::cout << a[i_][j_] << "\t" ;
				std::cout << "\n";
			}
                 //       exit(1); */
                } else sing = false ;

                if(i0 != k) {
                        det = -det;
                        for(int j =k;j<n;j++){
                                a[k][j] = a[k][j] + a[i0][j];
                                a[i0][j] = a[k][j] - a[i0][j];
                                a[k][j] = a[k][j] - a[i0][j];
                        }
                }

                for(int l=k+1;l<n;l++) {
                        a[l][k] = a[l][k]/a[k][k];  // Multiplier mij
                        for(int m=k+1;m<n;m++) {
                                a[l][m] = a[l][m] - a[l][k]*a[k][m];
                        }
                }
                det = det * a[k][k];
        }
        det = det*a[n-1][n-1];
        delete[] s;
        delete[] c;
}

void Solve(  double** a, int n,   double*& b, int*& npivot) {
          double temp;
        int i0;

        for(int k=0;k<n-1;k++) {
                i0 = npivot[k];
                if(i0 != k) {
                        b[i0] = b[i0] + b[k];
                        b[k] = b[i0] - b[k];
                        b[i0] = b[i0] - b[k];
                }
                for(int i=k+1;i<n;i++)
                        b[i] = b[i] - a[i][k]*b[k];
        }

        b[n-1] = b[n-1]/a[n-1][n-1];
        for(int p=n-2;p>=0;p--) {
                temp = 0;
                for(int q=p+1;q<n;q++) {
                        temp += a[p][q]*b[q];
                }
                b[p] = (b[p]-temp)/a[p][p];
        }
}

// Find inverse of a matrix using LU decomposition
void LUInverse(  double**& a_inv,   double** a, int n) {
        bool sing;
        int i, j, *pivot;
        double *Solution, **Temp_Mat, det ;

        Solution = new   double[n] ; pivot = new int[n] ;
        Allocate_2D_R(Temp_Mat,n,n);
        for(i = 0 ; i < n ; i++) {
                for(j = 0 ; j < n ; j++) Temp_Mat[i][j] = a[i][j] ;
        }

        Factor(Temp_Mat,n,pivot,det,sing);
        if(sing) {
                std::cout << "Exiting due to singular matrix in QR Inverse function";
                exit(1);
        } else {
                for(i = 0 ; i < n ; i++) {
                        for(j = 0 ; j < n ; j++) {
                                if(j == i) Solution[j] = 1.0;
                                else Solution[j] = 0.0 ;
                        }
                        Solve(Temp_Mat,n,Solution,pivot) ;
                        for(j = 0 ; j < n ; j++) a_inv[j][i] = Solution[j] ;
                }
        }

        delete [] pivot ; delete [] Solution ;
        for(i = 0 ; i < n ; i++ ) delete [] Temp_Mat[i] ;
        delete [] Temp_Mat ;
}



void LUOptimizedDiagonalSystem(double*& u, int m, int N) {
	int i0, i, j1;
	i0 = m*(3*m+1)/2 + (2*m+1)*(N-2*m) ;
	for(int n = 1 ; n <= N ; n++) {
		if(n <= m) {
			// diagonal element should be non zero
			if(fabs(u[(n-1)*m + (n*(n+1))/2 -1]) < 1E-15) { std::cout << "Trouble: Singularity at "<< n ; exit(1); }
			for(int j = n+1 ; j < n+1+m ; j++) {
				if(j <= m) {
					u[(j-1)*m+j*(j+1)/2-(j-n)-1] /= u[(n-1)*m+(n*(n+1))/2-1] ;
					for(int k = 1 ; k < m+1 ; k++) u[(j-1)*m+j*(j+1)/2-(j-n)+k-1] -= u[(j-1)*m+j*(j+1)/2-(j-n)-1]*u[(n-1)*m+(n*(n+1))/2-1+k];
				} else if((j > m) && (j <= (N-m))) {
					u[m*(3*m+3)/2+(2*m+1)*(j-m-1)-(j-n)] /= u[(n-1)*m+(n*(n+1))/2-1] ;
					for(int k = 1 ; k < m+1 ; k++) u[m*(3*m+3)/2+(2*m+1)*(j-m-1)-(j-n)+k] -= u[m*(3*m+3)/2+(2*m+1)*(j-m-1)-(j-n)]*u[(n-1)*m+(n*(n+1))/2-1+k];
				} else {
					i = j - N + m;
					u[i0+2*m*(i-1)-(i-1)*(i-2)/2+m-(j-n)] /= u[(n-1)*m+(n*(n+1))/2-1] ;
					for(int k = 1 ; k < m+1 ; k++) u[i0 + 2*m*(i-1)-(i-1)*(i-2)/2+m-(j-n)+k] -= u[i0+2*m*(i-1)-(i-1)*(i-2)/2+m-(j-n)]*u[(n-1)*m+(n*(n+1))/2-1+k];
				}
			}
		} else if((n > m) && (n <= (N-m))) {
			// diagonal element should be non zero
			if( fabs( u[m*(3*m+3)/2+(2*m+1)*(n-m-1)]) < 1E-15 ) { std::cout << "Trouble: Singularity at " << n ; exit(1); }
			for(int j = n+1;j < n+1+m;j++) {
				if(j <= m) { std::cout << "Trouble" ; exit(1); } 
				else if((j > m) && (j <= (N-m))) {
					u[m*(3*m+3)/2+(2*m+1)*(j-m-1)-(j-n)] /= u[m*(3*m+3)/2+(2*m+1)*(n-m-1)] ;
					for(int k = 1 ; k < m+1 ; k++) u[m*(3*m+3)/2+(2*m+1)*(j-m-1)-(j-n)+k] -= u[m*(3*m+3)/2+(2*m+1)*(j-m-1)-(j-n)]*u[m*(3*m+3)/2 + (2*m+1)*(n-m-1)+k];
				} else {
					i = j - N + m;
					u[i0+2*m*(i-1)-(i-1)*(i-2)/2+m-(j-n)] /= u[m*(3*m+3)/2+(2*m+1)*(n-m-1)];
					for(int k = 1 ; k < m+1 ; k++) u[i0 + 2*m*(i-1)-(i-1)*(i-2)/2+m-(j-n)+k] -= u[i0+2*m*(i-1)-(i-1)*(i-2)/2+m-(j-n)]*u[m*(3*m+3)/2+(2*m+1)*(n-m-1)+k];
				}
			}
		} else {
			i = n + m - N ;
			// diagonal element should be non zero
			if(fabs( u[i0+2*m*(i-1)-(i-1)*(i-2)/2+m]) < 1E-15 ) { std::cout << "Trouble: Singularity at " << n ; exit(1); }
			for(int j = n+1;j < N+1;j++) {
				j1 = j - N + m;
				if(j <= N-m) { std::cout << "Trouble" ; exit(1); 
				} else {
					u[i0+2*m*(j1-1)-(j1-1)*(j1-2)/2+m-(j-n)] /= u[i0+2*m*(i-1)-(i-1)*(i-2)/2+m];
					for(int k = 1 ; k < (m+1)-i ; k++) u[i0+2*m*(j1-1)-(j1-1)*(j1-2)/2+m-(j-n)+k] -= u[i0+2*m*(j1-1)-(j1-1)*(j1-2)/2+m-(j-n)]*u[i0+2*m*(i-1)-(i-1)*(i-2)/2+m+k];
				}
			}
		}
	}
}

void SolveOptimizedDiagonalSystem(double*& Sol, double* a, int m, int N) {
	int i0, i, j1, n;
	i0 = m*(3*m+1)/2 + (2*m+1)*(N-2*m) ;
	for( n = 1 ; n <= N ; n++) {
		if(n <= m) {
			for(int j = n+1 ; j < n+1+m ; j++) {
				if(j <= m) Sol[j-1] -= Sol[n-1]*a[(j-1)*m+j*(j+1)/2-(j-n)-1] ;
				else if((j > m) && (j <= (N-m))) Sol[j-1] -= Sol[n-1]*a[m*(3*m+3)/2+(2*m+1)*(j-m-1)-(j-n)] ;
				else {
					i = j - N + m;
					Sol[j-1] -= Sol[n-1]*a[i0+2*m*(i-1)-(i-1)*(i-2)/2+m-(j-n)];
				}
			}
		} else if((n > m) && (n <= (N-m))) {
			for(int j = n+1;j < n+1+m;j++) {
				if(j <= m) { std::cout << "Trouble" ; exit(1); } 
				else if((j > m) && (j <= (N-m))) Sol[j-1] -= Sol[n-1]*a[m*(3*m+3)/2+(2*m+1)*(j-m-1)-(j-n)];
				else {
					i = j - N + m;
					Sol[j-1] -= Sol[n-1]*a[i0+2*m*(i-1)-(i-1)*(i-2)/2+m-(j-n)] ;
				}
			}
		} else {
			i = n + m - N ;
			for(int j = n+1;j < (N+1);j++) {
				j1 = j - N + m;
				if(j <= N-m) { std::cout << "Trouble" ; exit(1); } 
				else Sol[j-1] -= Sol[n-1]*a[i0+2*m*(j1-1)-(j1-1)*(j1-2)/2+m-(j-n)] ;
			}
		}
	}
	// Backsolve to get the solution, solution is stored in the right hand side vector Sol
	for( n = N ; n > 0 ; n--) {
		if(n > N-m) {
			i = n - N + m ;
			for(int j = 1 ; j < (m-i)+1 ; j++) {
				Sol[n-1] -= Sol[n-1+j]*a[i0 + 2*m*(i-1) - (i-1)*(i-2)/2+m+j];
			}
			Sol[n-1] /= a[i0 + 2*m*(i-1) - (i-1)*(i-2)/2+m];
		}
		else if((n > m) && (n <= N-m)) {
			for(int j = 1 ; j < m+1 ; j++) {
				Sol[n-1] -= Sol[n-1+j]*a[3*m*(m+1)/2+(2*m+1)*(n-m-1)+j];
			}
			Sol[n-1] /= a[3*m*(m+1)/2+(2*m+1)*(n-m-1)];
		}
		else {
			for(int j = 1 ; j < m+1 ; j++) {
				Sol[n-1] -= Sol[n-1+j]*a[(n-1)*m + n*(n+1)/2 - 1 + j] ;
			}
			Sol[n-1] /= a[(n-1)*m + n*(n+1)/2 - 1];
		}
	}
}

// LU Decomposition of a general diagonal system corresponding to compact schemes.
// note that in a domain of count points has (2n+1) and B has (6n + 4m +1)
// non zero elements placed symmetric wrt. diagonal
// RHS consisting of function values is passed in C and derivative is returned.
void SolveDiagonal(double**& A, double**& B, double*& C, double* u, int m, int n, int count) {
	int i,j,k;
	double factor;
	// First find right hand side and store it in C, then store the final solution also in C.
	for(i = 0; i < count ; i++) {
		if(i < 3*n + 2*m) {
			C[i] = 0.0 ;
			for(j = 0 ; j < (3*n + 2*m + 1) + i ; j++) {
				C[i] += B[i][3*n+2*m+j-i]*u[j];
			}
		}
		else if(i > count-1-(3*n+2*m)) {
			C[i] = 0.0 ;
			for(j = 0 ; j < (3*n + 2*m) + (count-i) ; j++) {
				C[i] += B[i][j]*u[i-(3*n+2*m)+j];
			}
		}
		else {
			C[i] = 0.0 ;
			for(j = 0 ; j < 6*n+4*m+1 ; j++) {
				C[i] += B[i][j]*u[i-(3*n+2*m)+j];
			}
		}
	}
	// This will suffice for tridiagonal.
	// Back substitution. only for a tridiagonal system.
	for(i = 0 ; i < count - n ; i++) {
		for(j = 1; j <= n ; j++) {
			factor = A[i+j][n-j]/A[i][n];
			for(k = 0 ; k <= n; k++) {
				A[i+j][n+k-j] = A[i+j][n+k-j] - factor*A[i][n+k];
			}
			C[i+j] = C[i+j] - factor*C[i];
		}
	}
	// For pentadiagonal or higher systems need to do more.
	for(i = count - n ; i < count-2 ; i++) {
		for(j = 1 ; j <= count - 1 - i ; j++) {
			// Figure this out later.
		}
	}
	C[count-1] = C[count-1]/A[count-1][n];
	for(i = count-2 ; i >= 0 ; i--) {
		C[i] = (C[i] - A[i][n+1]*C[i+1])/A[i][n];
	}
}

void PeriodicLUOptimizedDiagonalSystem(long double*& u, int n, int N) {
	int i0, iD, i, Jindex, j, k, Iindex;
	for(i = 0 ; i < N ; i++) {
		if(i < n) {
			// diagonal element should be non zero
			if(fabs(u[(2*n+1)*i + (i*(i+1))/2]) < 1E-15) { std::cout << "Trouble: Singularity at "<< i ; exit(1); }
			for(j = i+1 ; j <= i+n ; j++) {
				if(j < n) {
					u[(2*n+1)*j + (j*(j+1))/2-(j-i)] /= u[(2*n+1)*i + (i*(i+1))/2] ;
					for(k = 1 ; k < n+1 ; k++) u[(2*n+1)*j + (j*(j+1))/2-(j-i)+k] -= u[(2*n+1)*j + (j*(j+1))/2-(j-i)]*u[(2*n+1)*i + (i*(i+1))/2+k];
					for(k = 1 ; k < n+1 ; k++) u[(2*n+1)*(j+1) + (j*(j+1))/2-k] -= u[(2*n+1)*j + (j*(j+1))/2-(j-i)]*u[(2*n+1)*(i+1) + (i*(i+1))/2-k];
				} else if((j >= n) && (j <= (N-2*n-1))) {
					i0 = (2*n+1)*n + n*(n+1)/2;
					u[i0+(3*n+1)*(j-n)-(j-i)] /= u[(2*n+1)*i + (i*(i+1))/2] ;
					for(k = 1 ; k < n+1 ; k++) u[i0+(3*n+1)*(j-n)-(j-i)+k] -= u[i0+(3*n+1)*(j-n)-(j-i)]*u[(2*n+1)*i + (i*(i+1))/2 + k];
					for(k = 1 ; k < n+1 ; k++) u[i0+(3*n+1)*(j+1-n) -n - k] -= u[i0+(3*n+1)*(j-n)-(j-i)]*u[(2*n+1)*(i+1) + (i*(i+1))/2 - k];
				} else if((j > N-1-2*n) && (j <= N-1-n)) {  
					Jindex = j - N + 2*n + 1;
					i0 = (2*n+1)*n + n*(n+1)/2 + (3*n+1)*(N-3*n) ; 
					u[i0 + (3*n+1)*(Jindex-1)-Jindex*(Jindex-1)/2-(j-i)] /= u[(2*n+1)*i + (i*(i+1))/2] ;
					for(k = 1 ; k < n+1 ; k++) u[i0 + (3*n+1)*(Jindex-1)-Jindex*(Jindex-1)/2-(j-i)+k] -= u[i0 + (3*n+1)*(Jindex-1)-Jindex*(Jindex-1)/2-(j-i)]*u[(2*n+1)*i + (i*(i+1))/2+k];
					for(k = 1 ; k < n+1 ; k++) u[i0 + (3*n+1)*Jindex-Jindex*(Jindex+1)/2-n-k] -= u[i0 + (3*n+1)*(Jindex-1)-Jindex*(Jindex-1)/2-(j-i)]*u[(2*n+1)*(i+1) + (i*(i+1))/2 - k];
				} else {
					Jindex = j - N + n + 1;
					i0 = 2*n*n + (3*n+1)*(N-2*n) ;
					u[i0 + N*(Jindex-1) + j -(j-i)] /= u[(2*n+1)*i + (i*(i+1))/2] ;
					for(k = 1 ; k < n+1 ; k++) u[i0 + N*(Jindex-1)+j-(j-i)+k] -= u[i0 + N*(Jindex-1) + j -(j-i)]*u[(2*n+1)*i + (i*(i+1))/2+k];
					for(k = 1 ; k < n+1 ; k++) u[i0 + N*Jindex - k] -= u[i0+N*(Jindex-1)+j-(j-i)]*u[(2*n+1)*(i+1) + (i*(i+1))/2 - k];
				}
			}
			// reduce the last n rows.
			for(j = N-n ; j <= N-1 ; j++) {
				Jindex = j - N + n + 1;
				i0 = 2*n*n + (3*n+1)*(N-2*n) ;
				u[i0 + N*(Jindex-1) + j -(j-i)] /= u[(2*n+1)*i + (i*(i+1))/2] ;
				for(k = 1 ; k < n+1 ; k++) u[i0 + N*(Jindex-1)+j-(j-i)+k] -= u[i0 + N*(Jindex-1) + j -(j-i)]*u[(2*n+1)*i + (i*(i+1))/2+k];
				for(k = 1 ; k < n+1 ; k++) u[i0 + N*Jindex - k] -= u[i0+N*(Jindex-1)+j-(j-i)]*u[(2*n+1)*(i+1) + (i*(i+1))/2 - k];
			}
		} else if((i >= n) && (i <= (N-2*n-1))) {
			// diagonal element should be non zero
			iD = (2*n+1)*n + n*(n+1)/2;
			if( fabs( u[iD+(3*n+1)*(i-n)]) < 1E-15 ) { std::cout << "Trouble: Singularity at " << i ; exit(1); }
			for(j = i+1 ; j <= i+n ; j++) {
				if(j < n) { std::cout << "Trouble" ; exit(1); } 
				else if((j >= n) && (j <= (N-2*n-1))) {
					i0 = (2*n+1)*n + n*(n+1)/2;
					u[i0+(3*n+1)*(j-n)-(j-i)] /= u[iD+(3*n+1)*(i-n)] ;
					for(k = 1 ; k < n+1 ; k++) u[i0+(3*n+1)*(j-n)-(j-i)+k] -= u[i0+(3*n+1)*(j-n)-(j-i)]*u[iD+(3*n+1)*(i-n) + k];
					for(k = 1 ; k < n+1 ; k++) u[i0+(3*n+1)*(j+1-n)-n-k] -= u[i0+(3*n+1)*(j-n)-(j-i)]*u[iD+(3*n+1)*(i+1-n)-n-k];
				} else if((j > N-1-2*n) && (j <= N-1-n)) {  
					Jindex = j - N + 2*n + 1;
					i0 = (2*n+1)*n + n*(n+1)/2 + (3*n+1)*(N-3*n) ; 
					u[i0 + (3*n+1)*(Jindex-1)-Jindex*(Jindex-1)/2-(j-i)] /= u[iD+(3*n+1)*(i-n)] ;
					for(k = 1 ; k < n+1 ; k++) u[i0 + (3*n+1)*(Jindex-1)-Jindex*(Jindex-1)/2-(j-i)+k] -= u[i0 + (3*n+1)*(Jindex-1)-Jindex*(Jindex-1)/2-(j-i)]*u[iD+(3*n+1)*(i-n)+k];
					for(k = 1 ; k < n+1 ; k++) u[i0 + (3*n+1)*Jindex-Jindex*(Jindex+1)/2-n-k] -= u[i0 + (3*n+1)*(Jindex-1)-Jindex*(Jindex-1)/2-(j-i)]*u[iD+(3*n+1)*(i+1-n)-n-k];
				} else {
					Jindex = j - N + n + 1;
					i0 = 2*n*n + (3*n+1)*(N-2*n) ;
					u[i0 + N*(Jindex-1) + j -(j-i)] /= u[iD+(3*n+1)*(i-n)] ;
					for(k = 1 ; k < n+1 ; k++) u[i0 + N*(Jindex-1)+j-(j-i)+k] -= u[i0 + N*(Jindex-1) + j -(j-i)]*u[iD+(3*n+1)*(i-n)+k];
					for(k = 1 ; k < n+1 ; k++) u[i0 + N*Jindex - k] -= u[i0+N*(Jindex-1)+j-(j-i)]*u[iD+(3*n+1)*(i+1-n)-n-k];
				}
			}
			// reduce the last n rows.
			for(j = N-n ; j <= N-1 ; j++) {
				Jindex = j - N + n + 1;
				i0 = 2*n*n + (3*n+1)*(N-2*n) ;
				u[i0 + N*(Jindex-1) + j -(j-i)] /= u[iD+(3*n+1)*(i-n)] ;
				for(k = 1 ; k < n+1 ; k++) u[i0 + N*(Jindex-1)+j-(j-i)+k] -= u[i0 + N*(Jindex-1) + j -(j-i)]*u[iD+(3*n+1)*(i-n)+k];
				for(k = 1 ; k < n+1 ; k++) u[i0 + N*Jindex - k] -= u[i0+N*(Jindex-1)+j-(j-i)]*u[iD+(3*n+1)*(i+1-n)-n- k];
			}
		} else if( (i > N-1-2*n) && (i <= N-1-n) ) {
			// diagonal element should be non zero
			Iindex = i - N + 2*n + 1;
			iD = (2*n+1)*n + n*(n+1)/2 + (3*n+1)*(N-3*n);
			if( fabs( u[iD+(3*n+1)*(Iindex-1)-Iindex*(Iindex-1)/2]) < 1E-15 ) { std::cout << "Trouble: Singularity at " << i ; exit(1); }
			for(j = i+1 ; j <= N-n-1 ; j++) { // Note (N-n-1) instead of i+n
				if(j < n) { std::cout << "Trouble" ; exit(1); } 
				else if((j >= n) && (j <= (N-2*n-1))) { std::cout << "Trouble" ; exit(1); }
				else if((j > N-1-2*n) && (j <= N-1-n)) {  
					Jindex = j - N + 2*n + 1;
					i0 = (2*n+1)*n + n*(n+1)/2 + (3*n+1)*(N-3*n) ; 
					u[i0 + (3*n+1)*(Jindex-1)-Jindex*(Jindex-1)/2-(j-i)] /= u[iD+(3*n+1)*(Iindex-1)-Iindex*(Iindex-1)/2] ;
					for(k = 1 ; k < n+1 ; k++) u[i0 + (3*n+1)*(Jindex-1)-Jindex*(Jindex-1)/2-(j-i)+k] -= u[i0 + (3*n+1)*(Jindex-1)-Jindex*(Jindex-1)/2-(j-i)]*u[iD+(3*n+1)*(Iindex-1)-Iindex*(Iindex-1)/2+k];
					for(k = 1 ; k < n+1-Iindex ; k++) u[i0 + (3*n+1)*Jindex-Jindex*(Jindex+1)/2-n-k] -= u[i0 + (3*n+1)*(Jindex-1)-Jindex*(Jindex-1)/2-(j-i)]*u[iD+(3*n+1)*(Iindex)-Iindex*(Iindex+1)/2-n-k];
				} else { std::cout << "Trouble" ; exit(1); }
			}
			// reduce the last n rows.
			for(j = N-n ; j <= N-1 ; j++) {
				Jindex = j - N + n + 1;
				i0 = 2*n*n + (3*n+1)*(N-2*n) ;
				u[i0 + N*(Jindex-1) + j -(j-i)] /= u[iD+(3*n+1)*(Iindex-1)-Iindex*(Iindex-1)/2] ;
				for(k = 1 ; k < n+1 ; k++) u[i0 + N*(Jindex-1)+j-(j-i)+k] -= u[i0 + N*(Jindex-1) + j -(j-i)]*u[iD+(3*n+1)*(Iindex-1)-Iindex*(Iindex-1)/2+k];
				for(k = 1 ; k < n+1-Iindex ; k++) u[i0 + N*Jindex - k] -= u[i0+N*(Jindex-1)+j-(j-i)]*u[iD+(3*n+1)*(Iindex)-Iindex*(Iindex+1)/2-n- k];
			}
		} else {
			iD = 2*n*n + (3*n+1)*(N-2*n) ;
			Iindex = i - N + n + 1 ; 
			// diagonal element should be non zero
			if( fabs( u[iD+N*(Iindex-1)+i]) < 1E-15 ) { std::cout << "Trouble: Singularity at " << i ; exit(1); }
			for(j = i+1 ; j <= N-1 ; j++) {
				Jindex = j - N + n + 1 ;
				i0 = 2*n*n + (3*n+1)*(N-2*n) ;
				u[i0 + N*(Jindex-1) + j - (j-i)] /= u[iD+N*(Iindex-1)+i] ;
				// Reduce all the rows simultaneously.
				for(k = 1 ; k < N-i ; k++) u[i0 + N*(Jindex-1)+j-(j-i)+k] -= u[i0 + N*(Jindex-1) + j -(j-i)]*u[iD+N*(Iindex-1)+i+k];
			}
		}
	}
}

void PeriodicSolveOptimizedDiagonalSystem(long double* Sol, long double* u, int n, int N) {
	int i0, i, Jindex, j, Iindex;
	for(i = 0 ; i < N ; i++) {
		if(i < n) {
			for(j = i+1 ; j <= i+n ; j++) {
				if(j < n) { 
					Sol[j] -= Sol[i]*u[(2*n+1)*j + (j*(j+1))/2-(j-i)] ; 
				} else if((j >= n) && (j <= (N-2*n-1))) { 
					i0 = (2*n+1)*n + n*(n+1)/2 ; 
					Sol[j] -= Sol[i]*u[i0+(3*n+1)*(j-n)-(j-i)] ; 
				} else if((j > N-1-2*n) && (j <= N-1-n)) { 
					Jindex = j - N + 2*n + 1 ; 
					i0 = (2*n+1)*n + n*(n+1)/2 + (3*n+1)*(N-3*n) ; 
					Sol[j] -= Sol[i]*u[i0 + (3*n+1)*(Jindex-1)-Jindex*(Jindex-1)/2-(j-i)] ;
				} else { 
					Jindex = j - N + n + 1 ; 
					i0 = 2*n*n + (3*n+1)*(N-2*n) ; 
					Sol[j] -= Sol[i]*u[i0 + N*(Jindex-1) + j -(j-i)] ; 
				}
			}
			// reduce the last n rows.
			for(j = N-n ; j <= N-1 ; j++) {
				Jindex = j - N + n + 1;
				i0 = 2*n*n + (3*n+1)*(N-2*n) ;
				Sol[j] -= Sol[i]*u[i0 + N*(Jindex-1) + j -(j-i)] ;
			}
		} else if((i >= n) && (i <= (N-2*n-1))) {
			for(j = i+1 ; j <= i+n ; j++) {
				if(j < n) { std::cout << "Trouble" ; exit(1); } 
				else if((j >= n) && (j <= (N-2*n-1))) {
					i0 = (2*n+1)*n + n*(n+1)/2;
					Sol[j] -= Sol[i]*u[i0+(3*n+1)*(j-n)-(j-i)] ;
				} else if((j > N-1-2*n) && (j <= N-1-n)) {  
					Jindex = j - N + 2*n + 1;
					i0 = (2*n+1)*n + n*(n+1)/2 + (3*n+1)*(N-3*n) ; 
					Sol[j] -= Sol[i]*u[i0 + (3*n+1)*(Jindex-1)-Jindex*(Jindex-1)/2-(j-i)] ;
				} else {
					Jindex = j - N + n + 1;
					i0 = 2*n*n + (3*n+1)*(N-2*n) ;
					Sol[j] -= Sol[i]*u[i0 + N*(Jindex-1) + j -(j-i)] ;
				}
			}
			// reduce the last n rows.
			for(j = N-n ; j <= N-1 ; j++) {
				Jindex = j - N + n + 1;
				i0 = 2*n*n + (3*n+1)*(N-2*n) ;
				Sol[j] -= Sol[i]*u[i0 + N*(Jindex-1) + j -(j-i)] ;
			}
		} else if( (i > N-1-2*n) && (i <= N-1-n) ) {
			Iindex = i - N + 2*n + 1;
			for(j = i+1 ; j <= N-n-1 ; j++) { // Note (N-n-1) instead of i+n
				if(j < n) { std::cout << "Trouble" ; exit(1); } 
				else if((j >= n) && (j <= (N-2*n-1))) { std::cout << "Trouble" ; exit(1); }
				else if((j > N-1-2*n) && (j <= N-1-n)) {  
					Jindex = j - N + 2*n + 1;
					i0 = (2*n+1)*n + n*(n+1)/2 + (3*n+1)*(N-3*n) ; 
					Sol[j] -= Sol[i]*u[i0 + (3*n+1)*(Jindex-1)-Jindex*(Jindex-1)/2-(j-i)] ;
				} else { std::cout << "Trouble" ; exit(1); }
			}
			// reduce the last n rows.
			for(j = N-n ; j <= N-1 ; j++) {
				Jindex = j - N + n + 1;
				i0 = 2*n*n + (3*n+1)*(N-2*n) ;
				Sol[j] -= Sol[i]*u[i0 + N*(Jindex-1) + j -(j-i)] ;
			}
		} else {
			Iindex = i - N + n + 1 ; 
			for(j = i+1 ; j <= N-1 ; j++) {
				Jindex = j - N + n + 1 ;
				i0 = 2*n*n + (3*n+1)*(N-2*n) ;
				Sol[j] -= Sol[i]*u[i0 + N*(Jindex-1) + j - (j-i)] ;
			}
		}
	}

	for(i = N-1 ; i >= 0 ; i--) {
		if(i < n) {
			for(j = 1 ; j < n+1 ; j++) Sol[i] -= Sol[i+j]*u[(2*n+1)*i + i*(i+1)/2+j] ;
			for(j = n+1 ; j < 2*n+1 ; j++) 	Sol[i] -= Sol[j+(N-2*n-1)]*u[(2*n+1)*i + i*(i+1)/2+j] ;
			Sol[i] /= u[(2*n+1)*i + i*(i+1)/2] ;
		} else if( (i >= n) && (i <= N-1-2*n) ) {
			i0 = (2*n+1)*n + n*(n+1)/2 ;
			for(j = 1 ; j < n+1 ; j++) Sol[i] -= Sol[i+j]*u[i0+(3*n+1)*(i-n)+j] ;
			for(j = n+1 ; j < 2*n+1 ; j++) Sol[i] -= Sol[j+(N-2*n-1)]*u[i0+(3*n+1)*(i-n)+j] ;
			Sol[i] /= u[i0+(3*n+1)*(i-n)] ;
		} else if( (i > N-2*n-1) && ( i <= N-n-1) ) {
			Iindex = i - N + 2*n + 1 ;
			for(j = i+1 ; j <= N-1 ; j++) {
				i0 = (2*n+1)*n + n*(n+1)/2 + (3*n+1)*(N-3*n) ; 
				Sol[i] -= Sol[j]*u[i0 + (3*n+1)*(Iindex-1)-Iindex*(Iindex-1)/2 +j-i];
			}
			Sol[i] /= u[i0 + (3*n+1)*(Iindex-1)-Iindex*(Iindex-1)/2] ;
		} else {
			Iindex = i - N + n + 1 ;
			for(j = i+1 ; j <= N-1 ; j++) {
				i0 = 2*n*n + (3*n+1)*(N-2*n) ;
				Sol[i] -= Sol[j]*u[i0 + N*(Iindex-1)  + j] ;
			}
			i0 = 2*n*n + (3*n+1)*(N-2*n) ;
			Sol[i] /= u[i0 + N*(Iindex-1)+i];
		}	
	}
}