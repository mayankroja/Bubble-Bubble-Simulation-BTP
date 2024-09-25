#ifndef _Weighted_ENO_H_
#define _Weighted_ENO_H_

typedef double* MatrixPtr;

MatrixPtr *AllocMatR(int size1, int size2);

void Allocate_2D_R(double **&m, int d1, int d2);

void Allocate_3D_R(double ***&m, int d1, int d2, int d3);

void Allocate_4D_R(double ****&m, int d1, int d2, int d3, int d4);

void Allocate_2D_I(int **&m, int d1, int d2);

void Allocate_3D_I(int ***&m, int d1, int d2, int d3);

void Allocate_4D_I(int ****&m, int d1, int d2, int d3, int d4);

int ISign(int x);

double Max(double a, double b);

double SIGN(double a, double b);

// Sign function
int I_DSign(double x);

double Sign(double x);

double Minimum2(double x, double y);

double Maximum2(double x, double y);

double Minimum3(double x, double y, double z);

double Maximum3(double x, double y, double z);

double Maximum4(double a0, double a1, double a2, double a3);

double Maximum5(double a0, double a1, double a2, double a3, double a4);

double Maximum6(double a0, double a1, double a2, double a3, double a4, double a5);

double Minimum4(double a0, double a1, double a2, double a3);

double Minimum5(double a0, double a1, double a2, double a3, double a4);

double Minimum6(double a0, double a1, double a2, double a3, double a4, double a5);

// Minmod function with two variables
double MinMod2(double x, double y);

//Minmod function with three variables
double MinMod3(double x, double y, double z);

// Minmod function with four variables
double MinMod4(double a0, double a1, double a2, double a3);

double MinMod6(double a0, double a1, double a2, double a3, double a4, double a5);

// Median function
double Median(double x, double y, double z);

// Henricks mapping function for mapped WENO...
double MappedWENOFunction(double Omega, double C);

void MappedWeights(double *&Omega_, double *d_, int m, int flag);

// The notation etc. is followed from the NASA/CR-97-206253 Report by Chi-Wang Shu
// k denotes the order of accuracy. C_{rj} : r denotes the stencil and j is the summation index.
// The coefficients have been tabulated here for convenience, they can also be computed.
// Smoothness indicators should be computed in the code itself. Coefficients from Balsara and Shu.

void GetTabulatedWENOCoefficientsUniform(MatrixPtr *&Crj, double *&d, int k);

// size of Beta is k, and that of V is 2k-1.
void GetTabulatedWENOSmoothnessIndicatorUniform(double *&Beta, double *V, int k);

// Gaussian quadrature weights for the WENO reconstruction in two dimensions.
void GetWENOGaussianQuadratureWeights(double **&GWMj, double **&GWPj, double *&dM, double *&dP, int k);

double DenominatorProd(int m, int k);

double NumeratorProd(int m, int l, int r, int k);

void ComputeWENOCoefficientsUniform(MatrixPtr *&Crj, double *&d, int k);

double IrrationalPlusNumeratorProd(int m, int l, int r, int k);

double IrrationalMinusNumeratorProd(int m, int l, int r, int k);

void ComputeWENOGaussianQuadratureWeights(MatrixPtr *&GWMj, MatrixPtr *&GWPj, double *&dM, double *&dP, int k);

// Construct QR decomposition of matrix a which is n X n.
// Upper Triangular matrix R is returned in upper triangle of a, except for diagonals of R which are returned in d.
// sing returns 1 if singularity is encountered , else 0.
void QRdecompose(double **a, int n, double *&c, double *&d, bool &sing);

// Solves Rx = B, R is upper triangular stored in a and d.
void rsolv(double **a, int n, double *d, double *&b);

// Solves AX = B. a[][], c[], d[] are input from output of QRdecompose and are not modified.
// b is right hand side known vector and is overwritten by solution.
void qrsolv(double **a, int n, double *c, double *d, double *&b);

// Find inverse of a matrix using QR factorization
void QRInverse(double **&a_inv, double **a, int n);

// Solution of a linear system AX = B using Gaussian elimination
void Factor(double **a, int n, int *&npivot, double &det, bool &sing);

void Solve(double **a, int n, double *&b, int *&npivot);

// Find inverse of a matrix using LU decomposition
void LUInverse(double **&a_inv, double **a, int n);

void LUOptimizedDiagonalSystem(double *&u, int m, int N);

void SolveOptimizedDiagonalSystem(double *&Sol, double *a, int m, int N);

// LU Decomposition of a general diagonal system corresponding to compact schemes.
// note that in a domain of count points has (2n+1) and B has (6n + 4m +1)
// non zero elements placed symmetric wrt. diagonal
// RHS consisting of function values is passed in C and derivative is returned.
void SolveDiagonal(double **&A, double **&B, double *&C, double *u, int m, int n, int count);

void PeriodicLUOptimizedDiagonalSystem(long double *&u, int n, int N);

void PeriodicSolveOptimizedDiagonalSystem(long double *Sol, long double *u, int n, int N);

#endif