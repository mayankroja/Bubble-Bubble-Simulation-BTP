#include <AdvectLS_helper.H>

namespace mycode
{

void WENO5_LS
(
    double &Psix_L, double &Psix_R, 
    double &Psiy_L, double &Psiy_R, 
    int i, int j,
    const amrex::Real* dx,
    amrex::Array4<amrex::Real const> const &Psi)
{
    double Dxp_im2, Dxp_im, Dxp_i, Dxp_ip, Dxc_im2, Dxc_im, Dxc_i, Dxc_ip, Dxc_ip2;
    double Dyp_jm2, Dyp_jm, Dyp_j, Dyp_jp, Dyc_jm2, Dyc_jm, Dyc_j, Dyc_jp, Dyc_jp2;
    //	double Dzp_km2, Dzp_km, Dzp_k, Dzp_kp, Dzc_km2, Dzc_km, Dzc_k, Dzc_kp, Dzc_kp2 ;
    double w0, w2, a0, a1, a2, a, b, c, d, IS0, IS1, IS2;

    int k = 0;
    const double TOL = 1e-10;

     //  ------------------------   x  direction -------------------------------------------

    Dxp_im2 = (Psi(i - 1, j, k) - Psi(i - 2, j, k)) / dx[0];
    Dxp_im = (Psi(i, j, k) - Psi(i - 1, j, k)) / dx[0];
    Dxp_i = (Psi(i + 1, j, k) - Psi(i, j, k)) / dx[0];
    Dxp_ip = (Psi(i + 2, j, k) - Psi(i + 1, j, k)) / dx[0];

    Dxc_im2 = (Psi(i - 3, j, k) - 2.0 * Psi(i - 2, j, k) + Psi(i - 1, j, k)) / dx[0];
    Dxc_im = (Psi(i - 2, j, k) - 2.0 * Psi(i - 1, j, k) + Psi(i, j, k)) / dx[0];
    Dxc_i = (Psi(i - 1, j, k) - 2.0 * Psi(i, j, k) + Psi(i + 1, j, k)) / dx[0];
    Dxc_ip = (Psi(i, j, k) - 2.0 * Psi(i + 1, j, k) + Psi(i + 2, j, k)) / dx[0];
    Dxc_ip2 = (Psi(i + 1, j, k) - 2.0 * Psi(i + 2, j, k) + Psi(i + 3, j, k)) / dx[0];

    Psix_L = Psix_R = (1.0 / 12.0) * (-Dxp_im2 + 7.0 * Dxp_im + 7.0 * Dxp_i - Dxp_ip);

    a = Dxc_im2;
    b = Dxc_im;
    c = Dxc_i;
    d = Dxc_ip;

    IS0 = 13.0 * (a - b) * (a - b) + 3.0 * (a - 3.0 * b) * (a - 3.0 * b);
    IS1 = 13.0 * (b - c) * (b - c) + 3.0 * (b + c) * (b + c);
    IS2 = 13.0 * (c - d) * (c - d) + 3.0 * (3.0 * c - d) * (3.0 * c - d);

    a0 = 1.0 / ((TOL + IS0) * (TOL + IS0));
    a1 = 6.0 / ((TOL + IS1) * (TOL + IS1));
    a2 = 3.0 / ((TOL + IS2) * (TOL + IS2));
    w0 = a0 / (a0 + a1 + a2);
    w2 = a2 / (a0 + a1 + a2);

    Psix_L -= (w0 / 3.0) * (a - 2.0 * b + c) + ((w2 - 0.5) / 6.0) * (b - 2.0 * c + d);

    a = Dxc_ip2;
    b = Dxc_ip;
    c = Dxc_i;
    d = Dxc_im;

    IS0 = 13.0 * (a - b) * (a - b) + 3.0 * (a - 3.0 * b) * (a - 3.0 * b);
    IS1 = 13.0 * (b - c) * (b - c) + 3.0 * (b + c) * (b + c);
    IS2 = 13.0 * (c - d) * (c - d) + 3.0 * (3.0 * c - d) * (3.0 * c - d);

    a0 = 1.0 / ((TOL + IS0) * (TOL + IS0));
    a1 = 6.0 / ((TOL + IS1) * (TOL + IS1));
    a2 = 3.0 / ((TOL + IS2) * (TOL + IS2));
    w0 = a0 / (a0 + a1 + a2);
    w2 = a2 / (a0 + a1 + a2);

    Psix_R += (w0 / 3.0) * (a - 2.0 * b + c) + ((w2 - 0.5) / 6.0) * (b - 2.0 * c + d);

    //  ------------------------   y  direction -------------------------------------------

    Dyp_jm2 = (Psi(i, j - 1, k) - Psi(i, j - 2, k)) / dx[1];
    Dyp_jm = (Psi(i, j, k) - Psi(i, j - 1, k)) / dx[1];
    Dyp_j = (Psi(i, j + 1, k) - Psi(i, j, k)) / dx[1];
    Dyp_jp = (Psi(i, j + 2, k) - Psi(i, j + 1, k)) / dx[1];

    Dyc_jm2 = (Psi(i, j - 3, k) - 2.0 * Psi(i, j - 2, k) + Psi(i, j - 1, k)) / dx[1];
    Dyc_jm = (Psi(i, j - 2, k) - 2.0 * Psi(i, j - 1, k) + Psi(i, j, k)) / dx[1];
    Dyc_j = (Psi(i, j - 1, k) - 2.0 * Psi(i, j, k) + Psi(i, j + 1, k)) / dx[1];
    Dyc_jp = (Psi(i, j, k) - 2.0 * Psi(i, j + 1, k) + Psi(i, j + 2, k)) / dx[1];
    Dyc_jp2 = (Psi(i, j + 1, k) - 2.0 * Psi(i, j + 2, k) + Psi(i, j + 3, k)) / dx[1];

    Psiy_L = Psiy_R = (1.0 / 12.0) * (-Dyp_jm2 + 7.0 * Dyp_jm + 7.0 * Dyp_j - Dyp_jp);

    a = Dyc_jm2;
    b = Dyc_jm;
    c = Dyc_j;
    d = Dyc_jp;

    IS0 = 13.0 * (a - b) * (a - b) + 3.0 * (a - 3.0 * b) * (a - 3.0 * b);
    IS1 = 13.0 * (b - c) * (b - c) + 3.0 * (b + c) * (b + c);
    IS2 = 13.0 * (c - d) * (c - d) + 3.0 * (3.0 * c - d) * (3.0 * c - d);

    a0 = 1.0 / ((TOL + IS0) * (TOL + IS0));
    a1 = 6.0 / ((TOL + IS1) * (TOL + IS1));
    a2 = 3.0 / ((TOL + IS2) * (TOL + IS2));
    w0 = a0 / (a0 + a1 + a2);
    w2 = a2 / (a0 + a1 + a2);

    Psiy_L -= (w0 / 3.0) * (a - 2.0 * b + c) + ((w2 - 0.5) / 6.0) * (b - 2.0 * c + d);

    a = Dyc_jp2;
    b = Dyc_jp;
    c = Dyc_j;
    d = Dyc_jm;

    IS0 = 13.0 * (a - b) * (a - b) + 3.0 * (a - 3.0 * b) * (a - 3.0 * b);
    IS1 = 13.0 * (b - c) * (b - c) + 3.0 * (b + c) * (b + c);
    IS2 = 13.0 * (c - d) * (c - d) + 3.0 * (3.0 * c - d) * (3.0 * c - d);

    a0 = 1.0 / ((TOL + IS0) * (TOL + IS0));
    a1 = 6.0 / ((TOL + IS1) * (TOL + IS1));
    a2 = 3.0 / ((TOL + IS2) * (TOL + IS2));
    w0 = a0 / (a0 + a1 + a2);
    w2 = a2 / (a0 + a1 + a2);

    Psiy_R += (w0 / 3.0) * (a - 2.0 * b + c) + ((w2 - 0.5) / 6.0) * (b - 2.0 * c + d);

    /*	//  ------------------------   z  direction -------------------------------------------	

Dzp_km2 = (Psi[i][j][k-1][p] - Psi[i][j][k-2][p])/deltaz ;
Dzp_km = (Psi[i][j][k][p] - Psi[i][j][k-1][p])/deltaz ;
Dzp_k = (Psi[i][j][k+1][p] - Psi[i][j][k][p])/deltaz ; 
Dzp_kp = (Psi[i][j][k+2][p] - Psi[i][j][k+1][p])/deltaz ; 

Dzc_km2 = (Psi[i][j][k-3][p] - 2.0*Psi[i][j][k-2][p] + Psi[i][j][k-1][p])/deltaz ;
Dzc_km = (Psi[i][j][k-2][p] - 2.0*Psi[i][j][k-1][p] + Psi[i][j][k][p])/deltaz ;
Dzc_k = (Psi[i][j][k-1][p] - 2.0*Psi[i][j][k][p] + Psi[i][j][k+1][p])/deltaz ;
Dzc_kp = (Psi[i][j][k][p] - 2.0*Psi[i][j][k+1][p] + Psi[i][j][k+2][p])/deltaz ;
Dzc_kp2 = (Psi[i][j][k+1][p] - 2.0*Psi[i][j][k+2][p] + Psi[i][j][k+3][p])/deltaz ;

Psiz_L = Psiz_R = (1.0/12.0)*( -Dzp_km2 + 7.0*Dzp_km + 7.0*Dzp_k - Dzp_kp) ;

a = Dzc_km2 ; b = Dzc_km ; c = Dzc_k ; d = Dzc_kp ;

IS0 = 13.0*(a - b)*(a - b) + 3.0*(a - 3.0*b)*(a - 3.0*b) ;
IS1 = 13.0*(b - c)*(b - c) + 3.0*(b + c)*(b + c) ;
IS2 = 13.0*(c - d)*(c - d) + 3.0*(3.0*c - d)*(3.0*c - d) ;

a0 = 1.0/( (TOL + IS0)*(TOL + IS0) ) ; a1 = 6.0/( (TOL + IS1)*(TOL + IS1) ) ; a2 = 3.0/( (TOL + IS2)*(TOL + IS2) ) ;
w0 = a0/(a0 + a1 + a2) ; w2 = a2/(a0 + a1 + a2) ;

Psiz_L -= (w0/3.0)*(a - 2.0*b + c) + ( (w2 - 0.5)/6.0 )*(b - 2.0*c + d) ;

a = Dzc_kp2 ; b = Dzc_kp ; c = Dzc_k ; d = Dzc_km ;

IS0 = 13.0*(a - b)*(a - b) + 3.0*(a - 3.0*b)*(a - 3.0*b) ;
IS1 = 13.0*(b - c)*(b - c) + 3.0*(b + c)*(b + c) ;
IS2 = 13.0*(c - d)*(c - d) + 3.0*(3.0*c - d)*(3.0*c - d) ;

a0 = 1.0/( (TOL + IS0)*(TOL + IS0) ) ; a1 = 6.0/( (TOL + IS1)*(TOL + IS1) ) ; a2 = 3.0/( (TOL + IS2)*(TOL + IS2) ) ;
w0 = a0/(a0 + a1 + a2) ; w2 = a2/(a0 + a1 + a2) ;

Psiz_R += (w0/3.0)*(a - 2.0*b + c) + ( (w2 - 0.5)/6.0 )*(b - 2.0*c + d) ; */
}

void WENO5_Fij
(
    double &Psix_L, double &Psix_R, 
    double &Psiy_L, double &Psiy_R, 
    int i, int j,
    const amrex::Real* dx,
    amrex::Array4<amrex::Real const> const &Psi)
{
    double Dxp_im2, Dxp_im, Dxp_i, Dxp_ip, Dxc_im2, Dxc_im, Dxc_i, Dxc_ip, Dxc_ip2;
    double Dyp_jm2, Dyp_jm, Dyp_j, Dyp_jp, Dyc_jm2, Dyc_jm, Dyc_j, Dyc_jp, Dyc_jp2;
    //	double Dzp_km2, Dzp_km, Dzp_k, Dzp_kp, Dzc_km2, Dzc_km, Dzc_k, Dzc_kp, Dzc_kp2 ;
    double w0, w1, w2, a0, a1, a2, a, b, c, d, IS0, IS1, IS2;
    double q0, q1, q2, q3;
    double q0p, q1p, q2p;
    double dm0,dm1,dm2;

    double dphim2,dphim1,dphi0,dphip1,dphip2;

    int k = 0;
    const double TOL = 1e-10;

     //  ------------------------   x  direction -------------------------------------------

    q0 = (1.0/3.0)*Psi(i - 2, j, k) - (7.0/6.0)*Psi(i - 1, j, k) + (11.0/6.0)*Psi(i, j, k);
    q1 = (-1.0/6.0)*Psi(i - 1, j, k) + (5.0/6.0)*Psi(i, j, k) + (1.0/3.0)*Psi(i + 1, j, k);
    q2 = (1.0/3.0)*Psi(i , j, k) + (5.0/6.0)*Psi(i + 1, j, k) - (1.0/6.0)*Psi(i + 2, j, k);
    q3 = (1.0/3.0)*Psi(i + 3, j, k) - (7.0/6.0)*Psi(i + 2, j, k) + (11.0/6.0)*Psi(i + 1, j, k);

    dphim2 = Psi(i - 2, j, k) - Psi(i - 1, j, k);
    dphim1 = Psi(i - 1, j, k) - Psi(i, j, k);
    dphi0 = Psi(i, j, k) - Psi(i + 1, j, k);
    dphip1 = Psi(i + 1, j, k) - Psi(i + 2, j, k);
    dphip2 = Psi(i + 2, j, k) - Psi(i + 3, j, k);

    dm0 = 1.0/10.0;
    dm1 = 6.0/10.0;
    dm2 = 3.0/10.0;

    a = dphim2;
    b = dphim1;
    c = dphi0;
    d = dphip1;

    IS0 = (13.0/12.0) * (a - b) * (a - b) + (1.0/4.0) * (a - 3.0 * b) * (a - 3.0 * b);
    IS1 = (13.0/12.0) * (b - c) * (b - c) + (1.0/4.0) * (b + c) * (b + c);
    IS2 = (13.0/12.0) * (c - d) * (c - d) + (1.0/4.0) * (3.0 * c - d) * (3.0 * c - d);

    a0 = 1.0 / (10.0 * (TOL + IS0) * (TOL + IS0));
    a1 = 6.0 / (10.0 * (TOL + IS1) * (TOL + IS1));
    a2 = 3.0 / (10.0 * (TOL + IS2) * (TOL + IS2));
    w0 = a0 / (a0 + a1 + a2);
    w1 = a1 / (a0 + a1 + a2);
    w2 = a2 / (a0 + a1 + a2);

    Psix_L = w0*q0 + w1*q1 + w2*q2;


    a = dphip2;
    b = dphip1;
    c = dphi0;
    d = dphim1;

    IS0 = (13.0/12.0) * (a - b) * (a - b) + (1.0/4.0) * (a - 3.0 * b) * (a - 3.0 * b);
    IS1 = (13.0/12.0) * (b - c) * (b - c) + (1.0/4.0) * (b + c) * (b + c);
    IS2 = (13.0/12.0) * (c - d) * (c - d) + (1.0/4.0) * (3.0 * c - d) * (3.0 * c - d);

    a0 = 1.0 / (10.0 * (TOL + IS0) * (TOL + IS0));
    a1 = 6.0 / (10.0 * (TOL + IS1) * (TOL + IS1));
    a2 = 3.0 / (10.0 * (TOL + IS2) * (TOL + IS2));
    w0 = a0 / (a0 + a1 + a2);
    w1 = a1 / (a0 + a1 + a2);
    w2 = a2 / (a0 + a1 + a2);

    Psix_R = w0*q3 + w1*q2 + w2*q1;


    //  ------------------------   y  direction -------------------------------------------

    q0 = (1.0/3.0)*Psi(i, j - 2, k) - (7.0/6.0)*Psi(i, j - 1, k) + (11.0/6.0)*Psi(i, j, k);
    q1 = (-1.0/6.0)*Psi(i, j - 1, k) + (5.0/6.0)*Psi(i, j, k) + (1.0/3.0)*Psi(i, j + 1, k);
    q2 = (1.0/3.0)*Psi(i , j, k) + (5.0/6.0)*Psi(i, j + 1, k) - (1.0/6.0)*Psi(i, j + 2, k);
    q3 = (1.0/3.0)*Psi(i, j + 3, k) - (7.0/6.0)*Psi(i, j + 2, k) + (11.0/6.0)*Psi(i, j + 1, k);

    dphim2 = Psi(i, j - 2, k) - Psi(i, j - 1, k);
    dphim1 = Psi(i, j - 1, k) - Psi(i, j, k);
    dphi0 = Psi(i, j, k) - Psi(i, j + 1, k);
    dphip1 = Psi(i, j + 1, k) - Psi(i, j + 2, k);
    dphip2 = Psi(i, j + 2, k) - Psi(i, j + 3, k);

    dm0 = 1.0/10.0;
    dm1 = 6.0/10.0;
    dm2 = 3.0/10.0;

    a = dphim2;
    b = dphim1;
    c = dphi0;
    d = dphip1;

    IS0 = (13.0/12.0) * (a - b) * (a - b) + (1.0/4.0) * (a - 3.0 * b) * (a - 3.0 * b);
    IS1 = (13.0/12.0) * (b - c) * (b - c) + (1.0/4.0) * (b + c) * (b + c);
    IS2 = (13.0/12.0) * (c - d) * (c - d) + (1.0/4.0) * (3.0 * c - d) * (3.0 * c - d);

    a0 = 1.0 / (10.0 * (TOL + IS0) * (TOL + IS0));
    a1 = 6.0 / (10.0 * (TOL + IS1) * (TOL + IS1));
    a2 = 3.0 / (10.0 * (TOL + IS2) * (TOL + IS2));
    w0 = a0 / (a0 + a1 + a2);
    w1 = a1 / (a0 + a1 + a2);
    w2 = a2 / (a0 + a1 + a2);

    Psiy_L = w0*q0 + w1*q1 + w2*q2;


    a = dphip2;
    b = dphip1;
    c = dphi0;
    d = dphim1;

    IS0 = (13.0/12.0) * (a - b) * (a - b) + (1.0/4.0) * (a - 3.0 * b) * (a - 3.0 * b);
    IS1 = (13.0/12.0) * (b - c) * (b - c) + (1.0/4.0) * (b + c) * (b + c);
    IS2 = (13.0/12.0) * (c - d) * (c - d) + (1.0/4.0) * (3.0 * c - d) * (3.0 * c - d);

    a0 = 1.0 / (10.0 * (TOL + IS0) * (TOL + IS0));
    a1 = 6.0 / (10.0 * (TOL + IS1) * (TOL + IS1));
    a2 = 3.0 / (10.0 * (TOL + IS2) * (TOL + IS2));
    w0 = a0 / (a0 + a1 + a2);
    w1 = a1 / (a0 + a1 + a2);
    w2 = a2 / (a0 + a1 + a2);

    Psiy_R = w0*q3 + w1*q2 + w2*q1;
    /*	//  ------------------------   z  direction -------------------------------------------	

Dzp_km2 = (Psi[i][j][k-1][p] - Psi[i][j][k-2][p])/deltaz ;
Dzp_km = (Psi[i][j][k][p] - Psi[i][j][k-1][p])/deltaz ;
Dzp_k = (Psi[i][j][k+1][p] - Psi[i][j][k][p])/deltaz ; 
Dzp_kp = (Psi[i][j][k+2][p] - Psi[i][j][k+1][p])/deltaz ; 

Dzc_km2 = (Psi[i][j][k-3][p] - 2.0*Psi[i][j][k-2][p] + Psi[i][j][k-1][p])/deltaz ;
Dzc_km = (Psi[i][j][k-2][p] - 2.0*Psi[i][j][k-1][p] + Psi[i][j][k][p])/deltaz ;
Dzc_k = (Psi[i][j][k-1][p] - 2.0*Psi[i][j][k][p] + Psi[i][j][k+1][p])/deltaz ;
Dzc_kp = (Psi[i][j][k][p] - 2.0*Psi[i][j][k+1][p] + Psi[i][j][k+2][p])/deltaz ;
Dzc_kp2 = (Psi[i][j][k+1][p] - 2.0*Psi[i][j][k+2][p] + Psi[i][j][k+3][p])/deltaz ;

Psiz_L = Psiz_R = (1.0/12.0)*( -Dzp_km2 + 7.0*Dzp_km + 7.0*Dzp_k - Dzp_kp) ;

a = Dzc_km2 ; b = Dzc_km ; c = Dzc_k ; d = Dzc_kp ;

IS0 = 13.0*(a - b)*(a - b) + 3.0*(a - 3.0*b)*(a - 3.0*b) ;
IS1 = 13.0*(b - c)*(b - c) + 3.0*(b + c)*(b + c) ;
IS2 = 13.0*(c - d)*(c - d) + 3.0*(3.0*c - d)*(3.0*c - d) ;

a0 = 1.0/( (TOL + IS0)*(TOL + IS0) ) ; a1 = 6.0/( (TOL + IS1)*(TOL + IS1) ) ; a2 = 3.0/( (TOL + IS2)*(TOL + IS2) ) ;
w0 = a0/(a0 + a1 + a2) ; w2 = a2/(a0 + a1 + a2) ;

Psiz_L -= (w0/3.0)*(a - 2.0*b + c) + ( (w2 - 0.5)/6.0 )*(b - 2.0*c + d) ;

a = Dzc_kp2 ; b = Dzc_kp ; c = Dzc_k ; d = Dzc_km ;

IS0 = 13.0*(a - b)*(a - b) + 3.0*(a - 3.0*b)*(a - 3.0*b) ;
IS1 = 13.0*(b - c)*(b - c) + 3.0*(b + c)*(b + c) ;
IS2 = 13.0*(c - d)*(c - d) + 3.0*(3.0*c - d)*(3.0*c - d) ;

a0 = 1.0/( (TOL + IS0)*(TOL + IS0) ) ; a1 = 6.0/( (TOL + IS1)*(TOL + IS1) ) ; a2 = 3.0/( (TOL + IS2)*(TOL + IS2) ) ;
w0 = a0/(a0 + a1 + a2) ; w2 = a2/(a0 + a1 + a2) ;

Psiz_R += (w0/3.0)*(a - 2.0*b + c) + ( (w2 - 0.5)/6.0 )*(b - 2.0*c + d) ; */
}

} // namespace mycode
