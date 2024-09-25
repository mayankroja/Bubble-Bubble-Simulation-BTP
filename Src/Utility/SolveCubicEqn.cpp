#include <cmath>
#include <complex>
#include<SolveCubicEqn.h>
        

void SolveCubicEqn(double &root1, double &root2, double &root3, double first_inv, double second_inv, double third_inv)
{
     double p;
     double q;
     double s;
     p = -1.0*first_inv;
     q = second_inv;
     s = -1.0*third_inv;

     double a = (1.0/3.0)*(3.0*q - p*p);
     double b = (-1.0/27.0)*(2*p*p*p - 9.0*p*q + 27.0*s);
     double Q = a/3.0;
     //double R = b/2.0;
     //double D = Q*Q*Q + R*R;
     //std::complex<double> S = std::pow( (R + std::sqrt(D)), 1.0/3.0 );
     //std::complex<double> T = std::pow( (R - std::sqrt(D)), 1.0/3.0 );
     //Now solve x^3 + ax + b = 0
     double rr = (-1.0/2.0)*b;
     double ii = std::sqrt(
                     (-1.0/27)*a*a*a
                    +(-1.0/4.0)*b*b
                     );
     std::complex<double> DD(rr, ii);
     std::complex<double> z = std::pow( DD, 1.0/3.0 );
     double RR = -1.0*z.real();
     double II = -1.0*z.imag();

     double root1p = 2.0*RR;
     double root2p = -1.0*RR + std::sqrt(3)*II;
     double root3p = -1.0*RR - std::sqrt(3)*II;

     root1 = root1p - p/3.0;
     root2 = root2p - p/3.0;
     root3 = root3p - p/3.0;
}
void SolveQuadraticEqn(double &root1, double &root2, double first_inv, double second_inv)
{
     double p = -1.0*first_inv;
     double q = second_inv;

     root1 = 0.5*(-1.0*p + std::sqrt(p*p - 4.0*q));
     root2 = 0.5*(-1.0*p - std::sqrt(p*p - 4.0*q));
}

