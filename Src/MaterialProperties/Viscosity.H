#ifndef _Viscosity_H_
#define _Viscosity_H_
#include <cmath>


/************************************************************************

Author: Mahesh K Sawardekar
Date: Sun Mar 21 21:25:55 2021

	* compute variable viscosity
	* created a class so different models can also be added in future

************************************************************************/

namespace mycode
{

class Viscosity
{
// public members
public:

// default constructor:
	Viscosity();

// destructor:
	~Viscosity();

	/**************************** functions *********************************/

	double MU() { return Mu_; }
	double MUMax() { return Mu_max; }
	double MUMin() { return Mu_min; }
	double Lambda_f() { return lambda_f; }
	void setMuMin(const double &m) { Mu_min = m; }
	void setMuMax(const double &m) { Mu_max = m; }
	double GetViscosity(const double &S11, const double &S12, const double &S22);
	double GetViscosity(const double &);
	double ComputeGammaDot(const double &S11, const double &S12, const double &S22);

        double GetViscosity(const double &S11, const double &S12, const double &S22, const double &phi);
	double GetViscosity(const double &S11, const double &S12, const double &S22, const double &phi, const double &dmg);
        double GetViscosity(const double &Gamma_Dot, const double &phi);
	double GetViscosityFromGammaDot(const double &Gamma_Dot, const double &phi, const double &dmg);
        double ShearModulus(const double &phi);
 
        /*********************Thermal conductivity********************************/
        double GetConductivity(double Mu,double Prandtl_no);
        


	// private members
private:
	/**************************** functions *********************************/


	/************************************************************************/


	/****************************   data    *********************************/
	
	double Mu_;
	double Mu_min, Mu_max;
	double lambda;
	double N;
        double eta = 0.0;

        double Mu1_;
        double Mu1_min, Mu1_max;
        double lambda1;
        double N1;
        double eta1 = 0.0;
        double lambda_f = std::sqrt(5.0);


	/************************************************************************/

}; /*End class Viscosity */

} /*End namespace mycode */

#endif
