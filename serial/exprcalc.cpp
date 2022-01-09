#include "exprcalc.h"

using namespace std;

exprcalc::exprcalc(parameters _parameters) {
    this->params = _parameters;
}

exprcalc::~exprcalc() 
{
}

double exprcalc::freePW(int l, int i, int j, double tau) {
    double dm;
	double r12;
	double besvar;
	double arg_exp;
    double lambda;
    double r1, r2;
    
    lambda = this->params.lambda;
    r1 = i*this->params.dr;
    r2 = j*this->params.dr;
	r12 = r1*r2;
	besvar = r12 / (2. * lambda * tau);

	switch(this->params.nDim) {
		case 2 :
			if (besvar < 600.)
			{
				dm = gsl_sf_bessel_In(l, besvar);	// integer order I_n
				dm = dm * sqrt(r12) / (2. * lambda * tau); 
				dm = dm * exp( -(r1*r1 + r2*r2) / (4. * lambda * tau) );
			}
			else	// use asymptotic expression
			{
				arg_exp = besvar - ((r1*r1 + r2*r2) / (4. * lambda * tau));
				dm = exp(arg_exp);
				dm = dm * this->asympBesselPF(l, besvar); 
				dm = dm  / sqrt(4. * PI * lambda * tau);
					
			} 
			return dm;
		case 3 :
			if (besvar < 600.)
			{
				if (l > 50) {
					gsl_set_error_handler_off();			// avoid underflow here
					dm = gsl_sf_bessel_Inu(l+0.5, besvar); 	// fractional order I_{n+0.5}
					gsl_set_error_handler(NULL);
				}
				else {
					dm = gsl_sf_bessel_Inu(l+0.5, besvar); 	// fractional order I_{n+0.5}
				}
				dm = dm *  sqrt(PI / (2. * besvar));
				dm = dm * 4. * PI * r12 * pow(4. * PI * lambda * tau, -1.5); 
				dm = dm * exp( -(r1*r1 + r2*r2) / (4. * lambda * tau) );

			}
			else	// use asymptotic expression
			{
				arg_exp = besvar - (l*(l+1.0)/(2.0*besvar));
				arg_exp = arg_exp - ((r1*r1 + r2*r2) / (4. * lambda * tau));
				dm = 4. * PI * r12 * pow((4. * PI * lambda * tau),-1.5);
				dm = dm * exp(arg_exp);
				dm = dm * this->asympBesselPF(l, besvar);
			} 
			return dm;

		default:
			return -1;

	}
}

double exprcalc::asympBesselPF(double n, double x)
{
	double asbes;
	asbes = (1./(2.*x)) * (1. - n*(n+1.)/(pow(2.*x, 2.)) + n*(n+1.)*(n-2.)*(n+3.)/(3.*pow(2.*x, 3.)) + n*(n+1.)*(5.*n*n + 5.*n - 12.)/(2.*pow(2.*x,4.)));
	return asbes;
}