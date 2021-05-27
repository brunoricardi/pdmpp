#include <stdio.h>
#include <iostream>
#include <fstream>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <cmath>

using namespace std;


int main () {
	double lambda = 1.855;
	double tau;//=1.0/800.0;
	double dr = 0.1;
	double PI = 3.14159;
	int i;

	double r;
	double besvar;
	double y1, y2;
	double arg_exp;
	int n;
	double x;
	bool underflow_flag=false;
	double underflow_control;
	underflow_control = pow(10.0,-314.0);
	int status;

/*
	for (i=1; i<=50; i++) {
		r = i*dr;
		besvar = r*r / (2.0*lambda*tau);

		gsl_set_error_handler_off();
		y1 = gsl_sf_bessel_Inu(l + 0.5, besvar);
		gsl_set_error_handler(NULL);
		//cout << status << endl;
		//if (status == GSL_EUNDRFLW) {
		//	y1 = 0;
		//	cout << "underflow" << endl;
		//}
		//else {
		//	y1 = status;
		//}
		y1 = y1 * sqrt(PI / (2.0 * besvar));
		y1 = y1 * 4. * PI * r*r * pow(4. * PI * lambda * tau, -1.5);
		y1 = y1 * exp( -2.*r*r / (4. * lambda * tau) );
	
		arg_exp = besvar - ( l*(l+1.)/(2.*besvar) );
		arg_exp = arg_exp - (2.*r*r / (4. * lambda * tau));
		y2 = 4. * PI * r*r * pow(4. * PI * lambda * tau, -1.5);
		y2 = y2 * exp(arg_exp);
		n = l;
		x = besvar;
		y2 = y2 * (1./(2.*x)) * (1. - n*(n+1.)/(pow(2.*x, 2.)) + n*(n+1.)*(n-2.)*(n+3.)/(3.*pow(2.*x, 3.)) + n*(n+1.)*(5.*n*n + 5.*n - 12.)/(2.*pow(2.*x,4.)));

		cout << r << "\t" << besvar << "\t" << y1 << "\t" << y2 << "\t" << y1-y2 <<endl;
	}
*/

	int l,t;
	double pwc;
	pwc = pow(10.0,-5.0);
	r = 5.0;
	for(t=0; t<20; t++) {
		tau = 1.0 / pow(2.0, double(t));
		besvar = r*r / (2.0*lambda*tau);
		if (besvar < 600.0) {
			for(l=0; l<1000; l++) {
				gsl_set_error_handler_off();
				y1 = gsl_sf_bessel_Inu(l + 0.5, besvar);
				gsl_set_error_handler(NULL);
				y1 = y1 * sqrt(PI / (2.0 * besvar));
				y1 = y1 * 4. * PI * r*r * pow(4. * PI * lambda * tau, -1.5);
				y1 = y1 * exp( -2.*r*r / (4. * lambda * tau) );
				if (y1 < pwc) {
					cout << "Tau = " << tau << " stopped at l = " << l << endl;
					break;
				}
			}
		}
		else {
			for(l=0;l<1000;l++) {
				arg_exp = besvar - ( l*(l+1.)/(2.*besvar) );
				arg_exp = arg_exp - (2.*r*r / (4. * lambda * tau));
				y2 = 4. * PI * r*r * pow(4. * PI * lambda * tau, -1.5);
				y2 = y2 * exp(arg_exp);
				n = l;
				x = besvar;
				y2 = y2 * (1./(2.*x)) * (1. - n*(n+1.)/(pow(2.*x, 2.)) + n*(n+1.)*(n-2.)*(n+3.)/(3.*pow(2.*x, 3.)) + n*(n+1.)*(5.*n*n + 5.*n - 12.)/(2.*pow(2.*x,4.)));
				if (y2 < pwc) {
					cout << "Tau = " << tau << " stopped at l = " << l << endl;
					break;
				}
			}
		}
	}


	return 0;
}
