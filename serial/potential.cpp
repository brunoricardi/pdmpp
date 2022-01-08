#include "potential.h"

using namespace std;

potential::potential(parameters _parameters) {
    this->params = _parameters;
}

potential::~potential() 
{
}

double potential::value(double r) {
	double pot;
	vector<double> pars(20);
	switch(this->params.potId) {
		case 0 :		// free particle
			pot = 0.;
			return pot;		


		case 1 :        // core-corona potential
			pars[0] = 1.0;	        // hard-core length scale
			pars[1] = 2.5;	        // soft-core length-scale
			pars[2] = 1000000.0; 	// hard-core energy (infite)
			pars[3] = 7.0;	        // soft-core energy
			if (r < pars[0])
				pot = pars[2];
			else if (r < pars[1])
				pot = pars[3];
			else
				pot = 0.0;
			return pot;

		case 2 : 		// Aziz1995 for He3
			{
				double a=186924.404;
				double alpha=10.5717543;
				double beta=-2.07758779;
				double d=1.438;
				double rm=0.29683;
				double eps=10.956;
				double c6=1.35186623;
				double c8=0.41495143;
				double c10=0.17151143;
				double sgmnm=0.2556;

				double rm2=pow(rm,2);
				double rm3=rm*rm2;
				double rm6=rm3*rm3;
				double rm8=rm6*rm2;
				double rm10=rm8*rm2;
				double epsa=eps*a;
				double ac=alpha/rm;
				double b2=-beta/rm2;
				double drm=d*rm;
				double epsc6=eps*c6*rm6;
				double epsc8=eps*c8*rm8;
				double epsc10=eps*c10*rm10;

				double dist;
				double oor,oor2,oor6;

				dist = r*sgmnm;
				oor = 1.0/dist;
				oor2 = pow(oor,2);
				oor6 = pow(oor2,3);
				pot = epsa * exp( (-ac - b2*dist) * dist );
				if (dist < drm)
					pot = pot - ( epsc6 + (epsc8+epsc10*oor2)*oor2 ) * oor6 * exp( -pow(drm*oor -1.0,2) );
				else
					pot = pot - ( epsc6 + (epsc8+epsc10*oor2)*oor2 ) * oor6;
				return pot;
				
			}

		case 3 :		// Lennard-Jones
			{
				double sig=2.556;
				double eps=10.22;
				double rinv, a;

				rinv = sig / r;
				a = pow(rinv, 6.0);
				pot = 4.0 * eps * a * (a - 1.0);
				return pot;
			}

		default :
			cout << "Invalid POTID!" << endl;
			return -1;

	}	
}