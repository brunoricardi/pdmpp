#include "pwaves.h"

using namespace std;

pwaves::pwaves(parameters _parameters) {
    this->params = _parameters;
}

pwaves::~pwaves() 
{
}

int pwaves::allocate() {
    for(int i = 0; i < this->params.nw; i++) {
        this->waves.push_back(gsl_matrix_alloc(this->params.nGrid, this->params.nGrid));
    }
    return 0;
}

int pwaves::initialize() {
    int check;
    check = this->allocate();
    if (check != 0) {
        cout << "Allocation of partial waves failed" << endl;
        return -1;
    }
    for(int i = 0; i < this->waves.size(); i++) {
        for(int j = 0; i < this->params.nGrid; j++) {
            for(int k = 0; k < j; k++) {
                gsl_matrix_set(this->waves[i], j, k, freepw(i,j,k))
            }
        }
    }
    return 0;
}

double pwaves::freepw(int l, int i, int j, double tau) {
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
				dm = dm * exp( -(r1*r1 + r2*r2)*this->nparams.dr / (4. * lambda * tau) );
			}
			else	// use asymptotic expression
			{
				arg_exp = besvar - ((r1*r1 + r2*r2) / (4. * lambda * tau));
				dm = exp(arg_exp);
				dm = dm * asympbessel_pf(l, besvar); 
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
				dm = dm * asympbessel_pf(l, besvar);
			} 
			return dm;

		default:
			return -1;

	}
}