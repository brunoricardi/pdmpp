#include "pwaves.h"

using namespace std;

PartWaves::PartWaves(Parameters _parameters) : exprCalc(_parameters), potential(_parameters) {
    this->parameters = _parameters;
}

PartWaves::~PartWaves() 
{
}

int PartWaves::allocate() {
    for(int i = 0; i < this->parameters.nw; i++) {
        this->waves.push_back(gsl_matrix_alloc(this->parameters.nGrid, this->parameters.nGrid));
    }
    return 0;
}

int PartWaves::initialize() {
    int check;
    check = this->allocate();
    if (check != 0) {
        cout << "Allocation of partial waves failed" << endl;
        return -1;
    }

    double dm;
    double it;
    for(int i = 0; i < this->parameters.nGrid; i++) {
        for(int j = 0; j <= i; j++) {
            // interaction does not depend on angular momentum
            it = this->potential.integrate(i*this->parameters.dr, j*this->parameters.dr);
            for(int k = 0; k < int(this->waves.size()); k++) {
                dm = it * this->exprCalc.freePW(k,i,j,this->parameters.tau);
                gsl_matrix_set(this->waves[k], i, j, dm);
            }
        }
    }
    // copy lower triangular block
    for(int i = 0; i < this->parameters.nGrid; i++) {
        for(int j = 0; j < i; j++) {
            for(int k = 0; k < int(this->waves.size()); k++) {
                gsl_matrix_set(this->waves[k], j, i, gsl_matrix_get(this->waves[k], i, j));
            }
        }
    }
    return 0;
}
