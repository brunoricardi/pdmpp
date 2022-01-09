#include "pwaves.h"

using namespace std;

pwaves::pwaves(parameters _parameters) : exprEngine(_parameters), pot(_parameters) {
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

    double dm;
    double it;
    for(int i = 0; i < this->params.nGrid; i++) {
        for(int j = 0; j <= i; j++) {
            // interaction does not depend on angular momentum
            it = this->pot.integrate(i*this->params.dr, j*this->params.dr);
            for(int k = 0; k < int(this->waves.size()); k++) {
                dm = it * this->exprEngine.freePW(k,i,j,this->params.tauMax);
                gsl_matrix_set(this->waves[k], i, j, dm);
            }
        }
    }
    // copy lower triangular block
    for(int i = 0; i < this->params.nGrid; i++) {
        for(int j = 0; j < i; j++) {
            for(int k = 0; k < int(this->waves.size()); k++) {
                gsl_matrix_set(this->waves[k], j, i, gsl_matrix_get(this->waves[k], i, j));
            }
        }
    }
    return 0;
}
