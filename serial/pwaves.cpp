#include "pwaves.h"

using namespace std;

pwaves::pwaves(parameters _parameters) : exprEngine(_parameters) {
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
    for(int i = 0; i < int(this->waves.size()); i++) {
        for(int j = 0; i < this->params.nGrid; j++) {
            for(int k = 0; k < j; k++) {
                gsl_matrix_set(this->waves[i], j, k, this->exprEngine.freePW(i,j,k,this->params.tauMax));
            }
        }
    }
    return 0;
}
