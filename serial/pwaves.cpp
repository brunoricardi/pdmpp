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