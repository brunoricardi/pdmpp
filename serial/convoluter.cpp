#include "convoluter.h"

using namespace std;

Convoluter::Convoluter() {
}

Convoluter::~Convoluter() 
{
}

gsl_matrix* Convoluter::run(gsl_matrix *A, gsl_matrix*B) {
    int check;
    gsl_matrix *C;
    C = gsl_matrix_alloc(A->size1, A->size2);
    check = gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, A, B, 0.0, C);
    if (check == 0) {
        return C;
    }
    else {
        cout << "Error multiplying matrices in convoluter.run()" << endl;
        return C;
    }
}

PartWaves Convoluter::runOverWaves(PartWaves setA, PartWaves setB) {
    PartWaves setC(setA.parameters);
    setC.allocate();
    setC.parameters.tau = setA.parameters.tau + setB.parameters.tau;
    setC.parameters.temp = 1.0 / setC.parameters.tau;
    for (int i=0; i < int(setA.waves.size()); i++) {
        setC.waves[i] = this->run(setA.waves[i], setB.waves[i]);
    }
    return setC;
}

PartWaves Convoluter::runSquarer(PartWaves setA) {
    PartWaves setB(setA.parameters);
    setB.waves = setA.waves;
    cout << "Initialized Squarer with nSq = " << setA.parameters.nSq << " convolutions\n\n" << endl;
    cout << flush;
    for (int n=0; n < setA.parameters.nSq; n++) {
        setB = this->runOverWaves(setB, setB);
        cout << "Squaring n = " << n << " completed" << endl;
        cout << flush;
    }
    return setB;
}