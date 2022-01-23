#include "convoluter.h"

using namespace std;

Convoluter::Convoluter() {
}

Convoluter::~Convoluter() 
{
}

gsl_matrix* Convoluter::run(gsl_matrix *A, gsl_matrix*B, double dr) {
    int check;
    gsl_matrix *C;
    C = gsl_matrix_alloc(A->size1, A->size2);
    check = gsl_blas_dsymm(CblasLeft, CblasUpper, dr, A, B, 0.0, C);
    if (check == 0) {
        return C;
    }
    else {
        cout << "Error multiplying matrices in convoluter.run()" << endl;
        return C;
    }
}

// YOU ARE FORGETTING TO MULTIPLY THE MATRICES BY DR

PartWaves Convoluter::runOverWaves(PartWaves setA, PartWaves setB) {
    PartWaves setC(setA.parameters);
    setC.allocate();
    setC.parameters.tau = setA.parameters.tau + setB.parameters.tau;
    setC.parameters.temp = 1.0 / setC.parameters.tau;
    for (int i=0; i < int(setA.waves.size()); i++) {
        setC.waves[i] = this->run(setA.waves[i], setB.waves[i], setC.parameters.dr);
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
    setB = this->truncCorrector(setB);
    return setB;
}


PartWaves Convoluter::truncCorrector(PartWaves setA) {
    PartWaves setB(setA.parameters);
    setB.waves = setA.waves;
    
    double twl;        // thermal wave length
    twl = 2.0 * sqrt(2.0 * setA.parameters.lambda * setA.parameters.tau);
    int ncut;
    ncut = setA.parameters.nGrid - int(twl / setA.parameters.dr);
    double dm, it;
    double r1, r2;
    double pot1;

	for(int i = ncut; i < setA.parameters.nGrid; i++)
	{
        r1 = (i+1)*setB.parameters.dr;
        pot1 = setB.potential.value(r1);
		for(int j = ncut; j < setA.parameters.nGrid; j++)
		{
            // primitive approximation for interaction
            r2 = (j+1)*setB.parameters.dr;
            it = pot1 + setB.potential.value(r2);
            it = it * setB.parameters.tau;
            it = it * 0.5;
            it = exp(-it);
            for(int k = 0; k < int(setB.waves.size()); k++) {
			    dm = this->exprCalc.freePW(k,i,j,setB.parameters);
			    dm = dm * it;  
			    gsl_matrix_set(setB.waves[k], i, j, dm);
            }
		}
	}

    return setB;
}