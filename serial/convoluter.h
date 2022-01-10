#ifndef CONVOLUTER_H
#define CONVOLUTER_H

#include "parameters.h"
#include "pwaves.h"
#include "gsl.h"

class Convoluter {
    // constructor, destructor
    public:
        Convoluter();
        ~Convoluter();

    // attributes


    // methods
    public:
        gsl_matrix *run(gsl_matrix *A, gsl_matrix *B);
        PartWaves runOverWaves(PartWaves setA, PartWaves setB);
        PartWaves runSquarer(PartWaves setA);
};


#endif /* CONVOLUTER_H */