#ifndef CONVOLUTER_H
#define CONVOLUTER_H

#include "parameters.h"
#include "pwaves.h"
#include "gsl.h"
#include "exprcalc.h"

class Convoluter {
    // constructor, destructor
    public:
        Convoluter();
        ~Convoluter();

    // attributes
    private:
        ExprCalc exprCalc;

    // methods
    public:
        PartWaves runSquarer(PartWaves setA);                   // run squarer over partial waves
    private:
        gsl_matrix *run(gsl_matrix *A, gsl_matrix *B, double dr);          // multiply two matrices
        PartWaves runOverWaves(PartWaves setA, PartWaves setB); // multiply entire set wave by wave
        PartWaves truncCorrector(PartWaves setA);               // truncation corrector after multiplication
};


#endif /* CONVOLUTER_H */