#ifndef PWAVES_H
#define PWAVES_H

#include "math.h"
#include "parameters.h"
#include "gsl.h"
#include "data.h"
#include "exprcalc.h"
#include "potential.h"

class PartWaves {
    // constructor, destructor
    public:
        PartWaves(Parameters _parameters);
        ~PartWaves();

    // ATTRIBUTES
    public:
        Parameters parameters;
        std::vector<gsl_matrix*> waves;  
    private:
        ExprCalc exprCalc;
        Potential potential;
    
    // METHODS
    public:
        int initialize();       // allocate and initialize with high temp expressions
    private:
        int allocate();         // allocates space for GSL matrices

};

#endif /* PWAVES_H */