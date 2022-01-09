#ifndef PWAVES_H
#define PWAVES_H

#include "math.h"
#include "parameters.h"
#include "gsl.h"
#include "data.h"
#include "exprcalc.h"

class pwaves {
    // constructor, destructor
    public:
        pwaves(parameters _parameters);
        ~pwaves();

    // ATTRIBUTES
    public:
        parameters params;
        std::vector<gsl_matrix*> waves;  
    private:
        exprcalc exprEngine;
    
    // METHODS
    public:
        int initialize();       // allocate and initialize with high temp expressions
    private:
        int allocate();         // allocates space for GSL matrices

};

#endif /* PWAVES_H */