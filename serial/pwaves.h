#ifndef PWAVES_H
#define PWAVES_H

#include "math.h"
#include "parameters.h"
#include "gsl.h"
#include "data.h"

class pwaves {
    public:
        pwaves(parameters _parameters);
        ~pwaves();
        parameters params;

    private:
        int allocate();

    public:
        std::vector<gsl_matrix*> waves;                 // the waves
        int initialize();                               // allocate and initialize with high temp expressions
        double freepw(int l, int i, int j, double tau); // calculate interactionless expressions

};

#endif /* PWAVES_H */