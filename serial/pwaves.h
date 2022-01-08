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
        std::vector<gsl_matrix*> waves;

};

#endif /* PWAVES_H */