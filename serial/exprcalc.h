#ifndef EXPRCALC_H
#define EXPRCALC_H

#include "gsl.h"
#include "math.h"
#include "parameters.h"

class exprcalc {
    public:
        exprcalc(parameters _parameters);
        ~exprcalc();
        parameters params;

    public:
        double freePW(int l, int i, int j, double tau);     // free partial wave expression
        double asympBesselPF(double n, double x);           // pre-factor of the asymptotic expression for Bessel functions
};

#endif /* EXPRCALC_H */