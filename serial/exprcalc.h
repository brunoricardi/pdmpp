#ifndef EXPRCALC_H
#define EXPRCALC_H

#include "gsl.h"
#include "math.h"
#include "parameters.h"

class ExprCalc {
    // constructor, destructor
    public:
        ExprCalc();
        ~ExprCalc();

    // attributes

    // methods
    public:
        double freePW(int l, int i, int j, Parameters parameters);  // free partial wave expression
        double asympBesselPF(double n, double x);                   // pre-factor of the asymptotic expression for Bessel functions
};

#endif /* EXPRCALC_H */