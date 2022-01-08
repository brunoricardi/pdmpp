#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "parameters.h"
#include "data.h"
#include "math.h"

class potential {
    public:
        potential(parameters _parameters);
        ~potential();
        parameters params;

    public:
        double value(double r);                 // value of potential at a certain relative distance          
        double integrate(double r1, double r2); // integrate potential over a straight line between r1 and r2

};
#endif /* POTENTIAL_H */