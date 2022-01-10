#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "parameters.h"
#include "data.h"
#include "math.h"

class Potential {
    // constructor, destructor
    public:
        Potential(Parameters _parameters);
        ~Potential();

    // attributes
    public:
        Parameters parameters;

    // methods
    public:
        double value(double r);                 // value of potential at a certain relative distance          
        double integrate(double r1, double r2); // integrate potential over a straight line between r1 and r2

};
#endif /* POTENTIAL_H */