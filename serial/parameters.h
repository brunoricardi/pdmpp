#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "inout.h"
#include "data.h"
#include "constants.h"

class parameters {
    public:
        parameters();
        ~parameters();

    public:
        double lambda;	// hbar^2/2m
        int potId;		// interaction to be used
        int nGrid;		// number of grid points in r, r'
        double dr;		// r,r' grid spacing
        int na;		    // number of grid points in theta
        double da;		// theta grid spacing
        int nw;		    // number of partial waves
        int nSq;		// number of convolutions
        double tHigh;	// initial temperatute (final=tHigh/2^nsq)
        double tauMax;	// initial time step
        int nDim;		// number of spatial dimensions
        bool printPw;	// print partial waves?
        int printRn;	// renormalization constant so output files are not too large

    public:
        int readInput();            // read inputs from file pdm.inp
    
    private:
        void checkValues();    // set defaults for values not provided by client

};
#endif /* PARAMETERS_H */