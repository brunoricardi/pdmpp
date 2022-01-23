#ifndef DENSMAT_H
#define DENSMAT_H

#include "data.h"
#include "parameters.h"

using namespace std;

class DensMat {
    // constructor, destructor
    public:
        DensMat(Parameters _parameters);
        ~DensMat();

    // attributes
    public:
        Parameters parameters;
        vector< vector < vector <double> >> value;
        vector<double> rGrid;
        vector<double> aGrid;
};



#endif /* EXPRCALC_H */

