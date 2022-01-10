#include "includes.h"

int main() {
    Parameters parameters;
    parameters.readInput();

    PartWaves waves(parameters);
    Convoluter convoluter;

    int check;
    check = waves.initialize();

    double mij;

    if(check == 0) {
        mij = gsl_matrix_get(waves.waves[0], 0, 0);
        std::cout << mij << std::endl;
        waves = convoluter.runSquarer(waves);
        mij = gsl_matrix_get(waves.waves[0], 0, 0);
        std::cout << mij << std::endl;
    }

    return 0;
}