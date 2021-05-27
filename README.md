# pdmpp

## Introduction
Welcome to PDM++. This is a scientific code written by Bruno Abreu (babreu.scientist@gmail.com).
PDM++ calculates pair density-matrices to be used in quantum many-body stochastic simulations, such as Path Integral Monte Carlo (PIMC) [cite].
The low temperature (short imaginary-time) expressions are obtained by numerical convolutions of their high-temperature expressions, which are calculated using a semi-classical approach [cite,cite]. 
Currently, a few interaction potentials are implemented:
- Free particle
- 6-12 Lennard-Jones
- Aziz1995 for He-3

Additional interactions can be easily implemented in one of the subroutines.

## How to cite this
As the *LICENSE* says, this is an open-source scientific code. You are free to modify, share, and use it however it suits your needs. However, if you do publish something using the code either on parts or as a whole, I largely appreciate if you can cite it (and I have the feeling you would understand why). Here's how you can do it:

Bruno R. de Abreu, *PDM++*, available at https://github.com/brunoricardi/pdmpp/ (2021)

I provided a *pdmpp.bib* file with this information.


## Overview of the code
Info on the functions

## Input/Output
pdm.inp, pdm.csv

## Details on the serial version
The file pdm.cpp contains the original code that has been used by me since 2012. Recently (2021), I took some time to work on optmizations. The source file /opt/pdm_opt.cpp contains a version that has loops and tensors adjusted for efficient use of the cache. This should be your choice for serial calculations.

## Details on the parallel version
I parallelized a large part of the code with OpenMP. If you can use all the cores in a node, this will greatly improve performance, up to about a 20x factor (see more details below). If you want to use the parallel version, your source file is /omp/pdm_omp.cpp. This has also been optimized for efficient cache usage.




## References
[1]
[2]
[3]
[4]
