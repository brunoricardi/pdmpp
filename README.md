# pdmpp

## Introduction
Welcome to PDM++. This is a scientific code written by Bruno Abreu (babreu.scientist@gmail.com).
PDM++ calculates pair density-matrices to be used in quantum many-body stochastic simulations, such as Path Integral Monte Carlo (PIMC) [1].
The low temperature (short imaginary-time) expressions are obtained by numerical convolutions of their high-temperature expressions, which are calculated using a semi-classical approach [2,3]. 
Currently, a few interaction potentials are implemented:
- Free particle
- 6-12 Lennard-Jones
- Aziz1995 for He-3

Additional interactions can be easily implemented in one of the subroutines.

## How to cite this
As the *LICENSE* says, this is an open-source scientific code. You are free to modify, share, and use it however it suits your needs. However, if you do publish something using the code either on parts or as a whole, I largely appreciate if you can cite it (and I have the feeling you would understand why). Here's how you can do it:

Bruno R. de Abreu, *PDM++*, available at https://github.com/brunoricardi/pdmpp (2021)

I provided a *pdmpp.bib* file with this information.


## Overview of the code
The code is composed of 7 functions and a main part. What each of them do is commented in the source files as well:
- read_parameters(): this looks for a *pmd.inp* file that will direct how to calcualtion is performed (see details below)
- check_global_variables(): this is just a quick check to make sure there's nothing absurd (e.g 0 dimension)
- asympbessel_pf(double, double): calculates the pre-factor of asymptotic expressions for Bessel functions
- pot(double): this is where the physics come in, with the specification of the interaction potential
- interaction(double, double): this integrates the potential between to spatial points, which originates from the WKB approach
- free_pw(double, double, int, double, int): calculates free-particle expressions for partial waves
- correction_integral(double, double, double): caclulates a correction for the truncation of the grid

In general terms, the code starts by reading the input file and checking to see if the variables are okay to go. 
It then creates all of the tensors that are needed, and populates them with the high-temprature expressions. After that, the numerical convolutions are performed. Finally, the partial wave sum is done, and the result is printed to the output file.

## Input/Output
Input must be in the file *pdm.inp* (sample provided here), which will be searched at execution time in the directory where you are running it. This file consists of a series of commands (five capital letters) followed by their values, separated by blank spaces. Each line corresponds to one command. They are:
- **PHYSI** double *lambda=hbar/2m* int *spatial_dimensions* int *potential_id*
  
  The potential_id can be found in the function *pot* of the source code.
- **RGRID** int *number_of_steps* double *step_size*
  
  This is the linear grid in the spatial variables. 
- **AGRID** int *number_of_steps*
  
  This is the linear grid in the angular variable. Note that the range is fixed: [0,pi).
- **NWAVE** int *number_of_waves*
  
  Total number of partial waves to be considered.
- **SQUAR** int *nsquares* double *initial_temp*
  
  Number of squarings to be performed. Each convolution reduces the initial temperature by a factor of two, so your final temperature is going to be *initial_temp* / 2^*nsquares*
- **PRWAV**
  
  This is the only optional command (no arguments to it). If presented, each partial wave will be printed in a different file.

The output is a *pdm.csv* that has four columns: *r1* (initial relative distance), *r2* (final relative distance), *theta* (angle between initial and final relative distances) and *pdm* (value of the pair density matrix). It can be a pretty large file (few GBs) depending on the size of your grids.



## Details on the serial version
The source file */serial/pdm_serial.cpp* contains a version that has loops and tensors adjusted for efficient use of the cache. This should be your choice for serial calculations. 


## Details on the parallel version
I parallelized a large part of the code with OpenMP. If you can use all the cores in a node, this will greatly improve performance, up to about a 20x factor. Some loops are parallelized in the number of partial waves, some others are parallelized in the spartial grid. A general rule for good performance is that you choose RGRID, AGRID and NWAVE values that are multiples of the number of threads you are goilg to call. If you want to use the parallel version, your source file is */omp/pdm_omp.cpp*. This has also been optimized for efficient cache usage.

## Additional documentation
I'm working on a PDF with the theoretical grounds for this code. It will be listed here at some point (I hope).

## Prospects
I have many ideas for improving this code, but very little time to implement them. If you would like to collaborate, do not hesitate and contact me!


## References
[1] Herman, M. F., E. J. Bruskin, and B. J. Berne. "On path integral Monte Carlo simulations." The Journal of Chemical Physics 76, no. 10 (1982): 5150-5155.

[2] Ceperley, David M. "Path integrals in the theory of condensed helium." Reviews of Modern Physics 67, no. 2 (1995): 279.

[3] de Abreu, Bruno R., Fabio Cinti, and Tommaso Macrì. "Superstripes and quasicrystals in bosonic systems with hard-soft corona interactions." arXiv preprint arXiv:2009.10203 (2020).
