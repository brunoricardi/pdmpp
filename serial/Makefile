# C++ compiler
CC=g++
COPTS=-Wall -std=gnu++2a

# GSL
GSLROOT=/home/babreu/NumericalLibraries/gsl-2.7
GSLINCLUDE=${GSLROOT}/include/gsl
GSLLIB=${GSLROOT}/lib

# LINKS
# gslcblas: GSL with BLAS
# gsl: GNU Scientific Library
LDFLAGS=-lblas -llapack -lgslcblas -lgsl 


objects: parameters.o potential.o

parameters.o: parameters.cpp
	$(CC) ${COPTS} -c $^

potential.o: potential.cpp
	$(CC) ${COPTS} -c $^


clean:
	rm ./*.o