# C++ compiler
CC=g++
COPTS=-Wall -std=gnu++2a

# GSL
GSLROOT=/home/babreu/NumericalLibraries/gsl-2.7
GSLINCLUDE=${GSLROOT}/include
GSLLIB=${GSLROOT}/lib

# LINKS
# gslcblas: GSL with BLAS
# gsl: GNU Scientific Library
LDFLAGS=-lgsl -lgslcblas -lm


all: objects
	$(CC) -L${GSLLIB} -o pdm.exe *.o ${LDFLAGS}

objects: convoluter.o exprcalc.o parameters.o pdm.o potential.o pwaves.o

convoluter.o: convoluter.cpp
	$(CC) ${COPTS} -I${GSLINCLUDE} -c $^

exprcalc.o: exprcalc.cpp
	$(CC) ${COPTS} -I${GSLINCLUDE} -c $^

parameters.o: parameters.cpp
	$(CC) ${COPTS} -c $^

pdm.o: pdm.cpp
	$(CC) ${COPTS} -I${GSLINCLUDE} -c $^

potential.o: potential.cpp
	$(CC) ${COPTS} -c $^

pwaves.o: pwaves.cpp
	$(CC) ${COPTS} -I${GSLINCLUDE} -c $^


clean:
	rm ./*.o