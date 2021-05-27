# C++ compiler
CPP=g++

# compilation flags
# blas: Basic Linear Algebra Subprograms
# lapack: Linear Algebra PACKage
# gslcblas: GSL with BLAS
# gsl: GNU Scientific Library
# std=gnu++14: standard GNU support for C++14
LDFLAGS=-lblas -llapack -lgslcblas -lgsl -std=gnu++14

# OpenMP flag
OMP=-fopenmp

# source files
SOURCE_SERIAL=./serial/pdm_serial.cpp
SOURCE_PARALLEL=./omp/pdm_omp.cpp

# path to where binary is going
PATH_BIN=./bin

# binary names
EXEC_SERIAL=pdm_serial
EXEC_PARALLEL=pdm_omp

serial:
	$(CPP) $(LDFLAGS) $(SOURCE_SERIAL) -o $(PATH_BIN)/$(EXEC_SERIAL)
	@echo -e "----- SERIAL COMPILATION DONE -----"


parallel:
	$(CPP) $(LDFLAGS) $(SOURCE_PARALLEL) -o $(PATH_BIN)/$(EXEC_PARALLEL)
	@echo -e "----- PARALLEL COMPILATION DONE -----"


clean:
	rm -r $(PATH_BIN)/pdm*


