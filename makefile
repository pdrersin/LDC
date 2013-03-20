FC=gfortran
FCFLAGS=-O3 -fopenmp -finline-functions #-traceback -openmp -parallel -fpp
LDC_INCLUDE=-I functions/


all: LDC

LDC: set_precision.o set_constants.o setup.o fileio.o matrix_manip.o solvers.o ldc.o
	$(FC) $(LDC_INCLUDE) $(FCFLAGS) set_precision.o set_constants.o setup.o fileio.o matrix_manip.o solvers.o ldc.o -o ldc

set_precision.o: set_precision.f90
	$(FC) $(LDC_INCLUDE) $(FCFLAGS) -c set_precision.f90

set_constants.o: set_precision.o set_constants.f90
	$(FC) $(LDC_INCLUDE) $(FCFLAGS) -c set_constants.f90

setup.o: set_precision.o set_constants.o setup.f90
	$(FC) $(LDC_INCLUDE) $(FCFLAGS) -c setup.f90

fileio.o: set_precision.o setup.o fileio.f90
	$(FC) $(LDC_INCLUDE) $(FCFLAGS) -c fileio.f90

matrix_manip.o: set_precision.o matrix_manip.f90
	$(FC) $(LDC_INCLUDE) $(FCFLAGS) -c matrix_manip.f90

solvers.o: set_precision.o set_constants.o setup.o matrix_manip.o solvers.f90
	$(FC) $(LDC_INCLUDE) $(FCFLAGS) -c solvers.f90

ldc.o: fileio.o setup.o solvers.o ldc.f90
	$(FC) $(LDC_INCLUDE) $(FCFLAGS) -c ldc.f90

clean:
	rm *.o *.mod