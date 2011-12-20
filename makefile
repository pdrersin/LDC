FC=gfortran
FCFLAGS=-O #-traceback -openmp -parallel -fpp

all: LDC

LDC: set_precision.o set_constants.o functions.o setup.o fileio.o ludcmp.o triblocksolve.o solvers.o ldc.o
	$(FC) $(FCFLAGS) set_precision.o set_constants.o functions.o setup.o fileio.o ludcmp.o triblocksolve.o solvers.o ldc.o -o ldc

set_precision.o: set_precision.f90
	$(FC) $(FCFLAGS) -c set_precision.f90

set_constants.o: set_constants.f90
	$(FC) $(FCFLAGS) -c set_constants.f90

functions.o: functions.f90
	$(FC) $(FCFLAGS) -c functions.f90

setup.o: setup.f90
	$(FC) $(FCFLAGS) -c setup.f90

fileio.o: fileio.f90
	$(FC) $(FCFLAGS) -c fileio.f90

ludcmp.o: ludcmp.f90
	$(FC) $(FCFLAGS) -c ludcmp.f90

triblocksolve.o: triblocksolve.f90
	$(FC) $(FCFLAGS) -c triblocksolve.f90

solvers.o: solvers.f90
	$(FC) $(FCFLAGS) -c solvers.f90

ldc.o: ldc.f90
	$(FC) $(FCFLAGS) -c ldc.f90

clean:
	rm *.o *.mod