FC=gfortran
FCFLAGS=-march=native -funroll-loops -O2

all: LDC

LDC: set_precision.o ludcmp.o triblocksolve.o LidDrivenCavityImplicit.o
	$(FC) $(FCFLAGS) set_precision.o ludcmp.o triblocksolve.o LidDrivenCavityImplicit.o -o LDC

LidDrivenCavityImplicit.o: LidDrivenCavityImplicit.f95
	$(FC) $(FCFLAGS) -c LidDrivenCavityImplicit.f95

set_precision.o: set_precision.f90
	$(FC) $(FCFLAGS) -c set_precision.f90

ludcmp.o: ludcmp.f90
	$(FC) $(FCFLAGS) -c ludcmp.f90

triblocksolve.o: triblocksolve.f90
	$(FC) $(FCFLAGS) -c triblocksolve.f90