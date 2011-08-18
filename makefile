FC=ifort
FCFLAGS=-O -g -pg #-traceback -openmp -parallel -fpp

all: LDC

LDC: set_precision.o ludcmp.o triblocksolve.o LidDrivenCavityImplicit.o
	$(FC) $(FCFLAGS) set_precision.o ludcmp.o triblocksolve.o LidDrivenCavityImplicit.o -o LDC

LidDrivenCavityImplicit.o: LidDrivenCavityImplicit.f90
	$(FC) $(FCFLAGS) -c LidDrivenCavityImplicit.f90

set_precision.o: set_precision.f90
	$(FC) $(FCFLAGS) -c set_precision.f90

ludcmp.o: ludcmp.f90
	$(FC) $(FCFLAGS) -c ludcmp.f90

triblocksolve.o: triblocksolve.f90
	$(FC) $(FCFLAGS) -c triblocksolve.f90

clean:
	rm *.o *.mod