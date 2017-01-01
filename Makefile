
# choose fortran compiler
#FC = gfortran
#LD = gfortran

FC = ifort
LD = ifort


# Fortran compiler flags
# Fortran compiler flags
#intel fortran compiler debug mode: 
# FFLAGS  = -g -O0 -fpe:0 -warn declarations -warn unused -warn ignore_loc -warn truncated_source -traceback -check all -implicitnone -openmp
# LDFLAGS = -mkl
#intel fortran compiler optimized mode: 
FFLAGS  = -O3 -fpe:0 -implicitnone -fast -ipo -xHost -parallel -openmp
LDFLAGS = -mkl

#gfortran compiler debug mode: 
#FFLAGS  = -Wall -g -fbounds-check
#LDFLAGS = 
#gfortran compiler optimized mode: 
#FFLAGS  = -O2 -s -fomit-frame-pointer -fexpensive-optimizations -ffast-math
#LDFLAGS = 


.DEFAULT:
	-touch $@
all: a.out
basic_io.o: ./basic_io.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./basic_io.f90
boundary.o: ./boundary.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./boundary.f90
fluxes.o: ./fluxes.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluxes.f90
harten.o: ./harten.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./harten.f90
main.o: ./main.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./main.f90
preproc.o: ./preproc.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./preproc.f90
shared_vars.o: ./shared_vars.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./shared_vars.f90
time_advance.o: ./time_advance.f90 shared_vars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./time_advance.f90
SRC = ./shared_vars.f90 ./fluxes.f90 ./preproc.f90 ./time_advance.f90 ./boundary.f90 ./main.f90 ./harten.f90 ./basic_io.f90
OBJ = shared_vars.o fluxes.o preproc.o time_advance.o boundary.o main.o harten.o basic_io.o
clean: neat
	-rm -f .cppdefs $(OBJ) *.mod a.out *.pdf *.dat
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
a.out: $(OBJ) 
	$(LD) $(OBJ) -o a.out  $(LDFLAGS)
