# Note: we can compare the times with option -O2 and -O3 and if they are huge
# different that means the code can be optimized very much
#
#FCC=ifort
FCC=gfortran
#FCC=mpifort

ID=Ising  # the dynamics of QSI with pairing

THERMALIZATION=
RESTART=
ACF=
ANNEALING=
HISTOGRAM=

ifeq ($(FCC),ifort)
	EXEC=$(ID)
	# FFTW
	MY_INCLUDE_PATH=-I$(HOME)/opt/fftw-icc/include
	MY_LIBRARY_PATH=-L$(HOME)/opt/fftw-icc/lib
	LIBS=-lfftw3

	FC_debug=-g -O0 -CU -CB
	FC=-O2 -cpp -DIFORT -DACF
	#FC=-O0 -g -CU -CB -cpp -DIFORT -DACF
else ifeq ($(FCC),gfortran)
	EXEC=$(ID)
	# FFTW
	MY_INCLUDE_PATH=-I$(HOME)/opt/fftw-gcc/include
	MY_LIBRARY_PATH=-L$(HOME)/opt/fftw-gcc/lib
	LIBS=-lfftw3

	FC=-O1 -cpp -DGFORT $(THERMALIZATION) $(RESTART) $(ACF) $(ANNEALING) $(HISTOGRAM)
	#FC=-O0 -cpp -DGFORT -fbounds-check -Wuninitialized $(THERMALIZATION) $(RESTART) $(ACF)
	#FC=-O0 -g -cpp -DGFORT -DACF -fbounds-check -Wuninitialized
else ifeq ($(FCC),mpifort)
	#EXEC=$(ID)-mpi
	EXEC=$(ID)
	# FFTW
	MY_INCLUDE_PATH=-I$(HOME)/opt/fftw-gcc/include
	MY_LIBRARY_PATH=-L$(HOME)/opt/fftw-gcc/lib
	#MY_INCLUDE_PATH=-I$(HOME)/opt/fftw-icc/include
	#MY_LIBRARY_PATH=-L$(HOME)/opt/fftw-icc/lib
	LIBS=-lfftw3

	#FC=-O2 -fpp -DMPI -DIFORT
	#FC=-O0 -g -CU -CB -fpp -DMPI -DIFORT
	FC=-O1 -cpp -DMPI -DGFORT
	FC=-O0 -g -cpp -DMPI -DGFORT -fbounds-check -Wuninitialized
endif


F90=mpi_run.f90


#PROF=prof
ifdef PROF
ifeq ($(FCC),ifort)
	FC=-p -O0
	#FC=-p -O2
else
	FC=-pg -O0
	#FC=-pg -O2
endif
EXEC=$(PROF)
endif

all: $(F90)
	$(FCC) $(FC) $(F90) -o $(EXEC) $(MY_INCLUDE_PATH) $(MY_LIBRARY_PATH) $(LIBS)
	rm -f *.mod

run:
	./$(EXEC)

debug: $(F90)
	$(FCC) $(FC_debug) $(F90) $(Oper) $(FFTW) $(LIB) -o debug
	idbc debug

clean:
	rm -f *.mod $(EXEC)
