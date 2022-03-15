# Project: FLATModel


#Variables for compiling

#Parallel
MPI_F77=mpif77
MPI_F90=mpif90
MPI_CC=mpicc

#serial
F77=xlf
F90=xlf90
CC=xlc

#OPENMP
F90_r=xlf90_r
export MP_FC=xlf90_r

#LIBRERIAS
LIBS= -lpthread      -lrt

#Linux
FC = gfortran

# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -g -fbounds-check
FCFLAGS = -O2

# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include

#Flags
#F90FLAGS= -qfree=f90 -qsuffix=f=f90
F90FLAGS= -ffree-line-length-none
# -ffree-form -fopenmp
#OPT_FLAGS= -O3 -qstrict -qtune=ppc970 -qcache=auto -qarch=ppc970 -qsmp=omp
OPT_FLAGS= -O3


#Compilation elements
OBJ  = GridOut_ASCII.o \
       DebrisDiceGlobals.mod \
       MaxPath.o \
       Stat.o \
       Normalize.o \
       Path.o \
       DebrisDice.o

LINKOBJ = GridOut_ASCII.o \
       DebrisDiceGlobals.mod \
       MaxPath.o \
       Stat.o \
       Normalize.o \
       Path.o \
       DebrisDice.o


#Default binary name
BIN  = DebrisDice
RM = /bin/rm -f


#Compilation alternatives
default: $(BIN)

all: $(BIN)

clean: clean-custom
	${RM} $(OBJ) $(BIN)



$(BIN): $(OBJ)
	$(FC) $(F90FLAGS) $(OPT_FLAGS) $(LINKOBJ) -o $(BIN) $(LIBS)

GridOut_ASCII.o: GridOut_ASCII.F90
	$(FC) $(F90FLAGS) $(OPT_FLAGS) -c GridOut_ASCII.F90 -o GridOut_ASCII.o

MontecarloGlobals.mod: MontecarloGlobals.F90
	$(FC) $(F90FLAGS) $(OPT_FLAGS) -c DebrisDiceGlobals.F90 -o DebrisDiceGlobals.mod

Stat.o: Stat.F90
	$(FC) $(F90FLAGS) $(OPT_FLAGS) -c Stat.F90 -o Stat.o

MaxPath.o: MaxPath.F90
	$(FC) $(F90FLAGS) $(OPT_FLAGS) -c MaxPath.F90 -o MaxPath.o

Montecarlo.o: Montecarlo.F90
	$(FC) $(F90FLAGS) $(OPT_FLAGS) -c DebrisDice.F90 -o DebrisDice.o

Normalize.o: Normalize.F90
	$(FC) $(F90FLAGS) $(OPT_FLAGS) -c Normalize.F90 -o Normalize.o

Path.o: Path.F90
	$(FC) $(F90FLAGS) $(OPT_FLAGS) -c Path.F90 -o Path.o





