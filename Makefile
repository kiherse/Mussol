################################################################################
# RATPENAT Makefile
################################################################################

#GCC
H5FC = h5pfc
FC = gfortran
H5PFC = h5pfc
MPIFC = mpif90

#inputs para hacer pruebas
HDF=1
OPENMP=1
MPI=1
DBG=0
MN=0

# Options

ifeq ($(HDF),1)
	UHDF = HDF
        ifeq ($(MPI),1)
                FC = $(H5PFC)
		UPARALELO=PARALELO
	else
                FC = $(H5FC)
		UPARALELO=NOPARALELO	
        endif
else
	UHDF = NOHDF
        ifeq ($(MPI),1)
                FC = $(MPIFC)
		UPARALELO=PARALELO
        else
		ifeq ($(IBM),1)
			FC = $(IFC)
		endif
                UPARALELO=NOPARALELO
        endif
endif


ifeq ($(OPENMP),1)
        UOPENMP = OPENMP
else
        UOPENMP = NOOPENMP
endif

	UCASERO=CASERO
	ifeq ($(OPENMP),1)
           OMPF = -fopenmp
	endif

        ifeq ($(DBG),1)
            FCOPT = -ffree-form -ffpe-trap=invalid,zero,overflow -ffree-line-length-170 -m64 -march=native -fimplicit-none -O3 -fcheck=all -g -fimplicit-none  
            #-Wall #-Wextra -Wuninitialized
        else
            FCOPT = -ffree-form -ffree-line-length-170 -m64 -march=native -fimplicit-none -O3 #-xHost #-r8
        endif 

        DEFINES =  -D$(UOPENMP) -D$(UHDF) -D$(UPARALELO) -D$(UCASERO) # -DNONSTCPP
        FLAGS = -free #-x f95-cpp-input
        FFLAGS = $(FCOPT) $(OMPF) $(FLAGS) $(DEFINES)


################################################################################

# DIRECTORIES
SRCDIR     := src
BINDIR     := bin
VPATH := $(SRCDIR) $(dir $(wildcard $(SRCDIR)/*/.))

MODULE_SRC := modulos.F mussol.F integrar.F paralelo.F rst_hdf5.F input.F
HEADERS = types.h

# EXECUTABLES
MAIN_EXE := RATPENAT
FFLAGS+= -J $(BINDIR)

# CREATE OBJECTS
MODULE_SRC := $(foreach file,$(MODULE_SRC),$(SRCDIR)/$(file))
MODULE_OBJ := $(patsubst %,$(BINDIR)/%,$(notdir $(MODULE_SRC:.F=.o)))

ifneq ($(BINDIR),)
$(shell test -d $(BINDIR) || mkdir -p $(BINDIR))
endif

################################################################################

RATPENAT: $(MODULE_OBJ)
	$(FC) $(FFLAGS) $^ -o $@

$(BINDIR)/%.o $(BINDIR)/%.mod : %.F $(HEADERS)
	$(FC) $(FFLAGS) -c -o $@ $<

################################################################################
clean:
	rm -rf *.bck
	rm -rf *.lst
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.f
	rm -rf *.err *.out
	rm -rf RATPENAT