# -*- mode: makefile -*-

##############################################################################
# GNU Compiler suite 4.x
#
##############################################################################
COMPILERNAME = GFORTRAN

# Default: No compiler wrapper commands
# This works both for builds of serial and parallel applications
F77       = openf90
F90       = openf90
CC        = opencc
CXX       = openCC
LD        = openf90

# If preprocessor switch -DENABLE_SERIAL_BUILD does not occur in compiler flags,
# a build for parallel execution is requested.
ifeq (,$(findstring -DENABLE_SERIAL_BUILD ,$(APPONLYFLAGS) $(CFLAGSF90) ))
ifeq ($(strip $(MPIWRAPPERS)), YES)
F77       = mpif77
F90       = mpif90
CC        = mpicc
CXX	  = mpic++
LD        = mpif90
endif
endif


##############################################################################
# Commands to get version information from compiler
##############################################################################
F77VERSION = $(F77) -version 2>&1 | head -n 1
F90VERSION = $(F90) -version 2>&1 | head -n 1
CCVERSION  = $(CC)  -version 2>&1 | head -n 1
CXXVERSION = $(CXX) -version 2>&1 | head -n 1


##############################################################################
# compiler flags 
# (including non-architecture specific optimisation flags)
##############################################################################

# Set default type of integer variables explicitly
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -i8
CFLAGSF90     := $(CFLAGSF90) -i8
endif
# $(CC) and $(CXX) do not have such a corresponding option, so we have to 
# pray that they default the 'int' type properly.



# Specify -fopenmp for all gcc compilers
ifeq ($(strip $(OPENMP)), YES)
CFLAGSF77LIBS := -DUSE_OPENMP -openmp $(CFLAGSF77LIBS)
CFLAGSF77     := -DUSE_OPENMP -openmp $(CFLAGSF77)
CFLAGSF90     := -DUSE_OPENMP -openmp $(CFLAGSF90)
CFLAGSC       := -DUSE_OPENMP -openmp $(CFLAGSC)
CFLAGSCXX     := -DUSE_OPENMP -openmp $(CFLAGSCXX)
LDFLAGS       := -DUSE_OPENMP -openmp $(LDFLAGS)
endif



# Set default compile flags
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -DUSE_COMPILER_OPEN64 $(CFLAGSF77LIBS) -O2 -LNO -mso \
	          -ffast-math -ffast-stdlib -CG:compute_to=on \
		  -finline-functions -inline -fno-second-underscore
CFLAGSF77     := $(CFLAGSF77LIBS) $(CFLAGSF77)
CFLAGSF90     := -DENABLE_USE_ONLY -DHAS_INTRINSIC_FLUSH \
	         $(CFLAGSF90) $(CFLAGSF77LIBS) \
		 -module $(OBJDIR)
CFLAGSC       := -DUSE_COMPILER_OPEN64 $(CFLAGSC) -O2 -LNO -mso \
		 -ffast-math -ffast-stdlib -CG:compute_to=on \
		 -finline-functions -inline -static
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
else
CFLAGSF77LIBS := -DUSE_COMPILER_OPEN64 $(CFLAGSF77LIBS) -O0 -g3 \
	         -fno-second-underscore
CFLAGSF77     := $(CFLAGSF77LIBS) $(CFLAGSF77)
CFLAGSF90     := -DENABLE_USE_ONLY -DHAS_INTRINSIC_FLUSH \
		 $(CFLAGSF90) $(CFLAGSF77LIBS) \
		 -module $(OBJDIR) -ffortran-bounds-check -fullwarn
CFLAGSC       := -DUSE_COMPILER_OPEN64 $(CFLAGSC) -O0 -g3 -trapuv
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
endif



ifeq ($(strip $(OPT)), EXPENSIVE)
CFLAGSF77LIBS := -Ofast -ipa $(CFLAGSF77LIBS)
CFLAGSF77     := -Ofast -ipa $(CFLAGSF77)
CFLAGSF90     := -Ofast -ipa $(CFLAGSF90)
CFLAGSC       := -O2    -ipa $(CFLAGSC)
CFLAGSCXX     := -Ofast -ipa $(CFLAGSCXX)
LDFLAGS       := -Ofast -ipa $(LDFLAGS)
endif



# Detect compiler version
OPEN64VERSION  := $(shell eval $(F90VERSION) )



##############################################################################
# Extension compiler uses for name of generated module information files
# (some compilers generate .mod files, others .MOD files)
##############################################################################
MODEXTENSION = mod


##############################################################################
# Manual moving of generated module information files to 
# object directory needed?
##############################################################################
MOVEMOD   = NO


##############################################################################
# Commands needed by the Sparse Banded Blas benchmark
##############################################################################
SBB_CVERSIONCMD  = $(F77) -V  2>&1 | sed 's|(R)||g; 1!d;'


# The settings needed to compile a FEAT2 application are "wildly" distributed
# over several files ((Makefile.inc and templates/*.mk) and if-branches 
# (in an attempt to reduce duplicate code and inconsistencies among all build 
# IDs that e.g. use the same MPI environment). Not having all settings at 
# *one* place entails the risk (especially in the event of setting up a new 
# build ID) that settings are incompletely defined. A simple typo in a matching 
# rule in Makefile.inc may prevent that the compiler and compiler command line
# flags are set. Compilation would fail with the most peculiar errors - if not
# the Makefile had been set up to catch such a case.
# Each build ID in FEAT2 has 6 tokens: architecture, cpu, operating system,
# compiler family, BLAS implementation, MPI environment. Whenever setting
# one of these, an according flag is set. They are named TOKEN1 up to TOKEN6.
# Before starting to actually compile a FEAT2 application, every Makefile
# generated by bin/configure checks whether all six tokens are set *for the 
# choosen build ID*. If not, the user gets an error message describing exactly 
# what information is missing, e.g. token 5 not set which means there is no
# information available which BLAS implementation to use and where to find the
# library.
#
# In this file, the fourth token has been set: compiler and compiler command
# line options. Set the flag accordingly.
TOKEN4 := 1
