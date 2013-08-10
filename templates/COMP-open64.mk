# -*- mode: makefile -*-

##############################################################################
# Open64 Compiler suite
#
##############################################################################
COMPILERNAME = OPEN64

# Default: No compiler wrapper commands
# This works both for builds of serial and parallel applications
F77       = openf90
F90       = openf90
CC        = opencc
CXX       = openCC
LD        = openf90

# Compiler flag to specify the directory where module files should be
# placed when created and where they should be searched for.
# Note: Do not remove ticks and whitespace!
MODOPTION = '-module '

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
F77VERSION = $(F77) -v 2>&1 | cat
F90VERSION = $(F90) -v 2>&1 | cat
CCVERSION  = $(CC)  -v 2>&1 | cat
CXXVERSION = $(CXX) -v 2>&1 | cat

# Detect compiler version
OPEN64VERSION  := $(shell eval $(F90VERSION) | \
		    sed -n -e 's/^.*Version \([0-9]*\.[0-9]*\.*[0-9]*\).*$$/\1/p; q;')
ifneq ($(OPEN64VERSION),)
OPEN64VERSION_MAJOR := $(shell echo $(OPEN64VERSION) | cut -d. -f1)
OPEN64VERSION_MINOR := $(shell echo $(OPEN64VERSION) | cut -d. -f2)
else
OPEN64VERSION_MAJOR := 0
OPEN64VERSION_MINOR := 0
endif

# Functions to detect minimal compiler version
open64minversion = $(shell if [ $(OPEN64VERSION_MAJOR) -gt $(1) ] || \
	                     ([ $(OPEN64VERSION_MAJOR) -ge $(1) ] && \
			      [ $(OPEN64VERSION_MINOR) -ge $(2) ]) ; then echo yes ; else echo no ; fi)

# Functions to detect maximal compiler version
open64maxversion = $(shell if [ $(OPEN64VERSION_MAJOR) -lt $(1) ] || \
	                     ([ $(OPEN64VERSION_MAJOR) -le $(1) ] && \
			      [ $(OPEN64VERSION_MINOR) -le $(2) ]) ; then echo yes ; else echo no ; fi)


##############################################################################
# Compiler flags
# (including non-architecture specific optimisation flags)
##############################################################################

# Set default type of integer variables explicitly
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -DUSE_LARGEINT -i8
endif
# $(CC) and $(CXX) do not have such a corresponding option, so we have to
# pray that they default the 'int' type properly.



# Specify -openmp for all Open64 compilers
ifeq ($(strip $(OPENMP)), YES)
CFLAGSF77     := -DUSE_OPENMP -openmp $(CFLAGSF77)
CFLAGSC       := -DUSE_OPENMP -openmp $(CFLAGSC)
LDFLAGS       := -DUSE_OPENMP -openmp $(LDFLAGS)
endif



ifneq (,$(findstring EXPENSIVE ,$(OPT)))
# Specify -Ofast -ipa for all Open64 compilers
CFLAGSF77     := -Ofast -ipa $(CFLAGSF77)
CFLAGSC       := -Ofast -ipa $(CFLAGSC)
LDFLAGS       := -Ofast -ipa $(LDFLAGS)
endif



# Check if Fortran source code files are preprocessed by F90CPP script.
# Otherwise, set additional flags to enforce use of CPP preprocessor.
ifeq ($(strip $(F90CPP)),)
CFLAGSF77 := $(CFLAGSF77) -cpp -I$(FEAT2BASEDIR)
# Note: Do not remove trailing whitespace!
MODOPTION = -module 
endif



# Set default compile flags
ifeq ($(call optimise), YES)
CFLAGSF77     := -DUSE_COMPILER_OPEN64 $(CFLAGSF77) -O2 -LNO -mso \
	          -ffast-math -ffast-stdlib -CG:compute_to=on \
		  -finline-functions -inline -fno-second-underscore
CFLAGSF90     := $(CFLAGSF90) $(CFLAGSF77)
CFLAGSC       := -DUSE_COMPILER_OPEN64 $(CFLAGSC) -O2 -LNO -mso \
		 -ffast-math -ffast-stdlib -CG:compute_to=on \
		 -finline-functions -inline -static
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
else
CFLAGSF77     := -DUSE_COMPILER_OPEN64 $(CFLAGSF77) -O0 -g3 \
	         -fno-second-underscore
CFLAGSF90     := $(CFLAGSF90) $(CFLAGSF77) -ffortran-bounds-check -fullwarn
CFLAGSC       := -DUSE_COMPILER_OPEN64 $(CFLAGSC) -O0 -g3 -trapuv
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
endif



##############################################################################
# Non-standard features supported by compiler
##############################################################################
CFLAGSF90 := -DHAS_INTRINSIC_FLUSH \
	     -DHAS_INTRINSIC_ISNAN \
	     -DHAS_ISO_C_BINDING $(CFLAGSF90)



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
