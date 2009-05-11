# -*- mode: makefile -*-

##############################################################################
# GNU Compiler suite 4.x
#
##############################################################################
COMPILERNAME = GFORTRAN

# Default: No compiler wrapper commands
# This works both for builds of serial and parallel applications
F77       = gfortran
F90       = gfortran
CC        = gcc
CXX       = g++
LD        = gfortran

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
F77VERSION = $(F77) --version | head -n 1
F90VERSION = $(F90) --version | head -n 1
CCVERSION  = $(CC)  --version | head -n 1
CXXVERSION = $(CXX) --version | head -n 1


##############################################################################
# compiler flags 
# (including non-architecture specific optimisation flags)
##############################################################################

# WARNING WARNING WARNING
# All integer variables in FEAT2 are explicitly typed to either 32 or 64 bits.
# The only native integers are literals and code written in F77 (blas, lapack, 
# sbblas) and C (metis, coproc backend). FEAT2 assumes that all these (native) 
# integers are 32-bit!!!
# So, to get things running with compilers that do not default native integers
# to 32 bits, we need to add an appropriate compiler flag to
# CFLAGSF77LIBS:. 
# This also applies when changing the kind-values in kernel/fsystem.f90.

# $(CC) and $(CXX) do not have such a corresponding option, so we have to 
# pray that they default the 'int' type properly.

# detect compiler version
GFORTRANVERSION  := $(shell eval $(F90VERSION) )

# command line options for gfortran 4.3.x
ifneq (,$(findstring 4.3.,$(GFORTRANVERSION)))
CFLAGSF77LIBS := $(CFLAGSF77LIBS)
endif

# command line options for gfortran 4.4.x
ifneq (,$(findstring 4.4.,$(GFORTRANVERSION)))
endif


ifeq ($(call optimise), YES)
# ****************************************** #
# First version of optimised compiler flags! #
# ****************************************** #
CFLAGSF77LIBS := -DUSE_COMPILER_GCC $(CFLAGSF77LIBS) -O3  \
		 -ffast-math -foptimize-register-move -fprefetch-loop-arrays \
		 -funroll-loops -static -fno-second-underscore
CFLAGSF77     := $(CFLAGSF77LIBS) $(CFLAGSF77)
CFLAGSF90     := -DENABLE_USE_ONLY -DHAS_INTRINSIC_FLUSH -DHAS_INTRINSIC_IARGC \
	         -DHAS_INTRINSIC_ISATTY $(CFLAGSF90) $(CFLAGSF77LIBS) \
		 -J$(OBJDIR) -Wuninitialized
CFLAGSC       := -DUSE_COMPILER_GCC $(CFLAGSC) -O3 \
		 -ffast-math -foptimize-register-move -fprefetch-loop-arrays \
		 -funroll-loops -static
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
else
CFLAGSF77LIBS := -DUSE_COMPILER_GCC $(CFLAGSF77LIBS) -O0 -g -fno-second-underscore \
		 -fbacktrace -fexternal-blas #-pg
CFLAGSF77     := $(CFLAGSF77LIBS) $(CFLAGSF77)
CFLAGSF90     := -DENABLE_USE_ONLY -DHAS_INTRINSIC_FLUSH -DHAS_INTRINSIC_IARGC \
		 -DHAS_INTRINSIC_ISATTY $(CFLAGSF90) $(CFLAGSF77LIBS) \
		 -J$(OBJDIR) -fbounds-check \
		 -Wcharacter-truncation -Winline \
		 -Wline-truncation -Wsurprising  \
		 -Wunreachable-code -Wunused-label -Wunused-variable
		 # -Wnonstd-intrinsics: not available in 4.4.0
		 # -Wuninitialized -Wimplicit-interface -Wunused-variable
# Do not specify:
# * -std=f95 
#   as it gives conflicts with LAM/MPI's mpif.h
# * -Wimplicit-interface
#   as it gives warnings for every MPI call, BLAS call etc.

CFLAGSC       := -DUSE_COMPILER_GCC $(CFLAGSC) -g -fbounds-check #-pg
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS) #-pg
endif



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
SBB_CVERSIONCMD  = $(F77) --version  | sed 's|[(|)]||g; 1!d;'


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
