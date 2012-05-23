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

# Set default type of integer variables explicitly
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -DUSE_LARGEINT -fdefault-integer-8
CFLAGSF90     := $(CFLAGSF90) -DUSE_LARGEINT -fdefault-integer-8
endif
# $(CC) and $(CXX) do not have such a corresponding option, so we have to 
# pray that they default the 'int' type properly.



# Specify -fopenmp for all gcc compilers
ifeq ($(strip $(OPENMP)), YES)
CFLAGSF77LIBS := -DUSE_OPENMP -fopenmp $(CFLAGSF77LIBS)
CFLAGSC       := -DUSE_OPENMP -fopenmp $(CFLAGSC)
LDFLAGS       := -DUSE_OPENMP -fopenmp $(LDFLAGS)
endif



# Set default compile flags
ifeq ($(call optimise), YES)
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



# Detect compiler version
GFORTRANVERSION  := $(shell eval $(F90VERSION) )

# Functions to detect minimal compiler version
gfortranminversion_4_7=\
	$(if $(findstring 4.7.,$(GFORTRANVERSION)),yes,no)
gfortranminversion_4_6=\
	$(if $(findstring yes,\
	$(call gfortranminversion_4_7) \
	$(if $(findstring 4.6.,$(GFORTRANVERSION)),yes,no)),yes,no)
gfortranminversion_4_5=\
	$(if $(findstring yes,\
	$(call gfortranminversion_4_6) \
	$(if $(findstring 4.5.,$(GFORTRANVERSION)),yes,no)),yes,no)
gfortranminversion_4_4=\
	$(if $(findstring yes,\
	$(call gfortranminversion_4_5) \
	$(if $(findstring 4.4.,$(GFORTRANVERSION)),yes,no)),yes,no)
gfortranminversion_4_3=\
	$(if $(findstring yes,\
	$(call gfortranminversion_4_4) \
	$(if $(findstring 4.3.,$(GFORTRANVERSION)),yes,no)),yes,no)
gfortranminversion_4_2=\
	$(if $(findstring yes,\
	$(call gfortranminversion_4_3) \
	$(if $(findstring 4.2.,$(GFORTRANVERSION)),yes,no)),yes,no)
gfortranminversion_4_1=\
	$(if $(findstring yes,\
	$(call gfortranminversion_4_2) \
	$(if $(findstring 4.1.,$(GFORTRANVERSION)),yes,no)),yes,no)

# Functions to detect maximal compiler version
gfortranmaxversion_4_7=\
	$(if $(findstring yes,\
	$(call gfortranmaxversion_4_6) \
	$(if $(findstring 4.7.,$(GFORTRANVERSION)),yes,no)),yes,no)
gfortranmaxversion_4_6=\
	$(if $(findstring yes,\
	$(call gfortranmaxversion_4_5) \
	$(if $(findstring 4.6.,$(GFORTRANVERSION)),yes,no)),yes,no)
gfortranmaxversion_4_5=\
	$(if $(findstring yes,\
	$(call gfortranmaxversion_4_4) \
	$(if $(findstring 4.5.,$(GFORTRANVERSION)),yes,no)),yes,no)
gfortranmaxversion_4_4=\
	$(if $(findstring yes,\
	$(call gfortranmaxversion_4_3) \
	$(if $(findstring 4.4.,$(GFORTRANVERSION)),yes,no)),yes,no)
gfortranmaxversion_4_3=\
	$(if $(findstring yes,\
	$(call gfortranmaxversion_4_2) \
	$(if $(findstring 4.3.,$(GFORTRANVERSION)),yes,no)),yes,no)
gfortranmaxversion_4_2=\
	$(if $(findstring yes,\
	$(call gfortranmaxversion_4_1) \
	$(if $(findstring 4.2.,$(GFORTRANVERSION)),yes,no)),yes,no)
gfortranmaxversion_4_1=\
	$(if $(findstring 4.1.,$(GFORTRANVERSION)),yes,no)



# Command line options for gfortran 4.3.x
ifneq (,$(findstring 4.3.,$(GFORTRANVERSION)))
CFLAGSF77LIBS := $(CFLAGSF77LIBS)
endif

# The gfortran compiler 4.1 and above supports ISO_C_BINDING 
ifeq ($(call gfortranminversion_4_1),yes)
CFLAGSF90     := -DHAS_ISO_C_BINDING $(CFLAGSF90)
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
