# -*- mode: makefile -*-

##############################################################################
# IBM XL Compiler suite
#
##############################################################################
COMPILERNAME = IBMXL

# Default: No compiler wrapper commands
# This works both for builds of serial and parallel applications
F77       = xlf_r
F90       = xlf90_r
CC	  = xlc_r
CXX       = xlC_r
LD        = xlf90_r

# Compiler flag to specify the directory where module files should be
# placed when created and where they should be searched for.
# Note: Do not remove ticks and whitespace!
MODOPTION = '-qmoddir='

# If preprocessor switch -DENABLE_SERIAL_BUILD does not occur in compiler flags,
# a build for parallel execution is requested.
ifeq (,$(findstring -DENABLE_SERIAL_BUILD ,$(APPONLYFLAGS) $(CFLAGSF90) ))
ifeq ($(strip $(MPIWRAPPERS)), YES)
F77       = mpxlf
F90       = mpxlf95
CC	  = xlc
CXX       = xlC
LD        = mpxlf95
endif
endif


##############################################################################
# Commands to get version information from compiler
##############################################################################
F77VERSION = $(F77) -qversion | head -n 1
F90VERSION = $(F90) -qversion | head -n 1
CCVERSION  = $(CC)  -qversion | head -n 1
CXXVERSION = $(CXX) -qversion | head -n 1


##############################################################################
# compiler flags
# (including non-architecture specific optimisation flags)
##############################################################################

# Set default compile flags
ifeq ($(call optimise), YES)
# !!! WARNING WARNING WARNING !!!
# You must not pass any -D option to xlf compilers, because it is *not*
# interpreted by the compiler like C-style -D macro definitions. You won't
# be able to link those object files together!
# You can, though, pass C-style -D macro definitions to CFLAGSF90 if and if
# only they are interpreted by $(F90CPP) and not passed to the underlying xlf
# compiler.
# !!! WARNING WARNING WARNING !!!
#
# Info:
# -qipa defaults to (among other options) -qipa=partition=medium.
# For stokes and disk application, that gives a linking error.
# -qipa=partition=large as well.
# -qipa=partition=small increases a bit memory consumption
# (in comparison to just omitting it) which leads to some tests
# on high levels failing with -qipa and finishing without -qipa.
CFLAGSF77     := $(CFLAGSF77) \
		 -qfullpath -q64 -qstrict \
		 -qalign=struct=natural -qarch=auto -qcache=auto \
		 -qtune=auto -qunroll=auto -qsclk=micro -qhot \
		 -qipa=partition=small -Q -O5
CFLAGSF90     := -DUSE_COMPILER_XLF $(CFLAGSF90) $(CFLAGSF77) -O5 -qfree=f90 \
		 -qsuffix=f=f90 -qmaxmem=-1
CFLAGSC       := $(CFLAGSC) -qcache=auto -qtune=auto -O3 -qstrict \
		 -qunroll=auto -q64
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS) -q64 -qmaxmem=-1 -qsclk=micro
else
CFLAGSF77     := $(CFLAGSF77) -DDEBUG \
		 -g -qfullpath -q64 -O0
CFLAGSF90     := -DUSE_COMPILER_XLF $(CFLAGSF90) $(CFLAGSF77) -qfree=f90 \
                 -qsuffix=f=f90
CFLAGSC       := $(CFLAGSF77) $(CFLAGSC)
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS) -g -qfullpath -q64
# !!! WARNING WARNING WARNING !!!
# You must not pass any -D option to xlf compilers, because it is *not*
# interpreted by the compiler like C-style -D macro definitions. You won't
# be able to link those object files together!
# You can, though, pass C-style -D macro definitions to CFLAGSF90 if and if
# only they are interpreted by $(F90CPP) and not passed to the underlying xlf
# compiler.
# !!! WARNING WARNING WARNING !!!
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
MOVEMOD   = YES


##############################################################################
# Commands needed by the Sparse Banded Blas benchmark
##############################################################################
SBB_CVERSIONCMD  = $(F77) -qversion | sed '1!d;'


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
