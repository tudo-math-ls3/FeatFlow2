# -*- mode: makefile -*-

##############################################################################
# G95 (www.g95.org) and GNU C/C++ compiler 3.x or 4.x
#
##############################################################################
COMPILERNAME = G95

# Default: No compiler wrapper commands
# This works both for builds of serial and parallel applications
F77       = g95
F90       = g95
CC        = gcc
CXX       = g++
LD        = g95

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
# All integer variables in FEAST are explicitly typed to either 32 or 64 bits.
# The only native integers are literals and code written in F77 (blas, lapack, 
# sbblas) and C (metis, coproc backend). FEAST assumes that all these (native) 
# integers are 32-bit!!!
# So, to get things running with compilers that do not default native integers
# to 32 bits, we need to add an appropriate compiler flag to
# CFLAGSF77LIBS: -i4. 
# This also applies when changing the kind-values in kernel/fsystem.f90.

# $(CC) and $(CXX) do not have such a corresponding option, so we have to 
# pray that they default the 'int' type properly.

ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -DUSE_COMPILER_G95 $(CFLAGSF77LIBS) -O3 \
		 -ffast-math -foptimize-register-move \
		 -fprefetch-loop-arrays -funroll-loops -static \
		 -fno-second-underscore -i4
CFLAGSF77     := $(CFLAGSF77LIBS) $(CFLAGSF77)
CFLAGSF90     := -DENABLE_USE_ONLY -DHAS_INTRINSIC_FLUSH \
	         $(CFLAGSF90) $(CFLAGSF77LIBS) -fmod=$(OBJDIR)
CFLAGSC       := -DUSE_COMPILER_G95 $(CFLAGSC) -O3 \
		 -ffast-math -foptimize-register-move -fprefetch-loop-arrays \
		 -funroll-loops -static 
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
else
CFLAGSF77LIBS := -DUSE_COMPILER_G95 $(CFLAGSF77LIBS) -g -fno-second-underscore \
		-i4 #-pg
# Don't include "-fbounds-check" in CFLAGSF77 as the SBBLAS routines
# are peppered with out-of-bounds accesses (i[0] and i[n]). That's 
# Turek-style code and nobody volunteered so far to fix it.
CFLAGSF77     := $(CFLAGSF77LIBS) $(CFLAGSF77) # -fbounds-check
CFLAGSF90     := -DENABLE_USE_ONLY -DHAS_INTRINSIC_FLUSH \
	         $(CFLAGSF90) $(CFLAGSF77LIBS) -fmod=$(OBJDIR) \
	         -fbounds-check -ftrace=full
CFLAGSC       := -DUSE_COMPILER_G95 $(CFLAGSC) -g -fbounds-check #-pg
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS) #-pg
endif


# Workarounds for known problems:

# Workaround 1:
# * g95 v0.9 does not provide an intrinsic isatty function, 
#   every daily snapshot since january 2006 does.
#
#   Workaround: provide own interface to C backend isatty 
#               for g95 versions that do not provide an isatty symbol.
#   Only test, if g95 is found in path.
ifneq ($(shell (which 2>/dev/null $(F90))),)
ifneq ($(shell (nm `$(F90) -print-libgcc-file-name | sed 's|/libgcc.a|libf95.a|;'` | grep 'T _g95_isatty')),)
CFLAGSF90     := -DHAS_INTRINSIC_ISATTY $(CFLAGSF90)
endif
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


# The settings needed to compile a FEAST application are "wildly" distributed
# over several files ((Makefile.inc and templates/*.mk) and if-branches 
# (in an attempt to reduce duplicate code and inconsistencies among all build 
# IDs that e.g. use the same MPI environment). Not having all settings at 
# *one* place entails the risk (especially in the event of setting up a new 
# build ID) that settings are incompletely defined. A simple typo in a matching 
# rule in Makefile.inc may prevent that the compiler and compiler command line
# flags are set. Compilation would fail with the most peculiar errors - if not
# the Makefile had been set up to catch such a case.
# Each build ID in FEAST has 6 tokens: architecture, cpu, operating system,
# compiler family, BLAS implementation, MPI environment. Whenever setting
# one of these, an according flag is set. They are named TOKEN1 up to TOKEN6.
# Before starting to actually compile a FEAST application, every Makefile
# generated by bin/configure checks whether all six tokens are set *for the 
# choosen build ID*. If not, the user gets an error message describing exactly 
# what information is missing, e.g. token 5 not set which means there is no
# information available which BLAS implementation to use and where to find the
# library.
#
# In this file, the fourth token has been set: compiler and compiler command
# line options. Set the flag accordingly.
TOKEN4 := 1
