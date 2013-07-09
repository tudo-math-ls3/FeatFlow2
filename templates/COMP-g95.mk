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

# Compiler flag to specify the directory where module files should be
# placed when created and where they should be searched for.
# Note: Do not remove ticks!
MODOPTION = '-fmod='

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
F77VERSION := $(F77) -v 2>&1
F90VERSION := $(F90) -v 2>&1
CCVERSION  := $(CC)  -v 2>&1
CXXVERSION := $(CXX) -v 2>&1

# Detect C/C++ compiler version
GCCVERSION := $(shell eval $(CCVERSION) | \
		sed -n -e 'y/GCV-/gcv /; /^gcc.*version/h;' -e 's/^.* \([0-9]*\.[0-9]\.[0-9]\) .*$$/\1/p')
ifneq ($(GCCVERSION),)
GCCVERSION_MAJOR := $(shell echo $(GCCVERSION) | cut -d. -f1)
GCCVERSION_MINOR := $(shell echo $(GCCVERSION) | cut -d. -f2)
else
GCCVERSION_MAJOR := 0
GCCVERSION_MINOR := 0
endif

# Functions to detect minimal Fortran compiler version
gccminversion = $(shell if [ $(GCCVERSION_MAJOR) -gt $(1) ] || \
	                  ([ $(GCCVERSION_MAJOR) -ge $(1) ] && \
			   [ $(GCCVERSION_MINOR) -ge $(2) ]) ; then echo yes ; else echo no ; fi)

# Functions to detect maximal Fortran compiler version
gccmaxversion = $(shell if [ $(GCCVERSION_MAJOR) -lt $(1) ] || \
	                  ([ $(GCCVERSION_MAJOR) -le $(1) ] && \
		   	   [ $(GCCVERSION_MINOR) -le $(2) ]) ; then echo yes ; else echo no ; fi)


# Detect compiler version
G95VERSION := $(shell eval $(F90VERSION) | \
		sed -n -e 'y/GCV-/gcv /; /^gcc.*version/h;' \
	               -e 's/^.*g95 \([0-9]*\.[0-9]*\.*[0-9]*\)[^0-9a-z.].*$$/\1/p')
ifneq ($(G95VERSION),)
G95VERSION_MAJOR := $(shell echo $(G95VERSION) | cut -d. -f1)
G95VERSION_MINOR := $(shell echo $(G95VERSION) | cut -d. -f2)
else
G95VERSION_MAJOR := 0
G95VERSION_MINOR := 0
endif

# Functions to detect minimal compiler version
g95minversion = $(shell if [ $(G95VERSION_MAJOR) -gt $(1) ] || \
	                  ([ $(G95VERSION_MAJOR) -ge $(1) ] && [ $(G95VERSION_MINOR) -ge $(2) ]) ; then echo yes ; else echo no ; fi)

# Functions to detect maximal compiler version
g95maxversion = $(shell if [ $(G95VERSION_MAJOR) -lt $(1) ] || \
	                  ([ $(G95VERSION_MAJOR) -le $(1) ] && [ $(G95VERSION_MINOR) -le $(2) ]) ; then echo yes ; else echo no ; fi)



##############################################################################
# compiler flags
# (including non-architecture specific optimisation flags)
##############################################################################

# Set default type of integer variables explicitly
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -DUSE_LARGEINT -i8
endif
# $(CC) and $(CXX) do not have such a corresponding option, so we have to
# pray that they default the 'int' type properly.



# Fortran compiler does not support OpenMP
ifeq ($(strip $(OPENMP)), YES)
MESSAGE  := $(MESSAGE) \
            echo; \
            echo '*** Warning: OpenMP is not supported by G95 compiler'; \
            echo;
CFLAGSC       := -DUSE_OPENMP -fopenmp $(CFLAGSC)
LDFLAGS       := -DUSE_OPENMP -fopenmp $(LDFLAGS)
endif



# Set default compile flags
ifeq ($(call optimise), YES)
CFLAGSF77     := -DUSE_COMPILER_G95 $(CFLAGSF77) -O3 \
		 -ffast-math -foptimize-register-move \
		 -fprefetch-loop-arrays -funroll-loops -static \
		 -fno-second-underscore
CFLAGSF90     := $(CFLAGSF90) $(CFLAGSF77)
CFLAGSC       := -DUSE_COMPILER_G95 $(CFLAGSC) -O3 \
		 -ffast-math -foptimize-register-move -fprefetch-loop-arrays \
		 -funroll-loops -static
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
else
CFLAGSF77     := -DUSE_COMPILER_G95 $(CFLAGSF77) -g -fno-second-underscore #-pg
# Don't include "-fbounds-check" in CFLAGSF77 as the SBBLAS routines
# are peppered with out-of-bounds accesses (i[0] and i[n]). That's
# Turek-style code and nobody volunteered so far to fix it.
CFLAGSF90     := $(CFLAGSF90) $(CFLAGSF77) -fbounds-check -ftrace=full
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
ifneq ($(shell (nm `$(F90) -print-libgcc-file-name | sed 's|/libgcc.a|libf95.a|;'` | grep 'T _g95_isnan')),)
CFLAGSF90     := -DHAS_INTRINSIC_ISNAN $(CFLAGSF90)
endif
ifneq ($(shell (nm `$(F90) -print-libgcc-file-name | sed 's|/libgcc.a|libf95.a|;'` | grep 'T _g95_iargc')),)
CFLAGSF90     := -DHAS_INTRINSIC_IARGC $(CFLAGSF90)
endif
endif


##############################################################################
# Non-standard features supported by compiler
##############################################################################
CFLAGSF90 := -DHAS_INTRINSIC_FLUSH $(CFLAGSF90)

# The gfortran compiler 4.1 and above supports ISO_C_BINDING
ifeq ($(call gfortranminversion,4,1),yes)
CFLAGSF90 := -DHAS_ISO_C_BINDING $(CFLAGSF90)
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
