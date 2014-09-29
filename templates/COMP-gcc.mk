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

# Compiler flag to specify the directory where module files should be
# placed when created and where they should be searched for.
# Note: Do not remove ticks!
MODOPTION = '-J'

# If preprocessor switch -DENABLE_SERIAL_BUILD does not occur in compiler flags,
# a build for parallel execution is requested.
ifeq (,$(findstring -DENABLE_SERIAL_BUILD ,$(APPONLYFLAGS) $(CFLAGSF90) ))
ifeq ($(strip $(MPIWRAPPERS)), YES)
F77       = mpif77
F90       = mpif90
CC        = mpicc
CXX       = mpic++
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

# Detect C/C++ compiler version
GCCVERSION := $(shell eval $(CCVERSION) | \
		sed -n -e 'y/GCV-/gcv /; /^gcc.*version/h;' -e 's/^gcc.* \([0-9]*\.[0-9]\.[0-9]\) .*$$/\1/p')
ifneq ($(GCCVERSION),)
GCCVERSION_MAJOR := $(shell echo $(GCCVERSION) | cut -d. -f1)
GCCVERSION_MINOR := $(shell echo $(GCCVERSION) | cut -d. -f2)
else
GCCVERSION_MAJOR := 0
GCCVERSION_MINOR := 0
endif

# Function to detect minimal GCC compiler version
gccminversion = $(shell if [ $(GCCVERSION_MAJOR) -gt $(1) ] || \
	                  ([ $(GCCVERSION_MAJOR) -ge $(1) ] && \
			   [ $(GCCVERSION_MINOR) -ge $(2) ]) ; then echo yes ; else echo no ; fi)

# Function to detect maximal GCC compiler version
gccmaxversion = $(shell if [ $(GCCVERSION_MAJOR) -lt $(1) ] || \
	                  ([ $(GCCVERSION_MAJOR) -le $(1) ] && \
		  	   [ $(GCCVERSION_MINOR) -le $(2) ]) ; then echo yes ; else echo no ; fi)


# Detect Fortran compiler version
GFORTRANVERSION := $(shell eval $(F90VERSION) | \
		     sed -n -e 'y/GCV-/gcv /; /^gcc.*version/h;' -e 's/^gcc.* \([0-9]*\.[0-9]\.[0-9]\) .*$$/\1/p')

ifneq ($(GFORTRANVERSION),)
GFORTRANVERSION_MAJOR := $(shell echo $(GFORTRANVERSION) | cut -d. -f1)
GFORTRANVERSION_MINOR := $(shell echo $(GFORTRANVERSION) | cut -d. -f2)
else
GFORTRANVERSION_MAJOR := 0
GFORTRANVERSION_MINOR := 0
endif

# Functions to detect minimal Fortran compiler version
gfortranminversion = $(shell if [ $(GFORTRANVERSION_MAJOR) -gt $(1) ] || \
	                       ([ $(GFORTRANVERSION_MAJOR) -ge $(1) ] && \
				[ $(GFORTRANVERSION_MINOR) -ge $(2) ]) ; then echo yes ; else echo no ; fi)

# Functions to detect maximal Fortran compiler version
gfortranmaxversion = $(shell if [ $(GFORTRANVERSION_MAJOR) -lt $(1) ] || \
	                       ([ $(GFORTRANVERSION_MAJOR) -le $(1) ] && \
				[ $(GFORTRANVERSION_MINOR) -le $(2) ]) ; then echo yes ; else echo no ; fi)



##############################################################################
# Compiler flags
# (including non-architecture specific optimisation flags)
##############################################################################

# Set default type of integer variables explicitly
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -DUSE_LARGEINT -fdefault-integer-8
endif
# $(CC) and $(CXX) do not have such a corresponding option, so we have to
# pray that they default the 'int' type properly.



# Specify -fPIC for shared builds
ifneq ($(strip $(SHARED)), NO)
CFLAGSF77     := $(CFLAGSF77) -fPIC
CFLAGSC       := $(CFLAGSC) -fPIC
LDFLAGS_LIB   := -shared
endif



# Specify -fopenmp for all gcc compilers
ifeq ($(strip $(OPENMP)), YES)
CFLAGSF77     := -DUSE_OPENMP -fopenmp $(CFLAGSF77)
CFLAGSC       := -DUSE_OPENMP -fopenmp $(CFLAGSC)
LDFLAGS       := -DUSE_OPENMP -fopenmp $(LDFLAGS)
ifeq ($(call gccminversion,4,2),yes)
CFLAGSF77 := $(CFLAGSF77) -DHAS_OPENMP25
CFLAGSC   := $(CFLAGSC) -DHAS_OPENMP25
LDFLAGS   := $(LDFLAGS) -DHAS_OPENMP25
endif
ifeq ($(call gccminversion,4,4),yes)
CFLAGSF77 := $(CFLAGSF77) -DHAS_OPENMP30
CFLAGSC   := $(CFLAGSC) -DHAS_OPENMP30
LDFLAGS   := $(LDFLAGS) -DHAS_OPENMP30
endif
ifeq ($(call gccminversion,4,7),yes)
CFLAGSF77 := $(CFLAGSF77) -DHAS_OPENMP31
CFLAGSC   := $(CFLAGSC) -DHAS_OPENMP31
LDFLAGS   := $(LDFLAGS) -DHAS_OPENMP31
endif
ifeq ($(call gccminversion,4,9),yes)
# Support of OpenMP 4.0 in gfortran comes with GCC 4.10
CFLAGSC   := $(CFLAGSC) -DHAS_OPENMP40
LDFLAGS   := $(LDFLAGS) -DHAS_OPENMP40
endif
ifeq ($(call gccminversion,4,10),yes)
CFLAGSF77 := $(CFLAGSF77) -DHAS_OPENMP40
CFLAGSC   := $(CFLAGSC) -DHAS_OPENMP40
LDFLAGS   := $(LDFLAGS) -DHAS_OPENMP40
endif
endif



# Specify -flto for all gcc compilers
ifneq (,$(findstring EXPENSIVE,$(OPT)))
ifeq ($(call gfortranminversion,4,1),yes)
CFLAGSF77     := -flto $(CFLAGSF77)
CFLAGSC       := -flto $(CFLAGSC)
LDFLAGS       := -flto $(LDFLAGS)
endif
endif



# Check if Fortran source code files are preprocessed by F90CPP script.
# Otherwise, set additional flags to enforce use of CPP preprocessor.
ifeq ($(strip $(F90CPP)),)
CFLAGSF77 := $(CFLAGSF77) -cpp -ffree-line-length-none -I$(FEAT2BASEDIR)
MODOPTION  = -J
endif



# Set default compile flags
ifeq ($(call optimise), YES)
# Don't specify -flto here. LTO optimisation is enabled by the flag opt=expensive.
# Description of compiler flags:
#  -O3                      : enables aggressive optimization
#  -ffast-math              : turns on fast math operations
#  -foptimize-register-move : attempt to reassign register numbers in move instructions and as
#                             operands of other simple instructions in order to maximize the
#                             amount of register tying. This is especially helpful on machines
#                             with two-operand instructions.
#  -fprefetch-loop-arrays   : If supported by the target machine, generate instructions to
#                             prefetch memory to improve the performance of loops that access
#                             large arrays.
#  -funroll-loops           : Unroll loops whose number of iterations can be determined at
#                             compile time or upon entry to the loop.
#  -static                  : static binary
#  -fno-second-underscore   : do not prepend second underscore
# or
#  -fno-underscoring        : do not prepend any underscore
#  -Wuninitialized          : warn uninitialised variables (from GCC 4.7 onwards warn even about
#                             possibly uninitialised variables, yielding a lot of false positives)
CFLAGSF77     := -DUSE_COMPILER_GCC $(CFLAGSF77) -O3 \
		 -ffast-math -foptimize-register-move -fprefetch-loop-arrays \
		 -funroll-loops -static
ifeq ($(strip $(UNDERSCORE)), YES)
CFLAGSF77     := $(CFLAGSF77) -fno-second-underscore
else
CFLAGSF77     := $(CFLAGSF77) -fno-underscoring
endif
CFLAGSF90     := $(CFLAGSF90) $(CFLAGSF77) -Wuninitialized
ifeq ($(call gfortranminversion,4,7),yes)
# Reduce number of false positives with GCC 4.7 and 4.8
CFLAGSF90     := $(CFLAGSF90) -Wno-maybe-uninitialized
endif
CFLAGSC       := -DUSE_COMPILER_GCC $(CFLAGSC) -O3 \
		 -ffast-math -foptimize-register-move -fprefetch-loop-arrays \
		 -funroll-loops -static
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
else
CFLAGSF77     := -DUSE_COMPILER_GCC $(CFLAGSF77) -DDEBUG -O0 -g \
		 -fbacktrace -fexternal-blas #-pg
ifeq ($(strip $(UNDERSCORE)), YES)
CFLAGSF77     := $(CFLAGSF77) -fno-second-underscore
else
CFLAGSF77     := $(CFLAGSF77) -fno-underscoring
endif
CFLAGSF90     := $(CFLAGSF90) $(CFLAGSF77) -fbounds-check \
		 -Wcharacter-truncation -Winline \
		 -Wline-truncation -Wsurprising  \
		 -Wunreachable-code -Wunused-label -Wunused-variable \
		 -Wuninitialized
		 # -Wnonstd-intrinsics: not available in 4.4.0
		 # -Wimplicit-interface -Wunused-variable
ifeq ($(call gfortranminversion,4,7),yes)
# Reduce number of false positives with GCC 4.7 and 4.8
CFLAGSF90     := $(CFLAGSF90) -Wno-maybe-uninitialized
endif
# Do not specify:
# * -std=f95
#   as it gives conflicts with LAM/MPI's mpif.h
# * -Wimplicit-interface
#   as it gives warnings for every MPI call, BLAS call etc.
CFLAGSC       := -DUSE_COMPILER_GCC $(CFLAGSC) -DDEBUG -g -fbounds-check #-pg
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS) #-pg
endif



##############################################################################
# Non-standard features supported by compiler
##############################################################################
CFLAGSF90     := -DHAS_INTRINSIC_FLUSH \
	         -DHAS_INTRINSIC_IARGC \
	         -DHAS_INTRINSIC_ISATTY \
	         -DHAS_INTRINSIC_ISNAN $(CFLAGSF90)

# The gfortran compiler 4.1 and above supports ISO_C_BINDING
ifeq ($(call gfortranminversion,4,1),yes)
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
