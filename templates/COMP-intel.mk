# -*- mode: makefile -*-

##############################################################################
# Intel Compiler suite
#
##############################################################################
COMPILERNAME = INTEL

# Default: No compiler wrapper commands
# This works both for builds of serial and parallel applications
F77       = ifort
F90       = ifort
CC        = icc
CXX       = icpc
LD        = ifort

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
CXX       = mpiCC
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
INTELVERSION := $(shell eval $(F90VERSION) | \
		  sed -n -e 'y/V/v/; /^ifort.*version/h;' -e 's/^.* \([0-9]*\.[0-9]*\.*[0-9]*\).*$$/\1/p')
ifneq ($(INTELVERSION),)
INTELVERSION_MAJOR := $(shell echo $(INTELVERSION) | cut -d. -f1)
INTELVERSION_MINOR := $(shell echo $(INTELVERSION) | cut -d. -f2)
else
INTELVERSION_MAJOR := 0
INTELVERSION_MINOR := 0
endif

# Functions to detect minimal compiler version
intelminversion = $(shell if [ $(INTELVERSION_MAJOR) -gt $(1) ] || \
	                    ([ $(INTELVERSION_MAJOR) -ge $(1) ] && \
			     [ $(INTELVERSION_MINOR) -ge $(2) ]) ; then echo yes ; else echo no ; fi)

# Functions to detect maximal compiler version
intelmaxversion = $(shell if [ $(INTELVERSION_MAJOR) -lt $(1) ] || \
	                    ([ $(INTELVERSION_MAJOR) -le $(1) ] && \
			     [ $(INTELVERSION_MINOR) -le $(2) ]) ; then echo yes ; else echo no ; fi)



##############################################################################
# Compiler flags
# (including non-architecture specific optimisation flags)
##############################################################################

# Set default type of integer variables explicitly
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -DUSE_LARGEINT -integer_size 64
endif
# $(CC) and $(CXX) do not have such a corresponding option, so we have to
# pray that they default the 'int' type properly.



# Specify -fpic for shared builds
ifneq ($(strip $(SHARED)), NO)
CFLAGSF77     := $(CFLAGSF77) -fpic
CFLAGSC       := $(CFLAGSC) -fpic
LDFLAGS_LIB   := -shared
endif



# Specify -openmp for all Intel compilers
ifeq ($(strip $(OPENMP)), YES)
CFLAGSF77     := -DUSE_OPENMP -openmp $(CFLAGSF77)
CFLAGSC       := -DUSE_OPENMP -openmp $(CFLAGSC)
LDFLAGS       := -DUSE_OPENMP -openmp $(LDFLAGS)
ifeq ($(call intelminversion,9,1),yes)
CFLAGSF77 := $(CFLAGSF77) -DHAS_OPENMP25
CFLAGSC   := $(CFLAGSC) -DHAS_OPENMP25
LDFLAGS   := $(LDFLAGS) -DHAS_OPENMP25
endif
ifeq ($(call intelminversion,11,0),yes)
CFLAGSF77 := $(CFLAGSF77) -DHAS_OPENMP30
CFLAGSC   := $(CFLAGSC) -DHAS_OPENMP30
LDFLAGS   := $(LDFLAGS) -DHAS_OPENMP30
endif
ifeq ($(call intelminversion,12,1),yes)
CFLAGSF77 := $(CFLAGSF77) -DHAS_OPENMP31
CFLAGSC   := $(CFLAGSC) -DHAS_OPENMP31
LDFLAGS   := $(LDFLAGS) -DHAS_OPENMP31
endif
ifeq ($(call intelminversion,14,0),yes)
CFLAGSF77 := $(CFLAGSF77) -DHAS_OPENMP40
CFLAGSC   := $(CFLAGSC) -DHAS_OPENMP40
LDFLAGS   := $(LDFLAGS) -DHAS_OPENMP40
endif
endif



ifneq (,$(findstring EXPENSIVE ,$(OPT)))
# Specifying -ipo for all Intel compilers is only feasible
# if all Intel compilers have the same build date!
CFLAGSF77     := -ipo $(CFLAGSF77)
CFLAGSC       := -ipo $(CFLAGSC)
LDFLAGS       := -ipo $(LDFLAGS)
endif



# Check if Fortran source code files are preprocessed by F90CPP script.
# Otherwise, set additional flags to enforce use of CPP preprocessor.
ifeq ($(strip $(F90CPP)),)
CFLAGSF77 := $(CFLAGSF77) -fpp -I$(FEAT2BASEDIR)
# Note: Do not remove trailing whitespace!
MODOPTION = -module
endif



# Treatment of trailing underscores
ifeq ($(strip $(UNDERSCORE)), YES)
CFLAGSF77     := $(CFLAGSF77) -assume underscore
else
CFLAGSF77     := $(CFLAGSF77) -DUSE_NO_UNDERSCORE -assume nounderscore
CFLAGSC       := $(CFLAGSC) -DUSE_NO_UNDERSCORE
endif




# Set default compile flags
ifeq ($(call optimise), YES)
# Don't specify -ipo here. IPO optimisation is enabled by the flag opt=expensive.
# Description of compiler flags:
#  -O3                   : enables aggressive optimization
#  -align records        : aligns derived-type components and record structure fields on default natural boundaries
#  -assume buffered_io   : data is buffered before written to disk
#  -assume underscore    : append an underscore character to external user-defined names
# or
#  -assume nounderscore  : do not append an underscore character to external user-defined names
#  -fp-model precise     : strictly adhere to value-safe  optimizations  when implementing floaing-point calculations
#  -funroll-loops        : enables loop unrolling
#  -ip                   : enables interprocedural optimizations for single-file compilation
#  -no-prec-div          : improves precision of floating-point divides
#  -opt-malloc-options=3 : enables optimal malloc algorithm
#  -pad                  : enables the changing of the variable and array memory layout.
#  -unroll               : enables aggressive loop unrolling
CFLAGSF77     := -DUSE_COMPILER_INTEL $(CFLAGSF77) -O3 \
		 -unroll-aggressive -ip -fp-model precise \
		 -no-prec-div -pad -opt-malloc-options=3
CFLAGSF90     := $(CFLAGSF90) $(CFLAGSF77) \
		 -align records -assume buffered_io
CFLAGSC       := -DUSE_COMPILER_INTEL $(CFLAGSC) -O3 -unroll-aggressive -ip -fp-model precise
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
#
# no optimisations
#
else
#
CFLAGSF77     := -DUSE_COMPILER_INTEL $(CFLAGSF77) -DDEBUG -O0 -g -fpe0
# the following (additional) flags make ifort super-strict and allow e.g.
# the detection of out-of-bounds accesses. Unfortunately, they already complain
# about such errors deep inside the numerical factorisation in UMFPACK (at
# least in version 5.5.2, which is the only one I tested), so the use of
# an iterative/alternative coarse grid solver is recommended.
#CFLAGSF77     := $(CFLAGSF77) -check all -debug -fp-stack-check -traceback -ftrapuv
CFLAGSF90     := $(CFLAGSF90) $(CFLAGSF77) \
		 -C -check bounds -traceback -warn all -assume buffered_io
CFLAGSC       := -DUSE_COMPILER_INTEL $(CFLAGSC) -DDEBUG -O0 -g
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
endif


##############################################################################
# Non-standard features supported by compiler
##############################################################################
CFLAGSF90     := -DHAS_INTRINSIC_FLUSH \
		 -DHAS_INTRINSIC_IARGC \
		 -DHAS_INTRINSIC_ISATTY \
		 -DHAS_INTRINSIC_ISNAN $(CFLAGSF90)

# The Intel compiler 10.1 and above supports ISO_C_BINDING
ifeq ($(call intelminversion,10,1),yes)
CFLAGSF90     := -DHAS_ISO_C_BINDING $(CFLAGSF90)
endif

# The Intel compiler 10.1 and above support IEEE_ARITHMETIC
# However, only the Intel compiler 11.0 and above implement intrinsic
# ieee_xxx functions as elemental. It is therefore advisable to
# require 11.0 is minimal compiler version here
ifeq ($(call intelminversion,11,0),yes)
CFLAGSF90     := -DHAS_INTRINSIC_IEEE_ARITHMETIC $(CFLAGSF90)
endif

# Enable workarounds for Intel 10.1.0[0-1][0-9] compiler releases,
# (not necessary for Intel 10.1.021)
ifneq (,$(findstring 10.1.00,$(INTELVERSION)))
CFLAGSF90     := -DUSE_COMPILER_INTEL_EARLY_10_1_WORKAROUNDS $(CFLAGSF90)
endif
ifneq (,$(findstring 10.1.01,$(INTELVERSION)))
CFLAGSF90     := -DUSE_COMPILER_INTEL_EARLY_10_1_WORKAROUNDS $(CFLAGSF90)
endif



##############################################################################
# command to create archive from object files
##############################################################################
AR        = xiar -rv


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
SBB_CVERSIONCMD = $(F77) -V  2>&1 | sed 's|(R)||g; 1!d;'


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
