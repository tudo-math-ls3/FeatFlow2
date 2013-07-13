# -*- mode: makefile -*-

##############################################################################
# Pathscale compiler suite
#
##############################################################################
COMPILERNAME = PATHSCALE

# Default: No compiler wrapper commands
# This works both for builds of serial and parallel applications
F77       = pathf90
F90       = pathf90
CC        = pathcc
CXX       = pathCC
LD        = pathf90

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
CXX	  = mpiCC
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
PATHSCALEVERSION := $(shell eval $(F90VERSION) | \
		      sed -n -e 's/^.*Version \([0-9]*\.[0-9]*\.*[0-9]*\).*$$/\1/p; q;')
ifneq ($(PATHSCALEVERSION),)
PATHSCALEVERSION_MAJOR := $(shell echo $(PATHSCALEVERSION) | cut -d. -f1)
PATHSCALEVERSION_MINOR := $(shell echo $(PATHSCALEVERSION) | cut -d. -f2)
else
PATHSCALEVERSION_MAJOR := 0
PATHSCALEVERSION_MINOR := 0
endif

# Functions to detect minimal compiler version
pathscaleminversion = $(shell if [ $(PATHSCALEVERSION_MAJOR) -gt $(1) ] || \
	                        ([ $(PATHSCALEVERSION_MAJOR) -ge $(1) ] && \
				 [ $(PATHSCALEVERSION_MINOR) -ge $(2) ]) ; then echo yes ; else echo no ; fi)

# Functions to detect maximal compiler version
pathscalemaxversion = $(shell if [ $(PATHSCALEVERSION_MAJOR) -lt $(1) ] || \
	                        ([ $(PATHSCALEVERSION_MAJOR) -le $(1) ] &&\
				 [ $(PATHSCALEVERSION_MINOR) -le $(2) ]) ; then echo yes ; else echo no ; fi)


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



# Specify -fopenmp for all Pathscale compilers
ifeq ($(strip $(OPENMP)), YES)
CFLAGSF77     := -DUSE_OPENMP -mp $(CFLAGSF77)
CFLAGSC       := -DUSE_OPENMP -mp $(CFLAGSC)
LDFLAGS       := -DUSE_OPENMP -mp $(LDFLAGS)
endif



ifeq ($(strip $(OPT)), EXPENSIVE)
# increases link time to 2 minutes per application, but is well worth it.
CFLAGSF77     := -ipa $(CFLAGSF77)
CFLAGSC       := -ipa $(CFLAGSC)
LDFLAGS       := -ipa $(LDFLAGS)
endif



# Set default compile flags
ifeq ($(call optimise), YES)
CFLAGSF77     := -DUSE_COMPILER_PATHSCALE $(CFLAGSF77) -O3 -OPT:Ofast \
		 -fno-math-errno #-Wuninitialized
CFLAGSF90     := $(CFLAGSF90) $(CFLAGSF77)
CFLAGSC       := -DUSE_COMPILER_PATHSCALE $(CFLAGSC) -O3 -OPT:Ofast \
		 -fno-math-errno
LDFLAGS       := $(LDFLAGS)
else
CFLAGSF77     := -DUSE_COMPILER_PATHSCALE $(CFLAGSF77) -g \
	         -Wuninitialized -ffortran-bounds-check
CFLAGSF90     := $(CFLAGSF90) $(CFLAGSF77)
CFLAGSC       := -DUSE_COMPILER_PATHSCALE $(CFLAGSC) -g
LDFLAGS       := $(LDFLAGS)
endif
# The option -Wuninitialized is highly recommended to let the compiler
# warn about variables that are accessed before being uninitialised.
# As the options doubles compilation time, it is not turned on by default.

# Pathscale F90 Compiler (v. 2.4) needs the -g flag (even for -O0
# optimisation level), otherwise FEAT2 crashes as soon as it tries to
# start solving something. This is no longer an issue with v. 3.1.



##############################################################################
# Non-standard features supported by compiler
##############################################################################
CFLAGSF90     := -DHAS_INTRINSIC_FLUSH \
	         -DHAS_INTRINSIC_ISATTY \
	         -DHAS_INTRINSIC_ISNAN \
		 -DHAS_INTRINSIC_IEEE_ARITHMETIC $(CFLAGSF90)

# The PathScale compiler 3.2 and above supports ISO_C_BINDING
ifeq ($(call pathscaleminversion,3,2),yes)
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
