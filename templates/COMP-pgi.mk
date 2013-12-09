# -*- mode: makefile -*-

##############################################################################
# Portland Group compiler suite
#
##############################################################################
COMPILERNAME = PGI

# Default: No compiler wrapper commands
# This works both for builds of serial and parallel applications
F77       = pgf77
F90       = pgf95
CC        = pgcc
CXX	  = pgCC
LD        = pgf95

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
F77VERSION = $(F77) -V
F90VERSION = $(F90) -V
CCVERSION  = $(CC)  -V
CXXVERSION = $(CXX) -V

# Detect compiler version
PGIVERSION := $(shell eval $(F90VERSION) | \
		sed -n -e '/^pgf90 .*target/h;' -e 's/^.* \([0-9]*\.[0-9]*-[0-9]*\) .*$$/\1/p')
ifneq ($(PGIVERSION),)
PGIVERSION_MAJOR := $(shell echo $(PGIVERSION) | cut -d. -f1)
PGIVERSION_MINOR := $(shell echo $(PGIVERSION) | cut -d. -f2 | cut -d- -f1)
else
PGIVERSION_MAJOR := 0
PGIVERSION_MINOR := 0
endif

# Functions to detect minimal compiler version
pgiminversion = $(shell if [ $(PGIVERSION_MAJOR) -gt $(1) ] || \
	                  ([ $(PGIVERSION_MAJOR) -ge $(1) ] && \
			   [ $(PGIVERSION_MINOR) -ge $(2) ]) ; then echo yes ; else echo no ; fi)

# Functions to detect maximal compiler version
pgimaxversion = $(shell if [ $(PGIVERSION_MAJOR) -lt $(1) ] || \
	                  ([ $(PGIVERSION_MAJOR) -le $(1) ] && \
			   [ $(PGIVERSION_MINOR) -le $(2) ]) ; then echo yes ; else echo no ; fi)


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



# Specify -openmp for all PGI compilers
ifeq ($(strip $(OPENMP)), YES)
CFLAGSF77     := -DUSE_OPENMP -mp $(CFLAGSF77)
CFLAGSC       := -DUSE_OPENMP -mp $(CFLAGSC)
LDFLAGS       := -DUSE_OPENMP -mp $(LDFLAGS)
endif



ifneq (,$(findstring EXPENSIVE ,$(OPT)))
# Specify -Mipa for all PGI compilers
CFLAGSF77     := -Mipa $(CFLAGSF77)
CFLAGSC       := -Mipa $(CFLAGSC)
LDFLAGS       := -Mipa $(LDFLAGS)
endif



# Check if Fortran source code files are preprocessed by F90CPP script.
# Otherwise, set additional flags to enforce use of CPP preprocessor.
ifeq ($(strip $(F90CPP)),)
CFLAGSF77 := $(CFLAGSF77) -Mpreprocess -I$(FEAT2BASEDIR)
# Note: Do not remove trailing whitespace!
MODOPTION = -module 
endif



# Set default compile flags
ifeq ($(call optimise), YES)
# -Mcache_align is important when using ACML.
CFLAGSF77     := -DUSE_COMPILER_PGI $(CFLAGSF77) -O4 -fastsse \
		 -Mcray=pointer -Mcache_align -Minline=size:32 -Munroll=c:4 \
		 -Mvect=assoc,prefetch,sse
# PGI F90 Compiler v6.1.x most likely needs "-g -Msave" otherwise FEAT2 used to
# crash as soon as it tries to start solving something! This might have been
# fixed, though, with revision 2.5 of parallel.f90.
CFLAGSF90     := $(CFLAGSF90) $(CFLAGSF77)
CFLAGSC       := -DUSE_COMPILER_PGI $(CFLAGSC) -O4 -fastsse \
		 -Mcache_align -Minline=size:32 -Munroll=c:4 \
		 -Mvect=assoc,prefetch,sse
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
else
CFLAGSF77     := -DUSE_COMPILER_PGI $(CFLAGSF77) -DDEBUG -O0 -g -Mbounds
# PGI F90 Compiler (at least 6.1.x) needs
# * -g flag (even for -O0 optimisation level)
# otherwise FEAT2 crashes as soon as it tries to start solving something!
CFLAGSF90     := $(CFLAGSF90) $(CFLAGSF77)
CFLAGSC       := -DUSE_COMPILER_PGI $(CFLAGSC) -DDEBUG -O0 -g -B -Mbounds
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
endif



##############################################################################
# Non-standard features supported by compiler
##############################################################################
CFLAGSF90     := -DHAS_INTRINSIC_FLUSH \
	         -DHAS_INTRINSIC_ISATTY $(CFLAGSF90)

# The PGI compiler 7.2 and above supports ISO_C_BINDING
ifeq ($(call pgiminversion,7,2),yes)
CFLAGSF90     := -DHAS_ISO_C_BINDING $(CFLAGSF90)
endif

# The PGI compiler 10.0 and above supports IEEE_ARITHMETIC
ifeq ($(call pgiminversion,10,0),yes)
CFLAGSF90     := -DHAS_INTRINSIC_IEEE_ARITHMETIC $(CFLAGSF90)
endif

# Enable workarounds for PGI 6.1 compiler
ifneq (,$(findstring pgf90 6.1-,$(PGIVERSION)))
CFLAGSF90     := -DUSE_COMPILER_PGI_6_1 $(CFLAGSF90)
endif

# Enable workarounds for PGI 6.2 compiler
ifneq (,$(findstring pgf90 6.2-,$(PGIVERSION)))
CFLAGSF90     := -DUSE_COMPILER_PGI_6_2 $(CFLAGSF90)
endif

# Enable workarounds for PGI 7.0 compiler
ifneq (,$(findstring pgf95 7.0-,$(PGIVERSION)))
CFLAGSF90     := -DUSE_COMPILER_PGI_7_0 $(CFLAGSF90)
endif

# Enable workarounds for PGI 7.2 compiler
ifneq (,$(findstring pgf95 7.2-,$(PGIVERSION)))
CFLAGSF90     := -DUSE_COMPILER_PGI_7_2 $(CFLAGSF90)
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
SBB_CVERSIONCMD = $(F77) -V  2>&1 | sed 's|(R)||g; 2!d;'


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
