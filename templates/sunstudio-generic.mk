# -*- mode: makefile -*-

##############################################################################
# Sun Studio Compiler suite
#
##############################################################################
COMPILERNAME = SUNSTUDIO

# Default: No compiler wrapper commands
# This works both for builds of serial and parallel applications
F77       = sunf77
F90       = sunf90
CC        = suncc
CXX       = sunCC
LD        = sunf90

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
F77VERSION = $(F77) -V 2>&1 1>/dev/null | head -n 1
F90VERSION = $(F90) -V 2>&1 1>/dev/null | head -n 1
CCVERSION  = $(CC)  -V 2>&1 1>/dev/null | head -n 1
CXXVERSION = $(CXX) -V 2>&1 1>/dev/null | head -n 1


##############################################################################
# compiler flags 
# (including non-architecture specific optimisation flags)
##############################################################################

# Set default type of integer variables explicitly
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -DUSE_LARGEINT -xtypemap=integer:64
CFLAGSF90     := $(CFLAGSF90) -DUSE_LARGEINT -xtypemap=integer:64
endif
# $(CC) and $(CXX) do not have such a corresponding option, so we have to 
# pray that they default the 'int' type properly.



# Specify -xopenmp for all Sun compilers
ifeq ($(strip $(OPENMP)), YES)
CFLAGSF77LIBS := -DUSE_OPENMP -xopenmp $(CFLAGSF77LIBS)
CFLAGSC       := -DUSE_OPENMP -xopenmp $(CFLAGSC)
LDFLAGS       := -DUSE_OPENMP -xopenmp $(LDFLAGS)
endif



# Set default compile flags
ifeq ($(call optimise), YES)

# MM: -fast flag produces internal compiler errors; this is a reduced
#     selection of compiler flags which would be set by the -fast flag
FAST := -libmil -dalign -xlibmopt -xdepend -pad=local -fround=nearest -xregs=frameptr -xprefetch -xvector

CFLAGSF77LIBS := -DUSE_COMPILER_SUNSTUDIO $(CFLAGSF77LIBS) $(FAST) -xtypemap=integer:32 
CFLAGSF77     := $(CFLAGSF77LIBS) $(CFLAGSF77)
CFLAGSF90     := -DHAS_INTRINSIC_FLUSH -DUSE_COMPILER_SUNSTUDIO $(CFLAGSF90) \
	         $(FAST) -moddir=$(OBJDIR)
CFLAGSC       := -DUSE_COMPILER_SUNSTUDIO $(CFLAGSC) $(FAST)
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS) $(FAST)
else
CFLAGSF77LIBS := -DUSE_COMPILER_SUNSTUDIO $(CFLAGSF77LIBS) -xtypemap=integer:32 -g -nolibmil 
CFLAGSF77     := $(CFLAGSF77LIBS) $(CFLAGSF77)
CFLAGSF90     := -DHAS_INTRINSIC_FLUSH -DUSE_COMPILER_SUNSTUDIO $(CFLAGSF90) \
		 -g -moddir=$(OBJDIR) -nolibmil \
		 #-C -xdebugformat=stabs
CFLAGSC       := -DUSE_COMPILER_SUNSTUDIO $(CFLAGSC) -g #-xdebugformat=stabs
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
endif


ifeq ($(strip $(OPT)), EXPENSIVE)
# Specifying -xipo for interprocedural optimizations
CFLAGSF77LIBS := -xipo $(CFLAGSF77LIBS)
CFLAGSF77     := -xipo $(CFLAGSF77)
CFLAGSF90     := -xipo $(CFLAGSF90)
CFLAGSC       := -xipo $(CFLAGSC)
CFLAGSCXX     := -xipo $(CFLAGSCXX)
LDFLAGS       := -xipo $(LDFLAGS)
endif


# SunStudio 10 Fortran compiler benefits when setting -DENABLE_USE_ONLY, 
# SunStudio 11 and 12 Fortran compiler, however, crash with internal compiler 
# errors with this setting, both for unoptimised and optimised builds.
#CFLAGSF90 := -DENABLE_USE_ONLY $(CFLAGSF90)

# Detect compiler version
SUNSTUDIOVERSION  := $(shell eval $(CXXVERSION) )

# Enable workarounds for Intel 10.1.0[0-1][0-9] compiler releases, 
# (not necessary for Intel 10.1.021)
ifneq (,$(findstring Sun C++ 5.10,$(SUNSTUDIOVERSION)))
ifneq (,$(findstring 2009/03/06,$(SUNSTUDIOVERSION)))
CFLAGSF90     := -DUSE_COMPILER_SUNSTUDIO_12_1_OR_PRERELEASE $(CFLAGSF90)
endif
ifneq (,$(findstring 2009/06/03,$(SUNSTUDIOVERSION)))
CFLAGSF90     := -DUSE_COMPILER_SUNSTUDIO_12_1_OR_PRERELEASE $(CFLAGSF90)
endif
endif
ifneq (,$(findstring Sun C++ 5.11,$(SUNSTUDIOVERSION)))
ifneq (,$(findstring 2010/05/10,$(SUNSTUDIOVERSION)))
CFLAGSF90     := -DUSE_COMPILER_SUNSTUDIO_12_2_OR_PRERELEASE $(CFLAGSF90)
endif
ifneq (,$(findstring 2010/08/13,$(SUNSTUDIOVERSION)))
CFLAGSF90     := -DUSE_COMPILER_SUNSTUDIO_12_2_OR_PRERELEASE $(CFLAGSF90)
endif
endif

# The Sun fortran compiler supports ISO_C_BINDING
CFLAGSF90     := -DHAS_ISO_C_BINDING $(CFLAGSF90)



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
SBB_CVERSIONCMD = $(F77) -V 2>&1 | sed '2!d;'


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
