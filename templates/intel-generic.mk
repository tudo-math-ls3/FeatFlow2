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
F77VERSION = $(F77) -V 2>&1 | head -n 1
F90VERSION = $(F90) -V 2>&1 | head -n 1
CCVERSION  = $(CC)  -V 2>&1 | head -n 1
CXXVERSION = $(CXX) -V 2>&1 | head -n 1


##############################################################################
# compiler flags
# (including non-architecture specific optimisation flags)
##############################################################################

# Set default type of integer variables explicitly
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -integer_size 64
CFLAGSF90     := $(CFLAGSF90) -integer_size 64
endif
# $(CC) and $(CXX) do not have such a corresponding option, so we have to 
# pray that they default the 'int' type properly.



# Specify -openmp for all Intel compilers
ifeq ($(strip $(OPENMP)), YES)
CFLAGSF77LIBS := -DUSE_OPENMP -openmp $(CFLAGSF77LIBS)
CFLAGSC       := -DUSE_OPENMP -openmp $(CFLAGSC)
LDFLAGS       := -DUSE_OPENMP -openmp $(LDFLAGS)
endif



# Set default compile flags
ifeq ($(call optimise), YES)
# Don't specify -ipo here. IPO optimisation is enabled by the flag opt=expensive.
CFLAGSF77LIBS := -DUSE_COMPILER_INTEL $(CFLAGSF77LIBS) -O3 \
		 -funroll-loops -ip -assume underscore \
                 -fp-model precise \
		 -no-prec-div -pad 
CFLAGSF77     := $(CFLAGSF77LIBS) $(CFLAGSF77)
CFLAGSF90     := -DHAS_INTRINSIC_FLUSH $(CFLAGSF90) $(CFLAGSF77LIBS) \
		 -module $(OBJDIR) -align records -assume buffered_io
CFLAGSC       := -DUSE_COMPILER_INTEL $(CFLAGSC) -O3 -unroll -ip -fp-model precise
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS) 
else
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -DUSE_COMPILER_INTEL -O0 -g -fpe0 -assume underscore
CFLAGSF77     := $(CFLAGSF77LIBS) $(CFLAGSF77)
CFLAGSF90     := -DHAS_INTRINSIC_FLUSH $(CFLAGSF90) $(CFLAGSF77LIBS) \
		 -module $(OBJDIR) -C -check bounds -traceback -warn all 
CFLAGSC       := -DUSE_COMPILER_INTEL $(CFLAGSC) -O0 -g
CFLAGSCXX     := $(CFLAGSC) $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS) 
endif



ifeq ($(strip $(OPT)), EXPENSIVE)
# Specifying -ipo for all Intel compilers is only feasible
# if all Intel compilers have the same build date!
CFLAGSF77LIBS := -ipo $(CFLAGSF77LIBS)
CFLAGSF77     := -ipo $(CFLAGSF77)
CFLAGSF90     := -ipo $(CFLAGSF90)
CFLAGSC       := -ipo $(CFLAGSC)
CFLAGSCXX     := -ipo $(CFLAGSCXX)
LDFLAGS       := -ipo $(LDFLAGS)
endif



# Detect compiler version
INTELVERSION  := $(shell eval $(CXXVERSION) )

# Enable workarounds for Intel 10.1.0[0-1][0-9] compiler releases, 
# (not necessary for Intel 10.1.021)
ifneq (,$(findstring 10.1.00,$(INTELVERSION)))
CFLAGSF90     := -DUSE_COMPILER_INTEL_EARLY_10_1_WORKAROUNDS $(CFLAGSF90)
endif
ifneq (,$(findstring 10.1.01,$(INTELVERSION)))
CFLAGSF90     := -DUSE_COMPILER_INTEL_EARLY_10_1_WORKAROUNDS $(CFLAGSF90)
endif



# Functions to detect minimal compiler version
intelminversion_12_0=\
	$(if $(findstring 12.0.,$(INTELVERSION)),yes,no)
intelminversion_11_1=\
	$(if $(findstring yes,\
	   $(call intelminversion_12_0) \
	      $(if $(findstring 11.1.,$(INTELVERSION)),yes,no)),yes,no)
intelminversion_11_0=\
	$(if $(findstring yes,\
	   $(call intelminversion_11_1) \
	      $(if $(findstring 11.0.,$(INTELVERSION)),yes,no)),yes,no)
intelminversion_10_1=\
	$(if $(findstring yes,\
	   $(call intelminversion_11_0) \
	      $(if $(findstring 10.1.,$(INTELVERSION)),yes,no)),yes,no)
intelminversion_10_0=\
	$(if $(findstring yes,\
	   $(call intelminversion_10_1) \
	      $(if $(findstring 10.0.,$(INTELVERSION)),yes,no)),yes,no)
intelminversion_9_1=\
	$(if $(findstring yes,\
	   $(call intelminversion_10_0) \
	      $(if $(findstring 9.1.,$(INTELVERSION)),yes,no)),yes,no)
intelminversion_9_0=\
	$(if $(findstring yes,\
	   $(call intelminversion_9_1) \
	      $(if $(findstring 9.0.,$(INTELVERSION)),yes,no)),yes,no)
intelminversion_8_1=\
	$(if $(findstring yes,\
	   $(call intelminversion_9_0) \
	      $(if $(findstring 8.1.,$(INTELVERSION)),yes,no)),yes,no)
intelminversion_8_0=\
	$(if $(findstring yes,\
	   $(call intelminversion_8_1) \
	      $(if $(findstring 8.0.,$(INTELVERSION)),yes,no)),yes,no)
intelminversion_7_1=\
	$(if $(findstring yes,\
	   $(call intelminversion_8_0) \
	      $(if $(findstring 7.1.,$(INTELVERSION)),yes,no)),yes,no)
intelminversion_7_0=\
	$(if $(findstring yes,\
	   $(call intelminversion_7_1) \
	      $(if $(findstring 7.0.,$(INTELVERSION)),yes,no)),yes,no)
intelminversion_6_0=\
	$(if $(findstring yes,\
	   $(call intelminversion_7_0) \
	      $(if $(findstring 6.0.,$(INTELVERSION)),yes,no)),yes,no)

# Functions to detect maximal compiler version
intelmaxversion_12_0=\
	$(if $(findstring yes,\
	   $(call intelmaxversion_11_1) \
	      $(if $(findstring 12.0.,$(INTELVERSION)),yes,no)),yes,no)
intelmaxversion_11_1=\
	$(if $(findstring yes,\
	   $(call intelmaxversion_11_0) \
	      $(if $(findstring 11.1.,$(INTELVERSION)),yes,no)),yes,no)
intelmaxversion_11_0=\
	$(if $(findstring yes,\
	   $(call intelmaxversion_10_1) \
	      $(if $(findstring 11.0.,$(INTELVERSION)),yes,no)),yes,no)
intelmaxversion_10_1=\
	$(if $(findstring yes,\
	   $(call intelmaxversion_10_0) \
	      $(if $(findstring 10.1.,$(INTELVERSION)),yes,no)),yes,no)
intelmaxversion_10_0=\
	$(if $(findstring yes,\
	   $(call intelmaxversion_9_1) \
	      $(if $(findstring 10.0.,$(INTELVERSION)),yes,no)),yes,no)
intelmaxversion_9_1=\
	$(if $(findstring yes,\
	   $(call intelmaxversion_9_0) \
	      $(if $(findstring 9.1.,$(INTELVERSION)),yes,no)),yes,no)
intelmaxversion_9_0=\
	$(if $(findstring yes,\
	   $(call intelmaxversion_8_1) \
	      $(if $(findstring 9.0.,$(INTELVERSION)),yes,no)),yes,no)
intelmaxversion_8_1=\
	$(if $(findstring yes,\
	   $(call intelmaxversion_8_0) \
	      $(if $(findstring 8.1.,$(INTELVERSION)),yes,no)),yes,no)
intelmaxversion_8_0=\
	$(if $(findstring yes,\
	   $(call intelmaxversion_7_1) \
	      $(if $(findstring 8.0.,$(INTELVERSION)),yes,no)),yes,no)
intelmaxversion_7_1=\
	$(if $(findstring yes,\
	   $(call intelmaxversion_7_0) \
	      $(if $(findstring 7.1.,$(INTELVERSION)),yes,no)),yes,no)
intelmaxversion_7_0=\
	$(if $(findstring yes,\
	   $(call intelmaxversion_6_0) \
	      $(if $(findstring 7.0.,$(INTELVERSION)),yes,no)),yes,no)
intelmaxversion_6_0=\
	$(if $(findstring 6.0.,$(INTELVERSION)),yes,no)



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
