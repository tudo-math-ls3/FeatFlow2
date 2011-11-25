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
F77VERSION = $(F77) -V | head -n 2
F90VERSION = $(F90) -V | head -n 2
CCVERSION  = $(CC)  -V | head -n 2
CXXVERSION = $(CXX) -V | head -n 2


##############################################################################
# compiler flags 
# (including non-architecture specific optimisation flags)
##############################################################################

# Set default type of integer variables explicitly
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -DUSE_LARGEINT -i8
CFLAGSF90     := $(CFLAGSF90) -DUSE_LARGEINT -i8
endif
# $(CC) and $(CXX) do not have such a corresponding option, so we have to 
# pray that they default the 'int' type properly.



# Specify -openmp for all PGI compilers
ifeq ($(strip $(OPENMP)), YES)
CFLAGSF77LIBS := -DUSE_OPENMP -mp $(CFLAGSF77LIBS)
CFLAGSC       := -DUSE_OPENMP -mp $(CFLAGSC)
LDFLAGS       := -DUSE_OPENMP -mp $(LDFLAGS)
endif



# Set default compile flags
ifeq ($(call optimise), YES)
# -Mcache_align is important when using ACML.
CFLAGSF77LIBS := -DUSE_COMPILER_PGI $(CFLAGSF77LIBS) -O4 -fastsse \
		 -Mcray=pointer -Mcache_align -Minline=size:32 -Munroll=c:4 \
		 -Mvect=assoc,cachesize:1048576,prefetch,sse
CFLAGSF77     := $(CFLAGSF77LIBS) $(CFLAGSF77)
# PGI F90 Compiler v6.1.x most likely needs "-g -Msave" otherwise FEAT2 used to
# crash as soon as it tries to start solving something! This might have been
# fixed, though, with revision 2.5 of parallel.f90.
CFLAGSF90     := -DENABLE_USE_ONLY -DHAS_INTRINSIC_FLUSH $(CFLAGSF90) \
		 $(CFLAGSF77LIBS) -module $(OBJDIR) -I $(OBJDIR)
CFLAGSC       := -DUSE_COMPILER_PGI $(CFLAGSC) -O4 -fastsse \
		 -Mcache_align -Minline=size:32 -Munroll=c:4 \
		 -Mvect=assoc,cachesize:1048576,prefetch,sse
LDFLAGS       := $(LDFLAGS)
else
CFLAGSF77LIBS := -DUSE_COMPILER_PGI $(CFLAGSF77LIBS) -O0 -g -Mbounds
CFLAGSF77     := $(CFLAGSF77LIBS) $(CFLAGSF77)
# PGI F90 Compiler (at least 6.1.x) needs 
# * -g flag (even for -O0 optimisation level)
# otherwise FEAT2 crashes as soon as it tries to start solving something!
CFLAGSF90     := -DENABLE_USE_ONLY -DHAS_INTRINSIC_FLUSH $(CFLAGSF90) \
	         $(CFLAGSF77LIBS) -module $(OBJDIR) -I $(OBJDIR)
CFLAGSC       := -DUSE_COMPILER_PGI $(CFLAGSC) -O0 -g -B -Mbounds
LDFLAGS       := $(LDFLAGS)
endif



ifeq ($(strip $(OPT)), EXPENSIVE)
# Specify -Mipa for all PGI compilers
CFLAGSF77LIBS := -Mipa $(CFLAGSF77LIBS)
CFLAGSF77     := -Mipa $(CFLAGSF77)
CFLAGSF90     := -Mipa $(CFLAGSF90)
CFLAGSC       := -Mipa $(CFLAGSC)
CFLAGSCXX     := -Mipa $(CFLAGSCXX)
LDFLAGS       := -Mipa $(LDFLAGS)
endif



# Detect compiler version
PGIVERSION  := $(shell eval $(F90VERSION) )

# Functions to detect minimal compiler version
pgiminversion_11_10=\
	$(if $(findstring 11.10-,$(PGIVERSION)),yes,no)
pgiminversion_11_9=\
	$(if $(findstring yes,\
	$(call pgiminversion_11_10) \
	$(if $(findstring 11.9-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_11_8=\
	$(if $(findstring yes,\
	$(call pgiminversion_11_9) \
	$(if $(findstring 11.8-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_11_7=\
	$(if $(findstring yes,\
	$(call pgiminversion_11_8) \
	$(if $(findstring 11.7-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_11_6=\
	$(if $(findstring yes,\
	$(call pgiminversion_11_7) \
	$(if $(findstring 11.6-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_11_5=\
	$(if $(findstring yes,\
	$(call pgiminversion_11_6) \
	$(if $(findstring 11.5-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_11_4=\
	$(if $(findstring yes,\
	$(call pgiminversion_11_5) \
	$(if $(findstring 11.4-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_11_3=\
	$(if $(findstring yes,\
	$(call pgiminversion_11_4) \
	$(if $(findstring 11.3-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_11_2=\
	$(if $(findstring yes,\
	$(call pgiminversion_11_3) \
	$(if $(findstring 11.2-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_11_1=\
	$(if $(findstring yes,\
	$(call pgiminversion_11_2) \
	$(if $(findstring 11.1-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_11_0=\
	$(if $(findstring yes,\
	$(call pgiminversion_11_1) \
	$(if $(findstring 11.0-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_10_9=\
	$(if $(findstring yes,\
	$(call pgiminversion_11_0) \
	$(if $(findstring 10.9-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_10_8=\
	$(if $(findstring yes,\
	$(call pgiminversion_10_9) \
	$(if $(findstring 10.8-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_10_7=\
	$(if $(findstring yes,\
	$(call pgiminversion_10_8) \
	$(if $(findstring 10.7-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_10_6=\
	$(if $(findstring yes,\
	$(call pgiminversion_10_7) \
	$(if $(findstring 10.6-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_10_5=\
	$(if $(findstring yes,\
	$(call pgiminversion_10_6) \
	$(if $(findstring 10.5-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_10_4=\
	$(if $(findstring yes,\
	$(call pgiminversion_10_5) \
	$(if $(findstring 10.4-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_10_3=\
	$(if $(findstring yes,\
	$(call pgiminversion_10_4) \
	$(if $(findstring 10.3-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_10_2=\
	$(if $(findstring yes,\
	$(call pgiminversion_10_3) \
	$(if $(findstring 10.2-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_10_1=\
	$(if $(findstring yes,\
	$(call pgiminversion_10_2) \
	$(if $(findstring 10.1-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_10_0=\
	$(if $(findstring yes,\
	$(call pgiminversion_10_1) \
	$(if $(findstring 10.0-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_9_0=\
	$(if $(findstring yes,\
	$(call pgiminversion_10_0) \
	$(if $(findstring 9.0-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_8_0=\
	$(if $(findstring yes,\
	$(call pgiminversion_8_0) \
	$(if $(findstring 8.0-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_7_2=\
	$(if $(findstring yes,\
	$(call pgiminversion_8_0) \
	$(if $(findstring 7.2-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_7_1=\
	$(if $(findstring yes,\
	$(call pgiminversion_7_2) \
	$(if $(findstring 7.1-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_6_2=\
	$(if $(findstring yes,\
	$(call pgiminversion_7_1) \
	$(if $(findstring 6.2-,$(PGIVERSION)),yes,no)),yes,no)
pgiminversion_6_1=\
	$(if $(findstring yes,\
	$(call pgiminversion_6_2) \
	$(if $(findstring 6.1-,$(PGIVERSION)),yes,no)),yes,no)

# Functions to detect maximal compiler version
pgimaxversion_11_10=\
	$(if $(findstring yes,\
	$(call pgimaxversion_11_9) \
	$(if $(findstring 11.10-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_11_9=\
	$(if $(findstring yes,\
	$(call pgimaxversion_11_8) \
	$(if $(findstring 11.9-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_11_8=\
	$(if $(findstring yes,\
	$(call pgimaxversion_11_7) \
	$(if $(findstring 11.8-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_11_7=\
	$(if $(findstring yes,\
	$(call pgimaxversion_11_6) \
	$(if $(findstring 11.7-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_11_6=\
	$(if $(findstring yes,\
	$(call pgimaxversion_11_5) \
	$(if $(findstring 11.6-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_11_5=\
	$(if $(findstring yes,\
	$(call pgimaxversion_11_4) \
	$(if $(findstring 11.5-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_11_4=\
	$(if $(findstring yes,\
	$(call pgimaxversion_11_3) \
	$(if $(findstring 11.4-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_11_3=\
	$(if $(findstring yes,\
	$(call pgimaxversion_11_2) \
	$(if $(findstring 11.3-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_11_2=\
	$(if $(findstring yes,\
	$(call pgimaxversion_11_1) \
	$(if $(findstring 11.2-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_11_1=\
	$(if $(findstring yes,\
	$(call pgimaxversion_11_0) \
	$(if $(findstring 11.1-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_11_0=\
	$(if $(findstring yes,\
	$(call pgimaxversion_10_9) \
	$(if $(findstring 10.0-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_10_9=\
	$(if $(findstring yes,\
	$(call pgimaxversion_10_8) \
	$(if $(findstring 10.9-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_10_8=\
	$(if $(findstring yes,\
	$(call pgimaxversion_10_7) \
	$(if $(findstring 10.8-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_10_7=\
	$(if $(findstring yes,\
	$(call pgimaxversion_10_6) \
	$(if $(findstring 10.7-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_10_6=\
	$(if $(findstring yes,\
	$(call pgimaxversion_10_5) \
	$(if $(findstring 10.6-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_10_5=\
	$(if $(findstring yes,\
	$(call pgimaxversion_10_4) \
	$(if $(findstring 10.5-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_10_4=\
	$(if $(findstring yes,\
	$(call pgimaxversion_10_3) \
	$(if $(findstring 10.4-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_10_3=\
	$(if $(findstring yes,\
	$(call pgimaxversion_10_2) \
	$(if $(findstring 10.3-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_10_2=\
	$(if $(findstring yes,\
	$(call pgimaxversion_10_1) \
	$(if $(findstring 10.2-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_10_1=\
	$(if $(findstring yes,\
	$(call pgimaxversion_10_0) \
	$(if $(findstring 10.1-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_10_0=\
	$(if $(findstring yes,\
	$(call pgimaxversion_9_0) \
	$(if $(findstring 10.0-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_9_0=\
	$(if $(findstring yes,\
	$(call pgimaxversion_8_0) \
	$(if $(findstring 9.0-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_8_0=\
	$(if $(findstring yes,\
	$(call pgimaxversion_7_2) \
	$(if $(findstring 8.0-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_7_2=\
	$(if $(findstring yes,\
	$(call pgimaxversion_7_1) \
	$(if $(findstring 7.2-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_7_1=\
	$(if $(findstring yes,\
	$(call pgimaxversion_6_2) \
	$(if $(findstring 7.1-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_6_2=\
	$(if $(findstring yes,\
	$(call pgimaxversion_6_1) \
	$(if $(findstring 6.2-,$(PGIVERSION)),yes,no)),yes,no)
pgimaxversion_6_1=\
	$(if $(findstring 6.1-,$(PGIVERSION)),yes,no)



# The PGI compiler 7.2 and above supports ISO_C_BINDING 
ifeq ($(call pgiminversion_7_2),yes)
CFLAGSF90     := -DHAS_ISO_C_BINDING $(CFLAGSF90)
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
