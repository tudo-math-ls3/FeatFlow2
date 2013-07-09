# -*- mode: makefile -*-

##############################################################################
# NEC Compiler suite
#
##############################################################################
COMPILERNAME = NEC

# Default: No compiler wrapper commands
# This works both for builds of serial and parallel applications
F77       = sxf90
F90       = sxf90
CC        = sxcc
CXX       = echo
LD        = sxf90

# Compiler flag to specify the directory where module files should be
# placed when created and where they should be searched for.
# Note: Do not remove ticks and whitespace!
MODOPTION = '-to '

# If preprocessor switch -DENABLE_SERIAL_BUILD does not occur in compiler flags,
# a build for parallel execution is requested.
ifeq (,$(findstring -DENABLE_SERIAL_BUILD ,$(APPONLYFLAGS) $(CFLAGSF90) ))
ifeq ($(strip $(MPIWRAPPERS)), YES)
F77       = sxmpif90
F90       = sxmpif90
CC        = sxmpicc
CXX	  = echo
LD        = sxmpif90
endif
endif


##############################################################################
# Commands to get version information from compiler
##############################################################################

# Don't use 'sxmpicc' and 'sxmpic++' here to retrieve version information!
# These MPI wrapper commands call the linker even in case '-V' is passed
# which leads to a message on standard error. As such a message is annoying
# and deranges our compilation screen output, we used to redirect it to
# standard output and cut it of with 'head -n 2', but redirecting had the
# (unnoticed) side effect of creating a temporary directory named
# /tmp/ccomXXXXX.dir (XXXXX the current pid) - a directory every compilation
# process creates.
# Subsequent compilations could then fail with an error message of kind
#
#   sxcc: Could not create working directory: mkdir (errno=17) : /tmp/ccom11532.dir
#
# Using the cross compiler sxcc and sxc++ directly (not the MPI wrapper commands)
# and passing '-V' to them, the linker is not called. Redirecting to a) standard error
# and b) piping through head is still needed, though, as a) version information is
# reported to standard error, not standard output and b) the redirection somehow
# gets lost in the chain of commands used to store and compare compiler versions.
#
# sxmpif90 is not affected.
F77VERSION = $(F77) -V 2>&1 | head -n 3
F90VERSION = $(F90) -V 2>&1 | head -n 3
CCVERSION  = sxcc  -V 2>&1 | head -n 3
CXXVERSION = sxc++ -V 2>&1 | head -n 3


#		 -pi expin=$(OBJDIR) -pi infomsg \
#		 -pi expin=$(OBJDIR)/fsystem.f90,$(OBJDIR)/element.f90,$(OBJDIR)/grid.f90,$(OBJDIR)/statistics.f90,$(OBJDIR)/storage.f90,$(OBJDIR)/hlayer.f90,$(OBJDIR)/parallel.f90,$(OBJDIR)/auxiliary.f90,$(OBJDIR)/userdef.f90
##############################################################################
# compiler flags
# (including non-architecture specific optimisation flags)
##############################################################################

# Set default compile flags
ifeq ($(call optimise), YES)
CFLAGSF77     := -DUSE_COMPILER_NEC $(CFLAGSF77) -NE -C hopt -pi auto
#CFLAGSF77     = -C debug -e C
CFLAGSF90     := -DENABLE_CPPMACRO_FOR_STORAGEGETBASE -DDISABLE_ERRORCONTROL \
		 $(CFLAGSF90) $(filter-out -NE, $(CFLAGSF77)) \
		 -pi line=512 #-pi nest=4 -ftrace
CFLAGSC       := -DUSE_COMPILER_NEC $(CFLAGSC) -C hopt # -ftrace
LDFLAGS       := $(LDFLAGS) # -ftrace
else
CFLAGSF77     := -DUSE_COMPILER_NEC $(CFLAGSF77) -g -ftrace -e CPR
CFLAGSF90     := -DENABLE_CPPMACRO_FOR_STORAGEGETBASE -DDISABLE_ERRORCONTROL \
		 $(CFLAGSF90) $(CFLAGSF77)
CFLAGSC       := -DUSE_COMPILER_NEC $(CFLAGSC) -g -ftrace
LDFLAGS       := $(LDFLAGS) -ftrace
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
SBB_CVERSIONCMD = $(F77) -V 2>&1 | awk '{if ( NR == 1 ) print $1 $3}'


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
