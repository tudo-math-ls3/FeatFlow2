#!/usr/bin/env make
########################################################################
#                                                                      #
#         FINITE ELEMENT ANALYSIS & SOLUTION TOOLS  F E A S T          #
#                                                                      #
# Authors: Ch.Becker,                                                  #
#          S.Buijssen, D.Goeddeke, M.Grajewski, H.Wobker,              #
#          S.Kilian, S.Turek                                           #
#                                                                      #
# Contact: Applied Mathematics, Dortmund University of Technology      #
#          Vogelpothsweg 87, 44227 Dortmund                            #
#          Germany                                                     #
#                                                                      #
# Web:     http://www.feast.tu-dortmund.de                             #
#          mailto:feast@math.tu-dortmund.de                            #
#                                                                      #
########################################################################
#                                                                      #
# Compile settings, variables and targets to be included in every      #
# Makefile for FEAT2 applications                                      #
#                                                                      #
# Author    : Sven H.M. Buijssen                                       #
# Maintainer: Chr. Becker, Sven H.M. Buijssen                          #
# Version   : $Id: Makefile.cpu.inc,v 2.12 2009/01/28 09:49:25 buijssen Exp $
########################################################################


# The function match is used to match ID agains wildcards.
#
# variant 1 (slower, especially on old machines, but possibly more portable)
#match=$(shell echo $(1) | awk '/$(2)/ { print "yes"; }')
#
# variant 2 (fastest so far)
match=$(shell echo $(1) | sed -e 's/$(2)/yes/')

# Ensure that a request for full optimisation is written in capital letters
OPT := $(shell echo $(OPT) | tr '[a-z]' '[A-Z]')

# Every Makefile should contain this line to avoid trouble on systems
# where the SHELL variable might be inherited from the environment.
# (This is never a problem with GNU make.)
SHELL = /bin/sh

# Set location of buildconf.h (i.e. file containing build ID
# for automatic choice of best SBBLAS version within FEAT2)
BUILDCONF_H = "buildconf.h"

# Files to store compiler modification dates and compiler settings etc.
# for current build ID. This information is used to check
# whether the entire project needs to be recompiled.
FILE_F77_INFO     = env_f77_compiler
FILE_F90_INFO     = env_f90_compiler
FILE_CC_INFO      = env_c_compiler
FILE_CXX_INFO     = env_cxx_compiler
COMPILER_SETTINGS = env_compiler_settings

# Directory where Makefile template files are stored
TEMPLATESDIR = $(FEAT2BASEDIR)/templates

##############################################################################
# functions
#
##############################################################################

# Map different levels of optimisation to a single string
#
# variant 1 (fast, but not extendible)
#optimise = $(subst EXPENSIVE, YES, $(OPT))

# variant 2 (more overhead, but extendible).
ifeq ($(strip $(OPT)), EXPENSIVE)
optimise = YES
endif
ifeq ($(strip $(OPT)), YES)
optimise = YES
endif



##############################################################################
# clean up environment
#
##############################################################################

# Don't let the script get confused by non-english messages
# from system information programs.
# (LC_ALL overrides the value of the LANG environment variable
# and the values of any other LC_* environment variables.)
LC_ALL=C

# Unset CDPATH to prevent problems when changing directories
CDPATH =



##############################################################################
# compiler settings for various machine types, list to be extended
#
##############################################################################
# initial values for all architectures.
# Be careful what you specify here! All Compiler and linker settings are
# overwritten later on, but values for CFLAGS*, SRCEXTRA, BUILDLIB etc.
# are inherited!

# C preprocessor
CPP           = cpp

# FEAT2's preprocessor
F90CPP        = $(FEAT2BASEDIR)/bin/f90cpp

# Fortran 77 compiler
F77           = echo; echo "No Fortran 77 compiler specified."; echo; exit 1; echo
# Fortran 77 compiler options
CFLAGSF77LIBS =
CFLAGSF77     =

# Fortran 90 compiler
F90           = @echo; echo "No Fortran 90 compiler specified."; echo; exit 1; echo
# Fortran 90 compiler options
CFLAGSF90     =
# Unoptimised builds perform additional (possibly run-time expensive) parameter checks
ifeq ($(strip $(OPT)), NO)
CFLAGSF90    := $(CFLAGSF90) -DENABLE_PARAMETER_CHECK -DENABLE_ERROR_TRACEBACK
endif
# Append all active flags found in feastcppswitches.h
CFLAGSF90    := $(CFLAGSF90) $(foreach preprocessordirective, \
			$(shell sed -e '/^\#/d; /^[ ]*$$/d;' $(FEAT2BASEDIR)/feat2cppswitches.h), \
			-D$(preprocessordirective))

# extension of files the compiler creates containing module information
MODEXTENSION  =
# Some compilers expect module information in the source directory,
# others in the object directory. Set to "YES", if newly created
# module files should be moved to $(OBJDIR).
MOVEMOD       =

# C compiler
CC            = @echo; echo "No C compiler specified."; echo; exit 1; echo
# C compiler options
CFLAGSC       =

# C++ compiler
CXX           = @echo; echo "No C++ compiler specified."; echo; exit 1; echo
# C++ compiler options
CFLAGSCXX     =

# Linker
LD            = @echo; echo "No linker specified."; echo; exit 1; echo

# Linker flags
LDFLAGS       =

# If you have libraries installed in a directories that are not searched
# by default, add them to LD_LIBRARY_PATH and uncomment the following
# three lines.
#ifneq ($(strip $(LD_LIBRARY_PATH)),)
#LDFLAGS      := $(LDFLAGS) -L$(subst :, -L,$(LD_LIBRARY_PATH))
#endif

# include paths
INC           =

# List of files only needed for a specific build ID
# (please include path information, at best using Makefile variables)
SRCEXTRA      = $(FEAT2BASEDIR)/kernel/Postprocessing/gmvwritef.c \
		$(FEAT2BASEDIR)/kernel/System/signal_ccode.c

# List of all libraries available (queried for purge* targets)
ALLLIB        = amd blas lapack splib umfpack

# libraries to be built (e.g. 'blas' if there is no system-wide blas;
# on a typical system set this to 'metis umfpack' only)
# (Non-alphabetic order. Order should be like in linker step)
BUILDLIB      = umfpack amd splib

# directories to look in for library archives
LIBDIR        = -L$(OBJDIR)
ifneq ($(OBJDIR_LIB), $(OBJDIR))
LIBDIR        = -L$(OBJDIR) -L$(OBJDIR_LIB)
endif

# which libraries to link the application against
ifneq ($(strip $(BUILDLIB)),)
# Split up string if multiple libraries are given
LIBS         := $(patsubst %,-l%,$(BUILDLIB))
endif

# In case MPI wrapper commands for the compilers are unavailable,
# use the following variables to specify MPI include paths, directories
# where to look up MPI library archives and which MPI libraries to link
# the application against
MPIINC        =
MPILIBDIR     =
MPILIBS       =

# path and arguments to 'ar'
AR            =

# path and arguments to 'ranlib'
RANLIB        =

# libraries to be built for SparseBanded BLAS benchmark
# (e.g. 'blas' if there is no system-wide blas)
SBB_BUILDLIB  =
SBB_LIBS      =
SBB_SRCEXTRA  = ztime.f



##############################################################################
# FEAT2 supports both parallel and serial builds from the same
# application's source code. configure command line option '--mode=serial'
# results in the Makefile variable MODE being set to SERIAL. Translate this
# setting into one the C preprocessor understands.
#
##############################################################################

ifeq ($(MPI), NO)
CFLAGSF77 := -DENABLE_SERIAL_BUILD $(CFLAGSF77)
CFLAGSF90 := -DENABLE_SERIAL_BUILD $(CFLAGSF90)
CFLAGSC   := -DENABLE_SERIAL_BUILD $(CFLAGSC)
CFLAGSCXX := -DENABLE_SERIAL_BUILD $(CFLAGSCXX)
endif



##############################################################################

# This block extends the compiler options by machine dependent settings
# and libraries (optimised BLAS/LAPack versions for example).
# If external BLAS/LAPack is used in BLASLIB then BLAS/LAPack can be
# omitted in BUILDLIB to save some compilation time.

# Syntax rules:
# Whenever adding a rule, make sure it has the following syntax:
#   ifeq ($(call match,$(ID),foo-bar-baz-AA-BB-CC),yes)
# where <foo>, <bar>, <baz> may be arbitrary strings, but no regular expressions,
# The combination of <foo>-<bar>-<baz> should be a string returned by the script
# bin/guess_id. In order to be able to present the user with a list of possible
# build ID based on the result of a call of bin/guess_id, it is forbidden
# to specify a catch-all rule like ".*-.*-.*-.*-goto.*". The first three tokens
# must be specified for every rule.

# In order to check whether all six tokens that make up a valid build ID
# set a flag for every token when it get matched (and set).
#
TOKEN1 := 0
TOKEN2 := 0
TOKEN3 := 0
TOKEN4 := 0
TOKEN5 := 0
TOKEN6 := 0


################################### Alpha ####################################
#if (m/^ifeq \(\$\(call match,\$\(ID\),(.*)\),yes\)\s*$/) {
ifeq ($(call match,$(ID),alpha-ev6-osf1-cf90-.*),yes)
# invokes the Compaq Fortran 77/90 (formerly DIGITAL Fortran 77/90)
# compiler, local umfpack
# No MPI wrappers installed on our alphas
MPIWRAPPERS = NO
include $(TEMPLATESDIR)/alpha-ev6-osf1-generic.mk
include $(TEMPLATESDIR)/decf90-generic.mk
endif


ifeq ($(call match,$(ID),alpha-ev6-osf1-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),alpha-ev6-osf1-.*-atlas.*),yes)
# Use ATLAS BLAS
include $(TEMPLATESDIR)/atlasblas-generic.mk
endif


ifeq ($(call match,$(ID),alpha-ev6-osf1-.*-dxml.*),yes)
# Use Digital BLAS
TOKEN5 := 1
LIBS      := $(LIBS) -ldxml
endif


ifeq ($(call match,$(ID),alpha-ev6-osf1-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),alpha-ev6-osf1-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


ifeq ($(call match,$(ID),alpha-ev6-osf1-.*-.*-dmpi),yes)
# Use Digital MPI
include $(TEMPLATESDIR)/dmpi-generic.mk
endif


ifeq ($(call match,$(ID),alpha-ev6-osf1-.*-.*-lammpi),yes)
# Use LAM/MPI
include $(TEMPLATESDIR)/lammpi-generic.mk
# On Alpha we had name conflicts with libmpi.a (Digital vs. LAM/MPI)
# and hence renamed the LAM/MPI libmpi.a to libmpilam.a. So, we need to
# override the generic LAM/MPI setting.
ifeq ($(strip $(MPIWRAPPERS)), NO)
MPILIBS = -llamf77mpi -lmpilam -llam -lpthread
endif
endif



################################ IBM Regatta #################################

ifeq ($(call match,$(ID),ibm-powerpc_power4-aix-xlf-.*-poempi),yes)
# invokes the IBM XL Fortran compiler, uses POE MPI
# local umfpack
TOKEN6 := 1
include $(TEMPLATESDIR)/xlf-generic.mk
include $(TEMPLATESDIR)/ibm-powerpc_power4-aix-generic.mk
# The system provides a METIS library, but we have problems using it.
# Ensure that we use our own METIS library by hard-coding the path.
LIBS      = $(OBJDIR_LIB)/libmetis.a -lumfpack -lamd
#SBB_LIBS := $(SBB_LIBS)
endif

ifeq ($(call match,$(ID),ibm-powerpc_power4-aix-.*-essl.*),yes)
# Use ESSL BLAS and self-compiled LAPACK
TOKEN5 := 1
BUILDLIB := $(BUILDLIB) lapack
LIBS     := $(LIBS) -llapack -lessl 
SBB_BUILDLIB := $(SBB_BUILDLIB) lapack
SBB_LIBS := $(SBB_LIBS) -llapack-lessl
endif

ifeq ($(call match,$(ID),ibm-powerpc_power4-aix-.*-goto.*),yes)
# Use Goto BLAS and self-compiled LAPACK
TOKEN5 := 1
BUILDLIB := $(BUILDLIB) lapack
LIBS     := $(LIBS) -lblas_Goto -llapack
SBB_BUILDLIB := $(SBB_BUILDLIB) lapack
SBB_LIBS := $(SBB_LIBS) -lblas_Goto -llapack
endif



################################ Itanium 2 #################################

ifeq ($(call match,$(ID),ia64-itanium2-linux-intel-blas-mpi),yes)
# invokes the intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/intel64-generic.mk
include $(TEMPLATESDIR)/blas-generic.mk
# Add preprocessor flag to enable 64bit UMFPACK support.
# UMFPACK documentation states to also add -DLP64, but there is no
# reference to it in the UMFPACK 5.2 sources.
CFLAGSF90     := -DENABLE_64BIT_UMFPACK $(CFLAGSF90)
CFLAGSC       := -DENABLE_64BIT_UMFPACK $(CFLAGSC)
TOKEN6 := 1
endif

ifeq ($(call match,$(ID),ia64-itanium2-linux-intel-mkl-mpi),yes)
# invokes the intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/intel64-generic.mk
# Add preprocessor flag to enable 64bit UMFPACK support.
# UMFPACK documentation states to also add -DLP64, but there is no
# reference to it in the UMFPACK 5.2 sources.
CFLAGSF90     := -DENABLE_64BIT_UMFPACK $(CFLAGSF90)
CFLAGSC       := -DENABLE_64BIT_UMFPACK $(CFLAGSC)
LIBS     := $(LIBS) -L/opt/MathKeisan/MKL/lib/64 -lguide -lmkl_lapack64 -lmkl
#LIBS     := $(LIBS) -L/nfs/home11/HLRS/xam/xamsbuij/intel/Compiler/11.0/074/mkl/lib/64 -lguide -lmkl_lapack -lmkl
SBB_LIBS := $(SBB_LIBS) -L/opt/MathKeisan/MKL/lib/64 -lguide -lmkl_lapack64 -lmkl
TOKEN5 := 1
TOKEN6 := 1
endif



############################ Itanium 2 Dual Core #############################

ifeq ($(call match,$(ID),ia64-itanium2x2-cygwin_nt.*-g95-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
#CFLAGSF77:=$(CFLAGSF77) -fno-underscoring
#CFLAGSF90:=$(CFLAGSF90) -fno-underscoring
endif

ifeq ($(call match,$(ID),ia64-itanium2x2-linux-g95-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
#CFLAGSF77:=$(CFLAGSF77) -fno-underscoring
#CFLAGSF90:=$(CFLAGSF90) -fno-underscoring
endif


ifeq ($(call match,$(ID),ia64-itanium2x2-linux-gcc-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -mtune=itanium2 $(CFLAGSF77LIBS)
CFLAGSF77     := -mtune=itanium2 $(CFLAGSF77)
CFLAGSF90     := -mtune=itanium2 $(CFLAGSF90)
CFLAGSC       := -mtune=itanium2 $(CFLAGSC)
CFLAGSCXX     := -mtune=itanium2 $(CFLAGSCXX)
LDFLAGS       := -mtune=itanium2 $(LDFLAGS)
endif
#CFLAGSF77:=$(CFLAGSF77) -fno-underscoring
#CFLAGSF90:=$(CFLAGSF90) -fno-underscoring
endif


ifeq ($(call match,$(ID),ia64-itanium2x2-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
endif


ifeq ($(call match,$(ID),ia64-itanium2x2-linux-.*-mkl.*),yes)
ifneq ($(strip $(INTEL_MKL_LIB)),)
# Split up string if multiple directories are given
# Note: Do not put whitespace between comma and the environment variable, because
#       if you do, a string like "-L /path/to/mkl" is the result and that string
#       won't make it into the command line.
LIBDIR   := $(LIBDIR) -L$(subst :, -L,$(INTEL_MKL_LIB))
endif
LIBS     := $(LIBS) -lmkl_lapack -lmkl_intel_thread -lmkl_core -lmkl_intel_lp64 -lguide -lpthread
SBB_LIBS := $(SBB_LIBS) -lmkl_lapack -lmkl_intel_thread -lmkl_core -lmkl_intel_lp64 -lguide -lpthread
TOKEN5 := 1
endif


ifeq ($(call match,$(ID),ia64-itanium2x2-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),ia64-itanium2x2-linux-.*-.*-ompi),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/openmpi-generic.mk
# Add preprocessor flag to enable 64bit UMFPACK support.
# UMFPACK documentation states to also add -DLP64, but there is no
# reference to it in the UMFPACK 5.2 sources.
CFLAGSF90     := -DENABLE_64BIT_UMFPACK $(CFLAGSF90)
CFLAGSC       := -DENABLE_64BIT_UMFPACK $(CFLAGSC)
endif



################################# Athlon 64 ##################################

ifeq ($(call match,$(ID),pc-athlon64-cygwin_nt.*-g95-.*),yes)
# invokes the G95 Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=athlon64 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=athlon64 $(CFLAGSF77)
CFLAGSF90     := -march=athlon64 $(CFLAGSF90)
CFLAGSC       := -march=athlon64 $(CFLAGSC)
CFLAGSCXX     := -march=athlon64 $(CFLAGSCXX)
LDFLAGS       := -march=athlon64 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-athlon64-linux-g95-.*),yes)
# invokes the G95 Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=athlon64 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=athlon64 $(CFLAGSF77)
CFLAGSF90     := -march=athlon64 $(CFLAGSF90)
CFLAGSC       := -march=athlon64 $(CFLAGSC)
CFLAGSCXX     := -march=athlon64 $(CFLAGSCXX)
LDFLAGS       := -march=athlon64 $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-athlon64-linux-gcc-.*),yes)
# invokes the gfortran Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=athlon64 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=athlon64 $(CFLAGSF77)
CFLAGSF90     := -march=athlon64 $(CFLAGSF90)
CFLAGSC       := -march=athlon64 $(CFLAGSC)
CFLAGSCXX     := -march=athlon64 $(CFLAGSCXX)
LDFLAGS       := -march=athlon64 $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-athlon64-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -xW $(CFLAGSF77LIBS)
CFLAGSF77     := -xW $(CFLAGSF77)
CFLAGSF90     := -xW $(CFLAGSF90)
CFLAGSC       := -xW $(CFLAGSC)
CFLAGSCXX     := -xW $(CFLAGSCXX)
LDFLAGS       := -xW $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-athlon64-linux-sunstudio-.*),yes)
# invokes the Sunstudio Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/sunstudio-generic.mk
# # # # # # # # # # #      WARNING !!!     # # # # # # # # # # #
# Architecture-specific optimisation missing!
# DOM: sunstudio seems to do this automatically.
# # # # # # # # # # #      WARNING !!!     # # # # # # # # # # #
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -xautopar  $(CFLAGSF77LIBS)
CFLAGSF77     := -xautopar  $(CFLAGSF77)
CFLAGSF90     := -xautopar  $(CFLAGSF90)
CFLAGSC       := -xautopar  $(CFLAGSC)
CFLAGSCXX     := -xautopar  $(CFLAGSCXX)
LDFLAGS       := -xautopar  $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-athlon64-linux-.*-.*-lammpi),yes)
# Use LAM/MPI
include $(TEMPLATESDIR)/lammpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-athlon64-linux-.*-.*-mpich2),yes)
# Use MPICH2
include $(TEMPLATESDIR)/mpich2-generic.mk
endif


ifeq ($(call match,$(ID),pc-athlon64-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-athlon64-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc-athlon64-linux-.*-goto),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc-athlon64-linux-.*-goto2),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


################################# Phenom II x4 ##################################
ifeq ($(call match,$(ID),pc64-phenomIIx4-linux-gcc-.*),yes)
# invokes the gfortran Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=amdfam10 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=amdfam10 $(CFLAGSF77)
CFLAGSF90     := -march=amdfam10 $(CFLAGSF90)
CFLAGSC       := -march=amdfam10 $(CFLAGSC)
CFLAGSCXX     := -march=amdfam10 $(CFLAGSCXX)
LDFLAGS       := -march=amdfam10 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-phenomIIx4-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif

################################# Turion X2 ##################################

ifeq ($(call match,$(ID),pc-turionx2-cygwin_nt.*-g95-*),yes)
# invokes the G95 Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
include $(TEMPLATESDIR)/openmpi-generic.mk
include $(TEMPLATESDIR)/gotoblas-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=athlon64 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=athlon64 $(CFLAGSF77)
CFLAGSF90     := -march=athlon64 $(CFLAGSF90)
CFLAGSC       := -march=athlon64 $(CFLAGSC)
CFLAGSCXX     := -march=athlon64 $(CFLAGSCXX)
LDFLAGS       := -march=athlon64 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-turionx2-linux-g95-goto-ompi),yes)
# invokes the G95 Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
include $(TEMPLATESDIR)/openmpi-generic.mk
include $(TEMPLATESDIR)/gotoblas-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=athlon64 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=athlon64 $(CFLAGSF77)
CFLAGSF90     := -march=athlon64 $(CFLAGSF90)
CFLAGSC       := -march=athlon64 $(CFLAGSC)
CFLAGSCXX     := -march=athlon64 $(CFLAGSCXX)
LDFLAGS       := -march=athlon64 $(LDFLAGS)
endif
endif



################################# Athlon XP ##################################

ifeq ($(call match,$(ID),pc-athlonxp-cygwin_nt.*-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
# -mmx -m3dnow -mfpmath=sse -msse are automatically set by -march=athlon-xp
CFLAGSF77LIBS := -march=athlon-xp $(CFLAGSF77LIBS)
CFLAGSF77     := -march=athlon-xp $(CFLAGSF77)
CFLAGSF90     := -march=athlon-xp $(CFLAGSF90)
CFLAGSC       := -march=athlon-xp $(CFLAGSC)
endif
endif

ifeq ($(call match,$(ID),pc-athlonxp-linux-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
# -mmx -m3dnow -mfpmath=sse -msse are automatically set by -march=athlon-xp
CFLAGSF77LIBS := -march=athlon-xp $(CFLAGSF77LIBS)
CFLAGSF77     := -march=athlon-xp $(CFLAGSF77)
CFLAGSF90     := -march=athlon-xp $(CFLAGSF90)
CFLAGSC       := -march=athlon-xp $(CFLAGSC)
endif
endif


ifeq ($(call match,$(ID),pc-athlonxp-linux-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
# -mmx -m3dnow -mfpmath=sse -msse are automatically set by -march=athlon-xp
CFLAGSF77LIBS := -march=athlon-xp $(CFLAGSF77LIBS)
CFLAGSF77     := -march=athlon-xp $(CFLAGSF77)
CFLAGSF90     := -march=athlon-xp $(CFLAGSF90)
CFLAGSC       := -march=athlon-xp $(CFLAGSC)
endif
endif


ifeq ($(call match,$(ID),pc-athlonxp-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -axK $(CFLAGSF77LIBS) #-tpp6
CFLAGSF77     := -axK $(CFLAGSF77)     #-tpp6
CFLAGSF90     := -axK $(CFLAGSF90)     #-tpp6
CFLAGSC       := -axK $(CFLAGSC)       #-tpp6
LDFLAGS       := -axK $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-athlonxp-linux-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp=athlonxp $(CFLAGSF77LIBS)
CFLAGSF77     := -tp=athlonxp $(CFLAGSF77)
CFLAGSF90     := -tp=athlonxp $(CFLAGSF90)
CFLAGSC       := -tp=athlonxp $(CFLAGSC)
endif
endif


ifeq ($(call match,$(ID),pc-athlonxp-linux-.*-.*-lammpi),yes)
# Use LAM/MPI
include $(TEMPLATESDIR)/lammpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-athlonxp-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-athlonxp-linux-.*-atlas.*),yes)
# Use ATLAS BLAS
include $(TEMPLATESDIR)/atlasblas-generic.mk
endif


ifeq ($(call match,$(ID),pc-athlonxp-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc-athlonxp-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc-athlonxp-linux-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


################################# CoreDuo 32bit ################################

ifeq ($(call match,$(ID),pc-coreduo-cygwin_nt.*-g95-.*),yes)
# invokes the GNU Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -msse3 -pipe -fomit-frame-pointer $(CFLAGSF77LIBS)
CFLAGSF77     := -msse3 -pipe -fomit-frame-pointer $(CFLAGSF77)
CFLAGSF90     := -msse3 -pipe -fomit-frame-pointer $(CFLAGSF90)
CFLAGSC       := -msse3 -pipe -fomit-frame-pointer $(CFLAGSC)
CFLAGSCXX     := -msse3 -pipe -fomit-frame-pointer $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-coreduo-cygwin_nt.*-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc-coreduo-linux-g95-.*),yes)
# invokes the GNU Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -msse3 -pipe -fomit-frame-pointer $(CFLAGSF77LIBS)
CFLAGSF77     := -msse3 -pipe -fomit-frame-pointer $(CFLAGSF77)
CFLAGSF90     := -msse3 -pipe -fomit-frame-pointer $(CFLAGSF90)
CFLAGSC       := -msse3 -pipe -fomit-frame-pointer $(CFLAGSC)
CFLAGSCXX     := -msse3 -pipe -fomit-frame-pointer $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-coreduo-linux-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
# optimal flags taken from script CPUFlags v0.7.2
#
# According to http://archives.gentoo.org/gentoo-amd64/msg_14402.xml
# H.J. Lu from Intel, very prominent on gcc 4.x bugzilla, advices
# different settings:
# * GCC 4.1:
#   - use -march=prescott for Intel Core Solo/Duo
#   - use -march=nocona (and an amd64 profile) for Core 2 Solo/Duo
# * GCC 4.2:
#   - use -march=prescott -mtune=generic for Intel Core Solo/Duo
#   - use -march=nocona -mtune=generic for Intel Core 2 Solo/Duo
CFLAGSF77LIBS := -march=nocona -msse3 -mfpmath=sse $(CFLAGSF77LIBS)
CFLAGSF77     := -march=nocona -msse3 -mfpmath=sse $(CFLAGSF77)
CFLAGSF90     := -march=nocona -msse3 -mfpmath=sse $(CFLAGSF90)
CFLAGSC       := -march=nocona -msse3 -mfpmath=sse $(CFLAGSC)
endif
endif


ifeq ($(call match,$(ID),pc-coreduo-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -xT $(CFLAGSF77LIBS)
CFLAGSF77     := -xT $(CFLAGSF77)
CFLAGSF90     := -xT $(CFLAGSF90)
CFLAGSC       := -xT $(CFLAGSC)
CFLAGSCXX     := -xT $(CFLAGSCXX)
LDFLAGS       := -xT $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-coreduo-linux-.*-.*-lammpi),yes)
# Use LAM/MPI
include $(TEMPLATESDIR)/lammpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-coreduo-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-coreduo-linux-.*-atlas.*),yes)
# Use ATLAS BLAS
include $(TEMPLATESDIR)/atlasblas-generic.mk
endif


ifeq ($(call match,$(ID),pc-coreduo-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc-coreduo-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc-coreduo-linux-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


ifeq ($(call match,$(ID),pc-coreduo-linux-.*-mkl.*),yes)
# Use MKL BLAS
include $(TEMPLATESDIR)/mkl-generic.mk
endif


################################# CoreDuo 64bit ################################

ifeq ($(call match,$(ID),pc64-coreduo-linux-g95-.*),yes)
# invokes the GNU Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -msse3 -pipe -fomit-frame-pointer $(CFLAGSF77LIBS)
CFLAGSF77     := -msse3 -pipe -fomit-frame-pointer $(CFLAGSF77)
CFLAGSF90     := -msse3 -pipe -fomit-frame-pointer $(CFLAGSF90)
CFLAGSC       := -msse3 -pipe -fomit-frame-pointer $(CFLAGSC)
CFLAGSCXX     := -msse3 -pipe -fomit-frame-pointer $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc64-coreduo-linux-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
# optimal flags taken from script CPUFlags v0.7.2
#
# According to http://archives.gentoo.org/gentoo-amd64/msg_14402.xml
# H.J. Lu from Intel, very prominent on gcc 4.x bugzilla, advices
# different settings:
# * GCC 4.1:
#   - use -march=prescott for Intel Core Solo/Duo
#   - use -march=nocona (and an amd64 profile) for Core 2 Solo/Duo
# * GCC 4.2:
#   - use -march=prescott -mtune=generic for Intel Core Solo/Duo
#   - use -march=nocona -mtune=generic for Intel Core 2 Solo/Duo
CFLAGSF77LIBS := -march=nocona -msse3 -mfpmath=sse $(CFLAGSF77LIBS)
CFLAGSF77     := -march=nocona -msse3 -mfpmath=sse $(CFLAGSF77)
CFLAGSF90     := -march=nocona -msse3 -mfpmath=sse $(CFLAGSF90)
CFLAGSC       := -march=nocona -msse3 -mfpmath=sse $(CFLAGSC)
endif
endif


ifeq ($(call match,$(ID),pc64-coreduo-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -xT $(CFLAGSF77LIBS)
CFLAGSF77     := -xT $(CFLAGSF77)
CFLAGSF90     := -xT $(CFLAGSF90)
CFLAGSC       := -xT $(CFLAGSC)
CFLAGSCXX     := -xT $(CFLAGSCXX)
LDFLAGS       := -xT $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc64-coreduo-linux-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: Not yet sure whether optimisation flags are optimal!';
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -tp=core2-64
CFLAGSF77     := $(CFLAGSF77) -tp=core2-64
CFLAGSF90     := $(CFLAGSF90) -tp=core2-64
CFLAGSC       := $(CFLAGSC) -tp=core2-64
LDFLAGS       := $(LDFLAGS) -tp=core2-64
endif
# Probably activate 64 bit integer support
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -i8
CFLAGSF90     := $(CFLAGSF90) -i8
endif
endif


ifeq ($(call match,$(ID),pc64-coreduo-linux-.*-.*-lammpi),yes)
# Use LAM/MPI
include $(TEMPLATESDIR)/lammpi-generic.mk
endif


ifeq ($(call match,$(ID),pc64-coreduo-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc64-coreduo-linux-.*-.*-mpich2),yes)
# Use MPICH2
include $(TEMPLATESDIR)/mpich2-generic.mk
endif


ifeq ($(call match,$(ID),pc64-coreduo-linux-.*-acml.*),yes)
# Use AMD Core Math Library
include $(TEMPLATESDIR)/acml-generic.mk
endif


ifeq ($(call match,$(ID),pc64-coreduo-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc64-coreduo-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc64-coreduo-linux-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


ifeq ($(call match,$(ID),pc64-coreduo-linux-.*-mkl.*),yes)
# Use MKL BLAS
include $(TEMPLATESDIR)/mkl-generic.mk
endif


ifeq ($(call match,$(ID),pc64-coreduo-linux-.*-woodymkl.*),yes)
# Use special BKL blas on Woodcrest cluster of the RRZE Erlangen
LIBS := $(LIBS) $(MKL_LIB)
TOKEN5 := 1
endif


################################# CoreSolo 32 bit ###############################

ifeq ($(call match,$(ID),pc-coresolo-cygwin_nt.*-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
# According to http://archives.gentoo.org/gentoo-amd64/msg_14402.xml
# H.J. Lu from Intel, very prominent on gcc 4.x bugzilla, advices
# different settings:
# * GCC 4.1:
#   - use -march=prescott for Intel Core Solo/Duo
#   - use -march=nocona (and an amd64 profile) for Core 2 Solo/Duo
# * GCC 4.2:
#   - use -march=prescott -mtune=generic for Intel Core Solo/Duo
#   - use -march=nocona -mtune=generic for Intel Core 2 Solo/Duo
CFLAGSF77LIBS := -march=pentium-m -mtune=prescott -msse $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium-m -mtune=prescott -msse $(CFLAGSF77)
CFLAGSF90     := -march=pentium-m -mtune=prescott -msse $(CFLAGSF90)
CFLAGSC       := -march=pentium-m -mtune=prescott -msse $(CFLAGSC)
CFLAGSCXX     := -march=pentium-m -mtune=prescott -msse $(CFLAGSCXX)
endif
endif


ifeq ($(call match,$(ID),pc-coresolo-linux-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
# According to http://archives.gentoo.org/gentoo-amd64/msg_14402.xml
# H.J. Lu from Intel, very prominent on gcc 4.x bugzilla, advices
# different settings:
# * GCC 4.1:
#   - use -march=prescott for Intel Core Solo/Duo
#   - use -march=nocona (and an amd64 profile) for Core 2 Solo/Duo
# * GCC 4.2:
#   - use -march=prescott -mtune=generic for Intel Core Solo/Duo
#   - use -march=nocona -mtune=generic for Intel Core 2 Solo/Duo
CFLAGSF77LIBS := -march=pentium-m -mtune=prescott -msse $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium-m -mtune=prescott -msse $(CFLAGSF77)
CFLAGSF90     := -march=pentium-m -mtune=prescott -msse $(CFLAGSF90)
CFLAGSC       := -march=pentium-m -mtune=prescott -msse $(CFLAGSC)
CFLAGSCXX     := -march=pentium-m -mtune=prescott -msse $(CFLAGSCXX)
endif
endif


ifeq ($(call match,$(ID),pc-coresolo-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -xP $(CFLAGSF77LIBS)
CFLAGSF77     := -xP $(CFLAGSF77)
CFLAGSF90     := -xP $(CFLAGSF90)
CFLAGSC       := -xP $(CFLAGSC)
CFLAGSCXX     := -xP $(CFLAGSCXX)
LDFLAGS       := -xP $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-coresolo-linux-gcc-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=prescott $(CFLAGSF77LIBS)
CFLAGSF77     := -march=prescott $(CFLAGSF77)
CFLAGSF90     := -march=prescott $(CFLAGSF90)
CFLAGSC       := -march=prescott $(CFLAGSC)
CFLAGSCXX     := -march=prescott $(CFLAGSCXX)
LDFLAGS       := -march=prescott $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-coresolo-linux-sunstudio-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/sunstudio-generic.mk
endif


ifeq ($(call match,$(ID),pc-coresolo-linux-.*-.*-lammpi),yes)
# Use LAM/MPI
include $(TEMPLATESDIR)/lammpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-coresolo-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-coresolo-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc-coresolo-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc-coresolo-linux-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


############################### Opteron 32bit ################################

ifeq ($(call match,$(ID),pc-opteron-cygwin_nt.*-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=opteron $(CFLAGSF77LIBS)
CFLAGSF77     := -march=opteron $(CFLAGSF77)
CFLAGSF90     := -march=opteron $(CFLAGSF90)
CFLAGSC       := -march=opteron $(CFLAGSC)
CFLAGSCXX     := -march=opteron $(CFLAGSCXX)
endif
endif


ifeq ($(call match,$(ID),pc-opteron-linux-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=opteron $(CFLAGSF77LIBS)
CFLAGSF77     := -march=opteron $(CFLAGSF77)
CFLAGSF90     := -march=opteron $(CFLAGSF90)
CFLAGSC       := -march=opteron $(CFLAGSC)
CFLAGSCXX     := -march=opteron $(CFLAGSCXX)
endif
endif


ifeq ($(call match,$(ID),pc-opteron-linux-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=opteron $(CFLAGSF77LIBS)
CFLAGSF77     := -march=opteron $(CFLAGSF77)
CFLAGSF90     := -march=opteron $(CFLAGSF90)
CFLAGSC       := -march=opteron $(CFLAGSC)
endif
endif


ifeq ($(call match,$(ID),pc-opteron-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -xW $(CFLAGSF77LIBS)
CFLAGSF77     := -xW $(CFLAGSF77)
CFLAGSF90     := -xW $(CFLAGSF90)
CFLAGSC       := -xW $(CFLAGSC)
CFLAGSCXX     := -xW $(CFLAGSCXX)
LDFLAGS       := -xW $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-opteron-linux-sunstudio-.*),yes)
# invokes the Sun Studio Compiler suite, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/sunstudio-generic.mk
# # # # # # # # # # #      WARNING !!!     # # # # # # # # # # #
# Architecture-specific optimisation missing!
# DOM: sunstudio seems to do this automatically.
# # # # # # # # # # #      WARNING !!!     # # # # # # # # # # #
#ifeq ($(call optimise), YES)
#CFLAGSF77LIBS := -xtarget=opteron -xarch=<DON'T KNOW> $(CFLAGSF77LIBS)
#CFLAGSF77     := -xtarget=opteron -xarch=<DON'T KNOW> $(CFLAGSF77)
#CFLAGSF90     := -xtarget=opteron -xarch=<DON'T KNOW> $(CFLAGSF90)
#CFLAGSC       := -xtarget=opteron -xarch=<DON'T KNOW> -xc99 $(CFLAGSC)
#endif
#CFLAGSC       := -xc99 $(CFLAGSC)
endif


ifeq ($(call match,$(ID),pc-opteron-linux-.*-.*-lammpi),yes)
# Use LAM/MPI
include $(TEMPLATESDIR)/lammpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-opteron-linux-.*-.*-mpich2),yes)
# Use MPICH2
include $(TEMPLATESDIR)/mpich2-generic.mk
endif


ifeq ($(call match,$(ID),pc-opteron-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-opteron-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc-opteron-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc-opteron-linux-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


########################### Opteron 32bit DualCore ##########################

ifeq ($(call match,$(ID),pc-opteronx2-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -xW #-tpp7
CFLAGSF77     := $(CFLAGSF77) -xW #-tpp7
CFLAGSF90     := $(CFLAGSF90) -xW #-tpp7
CFLAGSC       := $(CFLAGSC) -xW #-tpp7
LDFLAGS       := $(LDFLAGS) -xW #-tpp7
endif
endif


ifeq ($(call match,$(ID),pc-opteronx2-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-opteronx2-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc-opteronx2-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc-opteronx2-linux-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


########################### Opteron 64bit DualCore ##########################

ifeq ($(call match,$(ID),pc64-opteronx2-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -xW #-tpp7
CFLAGSF77     := $(CFLAGSF77) -xW #-tpp7
CFLAGSF90     := $(CFLAGSF90) -xW #-tpp7
CFLAGSC       := $(CFLAGSC) -xW #-tpp7
LDFLAGS       := $(LDFLAGS) -xW #-tpp7
endif
endif


ifeq ($(call match,$(ID),pc64-opteronx2-cygwin_nt.*-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -m64
CFLAGSF77     := $(CFLAGSF77) -m64
CFLAGSF90     := $(CFLAGSF90) -m64
CFLAGSC       := $(CFLAGSC) -m64
LDFLAGS       := $(LDFLAGS) -m64
ifeq ($(call optimise), YES)
# -mmx -m3dnow -mfpmath=sse -msse -msse2 are automatically set by -march=opteron
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -march=opteron
CFLAGSF77     := $(CFLAGSF77) -march=opteron
CFLAGSF90     := $(CFLAGSF90) -march=opteron
CFLAGSC       := $(CFLAGSC) -march=opteron
endif
endif


ifeq ($(call match,$(ID),pc64-opteronx2-linux-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -m64
CFLAGSF77     := $(CFLAGSF77) -m64
CFLAGSF90     := $(CFLAGSF90) -m64
CFLAGSC       := $(CFLAGSC) -m64
LDFLAGS       := $(LDFLAGS) -m64
ifeq ($(call optimise), YES)
# -mmx -m3dnow -mfpmath=sse -msse -msse2 are automatically set by -march=opteron
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -march=opteron
CFLAGSF77     := $(CFLAGSF77) -march=opteron
CFLAGSF90     := $(CFLAGSF90) -march=opteron
CFLAGSC       := $(CFLAGSC) -march=opteron
endif
endif


ifeq ($(call match,$(ID),pc64-opteronx2-linux-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -m64
CFLAGSF77     := $(CFLAGSF77) -m64
CFLAGSF90     := $(CFLAGSF90) -m64
CFLAGSC       := $(CFLAGSC) -m64
LDFLAGS       := $(LDFLAGS) -m64
ifeq ($(call optimise), YES)
# -mmx -m3dnow -mfpmath=sse -msse -msse2 are automatically set by -march=opteron
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -march=opteron
CFLAGSF77     := $(CFLAGSF77) -march=opteron
CFLAGSF90     := $(CFLAGSF90) -march=opteron
CFLAGSC       := $(CFLAGSC) -march=opteron
endif
endif


ifeq ($(call match,$(ID),pc64-opteronx2-linux-sunstudio-.*),yes)
# invokes the Sun Studio Compiler suite, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/sunstudio-generic.mk
# Sun Studio seems to optimise automatically for the detected CPU.
endif


ifeq ($(call match,$(ID),pc64-opteronx2-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc64-opteronx2-linux-.*-.*-lanlmpi),yes)
# Use LANL MPI
include $(TEMPLATESDIR)/lanlmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc64-opteronx2-linux-.*-.*-mpich2),yes)
# Use MPICH2
include $(TEMPLATESDIR)/mpich2-generic.mk
endif


ifeq ($(call match,$(ID),pc64-opteronx2-linux-.*-mkl.*),yes)
# Use Intel Math Kernel Library
include $(TEMPLATESDIR)/mkl-generic.mk
endif


ifeq ($(call match,$(ID),pc64-opteronx2-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc64-opteronx2-linux-.*-goto),yes)
# Use GotoBLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc64-opteronx2-linux-.*-goto2),yes)
# Use GotoBLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
# Workaround for linker problem:
# ld: hidden symbol `__svml_cosf4' in libsvml.a(svml_stub_scos4.o) is referenced by DSO
ifeq ($(call match,$(ID),pc64-opteronx2-linux-intel-goto2.*),yes)
LIBS := $(LIBS) -lsvml
endif
endif


ifeq ($(call match,$(ID),pc64-opteronx2-linux-.*-.*),yes)
# Add preprocessor flag to enable 64bit UMFPACK support.
# UMFPACK documentation states to also add -DLP64, but there is no
# reference to it in the UMFPACK 5.2 sources.
#CFLAGSF90     := -DENABLE_64BIT_UMFPACK $(CFLAGSF90)
#CFLAGSC       := -DENABLE_64BIT_UMFPACK $(CFLAGSC)
endif


############################### Opteron 64bit ################################

ifeq ($(call match,$(ID),pc64-opteron-cygwin_nt.*-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -m64
CFLAGSF77     := $(CFLAGSF77) -m64
CFLAGSF90     := $(CFLAGSF90) -m64
CFLAGSC       := $(CFLAGSC) -m64
LDFLAGS       := $(LDFLAGS) -m64
ifeq ($(call optimise), YES)
# -mmx -m3dnow -mfpmath=sse -msse -msse2 are automatically set by -march=opteron
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -march=opteron
CFLAGSF77     := $(CFLAGSF77) -march=opteron
CFLAGSF90     := $(CFLAGSF90) -march=opteron
CFLAGSC       := $(CFLAGSC) -march=opteron
endif
endif


ifeq ($(call match,$(ID),pc64-opteron-linux-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -m64
CFLAGSF77     := $(CFLAGSF77) -m64
CFLAGSF90     := $(CFLAGSF90) -m64
CFLAGSC       := $(CFLAGSC) -m64
LDFLAGS       := $(LDFLAGS) -m64
ifeq ($(call optimise), YES)
# -mmx -m3dnow -mfpmath=sse -msse -msse2 are automatically set by -march=opteron
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -march=opteron
CFLAGSF77     := $(CFLAGSF77) -march=opteron
CFLAGSF90     := $(CFLAGSF90) -march=opteron
CFLAGSC       := $(CFLAGSC) -march=opteron
endif
# Probably activate 64 bit integer support
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -i8
CFLAGSF90     := $(CFLAGSF90) -i8
endif
endif


ifeq ($(call match,$(ID),pc64-opteron-linux-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -m64
CFLAGSF77     := $(CFLAGSF77) -m64
CFLAGSF90     := $(CFLAGSF90) -m64
CFLAGSC       := $(CFLAGSC) -m64
LDFLAGS       := $(LDFLAGS) -m64
ifeq ($(call optimise), YES)
# -mmx -m3dnow -mfpmath=sse -msse -msse2 are automatically set by -march=opteron
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -march=opteron
CFLAGSF77     := $(CFLAGSF77) -march=opteron
CFLAGSF90     := $(CFLAGSF90) -march=opteron
CFLAGSC       := $(CFLAGSC) -march=opteron
endif
# Probably activate 64 bit integer support
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -fdefault-integer-8
CFLAGSF90     := $(CFLAGSF90) -fdefault-integer-8
endif
endif


ifeq ($(call match,$(ID),pc64-opteron-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -xW #-tpp7
CFLAGSF77     := $(CFLAGSF77) -xW     #-tpp7
CFLAGSF90     := $(CFLAGSF90) -xW     #-tpp7
CFLAGSC       := $(CFLAGSC) -xW       #-tpp7
LDFLAGS       := $(LDFLAGS) -xW       #-tpp7
endif
# Probably activate 64 bit integer support
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -integer_size 64 
CFLAGSF90     := $(CFLAGSF90) -integer_size 64 
endif
endif


ifeq ($(call match,$(ID),pc64-opteron-linux-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -tp=k8-64
CFLAGSF77     := $(CFLAGSF77) -tp=k8-64
CFLAGSF90     := $(CFLAGSF90) -tp=k8-64
CFLAGSC       := $(CFLAGSC) -tp=k8-64
LDFLAGS       := $(LDFLAGS) -tp=k8-64
endif
# Probably activate 64 bit integer support
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -i8
CFLAGSF90     := $(CFLAGSF90) -i8
endif
endif


ifeq ($(call match,$(ID),pc64-opteron-linux-psc-.*),yes)
# invokes the PathScale Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pathscale-generic.mk
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -m64
CFLAGSF77     := $(CFLAGSF77) -m64
CFLAGSF90     := $(CFLAGSF90) -m64
CFLAGSC       := $(CFLAGSC) -m64
LDFLAGS       := $(LDFLAGS) -m64
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -mtune=opteron -mcpu=opteron -march=opteron
CFLAGSF77     := $(CFLAGSF77) -mtune=opteron -mcpu=opteron -march=opteron
CFLAGSF90     := $(CFLAGSF90) -mtune=opteron -mcpu=opteron -march=opteron
CFLAGSC       := $(CFLAGSC) -mtune=opteron -mcpu=opteron -march=opteron
LDFLAGS       := $(LDFLAGS) -mtune=opteron -mcpu=opteron -march=opteron
endif
# Probably activate 64 bit integer support
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -i8
CFLAGSF90     := $(CFLAGSF90) -i8
endif
endif


ifeq ($(call match,$(ID),pc64-opteron-linux-sunstudio-.*),yes)
# invokes the Sun Studio Compiler suite, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/sunstudio-generic.mk
# # # # # # # # # # #      WARNING !!!     # # # # # # # # # # #
# Architecture-specific optimisation missing!
# DOM: sunstudio seems to do this automatically.
# # # # # # # # # # #      WARNING !!!     # # # # # # # # # # #
#CFLAGSC       := $(CFLAGSC) -xc99
CFLAGSF77     := $(CFLAGSF77) -m64 -xmodel=medium
CFLAGSF90     := $(CFLAGSF90) -m64 -xmodel=medium
CFLAGSC       := $(CFLAGSC) -m64 -xmodel=medium
LDFLAGS       := $(LDFLAGS) -m64 -xmodel=medium
# Probably activate 64 bit integer support
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -xtypemap=integer:64
CFLAGSF90     := $(CFLAGSF90) -xtypemap=integer:64
endif
endif


ifeq ($(call match,$(ID),pc64-opteron-linux-.*-.*-mpich),yes)
# Use MPICH
include $(TEMPLATESDIR)/mpich-generic.mk
endif

ifeq ($(call match,$(ID),pc64-opteron-linux-.*-.*-mpich2),yes)
# Use MPICH2
include $(TEMPLATESDIR)/mpich2-generic.mk
endif

ifeq ($(call match,$(ID),pc64-opteron-linux-.*-.*-mvapich),yes)
# Use MVAPICH
include $(TEMPLATESDIR)/mvapich-generic.mk
endif

ifeq ($(call match,$(ID),pc64-opteron-linux-.*-.*-lammpi),yes)
# Use LAM/MPI
include $(TEMPLATESDIR)/lammpi-generic.mk
endif

ifeq ($(call match,$(ID),pc64-opteron-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc64-opteron-linux-.*-.*-optmpich),yes)
# Use Allinea Opt MPI
include $(TEMPLATESDIR)/optmpich-generic.mk
# # # # # # # # # # #      WARNING !!!     # # # # # # # # # # #
#
# Code compiles, but does not run with OPT 1.0rc4. Dynamic linker
# can not find libopt.so. (Not a FEAT2 problem, though.)
#
# Additionally, OPT 1.0rc4 needs additional arguments in MPI calls.
# (Not sure about OPT 1.0 and above!!) parallel.f90 and parallelsys.f90
# are not prepared for them any more. (SB, 2006-06-02)
#
# # # # # # # # # # #      WARNING !!!     # # # # # # # # # # #
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -assume 2underscore
CFLAGSF77     := $(CFLAGSF77) -assume 2underscore
CFLAGSF90     := $(filter-out -DMPICALLS_NO_UNDERSCORE, $(CFLAGSF90)) \
		-DMPICALLS_NO_UNDERSCORE -assume 2underscore
SRCEXTRA      := $(SRCEXTRA) mpir_iargc.f

# Previous settings
# CFLAGSF77 = -g -us -assume 2underscore
# CFLAGSF90 = $(CFLAGSF77) -module $(OBJDIR) -check bounds -traceback \
#             -DHAS_INTRINSIC_IARGC -DHAS_FLUSH -DMPICALLS_NO_UNDERSCORE
# CFLAGSC   = -g -fpstkchk
# LDFLAGS   =
# SRCEXTRA  := $(SRCEXTRA) optmpi_wrapper.c
# MPIINC    = -I/usr/local/opt/mpich/include
# MPILIBDIR = -L/usr/local/opt/mpich/lib
# MPILIBS   = -L/usr/local/opt/opt/lib -lopt -lfmpich -lmpich -lpthread -lrt
# #MPILIBS   = -lpmpich -lfmpich -lmpich -lpthread -lrt
endif


ifeq ($(call match,$(ID),pc64-opteron-linux-.*-acml.*),yes)
# Use AMD Core Math Library
include $(TEMPLATESDIR)/acml-generic.mk
endif


ifeq ($(call match,$(ID),pc64-opteron-linux-.*-atlas.*),yes)
# Use ATLAS BLAS
include $(TEMPLATESDIR)/atlasblas-generic.mk
endif


ifeq ($(call match,$(ID),pc64-opteron-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc64-opteron-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc64-opteron-linux-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


ifeq ($(call match,$(ID),pc64-opteron-linux-.*-mkl.*),yes)
# Use Intel Math Kernel Library
include $(TEMPLATESDIR)/mkl-generic.mk
endif


ifeq ($(call match,$(ID),pc64-opteron-linux-.*-.*),yes)
# Add preprocessor flag to enable 64bit UMFPACK support. 
# UMFPACK documentation states to also add -DLP64, but there is no 
# reference to it in the UMFPACK 5.2 sources.
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF90     := -DENABLE_64BIT_UMFPACK $(CFLAGSF90)
CFLAGSC       := -DENABLE_64BIT_UMFPACK $(CFLAGSC)
endif
endif


###################### Opteron 64bit running 32bit app #######################

ifeq ($(call match,$(ID),pc64-opteron_ia32-cygwin_nt.*-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
CFLAGSF77LIBS := -m32 $(CFLAGSF77LIBS)
CFLAGSF77     := -m32 $(CFLAGSF77)
CFLAGSF90     := -m32 $(CFLAGSF90)
CFLAGSC       := -m32 $(CFLAGSC)
LDFLAGS       := -m32 $(LDFLAGS)
ifeq ($(call optimise), YES)
# -mmx -m3dnow -mfpmath=sse -msse -msse2 are automatically set by -march=opteron
CFLAGSF77LIBS := -march=opteron -mfpmath=sse -msse2 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=opteron -mfpmath=sse -msse2 $(CFLAGSF77)
CFLAGSF90     := -march=opteron -mfpmath=sse -msse2 $(CFLAGSF90)
CFLAGSC       := -march=opteron -mfpmath=sse -msse2 $(CFLAGSC)
endif
endif


ifeq ($(call match,$(ID),pc64-opteron_ia32-linux-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
CFLAGSF77LIBS := -m32 $(CFLAGSF77LIBS)
CFLAGSF77     := -m32 $(CFLAGSF77)
CFLAGSF90     := -m32 $(CFLAGSF90)
CFLAGSC       := -m32 $(CFLAGSC)
LDFLAGS       := -m32 $(LDFLAGS)
ifeq ($(call optimise), YES)
# -mmx -m3dnow -mfpmath=sse -msse -msse2 are automatically set by -march=opteron
CFLAGSF77LIBS := -march=opteron -mfpmath=sse -msse2 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=opteron -mfpmath=sse -msse2 $(CFLAGSF77)
CFLAGSF90     := -march=opteron -mfpmath=sse -msse2 $(CFLAGSF90)
CFLAGSC       := -march=opteron -mfpmath=sse -msse2 $(CFLAGSC)
endif
endif


ifeq ($(call match,$(ID),pc64-opteron_ia32-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc64-opteron_ia32-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc64-opteron_ia32-linux-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


################################# Pentium 3 ##################################

ifeq ($(call match,$(ID),pc-pentium3-linux-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium3 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium3 $(CFLAGSF77)
CFLAGSF90     := -march=pentium3 $(CFLAGSF90)
CFLAGSC       := -march=pentium3 $(CFLAGSC)
endif
endif


ifeq ($(call match,$(ID),pc-pentium3-cygwin_nt.*-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium3 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium3 $(CFLAGSF77)
CFLAGSF90     := -march=pentium3 $(CFLAGSF90)
CFLAGSC       := -march=pentium3 $(CFLAGSC)
LDFLAGS       := -march=pentium3 $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-pentium3-linux-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium3 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium3 $(CFLAGSF77)
CFLAGSF90     := -march=pentium3 $(CFLAGSF90)
CFLAGSC       := -march=pentium3 $(CFLAGSC)
LDFLAGS       := -march=pentium3 $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-pentium3-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -xW $(CFLAGSF77LIBS)
CFLAGSF77     := -xW $(CFLAGSF77)
CFLAGSF90     := -xW $(CFLAGSF90)
CFLAGSC       := -xW $(CFLAGSC)
CFLAGSCXX     := -xW $(CFLAGSCXX)
LDFLAGS       := -xW $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-pentium3-linux-.*-.*-lammpi),yes)
# Use LAM/MPI
include $(TEMPLATESDIR)/lammpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-pentium3-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-pentium3-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc-pentium3-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc-pentium3-linux-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


################################# Pentium M ##################################

ifeq ($(call match,$(ID),pc-pentiumm-linux-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium-m $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium-m $(CFLAGSF77)
CFLAGSF90     := -march=pentium-m $(CFLAGSF90)
CFLAGSC       := -march=pentium-m $(CFLAGSC)
endif
endif


ifeq ($(call match,$(ID),pc-pentiumm-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-pentiumm-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc-pentiumm-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc-pentiumm-linux-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


################################# Pentium 4 ##################################

ifeq ($(call match,$(ID),pc-pentium4-cygwin_nt.*-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium4 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium4 $(CFLAGSF77)
CFLAGSF90     := -march=pentium4 $(CFLAGSF90)
CFLAGSC       := -march=pentium4 $(CFLAGSC)
LDFLAGS       := -march=pentium4 $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-pentium4-linux-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium4 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium4 $(CFLAGSF77)
CFLAGSF90     := -march=pentium4 $(CFLAGSF90)
CFLAGSC       := -march=pentium4 $(CFLAGSC)
LDFLAGS       := -march=pentium4 $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-pentium4-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -xN $(CFLAGSF77LIBS) #-tpp6
CFLAGSF77     := -xN $(CFLAGSF77)     #-tpp6
CFLAGSF90     := -xN $(CFLAGSF90)     #-tpp6
CFLAGSC       := -xN $(CFLAGSC)       #-tpp6
LDFLAGS       := -xN $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-pentium4-linux-.*-.*-lammpi),yes)
# Use LAM/MPI
include $(TEMPLATESDIR)/lammpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-pentium4-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-pentium4-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc-pentium4-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc-pentium4-linux-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


############################# Pentium 4 (Nocona) #############################

ifeq ($(call match,$(ID),pc64-pentium4-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -xP $(CFLAGSF77LIBS)
CFLAGSF77     := -xP $(CFLAGSF77)
CFLAGSF90     := -xP $(CFLAGSF90)
CFLAGSC       := -xP $(CFLAGSC)
LDFLAGS       := -xP $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc64-pentium4-linux-.*-.*-lanlmpi),yes)
# Use LANL MPI
include $(TEMPLATESDIR)/lanlmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc64-pentium4-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),pc64-pentium4-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc64-pentium4-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc64-pentium4-linux-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


############################## Pentium 4 Mobile ##############################

ifeq ($(call match,$(ID),pc-pentium4m-cygwin_nt.*-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium4m $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium4m $(CFLAGSF77)
CFLAGSF90     := -march=pentium4m $(CFLAGSF90)
CFLAGSC       := -march=pentium4m $(CFLAGSC)
LDFLAGS       := -march=pentium4m $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-pentium4m-linux-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium4m $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium4m $(CFLAGSF77)
CFLAGSF90     := -march=pentium4m $(CFLAGSF90)
CFLAGSC       := -march=pentium4m $(CFLAGSC)
LDFLAGS       := -march=pentium4m $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-pentium4m-linux-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium4m $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium4m $(CFLAGSF77)
CFLAGSF90     := -march=pentium4m $(CFLAGSF90)
CFLAGSC       := -march=pentium4m $(CFLAGSC)
endif
endif


ifeq ($(call match,$(ID),pc-pentium4m-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -xB $(CFLAGSF77LIBS)
CFLAGSF77     := -xB $(CFLAGSF77)
CFLAGSF90     := -xB $(CFLAGSF90)
CFLAGSC       := -xB $(CFLAGSC)
LDFLAGS       := -xB $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-pentium4m-linux-sunstudio-.*),yes)
# invokes the Sun Studio Compiler suite, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/sunstudio-generic.mk
#CFLAGSC       := -xc99 $(CFLAGSC)
endif


ifeq ($(call match,$(ID),pc-pentium4m-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc-pentium4m-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc-pentium4m-linux-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


ifeq ($(call match,$(ID),pc-pentium4m-linux-.*-.*-lammpi),yes)
# Use LAM/MPI
include $(TEMPLATESDIR)/lammpi-generic.mk
endif


ifeq ($(call match,$(ID),pc-pentium4m-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif



######################### Nehalem 64 bit (Xeon X5560) #########################
#########################   Nehalem 64 bit (HP BETA)  #########################
#########################        Core i7 64 bit       #########################

ifeq ($(call match,$(ID),pc64-nehalem-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack, local metis
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(strip $(optimise)), YES)
# older intel versions
ifneq (,$(findstring 10.1.,$(INTELVERSION)))
CFLAGSF77LIBS := -xT $(CFLAGSF77LIBS)
CFLAGSF77     := -xT $(CFLAGSF77)
CFLAGSF90     := -xT $(CFLAGSF90)
CFLAGSC       := -xT $(CFLAGSC)
CFLAGSCXX     := -xT $(CFLAGSCXX)
CFLAGSCOPROC  := -xT $(CFLAGSCOPROC)
LDFLAGS       := -xT $(LDFLAGS)
endif
# intel 11.1
ifneq (,$(findstring 11.1,$(INTELVERSION)))
CFLAGSF77LIBS := -xSSE4.2 $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE4.2 $(CFLAGSF77)
CFLAGSF90     := -xSSE4.2 $(CFLAGSF90)
CFLAGSC       := -xSSE4.2 $(CFLAGSC)
CFLAGSCXX     := -xSSE4.2 $(CFLAGSCXX)
CFLAGSCOPROC  := -xSSE4.2 $(CFLAGSCOPROC)
LDFLAGS       := -xSSE4.2 $(LDFLAGS)
endif
endif
endif


ifeq ($(call match,$(ID),pc64-nehalem-linux-gcc-.*),yes)
# invokes the GCC Fortran 77/90 compiler, local umfpack, local metis
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(strip $(optimise)), YES)
CFLAGSF77LIBS := -O3 -m64 -march=core2 -mtune=core2 -msse4 -msse4.1 -msse4.2 -mfpmath=sse $(CFLAGSF77LIBS)
CFLAGSF77     := -O3 -m64 -march=core2 -mtune=core2 -msse4 -msse4.1 -msse4.2 -mfpmath=sse $(CFLAGSF77)
CFLAGSF90     := -O3 -m64 -march=core2 -mtune=core2 -msse4 -msse4.1 -msse4.2 -mfpmath=sse $(CFLAGSF90)
CFLAGSC       := -O3 -m64 -march=core2 -mtune=core2 -msse4 -msse4.1 -msse4.2 -mfpmath=sse $(CFLAGSC)
CFLAGSCXX     := -O3 -m64 -march=core2 -mtune=core2 -msse4 -msse4.1 -msse4.2 -mfpmath=sse $(CFLAGSCXX)
CFLAGSCOPROC  := -O3 -m64 -march=core2 -mtune=core2 -msse4 -msse4.1 -msse4.2 -mfpmath=sse $(CFLAGSCOPROC)
endif
endif


#sun (ideas, not tested yet):
#-fast -xtarget=nehalem -xarch=generic 


ifeq ($(call match,$(ID),pc64-nehalem-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc64-nehalem-linux-.*-mkl.*),yes)
# Use Intel Math Kernel Library
include $(TEMPLATESDIR)/mkl-generic.mk
endif


ifeq ($(call match,$(ID),pc64-nehalem-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif


ifeq ($(call match,$(ID),pc64-nehalem-linux-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


ifeq ($(call match,$(ID),pc64-nehalem-linux-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif

######################### Penryn 32 bit #########################

ifeq ($(call match,$(ID),pc-penryn-cygwin_nt.*-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack, local metis
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
# g95 is based on gcc 4.0.3. So, it supports less processors
# than gcc 4.3 and 4.4. For all available options see
# http://gcc.gnu.org/onlinedocs/gcc-4.0.3/gcc.pdf, page 145.
ifeq ($(strip $(optimise)), YES)
CFLAGSF77LIBS := -march=nocona $(CFLAGSF77LIBS)
CFLAGSF77     := -march=nocona $(CFLAGSF77)
CFLAGSF90     := -march=nocona $(CFLAGSF90)
CFLAGSC       := -march=nocona $(CFLAGSC)
CFLAGSCXX     := -march=nocona $(CFLAGSCXX)
CFLAGSCOPROC  := -march=nocona $(CFLAGSCOPROC)
LDFLAGS       := -march=nocona $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-penryn-darwin-g95-.*),yes)
# invokes the GNU Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -m32 -march=nocona $(CFLAGSF77LIBS)
CFLAGSF77     := -m32 -march=nocona $(CFLAGSF77)
CFLAGSF90     := -m32 -march=nocona $(CFLAGSF90)
CFLAGSC       := -m32 -march=nocona $(CFLAGSC)
CFLAGSCXX     := -m32 -march=nocona $(CFLAGSCXX)
LDFLAGS       := -m32 -march=nocona -read_only_relocs suppress $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-penryn-darwin-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
gcc43or44=NO
ifneq (,$(findstring 4.3.,$(GFORTRANVERSION)))
gcc43or44=YES
endif
ifneq (,$(findstring 4.4.,$(GFORTRANVERSION)))
gcc43or44=YES
endif
ifeq ($(strip $(gcc43or44)), YES)
ifeq ($(call optimise), YES)
# Options choosen based on 
# * the article at http://www.outsourcingi.com/articles_content.php?id=454
# * 'info gcc', section 3.17.14 Intel 386 and AMD x86-64 Options
# * and the options gcc itself passes when '-march=native --verbose' is 
#   invoked on a MacBookPro5,5.
# Note: It is absolutely mandatory to use a recent binutils package!
#       Otherwise SSE4.1 is not available and the compiler will bail out.
VERBOSEOPTFLAGS = \
         -m32 -march=core2 -msse4.1 -mcx16 -msahf --param l1-cache-size=32 \
         --param l1-cache-line-size=64 --param l2-cache-size=3072 \
	 -mtune=core2 -ftree-vectorize

CFLAGSF77LIBS := $(VERBOSEOPTFLAGS) $(CFLAGSF77LIBS)
CFLAGSF77     := $(VERBOSEOPTFLAGS) $(CFLAGSF77)
CFLAGSF90     := $(VERBOSEOPTFLAGS) $(CFLAGSF90)
CFLAGSC       := $(VERBOSEOPTFLAGS) $(CFLAGSC)
CFLAGSCXX     := $(VERBOSEOPTFLAGS) $(CFLAGSCXX)
CFLAGSCOPROC  := $(VERBOSEOPTFLAGS) $(CFLAGSCOPROC)
LDFLAGS       := $(VERBOSEOPTFLAGS) -read_only_relocs suppress $(LDFLAGS)
endif
else
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -m32 -march=nocona $(CFLAGSF77LIBS)
CFLAGSF77     := -m32 -march=nocona $(CFLAGSF77)
CFLAGSF90     := -m32 -march=nocona $(CFLAGSF90)
CFLAGSC       := -m32 -march=nocona $(CFLAGSC)
CFLAGSCXX     := -m32 -march=nocona $(CFLAGSCXX)
LDFLAGS       := -m32 -march=nocona -read_only_relocs suppress $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc-penryn-darwin-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif

ifeq ($(call match,$(ID),pc-penryn-darwin-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif

ifeq ($(call match,$(ID),pc-penryn-darwin-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif

######################### Penryn 64 bit (Xeon E5450) #########################
######################### Penryn 64 bit (LiDO 2)     #########################

ifeq ($(call match,$(ID),pc64-penryn-linux-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack, local metis
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
# g95 is based on gcc 4.0.3. So, it supports less processors
# than gcc 4.3 and 4.4. For all available options see
# http://gcc.gnu.org/onlinedocs/gcc-4.0.3/gcc.pdf, page 145.
ifeq ($(strip $(optimise)), YES)
CFLAGSF77LIBS := -march=nocona $(CFLAGSF77LIBS)
CFLAGSF77     := -march=nocona $(CFLAGSF77)
CFLAGSF90     := -march=nocona $(CFLAGSF90)
CFLAGSC       := -march=nocona $(CFLAGSC)
CFLAGSCXX     := -march=nocona $(CFLAGSCXX)
CFLAGSCOPROC  := -march=nocona $(CFLAGSCOPROC)
LDFLAGS       := -march=nocona $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc64-penryn-linux-gcc-.*),yes)
# invokes the GNU Compiler suite, local umfpack, local metis
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
gcc43or44=NO
ifneq (,$(findstring 4.3.,$(GFORTRANVERSION)))
gcc43or44=YES
endif
ifneq (,$(findstring 4.4.,$(GFORTRANVERSION)))
gcc43or44=YES
endif
ifeq ($(strip $(gcc43or44)), YES)
ifeq ($(strip $(optimise)), YES)
# Options choosen based on 
# * the article at http://www.outsourcingi.com/articles_content.php?id=454
# * 'info gcc', section 3.17.14 Intel 386 and AMD x86-64 Options
# * and the options gcc itself passes when '-march=native --verbose' is 
#   invoked on the LiDO-2 eth/ib compute nodes.
# Note: It is absolutely mandatory to use a recent binutils package!
#       Otherwise SSE4.1 is not available and the compiler will bail out.
LIDO2COMPUTENODESOPTFLAGS = \
         -march=core2 -msse4.1 -mcx16 -msahf --param l1-cache-size=32 \
         --param l1-cache-line-size=64 --param l2-cache-size=6144 \
         -ftree-vectorize
CFLAGSF77LIBS := $(LIDO2COMPUTENODESOPTFLAGS) $(CFLAGSF77LIBS)
CFLAGSF77     := $(LIDO2COMPUTENODESOPTFLAGS) $(CFLAGSF77)
CFLAGSF90     := $(LIDO2COMPUTENODESOPTFLAGS) $(CFLAGSF90)
CFLAGSC       := $(LIDO2COMPUTENODESOPTFLAGS) $(CFLAGSC)
CFLAGSCXX     := $(LIDO2COMPUTENODESOPTFLAGS) $(CFLAGSCXX)
CFLAGSCOPROC  := $(LIDO2COMPUTENODESOPTFLAGS) $(CFLAGSCOPROC)
LDFLAGS       := $(LIDO2COMPUTENODESOPTFLAGS) $(LDFLAGS)
endif
endif
endif


ifeq ($(call match,$(ID),pc64-penryn-linux-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifneq (,$(findstring 10.1.,$(INTELVERSION)))
CFLAGSF77LIBS := -xT $(CFLAGSF77LIBS)
CFLAGSF77     := -xT $(CFLAGSF77)
CFLAGSF90     := -xT $(CFLAGSF90)
CFLAGSC       := -xT $(CFLAGSC)
CFLAGSCXX     := -xT $(CFLAGSCXX)
LDFLAGS       := -xT $(LDFLAGS)
endif
ifneq (,$(findstring 11.0.0,$(INTELVERSION)))
CFLAGSF77LIBS := -xSSE4.1 $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE4.1 $(CFLAGSF77)
CFLAGSF90     := -xSSE4.1 $(CFLAGSF90)
CFLAGSC       := -xSSE4.1 $(CFLAGSC)
CFLAGSCXX     := -xSSE4.1 $(CFLAGSCXX)
LDFLAGS       := -xSSE4.1 $(LDFLAGS)
endif
endif
endif


ifeq ($(call match,$(ID),pc64-penryn-linux-psc-.*),yes)
# invokes the PathScale Fortran 77/90 compiler, local umfpack, local metis
MESSAGE  := $(MESSAGE) \
	    echo; \
	    echo '*** Warning: SB wants to state: Not yet sure whether optimisation flags are already optimal!'; \
	    echo;
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pathscale-generic.mk
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -m64 -march=auto
CFLAGSF77     := $(CFLAGSF77) -m64 -march=auto
CFLAGSF90     := $(CFLAGSF90) -m64 -march=auto
CFLAGSC       := $(CFLAGSC) -m64 -march=auto
CFLAGSCXX     := $(CFLAGSCXX) -m64 -march=auto
CFLAGSCOPROC  := $(CFLAGSCOPROC) -m64 -march=auto
LDFLAGS       := $(LDFLAGS) -m64 -march=auto
# default is -march=auto
ifeq ($(strip $(optimise)), YES)
# -msse4a
endif
endif


ifeq ($(call match,$(ID),pc64-penryn-linux-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
MESSAGE  := $(MESSAGE) \
	    echo; \
	    echo '*** Warning: SB wants to state: Not yet sure whether optimisation flags are already optimal!'; \
	    echo;
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -tp penryn-64
CFLAGSF77     := $(CFLAGSF77) -tp=penryn-64
CFLAGSF90     := $(CFLAGSF90) -tp=penryn-64
CFLAGSC       := $(CFLAGSC) -tp=penryn-64
LDFLAGS       := $(LDFLAGS) -tp=penryn-64
endif
# Probably activate 64 bit integer support
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSF77     := $(CFLAGSF77) -i8
CFLAGSF90     := $(CFLAGSF90) -i8
endif
endif


ifeq ($(call match,$(ID),pc64-penryn-linux-sunstudio-.*),yes)
# invokes the Sun Studio Compiler suite, local umfpack, local metis
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/sunstudio-generic.mk
ifeq ($(strip $(optimise)), YES)
# Options choosen based on 
# * the options the sun compiler themselves pass when '-fast -dryrun' is 
#   invoked on the LiDO-2 eth/ib compute nodes.
LIDO2COMPUTENODESOPTFLAGS = \
          -xO5 -xarch=sse4_1 -xcache=32/64/8:6144/64/24 -xchip=penryn \
          -fsimple=2 -fns=yes -ftrap=common -xlibmil -nofstore -xregs=frameptr
CFLAGSF77LIBS := $(CFLAGSF77LIBS) $(LIDO2COMPUTENODESOPTFLAGS) -dalign
CFLAGSF77     := $(CFLAGSF77) $(LIDO2COMPUTENODESOPTFLAGS) -dalign
CFLAGSF90     := $(CFLAGSF90) $(LIDO2COMPUTENODESOPTFLAGS) -dalign
CFLAGSC       := $(CFLAGSC) $(LIDO2COMPUTENODESOPTFLAGS)
CFLAGSCXX     := $(CFLAGSCXX) $(LIDO2COMPUTENODESOPTFLAGS)
CFLAGSCOPROC  := $(CFLAGSCOPROC) $(LIDO2COMPUTENODESOPTFLAGS)
LDFLAGS       := $(LDFLAGS) $(LIDO2COMPUTENODESOPTFLAGS) -dalign
endif
endif


ifeq ($(call match,$(ID),pc64-penryn-linux-.*-acml.*),yes)
# Use AMD Core Math Library
include $(TEMPLATESDIR)/acml-generic.mk
endif


ifeq ($(call match,$(ID),pc64-penryn-linux-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),pc64-penryn-linux-.*-mkl.*),yes)
# Use Intel Math Kernel Library
include $(TEMPLATESDIR)/mkl-generic.mk
endif


ifeq ($(call match,$(ID),pc64-penryn-linux-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif

ifeq ($(call match,$(ID),pc64-penryn-linux-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif

ifeq ($(call match,$(ID),pc64-penryn-darwin-g95-.*),yes)
# invokes the GNU Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -m64 -march=nocona $(CFLAGSF77LIBS)
CFLAGSF77     := -m64 -march=nocona $(CFLAGSF77)
CFLAGSF90     := -m64 -march=nocona $(CFLAGSF90)
CFLAGSC       := -m64 -march=nocona $(CFLAGSC)
CFLAGSCXX     := -m64 -march=nocona $(CFLAGSCXX)
LDFLAGS       := -m64 -march=nocona -read_only_relocs suppress $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-penryn-darwin-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
gcc43or44=NO
ifneq (,$(findstring 4.3.,$(GFORTRANVERSION)))
gcc43or44=YES
endif
ifneq (,$(findstring 4.4.,$(GFORTRANVERSION)))
gcc43or44=YES
endif
ifeq ($(strip $(gcc43or44)), YES)
ifeq ($(call optimise), YES)
# Options choosen based on 
# * the article at http://www.outsourcingi.com/articles_content.php?id=454
# * 'info gcc', section 3.17.14 Intel 386 and AMD x86-64 Options
# * and the options gcc itself passes when '-march=native --verbose' is 
#   invoked on a MacBookPro5,5.
# Note: It is absolutely mandatory to use a recent binutils package!
#       Otherwise SSE4.1 is not available and the compiler will bail out.
VERBOSEOPTFLAGS = \
         -m64 -march=core2 -msse4.1 -mcx16 -msahf --param l1-cache-size=32 \
         --param l1-cache-line-size=64 --param l2-cache-size=3072 \
	 -mtune=core2 -ftree-vectorize

CFLAGSF77LIBS := $(VERBOSEOPTFLAGS) $(CFLAGSF77LIBS)
CFLAGSF77     := $(VERBOSEOPTFLAGS) $(CFLAGSF77)
CFLAGSF90     := $(VERBOSEOPTFLAGS) $(CFLAGSF90)
CFLAGSC       := $(VERBOSEOPTFLAGS) $(CFLAGSC)
CFLAGSCXX     := $(VERBOSEOPTFLAGS) $(CFLAGSCXX)
CFLAGSCOPROC  := $(VERBOSEOPTFLAGS) $(CFLAGSCOPROC)
LDFLAGS       := $(VERBOSEOPTFLAGS) -read_only_relocs suppress $(LDFLAGS)
endif
else
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -m64 -march=nocona $(CFLAGSF77LIBS)
CFLAGSF77     := -m64 -march=nocona $(CFLAGSF77)
CFLAGSF90     := -m64 -march=nocona $(CFLAGSF90)
CFLAGSC       := -m64 -march=nocona $(CFLAGSC)
CFLAGSCXX     := -m64 -march=nocona $(CFLAGSCXX)
LDFLAGS       := -m64 -march=nocona -read_only_relocs suppress $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc64-penryn-darwin-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif

ifeq ($(call match,$(ID),pc64-penryn-darwin-.*-goto.*),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif

ifeq ($(call match,$(ID),pc64-penryn-darwin-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif


############################### Sun Ultrasparc v8 #############################

ifeq ($(call match,$(ID),sun4u-sparcv8-sunos-g95-blas.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
# Umfpack install instructions advice to set this flag on Sun.
CFLAGSC       := -DNSUNPERF $(CFLAGSC)
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -mcpu=ultrasparc $(CFLAGSF77LIBS)
CFLAGSF77     := -mcpu=ultrasparc $(CFLAGSF77)
CFLAGSF90     := -mcpu=ultrasparc $(CFLAGSF90)
CFLAGSC       := -mcpu=ultrasparc $(CFLAGSC)
LDFLAGS       := -mcpu=ultrasparc $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),sun4u-sparcv8-sunos-g95-perf.*),yes)
# invokes the GNU g95 compiler, uses LAM MPI, Perflib BLAS,
# local umfpack
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
include $(TEMPLATESDIR)/sunperf-generic.mk
CFLAGSC       := $(CFLAGSC) #-DNSUNPERF
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -mcpu=ultrasparc $(CFLAGSF77LIBS)
CFLAGSF77     := -mcpu=ultrasparc $(CFLAGSF77)
CFLAGSF90     := -mcpu=ultrasparc $(CFLAGSF90)
CFLAGSC       := -mcpu=ultrasparc $(CFLAGSC)
LDFLAGS       := -mcpu=ultrasparc $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),sun4u-sparcv8-sunos-.*-.*-lammpi),yes)
# Use LAM/MPI
include $(TEMPLATESDIR)/lammpi-generic.mk
endif


ifeq ($(call match,$(ID),sun4u-sparcv8-sunos-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),sun4u-sparcv8-sunos-sunstudio-.*),yes)
# invokes the Sun Studio Compiler suite, uses LAM MPI, local BLAS,
# local umfpack
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/sunstudio-generic.mk
# Umfpack install instructions advice to set this flag on Sun.
CFLAGSC       := -DNSUNPERF $(CFLAGSC)
# Have 32bit binaries built on Sun
CFLAGSF77LIBS := -xarch=native $(CFLAGSF77LIBS)
CFLAGSF77     := -xarch=native $(CFLAGSF77)
CFLAGSF90     := -xarch=native $(CFLAGSF90)
CFLAGSC       := -xarch=native $(CFLAGSC)
LDFLAGS       := -xarch=native $(LDFLAGS)
# -lmvec Needed? Performance benefit? I don't know (SB)
#
# If using -lmvec (some Vectorised Mathematical function) we also
# need socket and nsl libraries.
LIBS          := $(LIBS) -lmvec -lsocket -lnsl
# With static compiling, don't try to link against pthread library
# By the way, both dynamic and static compiling works! Just comment
# out the following three lines to get dynamic compiling again.
LDFLAGS       := $(LDFLAGS) -Bstatic
LIBS          := $(filter-out -lpthread, $(LIBS))
LIBS          := $(filter-out -ldl, $(LIBS))
MPILIBS       := $(filter-out -lpthread, $(MPILIBS))
MPILIBS       := $(filter-out -ldl, $(MPILIBS))
endif


ifeq ($(call match,$(ID),sun4u-sparcv8-sunos-sunstudio-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),sun4u-sparcv8-sunos-sunstudio-perf.*),yes)
# Use Perflib BLAS,
include $(TEMPLATESDIR)/sunperf-generic.mk
# Is it necessary to replace the way Sunperf Library is included?
# It works both ways! (SB)
LIBS          := $(filter-out -lsunperf, $(LIBS)) -xlic_lib=sunperf
endif



############################### Sun Ultrasparc v9 #############################

ifeq ($(call match,$(ID),sun4u-sparcv9-sunos-sunstudio-.*),yes)
# invokes the Sun Studio Compiler suite, local umfpack
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/sunstudio-generic.mk
# Umfpack install instructions advice to set this flag on Sun.
CFLAGSC       := -DNSUNPERF $(CFLAGSC)
# Have 64bit binaries built on Sun
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -xarch=native64
CFLAGSF77     := $(CFLAGSF77)     -xarch=native64
CFLAGSF90     := $(CFLAGSF90)     -xarch=native64
CFLAGSC       := $(CFLAGSC)       -xarch=native64
LDFLAGS       := $(LDFLAGS)       -xarch=native64
# If using -lmvec (some Vectorised Mathematical function) we also
# need socket and nsl libraries.
LIBS          := $(LIBS) -lmvec -lsocket -lnsl
endif


ifeq ($(call match,$(ID),sun4u-sparcv9-sunos-.*-.*-lammpi),yes)
# Use LAM/MPI
include $(TEMPLATESDIR)/lammpi-generic.mk
endif


ifeq ($(call match,$(ID),sun4u-sparcv9-sunos-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),sun4u-sparcv9-sunos-sunstudio-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif


ifeq ($(call match,$(ID),sun4u-sparcv9-sunos-sunstudio-perf.*),yes)
# Use Perflib BLAS,
include $(TEMPLATESDIR)/sunperf-generic.mk
# Is it necessary to replace the way Sunperf Library is included?
# It works both ways! (SB)
LIBS          := $(filter-out -lsunperf, $(LIBS)) -xlic_lib=sunperf
endif



########################## Sun Whatever Solaris 10 ###########################

ifeq ($(call match,$(ID),sun4v-sparcv9-sunos-sunstudio-.*),yes)
# invokes the Sun Studio Compiler suite, local umfpack
#
# not fully tested settings
MESSAGE := $(MESSAGE) \
	    echo; \
	    echo '*** Warning: The settings for this build ID have not been thoroughly tested yet'; \
	    echo;
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/sunstudio-generic.mk
# Umfpack install instructions advice to set this flag on Sun.
CFLAGSC       := -DNSUNPERF $(CFLAGSC)
# Have 64bit binaries built on Sun
CFLAGSF77LIBS := $(CFLAGSF77LIBS) -xarch=native64
CFLAGSF77     := $(CFLAGSF77)     -xarch=native64
CFLAGSF90     := $(CFLAGSF90)     -xarch=native64
CFLAGSC       := $(CFLAGSC)       -xarch=native64
LDFLAGS       := $(LDFLAGS)       -xarch=native64
# If using -lmvec (some Vectorised Mathematical function) we also
# need socket and nsl libraries.
LIBS          := $(LIBS) -lmvec -lsocket -lnsl
# Is it necessary to replace the way Sunperf Library is included?
# It works both ways! (SB)
LIBS          := $(filter-out -lsunperf, $(LIBS))
LIBS          := $(LIBS) -xlic_lib=sunperf
endif


ifeq ($(call match,$(ID),sun4v-sparcv9-sunos-gcc-.*),yes)
# invokes the GNU Compiler suite, local umfpack
#
# not fully tested settings
MESSAGE := $(MESSAGE) \
	    echo; \
	    echo '*** Warning: The settings for this build ID have not been thoroughly tested yet'; \
	    echo;
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
# Umfpack install instructions advice to set this flag on Sun.
CFLAGSC       := -DNSUNPERF $(CFLAGSC)
ifeq ($(call optimise), YES)
#CFLAGSF77LIBS := -march=athlon64 $(CFLAGSF77LIBS)
#CFLAGSF77     := -march=athlon64 $(CFLAGSF77)
#CFLAGSF90     := -march=athlon64 $(CFLAGSF90)
#CFLAGSC       := -march=athlon64 $(CFLAGSC)
#CFLAGSCXX     := -march=athlon64 $(CFLAGSCXX)
#LDFLAGS       := -march=athlon64 $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),sun4v-sparcv9-sunos-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif


ifeq ($(call match,$(ID),sun4v-sparcv9-sunos-.*-perf.*),yes)
# Use Perflib BLAS
include $(TEMPLATESDIR)/sunperf-generic.mk
endif



#################################### NEC #####################################

ifeq ($(call match,$(ID),sx6-none-superux-f90-keisan-mpi),yes)
TOKEN5 := 1
TOKEN6 := 1
include $(TEMPLATESDIR)/necsx-generic.mk
include $(TEMPLATESDIR)/necsxcompiler-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -sx6 $(CFLAGSF77LIBS)
CFLAGSF77     := -sx6 $(CFLAGSF77)
CFLAGSF90     := -sx6 $(CFLAGSF90)
CFLAGSC       := -sx6 $(CFLAGSC)
LDFLAGS       := -sx6 $(LDFLAGS)
endif
INC       = 
#BUILDLIB  = metis blas lapack umfpack amd
BUILDLIB  = metis umfpack amd
LIBS      = -lmetis -lumfpack -lamd -L/SX/opt/MathKeisan/lib -llapack -lblas
SBB_LIBS  = -L/SX/opt/MathKeisan/lib -llapack -lblas
MPILIBS   = -lmpi
endif


ifeq ($(call match,$(ID),sx8-none-superux-f90-.*-mpi),yes)
TOKEN6 := 1
include $(TEMPLATESDIR)/necsx-generic.mk
include $(TEMPLATESDIR)/necsxcompiler-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -sx8 $(CFLAGSF77LIBS) #-ftrace
CFLAGSF77     := -sx8 $(CFLAGSF77)  #-ftrace
CFLAGSF90     := -sx8 $(CFLAGSF90) #-Wf~-L fmtlist map summary transform~ -ftrace
CFLAGSC       := -sx8 $(CFLAGSC) #-ftrace
LDFLAGS       := -sx8 $(LDFLAGS) #-ftrace
endif
INC       =
SBB_LIBS  = -llapack -lblas
endif


ifeq ($(call match,$(ID),sx8-none-superux-f90-keisan-mpi),yes)
TOKEN5 := 1
#BUILDLIB  = metis blas lapack umfpack amd
BUILDLIB  = metis umfpack amd
#LIBDIR   := $(LIBDIR) -L/SX/opt/MathKeisan/lib
#LIBDIR   := $(LIBDIR) -L/SX/opt/mathkeisan/MK1_6/lib0/lib64
LIBS  := $(LIBS) -L/SX/opt/MathKeisan/lib -llapack -lblas
MPILIBS   = -lmpi
endif


ifeq ($(call match,$(ID),sx8-none-superux-f90-blas-mpi),yes)
TOKEN5 := 1
#BUILDLIB  = metis blas lapack umfpack amd
BUILDLIB  = metis umfpack amd blas lapack
#LIBDIR   := $(LIBDIR) -L/SX/usr/lib0/libp
MPILIBS   = -lmpi
MPILIBS   =
endif


##############################################################################
# auxiliary targets: .id, .idonly
#
##############################################################################
#

# hack to have this target in all Makefiles, the dot is to not
# consider it as a default rule when called without specific target
.PHONY: .id
.id:
	@echo 'Machine id:' $(ID)
	@echo
	@echo 'Compilers and preprocessors to be used:'
        # Catch cases where compiler/preprocessor is either not set or
        # compiler definition has been intentionally overwritten
        # to point the user to a problem.
	@if test -z "$(findstring exit 1, $(CPP))"; then \
	    echo '  C preprocessor     :' $(shell (which 2>/dev/null $(CPP) || \
		echo "\*\*\* $(CPP) not found ! \*\*\*")); \
	else \
	    echo '  C preprocessor     : *** not set ! ***'; \
	fi
	@if test -z "$(findstring exit 1, $(F90CPP))"; then \
	    echo '  FEAT2 preprocessor :' $(shell (which 2>/dev/null $(F90CPP) || \
		echo "\*\*\* $(F90CPP) not found ! \*\*\*")); \
	else \
	    echo '  FEAT2 preprocessor : *** not set ! ***'; \
	fi
	@if test -z "$(findstring exit 1, $(F77))"; then \
	    echo '  Fortran 77 compiler:' $(shell (which 2>/dev/null $(F77) || \
		echo "\*\*\* $(F77) not found ! \*\*\*")); \
	else \
	    echo '  Fortran 77 compiler: *** not set ! ***'; \
	fi
	@if test -z "$(findstring exit 1, $(F90))"; then \
	    echo '  Fortran 90 compiler:' $(shell (which 2>/dev/null $(F90) || \
		echo "\*\*\* $(F90) not found ! \*\*\*")); \
	else \
	    echo '  Fortran 90 compiler: *** not set ! ***'; \
	fi
	@if test -z "$(findstring exit 1, $(CC))"; then \
	    echo '  C compiler         :' $(shell (which 2>/dev/null $(CC)  || \
		echo "\*\*\* $(CC) not found ! \*\*\*")); \
	else \
	    echo '  C compiler         : *** not set ! ***'; \
	fi
	@if test -z "$(findstring exit 1, $(CXX))"; then \
	    echo '  C++ compiler       :' $(shell (which 2>/dev/null $(CXX) || \
		echo "\*\*\* $(CXX) not found ! \*\*\*")); \
	else \
	    echo '  C++ compiler       : *** not set ! ***'; \
	fi
	@echo
	@echo 'Linker to be used:'
	@if test -z "$(findstring exit 1, $(LD))"; then \
	    echo '  Linker             :' $(shell (which 2>/dev/null $(LD)  || \
		echo "\*\*\* $(LD) not found ! \*\*\*")); \
	else \
	    echo '  Linker             : *** not set ! ***'; \
	fi
	@echo
	@echo 'Flags to be used:'
	@echo '  OPT          =' $(OPT)
	@echo '  CFLAGSF77    =' $(CFLAGSF77)
	@echo '  CFLAGSF90    =' $(CFLAGSF90)
	@echo '  CFLAGSC      =' $(CFLAGSC)
	@echo '  CFLAGSCXX    =' $(CFLAGSCXX)
	@echo '  LDFLAGS      =' $(LDFLAGS)
	@echo '  APPONLYFLAGS =' $(APPONLYFLAGS)
	@echo '  MODEXTENSION =' $(MODEXTENSION)
	@echo '  MOVEMOD      =' $(MOVEMOD)
	@echo '  INC          =' $(INC)
	@echo '  INTSIZE      =' $(INTSIZE)
	@echo '  SRCEXTRA     =' $(SRCEXTRA)
	@echo '  BUILDLIB     =' $(BUILDLIB)
	@echo '  LIBDIR       =' $(LIBDIR)
	@echo '  LIBS         =' $(LIBS)
        # If preprocessor switch -DENABLE_SERIAL_BUILD does occur in compiler flags,
        # a build for serial execution is requested.
	@if test -z "$(findstring -DENABLE_SERIAL_BUILD ,$(APPONLYFLAGS) $(CFLAGSF90) )"; then \
	    echo '  MPIINC       =' $(MPIINC); \
	    echo '  MPILIBDIR    =' $(MPILIBDIR); \
	    echo '  MPILIBS      =' $(MPILIBS); \
	fi
	@echo '  AR           =' $(AR)
	@echo '  RANLIB       =' $(RANLIB)
	@echo '  OBJDIR       =' $(OBJDIR)
	@echo '  OBJDIR_LIB   =' $(OBJDIR_LIB)
ifneq ($(strip $(COPROCESSOR)),)
	@echo '  COPROCESSOR  =' $(COPROCESSOR)
endif
	@echo
        # If preprocessor switch -DENABLE_SERIAL_BUILD does occur in compiler flags,
        # a build for serial execution is requested.
	@if test -z "$(findstring -DENABLE_SERIAL_BUILD ,$(APPONLYFLAGS) $(CFLAGSF90) )"; then \
	    echo 'Build with MPI enabled.'; \
	else \
	    echo 'Serial build enabled.'; \
	fi
        # Catch cases where something requires a user's attention and
        # show the message we will get when trying to compile also here.
	@if test -n "$(MESSAGE)"; then \
	    echo; $(MESSAGE) \
	fi



# Note on options implied by composite optimization switches
# for GNU C/C++ compiler suite. This information is provided here
# in order to prevent having added unnecessary compile switches to
# CFLAGSC/CFLAGSCXX below.
#
# -O1 => -fcprop-registers -fcrossjumping -fdefer-pop -fdelayed-branch
#        -fdelayed-branch -fguess-branch-prob -fif-conversion
#        -fif-conversion2 -floop-optimize -fmerge-constants
#        -fomit-frame-pointer -fthread-jumps
#
# -O2 => -O1
#        -falign-functions -falign-jumps -falign-labels -falign-loops
#        -fcaller-saves -fcse-follow-jumps -fcse-skip-blocks
#        -fdelete-null-pointer-checks -fexpensive-optimizations -fforce-mem
#        -fgcse -foptimize-sibling-calls -fpeephole2 -fregmove -freorder-blocks
#        -frerun-cse-after-loop -frerun-loop-opt -fschedule-insns
#        -fschedule-insns-after-reload -fstrength-reduce -fstrict-aliasing
#        -funit-at-a-time
#
# -O3 => -O2
#        -finline-functions -frename-registers -fweb
#
# -march=pentium4 => -mmmx -msse -msse2
#                    -mieee-fp -malign-double -mno-push-args
#                    -maccumulate-outgoing-args -momit-leaf-frame-pointer
#                    -minline-all-stringops -mfpmath=
#

# Tell make to delete the target of a rule if it has changed and its commands
# exit with a nonzero exit status (just as it does when it receives a signal).
# By this rule, make will not create an empty .tex file when the java parser
# detects a syntax error in an input file. If we would allow such empty files
# being created, such an error would go undetected in subsequent invocations
# of make resulting in an incomplete documentation.
.DELETE_ON_ERROR:
