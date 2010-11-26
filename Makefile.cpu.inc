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
match=$(shell echo $(1) | awk '/$(2)/ { print "yes"; }')
#
# variant 2 (fastest so far, but does not work for rules like (pc|pc64)-)
#match=$(shell echo $(1) | sed -e 's/$(2)/yes/')

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
# C preprocessor flags
CPPFLAGS      = -x assembler-with-cpp \
		$(shell echo $(APPONLYFLAGS) $(CFLAGSF90) | tr ' ' '\n' | sed -e '/^-[DU]/!d')

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


# # # # # # # # #      IMPORTANT NOTICE !!!     # # # # # # # # #
#

# Writing consistent buildIDs for all combinations of systems,
# compilers and libraries is just not maintainable. Therefore, this
# file provides generic rules which should work for most systems. If
# thecompiler provides some auto-tune mechanism, then the generic
# rules even try to optimise the code for the underlying machine.
#
# There are some situation in which you want to add explicit buildIDs
#
# 1) if you have an exotic computing environment (e.g., NEC SX) which
#    need a detailed definition of compilers and libraries (e.g. with
#    absolute paths) which cannot be covered by generic rules
#
# 2) if you want to overwrite the default rule for a particular buildID;
#    keep in mind that some settings can be changed by using flags:
#
#    INTSIZE = [STD|LARGE] : forces the compiler to set the default
#    size of integer variables to 32|64 bit
#
#    BINARY = [32|64] : forces the compiler to build 32|64 bit binary
#    
# 3) some compilers do not provide an automised optimisation feature
#    (e.g., GCC 4.1.x and below, G95) so that you have to provide the
#    optimisation flags by hand
#
# 4) if you want to cross-compile your code
#
# # # # # # # # #      IMPORTANT NOTICE !!!     # # # # # # # # #



##############################################################################
# Alpha
#
##############################################################################

ifeq ($(call match,$(ID),alpha-ev6-osf1-cf90-.*),yes)
# invokes the Compaq Fortran 77/90 (formerly DIGITAL Fortran 77/90)
# compiler, local umfpack
# No MPI wrappers installed on our alphas
MPIWRAPPERS = NO
include $(TEMPLATESDIR)/alpha-ev6-osf1-generic.mk
include $(TEMPLATESDIR)/decf90-generic.mk
include $(TEMPLATESDIR)/osf1-generic.mk
endif

ifeq ($(call match,$(ID),alpha-ev6-osf1-.*-dxml.*),yes)
# Use Digital BLAS
TOKEN5 := 1
LIBS      := $(LIBS) -ldxml
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



##############################################################################
# AMD Athlon
#
##############################################################################

ifeq ($(call match,$(ID),pc-athlon-.*-g95-.*),yes)
# invokes the G95 Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=athlon $(CFLAGSF77LIBS)
CFLAGSF77     := -march=athlon $(CFLAGSF77)
CFLAGSF90     := -march=athlon $(CFLAGSF90)
CFLAGSC       := -march=athlon $(CFLAGSC)
CFLAGSCXX     := -march=athlon $(CFLAGSCXX)
LDFLAGS       := -march=athlon $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-athlon-.*-gcc-.*),yes)
# invokes the gfortran Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=athlon $(CFLAGSF77LIBS)
CFLAGSF77     := -march=athlon $(CFLAGSF77)
CFLAGSF90     := -march=athlon $(CFLAGSF90)
CFLAGSC       := -march=athlon $(CFLAGSC)
CFLAGSCXX     := -march=athlon $(CFLAGSCXX)
LDFLAGS       := -march=athlon $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-athlon-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE $(CFLAGSF77)
CFLAGSF90     := -xSSE $(CFLAGSF90)
CFLAGSC       := -xSSE $(CFLAGSC)
CFLAGSCXX     := -xSSE $(CFLAGSCXX)
LDFLAGS       := -xSSE $(LDFLAGS)
else
CFLAGSF77LIBS := -axK $(CFLAGSF77LIBS)
CFLAGSF77     := -axK $(CFLAGSF77)
CFLAGSF90     := -axK $(CFLAGSF90)
CFLAGSC       := -axK $(CFLAGSC)
CFLAGSCXX     := -axK $(CFLAGSCXX)
LDFLAGS       := -axK $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc-athlon-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp athlon $(CFLAGSF77LIBS)
CFLAGSF77     := -tp athlon $(CFLAGSF77)
CFLAGSF90     := -tp athlon $(CFLAGSF90)
CFLAGSC       := -tp athlon $(CFLAGSC)
CFLAGSCXX     := -tp athlon $(CFLAGSCXX)
LDFLAGS       := -tp athlon $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-athlon-.*-psc-.*),yes)
# invokes the PathScale Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pathscale-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=athlon $(CFLAGSF77LIBS)
CFLAGSF77     := -march=athlon $(CFLAGSF77)
CFLAGSF90     := -march=athlon $(CFLAGSF90)
CFLAGSC       := -march=athlon $(CFLAGSC)
CFLAGSCXX     := -march=athlon $(CFLAGSCXX)
LDFLAGS       := -march=athlon $(LDFLAGS)
endif
endif



##############################################################################
# AMD Athlon 64 (X2) / Turion 64 (X2)
#
##############################################################################

ifeq ($(call match,$(ID),(pc|pc64)-(athlon64|athlon64x2|turion64|turion64x2)-.*-g95-.*),yes)
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

ifeq ($(call match,$(ID),(pc|pc64)-(athlon64|athlon64x2|turion64|turion64x2)-.*-gcc-.*),yes)
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

ifeq ($(call match,$(ID),(pc|pc64)-athlon64-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE2 $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE2 $(CFLAGSF77)
CFLAGSF90     := -xSSE2 $(CFLAGSF90)
CFLAGSC       := -xSSE2 $(CFLAGSC)
CFLAGSCXX     := -xSSE2 $(CFLAGSCXX)
LDFLAGS       := -xSSE2 $(LDFLAGS)
else
CFLAGSF77LIBS := -axW $(CFLAGSF77LIBS)
CFLAGSF77     := -axW $(CFLAGSF77)
CFLAGSF90     := -axW $(CFLAGSF90)
CFLAGSC       := -axW $(CFLAGSC)
CFLAGSCXX     := -axW $(CFLAGSCXX)
LDFLAGS       := -axW $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-(athlon64x2|turion64|turion64x2)-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE3 $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE3 $(CFLAGSF77)
CFLAGSF90     := -xSSE3 $(CFLAGSF90)
CFLAGSC       := -xSSE3 $(CFLAGSC)
CFLAGSCXX     := -xSSE3 $(CFLAGSCXX)
LDFLAGS       := -xSSE3 $(LDFLAGS)
else
CFLAGSF77LIBS := -axP $(CFLAGSF77LIBS)
CFLAGSF77     := -axP $(CFLAGSF77)
CFLAGSF90     := -axP $(CFLAGSF90)
CFLAGSC       := -axP $(CFLAGSC)
CFLAGSCXX     := -axP $(CFLAGSCXX)
LDFLAGS       := -axP $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc-(athlon64|athlon64x2|turion64|turion64x2)-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp k8-32 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp k8-32 $(CFLAGSF77)
CFLAGSF90     := -tp k8-32 $(CFLAGSF90)
CFLAGSC       := -tp k8-32 $(CFLAGSC)
CFLAGSCXX     := -tp k8-32 $(CFLAGSCXX)
LDFLAGS       := -tp k8-32 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-(athlon64|athlon64x2|turion64|turion64x2)-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp k8-64 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp k8-64 $(CFLAGSF77)
CFLAGSF90     := -tp k8-64 $(CFLAGSF90)
CFLAGSC       := -tp k8-64 $(CFLAGSC)
CFLAGSCXX     := -tp k8-64 $(CFLAGSCXX)
LDFLAGS       := -tp k8-64 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-(athlon64|athlon64x2|turion64|turion64x2)-.*-psc-.*),yes)
# invokes the PathScale Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pathscale-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=athlon64 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=athlon64 $(CFLAGSF77)
CFLAGSF90     := -march=athlon64 $(CFLAGSF90)
CFLAGSC       := -march=athlon64 $(CFLAGSC)
CFLAGSCXX     := -march=athlon64 $(CFLAGSCXX)
LDFLAGS       := -march=athlon64 $(LDFLAGS)
endif
endif



##############################################################################
# AMD Athlon XP/MP
#
##############################################################################

ifeq ($(call match,$(ID),pc-athlonxp-.*-g95-.*),yes)
# invokes the G95 Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=athlon-xp $(CFLAGSF77LIBS)
CFLAGSF77     := -march=athlon-xp $(CFLAGSF77)
CFLAGSF90     := -march=athlon-xp $(CFLAGSF90)
CFLAGSC       := -march=athlon-xp $(CFLAGSC)
CFLAGSCXX     := -march=athlon-xp $(CFLAGSCXX)
LDFLAGS       := -march=athlon-xp $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-athlonxp-.*-gcc-.*),yes)
# invokes the gfortran Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=athlon-xp $(CFLAGSF77LIBS)
CFLAGSF77     := -march=athlon-xp $(CFLAGSF77)
CFLAGSF90     := -march=athlon-xp $(CFLAGSF90)
CFLAGSC       := -march=athlon-xp $(CFLAGSC)
CFLAGSCXX     := -march=athlon-xp $(CFLAGSCXX)
LDFLAGS       := -march=athlon-xp $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-athlonxp-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE $(CFLAGSF77)
CFLAGSF90     := -xSSE $(CFLAGSF90)
CFLAGSC       := -xSSE $(CFLAGSC)
CFLAGSCXX     := -xSSE $(CFLAGSCXX)
LDFLAGS       := -xSSE $(LDFLAGS)
else
CFLAGSF77LIBS := -axK $(CFLAGSF77LIBS)
CFLAGSF77     := -axK $(CFLAGSF77)
CFLAGSF90     := -axK $(CFLAGSF90)
CFLAGSC       := -axK $(CFLAGSC)
CFLAGSCXX     := -axK $(CFLAGSCXX)
LDFLAGS       := -axK $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc-athlonxp-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp athlonxp $(CFLAGSF77LIBS)
CFLAGSF77     := -tp athlonxp $(CFLAGSF77)
CFLAGSF90     := -tp athlonxp $(CFLAGSF90)
CFLAGSC       := -tp athlonxp $(CFLAGSC)
CFLAGSCXX     := -tp athlonxp $(CFLAGSCXX)
LDFLAGS       := -tp athlonxp $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-athlonxp-.*-psc-.*),yes)
# invokes the PathScale Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pathscale-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=athlon $(CFLAGSF77LIBS)
CFLAGSF77     := -march=athlon $(CFLAGSF77)
CFLAGSF90     := -march=athlon $(CFLAGSF90)
CFLAGSC       := -march=athlon $(CFLAGSC)
CFLAGSCXX     := -march=athlon $(CFLAGSCXX)
LDFLAGS       := -march=athlon $(LDFLAGS)
endif
endif



##############################################################################
# AMD Athlon II
#
##############################################################################

# MM: Don't know yet what type of processor this is



##############################################################################
# AMD Opteron (X2)
#
##############################################################################

ifeq ($(call match,$(ID),(pc|pc64)-(opteron|opteronx2)-.*-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=opteron $(CFLAGSF77LIBS)
CFLAGSF77     := -march=opteron $(CFLAGSF77)
CFLAGSF90     := -march=opteron $(CFLAGSF90)
CFLAGSC       := -march=opteron $(CFLAGSC)
CFLAGSCXX     := -march=opteron $(CFLAGSCXX)
LDFLAGS       := -march=opteron $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-(opteron|opteronx2)-.*-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=opteron $(CFLAGSF77LIBS)
CFLAGSF77     := -march=opteron $(CFLAGSF77)
CFLAGSF90     := -march=opteron $(CFLAGSF90)
CFLAGSC       := -march=opteron $(CFLAGSC)
CFLAGSCXX     := -march=opteron $(CFLAGSCXX)
LDFLAGS       := -march=opteron $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),(pc|pc64)-(opteron|opteronx2)-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE2 $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE2 $(CFLAGSF77)
CFLAGSF90     := -xSSE2 $(CFLAGSF90)
CFLAGSC       := -xSSE2 $(CFLAGSC)
CFLAGSCXX     := -xSSE2 $(CFLAGSCXX)
LDFLAGS       := -xSSE2 $(LDFLAGS)
else
CFLAGSF77LIBS := -axW $(CFLAGSF77LIBS)
CFLAGSF77     := -axW $(CFLAGSF77)
CFLAGSF90     := -axW $(CFLAGSF90)
CFLAGSC       := -axW $(CFLAGSC)
CFLAGSCXX     := -axW $(CFLAGSCXX)
LDFLAGS       := -axW $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc-phenom(opteron|opteronx2)-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp k8-32 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp k8-32 $(CFLAGSF77)
CFLAGSF90     := -tp k8-32 $(CFLAGSF90)
CFLAGSC       := -tp k8-32 $(CFLAGSC)
CFLAGSCXX     := -tp k8-32 $(CFLAGSCXX)
LDFLAGS       := -tp k8-32 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-(opteron|opteronx2)-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp k8-64 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp k8-64 $(CFLAGSF77)
CFLAGSF90     := -tp k8-64 $(CFLAGSF90)
CFLAGSC       := -tp k8-64 $(CFLAGSC)
CFLAGSCXX     := -tp k8-64 $(CFLAGSCXX)
LDFLAGS       := -tp k8-64 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-(opteron|opteronx2)-.*-psc-.*),yes)
# invokes the PathScale Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pathscale-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -mtune=opteron -mcpu=opteron -march=opteron $(CFLAGSF77LIBS)
CFLAGSF77     := -mtune=opteron -mcpu=opteron -march=opteron $(CFLAGSF77)
CFLAGSF90     := -mtune=opteron -mcpu=opteron -march=opteron $(CFLAGSF90)
CFLAGSC       := -mtune=opteron -mcpu=opteron -march=opteron $(CFLAGSC)
CFLAGSCCC     := -mtune=opteron -mcpu=opteron -march=opteron $(CFLAGSCCC)
LDFLAGS       := -mtune=opteron -mcpu=opteron -march=opteron $(LDFLAGS)
endif
endif


###
### SPECIAL RULES NOT CHECKED YET !!!
###

ifeq ($(call match,$(ID),pc64-opteronx2-linux-.*-goto2),yes)
# Use GotoBLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
include $(TEMPLATESDIR)/linux-generic.mk
# Workaround for linker problem:
# ld: hidden symbol `__svml_cosf4' in libsvml.a(svml_stub_scos4.o) is referenced by DSO
ifeq ($(call match,$(ID),pc64-opteronx2-linux-intel-goto2.*),yes)
LIBS := $(LIBS) -lsvml
endif
endif

ifeq ($(call match,$(ID),pc64-opteron-linux-.*-.*-optmpich),yes)
# Use Allinea Opt MPI
include $(TEMPLATESDIR)/optmpich-generic.mk
include $(TEMPLATESDIR)/linux-generic.mk
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



##############################################################################
# AMD Phenom (Barcelona)
#
##############################################################################

ifeq ($(call match,$(ID),(pc|pc64)-phenom-.*-g95-.*),yes)
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

ifeq ($(call match,$(ID),(pc|pc64)-phenom-.*-gcc-.*),yes)
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

ifeq ($(call match,$(ID),(pc|pc64)-phenom-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE4.1 $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE4.1 $(CFLAGSF77)
CFLAGSF90     := -xSSE4.1 $(CFLAGSF90)
CFLAGSC       := -xSSE4.1 $(CFLAGSC)
CFLAGSCXX     := -xSSE4.1 $(CFLAGSCXX)
LDFLAGS       := -xSSE4.1 $(LDFLAGS)
else
CFLAGSF77LIBS := -xT $(CFLAGSF77LIBS)
CFLAGSF77     := -xT $(CFLAGSF77)
CFLAGSF90     := -xT $(CFLAGSF90)
CFLAGSC       := -xT $(CFLAGSC)
CFLAGSCXX     := -xT $(CFLAGSCXX)
LDFLAGS       := -xT $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc-phenom-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp barcelona-32 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp barcelona-32 $(CFLAGSF77)
CFLAGSF90     := -tp barcelona-32 $(CFLAGSF90)
CFLAGSC       := -tp barcelona-32 $(CFLAGSC)
CFLAGSCXX     := -tp barcelona-32 $(CFLAGSCXX)
LDFLAGS       := -tp barcelona-32 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-phenom-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp barcelona-64 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp barcelona-64 $(CFLAGSF77)
CFLAGSF90     := -tp barcelona-64 $(CFLAGSF90)
CFLAGSC       := -tp barcelona-64 $(CFLAGSC)
CFLAGSCXX     := -tp barcelona-64 $(CFLAGSCXX)
LDFLAGS       := -tp barcelona-64 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-phenom-.*-psc-.*),yes)
# invokes the PathScale Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pathscale-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=barcelona $(CFLAGSF77LIBS)
CFLAGSF77     := -march=barcelona $(CFLAGSF77)
CFLAGSF90     := -march=barcelona $(CFLAGSF90)
CFLAGSC       := -march=barcelona $(CFLAGSC)
CFLAGSCXX     := -march=barcelona $(CFLAGSCXX)
LDFLAGS       := -march=barcelona $(LDFLAGS)
endif
endif



##############################################################################
# AMD Phenom II (Shanghai)
#
##############################################################################

ifeq ($(call match,$(ID),(pc|pc64)-phenomII-.*-g95-.*),yes)
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

ifeq ($(call match,$(ID),(pc|pc64)-phenomII-.*-gcc-.*),yes)
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

ifeq ($(call match,$(ID),(pc|pc64)-phenomII-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE4.1 $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE4.1 $(CFLAGSF77)
CFLAGSF90     := -xSSE4.1 $(CFLAGSF90)
CFLAGSC       := -xSSE4.1 $(CFLAGSC)
CFLAGSCXX     := -xSSE4.1 $(CFLAGSCXX)
LDFLAGS       := -xSSE4.1 $(LDFLAGS)
else
CFLAGSF77LIBS := -xT $(CFLAGSF77LIBS)
CFLAGSF77     := -xT $(CFLAGSF77)
CFLAGSF90     := -xT $(CFLAGSF90)
CFLAGSC       := -xT $(CFLAGSC)
CFLAGSCXX     := -xT $(CFLAGSCXX)
LDFLAGS       := -xT $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc-phenomII-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp shanghai-32 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp shanghai-32 $(CFLAGSF77)
CFLAGSF90     := -tp shanghai-32 $(CFLAGSF90)
CFLAGSC       := -tp shanghai-32 $(CFLAGSC)
CFLAGSCXX     := -tp shanghai-32 $(CFLAGSCXX)
LDFLAGS       := -tp shanghai-32 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-phenomII-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp shanghai-64 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp shanghai-64 $(CFLAGSF77)
CFLAGSF90     := -tp shanghai-64 $(CFLAGSF90)
CFLAGSC       := -tp shanghai-64 $(CFLAGSC)
CFLAGSCXX     := -tp shanghai-64 $(CFLAGSCXX)
LDFLAGS       := -tp shanghai-64 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-phenomII-.*-psc-.*),yes)
# invokes the PathScale Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pathscale-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=barcelona $(CFLAGSF77LIBS)
CFLAGSF77     := -march=barcelona $(CFLAGSF77)
CFLAGSF90     := -march=barcelona $(CFLAGSF90)
CFLAGSC       := -march=barcelona $(CFLAGSC)
CFLAGSCXX     := -march=barcelona $(CFLAGSCXX)
LDFLAGS       := -march=barcelona $(LDFLAGS)
endif
endif



##############################################################################
# AMD Opteron 4100 (Istanbul) / 6100 (Magny Cours)
#
##############################################################################

ifeq ($(call match,$(ID),(pc|pc64)-(opteron4100|opteron6100)-.*-g95-.*),yes)
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

ifeq ($(call match,$(ID),(pc|pc64)-(opteron4100|opteron6100)-.*-gcc-.*),yes)
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

ifeq ($(call match,$(ID),(pc|pc64)-(opteron4100|opteron6100)-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE4.1 $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE4.1 $(CFLAGSF77)
CFLAGSF90     := -xSSE4.1 $(CFLAGSF90)
CFLAGSC       := -xSSE4.1 $(CFLAGSC)
CFLAGSCXX     := -xSSE4.1 $(CFLAGSCXX)
LDFLAGS       := -xSSE4.1 $(LDFLAGS)
else
CFLAGSF77LIBS := -xT $(CFLAGSF77LIBS)
CFLAGSF77     := -xT $(CFLAGSF77)
CFLAGSF90     := -xT $(CFLAGSF90)
CFLAGSC       := -xT $(CFLAGSC)
CFLAGSCXX     := -xT $(CFLAGSCXX)
LDFLAGS       := -xT $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc-(opteron4100|opteron6100)-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp istanbul-32 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp istanbul-32 $(CFLAGSF77)
CFLAGSF90     := -tp istanbul-32 $(CFLAGSF90)
CFLAGSC       := -tp istanbul-32 $(CFLAGSC)
CFLAGSCXX     := -tp istanbul-32 $(CFLAGSCXX)
LDFLAGS       := -tp istanbul-32 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-(opteron4100|opteron6100)-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp istanbul-64 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp istanbul-64 $(CFLAGSF77)
CFLAGSF90     := -tp istanbul-64 $(CFLAGSF90)
CFLAGSC       := -tp istanbul-64 $(CFLAGSC)
CFLAGSCXX     := -tp istanbul-64 $(CFLAGSCXX)
LDFLAGS       := -tp istanbul-64 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-phenomII-.*-psc-.*),yes)
# invokes the PathScale Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pathscale-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=barcelona $(CFLAGSF77LIBS)
CFLAGSF77     := -march=barcelona $(CFLAGSF77)
CFLAGSF90     := -march=barcelona $(CFLAGSF90)
CFLAGSC       := -march=barcelona $(CFLAGSC)
CFLAGSCXX     := -march=barcelona $(CFLAGSCXX)
LDFLAGS       := -march=barcelona $(LDFLAGS)
endif
endif



##############################################################################
# IBM Regatta
#
##############################################################################

ifeq ($(call match,$(ID),ibm-powerpc_power4-aix-xlf-.*-poempi),yes)
# invokes the IBM XL Fortran compiler, uses POE MPI
# local umfpack
TOKEN6 := 1
include $(TEMPLATESDIR)/aix-generic.mk
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



##############################################################################
# Intel Pentium 3
#
##############################################################################

ifeq ($(call match,$(ID),pc-pentium3-.*-g95-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium3 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium3 $(CFLAGSF77)
CFLAGSF90     := -march=pentium3 $(CFLAGSF90)
CFLAGSC       := -march=pentium3 $(CFLAGSC)
CFLAGSCXX     := -march=pentium3 $(CFLAGSCXX)
LDFLAGS       := -march=pentium3 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-pentium3-.*-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium3 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium3 $(CFLAGSF77)
CFLAGSF90     := -march=pentium3 $(CFLAGSF90)
CFLAGSC       := -march=pentium3 $(CFLAGSC)
CFLAGSCXX     := -march=pentium3 $(CFLAGSCXX)
LDFLAGS       := -march=pentium3 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-pentium3-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE $(CFLAGSF77)
CFLAGSF90     := -xSSE $(CFLAGSF90)
CFLAGSC       := -xSSE $(CFLAGSC)
CFLAGSCXX     := -xSSE $(CFLAGSCXX)
LDFLAGS       := -xSSE $(LDFLAGS)
else
CFLAGSF77LIBS := -xK $(CFLAGSF77LIBS)
CFLAGSF77     := -xK $(CFLAGSF77)
CFLAGSF90     := -xK $(CFLAGSF90)
CFLAGSC       := -xK $(CFLAGSC)
CFLAGSCXX     := -xK $(CFLAGSCXX)
LDFLAGS       := -xK $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc-pentium3-.*-pgi-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp piii $(CFLAGSF77LIBS)
CFLAGSF77     := -tp piii $(CFLAGSF77)
CFLAGSF90     := -tp piii $(CFLAGSF90)
CFLAGSC       := -tp piii $(CFLAGSC)
CFLAGSCXX     := -tp piii $(CFLAGSCXX)
LDFLAGS       := -tp piii $(LDFLAGS)
endif
endif



##############################################################################
# Intel Pentium M
#
##############################################################################

ifeq ($(call match,$(ID),pc-pentiumm-.*-g95-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium-m $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium-m $(CFLAGSF77)
CFLAGSF90     := -march=pentium-m $(CFLAGSF90)
CFLAGSC       := -march=pentium-m $(CFLAGSC)
CFLAGSCXX     := -march=pentium-m $(CFLAGSCXX)
LDFLAGS       := -march=pentium-m $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-pentiumm-.*-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium-m $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium-m $(CFLAGSF77)
CFLAGSF90     := -march=pentium-m $(CFLAGSF90)
CFLAGSC       := -march=pentium-m $(CFLAGSC)
CFLAGSCXX     := -march=pentium-m $(CFLAGSCXX)
LDFLAGS       := -march=pentium-m $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-pentiumm-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE2 $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE2 $(CFLAGSF77)
CFLAGSF90     := -xSSE2 $(CFLAGSF90)
CFLAGSC       := -xSSE2 $(CFLAGSC)
CFLAGSCXX     := -xSSE2 $(CFLAGSCXX)
LDFLAGS       := -xSSE2 $(LDFLAGS)
else
CFLAGSF77LIBS := -xB $(CFLAGSF77LIBS)
CFLAGSF77     := -xB $(CFLAGSF77)
CFLAGSF90     := -xB $(CFLAGSF90)
CFLAGSC       := -xB $(CFLAGSC)
CFLAGSCXX     := -xB $(CFLAGSCXX)
LDFLAGS       := -xB $(LDFLAGS)
endif
endif
endif



##############################################################################
# Intel Pentium 4
#
##############################################################################

ifeq ($(call match,$(ID),pc-pentium4-.*-g95-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium4 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium4 $(CFLAGSF77)
CFLAGSF90     := -march=pentium4 $(CFLAGSF90)
CFLAGSC       := -march=pentium4 $(CFLAGSC)
CFLAGSCXX     := -march=pentium4 $(CFLAGSCXX)
LDFLAGS       := -march=pentium4 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-pentium4-.*-g95-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=nocona $(CFLAGSF77LIBS)
CFLAGSF77     := -march=nocona $(CFLAGSF77)
CFLAGSF90     := -march=nocona $(CFLAGSF90)
CFLAGSC       := -march=nocona $(CFLAGSC)
CFLAGSCXX     := -march=nocona $(CFLAGSCXX)
LDFLAGS       := -march=nocona $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-pentium4-.*-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium4 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium4 $(CFLAGSF77)
CFLAGSF90     := -march=pentium4 $(CFLAGSF90)
CFLAGSC       := -march=pentium4 $(CFLAGSC)
CFLAGSCXX     := -march=pentium4 $(CFLAGSCXX)
LDFLAGS       := -march=pentium4 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-pentium4-.*-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=nocona $(CFLAGSF77LIBS)
CFLAGSF77     := -march=nocona $(CFLAGSF77)
CFLAGSF90     := -march=nocona $(CFLAGSF90)
CFLAGSC       := -march=nocona $(CFLAGSC)
CFLAGSCXX     := -march=nocona $(CFLAGSCXX)
LDFLAGS       := -march=nocona $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-pentium4-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE2 $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE2 $(CFLAGSF77)
CFLAGSF90     := -xSSE2 $(CFLAGSF90)
CFLAGSC       := -xSSE2 $(CFLAGSC)
CFLAGSCXX     := -xSSE2 $(CFLAGSCXX)
LDFLAGS       := -xSSE2 $(LDFLAGS)
else
CFLAGSF77LIBS := -xN $(CFLAGSF77LIBS)
CFLAGSF77     := -xN $(CFLAGSF77)
CFLAGSF90     := -xN $(CFLAGSF90)
CFLAGSC       := -xN $(CFLAGSC)
CFLAGSCXX     := -xN $(CFLAGSCXX)
LDFLAGS       := -xN $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc64-pentium4-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE2 $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE2 $(CFLAGSF77)
CFLAGSF90     := -xSSE2 $(CFLAGSF90)
CFLAGSC       := -xSSE2 $(CFLAGSC)
CFLAGSCXX     := -xSSE2 $(CFLAGSCXX)
LDFLAGS       := -xSSE2 $(LDFLAGS)
else
CFLAGSF77LIBS := -xW $(CFLAGSF77LIBS)
CFLAGSF77     := -xW $(CFLAGSF77)
CFLAGSF90     := -xW $(CFLAGSF90)
CFLAGSC       := -xW $(CFLAGSC)
CFLAGSCXX     := -xW $(CFLAGSCXX)
LDFLAGS       := -xW $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-pentium4-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp piv $(CFLAGSF77LIBS)
CFLAGSF77     := -tp piv $(CFLAGSF77)
CFLAGSF90     := -tp piv $(CFLAGSF90)
CFLAGSC       := -tp piv $(CFLAGSC)
CFLAGSCXX     := -tp piv $(CFLAGSCXX)
LDFLAGS       := -tp piv $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-pentium4-.*-psc-.*),yes)
# invokes the PathScale Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pathscale-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=pentium4 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=pentium4 $(CFLAGSF77)
CFLAGSF90     := -march=pentium4 $(CFLAGSF90)
CFLAGSC       := -march=pentium4 $(CFLAGSC)
CFLAGSCXX     := -march=pentium4 $(CFLAGSCXX)
LDFLAGS       := -march=pentium4 $(LDFLAGS)
endif
endif



##############################################################################
# Intel CoreDuo 2
#
##############################################################################

ifeq ($(call match,$(ID),pc-coreduo-.*-g95-.*),yes)
# invokes the GNU Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=prescott $(CFLAGSF77LIBS)
CFLAGSF77     := -march=prescott $(CFLAGSF77)
CFLAGSF90     := -march=prescott $(CFLAGSF90)
CFLAGSC       := -march=prescott $(CFLAGSC)
CFLAGSCXX     := -march=prescott $(CFLAGSCXX)
LDFLAGS       := -march=prescott $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-coreduo-.*-g95-.*),yes)
# invokes the GNU Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=nocona $(CFLAGSF77LIBS)
CFLAGSF77     := -march=nocona $(CFLAGSF77)
CFLAGSF90     := -march=nocona $(CFLAGSF90)
CFLAGSC       := -march=nocona $(CFLAGSC)
CFLAGSCXX     := -march=nocona $(CFLAGSCXX)
LDFLAGS       := -march=nocona $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc-coreduo-.*-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=prescott -msse3 -mfpmath=sse $(CFLAGSF77LIBS)
CFLAGSF77     := -march=prescott -msse3 -mfpmath=sse $(CFLAGSF77)
CFLAGSF90     := -march=prescott -msse3 -mfpmath=sse $(CFLAGSF90)
CFLAGSC       := -march=prescott -msse3 -mfpmath=sse $(CFLAGSC)
CFLAGSCXX     := -march=prescott -msse3 -mfpmath=sse $(CFLAGSCXX)
LDFLAGS       := -march=prescott -msse3 -mfpmath=sse $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-coreduo-.*-gcc-.*),yes)
# invokes the GNU Fortran compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call gfortranminversion_4_3),yes)
CFLAGSF77LIBS := -march=core2 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=core2 $(CFLAGSF77)
CFLAGSF90     := -march=core2 $(CFLAGSF90)
CFLAGSC       := -march=core2 $(CFLAGSC)
CFLAGSCXX     := -march=core2 $(CFLAGSCXX)
LDFLAGS       := -march=core2 $(LDFLAGS)
else
CFLAGSF77LIBS := -march=nocona -msse3 -mfpmath=sse $(CFLAGSF77LIBS)
CFLAGSF77     := -march=nocona -msse3 -mfpmath=sse $(CFLAGSF77)
CFLAGSF90     := -march=nocona -msse3 -mfpmath=sse $(CFLAGSF90)
CFLAGSC       := -march=nocona -msse3 -mfpmath=sse $(CFLAGSC)
CFLAGSCXX     := -march=nocona -msse3 -mfpmath=sse $(CFLAGSCXX)
LDFLAGS       := -march=nocona -msse3 -mfpmath=sse $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc-coreduo-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE3 $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE3 $(CFLAGSF77)
CFLAGSF90     := -xSSE3 $(CFLAGSF90)
CFLAGSC       := -xSSE3 $(CFLAGSC)
CFLAGSCXX     := -xSSE3 $(CFLAGSCXX)
LDFLAGS       := -xSSE3 $(LDFLAGS)
else
CFLAGSF77LIBS := -xP $(CFLAGSF77LIBS)
CFLAGSF77     := -xP $(CFLAGSF77)
CFLAGSF90     := -xP $(CFLAGSF90)
CFLAGSC       := -xP $(CFLAGSC)
CFLAGSCXX     := -xP $(CFLAGSCXX)
LDFLAGS       := -xP $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc64-coreduo-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE3 $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE3 $(CFLAGSF77)
CFLAGSF90     := -xSSE3 $(CFLAGSF90)
CFLAGSC       := -xSSE3 $(CFLAGSC)
CFLAGSCXX     := -xSSE3 $(CFLAGSCXX)
LDFLAGS       := -xSSE3 $(LDFLAGS)
else
CFLAGSF77LIBS := -xT $(CFLAGSF77LIBS)
CFLAGSF77     := -xT $(CFLAGSF77)
CFLAGSF90     := -xT $(CFLAGSF90)
CFLAGSC       := -xT $(CFLAGSC)
CFLAGSCXX     := -xT $(CFLAGSCXX)
LDFLAGS       := -xT $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc-coreduo-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp core2-32 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp core2-32 $(CFLAGSF77)
CFLAGSF90     := -tp core2-32 $(CFLAGSF90)
CFLAGSC       := -tp core2-32 $(CFLAGSC)
CFLAGSCXX     := -tp core2-32 $(CFLAGSCXX)
LDFLAGS       := -tp core2-32 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-coreduo-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp core2-64 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp core2-64 $(CFLAGSF77)
CFLAGSF90     := -tp core2-64 $(CFLAGSF90)
CFLAGSC       := -tp core2-64 $(CFLAGSC)
CFLAGSCXX     := -tp core2-64 $(CFLAGSCXX)
LDFLAGS       := -tp core2-64 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-coreduo-.*-psc-.*),yes)
# invokes the PathScale Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pathscale-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=core $(CFLAGSF77LIBS)
CFLAGSF77     := -march=core $(CFLAGSF77)
CFLAGSF90     := -march=core $(CFLAGSF90)
CFLAGSC       := -march=core $(CFLAGSC)
CFLAGSCXX     := -march=core $(CFLAGSCXX)
LDFLAGS       := -march=core $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-coreduo-linux-.*-woodymkl.*),yes)
# Use special BKL blas on Woodcrest cluster of the RRZE Erlangen
include $(TEMPLATESDIR)/linux-generic.mk
LIBS := $(LIBS) $(MKL_LIB)
TOKEN5 := 1
endif



##############################################################################
# Intel CoreSolo
#
##############################################################################

ifeq ($(call match,$(ID),pc-coresolo-.*-g95-.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=prescott -mtune=prescott -msse $(CFLAGSF77LIBS)
CFLAGSF77     := -march=prescott -mtune=prescott -msse $(CFLAGSF77)
CFLAGSF90     := -march=prescott -mtune=prescott -msse $(CFLAGSF90)
CFLAGSC       := -march=prescott -mtune=prescott -msse $(CFLAGSC)
CFLAGSCXX     := -march=prescott -mtune=prescott -msse $(CFLAGSCXX)
endif
endif

ifeq ($(call match,$(ID),pc-coresolo-.*-gcc-.*),yes)
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

ifeq ($(call match,$(ID),pc-coresolo-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -msse3 $(CFLAGSF77LIBS)
CFLAGSF77     := -msse3 $(CFLAGSF77)
CFLAGSF90     := -msse3 $(CFLAGSF90)
CFLAGSC       := -msse3 $(CFLAGSC)
CFLAGSCXX     := -msse3 $(CFLAGSCXX)
LDFLAGS       := -msse3 $(LDFLAGS)
else
CFLAGSF77LIBS := -xO $(CFLAGSF77LIBS)
CFLAGSF77     := -xO $(CFLAGSF77)
CFLAGSF90     := -xO $(CFLAGSF90)
CFLAGSC       := -xO $(CFLAGSC)
CFLAGSCXX     := -xO $(CFLAGSCXX)
LDFLAGS       := -xO $(LDFLAGS)
endif
endif
endif



##############################################################################
# Intel Penryn
#
##############################################################################

ifeq ($(call match,$(ID),pc-penryn-.*-g95-.*),yes)
# invokes the GNU Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=prescott $(CFLAGSF77LIBS)
CFLAGSF77     := -march=prescott $(CFLAGSF77)
CFLAGSF90     := -march=prescott $(CFLAGSF90)
CFLAGSC       := -march=prescott $(CFLAGSC)
CFLAGSCXX     := -march=prescott $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-penryn-.*-g95-.*),yes)
# invokes the GNU Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=nocona $(CFLAGSF77LIBS)
CFLAGSF77     := -march=nocona $(CFLAGSF77)
CFLAGSF90     := -march=nocona $(CFLAGSF90)
CFLAGSC       := -march=nocona $(CFLAGSC)
CFLAGSCXX     := -march=nocona $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-penryn-.*-gcc-.*),yes)
# invokes the GNU Fortran 77/90 compiler, local umfpack
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

ifeq ($(call match,$(ID),pc64-penryn-.*-gcc-.*),yes)
# invokes the GNU Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call gfortranminversion_4_3),yes)
CFLAGSF77LIBS := -march=core2 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=core2 $(CFLAGSF77)
CFLAGSF90     := -march=core2 $(CFLAGSF90)
CFLAGSC       := -march=core2 $(CFLAGSC)
CFLAGSCXX     := -march=core2 $(CFLAGSCXX)
LDFLAGS       := -march=core2 $(LDFLAGS)
else
CFLAGSF77LIBS := -march=nocona $(CFLAGSF77LIBS)
CFLAGSF77     := -march=nocona $(CFLAGSF77)
CFLAGSF90     := -march=nocona $(CFLAGSF90)
CFLAGSC       := -march=nocona $(CFLAGSC)
CFLAGSCXX     := -march=nocona $(CFLAGSCXX)
LDFLAGS       := -march=nocona $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-penryn-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE4.1 $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE4.1 $(CFLAGSF77)
CFLAGSF90     := -xSSE4.1 $(CFLAGSF90)
CFLAGSC       := -xSSE4.1 $(CFLAGSC)
CFLAGSCXX     := -xSSE4.1 $(CFLAGSCXX)
LDFLAGS       := -xSSE4.1 $(LDFLAGS)
else
CFLAGSF77LIBS := -xT $(CFLAGSF77LIBS)
CFLAGSF77     := -xT $(CFLAGSF77)
CFLAGSF90     := -xT $(CFLAGSF90)
CFLAGSC       := -xT $(CFLAGSC)
CFLAGSCXX     := -xT $(CFLAGSCXX)
LDFLAGS       := -xT $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc-penryn-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp penryn-32 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp penryn-32 $(CFLAGSF77)
CFLAGSF90     := -tp penryn-32 $(CFLAGSF90)
CFLAGSC       := -tp penryn-32 $(CFLAGSC)
CFLAGSCXX     := -tp penryn-32 $(CFLAGSCXX)
LDFLAGS       := -tp penryn-32 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-penryn-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp penryn-64 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp penryn-64 $(CFLAGSF77)
CFLAGSF90     := -tp penryn-64 $(CFLAGSF90)
CFLAGSC       := -tp penryn-64 $(CFLAGSC)
CFLAGSCXX     := -tp penryn-64 $(CFLAGSCXX)
LDFLAGS       := -tp penryn-64 $(LDFLAGS)
endif
endif



##############################################################################
# Intel Nehalem
#
##############################################################################

ifeq ($(call match,$(ID),pc-nehalem-.*-g95-.*),yes)
# invokes the GNU Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=prescott $(CFLAGSF77LIBS)
CFLAGSF77     := -march=prescott $(CFLAGSF77)
CFLAGSF90     := -march=prescott $(CFLAGSF90)
CFLAGSC       := -march=prescott $(CFLAGSC)
CFLAGSCXX     := -march=prescott $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-nehalem-.*-g95-.*),yes)
# invokes the GNU Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=nocona $(CFLAGSF77LIBS)
CFLAGSF77     := -march=nocona $(CFLAGSF77)
CFLAGSF90     := -march=nocona $(CFLAGSF90)
CFLAGSC       := -march=nocona $(CFLAGSC)
CFLAGSCXX     := -march=nocona $(CFLAGSCXX)
LDFLAGS       := $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc-nehalem-.*-gcc-.*),yes)
# invokes the GNU Fortran 77/90 compiler, local umfpack
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

ifeq ($(call match,$(ID),pc64-nehalem-.*-gcc-.*),yes)
# invokes the GNU Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call gfortranminversion_4_3),yes)
CFLAGSF77LIBS := -march=core2 -mtune=core2 -msse4 -msse4.1 -msse4.2 -mfpmath=sse $(CFLAGSF77LIBS)
CFLAGSF77     := -march=core2 -mtune=core2 -msse4 -msse4.1 -msse4.2 -mfpmath=sse $(CFLAGSF77)
CFLAGSF90     := -march=core2 -mtune=core2 -msse4 -msse4.1 -msse4.2 -mfpmath=sse $(CFLAGSF90)
CFLAGSC       := -march=core2 -mtune=core2 -msse4 -msse4.1 -msse4.2 -mfpmath=sse $(CFLAGSC)
CFLAGSCXX     := -march=core2 -mtune=core2 -msse4 -msse4.1 -msse4.2 -mfpmath=sse $(CFLAGSCXX)
LDFLAGS       := -march=core2 -mtune=core2 -msse4 -msse4.1 -msse4.2 -mfpmath=sse $(LDFLAGS)
else
CFLAGSF77LIBS := -march=nocona $(CFLAGSF77LIBS)
CFLAGSF77     := -march=nocona $(CFLAGSF77)
CFLAGSF90     := -march=nocona $(CFLAGSF90)
CFLAGSC       := -march=nocona $(CFLAGSC)
CFLAGSCXX     := -march=nocona $(CFLAGSCXX)
LDFLAGS       := -march=nocona $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-nehalem-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xSSE4.2 $(CFLAGSF77LIBS)
CFLAGSF77     := -xSSE4.2 $(CFLAGSF77)
CFLAGSF90     := -xSSE4.2 $(CFLAGSF90)
CFLAGSC       := -xSSE4.2 $(CFLAGSC)
CFLAGSCXX     := -xSSE4.2 $(CFLAGSCXX)
LDFLAGS       := -xSSE4.2 $(LDFLAGS)
else
CFLAGSF77LIBS := -xT $(CFLAGSF77LIBS)
CFLAGSF77     := -xT $(CFLAGSF77)
CFLAGSF90     := -xT $(CFLAGSF90)
CFLAGSC       := -xT $(CFLAGSC)
CFLAGSCXX     := -xT $(CFLAGSCXX)
LDFLAGS       := -xT $(LDFLAGS)
endif
endif
endif

ifeq ($(call match,$(ID),pc-nehalem-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp nehalem-32 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp nehalem-32 $(CFLAGSF77)
CFLAGSF90     := -tp nehalem-32 $(CFLAGSF90)
CFLAGSC       := -tp nehalem-32 $(CFLAGSC)
CFLAGSCXX     := -tp nehalem-32 $(CFLAGSCXX)
LDFLAGS       := -tp nehalem-32 $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),pc64-nehalem-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -tp nehalem-64 $(CFLAGSF77LIBS)
CFLAGSF77     := -tp nehalem-64 $(CFLAGSF77)
CFLAGSF90     := -tp nehalem-64 $(CFLAGSF90)
CFLAGSC       := -tp nehalem-64 $(CFLAGSC)
CFLAGSCXX     := -tp nehalem-64 $(CFLAGSCXX)
LDFLAGS       := -tp nehalem-64 $(LDFLAGS)
endif
endif



##############################################################################
# Itanium
#
##############################################################################

ifeq ($(call match,$(ID),ia64-itanium-.*-gcc-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=itanium $(CFLAGSF77LIBS)
CFLAGSF77     := -march=itanium $(CFLAGSF77)
CFLAGSF90     := -march=itanium $(CFLAGSF90)
CFLAGSC       := -march=itanium $(CFLAGSC)
CFLAGSCXX     := -march=itanium $(CFLAGSCXX)
LDFLAGS       := -march=itanium $(LDFLAGS)
endif
#CFLAGSF77:=$(CFLAGSF77) -fno-underscoring
#CFLAGSF90:=$(CFLAGSF90) -fno-underscoring
endif



##############################################################################
# Itanium 2
#
##############################################################################

ifeq ($(call match,$(ID),ia64-itanium2-.*-gcc-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=itanium2 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=itanium2 $(CFLAGSF77)
CFLAGSF90     := -march=itanium2 $(CFLAGSF90)
CFLAGSC       := -march=itanium2 $(CFLAGSC)
CFLAGSCXX     := -march=itanium2 $(CFLAGSCXX)
LDFLAGS       := -march=itanium2 $(LDFLAGS)
endif
#CFLAGSF77:=$(CFLAGSF77) -fno-underscoring
#CFLAGSF90:=$(CFLAGSF90) -fno-underscoring
endif

ifeq ($(call match,$(ID),ia64-itanium2-linux-intel-blas-mpi),yes)
# invokes the intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/linux-generic.mk
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
include $(TEMPLATESDIR)/linux-generic.mk
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



##############################################################################
# Itanium 2 Dual Core
#
##############################################################################

ifeq ($(call match,$(ID),ia64-itanium2x2-.*-gcc-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
CFLAGSF77LIBS := -march=itanium2 $(CFLAGSF77LIBS)
CFLAGSF77     := -march=itanium2 $(CFLAGSF77)
CFLAGSF90     := -march=itanium2 $(CFLAGSF90)
CFLAGSC       := -march=itanium2 $(CFLAGSC)
CFLAGSCXX     := -march=itanium2 $(CFLAGSCXX)
LDFLAGS       := -march=itanium2 $(LDFLAGS)
endif
#CFLAGSF77:=$(CFLAGSF77) -fno-underscoring
#CFLAGSF90:=$(CFLAGSF90) -fno-underscoring
endif

ifeq ($(call match,$(ID),ia64-itanium2x2-linux-.*-mkl.*),yes)
include $(TEMPLATESDIR)/linux-generic.mk
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

ifeq ($(call match,$(ID),ia64-itanium2x2-linux-.*-.*-ompi),yes)
include $(TEMPLATESDIR)/linux-generic.mk
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/openmpi-generic.mk
# Add preprocessor flag to enable 64bit UMFPACK support.
# UMFPACK documentation states to also add -DLP64, but there is no
# reference to it in the UMFPACK 5.2 sources.
CFLAGSF90     := -DENABLE_64BIT_UMFPACK $(CFLAGSF90)
CFLAGSC       := -DENABLE_64BIT_UMFPACK $(CFLAGSC)
endif



##############################################################################
# NEC
#
##############################################################################

ifeq ($(call match,$(ID),sx6-none-superux-f90-keisan-mpi),yes)
TOKEN5 := 1
TOKEN6 := 1
include $(TEMPLATESDIR)/necsx-generic.mk
include $(TEMPLATESDIR)/necsxcompiler-generic.mk
include $(TEMPLATESDIR)/superux-generic.mk
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
include $(TEMPLATESDIR)/superux-generic.mk
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
include $(TEMPLATESDIR)/superux-generic.mk
#BUILDLIB  = metis blas lapack umfpack amd
BUILDLIB  = metis umfpack amd
#LIBDIR   := $(LIBDIR) -L/SX/opt/MathKeisan/lib
#LIBDIR   := $(LIBDIR) -L/SX/opt/mathkeisan/MK1_6/lib0/lib64
LIBS  := $(LIBS) -L/SX/opt/MathKeisan/lib -llapack -lblas
MPILIBS   = -lmpi
endif




ifeq ($(call match,$(ID),sx8-none-superux-f90-blas-mpi),yes)
TOKEN5 := 1
include $(TEMPLATESDIR)/superux-generic.mk
#BUILDLIB  = metis blas lapack umfpack amd
BUILDLIB  = metis umfpack amd blas lapack
#LIBDIR   := $(LIBDIR) -L/SX/usr/lib0/libp
MPILIBS   = -lmpi
MPILIBS   =
endif



##############################################################################
# SUN Ultrasparc v8
#
##############################################################################

ifeq ($(call match,$(ID),sun4u-sparcv8-sunos-g95-blas.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/sunos-generic.mk
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
include $(TEMPLATESDIR)/sunos-generic.mk
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

ifeq ($(call match,$(ID),sun4u-sparcv8-sunos-sunstudio-.*),yes)
# invokes the Sun Studio Compiler suite, uses LAM MPI, local BLAS,
# local umfpack
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/sunos-generic.mk
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

ifeq ($(call match,$(ID),sun4u-sparcv8-sunos-sunstudio-perf.*),yes)
# Use Perflib BLAS,
include $(TEMPLATESDIR)/sunperf-generic.mk
# Is it necessary to replace the way Sunperf Library is included?
# It works both ways! (SB)
LIBS          := $(filter-out -lsunperf, $(LIBS)) -xlic_lib=sunperf
endif



##############################################################################
# SUN Ultrasparc v9
#
##############################################################################

ifeq ($(call match,$(ID),sun4u-sparcv9-sunos-sunstudio-.*),yes)
# invokes the Sun Studio Compiler suite, local umfpack
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/sunos-generic.mk
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

ifeq ($(call match,$(ID),sun4u-sparcv9-sunos-sunstudio-perf.*),yes)
# Use Perflib BLAS,
include $(TEMPLATESDIR)/sunperf-generic.mk
# Is it necessary to replace the way Sunperf Library is included?
# It works both ways! (SB)
LIBS          := $(filter-out -lsunperf, $(LIBS)) -xlic_lib=sunperf
endif



##############################################################################
# SUN Whatever Solaris 10
#
##############################################################################

ifeq ($(call match,$(ID),sun4v-sparcv9-sunos-sunstudio-.*),yes)
# invokes the Sun Studio Compiler suite, local umfpack
#
# not fully tested settings
MESSAGE := $(MESSAGE) \
	    echo; \
	    echo '*** Warning: The settings for this build ID have not been thoroughly tested yet'; \
	    echo;
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/sunos-generic.mk
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
include $(TEMPLATESDIR)/sunos-generic.mk
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



# # # # # # # # #      IMPORTANT NOTICE !!!     # # # # # # # # #
#
# DO NOT EDIT THE LINES BELOW UNLESS YOU REALLY(!!!) KNOW WHAT
# YOU ARE DOING. AN ERROR IN THE GENERIC RULES WILL MOST PROPABLY
# STOP THE BUILD-SYSTEM TO WORK.
#
# # # # # # # # #      IMPORTANT NOTICE !!!     # # # # # # # # #

ifeq ($(TOKEN3), 0)

##############################################################################
# Operating systems
#
##############################################################################

ifeq ($(call match,$(ID),.*-.*-aix-.*-.*),yes)
# Use IBM Aix as operating system
include $(TEMPLATESDIR)/aix-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-cygwin_nt-.*-.*),yes)
# Use Cygwin NT as operating system
include $(TEMPLATESDIR)/cygwin_nt-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-darwin-.*-.*),yes)
# Use Darwin (Mac OS X) as operating system
include $(TEMPLATESDIR)/darwin-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-hpux-.*-.*),yes)
# Use HP UX as operating system
include $(TEMPLATESDIR)/hpux-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-linux-.*-.*),yes)
# Use Linux as operating system
include $(TEMPLATESDIR)/linux-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-osf1-.*-.*),yes)
# Use DEC OSF1 as operating system
include $(TEMPLATESDIR)/osf1-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-sunos-.*-.*),yes)
# Use SUN Solaris as operating system
include $(TEMPLATESDIR)/sunos-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-superux-.*-.*),yes)
# Use NEC Super-UX as operating system
include $(TEMPLATESDIR)/superux-generic.mk
endif

# End of generic operating system rules
endif



ifeq ($(TOKEN4), 0)

##############################################################################
# Alpha generic
#
##############################################################################

ifeq ($(call match,$(ID),alpha-.*-.*-.*-.*),yes)
MESSAGE  := $(MESSAGE) \
	    echo '*** Error: No generic rule available for Alpha-bases systems';
endif



##############################################################################
# IBM generic
#
##############################################################################

ifeq ($(call match,$(ID),ibm-.*-.*-.*-.*),yes)
MESSAGE  := $(MESSAGE) \
	    echo '*** Error: No generic rule available for IBM systems';
endif



##############################################################################
# x86 / x86_64 generic
#
# This rule tries to find the optimal compiler settings for 
# x86- and x86_64-based systems.
#
# At the time of writing, the following compilers support auto-tuning:
# - GNU GCC 4.2.0 and better
# - Portland Group Fortran Compiler
# - PathScale Fortran Compiler
# - SunStudio Fortran Compiler
#
##############################################################################

ifeq ($(call match,$(ID),(pc|pc64)-.*-.*-g95-.*),yes)
# invokes the G95 Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: No auto-tuning for G95 compiler available!';
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-.*-.*-gcc-.*),yes)
# invokes the gfortran Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call gfortranminversion_4_2),yes)
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by GNU GCC compiler!';
CFLAGSF77LIBS := -march=native $(CFLAGSF77LIBS)
CFLAGSF77     := -march=native $(CFLAGSF77)
CFLAGSF90     := -march=native $(CFLAGSF90)
CFLAGSC       := -march=native $(CFLAGSC)
CFLAGSCXX     := -march=native $(CFLAGSCXX)
LDFLAGS       := -march=native $(LDFLAGS)
else
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: No auto-tuning for GNU GCC below 4.2.0 compiler available!';
endif
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-.*-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xHOST $(CFLAGSF77LIBS)
CFLAGSF77     := -xHOST $(CFLAGSF77)
CFLAGSF90     := -xHOST $(CFLAGSF90)
CFLAGSC       := -xHOST $(CFLAGSC)
CFLAGSCXX     := -xHOST $(CFLAGSCXX)
LDFLAGS       := -xHOST $(LDFLAGS)
else
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: No auto-tuning for Intel compiler below 11.0 available!';
endif
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-.*-.*-open64-.*),yes)
# invokes the Open64 Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/open64-generic.mk
ifeq ($(call optimise), YES)
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: No auto-tuning for Open64 compiler available!';
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-.*-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by Portland Group Fortan compiler!';
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-.*-.*-psc-.*),yes)
# invokes the PathScale Fortran 77/90 compiler, local umfpack, local metis
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/pathscale-generic.mk
ifeq ($(call optimise), YES)
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by PathScale Fortran compiler!';
CFLAGSF77LIBS := -march=auto $(CFLAGSF77LIBS)
CFLAGSF77     := -march=auto $(CFLAGSF77)
CFLAGSF90     := -march=auto $(CFLAGSF90)
CFLAGSC       := -march=auto $(CFLAGSC)
CFLAGSCXX     := -march=auto $(CFLAGSCXX)
CFLAGSCOPROC  := -march=auto $(CFLAGSCOPROC)
LDFLAGS       := -march=auto $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-.*-.*-sunstudio-.*),yes)
# invokes the Sunstudio Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/sunstudio-generic.mk
ifeq ($(call optimise), YES)
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by SunStudio Fortran compiler!';
CFLAGSF77LIBS := -xautopar  $(CFLAGSF77LIBS)
CFLAGSF77     := -xautopar  $(CFLAGSF77)
CFLAGSF90     := -xautopar  $(CFLAGSF90)
CFLAGSC       := -xautopar  $(CFLAGSC)
CFLAGSCXX     := -xautopar  $(CFLAGSCXX)
LDFLAGS       := -xautopar  $(LDFLAGS)
endif
endif



##############################################################################
# i64 (Itanium/Itanium 2) generic
#
# This rule tries to find the optimal compiler settings for 
# i64-based systems.
#
# At the time of writing, the following compilers support auto-tuning:
# - GNU GCC 4.2.0 and better
# - Portland Group Fortran Compiler
# - PathScale Fortran Compiler
# - SunStudio Fortran Compiler
#
##############################################################################

ifeq ($(call match,$(ID),ia64-.*-.*.*-g95-.*),yes)
# invokes the G95 Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/g95-generic.mk
ifeq ($(call optimise), YES)
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: No auto-tuning for G95 compiler available!';
endif
endif

ifeq ($(call match,$(ID),ia64-.*-.*-gcc-.*),yes)
# invokes the gfortran Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/gcc-generic.mk
ifeq ($(call optimise), YES)
# detect compiler version
gcc42better=NO
ifneq (,$(findstring 4.2.,$(GFORTRANVERSION)))
gcc42better=YES
endif
ifneq (,$(findstring 4.3.,$(GFORTRANVERSION)))
gcc42better=YES
endif
ifneq (,$(findstring 4.4.,$(GFORTRANVERSION)))
gcc42better=YES
endif
ifneq (,$(findstring 4.5.,$(GFORTRANVERSION)))
gcc42better=YES
endif

ifeq ($(strip $(gcc42better)), YES)
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by GNU GCC compiler!';
CFLAGSF77LIBS := -march=native $(CFLAGSF77LIBS)
CFLAGSF77     := -march=native $(CFLAGSF77)
CFLAGSF90     := -march=native $(CFLAGSF90)
CFLAGSC       := -march=native $(CFLAGSC)
CFLAGSCXX     := -march=native $(CFLAGSCXX)
LDFLAGS       := -march=native $(LDFLAGS)
else
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: No auto-tuning for GNU GCC below 4.2.0 compiler available!';
endif
endif
endif

ifeq ($(call match,$(ID),ia64-.*-.*-intel-.*),yes)
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/intel-generic.mk
ifeq ($(call optimise), YES)
ifeq ($(call intelminversion_11_0),yes)
CFLAGSF77LIBS := -xHOST $(CFLAGSF77LIBS)
CFLAGSF77     := -xHOST $(CFLAGSF77)
CFLAGSF90     := -xHOST $(CFLAGSF90)
CFLAGSC       := -xHOST $(CFLAGSC)
CFLAGSCXX     := -xHOST $(CFLAGSCXX)
LDFLAGS       := -xHOST $(LDFLAGS)
else
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: No auto-tuning for Intel compiler below 11.0 available!';
endif
endif
endif

ifeq ($(call match,$(ID),ia64-.*-.*-pgi-.*),yes)
# invokes the Portland Group Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/pgi-generic.mk
ifeq ($(call optimise), YES)
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Auto-tuning by Portland Group Fortan compiler!';
endif
endif

ifeq ($(call match,$(ID),ia64-.*-.*-psc-.*),yes)
# invokes the PathScale Fortran 77/90 compiler, local umfpack, local metis
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/pathscale-generic.mk
ifeq ($(call optimise), YES)
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by PathScale Fortran compiler!';
CFLAGSF77LIBS := -march=auto $(CFLAGSF77LIBS)
CFLAGSF77     := -march=auto $(CFLAGSF77)
CFLAGSF90     := -march=auto $(CFLAGSF90)
CFLAGSC       := -march=auto $(CFLAGSC)
CFLAGSCXX     := -march=auto $(CFLAGSCXX)
CFLAGSCOPROC  := -march=auto $(CFLAGSCOPROC)
LDFLAGS       := -march=auto $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),ia64-.*-.*-sunstudio-.*),yes)
# invokes the Sunstudio Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/sunstudio-generic.mk
ifeq ($(call optimise), YES)
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by SunStudio Fortran compiler!';
CFLAGSF77LIBS := -xautopar  $(CFLAGSF77LIBS)
CFLAGSF77     := -xautopar  $(CFLAGSF77)
CFLAGSF90     := -xautopar  $(CFLAGSF90)
CFLAGSC       := -xautopar  $(CFLAGSC)
CFLAGSCXX     := -xautopar  $(CFLAGSCXX)
LDFLAGS       := -xautopar  $(LDFLAGS)
endif
endif



##############################################################################
# NEC generic
#
##############################################################################

ifeq ($(call match,$(ID),sx.*-.*-.*-.*-.*),yes)
MESSAGE  := $(MESSAGE) \
	    echo '*** Error: No generic rule available for NEC systems';
endif

# End of generic compiler rules
endif



##############################################################################
# BLAS/LAPACK implementations
#
##############################################################################

ifeq ($(TOKEN5), 0)

# Include BLAS/Lapack libraries the default way

ifeq ($(call match,$(ID),.*-.*-.*-.*-acml.*),yes)
# Use AMD Core Math Library
include $(TEMPLATESDIR)/acml-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-atlas.*),yes)
# Use ATLAS BLAS
include $(TEMPLATESDIR)/atlasblas-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/blas-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-goto\(\$\|-.*\))),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/gotoblas-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/gotoblas2-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-mkl.*),yes)
# Use MKL BLAS
include $(TEMPLATESDIR)/mkl-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-perf.*),yes)
# Use Perflib BLAS
include $(TEMPLATESDIR)/sunperf-generic.mk
endif

# End of generic BLAS/LAPACK rules
endif



##############################################################################
# MPI environments
#
##############################################################################

ifeq ($(TOKEN6), 0)

# Include MPI libraries the default way

ifeq ($(call match,$(ID),.*-.*-.*-.*-.*-lammpi),yes)
# Use LAM/MPI
include $(TEMPLATESDIR)/lammpi-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-.*-lanlmpi),yes)
# Use LANL MPI
include $(TEMPLATESDIR)/lanlmpi-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-.*-mpich),yes)
# Use MPICH
include $(TEMPLATESDIR)/mpich-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-.*-mpich2),yes)
# Use MPICH2
include $(TEMPLATESDIR)/mpich2-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-.*-mvapich),yes)
# Use MVAPICH
include $(TEMPLATESDIR)/mvapich-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/openmpi-generic.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-.*-optmpich),yes)
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

# End of generic MPI rules
endif



##############################################################################
# Check for forced 32bit/64bit binary
##############################################################################
ifeq ($(strip $(BINARY)), 32)
CFLAGSF77LIBS := -m32 $(CFLAGSF77LIBS)
CFLAGSF77     := -m32 $(CFLAGSF77)
CFLAGSF90     := -m32 $(CFLAGSF90)
CFLAGSC       := -m32 $(CFLAGSC)
CFLAGSCXX     := -m32 $(CFLAGSCXX)
LDFLAGS       := -m32 $(LDFLAGS)
endif

ifeq ($(strip $(BINARY)), 64)
CFLAGSF77LIBS := -m64 $(CFLAGSF77LIBS)
CFLAGSF77     := -m64 $(CFLAGSF77)
CFLAGSF90     := -m64 $(CFLAGSF90)
CFLAGSC       := -m64 $(CFLAGSC)
CFLAGSCXX     := -m64 $(CFLAGSCXX)
LDFLAGS       := -m64 $(LDFLAGS)
endif



##############################################################################
# auxiliary targets: .id, .idonly
#
##############################################################################


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
	@echo '  CPPFLAGS     =' $(CPPFLAGS)
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
