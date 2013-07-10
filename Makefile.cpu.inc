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
FILE_CUDA_INFO    = env_cuda_compiler
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
ifeq ($(strip $(OPT)), NATIVE-EXPENSIVE)
optimise = YES
endif
ifeq ($(strip $(OPT)), NATIVE)
optimise = YES
endif
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
CPPFLAGS      = -x assembler-with-cpp $(filter -D% -U% '-D% '-U%, $(APPONLYFLAGS) $(CFLAGSF90))

# FEAT2's preprocessor
F90CPP        = $(FEAT2BASEDIR)/bin/f90cpp

# Fortran 77 compiler
F77           = echo; echo "No Fortran 77 compiler specified."; echo; exit 1; echo
# Fortran 77 compiler options
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
INC           = -I$(FEAT2BASEDIR)/librariesNEW/coproc

# List of files only needed for a specific build ID
# (please include path information, at best using Makefile variables)
SRCEXTRA      = $(FEAT2BASEDIR)/kernel/Postprocessing/gmvwritef.c \
                $(FEAT2BASEDIR)/kernel/System/high_res_clock.c \
		$(FEAT2BASEDIR)/kernel/System/signal_ccode.c \
		$(FEAT2BASEDIR)/kernel/System/storage_malloc.c

# List of all libraries available (queried for purge* targets)
ALLLIB        = agmg amd blas coproc inshape2d lapack metis splib umfpack

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
# FEAT2 can be linked with the external AGMG library by Y.Notay

##############################################################################

ifeq ($(AGMG), YES)
BUILDLIB := agmg $(BUILDLIB)
LIBS     := -lagmg $(LIBS)
APPONLYFLAGS := -DENABLE_AGMG $(APPONLYFLAGS)
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
# Token 1 and 2: Architecture and CPU-type
#
##############################################################################

# Intel Itanium IA-64 architecture
ifeq ($(call match,$(ID),ia64-.*-.*-.*-.*),yes)
include $(TEMPLATESDIR)/ARCH_CPU-ia64.mk
endif

# Intel x86 / x86_64 architecture
ifeq ($(call match,$(ID),(pc|pc64)-.*-.*-.*-.*),yes)
include $(TEMPLATESDIR)/ARCH_CPU-pc.mk
endif

# Sun UltraSPARC I, II, III, IV and UltraSPARC T1, T2 CPUs
ifeq ($(call match,$(ID),(sun4u|sun4v)-.*-.*-.*-.*),yes)
include $(TEMPLATESDIR)/ARCH_CPU-sun4u.mk
endif

# NEC SX machines
ifeq ($(call match,$(ID),(sx6|sx8)-.*-.*-.*-.*),yes)
include $(TEMPLATESDIR)/ARCH_CPU-necsx.mk
endif

##############################################################################
# Token 3: Operating System
#
##############################################################################

# IBM AIX unix OS
ifeq ($(call match,$(ID),.*-.*-aix-.*-.*),yes)
include $(TEMPLATESDIR)/OS-aix.mk
endif

# Cygwin OS
ifeq ($(call match,$(ID),.*-.*-cygwin_nt-.*-.*),yes)
include $(TEMPLATESDIR)/OS-cygwin_nt.mk
endif

# Darwin OS (aka Mac OS X)
ifeq ($(call match,$(ID),.*-.*-darwin-.*-.*),yes)
include $(TEMPLATESDIR)/OS-darwin.mk
endif

# HP-UX OS
ifeq ($(call match,$(ID),.*-.*-hpux-.*-.*),yes)
include $(TEMPLATESDIR)/OS-hpux.mk
endif

# Linux OS
ifeq ($(call match,$(ID),.*-.*-linux-.*-.*),yes)
include $(TEMPLATESDIR)/OS-linux.mk
endif

# DEC OSF1 OS
ifeq ($(call match,$(ID),.*-.*-osf1-.*-.*),yes)
include $(TEMPLATESDIR)/OS-osf1.mk
endif

# Sun/Oracle SunOS
ifeq ($(call match,$(ID),.*-.*-sunos-.*-.*),yes)
include $(TEMPLATESDIR)/OS-sunos.mk
endif

# Super-UX (used on NEC SX)
ifeq ($(call match,$(ID),.*-.*-superux-.*-.*),yes)
include $(TEMPLATESDIR)/OS-superux.mk
endif



##############################################################################
# Token 4: Compiler
#
##############################################################################

# DEC compiler suite
ifeq ($(call match,$(ID),.*-.*-.*-decf90-.*),yes)
include $(TEMPLATESDIR)/COMP-decf90.mk
endif

# G95 Fortran compiler
ifeq ($(call match,$(ID),.*-.*-.*-g95-.*),yes)
include $(TEMPLATESDIR)/COMP-g95.mk
ifeq ($(call optimise), YES)
# Enable optimization flags for x86 / x86_64
ifeq ($(call match,$(ID),(pc|pc64)-*-.*-g95-.*),yes)
include $(TEMPLATESDIR)/CFLAGSOPT-pc-g95.mk
endif
# Enable optimization flags for IA-64
ifeq ($(call match,$(ID),ia64-*-.*-g95-.*),yes)
include $(TEMPLATESDIR)/CFLAGSOPT-ia64-g95.mk
endif
# Optimization flags for further architectures go here
# ...
CFLAGSF77 := $(CFLAGSOPTF77) $(CFLAGSF77)
CFLAGSF90 := $(CFLAGSOPTF90) $(CFLAGSF90)
CFLAGSC   := $(CFLAGSOPTC)   $(CFLAGSC)
CFLAGSCXX := $(CFLAGSOPTCXX) $(CFLAGSCXX)
LDFLAGS   := $(LDFLAGSOPT)   $(LDFLAGS)
endif
endif

# GNU compiler suite
ifeq ($(call match,$(ID),.*-.*-.*-gcc-.*),yes)
include $(TEMPLATESDIR)/COMP-gcc.mk
ifeq ($(call optimise), YES)
# Enable optimization flags for x86 / x86_64
ifeq ($(call match,$(ID),(pc|pc64)-*-.*-gcc-.*),yes)
include $(TEMPLATESDIR)/CFLAGSOPT-pc-gcc.mk
endif
# Enable optimization flags for IA-64
ifeq ($(call match,$(ID),ia64-*-.*-gcc-.*),yes)
include $(TEMPLATESDIR)/CFLAGSOPT-ia64-gcc.mk
endif
# Optimization flags for further architectures go here
# ...
CFLAGSF77 := $(CFLAGSOPTF77) $(CFLAGSF77)
CFLAGSF90 := $(CFLAGSOPTF90) $(CFLAGSF90)
CFLAGSC   := $(CFLAGSOPTC)   $(CFLAGSC)
CFLAGSCXX := $(CFLAGSOPTCXX) $(CFLAGSCXX)
LDFLAGS   := $(LDFLAGSOPT)   $(LDFLAGS)
endif
endif

# Intel compiler suite
ifeq ($(call match,$(ID),.*-.*-.*-intel-.*),yes)
include $(TEMPLATESDIR)/COMP-intel.mk
ifeq ($(call optimise), YES)
# Enable optimization flags for x86 / x86_64
ifeq ($(call match,$(ID),(pc|pc64)-*-.*-intel-.*),yes)
include $(TEMPLATESDIR)/CFLAGSOPT-pc-intel.mk
endif
# Optimization flags for further architectures go here
# ...
CFLAGSF77 := $(CFLAGSOPTF77) $(CFLAGSF77)
CFLAGSF90 := $(CFLAGSOPTF90) $(CFLAGSF90)
CFLAGSC   := $(CFLAGSOPTC)   $(CFLAGSC)
CFLAGSCXX := $(CFLAGSOPTCXX) $(CFLAGSCXX)
LDFLAGS   := $(LDFLAGSOPT)   $(LDFLAGS)
endif
endif

# Intel compiler suite for Itanium
ifeq ($(call match,$(ID),.*-.*-.*-intel64-.*),yes)
include $(TEMPLATESDIR)/COMP-intel64.mk
ifeq ($(call optimise), YES)
# Enable optimization flags for i64
ifeq ($(call match,$(ID),ia64-*-.*-intel64-.*),yes)
include $(TEMPLATESDIR)/CFLAGSOPT-ia-intel64.mk
CFLAGSF77 := $(CFLAGSOPTF77) $(CFLAGSF77)
CFLAGSF90 := $(CFLAGSOPTF90) $(CFLAGSF90)
CFLAGSC   := $(CFLAGSOPTC)   $(CFLAGSC)
CFLAGSCXX := $(CFLAGSOPTCXX) $(CFLAGSCXX)
LDFLAGS   := $(LDFLAGSOPT)   $(LDFLAGS)
endif
endif
endif

# NEC SX compiler suite
ifeq ($(call match,$(ID),.*-.*-.*-necsx-.*),yes)
include $(TEMPLATESDIR)/COMP-necsx.mk
endif

# Open64 compiler suite
ifeq ($(call match,$(ID),.*-.*-.*-open64-.*),yes)
include $(TEMPLATESDIR)/COMP-open64.mk
ifeq ($(call optimise), YES)
# Enable optimization flags for x86 / x86_64
ifeq ($(call match,$(ID),(pc|pc64)-*-.*-open64-.*),yes)
include $(TEMPLATESDIR)/CFLAGSOPT-pc-open64.mk
endif
# Optimization flags for further architectures go here
# ...
CFLAGSF77 := $(CFLAGSOPTF77) $(CFLAGSF77)
CFLAGSF90 := $(CFLAGSOPTF90) $(CFLAGSF90)
CFLAGSC   := $(CFLAGSOPTC)   $(CFLAGSC)
CFLAGSCXX := $(CFLAGSOPTCXX) $(CFLAGSCXX)
LDFLAGS   := $(LDFLAGSOPT)   $(LDFLAGS)
endif
endif

# PathScale compiler suite
ifeq ($(call match,$(ID),.*-.*-.*-psc-.*),yes)
include $(TEMPLATESDIR)/COMP-pathscale.mk
ifeq ($(call optimise), YES)
# Enable optimization flags for x86 / x86_64
ifeq ($(call match,$(ID),(pc|pc64)-*-.*-psc-.*),yes)
include $(TEMPLATESDIR)/CFLAGSOPT-pc-pathscale.mk
endif
# Optimization flags for further architectures go here
# ...
CFLAGSF77 := $(CFLAGSOPTF77) $(CFLAGSF77)
CFLAGSF90 := $(CFLAGSOPTF90) $(CFLAGSF90)
CFLAGSC   := $(CFLAGSOPTC)   $(CFLAGSC)
CFLAGSCXX := $(CFLAGSOPTCXX) $(CFLAGSCXX)
LDFLAGS   := $(LDFLAGSOPT)   $(LDFLAGS)
endif
endif

# PGI compiler suite
ifeq ($(call match,$(ID),.*-.*-.*-pgi-.*),yes)
include $(TEMPLATESDIR)/COMP-pgi.mk
ifeq ($(call optimise), YES)
# Enable optimization flags for x86 / x86_64
ifeq ($(call match,$(ID),(pc|pc64)-*-.*-pgi-.*),yes)
include $(TEMPLATESDIR)/CFLAGSOPT-pc-pgi.mk
endif
# Optimization flags for further architectures go here
# ...
CFLAGSF77 := $(CFLAGSOPTF77) $(CFLAGSF77)
CFLAGSF90 := $(CFLAGSOPTF90) $(CFLAGSF90)
CFLAGSC   := $(CFLAGSOPTC)   $(CFLAGSC)
CFLAGSCXX := $(CFLAGSOPTCXX) $(CFLAGSCXX)
LDFLAGS   := $(LDFLAGSOPT)   $(LDFLAGS)
endif
endif

# SunStudio/OracleStudio compiler suite
ifeq ($(call match,$(ID),.*-.*-.*-sunstudio-.*),yes)
include $(TEMPLATESDIR)/COMP-sunstudio.mk
ifeq ($(call optimise), YES)
# Enable optimization flags for sun4u / sun4v
ifeq ($(call match,$(ID),(sun4u|sun4v)-*-.*-sunstudio-.*),yes)
include $(TEMPLATESDIR)/CFLAGSOPT-sun4u-sunstudio.mk
endif
# Enable optimization flags for x86 / x86_64
ifeq ($(call match,$(ID),(pc|pc64)-*-.*-sunstudio-.*),yes)
include $(TEMPLATESDIR)/CFLAGSOPT-pc-sunstudio.mk
endif
CFLAGSF77 := $(CFLAGSOPTF77) $(CFLAGSF77)
CFLAGSF90 := $(CFLAGSOPTF90) $(CFLAGSF90)
CFLAGSC   := $(CFLAGSOPTC)   $(CFLAGSC)
CFLAGSCXX := $(CFLAGSOPTCXX) $(CFLAGSCXX)
LDFLAGS   := $(LDFLAGSOPT)   $(LDFLAGS)
endif
endif

# IBM Xlf compiler
ifeq ($(call match,$(ID),.*-.*-.*-xlf-.*),yes)
include $(TEMPLATESDIR)/COMP-xlf.mk
endif


##############################################################################
# Token 5: BLAS/LAPACK implementations
#
##############################################################################

ifeq ($(TOKEN5), 0)

# Include BLAS/Lapack libraries the default way

ifeq ($(call match,$(ID),.*-.*-.*-.*-acml.*),yes)
# Use AMD Core Math Library
include $(TEMPLATESDIR)/BLAS-acml.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-atlas.*),yes)
# Use ATLAS BLAS
include $(TEMPLATESDIR)/BLAS-atlasblas.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-blas.*),yes)
# Use self-compiled local BLAS
include $(TEMPLATESDIR)/BLAS-blas.mk
endif

# sed syntax for match variant 1:
ifeq ($(call match,$(ID),.*-.*-.*-.*-goto($$|-.*)),yes)
# sed syntax for match variant 2:
# ifeq ($(call match,$(ID),.*-.*-.*-.*-goto\(\$\|-.*\))),yes)
# Use Goto BLAS
include $(TEMPLATESDIR)/BLAS-gotoblas.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-goto2.*),yes)
# Use Goto BLAS 2
include $(TEMPLATESDIR)/BLAS-gotoblas2.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-mkl.*),yes)
# Use MKL BLAS
include $(TEMPLATESDIR)/BLAS-mkl.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-openblas.*),yes)
# Use Open BLAS
include $(TEMPLATESDIR)/BLAS-openblas.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-perf.*),yes)
# Use Perflib BLAS
include $(TEMPLATESDIR)/BLAS-sunperf.mk
endif

endif



##############################################################################
# Token 6: MPI environments
#
##############################################################################

ifeq ($(TOKEN6), 0)

# Include MPI libraries the default way

ifeq ($(call match,$(ID),.*-.*-.*-.*-.*-lammpi),yes)
# Use LAM/MPI
include $(TEMPLATESDIR)/MPI-lammpi.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-.*-lanlmpi),yes)
# Use LANL MPI
include $(TEMPLATESDIR)/MPI-lanlmpi.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-.*-mpich),yes)
# Use MPICH
include $(TEMPLATESDIR)/MPI-mpich.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-.*-mpich2),yes)
# Use MPICH2
include $(TEMPLATESDIR)/MPI-mpich2.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-.*-mvapich),yes)
# Use MVAPICH
include $(TEMPLATESDIR)/MPI-mvapich.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-.*-ompi),yes)
# Use OpenMPI
include $(TEMPLATESDIR)/MPI-openmpi.mk
endif

ifeq ($(call match,$(ID),.*-.*-.*-.*-.*-optmpich),yes)
# Use Allinea Opt MPI
include $(TEMPLATESDIR)/MPI-optmpich.mk
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
# Add CUDA support
##############################################################################
HAS_CUDA10 := $(shell echo $(HAS_CUDA10) | tr '[a-z]' '[A-Z]')
HAS_CUDA11 := $(shell echo $(HAS_CUDA11) | tr '[a-z]' '[A-Z]')
HAS_CUDA12 := $(shell echo $(HAS_CUDA12) | tr '[a-z]' '[A-Z]')
HAS_CUDA13 := $(shell echo $(HAS_CUDA13) | tr '[a-z]' '[A-Z]')
HAS_CUDA20 := $(shell echo $(HAS_CUDA20) | tr '[a-z]' '[A-Z]')
HAS_CUDA21 := $(shell echo $(HAS_CUDA21) | tr '[a-z]' '[A-Z]')
HAS_CUDA30 := $(shell echo $(HAS_CUDA30) | tr '[a-z]' '[A-Z]')
ifneq (,$(findstring YES,$(HAS_CUDA10)$(HAS_CUDA11)$(HAS_CUDA12)$(HAS_CUDA13)$(HAS_CUDA20)$(HAS_CUDA21)$(HAS_CUDA30)))
include $(TEMPLATESDIR)/COMP-nvcc.mk
endif


##############################################################################
# Check for forced 32bit/64bit binary
##############################################################################
ifeq ($(strip $(BINARY)), 32)
CFLAGSF77     := -m32 $(CFLAGSF77)
CFLAGSF90     := -m32 $(CFLAGSF90)
CFLAGSC       := -m32 $(CFLAGSC)
CFLAGSCXX     := -m32 $(CFLAGSCXX)
LDFLAGS       := -m32 $(LDFLAGS)
endif

ifeq ($(strip $(BINARY)), 64)
CFLAGSF77     := -m64 $(CFLAGSF77)
CFLAGSF90     := -m64 $(CFLAGSF90)
CFLAGSC       := -m64 $(CFLAGSC)
CFLAGSCXX     := -m64 $(CFLAGSCXX)
LDFLAGS       := -m64 $(LDFLAGS)
endif


##############################################################################
# Check whether to enable 32bit or 64bit UMFPACK interface
##############################################################################
ifeq ($(strip $(INTSIZE)), LARGE)
# Add preprocessor flag to enable 64bit UMFPACK support.
# UMFPACK documentation states to also add -DLP64, but there is no
# reference to it in the UMFPACK 5.2-5.6.2 sources.
CFLAGSF90     := -DENABLE_64BIT_UMFPACK $(CFLAGSF90)
CFLAGSC       := -DENABLE_64BIT_UMFPACK $(CFLAGSC)
endif


##############################################################################
# Check for enabled Warnings
##############################################################################
ifeq ($(strip $(WARN_DEBUG)),YES)
CFLAGSF77     := -DWARN_DEBUG $(CFLAGSF77)
CFLAGSF90     := -DWARN_DEBUG $(CFLAGSF90)
CFLAGSC       := -DWARN_DEBUG $(CFLAGSC)
CFLAGSCXX     := -DWARN_DEBUG $(CFLAGSCXX)
endif

ifeq ($(strip $(WARN_DEPRECATED)),YES)
CFLAGSF77     := -DWARN_DEPREC $(CFLAGSF77)
CFLAGSF90     := -DWARN_DEPREC $(CFLAGSF90)
CFLAGSC       := -DWARN_DEPREC $(CFLAGSC)
CFLAGSCXX     := -DWARN_DEPREC $(CFLAGSCXX)
endif

ifeq ($(strip $(WARN_UNSUPPORTED)),YES)
CFLAGSF77     := -DWARN_UNSUPP $(CFLAGSF77)
CFLAGSF90     := -DWARN_UNSUPP $(CFLAGSF90)
CFLAGSC       := -DWARN_UNSUPP $(CFLAGSC)
CFLAGSCXX     := -DWARN_UNSUPP $(CFLAGSCXX)
endif

ifeq ($(strip $(WARN_UNUSUAL)),YES)
CFLAGSF77     := -DWARN_UNUSUAL $(CFLAGSF77)
CFLAGSF90     := -DWARN_UNUSUAL $(CFLAGSF90)
CFLAGSC       := -DWARN_UNUSUAL $(CFLAGSC)
CFLAGSCXX     := -DWARN_UNUSUAL $(CFLAGSCXX)
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
	    echo '  Fortran 77 compiler:' $(PREF77) $(shell (which 2>/dev/null $(F77) || \
		echo "\*\*\* $(F77) not found ! \*\*\*")); \
	else \
	    echo '  Fortran 77 compiler: *** not set ! ***'; \
	fi
	@if test -z "$(findstring exit 1, $(F90))"; then \
	    echo '  Fortran 90 compiler:' $(PREF90) $(shell (which 2>/dev/null $(F90) || \
		echo "\*\*\* $(F90) not found ! \*\*\*")); \
	else \
	    echo '  Fortran 90 compiler: *** not set ! ***'; \
	fi
	@if test -z "$(findstring exit 1, $(CC))"; then \
	    echo '  C compiler         :' $(PRECC) $(shell (which 2>/dev/null $(CC)  || \
		echo "\*\*\* $(CC) not found ! \*\*\*")); \
	else \
	    echo '  C compiler         : *** not set ! ***'; \
	fi
	@if test -z "$(findstring exit 1, $(CXX))"; then \
	    echo '  C++ compiler       :' $(PRECXX) $(shell (which 2>/dev/null $(CXX) || \
		echo "\*\*\* $(CXX) not found ! \*\*\*")); \
	else \
	    echo '  C++ compiler       : *** not set ! ***'; \
	fi
        # If cuda switch -DHAS_CUDAxx does occur in compiler flags,
        # show CUDA compiler
ifneq (,$(findstring YES,$(HAS_CUDA10)$(HAS_CUDA11)$(HAS_CUDA12)$(HAS_CUDA13)$(HAS_CUDA20)$(HAS_CUDA21)$(HAS_CUDA30)))
	@if test -z "$(findstring exit 1, $(CUDA))"; then \
	    echo '  CUDA compiler      :' $(PRECUDA) $(shell (which 2>/dev/null $(CUDA) || \
		echo "\*\*\* $(CUDA) not found ! \*\*\*")); \
	else \
	    echo '  CUDA compiler      : *** not set ! ***'; \
	fi
endif
	@echo
	@echo 'Linker to be used:'
	@if test -z "$(findstring exit 1, $(LD))"; then \
	    echo '  Linker             :' $(PRELD) $(shell (which 2>/dev/null $(LD)  || \
		echo "\*\*\* $(LD) not found ! \*\*\*")); \
	else \
	    echo '  Linker             : *** not set ! ***'; \
	fi
	@echo
	@echo 'Flags to be used:'
	@echo '  AGMG         =' $(AGMG)
	@echo '  APPNAME      =' $(APPNAME)
	@echo '  APPONLYFLAGS =' $(APPONLYFLAGS)
	@echo '  AR           =' $(AR)
	@echo '  BINARY       =' $(BINARY)
	@echo '  BUILDTARGET  =' $(BUILDTARGET)
	@echo '  BUILDLIB     =' $(BUILDLIB)
	@echo '  CPPFLAGS     =' $(CPPFLAGS)
	@echo '  CFLAGSF77    =' $(CFLAGSF77)
	@echo '  CFLAGSF90    =' $(CFLAGSF90) $(MODOPTION)$(OBJDIR)
	@echo '  CFLAGSC      =' $(CFLAGSC)
	@echo '  CFLAGSCXX    =' $(CFLAGSCXX)
ifneq (,$(findstring YES,$(HAS_CUDA10)$(HAS_CUDA11)$(HAS_CUDA12)$(HAS_CUDA13)$(HAS_CUDA20)$(HAS_CUDA21)$(HAS_CUDA30)))
	@echo '  CFLAGSCUDA   =' $(CFLAGSCUDA)
endif
ifneq ($(strip $(COPROCLIST)),)
	@echo '  COPROCLIST   =' $()
	@echo '  CFLAGSCOPROC =' $(CFLAGSCOPROC)
endif
	@echo '  FORCE_ID     =' $(FORCE_ID)
	@echo '  INC          =' $(INC)
	@echo '  INTSIZE      =' $(INTSIZE)
	@echo '  LDFLAGS      =' $(LDFLAGS)
	@echo '  LIBDIR       =' $(LIBDIR)
	@echo '  LIBS         =' $(LIBS)
	@echo '  MODE         =' $(MODE)
	@echo '  MODEXTENSION =' $(MODEXTENSION)
	@echo '  MOVEMOD      =' $(MOVEMOD)
        # If preprocessor switch -DENABLE_SERIAL_BUILD does occur in compiler flags,
        # a build for serial execution is requested.
	@if test -z "$(findstring -DENABLE_SERIAL_BUILD ,$(APPONLYFLAGS) $(CFLAGSF90) )"; then \
	    echo '  MPI          =' $(MPI); \
	    echo '  MPIINC       =' $(MPIINC); \
	    echo '  MPILIBDIR    =' $(MPILIBDIR); \
	    echo '  MPILIBS      =' $(MPILIBS); \
	    echo '  MPIWRAPPERS  =' $(MPIWRAPPERS); \
	fi
	@echo '  OBJDIR       =' $(OBJDIR)
	@echo '  OBJDIR_LIB   =' $(OBJDIR_LIB)
	@echo '  OPENMP       =' $(OPENMP)
	@echo '  OPT          =' $(OPT)
	@echo '  RANLIB       =' $(RANLIB)
	@echo '  RELINK_ALWAYS=' $(RELINK_ALWAYS)
	@echo '  SRCEXTRA     =' $(SRCEXTRA)
	@echo '  SRCEXTRA_APP =' $(SRCEXTRA_APP)
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
