#!/usr/bin/env gmake
## -*- makefile -*-

########################################################################
# Every makefile that includes this Globals.mk-file has to define
# a variable FEATFLOW pointing to the path of the Featflow installation
# (i.e. the directory containing the Globals.mk file).
# The makefiles in this installation specify the path relatively
# to the current path (e.g. "../../") to allow the whole installation
# directory being moved anywhere. If the user copies an application
# directory to another location, the FEATFLOW variable in the
# makefile has to be modified properly!
########################################################################

# For version management we define a version-tag of the Featflow-version
# here - this is used mainly for output on screen.

FFVER:=2.0ALPHA

########################################################################
# machine and architecture info, machine ID in the form arch-cpu-os
#
# The machine-ID is guessed by the shell script "guess_id" in the
# bin/ subdirectory.
# The script match_id is used to match the ID agains wildcards.
########################################################################

#ID:=$(shell $(FEATFLOW)/bin/guess_id)
FFHOST:=$(shell uname -n)
match= $(shell $(FEATFLOW)/bin/match_id $(1) '$(2)')

ifndef FFARCH
FFARCH:=$(shell $(FEATFLOW)/bin/guess_id 1)
endif

ifndef FFCPU
FFCPU:=$(shell $(FEATFLOW)/bin/guess_id 2)
endif

ifndef FFCORE
FFCORE:=$(shell $(FEATFLOW)/bin/guess_id 3)
endif

FFOS:=$(shell $(FEATFLOW)/bin/guess_id 4)

ifndef ID
ID:=${FFARCH}-${FFCPU}-${FFOS}
endif

########################################################################
# There is a possibility to overide the autodetected ID.  If you know
# your PC while the script is not guessing it correctly, you can either
# modify this file or use this, e.g.
#
#       make build ID=pc-pentium4-linux
#
# If you have a Pentium 5 or similar ( :-) ) and you want to force the
# script to take the settings of the pentium 4 on linux
########################################################################

########################################################################
# Extend the id to arch-cpu-os-$(ALT) if variable ALT is defined.
#
# This way multiple compilers on one machine can be supported. The
# User has then to call MAKE with the appropriate alternative target,
# e.g.
#       make build ALT=ifc
#
# If the machine-ID is for example detected as pc-athlonxp-linux,
# the new machine-ID will be pc-athlonxp-linux-ifc, thus using the
# Intel Fortran Compiler with the settings defined below.
########################################################################

ifdef ALT
IDORIG:=${ID}
ID:=$(ID)-$(ALT)
endif

########################################################################
# Default target .log-file for benchmark results
########################################################################

BENCHLOGFILENAME:=bench-$(HOST)-$(ID).log
BENCHLOG:=$(FEATFLOW)/$(BENCHLOGFILENAME)

########################################################################
# Global path to libraries of current machine.
########################################################################

LIBDIR=$(FEATFLOW)/object/libraries/lib-$(ID)

# list of application modules to create by the top level make
APPS:= $(shell ls $(FEATFLOW)/applications)

# list of all library modules available at the top level
LIBS= feat3d feat2d sysutils umfpack2 amd umfpack4 minisplib lapack blas

ifeq ($(HDF5),YES)
LIBS:=$(LIBS) hdf5_fortran hdf5 zlib sz 
endif

########################################################################
# General name and location of compilers.
#
# These entries are machine-dependend and might be necessary to be
# modified by the user. The entries here are the general ones (!), used
# IF THE MACHINE COULD NOT BE DETECTED. If the script is able to
# detect a known machine, the machine-specific configuration branch
# (see below) is used.
#
# Featflow needs an f90-compiler as well as a C-compiler.
#
# The following entries are defined:
# FC = Fortran 77 and 90 compiler
# CC = ANSI C compiler
# CPP= C++ compiler. If set to "", this defaults to CC
# LD = linker (usually same as the Fortran Compiler) 
# AR = library creation tool (archiver)
# ARC=  for C libraries, if undifined AR is used
########################################################################

FC=f90     # Fortran 90 compiler
CC=cc      # ANSI C compiler
CPP=       # C++ compiler. Defaults to CC if not defined
LD=$(FC)   # linker (usually same as the Fortran Compiler) 
AR=ar      # library creation tool (archiver)
ARC=       # for C libraries, if undifined AR is used

MPIFC=mpif90 # Fortran 90 compiler with MPI support

########################################################################
# compiler settings for various machine types, list to be extended
# make sure to include some BLAS implementation like the one provided
########################################################################

# default values for general arch 

# special compiler optimization flags for both, Fortran and C compiler:

OPTFLAGS=

# special compiler optimization flags for Fortran compiler

OPTFLAGSF=

# special compiler optimization flags for C compiler

OPTFLAGSC=

# special compiler optimization flags for C++ compiler

OPTFLAGSCPP=

# standard flags applied when 'make debug' is executed,
# for both, Fortran and C compiler

OPTFLAGSDEBUG= -g

# standard flags applied when 'make debug' is executed,
# for Fortran compiler; additionally to OPTFLAGSDEBUG!

OPTFLAGSFDEBUG= 

# standard flags applied when 'make debug' is executed,
# for C compiler; additionally to OPTFLAGSDEBUG!

OPTFLAGSCDEBUG= 

# standard flags applied when 'make debug' is executed,
# for C++ compiler; additionally to OPTFLAGSDEBUG!

OPTFLAGSCPPDEBUG= 

# general compiler options for Fortran compiler:

FCFLAGS=

# general compiler options for C compiler:

CCFLAGS=

# general compiler options for C++ compiler:

CPPFLAGS=

# list of featflow included libs to be build by the top level make
# feat2d, feat3d and sysutils required, include lapack and blas if necessary

BUILDLIB = $(LIBS)

# BLAS library to be used, if left empty the included blas and lapack
# is used automatically; to use a single combined BLAS/LAPACK
# library (e.g. ATLAS), set BLASLIB=LAPACKLIB!

BLASLIB  =          # will use included blas
LAPACKLIB=          # will use included lapack

# additional linker options
LDFLAGS=

########################################################################
# additional, system specific libraries, e.g. POSIX,... if necessary
#
# Has to be specified on some systems to include libraries that contain
# standard Fortran intrinsic routines. E.g. if the Intel Fortran
# compiler is not installed properly, eventually the name+path to the
# library libPEPCF90.a has to be spcified here; otherwise the routine
# ETIME might not be found by the linker!
# Example:  LDLIBS=/opt/intel/compiler70/ia32/lib/libPEPCF90.a
########################################################################

LDLIBS=

########################################################################
# This block redefines machine dependent compiler options and
# libraries (optimized blas/lapack versions for example), depending 
# on the current machine-ID. 
#
# If external blas/lapack is used in BLASLIB/LAPACKLIB then blas/
# lapack can be omitted in BUILDLIB to save some compilation time.
# When BLASLIB=LAPACKLIB is set, one single common BLAS/LAPACK library
# is used as specified in both, BLASLIB and LAPACKLIB.
#
# You can detect the id of your machine by typing
#   "make id"
# in the Featflow installation directory. The detected machine-id of
# your computer will then be printed on the screen at the end of the
# line denoted by "Machine-ID (xxx) : ".
# Remember that you can modify/extend your machine-id by an alternative
# setting. E.g. "make id ALT=ifc" will add the string "-ifc" to the
# machine-id before evaluating it and printing it on screen.
#
# The GENERIC flag is used to detect if a generic compiler configuration
# is adopted or if some detailed specification for the employed cpu
# was found
########################################################################

GENERIC = 1

ifeq ($(call match,$(FFARCH),pc),yes)
include $(FEATFLOW)/Globals.x86
endif

ifeq ($(call match,$(FFARCH),pc64),yes)
include $(FEATFLOW)/Globals.x86_64
endif

ifeq ($(call match,$(FFARCH),(ia64|hp*)),yes)
include $(FEATFLOW)/Globals.x86_64
endif

ifeq ($(call match,$(FFARCH),sun),yes)
include $(FEATFLOW)/Globals.sparc
endif

ifeq ($(call match,$(FFARCH),alpha),yes)
include $(FEATFLOW)/Globals.alpha
endif

ifeq ($(call match,$(FFARCH),(power|ppc)),yes)
include $(FEATFLOW)/Globals.power
endif

########################################################################
# If no matching section could be found which matches the overwritten
# ID string then the original ID is adopted
########################################################################

ifdef ALT
ifeq ($(GENERIC),1)
ID:=$(IDORIG)
endif
endif

########################################################################
# Combine BLAS and LAPACK together into one compiler flag
########################################################################

# In case the content of BLASLIB/LAPACKLIB is the same, the user
# uses one shared library which contains both, BLAS and LAPACK routines.

BLASLAPACKLIB = $(BLASLIB) 

ifneq "$(LAPACKLIB)" "$(BLASLIB)"
BLASLAPACKLIB = $(LAPACKLIB) $(BLASLIB)
endif

########################################################################
# Set undefined variables to default values
########################################################################

ifeq "$(ARC)" ""
ARC := $(AR)
endif

# Default C++ compiler is the CC compiler

ifeq "$(CPP)" ""
CPP := $(CC)
OPTFLAGSCPP:=$(OPTFLAGSC)
CPPFLAGS:=$(CCFLAGS)
endif

########################################################################
# hack debug flags if 'make debug' is applied
########################################################################

debug: OPTFLAGS=$(OPTFLAGSDEBUG)
debug: OPTFLAGSC=$(OPTFLAGSCDEBUG)
debug: OPTFLAGSCPP=$(OPTFLAGSCPPDEBUG)
debug: OPTFLAGSF=$(OPTFLAGSFDEBUG)

########################################################################
# hack mpi flags if 'make mpi' is applied
########################################################################

# mpi: FC=$(MPIFC)

########################################################################
# hack to have this target in all Makefiles, the dot is to not
# consider it as a default rule when called without specific target
########################################################################

FCC:=$(shell (which $(CC) 2>/dev/null || echo "$(CC) not found !!"))
FCPP:=$(shell (which $(CPP) 2>/dev/null || echo "$(CPP) not found !!"))
FFC:=$(shell (which $(FC) 2>/dev/null || echo "$(FC) not found !!"))
ARF:=$(shell (which $(AR) 2>/dev/null || echo "$(AR) not found !!"))
ARC:=$(shell (which $(ARC) 2>/dev/null || echo "$(ARC) not found !!"))
.id:
	@echo
	@echo 'Machine-ID' "($(shell uname -n))" ':' $(ID) 
	@echo 
	@echo 'Compilers to be used:'
	@echo '  C compiler:        ' $(FCC)
	@echo '  C++ compiler:      ' $(FCPP)
	@echo '  Fortran compiler:  ' $(FFC)
	@echo '  F-Library archiver:' $(ARF)
	@echo '  C-Library archiver:' $(ARC)
	@echo
	@echo 'Flags to be used:'
	@echo '  OPTFLAGS        =' $(OPTFLAGS)
	@echo '  OPTFLAGSC       =' $(OPTFLAGSC)
	@echo '  OPTFLAGSCPP     =' $(OPTFLAGSCPP)
	@echo '  OPTFLAGSF       =' $(OPTFLAGSC)
	@echo '  OPTFLAGSDEBUG   =' $(OPTFLAGSDEBUG)
	@echo '  OPTFLAGSCDEBUG  =' $(OPTFLAGSCDEBUG)
	@echo '  OPTFLAGSCPPDEBUG=' $(OPTFLAGSCPPDEBUG)
	@echo '  OPTFLAGSFDEBUG  =' $(OPTFLAGSFDEBUG)
	@echo '  FCFLAGS         =' $(FCFLAGS)
	@echo '  CCFLAGS         =' $(CCFLAGS)
	@echo '  CPPFLAGS        =' $(CPPFLAGS)
	@echo '  BUILDLIB        =' $(BUILDLIB)
	@echo '  BLASLIB         =' $(if $(BLASLIB),$(BLASLIB),"(standard BLAS, included in installation package)")
ifeq "$(LAPACKLIB)" ""
	@echo '  LAPACKLIB       = (standard LAPACK, included in installation package)'
else
ifeq "$(LAPACKLIB)" "$(BLASLIB)"
			@echo '  LAPACKLIB= (Shared BLAS/LAPACK)'
else
			@echo '  LAPACKLIB= ' $(LAPACKLIB)
endif
endif
	@echo '  LDLIBS          =' $(LDLIBS)
	@echo '  LDFLAGS         =' $(LDFLAGS)
	@echo 
	@(if [ ! -x "$(FCC)" ] ; then echo 'Please edit Globals.mk to specify your C compiler' ; exit 1; fi)
	@(if [ ! -x "$(FFC)" ] ; then echo 'Please edit Globals.mk to specify your Fortran compiler' ; exit 1; fi)

.help:
	@echo 'Available global make targets:'
	@echo ' all           - compile library / application module'
	@echo ' debug         - like "all", but compile with debug infos'
	@echo ' test          - runs an application test'
	@echo ' id            - print out settings for current ID'
	@echo ' clean         - remove all for not needed for run (object files)'
	@echo ' purge         - remove all that can be removed'
	@echo 
	@echo 'Available make modifiers:'
	@echo ' ALT=xxx       - specification of alternative ID to use ID-xxx as a new ID.'
	@echo '                 (See Globals.mk for examples)'
	@echo ' ID=xxx        - overides the autodetected architecture ID by xxx'
	@echo '                 (See Globals.mk for details)'
	@echo ' FFARCH=xxx    - overwrites the autodetected value for FFARCH by xxx'
	@echo '                 (See Globals.mk for details)'
	@echo ' FFCPU=xxx     - overwrites the autodetected value for FFCPU by xxx'
	@echo '                 (See Globals.mk for details)'
	@echo ' FFCORE=xxx    - overwrites the autodetected value for FFCORE by xxx'
	@echo '                 (See Globals.mk for details)'
	@echo ' HDF5=YES      - compile and link HDF5 libraries'
	@echo '                 (See Globals.mk for details)'

