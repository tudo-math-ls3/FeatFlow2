#!/usr/bin/env gmake

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

ID:=$(shell $(FEATFLOW)/bin/guess_id)
HOST:=$(shell uname -n)
match= $(shell $(FEATFLOW)/bin/match_id $(1) '$(2)')

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
APPS= trigen2d trigen3d tr2to3 cc2d pp2d cc3d pp3d 

# list of all library modules available at the top level
LIBS= feat3d feat2d sysutils umfpack2 amd umfpack4 minisplib lapack blas

########################################################################
# General name and location of compilers.
#
# These entries are machine-dependend and might be necessary to be
# modified by the user. The entries here are the general ones (!), used
# IF THE MACHINE COULD NOT BE DETECTED. If the script is able to
# detect a known machine, the machine-specific configuration branch
# (see below) is used.
#
# Featflow needs an F77-compiler as well as a C-compiler. 
########################################################################

FC=f77     # Fortran 77 compiler
CC=cc      # ANSI C compiler
LD=$(FC)   # linker (usually same as the Fortran Compiler) 
AR=ar      # library creation tool (archiver)

########################################################################
# compiler settings for various machine types, list to be extended
# make sure to include some BLAS implementation like the one provided
########################################################################

# default values for general arch 

# special compiler optimization flags:

OPTFLAGS=

# general compiler options for Fortran compiler:

FCFLAGS=

# general compiler options for C compiler:

CCFLAGS=

# list of featflow included libs to be build by the top level make
# feat2d, feat3d and sysutils required, include lapack and blas if necessary

#BUILDLIB= feat3d feat2d sysutils #lapack blas 

BUILDLIB = $(LIBS)

# BLAS library to be used, if left empty the included blas and lapack
# is used automatically; to use a single combined BLAS/LAPACK
# library (e.g. ATLAS), set BLASLIB=LAPACKLIB!

BLASLIB  =          # will use included blas
LAPACKLIB=          # will use included lapack

# additional linker options
LDFLAGS=

# default directory for .mod files in Fortran 95.
MODDIRF90=

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
########################################################################


# This specifies the C compiler (CC) and fortran compiler (FC) and
# options for the system with ID=sun4u-sparcv?-sunos (i.e
# sun4u-sparcv7-sunos, sun4u-sparcv8-sunos or sun4u-sparcv9-sunos)
# The compilers CC/FC can be specified with full directory path if needed. 

ifeq ($(call match,$(ID),sun4u-sparcv[789]-sunos),yes)
CC=cc
FC=f77
OPTFLAGS  = -fast 
FCFLAGS   = -xarch=native -silent
CCFLAGS   = -xarch=native
BLASLIB   = -xlic_lib=sunperf
LAPACKLIB = -xlic_lib=sunperf
# BLASLIB=LAPACKLIB enforces to use a single combined library for
# both, BLAS and LAPACK.
endif

# This specification will be used on the machine with
# ID=sun4u-sparcv9-sunos and if the ALT variable is set to 64bit by
# make ALT=64bit ...
# the ALT=64bit has to be specified on every call of the make command in
# order to use these settings

ifeq ($(ID),sun4u-sparcv9-sunos-64bit)
CC=cc
FC=f77
OPTFLAGS  = -fast
FCFLAGS   = -xarch=native64 -silent
CCFLAGS   = -xarch=native64
BLASLIB   = -xlic_lib=sunperf
LAPACKLIB = -xlic_lib=sunperf
# BLASLIB=LAPACKLIB enforces to use a single combined library for
# both, BLAS and LAPACK.
endif

# This will use the gcc compiler on the sun machine and will apply if
# the ALT variable is set to gcc, i.e.
# make ALT=gcc ...

ifeq ($(ID),sun4u-sparcv9-sunos-gcc)
CC=gcc
FC=g77
OPTFLAGS  = -O3 -ffast-math -fexpensive-optimizations -fomit-frame-pointer -funroll-loops -fprefetch-loop-arrays
FCFLAGS   = -mcpu=v9 -m64 -fno-globals -Wno-globals -pipe
CCFLAGS   = -mcpu=v9 -m64 -fno-globals -Wno-globals -pipe
BLASLIB   = 
LAPACKLIB =
endif

# This specifies the C compiler (CC) and fortran compiler (FC) and
# options for the systems with ID matching regexp alpha-ev[4567]-osf1
# i.e. alpha-ev5-osf1, alpha-ev6-osf1 and alpha-ev7-osf1

ifeq ($(call match,$(ID),alpha-ev[4567]-osf1),yes)
CC=cc
FC=f77
OPTFLAGS  = -fast
FCFLAGS   =
CCFLAGS   = 
BLASLIB   = -ldxml
LAPACKLIB = -ldxml
# BLASLIB=LAPACKLIB enforces to use a single combined library for
# both, BLAS and LAPACK.
endif

# This specifies the C compiler (CC) and fortran compiler (FC) and
# options for the systems with ID matching regexp
# pc-unknown-(linux|cygwin_nt?.?)  
# i.e. pc-unknown-linux or pc-unknown-cygwin_nt?.?
# with cygwin_nt?.? matching any version of cygwin_ntX.X

ifeq ($(call match,$(ID),pc-unknown-(linux|cygwin_nt?.?)),yes)
CC=gcc
FC=g77
OPTFLAGS  = -O3 -ffast-math -fexpensive-optimizations -fomit-frame-pointer -funroll-loops
FCFLAGS   = -fno-globals -Wno-globals -pipe
CCFLAGS   = -fno-globals -Wno-globals -pipe
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(call match,$(ID),pc-athlon-(linux|cygwin_nt?.?)),yes)
CC=gcc
FC=g77
OPTFLAGS  = -O3 -ffast-math -fexpensive-optimizations -fomit-frame-pointer -funroll-loops -fprefetch-loop-arrays
FCFLAGS   = -march=athlon -fno-globals -Wno-globals -pipe
CCFLAGS   = -march=athlon -fno-globals -Wno-globals -pipe
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(call match,$(ID),pc-athlonxp-(linux|cygwin_nt?.?)),yes)
CC=gcc
FC=g77
OPTFLAGS  = -O3 -ffast-math -fexpensive-optimizations -fomit-frame-pointer -funroll-loops -fprefetch-loop-arrays
FCFLAGS   = -march=athlon-xp -fno-globals -Wno-globals -pipe
CCFLAGS   = -march=athlon-xp -fno-globals -Wno-globals -pipe
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(call match,$(ID),pc-pentium3-(linux|cygwin_nt?.?)),yes)
CC=gcc
FC=g77
OPTFLAGS  = -O3 -mfpmath=sse -ffast-math -fexpensive-optimizations -fomit-frame-pointer -funroll-loops -fprefetch-loop-arrays
FCFLAGS   = -march=pentium3 -fno-globals -Wno-globals -pipe
CCFLAGS   = -march=pentium3 -fno-globals -Wno-globals -pipe
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(call match,$(ID),pc-pentiumm-(linux|cygwin_nt?.?)),yes)
CC=gcc
FC=g95
#OPTFLAGS  = -O3 -mfpmath=sse -msse2 -ffast-math -fexpensive-optimizations -fomit-frame-pointer -funroll-loops -fprefetch-loop-arrays
#FCFLAGS   = -march=pentium3 -fno-globals -Wno-globals -pipe
#CCFLAGS   = -march=pentium3 -fno-globals -Wno-globals -pipe
OPTFLAGS  = -O3 -mfpmath=sse -msse2 -ffast-math -fexpensive-optimizations -fomit-frame-pointer -funroll-loops -fprefetch-loop-arrays
FCFLAGS   = -march=pentium3 -pipe -fmod=$(MODDIRF90)
CCFLAGS   = -march=pentium3 -fno-globals -Wno-globals -pipe
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(call match,$(ID),pc-pentium4-(linux|cygwin_nt?.?)),yes)
CC=gcc
FC=g95
OPTFLAGS  = -O3 -mfpmath=sse -ffast-math -fexpensive-optimizations -fomit-frame-pointer -funroll-loops -fprefetch-loop-arrays
FCFLAGS   = -march=pentium4 -pipe -fmod=$(MODDIRF90)
CCFLAGS   = -march=pentium4 -fno-globals -Wno-globals -pipe
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),ppc64-power5-linux)
CC=gcc
FC=gfortran
OPTFLAGS  = -mcpu=power5 -mtune=power5 -O3 -maltivec -mabi=altivec -ffast-math -fexpensive-optimizations -fomit-frame-pointer -funroll-loops -fprefetch-loop-arrays
FCFLAGS   = -pipe
CCFLAGS   = -pipe
BLASLIB   = 
LAPACKLIB = 
endif

# Apple mac/osx with gcc/gfort 4.0
ifeq ($(ID),power_macintosh-ppc_7450-darwin)
CC=gcc
FC=gfortran
OPTFLAGS  = -mtune=G4 -mcpu=G4 -O3 -maltivec -mabi=altivec -ffast-math -fexpensive-optimizations -fomit-frame-pointer -funroll-loops -fprefetch-loop-arrays
FCFLAGS   = -pipe
CCFLAGS   = -pipe
LDFLAGS   = -bind_at_load
BLASLIB   = -faltivec -framework Accelerate
LAPACKLIB = -faltivec -framework Accelerate
endif

# Apple mac/osx with absoft 8 compiler
ifeq ($(ID),power_macintosh-ppc_7450-darwin-absoft)
CC=cc
FC=/Applications/Absoft/bin/f77
OPTFLAGS  = 
FCFLAGS   = -s -N1 -N109
CCFLAGS   = 
LDFLAGS   = 
BLASLIB   = -unixlib
LAPACKLIB = -unixlib
endif

# Apple mac/osx with gcc/g95 3.x
ifeq ($(ID),power_macintosh-ppc_7450-darwin-gcc3)
CC=gcc
FC=g95
OPTFLAGS  = -mtune=G4 -mcpu=G4 -O3 -maltivec -mabi=altivec -ffast-math -fexpensive-optimizations -fomit-frame-pointer -funroll-loops -fprefetch-loop-arrays
FCFLAGS   = -pipe -fmod=$(MODDIRF90)
CCFLAGS   = -pipe
LDFLAGS   = 
BLASLIB   = 
LAPACKLIB = 
endif

# Apple mac/osx with IBM compilers
ifeq ($(ID),power_macintosh-ppc_7450-darwin-xlf)
CC=/opt/ibmcmp/vac/6.0/bin/xlc
FC=/opt/ibmcmp/xlf/8.1/bin/f77
OPTFLAGS  = -O4
FCFLAGS   = 
CCFLAGS   = 
BLASLIB   = 
LAPACKLIB = 
endif


ifeq ($(ID),pc-athlonxp-linux-ifc)
CC=icc
FC=ifort
AR=xiar
OPTFLAGS  = -O3 -ipo -tpp6
FCFLAGS   = -cm -vec_report0 -fpe0 -module $(MODDIRF90)
CCFLAGS   = -cm -vec_report0 -fpe0
LDFLAGS   = -L/usr/local/ifc/lib -lifcore 
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),pc-athlonxp-linux-ifc8)
# For Intel Fortran Compiler 8, which supports inter-file optimizations
# with libraries (-ipo, -ipo_obj). Version 9.xx can't link anymore :(
CC=icc
FC=ifort
OPTFLAGS  = -O3 -ipo -ipo_obj -tpp6
FCFLAGS   = -cm -vec_report0 -fpe0  -module $(MODDIRF90)
CCFLAGS   = -cm -vec_report0 -fpe0
LDFLAGS   = -L/usr/local/ifc/lib -lifcore 
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),pc-athlonxp-linux-pgi)
CC=gcc
FC=pgf77
OPTFLAGS  = -fastsse -O4 
FCFLAGS   = -tp athlonxp
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),pc-pentium4-linux-ifc)
CC=icc
FC=ifort
AR=xiar
OPTFLAGS  = -O3 -xN -ipo
FCFLAGS   = -f77rtl -cm -fpe0 -vec_report0 -module $(MODDIRF90)
CCFLAGS   = -cm -fpe0 -vec_report0
LDFLAGS   = 
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),pc-pentium4-linux-ifc8)
# For Intel Fortran Compiler 8, which supports inter-file optimizations
# with libraries (-ipo, -ipo_obj). Version 9.xx can't link anymore :(
CC=icc
FC=ifort
OPTFLAGS  = -O3 -xN -ipo -ipo_obj
FCFLAGS   = -f77rtl -cm -fpe0 -vec_report0 -module $(MODDIRF90)
CCFLAGS   = -cm -fpe0 -vec_report0
LDFLAGS   = 
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),pc-pentiumm-linux-ifc)
# Intel Fortran Compiler 9, Intel C Compiler 9, Pentium M
CC=icc
FC=ifort
AR=xiar
OPTFLAGS  = -O3 -xB -ipo
FCFLAGS   = -f77rtl -cm -fpe0 -vec_report0 -module $(MODDIRF90)
CCFLAGS   = -cm -fpe0 -vec_report0
LDFLAGS   = 
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),pc-pentiumm-linux-ifc8)
# Intel Fortran Compiler 8, Intel C Compiler 8, Pentium M
# For Intel Fortran Compiler 8, which supports inter-file optimizations
# with libraries (-ipo, -ipo_obj). Version 9.xx can't link anymore :(
CC=icc
FC=ifort
AR=xiar
OPTFLAGS  = -O3 -xB -ipo -ipo_obj
FCFLAGS   = -f77rtl -cm -fpe0 -vec_report0 -module $(MODDIRF90)
CCFLAGS   = -cm -fpe0 -vec_report0
LDFLAGS   = 
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),pc64-opteron-linux)
CC=gcc
FC=g95
OPTFLAGS  = -O3 -m64 -mmmx -msse -msse2 -m3dnow -mfpmath=sse \
            -ffast-math -fexpensive-optimizations -ffinite-math-only \
            -fgcse -floop-optimize -foptimize-register-move -foptimize-sibling-calls -frename-registers -freorder-blocks -fomit-frame-pointer -funroll-loops -fprefetch-loop-arrays -fsched-interblock -frerun-loop-opt -frerun-cse-after-loop -freorder-functions
FCFLAGS   = -march=opteron -pipe -fmod=$(MODDIRF90)
CCFLAGS   = -march=opteron -fno-globals -Wno-globals -pipe
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),pc64-opteron-linux-path)
CC=/opt/pathscale/bin/pathcc
FC=/opt/pathscale/bin/pathf90
OPTFLAGS  = -Ofast
FCFLAGS   = -march=opteron -pipe
CCFLAGS   = -march=opteron -pipe
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),pc64-athlon64-linux)
CC=gcc
FC=g95
OPTFLAGS  = -O3 -m64 -mmmx -msse -msse2 -m3dnow -mfpmath=sse -ffast-math -fexpensive-optimizations -ffinite-math-only -fgcse -floop-optimize -fmove-all-movables -foptimize-register-move -foptimize-sibling-calls -frename-registers -freorder-blocks -fomit-frame-pointer -funroll-loops -fprefetch-loop-arrays -fsched-interblock -frerun-loop-opt -frerun-cse-after-loop -freorder-functions
FCFLAGS   = -march=athlon64 -pipe -fmod=$(MODDIRF90)
CCFLAGS   = -march=athlon64 -fno-globals -Wno-globals -pipe
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(call match,$(ID),pc-opteron-(linux|cygwin_nt?.?)),yes)
CC=gcc
FC=g95
OPTFLAGS  = -O3 -m32 -mmmx -msse -msse2 -m3dnow -mfpmath=sse -ffast-math -fexpensive-optimizations -ffinite-math-only -fgcse -floop-optimize -fmove-all-movables -foptimize-register-move -foptimize-sibling-calls -frename-registers -freorder-blocks -fomit-frame-pointer -funroll-loops -fprefetch-loop-arrays -fsched-interblock -frerun-loop-opt -frerun-cse-after-loop -freorder-functions
FCFLAGS   = -march=opteron -pipe -fmod=$(MODDIRF90)
CCFLAGS   = -march=opteron -fno-globals -Wno-globals -pipe
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),pc64-opteron-linux-gcc4)
CC=gcc
FC=gfortran
OPTFLAGS  = -O3 -ffast-math -fexpensive-optimizations
FCFLAGS   = -march=opteron -pipe
CCFLAGS   = -march=opteron -pipe
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),pc64-opteron-linux-ifc)
CC=gcc
FC=ifort
AR=xiar
OPTFLAGS  = -O3 -ipo -msse2 -mtune=pentiumpro -march=pentium4
FCFLAGS   = -f77rtl -cm -vec_report0 -module $(MODDIRF90)
CCFLAGS   = -vec_report0
LDFLAGS   = 
BLASLIB   = 
LAPACKLIB = 
endif

# IFORT compiler, GOTO BLAS, large arrays > 2GB
# High precision: 128 bit double, 64 bit integer

ifeq ($(ID),pc64-opteron-linux-ifclargegoto)
CC=/usr/local/icce/bin/icc
FC=/usr/local/ifce/bin/ifort
OPTFLAGS  = -O2 -ip -no-prec-div
FCFLAGS   = -mcmodel=medium -module $(MODDIRF90)
CCFLAGS   = -mcmodel=medium 
LDFLAGS   = -i-dynamic -threads
BLASLIB   = ./libgoto_opteron64p-r1.00.so 
LAPACKLIB = 
endif

# IFORT compiler, standard BLAS, large arrays > 2GB
# High precision: 128 bit double, 64 bit integer

ifeq ($(ID),pc64-opteron-linux-ifclarge)
CC=/usr/local/icce/bin/icc
FC=/usr/local/ifce/bin/ifort
OPTFLAGS  = 
FCFLAGS   = -mcmodel=large -integer_size 64 -double_size 128 -real_size 64 -module $(MODDIRF90)
CCFLAGS   = -mcmodel=large -integer_size 64 -double_size 128 -real_size 64 
BLASLIB   = 
LAPACKLIB = 
endif

# GCC compiler, standard BLAS, large arrays > 2GB

ifeq ($(ID),pc64-opteron-linux-gcclarge)
CC=g++
FC=g95
OPTFLAGS  = 
FCFLAGS   = -mcmodel=medium -m64 -fmod=$(MODDIRF90)
CCFLAGS   = -mcmodel=medium -m64
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),pc64-opteron-linux-pgi)
CC=pgcc
FC=pgf77
OPTFLAGS  = -fastsse -O4 -tp k8-64 -Mipa -mcmodel=medium -Mlarge_arrays 
FCFLAGS   = 
CCFLAGS   = 
LDFLAGS   = -lpgftnrtl -lm
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),ia64-itanium2-linux)
CC=cc
FC=f77
OPTFLAGS  = -O3 -ffast-math -fexpensive-optimizations -fomit-frame-pointer -funroll-loops -fprefetch-loop-arrays
FCFLAGS   = -pipe
CCFLAGS   = -pipe
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),ia64-itanium2-linux-ifc)
CC=icc
FC=ifort
AR=xiar
OPTFLAGS  = -fast
FCFLAGS   = -module $(MODDIRF90)
CCFLAGS   = 
BLASLIB   = 
LAPACKLIB = 
endif

ifeq ($(ID),hpia64-hppa-hpux)
CC=cc
FC=f90
OPTFLAGS  = +O4 +Ofaster +U77
FCFLAGS   = 
CCFLAGS   = 
BLASLIB   = -llapack
LAPACKLIB = -llapack
endif

ifeq ($(ID),hp9000800-hppa-hpux)
CC=cc
FC=f90
OPTFLAGS  = +O4 +Ofaster +U77
FCFLAGS   = 
CCFLAGS   = 
BLASLIB   = -llapack
LAPACKLIB = -llapack
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
# hacked debug cflags
########################################################################

debug: OPTFLAGS= -g

########################################################################
# hack to have this target in all Makefiles, the dot is to not
# consider it as a default rule when called without specific target
########################################################################

FCC:=$(shell (which $(CC) 2>/dev/null || echo "$(CC) not found !!"))
FFC:=$(shell (which $(FC) 2>/dev/null || echo "$(FC) not found !!"))
ARC:=$(shell (which $(AR) 2>/dev/null || echo "$(AR) not found !!"))
.id:
	@echo
	@echo 'Machine-ID' "($(shell uname -n))" ':' $(ID) 
	@echo 
	@echo 'Compilers to be used:'
	@echo '  C compiler:      ' $(FCC)
	@echo '  Fortran compiler:' $(FFC)
	@echo '  Library archiver:' $(ARC)
	@echo
	@echo 'Flags to be used:'
	@echo '  OPTFLAGS =' $(OPTFLAGS)
	@echo '  FCFLAGS  =' $(FCFLAGS)
	@echo '  CCFLAGS  =' $(CCFLAGS)
	@echo '  BUILDLIB =' $(BUILDLIB)
	@echo '  BLASLIB  =' $(if $(BLASLIB),$(BLASLIB),"(standard BLAS, included in installation package)")
ifeq "$(LAPACKLIB)" ""
	@echo '  LAPACKLIB= (standard LAPACK, included in installation package)'
else
ifeq "$(LAPACKLIB)" "$(BLASLIB)"
			@echo '  LAPACKLIB= (Shared BLAS/LAPACK)'
else
			@echo '  LAPACKLIB= ' $(LAPACKLIB)
endif
endif
	@echo '  LDLIBS   =' $(LDLIBS)
	@echo '  LDFLAGS  =' $(LDFLAGS)
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

