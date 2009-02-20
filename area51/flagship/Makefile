#!/usr/bin/env gmake
## -*- makefile -*-

# define the root directory of the Featflow2 installation
FEATFLOW = ../..

# include compiler definitions
include $(FEATFLOW)/Globals.mk

# Overwrite compiler flags
BLASLIB = -L$(GOTOBLAS_LIB) -lgoto -lpthread

# application include dirs
INCDIR = -I./src -I./src/kernel

# include the standard kernel sources
include $(FEATFLOW)/kernel/kernel.mk
SRC = $(KERNELSRC)

# application kernel source files
SRC += boundaryfilter.f90 \
       problem.f90 \
       solver.f90 \
       solverlinear.f90 \
       solvernonlinear.f90 \
       thermodynamics.f90 \
       timestep.f90

# common application source files
SRC += flagship.f90 \
       flagship_basic.f90 \
       flagship_callback.f90

# application source files for CODIRE module
SRC += codire_application.f90 \
       codire_basic.f90 \
       codire_callback.f90

# application source files for EULER module
SRC += euler_application.f90 \
       euler_basic.f90 \
       euler_callback.f90 \
       euler_callback1d.f90 \
       euler_callback2d.f90 \
       euler_callback3d.f90

# application source files for MHD module
SRC += mhd_application.f90 \
       mhd_basic.f90 \
       mhd_callback.f90

# path for the make where to look for application source files
vpath %.f src src/kernel
vpath %.f90 src src/kernel
vpath %.inc src src/kernel

# The name of the final executable
EXEC=flagship-$(ID)

# used libraries
FEATLIB= minisplib umfpack4 amd

default: version all

# include make rules
include $(FEATFLOW)/Rules_apps.mk

# include automatic dependencies
-include Deps.mk

version: .version
	@echo "CHARACTER(LEN=5) :: VERSION = '"`date +%y`"/"`date +%m`"'" > .version
	@echo "CHARACTER(LEN=8) :: BUILD   = '"`date +%Y%m%d`"'" >> .version
