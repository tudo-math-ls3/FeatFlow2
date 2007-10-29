#!/usr/bin/env gmake

# define the root directory of the Featflow installation
FEATFLOW=../..

# include compiler definitions
include $(FEATFLOW)/Globals.mk

# application include dirs
INCDIR= -I./src/geometry -I./src/misc -I./src/parametrization \
        -I./src/triangulation -I./src -I./include

# The f90 source files
SRC +=  heatcond_callback.f90 \
	heatcond_basic.f90 heatcond_boundarycondition.f90 \
	heatcond_partridiscr.f90 heatcond_solver.f90 \
	heatcond_matvec.f90 heatcond_timeloop.f90 \
	heatcond_method1.f90 heatcond_method2.f90 heatcond_method5.f90 \
	heatcond.f90

# include the standart kernel sources
include $(FEATFLOW)/kernel/kernel.mk
SRC += $(KERNELSRC)

# path for the make where to look for application source files
vpath %.f src/geometry src/misc src/parametrization src/triangulation src

vpath %.f90 src/geometry src/misc src/parametrization src/triangulation src

vpath %.inc include \
      src/geometry src/misc src/parametrization src/triangulation src

# The name of the final executable
EXEC=heatcond-$(ID)

# used libraries
FEATLIB= minisplib feat2d sysutils umfpack4 amd

# include make rules
include $(FEATFLOW)/Rules_apps.mk

# include automatic dependencies
-include Deps.mk

test: clean_exec run

