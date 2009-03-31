#!/usr/bin/env gmake

# define the root directory of the Featflow installation
FEATFLOW=../..

# include compiler definitions
include $(FEATFLOW)/Globals.mk

# application include dirs
INCDIR= -I./src/geometry -I./src/misc -I./src/parametrization \
        -I./src/triangulation -I./src -I./include \
        -I$(MPI_HOME)/include

# The source files 
SRC =   elemdbg1d_callback.f90 elemdbg1d_test1.f90 \
        elemdbg2d_callback.f90 elemdbg2d_test1.f90 \
        elemdbg3d_callback.f90 elemdbg3d_test1.f90 \
        disto2d_aux.f90 disto3d_aux.f90 \
        elemdbg.f90

# include the standart kernel sources
include $(FEATFLOW)/kernel/kernel.mk
SRC += $(KERNELSRC)

# path for the make where to look for application source files
vpath %.f src/geometry src/misc src/parametrization src/triangulation src

vpath %.f90 src/geometry src/misc src/parametrization src/triangulation src

vpath %.inc include \
      src/geometry src/misc src/parametrization src/triangulation src

vpath %.h $(MPI_HOME)/include

# The name of the final executable
EXEC=elemdbg-$(ID)

# used libraries
FEATLIB= minisplib sysutils umfpack4 amd

# include make rules
include $(FEATFLOW)/Rules_apps.mk

# include automatic dependencies
-include Deps.mk

test: clean_exec run

