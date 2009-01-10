#!/usr/bin/env gmake

# define the root directory of the Featflow installation
FEATFLOW=../..

# include compiler definitions
include $(FEATFLOW)/Globals.mk

# application include dirs
INCDIR= -I./src/geometry -I./src/misc -I./src/parametrization \
        -I./src/triangulation -I./src -I./include

# The f90 source files 
SRC =  stokes2d_callback.f90 stokes2d_method0_simple.f90 stokes2d_method1_mg.f90 \
       stokes2d_method2_sv.f90 stokes2d_method2_gv.f90 navst2d_method1_mg.f90 \
       stokes3d_callback.f90 stokes3d_method0_simple.f90 stokes3d_method1_mg.f90 \
       dom3d_cube.f90 dom3d_c3d0.f90 navst3d_method1_mg.f90 stokes.f90

# include the standart kernel sources
include $(FEATFLOW)/kernel/kernel.mk
SRC += $(KERNELSRC)

# path for the make where to look for application source files
vpath %.f src/geometry src/misc src/parametrization src/triangulation src

vpath %.f90 src/geometry src/misc src/parametrization src/triangulation src

vpath %.inc include \
      src/geometry src/misc src/parametrization src/triangulation src

# The name of the final executable
EXEC=stokes-$(ID)

# used libraries
FEATLIB= minisplib sysutils umfpack4 amd

# include make rules
include $(FEATFLOW)/Rules_apps.mk

# include automatic dependencies
-include Deps.mk

test: clean_exec run

