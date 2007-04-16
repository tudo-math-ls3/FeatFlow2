#!/usr/bin/env gmake

# define the root directory of the Featflow installation
FEATFLOW=../..

# include compiler definitions
include $(FEATFLOW)/Globals.mk

# application include dirs
INCDIR= -I./src/geometry -I./src/misc -I./src/parametrization \
        -I./src/triangulation -I./src -I./include

# The source files 
SRC =   poisson_callback.f90 poisson_method1.f90 \
        poisson_method2.f90 poisson_method3.f90 poisson_method4.f90 \
        poisson_method5.f90 poisson_method6.f90 poisson_method7.f90 \
        poisson.f90

# include the standart kernel sources
include $(FEATFLOW)/kernel/kernel.mk
SRC += $(KERNELSRC)

# path for the make where to look for application source files
vpath %.f src/geometry src/misc src/parametrization src/triangulation src

vpath %.f90 src/geometry src/misc src/parametrization src/triangulation src

vpath %.inc include \
      src/geometry src/misc src/parametrization src/triangulation src

# The name of the final executable
EXEC=poisson-$(ID)

# used libraries
FEATLIB= minisplib sysutils umfpack4 amd

# include make rules
include $(FEATFLOW)/Rules_apps.mk

# include automatic dependencies
-include Deps.mk

test: clean_exec run

