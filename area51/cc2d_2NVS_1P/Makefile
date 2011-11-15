#!/usr/bin/env gmake

# define the root directory of the Featflow installation
FEATFLOW=../..

# include compiler definitions
include $(FEATFLOW)/Globals.mk

# application include dirs
INCDIR= -I./src/geometry -I./src/misc -I./src/parametrization \
        -I./src/triangulation -I./src -I./include

# The f90 source files 
SRC = adaptivetimestep.f90 \
      cc2d.f90 \
      ccbasic.f90 \
      ccboundarycondition.f90 \
      ccboundaryconditionparser.f90 \
      cccallback.f90 \
      ccgeneraldiscretisation.f90 \
      ccinitgeneralparameters.f90 \
      ccinitparamtriang.f90 \
      ccmainproblem.f90 \
      ccmatvecassembly.f90 \
      ccnonlinearcore.f90 \
      ccnonlinearcoreinit.f90 \
      ccnonstationary.f90 \
      ccpostprocessing.f90 \
      ccstationary.f90 \
      cctimeanalysis.f90


# include the standart kernel sources
include $(FEATFLOW)/kernel/kernel.mk
SRC += $(KERNELSRC)

# path for the make where to look for application source files
vpath %.f src/geometry src/misc src/parametrization src/triangulation src

vpath %.f90 src/geometry src/misc src/parametrization src/triangulation src

vpath %.inc include \
      src/geometry src/misc src/parametrization src/triangulation src

# The name of the final executable
EXEC=cc2d-$(ID)

# used libraries
FEATLIB= minisplib feat2d sysutils umfpack4 amd

# include make rules
include $(FEATFLOW)/Rules_apps.mk

# include automatic dependencies
-include Deps.mk

test: clean_exec run

