#!/usr/bin/env gmake

# define the root directory of the Featflow installation
FEATFLOW=../..

# include compiler definitions
include $(FEATFLOW)/Globals.mk

# application include dirs
INCDIR= -I./src/geometry -I./src/misc -I./src/parametrization \
        -I./src/triangulation -I./src -I./include

SRC = adaptivetimestep.f90 cc2dmediumm2basic.f90 \
  cc2dmediumm2boundary.f90 cc2dmediumm2boundarydef.f90 \
  cc2dmediumm2init.f90     cc2dmediumm2scriptfile.f90 \
  cc2dmediumoptc.f90       cc2dmedium_callback.f90 \
  cc2dmedium_method2.f90   spacetimelinearsolver.f90 \
  nonlinearoneshotspacetimesolver.f90 nonstationaryoptcsolver.f90\
  optcanalysis.f90         postprocessing.f90 \
  spacediscretisation.f90  spacematvecassembly.f90\
  spacepreconditioner.f90  spacepreconditionerinit.f90\
  spacetimediscretisation.f90 spacetimeinterlevelprojection.f90\
  spacetimelinearsystem.f90   spacetimevanca.f90 \
  spacetimevectors.f90        stationaryoptcsolver.f90 \
  timeanalysis.f90            timeboundaryconditions.f90 \
  timediscretisation.f90      timeevaluation.f90 \
  timerhsevaluation.f90

# include the standart kernel sources
include $(FEATFLOW)/kernel/kernel.mk
SRC += $(KERNELSRC)

# path for the make where to look for application source files
vpath %.f src/geometry src/misc src/parametrization src/triangulation src

vpath %.f90 src/geometry src/misc src/parametrization src/triangulation src

vpath %.inc include \
      src/geometry src/misc src/parametrization src/triangulation src

# The name of the final executable
EXEC=cc2dmediumoptc-$(ID)

# used libraries
FEATLIB= minisplib feat2d sysutils umfpack4 amd

# include make rules
include $(FEATFLOW)/Rules_apps.mk

# include automatic dependencies
-include Deps.mk

test: clean_exec run

