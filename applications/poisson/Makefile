#!/usr/bin/env gmake

FEATFLOW=../..

include $(FEATFLOW)/Globals.mk

INCDIR= -I./src/geometry -I./src/misc -I./src/parametrization \
        -I./src/triangulation -I./src -I./include

# The f77 source files

SRC = qtrafo.f \
      condoutput.f dstrings1.f dstrings2.f dstrings3.f dstringsaux.f \
      heapsort.f invert.f\
      inipar.f parq2d.f parq2domega.f parq2dwrap.f rdparm.f\
      gmv2doutput.f initria.f inmtri.f quad2tri.f structria.f triutil2d.f

# The F90 source files

SRCF90=
include $(FEATFLOW)/kernel/kernel.mk
SRCF90:=$(SRCF90) \
        poisson_callback.f90 poisson_method1.f90 \
        poisson_method2.f90 poisson_method3.f90 poisson.f90

# put all files together

# path for the make where to look for which files

vpath %.f \
      src/geometry src/misc src/parametrization src/triangulation src

vpath %.f90 \
      src/geometry src/misc src/parametrization src/triangulation src

vpath %.inc \
      include \
      src/geometry src/misc src/parametrization src/triangulation src

EXEC=poisson-$(ID)

FEATLIB= feat2d sysutils #umfpack4 amd

include $(FEATFLOW)/Rules_apps_f90.mk
-include Deps.mk

test: clean_exec run

