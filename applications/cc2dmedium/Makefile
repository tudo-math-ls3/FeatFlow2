#!/usr/bin/env gmake

# define the root directory of the Featflow installation
FEATFLOW=../..

# include compiler definitions
include $(FEATFLOW)/Globals.mk

# application include dirs
INCDIR= -I./src/geometry -I./src/misc -I./src/parametrization \
        -I./src/triangulation -I./src -I./include

# The f77 source files
SRC = qtrafo.f \
      condoutput.f dstrings1.f dstrings2.f dstrings3.f dstringsaux.f \
      heapsort.f invert.f\
      inipar.f parq2d.f parq2domega.f parq2dwrap.f rdparm.f\
      gmv2doutput.f initria.f inmtri.f quad2tri.f structria.f triutil2d.f

# The f90 source files 
SRC +=  cc2dmedium_callback.f90 cc2dmediumm2basic.f90 \
	cc2dmediumm2init.f90 cc2dmediumm2boundary.f90 \
	cc2dmediumm2discretisation.f90 cc2dmediumm2postprocessing.f90 \
	cc2dmediumm2stationary.f90 cc2dmedium_method2.f90\
        cc2dmedium.f90

# include the standart kernel sources
include $(FEATFLOW)/kernel/kernel.mk
SRC += $(KERNELSRC)

# path for the make where to look for application source files
vpath %.f src/geometry src/misc src/parametrization src/triangulation src

vpath %.f90 src/geometry src/misc src/parametrization src/triangulation src

vpath %.inc include \
      src/geometry src/misc src/parametrization src/triangulation src

# The name of the final executable
EXEC=cc2dmedium-$(ID)

# used libraries
FEATLIB= minisplib feat2d sysutils umfpack4 amd

# include make rules
include $(FEATFLOW)/Rules_apps.mk

# include automatic dependencies
-include Deps.mk

test: clean_exec run

