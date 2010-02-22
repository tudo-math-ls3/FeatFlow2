 
#!/usr/bin/env gmake

FEATFLOW=../..

include $(FEATFLOW)/Globals.mk

LIBNAME=inshape2d

INCDIR=-I./include -I./src -I../opennurbs/include 

SRCLIST=        aabb2.cpp                distops.cpp              nurbs.cpp\
		aabb3.cpp                dynamicarray.cpp         obb2.cpp\
		aabbtree.cpp             exception.cpp            ObjLoader.cpp\
		approxcurve.cpp          paramline.cpp\
		feat2ops.cpp             primops.cpp\
		basiccurve.cpp           fileParser.cpp\
		geomops.cpp              quaternion.cpp\
		ray3.cpp		 interface.cpp\
		bctree.cpp               rootfinder.cpp\
		graham.cpp               sbb.cpp\
		sbbtree.cpp\
		bspline.cpp              grid.cpp\
		torus.cpp\
		Triangle3.cpp\
		matrix2x2.cpp            vector2.cpp\
		matrix3x3.cpp            vector3.cpp\
		matrix4x4.cpp\
		colldetect.cpp
		


# OPENNURBSFLAGS is set in globals.* sccording to the svtusl compiler.
# E.g. if a GCC compiler is involved, this is set to -D_GNU_SOURCE.
# For other compilers, a different definition may be specified.
DEFS= $(OPENNURBSFLAGS)


include $(FEATFLOW)/Rules_libs.mk
