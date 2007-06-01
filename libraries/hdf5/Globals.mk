#!/usr/bin/env gmake

########################################################################
# This package is a little bit tricky to compile. For some compilers,
# it needs special options to compile properly. These options should
# not be included in the global Global.mk since they are only required
# for the HDF5 package.
########################################################################

########################################################################
# Special compiler flags for the C compiler
########################################################################

FFLAGS=

########################################################################
# Special compiler flags for the F90 compiler
########################################################################

CFLAGS=

########################################################################
# Compiler dependent flags
########################################################################

ifeq ($(FC),g95)
FFLAGS=-Wno-globals -fno-second-underscore
endif

