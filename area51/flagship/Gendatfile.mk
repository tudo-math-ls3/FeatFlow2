#!/usr/bin/env gmake
## -*- makefile -*-

# define the root directory of the Featflow2 installation
FEATFLOW = ../..

# include compiler definitions
include $(FEATFLOW)/Globals.mk

# Overwrite compiler flags
BLASLIB = -L$(GOTOBLAS_LIB) -lgoto -lpthread

# application include dirs
INCDIR = -I./src 

# include the standard kernel sources
include $(FEATFLOW)/kernel/kernel.mk
SRC = $(KERNELSRC)

# application source files
SRC += gendatfile.f90

# path for the make where to look for application source files
vpath %.f src
vpath %.f90 src
vpath %.inc src

# The name of the final executable
EXEC=gendatfile-$(ID)

# used libraries
FEATLIB= minisplib umfpack4 amd

default: version all

# include make rules
include $(FEATFLOW)/Rules_apps.mk

# include automatic dependencies
-include Deps.mk

version: .version
	@echo "CHARACTER(LEN=5) :: VERSION = '"`date +%y`"/"`date +%m`"'" > .version
	@echo "CHARACTER(LEN=8) :: BUILD   = '"`date +%Y%m%d`"'" >> .version
