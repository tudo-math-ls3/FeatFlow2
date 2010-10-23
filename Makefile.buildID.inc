#!/usr/bin/env make
########################################################################
#                                                                      #
#         FINITE ELEMENT ANALYSIS & SOLUTION TOOLS  F E A S T          #
#                                                                      #
# Authors: Ch.Becker,                                                  #
#          S.Buijssen, D.Goeddeke, M.Grajewski, H.Wobker,              #
#          S.Kilian, S.Turek                                           #
#                                                                      #
# Contact: Applied Mathematics, Dortmund University of Technology      #
#          Vogelpothsweg 87, 44227 Dortmund                            #
#          Germany                                                     #
#                                                                      #
# Web:     http://www.feast.tu-dortmund.de                             #
#          mailto:feast@math.tu-dortmund.de                            #
#                                                                      #
########################################################################
#                                                                      #
# Determine cpu, architecture and operating system and set a default   #
# build ID. Is to be included in every Makefile for FEAT2 applications #
#                                                                      #
# Author    : Sven H.M. Buijssen                                       #
# Maintainer: Chr. Becker, Sven H.M. Buijssen                          #
# Version   : $Id: Makefile.buildID.inc,v 2.4 2009/01/13 19:33:43 buijssen Exp $
########################################################################

# Path to this FEAT2 installation main directory
ifeq ($(strip $(FEAT2BASEDIR)),)
FEAT2BASEDIR:=../..
endif

# machine and architecture info, machine ID
REALID:=$(shell $(FEAT2BASEDIR)/bin/guess_id)
ifeq ($(strip $(ID)),)
ID:=$(REALID)
endif


# The function match is used to match ID agains wildcards.
#
# variant 1 (slower, especially on old machines, but possibly more portable)
#match=$(shell echo $(1) | awk '/$(2)/ { print "yes"; }')
#
# variant 2 (fastest so far)
match=$(shell echo $(1) | sed -e 's/$(2)/yes/')



##############################################################################
# clean up environment
#
##############################################################################

# Don't let the script get confused by non-english messages
# from system information programs.
# (LC_ALL overrides the value of the LANG environment variable
# and the values of any other LC_* environment variables.)
LC_ALL=C

# Unset CDPATH to prevent problems when changing directories
CDPATH=



##############################################################################
# Set default target platforms for compilation of serial applications.
#
# Refine the machine ID by specifying main fortran compiler and blas library
# implementation to come up with default build IDs for every platform. 
# List is to be extended for every new host FEAT is ported to.
##############################################################################

ifeq ($(call match,$(ID),pc-.*-cygwin_nt.*),yes)
ID:=$(strip $(ID))-g95-blas
endif

ifeq ($(strip $(ID)), alpha-ev6-osf1)
ID:=$(strip $(ID))-cf90-dxml
endif

ifeq ($(strip $(ID)), ibm-powerpc_power4-aix)
ID:=$(strip $(ID))-xlf-goto
endif

ifeq ($(strip $(ID)), pc-athlonxp-linux)
ID:=$(strip $(ID))-g95-goto
endif

ifeq ($(strip $(ID)), pc-coresolo-linux)
ID:=$(strip $(ID))-g95-blas
endif

ifeq ($(strip $(ID)), pc-opteron-linux)
ID:=$(strip $(ID))-g95-goto
endif

ifeq ($(strip $(ID)), pc-coreduo-linux)
ID:=$(strip $(ID))-intel-goto
endif

ifeq ($(strip $(ID)), pc64-coreduo-linux)
ID:=$(strip $(ID))-intel-goto
endif

ifeq ($(strip $(ID)), pc-nehalem-linux)
ID:=$(strip $(ID))-intel-goto
endif

ifeq ($(strip $(ID)), pc64-nehalem-linux)
ID:=$(strip $(ID))-intel-goto
endif

ifeq ($(strip $(ID)), pc-athlon64-linux)
ID:=$(strip $(ID))-gcc-goto
endif

ifeq ($(strip $(ID)), pc-turionx2-linux)
ID:=$(strip $(ID))-g95-goto
endif

ifeq ($(strip $(ID)), pc64-opteron-linux)
ID:=$(strip $(ID))-intel-goto
endif

ifeq ($(strip $(ID)), pc64-opteronx2-linux)
ID:=$(strip $(ID))-intel-goto
endif

ifeq ($(strip $(ID)), pc-opteronx2-linux)
ID:=$(strip $(ID))-intel-goto
endif

ifeq ($(strip $(ID)), pc64-penryn-linux)
ID:=$(strip $(ID))-intel-goto
endif

ifeq ($(strip $(ID)), pc-pentium3-linux)
ID:=$(strip $(ID))-gcc-blas
endif

ifeq ($(strip $(ID)), pc-pentiumm-linux)
ID:=$(strip $(ID))-gcc-goto
endif

ifeq ($(strip $(ID)), pc-pentium4-linux)
ID:=$(strip $(ID))-gcc-goto
endif

ifeq ($(strip $(ID)), pc-pentium4m-linux)
ID:=$(strip $(ID))-gcc-goto
endif

ifeq ($(strip $(ID)), pc64-pentium4-linux)
ID:=$(strip $(ID))-intel-goto
endif

ifeq ($(strip $(ID)), sun4u-sparcv9-sunos)
ID:=$(strip $(ID))-sunstudio-perf
endif

ifeq ($(strip $(ID)), sun4v-sparcv9-sunos)
ID:=$(strip $(ID))-sunstudio-perf
endif

ifeq ($(strip $(ID)), sun4u-sparcv8-sunos)
ID:=$(strip $(ID))-sunstudio-perf
endif

ifeq ($(strip $(ID)), sx6-none-superux)
ID:=$(strip $(ID))-f90-keisan
endif

ifeq ($(strip $(ID)), sx8-none-superux)
ID:=$(strip $(ID))-f90-keisan
endif

ifeq ($(strip $(ID)), ia64-itanium2-linux)
ID:=$(strip $(ID))-intel-mkl
endif

ifeq ($(strip $(ID)), ia64-itanium2x2-linux)
ID:=$(strip $(ID))-intel-mkl
endif

ifeq ($(strip $(ID)), pc-coresolo-darwin)
ID:=$(strip $(ID))-gcc-goto
endif

ifeq ($(strip $(ID)), pc-coreduo-darwin)
ID:=$(strip $(ID))-gcc-blas
endif

ifeq ($(strip $(ID)), pc64-coreduo-darwin)
ID:=$(strip $(ID))-gcc-blas
endif

ifeq ($(strip $(ID)), pc-penryn-darwin)
ID:=$(strip $(ID))-gcc-blas
endif

ifeq ($(strip $(ID)), pc64-penryn-darwin)
ID:=$(strip $(ID))-gcc-blas
endif

ifeq ($(strip $(ID)), pc64-nehalem-darwin)
ID:=$(strip $(ID))-gcc-blas
endif


# Set default target platforms for compilation of parallel applications.
# i.e. refine the machine ID even more by specifying the MPI environment.
ifeq ($(strip $(MPI)), YES)
ifeq ($(strip $(ID)), ibm-powerpc_power4-aix-xlf-goto)
ID:=$(strip $(ID))-poempi
endif

ifeq ($(strip $(ID)), pc-pentium4-linux64-intel-goto)
ID:=$(strip $(ID))-lanlmpi
endif

ifeq ($(strip $(ID)), sx6-none-superux-f90-keisan)
ID:=$(strip $(ID))-mpi
endif

ifeq ($(strip $(ID)), sx8-none-superux-f90-keisan)
ID:=$(strip $(ID))-mpi
endif

ifeq ($(strip $(ID)), ia64-itanium2-linux64-intel-mkl)
ID:=$(strip $(ID))-mpi
endif

ifeq ($(call match,$(ID),[^-][^-]*-[^-][^-]*-[^-][^-]*-[^-][^-]*-[^-][^-]*),yes)
ID:=$(ID)-ompi
endif
endif


# ID distinguishes machine and architecture, machine ID, Compiler, BLAS implementation
# and MPI environment. However, libraries like umfpack, metis, lapack, blas do not
# refer to MPI environment and/or BLAS implementation. So, any combination of MPI/BLAS
# can use the same libraries. So, before extending ID, set LIBID to the first four tokens
# of ID's value.
#
# Unfortunately, the complete ID can also be given on command line, in which
# case we have to cut the last two identifiers (specifying blas library implementation and
# MPI environment). So, use current ID setting and cutting it will in any case give us 
# what we need.
LIBID:=$(shell echo "$(ID)" | LC_ALL=C perl -pe 's/^([^-]+-[^-]+-[^-]+-[^-]+-[^-]+)(-[^-]+|)$$/$$1/;')


# hack to have this target in all Makefiles, the dot is to not
# consider it as a default rule when called without specific target
.PHONY: .idonly
.idonly:
	@echo $(ID)

.PHONY: .libidonly
.libidonly:
	@echo $(LIBID)

.PHONY: .mpi
.mpi:
	@echo $(MPI)
