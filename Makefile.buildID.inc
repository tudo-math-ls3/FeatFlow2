#!/usr/bin/env make
########################################################################
#                                                                      #
#         FINITE ELEMENT ANALYSIS TOOLS  2                             #
#                                                                      #
# Authors: M. Koester, M. Moeller, S. Turek                            #
#          S.Buijssen                                                  #
#                                                                      #
# Contact: Applied Mathematics, TU Dortmund                            #
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
# Maintainer: Sven H.M. Buijssen, Matthias Moeller, Michael Koester    #
# Version   : $Id: Makefile.buildID.inc,v 2.4 2009/01/13 19:33:43 buijssen Exp $
########################################################################

# Path to this FEAT2 installation main directory
ifeq ($(strip $(FEAT2BASEDIR)),)
FEAT2BASEDIR:=../..
endif

# machine and architecture info, machine ID
REALID:=$(shell $(FEAT2BASEDIR)/bin/guess_id)
REALID_SHORT := $(REALID)

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

ifeq ($(call match,$(REALID),pc-.*-cygwin_nt.*),yes)
REALID:=$(strip $(REALID))-g95-blas
endif

ifeq ($(strip $(REALID)), alpha-ev6-osf1)
REALID:=$(strip $(REALID))-cf90-dxml
endif

ifeq ($(strip $(REALID)), ibm-powerpc_power4-aix)
REALID:=$(strip $(REALID))-xlf-goto
endif

ifeq ($(strip $(REALID)), pc-athlonxp-linux)
REALID:=$(strip $(REALID))-g95-goto
endif

ifeq ($(strip $(REALID)), pc-coresolo-linux)
REALID:=$(strip $(REALID))-g95-blas
endif

ifeq ($(strip $(REALID)), pc-opteron-linux)
REALID:=$(strip $(REALID))-g95-goto
endif

ifeq ($(strip $(REALID)), pc-coreduo-linux)
REALID:=$(strip $(REALID))-intel-goto
endif

ifeq ($(strip $(REALID)), pc64-coreduo-linux)
REALID:=$(strip $(REALID))-intel-goto
endif

ifeq ($(strip $(REALID)), pc-nehalem-linux)
REALID:=$(strip $(REALID))-intel-goto
endif

ifeq ($(strip $(REALID)), pc64-nehalem-linux)
REALID:=$(strip $(REALID))-intel-goto
endif

ifeq ($(strip $(REALID)), pc-athlon64-linux)
REALID:=$(strip $(REALID))-gcc-goto
endif

ifeq ($(strip $(REALID)), pc64-athlon64x2-linux)
REALID:=$(strip $(REALID))-gcc-goto
endif

ifeq ($(strip $(REALID)), pc-turionx2-linux)
REALID:=$(strip $(REALID))-g95-goto
endif

ifeq ($(strip $(REALID)), pc64-opteron-linux)
REALID:=$(strip $(REALID))-intel-goto
endif

ifeq ($(strip $(REALID)), pc64-opteronx2-linux)
REALID:=$(strip $(REALID))-intel-goto
endif

ifeq ($(strip $(REALID)), pc-opteronx2-linux)
REALID:=$(strip $(REALID))-intel-goto
endif

ifeq ($(strip $(REALID)), pc64-penryn-linux)
REALID:=$(strip $(REALID))-intel-goto
endif

ifeq ($(strip $(REALID)), pc-pentium3-linux)
REALID:=$(strip $(REALID))-gcc-blas
endif

ifeq ($(strip $(REALID)), pc-pentiumm-linux)
REALID:=$(strip $(REALID))-gcc-goto
endif

ifeq ($(strip $(REALID)), pc-pentium4-linux)
REALID:=$(strip $(REALID))-gcc-goto
endif

ifeq ($(strip $(REALID)), pc-pentium4m-linux)
REALID:=$(strip $(REALID))-gcc-goto
endif

ifeq ($(strip $(REALID)), pc64-pentium4-linux)
REALID:=$(strip $(REALID))-intel-goto
endif

ifeq ($(strip $(REALID)), sun4u-sparcv9-sunos)
REALID:=$(strip $(REALID))-sunstudio-perf
endif

ifeq ($(strip $(REALID)), sun4v-sparcv9-sunos)
REALID:=$(strip $(REALID))-sunstudio-perf
endif

ifeq ($(strip $(REALID)), sun4u-sparcv8-sunos)
REALID:=$(strip $(REALID))-sunstudio-perf
endif

ifeq ($(strip $(REALID)), sx6-none-superux)
REALID:=$(strip $(REALID))-f90-keisan
endif

ifeq ($(strip $(REALID)), sx8-none-superux)
REALID:=$(strip $(REALID))-f90-keisan
endif

ifeq ($(strip $(REALID)), ia64-itanium2-linux)
REALID:=$(strip $(REALID))-intel-mkl
endif

ifeq ($(strip $(REALID)), ia64-itanium2x2-linux)
REALID:=$(strip $(REALID))-intel-mkl
endif

ifeq ($(strip $(REALID)), pc-coresolo-darwin)
REALID:=$(strip $(REALID))-gcc-goto
endif

ifeq ($(strip $(REALID)), pc-coreduo-darwin)
REALID:=$(strip $(REALID))-gcc-blas
endif

ifeq ($(strip $(REALID)), pc64-coreduo-darwin)
REALID:=$(strip $(REALID))-gcc-blas
endif

ifeq ($(strip $(REALID)), pc-penryn-darwin)
REALID:=$(strip $(REALID))-gcc-blas
endif

ifeq ($(strip $(REALID)), pc64-penryn-darwin)
REALID:=$(strip $(REALID))-gcc-blas
endif

ifeq ($(strip $(REALID)), pc64-nehalem-darwin)
REALID:=$(strip $(REALID))-gcc-blas
endif


# Set default target platforms for compilation of parallel applications.
# i.e. refine the machine REALID even more by specifying the MPI environment.
ifeq ($(strip $(MPI)), YES)
ifeq ($(strip $(REALID)), ibm-powerpc_power4-aix-xlf-goto)
REALID:=$(strip $(REALID))-poempi
endif

ifeq ($(strip $(REALID)), pc-pentium4-linux64-intel-goto)
REALID:=$(strip $(REALID))-lanlmpi
endif

ifeq ($(strip $(REALID)), sx6-none-superux-f90-keisan)
REALID:=$(strip $(REALID))-mpi
endif

ifeq ($(strip $(REALID)), sx8-none-superux-f90-keisan)
REALID:=$(strip $(REALID))-mpi
endif

ifeq ($(strip $(REALID)), ia64-itanium2-linux64-intel-mkl)
REALID:=$(strip $(REALID))-mpi
endif

ifeq ($(call match,$(REALID),[^-][^-]*-[^-][^-]*-[^-][^-]*-[^-][^-]*-[^-][^-]*),yes)
REALID:=$(REALID)-ompi
endif
endif


# Set default build ID if none set
ifeq ($(strip $(ID)),)
ID:=$(REALID)
endif

# Automatically fix an incompletely given build ID in case MPI=YES
ifeq ($(call match,$(ID),[^-][^-]*-[^-][^-]*-[^-][^-]*-[^-][^-]*-[^-][^-]*),yes)
ID:=$(ID)-ompi
endif




# ID distinguishes machine and architecture, machine ID, Compiler, BLAS implementation
# and MPI environment. However, libraries like umfpack, metis, lapack, blas do not refer
# to MPI environment and/or BLAS implementation. So, any combination of MPI/BLAS can use
# the same libraries. So, before extending ID, set LIBID to the first four tokens of
# ID's value.
#
# Unfortunately, the complete ID can also be given on command line, in which case we
# have to cut the last two identifiers (specifying blas library implementation and MPI
# environment). So, use current ID setting and cutting it will in any case give us what
# we need.
LIBID:=$(shell echo "$(ID)" | LC_ALL=C perl -pe 's/^([^-]+-[^-]+-[^-]+-[^-]+-[^-]+)(-[^-]+|)$$/$$1/;')

# Catch the possible error of configuring for a specific build ID and serial execution
# mode, but still giving 6 tokens, i.e. the error of unnecessarily adding the MPI token.
# Remove it.
ifeq ($(strip $(MPI)), NO)
ID:=$(shell echo "$(ID)" | LC_ALL=C perl -pe 's/^([^-]+-[^-]+-[^-]+-[^-]+-[^-]+)(-[^-]+|)$$/$$1/;')
endif


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
