#!/usr/bin/env make
########################################################################
#                                                                      #
#                   FINITE ELEMENT ANALYSIS TOOLS 2                    #
#                                                                      #
# Authors: M. Koester, M. Moeller, S. Turek, S. Buijssen               #
#                                                                      #
#                                                                      #
# Contact: Applied Mathematics, TU Dortmund University                 #
#          Vogelpothsweg 87, 44227 Dortmund                            #
#          Germany                                                     #
#                                                                      #
# Web:     http://www.featflow.de/en/software/featflow2.html           #
#          mailto:featflow@featflow.de                                 #
#                                                                      #
########################################################################
#                                                                      #
# Determine cpu, architecture and operating system and set a default   #
# build ID. Is to be included in every Makefile for FEAT2 applications #
#                                                                      #
# Author    : Sven H.M. Buijssen                                       #
# Maintainer: Sven H.M. Buijssen, Matthias Moeller, Michael Koester    #
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
# Refine the machine ID by specifying main Fortran compiler and blas library
# implementation to come up with default build IDs for every platform. 
# List is to be extended for every new host FEAT is ported to.
##############################################################################

# Use DEC Fortran compiler for DEC Alpha machines
ifeq ($(call match,$(REALID),alpha-.*-osf1.*),yes)
  REALID:=$(strip $(REALID))-cf90-dxml
endif

# Use xlf Fortran compiler on IBM machines running AIX
ifeq ($(call match,$(REALID),ibm-.*-aix.*),yes)
  REALID:=$(strip $(REALID))-xlf-blas
endif

# Use SunStudio Fortran compiler on Sun machines running SunOS
ifeq ($(call match,$(REALID),sun.*-.*-sunos.*),yes)
  REALID:=$(strip $(REALID))-sunstudio-perf
endif

# Use NEC Fortran compiler on NEC SX machines running SuperUX
ifeq ($(call match,$(REALID),sx.*-.*-superux.*),yes)
  REALID:=$(strip $(REALID))-f90-keisan
endif


# Default build ID extenstion if no match in the previous definitions
#
# part 1:
# Use Intel compiler if a related environment variable exists, use GCC
# compiler if not. With Intel, try to use MKL.
ifeq ($(call match,$(REALID),[^-]*-[^-]*-[^-]*),yes)
  ifneq ($(strip $(INTEL_LICENSE_FILE)$(INTEL_CC_HOME)$(INTEL_FC_HOME)),)
    REALID:=$(strip $(REALID))-intel
    ifneq ($(strip $(INTEL_MKL_LIB)$(MKLROOT)),)
      REALID:=$(strip $(REALID))-mkl
    endif
  else
    REALID:=$(strip $(REALID))-gcc
  endif
endif

# part 2:
# If no BLAS/LAPACK installation predescribed yet in build ID, try OpenBLAS,
# GotoBLAS2, then non-optimised, self-compiled local BLAS.
ifeq ($(call match,$(REALID),[^-]*-[^-]*-[^-]*-[^-]*),yes)
  ifneq ($(strip $(OPENBLAS_HOME)$(OPENBLAS_LIB)),)
    REALID:=$(strip $(REALID))-openblas
  else
    ifneq ($(strip $(GOTOBLAS2_HOME)$(GOTOBLAS2_LIB)),)
      REALID:=$(strip $(REALID))-goto2
    else
      REALID:=$(strip $(REALID))-blas
    endif
  endif
endif


# Default build id extenstion if no match in the previous definitions
# gcc compiler with non-optimised BLAS, only standard optimisation.
ifeq ($(call match,$(REALID),[^-]*-[^-]*-[^-]*),yes)
  REALID:=$(strip $(REALID))-gcc-blas
endif


# Set default target platforms for compilation of parallel applications.
# i.e. refine the machine REALID even more by specifying the MPI environment.
ifeq ($(strip $(MPI)), YES)
  ifeq ($(call match,$(REALID),ibm-.*-aix-xlf-goto),yes)
    REALID:=$(strip $(REALID))-poempi
  endif

  ifeq ($(call match,$(REALID),sx.*-.*-superux-f90-keisan),yes)
    REALID:=$(strip $(REALID))-mpi
  endif

  # If no MPI environment predescribed yet in build ID,
  # refine the machine REALID even more by specifying the MPI environment.
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

  # Automatically fix an incompletely given build ID in case MPI=YES
  # to use Open MPI.
  ifeq ($(call match,$(REALID),[^-][^-]*-[^-][^-]*-[^-][^-]*-[^-][^-]*-[^-][^-]*),yes)
    REALID:=$(REALID)-ompi
  endif

  # Automatically fix an incompletely given build ID in case MPI=YES
  # to use Open MPI.
  ifeq ($(call match,$(ID),[^-][^-]*-[^-][^-]*-[^-][^-]*-[^-][^-]*-[^-][^-]*),yes)
    ID:=$(ID)-ompi
  endif
endif   # end MPI=YES


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
