# -*- mode: makefile -*-

##############################################################################
# Intel Math Kernel Library (MKL) for Intel compiler suite
#
# providing both an implementation of the BLAS and LAPACK libraries
##############################################################################

ifneq ($(strip $(INTEL_MKL_LIB)),)
# Split up string if multiple directories are given
# Note: Do not put whitespace between comma and the environment variable, because
#       if you do, a string like "-L /path/to/mkl" is the result and that string
#       won't make it into the command line.
LIBDIR   := $(LIBDIR) -L$(subst :, -L,$(INTEL_MKL_LIB))
endif

# Settings working with older Intel MKL releases (that did provide a generic libmkl.so)
# LIBS     := $(LIBS) -lmkl -lguide -lpthread

# Settings working with Intel MKL 10.2 Update X and above
ifeq ($(strip $(OPENMP)), YES)
# Threaded implementation

# 32 bit
ifeq ($(call match,$(ID),pc-.*-linux-.*-mkl.*),yes)
LIBS := $(LIBS) -lmkl_intel -lmkl_intel_thread -lmkl_core
endif

# 64 bit
ifeq ($(call match,$(ID),pc64-.*-linux-.*-mkl.*),yes)
ifeq ($(strip $(INTSIZE)), LARGE)
LIBS := $(LIBS) -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core
else
LIBS := $(LIBS) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
endif
endif

else
# Sequential implementation

# 32 bit
ifeq ($(call match,$(ID),pc-.*-linux-.*-mkl.*),yes)
LIBS := $(LIBS) -lmkl_intel -lmkl_sequential -lmkl_core
endif

# 64 bit
ifeq ($(call match,$(ID),pc64-.*-linux-.*-mkl.*),yes)
ifeq ($(strip $(INTSIZE)), LARGE)
LIBS := $(LIBS) -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
else
LIBS := $(LIBS) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
endif
endif

endif


##############################################################################
# BLAS and LAPACK also needed by the Sparse Banded Blas benchmark
##############################################################################

# Settings working with older Intel MKL releases (that did provide a generic libmkl.so)
# SBB_LIBS     := $(SBB_LIBS) -lmkl -lguide -lpthread

# Settings working with Intel MKL 10.2 Update X and above
ifeq ($(strip $(OPENMP)), YES)
# Threaded implementation

# 32 bit
ifeq ($(call match,$(ID),pc-.*-linux-.*-mkl.*),yes)
SBB_LIBS := $(SBB_LIBS) -lmkl_intel -lmkl_intel_thread -lmkl_core
endif

# 64 bit
ifeq ($(call match,$(ID),pc64-.*-linux-.*-mkl.*),yes)
ifeq ($(strip $(INTSIZE)), LARGE)
SBB_LIBS := $(SBB_LIBS) -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core
else
SBB_LIBS := $(SBB_LIBS) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
endif
endif

else
# Sequential implementation

# 32 bit
ifeq ($(call match,$(ID),pc-.*-linux-.*-mkl.*),yes)
SBB_LIBS := $(SBB_LIBS) -lmkl_intel -lmkl_sequential -lmkl_core
endif

# 64 bit
ifeq ($(call match,$(ID),pc64-.*-linux-.*-mkl.*),yes)
ifeq ($(strip $(INTSIZE)), LARGE)
SBB_LIBS := $(SBB_LIBS) -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
else
SBB_LIBS := $(SBB_LIBS) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
endif
endif

endif


# The settings needed to compile a FEAT2 application are "wildly" distributed
# over several files ((Makefile.inc and templates/*.mk) and if-branches
# (in an attempt to reduce duplicate code and inconsistencies among all build
# IDs that e.g. use the same MPI environment). Not having all settings at
# *one* place entails the risk (especially in the event of setting up a new
# build ID) that settings are incompletely defined. A simple typo in a matching
# rule in Makefile.inc may prevent that the compiler and compiler command line
# flags are set. Compilation would fail with the most peculiar errors - if not
# the Makefile had been set up to catch such a case.
# Each build ID in FEAT2 has 6 tokens: architecture, cpu, operating system,
# compiler family, BLAS implementation, MPI environment. Whenever setting
# one of these, an according flag is set. They are named TOKEN1 up to TOKEN6.
# Before starting to actually compile a FEAT2 application, every Makefile
# generated by bin/configure checks whether all six tokens are set *for the
# choosen build ID*. If not, the user gets an error message describing exactly
# what information is missing, e.g. token 5 not set which means there is no
# information available which BLAS implementation to use and where to find the
# library.
#
# In this file, the fifth token has been set: BLAS and LAPACK implementation.
# Set the flag accordingly.
TOKEN5 := 1
