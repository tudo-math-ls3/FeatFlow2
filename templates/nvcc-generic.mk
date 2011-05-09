# -*- mode: makefile -*-

##############################################################################
# NVIDIA Cuda compiler 3.x
#
##############################################################################
COMPILERNAME = NVCC

# Default: Inherit flags that are prone to changes
CUDA         = nvcc
CFLAGSCUDA   = $(CFLAGSCOPROC_FROM_FEAT2) $(APPONLYFLAGS)


##############################################################################
# Commands to get version information from compiler
##############################################################################
CUDAVERSION  = $(CUDA)  --version | tail -n 1


##############################################################################
# compiler flags 
# (including non-architecture specific optimisation flags)
##############################################################################


# Add CUDA support
CFLAGSCUDA := $(CFLAGSCUDA) -DHAS_CUDA

# Set CUDA_INC_PATH if it is not empty
ifneq ($(strip $(CUDA_INC_PATH)),)
CFLAGSCUDA := $(CFLAGSCUDA) -I$(CUDA_INC_PATH)
endif


# Set default compile flags part 1
ifeq ($(call optimise), YES)
CFLAGSCUDA := $(CFLAGSCUDA) -O3 --ptxas-options=-v
else
CFLAGSCUDA := $(CFLAGSCUDA) -g -DENABLE_PARAMETER_CHECK --ptxas-options=-v
endif


# Set compile flags part 2
ifeq ($(strip $(HAS_CUDA10)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch sm_10 
endif
ifeq ($(strip $(HAS_CUDA11)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch sm_11
endif
ifeq ($(strip $(HAS_CUDA12)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch sm_12
endif
ifeq ($(strip $(HAS_CUDA13)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch=compute_13 -code=compute_13 -m64  -DHAS_CUDADOUBLEPREC
endif
ifeq ($(strip $(HAS_CUDA20)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch=compute_20 -code=compute_20 -m64  -DHAS_CUDADOUBLEPREC
endif
ifeq ($(strip $(HAS_CUDA21)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch=compute_21 -code=compute_21 -m64  -DHAS_CUDADOUBLEPREC
endif


# Detect compiler version
NVCCVERSION := $(shell eval $(CUDAVERSION) )
