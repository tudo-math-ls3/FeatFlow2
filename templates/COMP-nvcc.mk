# -*- mode: makefile -*-

##############################################################################
# NVIDIA Cuda compiler 3.x, 4.x, 5.x
#
##############################################################################
COMPILERNAME = NVCC

# Default: Inherit flags that are prone to changes
CUDA         = nvcc
CFLAGSCUDA   = $(CFLAGSCOPROC_FROM_FEAT2) $(APPONLYFLAGS)


##############################################################################
# Commands to get version information from compiler
##############################################################################
CUDAVERSION  = $(CUDA) --version | tail -n 1


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



# Set default type of integer variables explicitly
ifeq ($(strip $(INTSIZE)), LARGE)
CFLAGSCUDA := $(CFLAGSCUDA) -DUSE_LARGEINT
endif



# Set default compile flags part 1
ifeq ($(call optimise), YES)
CFLAGSCUDA := $(CFLAGSCUDA) -O3 --ptxas-options=-v
else
CFLAGSCUDA := $(CFLAGSCUDA) -DDEBUG -g -G -DENABLE_PARAMETER_CHECK \
	                    --ptxas-options=-v
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
CFLAGSCUDA := $(CFLAGSCUDA)  -arch=compute_13 -code=sm_13 -m64  -DHAS_CUDADOUBLEPREC
endif
ifeq ($(strip $(HAS_CUDA20)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch=compute_20 -code=sm_20 -m64  -DHAS_CUDADOUBLEPREC
endif
ifeq ($(strip $(HAS_CUDA21)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch=compute_21 -code=sm_21 -m64  -DHAS_CUDADOUBLEPREC
endif
ifeq ($(strip $(HAS_CUDA30)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch=compute_30 -code=sm_30 -m64  -DHAS_CUDADOUBLEPREC
endif
ifeq ($(strip $(HAS_CUDA32)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch=compute_32 -code=sm_32 -m64  -DHAS_CUDADOUBLEPREC
endif
ifeq ($(strip $(HAS_CUDA35)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch=compute_35 -code=sm_35 -m64  -DHAS_CUDADOUBLEPREC
endif
ifeq ($(strip $(HAS_CUDA37)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch=compute_37 -code=sm_37 -m64  -DHAS_CUDADOUBLEPREC
endif
ifeq ($(strip $(HAS_CUDA50)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch=compute_50 -code=sm_50 -m64  -DHAS_CUDADOUBLEPREC
endif
ifeq ($(strip $(HAS_CUDA52)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch=compute_52 -code=sm_52 -m64  -DHAS_CUDADOUBLEPREC
endif
ifeq ($(strip $(HAS_CUDA53)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch=compute_53 -code=sm_53 -m64  -DHAS_CUDADOUBLEPREC
endif



# Detect compiler version
NVCCVERSION := $(shell eval $(CUDAVERSION))
ifneq (,$(findstring Cuda,$(NVCCVERSION)))
NVCCVERSION_MAJOR := $(shell eval $(CUDAVERSION) | cut -d' ' -f5 | cut -d'.' -f1)
NVCCVERSION_MINOR := $(shell eval $(CUDAVERSION) | cut -d' ' -f5 | cut -d'.' -f2 | cut -d',' -f1)
else
NVCCVERSION_MAJOR := 0
NVCCVERSION_MINOR := 0
endif

# Functions to detect minimal compiler version
nvccminversion = $(shell if [ $(NVCCVERSION_MAJOR) -gt $(1) ] || \
	                   ([ $(NVCCVERSION_MAJOR) -ge $(1) ] && [ $(NVCCVERSION_MINOR) -ge $(2) ]) ; then echo yes ; else echo no ; fi)

# Functions to detect maximal compiler version
nvccmaxversion = $(shell if [ $(NVCCVERSION_MAJOR) -lt $(1) ] || \
	                   ([ $(NVCCVERSION_MAJOR) -le $(1) ] && [ $(NVCCVERSION_MINOR) -le $(2) ]) ; then echo yes ; else echo no ; fi)



# Set features which require special minimum CUDA version
ifeq ($(call nvccminversion,4,0),yes)
CFLAGSCUDA := $(CFLAGSCUDA) -DHAS_INLINE_PTX
endif
