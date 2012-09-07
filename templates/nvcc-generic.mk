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
ifeq ($(strip $(HAS_CUDA35)),YES)
CFLAGSCUDA := $(CFLAGSCUDA)  -arch=compute_35 -code=sm_35 -m64  -DHAS_CUDADOUBLEPREC
endif



# Detect compiler version
NVCCVERSION := $(shell eval $(CUDAVERSION) )

# Functions to detect minimal compiler version
nvccminversion_5_0=\
	$(if $(findstring 5.0,$(NVCCVERSION)),yes,no)
nvccminversion_4_2=\
	$(if $(findstring yes,\
	$(call nvccminversion_5_0) \
	$(if $(findstring 4.2,$(NVCCVERSION)),yes,no)),yes,no)
nvccminversion_4_1=\
	$(if $(findstring yes,\
	$(call nvccminversion_4_2) \
	$(if $(findstring 4.1,$(NVCCVERSION)),yes,no)),yes,no)
nvccminversion_4_0=\
	$(if $(findstring yes,\
	$(call nvccminversion_4_1) \
	$(if $(findstring 4.0,$(NVCCVERSION)),yes,no)),yes,no)
nvccminversion_3_2=\
	$(if $(findstring yes,\
	$(call nvccminversion_4_0) \
	$(if $(findstring 3.2,$(NVCCVERSION)),yes,no)),yes,no)
nvccminversion_3_1=\
	$(if $(findstring yes,\
	$(call nvccminversion_3_2) \
	$(if $(findstring 3.1,$(NVCCVERSION)),yes,no)),yes,no)
nvccminversion_3_0=\
	$(if $(findstring yes,\
	$(call nvccminversion_3_1) \
	$(if $(findstring 3.0,$(NVCCVERSION)),yes,no)),yes,no)
nvccminversion_2_3=\
	$(if $(findstring yes,\
	$(call nvccminversion_3_0) \
	$(if $(findstring 2.3,$(NVCCVERSION)),yes,no)),yes,no)
nvccminversion_2_2=\
	$(if $(findstring yes,\
	$(call nvccminversion_2_3) \
	$(if $(findstring 2.2,$(NVCCVERSION)),yes,no)),yes,no)
nvccminversion_2_1=\
	$(if $(findstring yes,\
	$(call nvccminversion_2_2) \
	$(if $(findstring 2.1,$(NVCCVERSION)),yes,no)),yes,no)
nvccminversion_2_0=\
	$(if $(findstring yes,\
	$(call nvccminversion_2_1) \
	$(if $(findstring 2.0,$(NVCCVERSION)),yes,no)),yes,no)
nvccminversion_1_1=\
	$(if $(findstring yes,\
	$(call nvccminversion_2_0) \
	$(if $(findstring 1.1,$(NVCCVERSION)),yes,no)),yes,no)
nvccminversion_1_0=\
	$(if $(findstring yes,\
	$(call nvccminversion_1_1) \
	$(if $(findstring 1.0,$(NVCCVERSION)),yes,no)),yes,no)

# Functions to detect maximal compiler version
nvccmaxversion_5_0=\
	$(if $(findstring yes,\
	$(call nvccmaxversion_4_2) \
	$(if $(findstring 5.0,$(NVCCVERSION)),yes,no)),yes,no)
nvccmaxversion_4_2=\
	$(if $(findstring yes,\
	$(call nvccmaxversion_4_1) \
	$(if $(findstring 4.2,$(NVCCVERSION)),yes,no)),yes,no)
nvccmaxversion_4_1=\
	$(if $(findstring yes,\
	$(call nvccmaxversion_4_0) \
	$(if $(findstring 4.1,$(NVCCVERSION)),yes,no)),yes,no)
nvccmaxversion_4_0=\
	$(if $(findstring yes,\
	$(call nvccmaxversion_3_2) \
	$(if $(findstring 4.0,$(NVCCVERSION)),yes,no)),yes,no)
nvccmaxversion_3_2=\
	$(if $(findstring yes,\
	$(call nvccmaxversion_3_1) \
	$(if $(findstring 3.2,$(NVCCVERSION)),yes,no)),yes,no)
nvccmaxversion_3_1=\
	$(if $(findstring yes,\
	$(call nvccmaxversion_3_0) \
	$(if $(findstring 3.1,$(NVCCVERSION)),yes,no)),yes,no)
nvccmaxversion_3_0=\
	$(if $(findstring yes,\
	$(call nvccmaxversion_2_3) \
	$(if $(findstring 3.0,$(NVCCVERSION)),yes,no)),yes,no)
nvccmaxversion_2_3=\
	$(if $(findstring yes,\
	$(call nvccmaxversion_2_2) \
	$(if $(findstring 2.3,$(NVCCVERSION)),yes,no)),yes,no)
nvccmaxversion_2_2=\
	$(if $(findstring yes,\
	$(call nvccmaxversion_2_1) \
	$(if $(findstring 2.2,$(NVCCVERSION)),yes,no)),yes,no)
nvccmaxversion_2_1=\
	$(if $(findstring yes,\
	$(call nvccmaxversion_2_0) \
	$(if $(findstring 2.1,$(NVCCVERSION)),yes,no)),yes,no)
nvccmaxversion_2_0=\
	$(if $(findstring yes,\
	$(call nvccmaxversion_1_1) \
	$(if $(findstring 2.0,$(NVCCVERSION)),yes,no)),yes,no)
nvccmaxversion_1_1=\
	$(if $(findstring yes,\
	$(call nvccmaxversion_1_0) \
	$(if $(findstring 1.1,$(NVCCVERSION)),yes,no)),yes,no)
nvccmaxversion_1_0=\
	$(if $(findstring 1.0,$(NVCCVERSION)),yes,no)



# Set features which require special minimum CUDA version
ifeq ($(call nvccminversion_4_0),yes)
CFLAGSCUDA := $(CFLAGSCUDA) -DHAS_INLINE_PTX
endif