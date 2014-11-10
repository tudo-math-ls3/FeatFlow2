# -*- mode: makefile -*-

##############################################################################
# GNU Compiler suite 4.x with individual settings for ARM CPU IDs
#
##############################################################################

# Option -march=cpu-type specifies that the code is compiled with
# optimization for the specified machine type "cpu-type". For GCC
# versions later than 4.2 -march=native tries to find the optimal
# configuration for the particular CPU the code is compiled on.  A
# detailed list of options enabled/disabled by this switch can be
# obtained by invoking the following command:
#
# gcc -march=native -E -v - </dev/null 2>&1 | grep cc1

##############################################################################
# ARM Cortex CPUs
#
##############################################################################

# List of supported cpu-types and minimally required compiler version
# CPU-TYPE     min. GCC version
# cortex-a5    4.5.0 
# cortex-a7    4.7.0
# cortex-a8    4.3.0
# cortex-a9    4.4.0
# cortex-a15   4.6.0
# cortex-a17   4.9.0
# cortex-a53   4.9.0
# cortex-a57   4.9.0
# cortex-r4    4.3.0
# cortex-r4f   4.4.0
# cortex-r5    4.7.0
# cortex-r7    4.9.0
# cortex-m0    4.5.0
# cortex-m1    4.4.0
# cortex-m3    4.3.0
# cortex-m4    4.4.0
# cortex-m7    ?.?.?

# ARM Cortex A5 CPU.
ifeq ($(call match,$(ID),arm-cortexa5-.*-gcc-.*),yes)
ifeq ($(call gfortranminversion,4,5),yes)
CFLAGSOPTF77 := -mcpu=cortex-a5 -mfpu=neon-fp16 -mfloat-abi=hard
CFLAGSOPTF90 := -mcpu=cortex-a5 -mfpu=neon-fp16 -mfloat-abi=hard
CFLAGSOPTC   := -mcpu=cortex-a5 -mfpu=neon-fp16 -mfloat-abi=hard
CFLAGSOPTCXX := -mcpu=cortex-a5 -mfpu=neon-fp16 -mfloat-abi=hard
LDFLAGSOPT   := -mcpu=cortex-a5 -mfpu=neon-fp16 -mfloat-abi=hard
else
CFLAGSOPTF77 := -mfloat-abi=soft
CFLAGSOPTF90 := -mfloat-abi=soft
CFLAGSOPTC   := -mfloat-abi=soft
CFLAGSOPTCXX := -mfloat-abi=soft
LDFLAGSOPT   := -mfloat-abi=soft
endif
endif

# ARM Cortex A7 CPU. 
ifeq ($(call match,$(ID),arm-cortexa7-.*-gcc-.*),yes)
ifeq ($(call gfortranminversion,4,7),yes)
CFLAGSOPTF77 := -mcpu=cortex-a6 -mfpu=neon-vfpv4 -mfloat-abi=hard
CFLAGSOPTF90 := -mcpu=cortex-a6 -mfpu=neon-vfpv4 -mfloat-abi=hard
CFLAGSOPTC   := -mcpu=cortex-a6 -mfpu=neon-vfpv4 -mfloat-abi=hard
CFLAGSOPTCXX := -mcpu=cortex-a6 -mfpu=neon-vfpv4 -mfloat-abi=hard
LDFLAGSOPT   := -mcpu=cortex-a6 -mfpu=neon-vfpv4 -mfloat-abi=hard
else
CFLAGSOPTF77 := -mfloat-abi=soft
CFLAGSOPTF90 := -mfloat-abi=soft
CFLAGSOPTC   := -mfloat-abi=soft
CFLAGSOPTCXX := -mfloat-abi=soft
LDFLAGSOPT   := -mfloat-abi=soft
endif
endif

# ARM Cortex A8 CPU.
ifeq ($(call match,$(ID),arm-cortexa8-.*-gcc-.*),yes)
ifeq ($(call gfortranminversion,4,3),yes)
CFLAGSOPTF77 := -mcpu=cortex-a8 -mfpu=neon -mfloat-abi=hard
CFLAGSOPTF90 := -mcpu=cortex-a8 -mfpu=neon -mfloat-abi=hard
CFLAGSOPTC   := -mcpu=cortex-a8 -mfpu=neon -mfloat-abi=hard
CFLAGSOPTCXX := -mcpu=cortex-a8 -mfpu=neon -mfloat-abi=hard
LDFLAGSOPT   := -mcpu=cortex-a8 -mfpu=neon -mfloat-abi=hard
else
CFLAGSOPTF77 := -mfloat-abi=soft
CFLAGSOPTF90 := -mfloat-abi=soft
CFLAGSOPTC   := -mfloat-abi=soft
CFLAGSOPTCXX := -mfloat-abi=soft
LDFLAGSOPT   := -mfloat-abi=soft
endif
endif

# ARM Cortex A9 CPU.
ifeq ($(call match,$(ID),arm-cortexa9-.*-gcc-.*),yes)
ifeq ($(call gfortranminversion,4,4),yes)
CFLAGSOPTF77 := -mcpu=cortex-a9 -mfpu=neon-fp16 -mfloat-abi=hard
CFLAGSOPTF90 := -mcpu=cortex-a9 -mfpu=neon-fp16 -mfloat-abi=hard
CFLAGSOPTC   := -mcpu=cortex-a9 -mfpu=neon-fp16 -mfloat-abi=hard
CFLAGSOPTCXX := -mcpu=cortex-a9 -mfpu=neon-fp16 -mfloat-abi=hard
LDFLAGSOPT   := -mcpu=cortex-a9 -mfpu=neon-fp16 -mfloat-abi=hard
else
CFLAGSOPTF77 := -mfloat-abi=soft
CFLAGSOPTF90 := -mfloat-abi=soft
CFLAGSOPTC   := -mfloat-abi=soft
CFLAGSOPTCXX := -mfloat-abi=soft
LDFLAGSOPT   := -mfloat-abi=soft
endif
endif

# ARM Cortex A15 CPU.
ifeq ($(call match,$(ID),arm-cortexa15-.*-gcc-.*),yes)
ifeq ($(call gfortranminversion,4,6),yes)
CFLAGSOPTF77 := -mcpu=cortex-a15 -mfpu=neon-vfpv4 -mfloat-abi=hard
CFLAGSOPTF90 := -mcpu=cortex-a15 -mfpu=neon-vfpv4 -mfloat-abi=hard
CFLAGSOPTC   := -mcpu=cortex-a15 -mfpu=neon-vfpv4 -mfloat-abi=hard
CFLAGSOPTCXX := -mcpu=cortex-a15 -mfpu=neon-vfpv4 -mfloat-abi=hard
LDFLAGSOPT   := -mcpu=cortex-a15 -mfpu=neon-vfpv4 -mfloat-abi=hard
else
CFLAGSOPTF77 := -mfloat-abi=soft
CFLAGSOPTF90 := -mfloat-abi=soft
CFLAGSOPTC   := -mfloat-abi=soft
CFLAGSOPTCXX := -mfloat-abi=soft
LDFLAGSOPT   := -mfloat-abi=soft
endif
endif

# ARM Cortex A17 CPU.
ifeq ($(call match,$(ID),arm-cortexa17-.*-gcc-.*),yes)
ifeq ($(call gfortranminversion,4,9),yes)
CFLAGSOPTF77 := -mcpu=cortex-a17 -mfpu=neon-vfpv4 -mfloat-abi=hard
CFLAGSOPTF90 := -mcpu=cortex-a17 -mfpu=neon-vfpv4 -mfloat-abi=hard
CFLAGSOPTC   := -mcpu=cortex-a17 -mfpu=neon-vfpv4 -mfloat-abi=hard
CFLAGSOPTCXX := -mcpu=cortex-a17 -mfpu=neon-vfpv4 -mfloat-abi=hard
LDFLAGSOPT   := -mcpu=cortex-a17 -mfpu=neon-vfpv4 -mfloat-abi=hard
else
CFLAGSOPTF77 := -mfloat-abi=soft
CFLAGSOPTF90 := -mfloat-abi=soft
CFLAGSOPTC   := -mfloat-abi=soft
CFLAGSOPTCXX := -mfloat-abi=soft
LDFLAGSOPT   := -mfloat-abi=soft
endif
endif

# ARM Cortex A53 CPU.
ifeq ($(call match,$(ID),arm-cortexa53-.*-gcc-.*),yes)
ifeq ($(call gfortranminversion,4,9),yes)
CFLAGSOPTF77 := -mcpu=cortex-a53 -mfpu=neon-vfpv4 -mfloat-abi=hard
CFLAGSOPTF90 := -mcpu=cortex-a53 -mfpu=neon-vfpv4 -mfloat-abi=hard
CFLAGSOPTC   := -mcpu=cortex-a53 -mfpu=neon-vfpv4 -mfloat-abi=hard
CFLAGSOPTCXX := -mcpu=cortex-a53 -mfpu=neon-vfpv4 -mfloat-abi=hard
LDFLAGSOPT   := -mcpu=cortex-a53 -mfpu=neon-vfpv4 -mfloat-abi=hard
else
CFLAGSOPTF77 := -mfloat-abi=soft
CFLAGSOPTF90 := -mfloat-abi=soft
CFLAGSOPTC   := -mfloat-abi=soft
CFLAGSOPTCXX := -mfloat-abi=soft
LDFLAGSOPT   := -mfloat-abi=soft
endif
endif

# ARM Cortex A57 CPU.
ifeq ($(call match,$(ID),arm-cortexa57-.*-gcc-.*),yes)
ifeq ($(call gfortranminversion,4,9),yes)
CFLAGSOPTF77 := -mcpu=cortex-a57 -mfpu=neon-vfpv4 -mfloat-abi=hard
CFLAGSOPTF90 := -mcpu=cortex-a57 -mfpu=neon-vfpv4 -mfloat-abi=hard
CFLAGSOPTC   := -mcpu=cortex-a57 -mfpu=neon-vfpv4 -mfloat-abi=hard
CFLAGSOPTCXX := -mcpu=cortex-a57 -mfpu=neon-vfpv4 -mfloat-abi=hard
LDFLAGSOPT   := -mcpu=cortex-a57 -mfpu=neon-vfpv4 -mfloat-abi=hard
else
CFLAGSOPTF77 := -mfloat-abi=soft
CFLAGSOPTF90 := -mfloat-abi=soft
CFLAGSOPTC   := -mfloat-abi=soft
CFLAGSOPTCXX := -mfloat-abi=soft
LDFLAGSOPT   := -mfloat-abi=soft
endif
endif

##############################################################################
# Automatic CPU type detection
#
##############################################################################

ifneq (,$(findstring NATIVE,$(OPT)))
ifeq ($(call gfortranminversion,4,2),yes)
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by GNU gfortran compiler!';
CFLAGSOPTF77 := -march=native
CFLAGSOPTF90 := -march=native
LDFLAGSOPT   := -march=native
else
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: No auto-tuning for GNU gfortran below 4.2.0 compiler available!';
endif
ifeq ($(call gccminversion,4,2),yes)
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by GNU GCC compiler!';
CFLAGSOPTC   := -march=native
CFLAGSOPTCXX := -march=native
else
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: No auto-tuning for GNU GCC below 4.2.0 compiler available!';
endif
endif
