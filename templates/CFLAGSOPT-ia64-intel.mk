# -*- mode: makefile -*-

##############################################################################
# Intel Compiler suite with individual settings for Itanium CPU IDs
#
##############################################################################

# The intel compiler suite frequently changed the naming convention
# for optimization flags so that it is necessary to carefully check
# for the compiler versions.
#
# Intel Compiler suite later than and including 11.0 support auto-tuning
#
# ifort -xHost -dryrun dummy
#
# Intel Compiler suite prior to and including 8.1
#
# -tpp1  Intel Itanium CPU
# -tpp2  Intel Itanium2 CPU
#
#
# Intel Compiler suite prior to and including 9.1
#
# -mcpu=<cpu>  optimize for a specific cpu
#       itanium - optimize for Itanium(R) processor
#       itanium2 - optimize for Itanium(R) 2 processor (DEFAULT)
#       itanium2-p9000 - optimize for Dual-Core Intel(R) Itanium(R) 2 Processor 9000 Sequence
# -mtune=<cpu> optimize for a specific cpu
#       itanium - optimize for Itanium(R) processor
#       itanium2 - optimize for Itanium(R) 2 processor (DEFAULT)
#       itanium2-p9000 - optimize for Dual-Core Intel(R) Itanium(R) 2 Processor 9000 Sequence

##############################################################################
# Intel IA-64 CPUs
#
##############################################################################

# Intel Itanium CPU with 64-bit extensions
ifeq ($(call match,$(ID),ia64-itanium-.*-intel-.*),yes)
ifeq ($(call intelminversion,9,1),yes)
CFLAGSOPTF77 := -mtune=itanium
CFLAGSOPTF90 := -mtune=itanium
CFLAGSOPTC   := -mtune=itanium
CFLAGSOPTCXX := -mtune=itanium
LDFLAGSOPT   := -mtune=itanium
else
CFLAGSOPTF77 := -tpp1
CFLAGSOPTF90 := -tpp1
CFLAGSOPTC   := -tpp1
CFLAGSOPTCXX := -tpp1
LDFLAGSOPT   := -tpp1
endif
endif

# Intel Itanium 2 CPU with 64-bit extensions
ifeq ($(call match,$(ID),ia64-itanium(2|2x2)-.*-intel-.*),yes)
ifeq ($(call intelminversion,9,1),yes)
CFLAGSOPTF77 := -mtune=itanium2
CFLAGSOPTF90 := -mtune=itanium2
CFLAGSOPTC   := -mtune=itanium2
CFLAGSOPTCXX := -mtune=itanium2
LDFLAGSOPT   := -mtune=itanium2
else
CFLAGSOPTF77 := -tpp2
CFLAGSOPTF90 := -tpp2
CFLAGSOPTC   := -tpp2
CFLAGSOPTCXX := -tpp2
LDFLAGSOPT   := -tpp2
endif
endif

##############################################################################
# Automatic CPU type detection
#
##############################################################################

ifneq (,$(findstring NATIVE,$(OPT)))
ifeq ($(call intelminversion,11,0),yes)
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by Intel compiler!';
CFLAGSOPTF77 := -xHost
CFLAGSOPTF90 := -xHost
CFLAGSOPTC   := -xHost
CFLAGSOPTCXX := -xHost
LDFLAGSOPT   := -xHost
else
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: No auto-tuning for Intel compiler below 11.0 available!';
endif
endif