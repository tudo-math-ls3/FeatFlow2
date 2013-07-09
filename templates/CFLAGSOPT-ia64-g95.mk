# -*- mode: makefile -*-

##############################################################################
# G95 (www.g95.org) with individual settings for Itanium CPU IDs
#
##############################################################################

# Option -march=cpu-type specifies that the code is compiled with
# optimization for the specified machine type "cpu-type". Since G95 is
# based on GCC 4.1.2, the handy option -march=native is not
# available. Moreover, optimization for newer CPUs is not possible.

##############################################################################
# Intel Itanium CPUs
#
##############################################################################

# List of supported cpu-types and minimally required compiler version
# CPU-TYPE     min. GCC version   min. G95 version
# Itanium      4.0.x
# Itanium1     4.0.x
# Merced       4.0.x
# Itanium2     4.0.x
# McKinley     4.0.x

# Intel Itanium CPU with 64-bit extensions
ifeq ($(call match,$(ID),ia64-itanium-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=itanium
CFLAGSOPTF90 := -march=itanium
CFLAGSOPTC   := -march=itanium
CFLAGSOPTCXX := -march=itanium
LDFLAGSOPT   := -march=itanium
endif

# Intel Itanium 2 CPU with 64-bit extensions
ifeq ($(call match,$(ID),ia64-itanium(2|2x2)-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=itanium2
CFLAGSOPTF90 := -march=itanium2
CFLAGSOPTC   := -march=itanium2
CFLAGSOPTCXX := -march=itanium2
LDFLAGSOPT   := -march=itanium2
endif

##############################################################################
# Automatic CPU type detection
#
##############################################################################

ifneq (,$(findstring NATIVE,$(OPT)))
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: No auto-tuning for g95 compiler available!';
endif
