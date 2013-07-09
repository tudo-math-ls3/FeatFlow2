# -*- mode: makefile -*-

##############################################################################
# GNU Compiler suite 4.x with individual settings for Itanium CPU IDs
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
# Intel IA-64 CPUs
#
##############################################################################

# List of supported cpu-types and minimally required compiler version
# CPU-TYPE     min. GCC version
# Itanium      4.0.x
# Itanium1     4.0.x
# Merced       4.0.x
# Itanium2     4.0.x
# McKinley     4.0.x

# Intel Itanium CPU with 64-bit extensions
ifeq ($(call match,$(ID),ia64-itanium-.*-gcc-.*),yes)
CFLAGSOPTF77 := -march=itanium
CFLAGSOPTF90 := -march=itanium
CFLAGSOPTC   := -march=itanium
CFLAGSOPTCXX := -march=itanium
LDFLAGSOPT   := -march=itanium
endif

# Intel Itanium 2 CPU with 64-bit extensions
ifeq ($(call match,$(ID),ia64-itanium(2|2x2)-.*-gcc-.*),yes)
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
ifeq ($(call gfortranminversion,4,2),yes)
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by GNU GCC compiler!';
CFLAGSOPTF77 := -march=native
CFLAGSOPTF90 := -march=native
LDFLAGSOPT   := -march=native
else
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: No auto-tuning for GNU GCC below 4.2.0 compiler available!';
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
