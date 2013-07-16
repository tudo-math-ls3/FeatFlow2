# -*- mode: makefile -*-

##############################################################################
# GNU Compiler suite 4.x with individual settings for x86/x86_64 CPU IDs
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

# List of supported cpu-types and minimally required compiler version
# CPU-TYPE     min. GCC version
# ultrasparc   4.0.x
# ultrasparc3  4.0.x
# niagara      4.2.x
# niagara2     4.3.x
# niagara3     4.7.x
# niagara4     4.7.x
# native       4.7.x

# Sun UltraSPARC-I CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCI-.*-gcc-.*),yes)
CFLAGSOPTF77 := -mcpu=ultrasparc -mtune=ultrasparc
CFLAGSOPTF90 := -mcpu=ultrasparc -mtune=ultrasparc
CFLAGSOPTC   := -mcpu=ultrasparc -mtune=ultrasparc
CFLAGSOPTCXX := -mcpu=ultrasparc -mtune=ultrasparc
LDFLAGSOPT   := -mcpu=ultrasparc -mtune=ultrasparc
endif

# Sun UltraSPARC-II CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCII-.*-gcc-.*),yes)
CFLAGSOPTF77 := -mcpu=ultrasparc -mtune=ultrasparc
CFLAGSOPTF90 := -mcpu=ultrasparc -mtune=ultrasparc
CFLAGSOPTC   := -mcpu=ultrasparc -mtune=ultrasparc
CFLAGSOPTCXX := -mcpu=ultrasparc -mtune=ultrasparc
LDFLAGSOPT   := -mcpu=ultrasparc -mtune=ultrasparc
endif

# Sun UltraSPARC-IIi CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCIIi-.*-gcc-.*),yes)
CFLAGSOPTF77 := -mcpu=ultrasparc -mtune=ultrasparc
CFLAGSOPTF90 := -mcpu=ultrasparc -mtune=ultrasparc
CFLAGSOPTC   := -mcpu=ultrasparc -mtune=ultrasparc
CFLAGSOPTCXX := -mcpu=ultrasparc -mtune=ultrasparc
LDFLAGSOPT   := -mcpu=ultrasparc -mtune=ultrasparc
endif

# Sun UltraSPARC-IIe CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCIIe-.*-gcc-.*),yes)
CFLAGSOPTF77 := -mcpu=ultrasparc -mtune=ultrasparc
CFLAGSOPTF90 := -mcpu=ultrasparc -mtune=ultrasparc
CFLAGSOPTC   := -mcpu=ultrasparc -mtune=ultrasparc
CFLAGSOPTCXX := -mcpu=ultrasparc -mtune=ultrasparc
LDFLAGSOPT   := -mcpu=ultrasparc -mtune=ultrasparc
endif

# Sun UltraSPARC-III CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCIII-.*-gcc-.*),yes)
CFLAGSOPTF77 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTF90 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTC   := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTCXX := -mcpu=ultrasparc3 -mtune=ultrasparc3
LDFLAGSOPT   := -mcpu=ultrasparc3 -mtune=ultrasparc3
endif

# Sun UltraSPARC-IIIi CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCIIIi-.*-gcc-.*),yes)
CFLAGSOPTF77 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTF90 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTC   := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTCXX := -mcpu=ultrasparc3 -mtune=ultrasparc3
LDFLAGSOPT   := -mcpu=ultrasparc3 -mtune=ultrasparc3
endif


# Sun UltraSPARC-IV CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCIV-.*-gcc-.*),yes)
CFLAGSOPTF77 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTF90 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTC   := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTCXX := -mcpu=ultrasparc3 -mtune=ultrasparc3
LDFLAGSOPT   := -mcpu=ultrasparc3 -mtune=ultrasparc3
endif

# Sun UltraSPARC-IV+ CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCIV+-.*-gcc-.*),yes)
CFLAGSOPTF77 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTF90 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTC   := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTCXX := -mcpu=ultrasparc3 -mtune=ultrasparc3
LDFLAGSOPT   := -mcpu=ultrasparc3 -mtune=ultrasparc3
endif

# Sun SPARC64-VI CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-SPARC64VI-.*-gcc-.*),yes)
CFLAGSOPTF77 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTF90 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTC   := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTCXX := -mcpu=ultrasparc3 -mtune=ultrasparc3
LDFLAGSOPT   := -mcpu=ultrasparc3 -mtune=ultrasparc3
endif

# Sun SPARC64-VII CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-SPARC64VII-.*-gcc-.*),yes)
CFLAGSOPTF77 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTF90 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTC   := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTCXX := -mcpu=ultrasparc3 -mtune=ultrasparc3
LDFLAGSOPT   := -mcpu=ultrasparc3 -mtune=ultrasparc3
endif

# Sun SPARC64-VII+ CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-SPARC64VII+-.*-gcc-.*),yes)
CFLAGSOPTF77 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTF90 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTC   := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTCXX := -mcpu=ultrasparc3 -mtune=ultrasparc3
LDFLAGSOPT   := -mcpu=ultrasparc3 -mtune=ultrasparc3
endif

# Sun UltraSPARC-T1 CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCT1-.*-gcc-.*),yes)
ifeq ($(call gfortranminversion,4,2),yes)
CFLAGSOPTF77 := -mcpu=niagara -mtune=niagara
CFLAGSOPTF90 := -mcpu=niagara -mtune=niagara
LDFLAGSOPT   := -mcpu=niagara -mtune=niagara
else
CFLAGSOPTF77 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTF90 := -mcpu=ultrasparc3 -mtune=ultrasparc3
LDFLAGSOPT   := -mcpu=ultrasparc3 -mtune=ultrasparc3
endif
ifeq ($(call gccminversion,4,2),yes)
CFLAGSOPTC   := -mcpu=niagara -mtune=niagara
CFLAGSOPTCXX := -mcpu=niagara -mtune=niagara
else
CFLAGSOPTC   := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTCXX := -mcpu=ultrasparc3 -mtune=ultrasparc3
endif
endif

# Sun UltraSPARC-T2 CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCT2-.*-gcc-.*),yes)
ifeq ($(call gfortranminversion,4,3),yes)
CFLAGSOPTF77 := -mcpu=niagara2 -mtune=niagara2
CFLAGSOPTF90 := -mcpu=niagara2 -mtune=niagara2
LDFLAGSOPT   := -mcpu=niagara2 -mtune=niagara2
else
CFLAGSOPTF77 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTF90 := -mcpu=ultrasparc3 -mtune=ultrasparc3
LDFLAGSOPT   := -mcpu=ultrasparc3 -mtune=ultrasparc3
endif
ifeq ($(call gccminversion,4,3),yes)
CFLAGSOPTC   := -mcpu=niagara2 -mtune=niagara2
CFLAGSOPTCXX := -mcpu=niagara2 -mtune=niagara2
else
CFLAGSOPTC   := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTCXX := -mcpu=ultrasparc3 -mtune=ultrasparc3
endif
endif

# Sun UltraSPARC-T3 CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCT3-.*-gcc-.*),yes)
ifeq ($(call gfortranminversion,4,7),yes)
CFLAGSOPTF77 := -mcpu=niagara3 -mtune=niagara3
CFLAGSOPTF90 := -mcpu=niagara3 -mtune=niagara3
LDFLAGSOPT   := -mcpu=niagara3 -mtune=niagara3
else
ifeq ($(call gfortranminversion,4,3),yes)
CFLAGSOPTF77 := -mcpu=niagara2 -mtune=niagara2
CFLAGSOPTF90 := -mcpu=niagara2 -mtune=niagara2
LDFLAGSOPT   := -mcpu=niagara2 -mtune=niagara2
else
CFLAGSOPTF77 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTF90 := -mcpu=ultrasparc3 -mtune=ultrasparc3
LDFLAGSOPT   := -mcpu=ultrasparc3 -mtune=ultrasparc3
endif
endif
ifeq ($(call gccminversion,4,7),yes)
CFLAGSOPTC   := -mcpu=niagara3 -mtune=niagara3
CFLAGSOPTCXX := -mcpu=niagara3 -mtune=niagara3
else
ifeq ($(call gccminversion,4,3),yes)
CFLAGSOPTC   := -mcpu=niagara2 -mtune=niagara2
CFLAGSOPTCXX := -mcpu=niagara2 -mtune=niagara2
else
CFLAGSOPTC   := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTCXX := -mcpu=ultrasparc3 -mtune=ultrasparc3
endif
endif
endif

# Sun UltraSPARC-T4 CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCT4-.*-gcc-.*),yes)
ifeq ($(call gfortranminversion,4,7),yes)
CFLAGSOPTF77 := -mcpu=niagara4 -mtune=niagara4
CFLAGSOPTF90 := -mcpu=niagara4 -mtune=niagara4
LDFLAGSOPT   := -mcpu=niagara4 -mtune=niagara4
ifeq ($(call gfortranminversion,4,3),yes)
CFLAGSOPTF77 := -mcpu=niagara2 -mtune=niagara2
CFLAGSOPTF90 := -mcpu=niagara2 -mtune=niagara2
LDFLAGSOPT   := -mcpu=niagara2 -mtune=niagara2
else
CFLAGSOPTF77 := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTF90 := -mcpu=ultrasparc3 -mtune=ultrasparc3
LDFLAGSOPT   := -mcpu=ultrasparc3 -mtune=ultrasparc3
endif
endif
ifeq ($(call gccminversion,4,7),yes)
CFLAGSOPTC   := -mcpu=niagara4 -mtune=niagara4
CFLAGSOPTCXX := -mcpu=niagara4 -mtune=niagara4
ifeq ($(call gccminversion,4,3),yes)
CFLAGSOPTC   := -mcpu=niagara2 -mtune=niagara2
CFLAGSOPTCXX := -mcpu=niagara2 -mtune=niagara2
else
CFLAGSOPTC   := -mcpu=ultrasparc3 -mtune=ultrasparc3
CFLAGSOPTCXX := -mcpu=ultrasparc3 -mtune=ultrasparc3
endif
endif
endif

##############################################################################
# Automatic CPU type detection
#
##############################################################################

ifneq (,$(findstring NATIVE,$(OPT)))
ifeq ($(call gfortranminversion,4,7),yes)
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by GCC compiler!';
CFLAGSOPTF77 := -mcpu=native -mtune=native
CFLAGSOPTF90 := -mcpu=native -mtune=native
LDFLAGSOPT   := -mcpu=native -mtune=native
else
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: No auto-tuning for GNU gfortran below 4.7.0 compiler available!';
endif
endif
ifeq ($(call gccminversion,4,7),yes)
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by GCC compiler!';
CFLAGSOPTC   := -mcpu=native -mtune=native
CFLAGSOPTCXX := -mcpu=native -mtune=native
else
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: No auto-tuning for GNU gcc below 4.7.0 compiler available!';
endif
endif
endif
