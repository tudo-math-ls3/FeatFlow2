# -*- mode: makefile -*-

##############################################################################
# SunStudio Compiler suite with individual settings for x86/x86_64 CPU IDs
#
##############################################################################

# Option -xtarget=cpu-type specifies that the code is compiled with
# optimization for the specified machine type "cpu-type".

# Sun UltraSPARC-I CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCI-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=ultra
CFLAGSOPTF90 := -xtarget=ultra
CFLAGSOPTC   := -xtarget=ultra
CFLAGSOPTCXX := -xtarget=ultra
LDFLAGSOPT   := -xtarget=ultra
endif

# Sun UltraSPARC-II CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCII-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=ultra2
CFLAGSOPTF90 := -xtarget=ultra2
CFLAGSOPTC   := -xtarget=ultra2
CFLAGSOPTCXX := -xtarget=ultra2
LDFLAGSOPT   := -xtarget=ultra2
endif

# Sun UltraSPARC-IIi CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCIIi-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=ultra2i
CFLAGSOPTF90 := -xtarget=ultra2i
CFLAGSOPTC   := -xtarget=ultra2i
CFLAGSOPTCXX := -xtarget=ultra2i
LDFLAGSOPT   := -xtarget=ultra2i
endif

# Sun UltraSPARC-IIe CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCIIe-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=ultra2e
CFLAGSOPTF90 := -xtarget=ultra2e
CFLAGSOPTC   := -xtarget=ultra2e
CFLAGSOPTCXX := -xtarget=ultra2e
LDFLAGSOPT   := -xtarget=ultra2e
endif

# Sun UltraSPARC-III CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCIII-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=ultra3
CFLAGSOPTF90 := -xtarget=ultra3
CFLAGSOPTC   := -xtarget=ultra3
CFLAGSOPTCXX := -xtarget=ultra3
LDFLAGSOPT   := -xtarget=ultra3
endif

# Sun UltraSPARC-IIIi CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCIIIi-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=ultra3i
CFLAGSOPTF90 := -xtarget=ultra3i
CFLAGSOPTC   := -xtarget=ultra3i
CFLAGSOPTCXX := -xtarget=ultra3i
LDFLAGSOPT   := -xtarget=ultra3i
endif


# Sun UltraSPARC-IV CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCIV-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=ultra4
CFLAGSOPTF90 := -xtarget=ultra4
CFLAGSOPTC   := -xtarget=ultra4
CFLAGSOPTCXX := -xtarget=ultra4
LDFLAGSOPT   := -xtarget=ultra4
endif

# Sun UltraSPARC-IV+ CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCIV+-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=ultra4plus
CFLAGSOPTF90 := -xtarget=ultra4plus
CFLAGSOPTC   := -xtarget=ultra4plus
CFLAGSOPTCXX := -xtarget=ultra4plus
LDFLAGSOPT   := -xtarget=ultra4plus
endif

# Sun SPARC64-VI CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-SPARC64VI-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=sparc64vi
CFLAGSOPTF90 := -xtarget=sparc64vi
CFLAGSOPTC   := -xtarget=sparc64vi
CFLAGSOPTCXX := -xtarget=sparc64vi
LDFLAGSOPT   := -xtarget=sparc64vi
endif

# Sun SPARC64-VII CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-SPARC64VII-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=sparc64vii
CFLAGSOPTF90 := -xtarget=sparc64vii
CFLAGSOPTC   := -xtarget=sparc64vii
CFLAGSOPTCXX := -xtarget=sparc64vii
LDFLAGSOPT   := -xtarget=sparc64vii
endif

# Sun SPARC64-VII+ CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-SPARC64VII+-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=sparc64viiplus
CFLAGSOPTF90 := -xtarget=sparc64viiplus
CFLAGSOPTC   := -xtarget=sparc64viiplus
CFLAGSOPTCXX := -xtarget=sparc64viiplus
LDFLAGSOPT   := -xtarget=sparc64viiplus
endif

# Sun UltraSPARC-T1 CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCT1-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=ultraT1
CFLAGSOPTF90 := -xtarget=ultraT1
CFLAGSOPTC   := -xtarget=ultraT1
CFLAGSOPTCXX := -xtarget=ultraT1
LDFLAGSOPT   := -xtarget=ultraT1
endif

# Sun UltraSPARC-T2 CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCT2-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=ultraT2
CFLAGSOPTF90 := -xtarget=ultraT2
CFLAGSOPTC   := -xtarget=ultraT2
CFLAGSOPTCXX := -xtarget=ultraT2
LDFLAGSOPT   := -xtarget=ultraT2
endif

# Sun UltraSPARC-T3 CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCT3-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=T3
CFLAGSOPTF90 := -xtarget=T3
CFLAGSOPTC   := -xtarget=T3
CFLAGSOPTCXX := -xtarget=T3
LDFLAGSOPT   := -xtarget=T3
endif

# Sun UltraSPARC-T4 CPU.
ifeq ($(call match,$(ID),(sun4u|sun4v)-UltraSPARCT4-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=T4
CFLAGSOPTF90 := -xtarget=T4
CFLAGSOPTC   := -xtarget=T4
CFLAGSOPTCXX := -xtarget=T4
LDFLAGSOPT   := -xtarget=T4
endif

##############################################################################
# Automatic CPU type detection
#
##############################################################################

ifneq (,$(findstring NATIVE,$(OPT)))
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by SunStudio compiler!';
CFLAGSOPTF77 := -xtarget=native
CFLAGSOPTF90 := -xtarget=native
CFLAGSOPTC   := -xtarget=native
CFLAGSOPTCXX := -xtarget=native
LDFLAGSOPT   := -xtarget=native
endif
