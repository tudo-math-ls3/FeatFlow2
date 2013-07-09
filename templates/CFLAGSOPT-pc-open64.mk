# -*- mode: makefile -*-

##############################################################################
# Open64 Compiler suite with individual settings for x86/x86_64 CPU IDs
#
##############################################################################

# Option -march=cpu-type specifies that the code is compiled with
# optimization for the specified machine type "cpu-type". A detailed
# list of options enabled/disabled by this switch can be obtained by
# invoking the following command:
#
# opencc -show -march=auto file.c

##############################################################################
# Intel x86/x86_64 CPUs
#
##############################################################################

# Original Intel i386 CPU.
ifeq ($(call match,$(ID),pc-i386-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=anyx86
CFLAGSOPTF90 := -march=anyx86
CFLAGSOPTC   := -march=anyx86
CFLAGSOPTCXX := -march=anyx86
LDFLAGSOPT   := -march=anyx86
endif

# Intel i486 CPU. (No scheduling is implemented for this chip.)
ifeq ($(call match,$(ID),pc-i486-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=anyx86
CFLAGSOPTF90 := -march=anyx86
CFLAGSOPTC   := -march=anyx86
CFLAGSOPTCXX := -march=anyx86
LDFLAGSOPT   := -march=anyx86
endif

# Intel Pentium CPU with or without MMX support.
ifeq ($(call match,$(ID),pc-(i568|pentium)-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=anyx86
CFLAGSOPTF90 := -march=anyx86
CFLAGSOPTC   := -march=anyx86
CFLAGSOPTCXX := -march=anyx86
LDFLAGSOPT   := -march=anyx86
endif

# Intel Pentium Pro CPU.
ifeq ($(call match,$(ID),pc-pentiumpro-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=anyx86
CFLAGSOPTF90 := -march=anyx86
CFLAGSOPTC   := -march=anyx86
CFLAGSOPTCXX := -march=anyx86
LDFLAGSOPT   := -march=anyx86
endif

# Intel Pentium II CPU, based on Pentium Pro core with MMX instruction
# set support.
ifeq ($(call match,$(ID),pc-pentium2-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=anyx86
CFLAGSOPTF90 := -march=anyx86
CFLAGSOPTC   := -march=anyx86
CFLAGSOPTCXX := -march=anyx86
LDFLAGSOPT   := -march=anyx86
endif

# Intel Pentium III CPU, based on Pentium Pro core with MMX and SSE
# instruction set support.
ifeq ($(call match,$(ID),pc-pentium3-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=anyx86 -msse
CFLAGSOPTF90 := -march=anyx86 -msse
CFLAGSOPTC   := -march=anyx86 -msse
CFLAGSOPTCXX := -march=anyx86 -msse
LDFLAGSOPT   := -march=anyx86 -msse
endif

# Intel Pentium M; low-power version of Intel Pentium III CPU with MMX, 
# SSE and SSE2 instruction set support. Used by Centrino notebooks.
ifeq ($(call match,$(ID),pc-pentiumm-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=anyx86 -msse -msse2
CFLAGSOPTF90 := -march=anyx86 -msse -msse2
CFLAGSOPTC   := -march=anyx86 -msse -msse2
CFLAGSOPTCXX := -march=anyx86 -msse -msse2
LDFLAGSOPT   := -march=anyx86 -msse -msse2
endif

# Intel Pentium 4 CPU with MMX, SSE and SSE2 instruction set support.
ifeq ($(call match,$(ID),pc-pentium4-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=pentium4 -msse -msse2
CFLAGSOPTF90 := -march=pentium4 -msse -msse2
CFLAGSOPTC   := -march=pentium4 -msse -msse2
CFLAGSOPTCXX := -march=pentium4 -msse -msse2
LDFLAGSOPT   := -march=pentium4 -msse -msse2
endif

# Intel CoreSolo/Duo. Improved version of Intel Pentium 4 CPU with MMX,
# SSE, SSE2 and SSE3 instruction set support.
ifeq ($(call match,$(ID),pc-coresolo-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=anyx86 -msse -msse2 -msse3
CFLAGSOPTF90 := -march=anyx86 -msse -msse2 -msse3
CFLAGSOPTC   := -march=anyx86 -msse -msse2 -msse3
CFLAGSOPTCXX := -march=anyx86 -msse -msse2 -msse3
LDFLAGSOPT   := -march=anyx86 -msse -msse2 -msse3
endif

# Intel Core 2 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-(coreduo|penryn)-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=core -msse -msse2 -msse3 -mssse3
CFLAGSOPTF90 := -march=core -msse -msse2 -msse3 -mssse3
CFLAGSOPTC   := -march=core -msse -msse2 -msse3 -mssse3
CFLAGSOPTCXX := -march=core -msse -msse2 -msse3 -mssse3
LDFLAGSOPT   := -march=core -msse -msse2 -msse3 -mssse3
endif

# Intel Nehalem CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-nehalem-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=core -msse -msse2 -msse3 -mssse3
CFLAGSOPTF90 := -march=core -msse -msse2 -msse3 -mssse3
CFLAGSOPTC   := -march=core -msse -msse2 -msse3 -mssse3
CFLAGSOPTCXX := -march=core -msse -msse2 -msse3 -mssse3
LDFLAGSOPT   := -march=core -msse -msse2 -msse3 -mssse3
endif

# Intel Core i7 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3,
# SSSE3, SSE4.1, SSE4.2, AVX, AES and PCLMUL instruction set
# support.
ifeq ($(call match,$(ID),(pc|pc64)-sandybridge-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=core -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpclmul
CFLAGSOPTF90 := -march=core -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpclmul
CFLAGSOPTC   := -march=core -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpclmul
CFLAGSOPTCXX := -march=core -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpclmul
LDFLAGSOPT   := -march=core -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpclmul
endif

# Intel Core CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3, SSSE3,
# SSE4.1, SSE4.2, AVX, AES, PCLMUL, FSGSBASE, RDRND and F16C
# instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-ivybridge-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=core -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpclmul
CFLAGSOPTF90 := -march=core -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpclmul
CFLAGSOPTC   := -march=core -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpclmul
CFLAGSOPTCXX := -march=core -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpclmul
LDFLAGSOPT   := -march=core -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpclmul
endif

# Intel Core CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3, SSSE3,
# SSE4.1, SSE4.2, AVX, AES, PCLMUL, FSGSBASE, RDRND and F16C
# instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-haswell-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=core -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpclmul
CFLAGSOPTF90 := -march=core -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpclmul
CFLAGSOPTC   := -march=core -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpclmul
CFLAGSOPTCXX := -march=core -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpclmul
LDFLAGSOPT   := -march=core -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpclmul
endif

# Intel Atom CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-atom-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=anyx86 -msse -msse2 -msse3 -mssse3
CFLAGSOPTF90 := -march=anyx86 -msse -msse2 -msse3 -mssse3
CFLAGSOPTC   := -march=anyx86 -msse -msse2 -msse3 -mssse3
CFLAGSOPTCXX := -march=anyx86 -msse -msse2 -msse3 -mssse3
LDFLAGSOPT   := -march=anyx86 -msse -msse2 -msse3 -mssse3
endif

##############################################################################
# AMD x86/x86_64 CPUs
#
##############################################################################

# AMD K6 CPU with MMX instruction set support and improved versions of
# AMD K6 CPU with MMX and 3DNow! instruction set support,
# respectively.
ifeq ($(call match,$(ID),pc-k6-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=anyx86 -m3dnow
CFLAGSOPTF90 := -march=anyx86 -m3dnow
CFLAGSOPTC   := -march=anyx86 -m3dnow
CFLAGSOPTCXX := -march=anyx86 -m3dnow
LDFLAGSOPT   := -march=anyx86 -m3dnow
endif

# AMD Athlon CPU with MMX, 3dNOW!, enhanced 3DNow! and SSE prefetch
# instructions support.
ifeq ($(call match,$(ID),pc-athlon-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=athlon
CFLAGSOPTF90 := -march=athlon
CFLAGSOPTC   := -march=athlon
CFLAGSOPTCXX := -march=athlon
LDFLAGSOPT   := -march=athlon
endif

# Improved AMD Athlon CPU with MMX, 3DNow!, enhanced 3DNow! and full
# SSE instruction set support.
ifeq ($(call match,$(ID),pc-athlonxp-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=athlon
CFLAGSOPTF90 := -march=athlon
CFLAGSOPTC   := -march=athlon
CFLAGSOPTCXX := -march=athlon
LDFLAGSOPT   := -march=athlon
endif

# Processors based on the AMD K8 core with x86-64 instruction set
# support, including the AMD Opteron, Athlon 64, and Athlon 64 FX
# processors. (This supersets MMX, SSE, SSE2, 3DNow!, enhanced 3DNow!
# and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(athlon64|athlon64x2|turion64|turion64x2)-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=athlon64
CFLAGSOPTF90 := -march=athlon64
CFLAGSOPTC   := -march=athlon64
CFLAGSOPTCXX := -march=athlon64
LDFLAGSOPT   := -march=athlon64
endif

#Processors based on the AMD K8 core with x86-64 instruction set
#support, including the AMD Opteron, Athlon 64, and Athlon 64 FX
#processors. (This supersets MMX, SSE, SSE2, 3DNow!, enhanced 3DNow!
#and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(opteron|opteronx2)-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=athlon64
CFLAGSOPTF90 := -march=athlon64
CFLAGSOPTC   := -march=athlon64
CFLAGSOPTCXX := -march=athlon64
LDFLAGSOPT   := -march=athlon64
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-phenom-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=barcelona
CFLAGSOPTF90 := -march=barcelona
CFLAGSOPTC   := -march=barcelona
CFLAGSOPTCXX := -march=barcelona
LDFLAGSOPT   := -march=barcelona
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-phenomII-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=barcelona
CFLAGSOPTF90 := -march=barcelona
CFLAGSOPTC   := -march=barcelona
CFLAGSOPTCXX := -march=barcelona
LDFLAGSOPT   := -march=barcelona
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(opteron4100|opteron6100)-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=barcelona
CFLAGSOPTF90 := -march=barcelona
CFLAGSOPTC   := -march=barcelona
CFLAGSOPTCXX := -march=barcelona
LDFLAGSOPT   := -march=barcelona
endif

# AMD Family 15: Bulldozer core
# CPUs based on AMD Family 15h cores with x86-64 instruction set
# support. (This supersets FMA4, AVX, XOP, LWP, AES, PCL_MUL, CX16,
# MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2, ABM and 64-bit
# instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-bulldozer-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=bdver1
CFLAGSOPTF90 := -march=bdver1
CFLAGSOPTC   := -march=bdver1
CFLAGSOPTCXX := -march=bdver1
LDFLAGSOPT   := -march=bdver1
endif

# AMD Family 15h: Piledriver core
# AMD Family 15h core based CPUs with x86-64 instruction set
# support. (This supersets BMI, TBM, F16C, FMA, AVX, XOP, LWP, AES,
# PCL_MUL, CX16, MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2,
# ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-piledriver-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=bdver2
CFLAGSOPTF90 := -march=bdver2
CFLAGSOPTC   := -march=bdver2
CFLAGSOPTCXX := -march=bdver2
LDFLAGSOPT   := -march=bdver2
endif

# AMD Family 15h: Steamroller core
# AMD Family 15h core based CPUs with x86-64 instruction set
# support. (This supersets BMI, TBM, F16C, FMA, AVX, XOP, LWP, AES,
# PCL_MUL, CX16, MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2,
# ABM and 64-bit instruction set extensions.
ifeq ($(call match,$(ID),(pc|pc64)-steamroller-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=bdver2
CFLAGSOPTF90 := -march=bdver2
CFLAGSOPTC   := -march=bdver2
CFLAGSOPTCXX := -march=bdver2
LDFLAGSOPT   := -march=bdver2
endif

# AMD Family 14h: Bobcat core
# CPUs based on AMD Family 14h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSSE3, SSE4A, CX16,
# ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-bobcat-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=bdver2
CFLAGSOPTF90 := -march=bdver2
CFLAGSOPTC   := -march=bdver2
CFLAGSOPTCXX := -march=bdver2
LDFLAGSOPT   := -march=bdver2
endif


# AMD Family 16h: Jaguar core
# CPUs based on AMD Family 16h cores with x86-64 instruction set
# support. This includes MOVBE, F16C, BMI, AVX, PCL_MUL, AES, SSE4.2,
# SSE4.1, CX16, ABM, SSE4A, SSSE3, SSE3, SSE2, SSE, MMX and 64-bit
# instruction set extensions.
ifeq ($(call match,$(ID),(pc|pc64)-jaguar-.*-open64-.*),yes)
CFLAGSOPTF77 := -march=bdver2
CFLAGSOPTF90 := -march=bdver2
CFLAGSOPTC   := -march=bdver2
CFLAGSOPTCXX := -march=bdver2
LDFLAGSOPT   := -march=bdver2
endif

##############################################################################
# Automatic CPU type detection
#
##############################################################################

ifneq (,$(findstring NATIVE,$(OPT)))
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by Open64 compiler!';
CFLAGSOPTF77 := -march=auto
CFLAGSOPTF90 := -march=auto
CFLAGSOPTC   := -march=auto
CFLAGSOPTCXX := -march=auto
LDFLAGSOPT   := -march=auto
endif
