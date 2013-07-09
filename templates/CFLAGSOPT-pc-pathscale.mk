# -*- mode: makefile -*-

##############################################################################
# PathScale Compiler suite with individual settings for x86/x86_64 CPU IDs
#
##############################################################################

##############################################################################
# Intel x86/x86_64 CPUs
#
##############################################################################

# Original Intel i386 CPU.
ifeq ($(call match,$(ID),pc-i386-.*-psc-.*),yes)
CFLAGSOPTF77 := 
CFLAGSOPTF90 := 
CFLAGSOPTC   := 
CFLAGSOPTCXX := 
LDFLAGSOPT   := 
endif

# Intel i486 CPU. (No scheduling is implemented for this chip.)
ifeq ($(call match,$(ID),pc-i486-.*-psc-.*),yes)
CFLAGSOPTF77 := 
CFLAGSOPTF90 := 
CFLAGSOPTC   := 
CFLAGSOPTCXX := 
LDFLAGSOPT   := 
endif

# Intel Pentium CPU with or without MMX support.
ifeq ($(call match,$(ID),pc-(i568|pentium)-.*-psc-.*),yes)
CFLAGSOPTF77 := 
CFLAGSOPTF90 := 
CFLAGSOPTC   := 
CFLAGSOPTCXX := 
LDFLAGSOPT   := 
endif

# Intel Pentium Pro CPU.
ifeq ($(call match,$(ID),pc-pentiumpro-.*-psc-.*),yes)
CFLAGSOPTF77 := 
CFLAGSOPTF90 := 
CFLAGSOPTC   := 
CFLAGSOPTCXX := 
LDFLAGSOPT   := 
endif

# Intel Pentium II CPU, based on Pentium Pro core with MMX instruction
# set support.
ifeq ($(call match,$(ID),pc-pentium2-.*-psc-.*),yes)
CFLAGSOPTF77 := 
CFLAGSOPTF90 := 
CFLAGSOPTC   := 
CFLAGSOPTCXX := 
LDFLAGSOPT   := 
endif

# Intel Pentium III CPU, based on Pentium Pro core with MMX and SSE
# instruction set support.
ifeq ($(call match,$(ID),pc-pentium3-.*-psc-.*),yes)
CFLAGSOPTF77  := 
CFLAGSOPTF90  := 
CFLAGSOPTC    := 
CFLAGSOPTCXX  := 
LDFLAGSOPT    := 
endif

# Intel Pentium M; low-power version of Intel Pentium III CPU with MMX, 
# SSE and SSE2 instruction set support. Used by Centrino notebooks.
ifeq ($(call match,$(ID),pc-pentiumm-.*-psc-.*),yes)
CFLAGSOPTF77  := -msse2
CFLAGSOPTF90  := -msse2
CFLAGSOPTC    := -msse2
CFLAGSOPTCXX  := -msse2
LDFLAGSOPT    := -msse2
endif

# Intel Pentium 4 CPU with MMX, SSE and SSE2 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-pentium4-.*-psc-.*),yes)
ifeq ($(call match,$(core),"prescott"),yes)
CFLAGSOPTF77  := -march=prescott
CFLAGSOPTF90  := -march=prescott
CFLAGSOPTC    := -march=prescott
CFLAGSOPTCXX  := -march=prescott
LDFLAGSOPT    := -march=prescott
else
ifeq ($(call match,$(core),"nocona"),yes)
CFLAGSOPTF77  := -march=nocona
CFLAGSOPTF90  := -march=nocona
CFLAGSOPTC    := -march=nocona
CFLAGSOPTCXX  := -march=nocona
LDFLAGSOPT    := -march=nocona
else
CFLAGSOPTF77  := -msse2
CFLAGSOPTF90  := -msse2
CFLAGSOPTC    := -msse2
CFLAGSOPTCXX  := -msse2
LDFLAGSOPT    := -msse2
endif
endif
endif

# Intel CoreSolo/Duo. Improved version of Intel Pentium 4 CPU with MMX,
# SSE, SSE2 and SSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-coresolo-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=core
CFLAGSOPTF90  := -march=core
CFLAGSOPTC    := -march=core
CFLAGSOPTCXX  := -march=core
LDFLAGSOPT    := -march=core
endif

# Intel Core 2 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-coreduo-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=core2
CFLAGSOPTF90  := -march=core2
CFLAGSOPTC    := -march=core2
CFLAGSOPTCXX  := -march=core2
LDFLAGSOPT    := -march=core2
endif

# Intel Core 2 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-penryn-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=core2
CFLAGSOPTF90  := -march=core2
CFLAGSOPTC    := -march=core2
CFLAGSOPTCXX  := -march=core2
LDFLAGSOPT    := -march=core2
endif

# Intel Nehalem CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-nehalem-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=nehalem
CFLAGSOPTF90  := -march=nehalem
CFLAGSOPTC    := -march=nehalem
CFLAGSOPTCXX  := -march=nehalem
LDFLAGSOPT    := -march=nehalem
endif

# Intel Core i7 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3,
# SSSE3, SSE4.1, SSE4.2, AVX, AES and PCLMUL instruction set
# support.
ifeq ($(call match,$(ID),(pc|pc64)-sandybridge-.*-psc-.*),yes)
CFLAGSOPTF77 := -march=sandy
CFLAGSOPTF90 := -march=sandy
CFLAGSOPTC   := -march=sandy
CFLAGSOPTCXX := -march=sandy
LDFLAGSOPT   := -march=sandy
endif

# Intel Core CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3, SSSE3,
# SSE4.1, SSE4.2, AVX, AES, PCLMUL, FSGSBASE, RDRND and F16C
# instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-ivybridge-.*-psc-.*),yes)
CFLAGSOPTF77 := -march=sandy
CFLAGSOPTF90 := -march=sandy
CFLAGSOPTC   := -march=sandy
CFLAGSOPTCXX := -march=sandy
LDFLAGSOPT   := -march=sandy
endif

# Intel Core CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3, SSSE3,
# SSE4.1, SSE4.2, AVX, AES, PCLMUL, FSGSBASE, RDRND and F16C
# instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-haswell-.*-psc-.*),yes)
CFLAGSOPTF77 := -march=sandy
CFLAGSOPTF90 := -march=sandy
CFLAGSOPTC   := -march=sandy
CFLAGSOPTCXX := -march=sandy
LDFLAGSOPT   := -march=sandy
endif

# Intel Atom CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-atom-.*-psc-.*),yes)
CFLAGSOPTF77 := -march=core2
CFLAGSOPTF90 := -march=core2
CFLAGSOPTC   := -march=core2
CFLAGSOPTCXX := -march=core2
LDFLAGSOPT   := -march=core2
endif

##############################################################################
# AMD x86/x86_64 CPUs
#
##############################################################################

# AMD K6 CPU with MMX instruction set support and improved versions of
# AMD K6 CPU with MMX and 3DNow! instruction set support,
# respectively.
ifeq ($(call match,$(ID),pc-k6-.*-psc-.*),yes)
CFLAGSOPTF77 := 
CFLAGSOPTF90 := 
CFLAGSOPTC   := 
CFLAGSOPTCXX := 
LDFLAGSOPT   := 
endif

# AMD Athlon CPU with MMX, 3dNOW!, enhanced 3DNow! and SSE prefetch
# instructions support.
ifeq ($(call match,$(ID),pc-athlon-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=athlon
CFLAGSOPTF90  := -march=athlon
CFLAGSOPTC    := -march=athlon
CFLAGSOPTCXX  := -march=athlon
LDFLAGSOPT    := -march=athlon
endif

# Improved AMD Athlon CPU with MMX, 3DNow!, enhanced 3DNow! and full
# SSE instruction set support.
ifeq ($(call match,$(ID),pc-athlonxp-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=athlon
CFLAGSOPTF90  := -march=athlon
CFLAGSOPTC    := -march=athlon
CFLAGSOPTCXX  := -march=athlon
LDFLAGSOPT    := -march=athlon
endif

# Processors based on the AMD K8 core with x86-64 instruction set
# support, including the AMD Opteron, Athlon 64, and Athlon 64 FX
# processors. (This supersets MMX, SSE, SSE2, 3DNow!, enhanced 3DNow!
# and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(athlon64|athlon64x2|turion64|turion64x2)-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=athlon
CFLAGSOPTF90  := -march=athlon
CFLAGSOPTC    := -march=athlon
CFLAGSOPTCXX  := -march=athlon
LDFLAGSOPT    := -march=athlon
endif

# Processors based on the AMD K8 core with x86-64 instruction set
# support, including the AMD Opteron, Athlon 64, and Athlon 64 FX
# processors. (This supersets MMX, SSE, SSE2, 3DNow!, enhanced 3DNow!
# and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(opteron|opteronx2)-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=opteron
CFLAGSOPTF90  := -march=opteron
CFLAGSOPTC    := -march=opteron
CFLAGSOPTCXX  := -march=opteron
LDFLAGSOPT    := -march=opteron
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-phenom-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=barcelona
CFLAGSOPTF90  := -march=barcelona
CFLAGSOPTC    := -march=barcelona
CFLAGSOPTCXX  := -march=barcelona
LDFLAGSOPT    := -march=barcelona
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-phenomII-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=barcelona
CFLAGSOPTF90  := -march=barcelona
CFLAGSOPTC    := -march=barcelona
CFLAGSOPTCXX  := -march=barcelona
LDFLAGSOPT    := -march=barcelona
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(opteron4100|opteron6100)-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=barcelona
CFLAGSOPTF90  := -march=barcelona
CFLAGSOPTC    := -march=barcelona
CFLAGSOPTCXX  := -march=barcelona
LDFLAGSOPT    := -march=barcelona
endif

# AMD Family 15: Bulldozer core
# CPUs based on AMD Family 15h cores with x86-64 instruction set
# support. (This supersets FMA4, AVX, XOP, LWP, AES, PCL_MUL, CX16,
# MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2, ABM and 64-bit
# instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-bulldozer-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=bdver1
CFLAGSOPTF90  := -march=bdver1
CFLAGSOPTC    := -march=bdver1
CFLAGSOPTCXX  := -march=bdver1
LDFLAGSOPT    := -march=bdver1
endif

# AMD Family 15h: Piledriver core
# AMD Family 15h core based CPUs with x86-64 instruction set
# support. (This supersets BMI, TBM, F16C, FMA, AVX, XOP, LWP, AES,
# PCL_MUL, CX16, MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2,
# ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-piledriver-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=bdver1
CFLAGSOPTF90  := -march=bdver1
CFLAGSOPTC    := -march=bdver1
CFLAGSOPTCXX  := -march=bdver1
LDFLAGSOPT    := -march=bdver1
endif

# AMD Family 15h: Steamroller core
# AMD Family 15h core based CPUs with x86-64 instruction set
# support. (This supersets BMI, TBM, F16C, FMA, AVX, XOP, LWP, AES,
# PCL_MUL, CX16, MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2,
# ABM and 64-bit instruction set extensions.
ifeq ($(call match,$(ID),(pc|pc64)-steamroller-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=bdver1
CFLAGSOPTF90  := -march=bdver1
CFLAGSOPTC    := -march=bdver1
CFLAGSOPTCXX  := -march=bdver1
LDFLAGSOPT    := -march=bdver1
endif

# AMD Family 14h: Bobcat core
# CPUs based on AMD Family 14h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSSE3, SSE4A, CX16,
# ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-bobcat-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=bdver1
CFLAGSOPTF90  := -march=bdver1
CFLAGSOPTC    := -march=bdver1
CFLAGSOPTCXX  := -march=bdver1
LDFLAGSOPT    := -march=bdver1
endif

# AMD Family 16h: Jaguar core
# CPUs based on AMD Family 16h cores with x86-64 instruction set
# support. This includes MOVBE, F16C, BMI, AVX, PCL_MUL, AES, SSE4.2,
# SSE4.1, CX16, ABM, SSE4A, SSSE3, SSE3, SSE2, SSE, MMX and 64-bit
# instruction set extensions.
ifeq ($(call match,$(ID),(pc|pc64)-jaguar-.*-psc-.*),yes)
CFLAGSOPTF77  := -march=bdver1
CFLAGSOPTF90  := -march=bdver1
CFLAGSOPTC    := -march=bdver1
CFLAGSOPTCXX  := -march=bdver1
LDFLAGSOPT    := -march=bdver1
endif

##############################################################################
# Automatic CPU type detection
#
##############################################################################

ifneq (,$(findstring NATIVE,$(OPT)))
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by Pathscale compiler!';
CFLAGSOPTF77 := -march=auto
CFLAGSOPTF90 := -march=auto
CFLAGSOPTC   := -march=auto
CFLAGSOPTCXX := -march=auto
LDFLAGSOPT   := -march=auto
endif
