# -*- mode: makefile -*-

##############################################################################
# PGI Compiler suite with individual settings for x86/x86_64 CPU IDs
#
##############################################################################

##############################################################################
# Intel x86/x86_64 CPUs
#
##############################################################################

# Original Intel i386 CPU.
ifeq ($(call match,$(ID),pc-i386-.*-pgi-.*),yes)
CFLAGSOPTF77 := -tp px
CFLAGSOPTF90 := -tp px
CFLAGSOPTC   := -tp px
CFLAGSOPTCXX := -tp px
LDFLAGSOPT   := -tp px
endif

# Intel i486 CPU. (No scheduling is implemented for this chip.)
ifeq ($(call match,$(ID),pc-i486-.*-pgi-.*),yes)
CFLAGSOPTF77 := -tp px
CFLAGSOPTF90 := -tp px
CFLAGSOPTC   := -tp px
CFLAGSOPTCXX := -tp px
LDFLAGSOPT   := -tp px
endif

# Intel Pentium CPU with or without MMX support.
ifeq ($(call match,$(ID),pc-(i568|pentium)-.*-pgi-.*),yes)
CFLAGSOPTF77 := -tp p5
CFLAGSOPTF90 := -tp p5
CFLAGSOPTC   := -tp p5
CFLAGSOPTCXX := -tp p5
LDFLAGSOPT   := -tp p5
endif

# Intel Pentium Pro CPU.
ifeq ($(call match,$(ID),pc-pentiumpro-.*-pgi-.*),yes)
CFLAGSOPTF77 := -tp p6
CFLAGSOPTF90 := -tp p6
CFLAGSOPTC   := -tp p6
CFLAGSOPTCXX := -tp p6
LDFLAGSOPT   := -tp p6
endif

# Intel Pentium II CPU, based on Pentium Pro core with MMX instruction
# set support.
ifeq ($(call match,$(ID),pc-pentium2-.*-pgi-.*),yes)
CFLAGSOPTF77 := -tp p6
CFLAGSOPTF90 := -tp p6
CFLAGSOPTC   := -tp p6
CFLAGSOPTCXX := -tp p6
LDFLAGSOPT   := -tp p6
endif

# Intel Pentium III CPU, based on Pentium Pro core with MMX and SSE
# instruction set support.
ifeq ($(call match,$(ID),pc-pentium3-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp piii
CFLAGSOPTF90  := -tp piii
CFLAGSOPTC    := -tp piii
CFLAGSOPTCXX  := -tp piii
LDFLAGSOPT    := -tp piii
endif

# Intel Pentium M; low-power version of Intel Pentium III CPU with MMX, 
# SSE and SSE2 instruction set support. Used by Centrino notebooks.
ifeq ($(call match,$(ID),pc-pentiumm-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp piii
CFLAGSOPTF90  := -tp piii
CFLAGSOPTC    := -tp piii
CFLAGSOPTCXX  := -tp piii
LDFLAGSOPT    := -tp piii
endif

# Intel Pentium 4 CPU with MMX, SSE and SSE2 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-pentium4-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp pvi
CFLAGSOPTF90  := -tp pvi
CFLAGSOPTC    := -tp pvi
CFLAGSOPTCXX  := -tp pvi
LDFLAGSOPT    := -tp pvi
endif

# Intel CoreSolo/Duo. Improved version of Intel Pentium 4 CPU with MMX,
# SSE, SSE2 and SSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-coresolo-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp pvi
CFLAGSOPTF90  := -tp pvi
CFLAGSOPTC    := -tp pvi
CFLAGSOPTCXX  := -tp pvi
LDFLAGSOPT    := -tp pvi
endif

# Intel Core 2 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-coreduo-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp core2
CFLAGSOPTF90  := -tp core2
CFLAGSOPTC    := -tp core2
CFLAGSOPTCXX  := -tp core2
LDFLAGSOPT    := -tp core2
endif

# Intel Core 2 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-penryn-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp penryn
CFLAGSOPTF90  := -tp penryn
CFLAGSOPTC    := -tp penryn
CFLAGSOPTCXX  := -tp penryn
LDFLAGSOPT    := -tp penryn
endif

# Intel Nehalem CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-nehalem-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp nehalem
CFLAGSOPTF90  := -tp nehalem
CFLAGSOPTC    := -tp nehalem
CFLAGSOPTCXX  := -tp nehalem
LDFLAGSOPT    := -tp nehalem
endif

# Intel Core i7 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3,
# SSSE3, SSE4.1, SSE4.2, AVX, AES and PCLMUL instruction set
# support.
ifeq ($(call match,$(ID),(pc|pc64)-sandybridge-.*-pgi-.*),yes)
CFLAGSOPTF77 := -tp sandybridge
CFLAGSOPTF90 := -tp sandybridge
CFLAGSOPTC   := -tp sandybridge
CFLAGSOPTCXX := -tp sandybridge
LDFLAGSOPT   := -tp sandybridge
endif

# Intel Core CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3, SSSE3,
# SSE4.1, SSE4.2, AVX, AES, PCLMUL, FSGSBASE, RDRND and F16C
# instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-ivybridge-.*-pgi-.*),yes)
CFLAGSOPTF77 := -tp ivybridge
CFLAGSOPTF90 := -tp ivybridge
CFLAGSOPTC   := -tp ivybridge
CFLAGSOPTCXX := -tp ivybridge
LDFLAGSOPT   := -tp ivybridge
endif

# Intel Core CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3, SSSE3,
# SSE4.1, SSE4.2, AVX, AES, PCLMUL, FSGSBASE, RDRND and F16C
# instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-haswell-.*-pgi-.*),yes)
CFLAGSOPTF77 := -tp ivybridge
CFLAGSOPTF90 := -tp ivybridge
CFLAGSOPTC   := -tp ivybridge
CFLAGSOPTCXX := -tp ivybridge
LDFLAGSOPT   := -tp ivybridge
endif

# Intel Atom CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-atom-.*-pgi-.*),yes)
CFLAGSOPTF77 := -tp core2
CFLAGSOPTF90 := -tp core2
CFLAGSOPTC   := -tp core2
CFLAGSOPTCXX := -tp core2
LDFLAGSOPT   := -tp core2
endif

##############################################################################
# AMD x86/x86_64 CPUs
#
##############################################################################

# AMD K6 CPU with MMX instruction set support and improved versions of
# AMD K6 CPU with MMX and 3DNow! instruction set support,
# respectively.
ifeq ($(call match,$(ID),pc-k6-.*-pgi-.*),yes)
CFLAGSOPTF77 := -tp px
CFLAGSOPTF90 := -tp px
CFLAGSOPTC   := -tp px
CFLAGSOPTCXX := -tp px
LDFLAGSOPT   := -tp px
endif

# AMD Athlon CPU with MMX, 3dNOW!, enhanced 3DNow! and SSE prefetch
# instructions support.
ifeq ($(call match,$(ID),pc-athlon-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp athlon
CFLAGSOPTF90  := -tp athlon
CFLAGSOPTC    := -tp athlon
CFLAGSOPTCXX  := -tp athlon
LDFLAGSOPT    := -tp athlon
endif

# Improved AMD Athlon CPU with MMX, 3DNow!, enhanced 3DNow! and full
# SSE instruction set support.
ifeq ($(call match,$(ID),pc-athlonxp-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp athlonxp
CFLAGSOPTF90  := -tp athlonxp
CFLAGSOPTC    := -tp athlonxp
CFLAGSOPTCXX  := -tp athlonxp
LDFLAGSOPT    := -tp athlonxp
endif

# Processors based on the AMD K8 core with x86-64 instruction set
# support, including the AMD Opteron, Athlon 64, and Athlon 64 FX
# processors. (This supersets MMX, SSE, SSE2, 3DNow!, enhanced 3DNow!
# and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(athlon64|athlon64x2|turion64|turion64x2)-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp k8
CFLAGSOPTF90  := -tp k8
CFLAGSOPTC    := -tp k8
CFLAGSOPTCXX  := -tp k8
LDFLAGSOPT    := -tp k8
endif

# Processors based on the AMD K8 core with x86-64 instruction set
# support, including the AMD Opteron, Athlon 64, and Athlon 64 FX
# processors. (This supersets MMX, SSE, SSE2, 3DNow!, enhanced 3DNow!
# and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(opteron|opteronx2)-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp k8
CFLAGSOPTF90  := -tp k8
CFLAGSOPTC    := -tp k8
CFLAGSOPTCXX  := -tp k8
LDFLAGSOPT    := -tp k8
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-phenom-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp barcelona
CFLAGSOPTF90  := -tp barcelona
CFLAGSOPTC    := -tp barcelona
CFLAGSOPTCXX  := -tp barcelona
LDFLAGSOPT    := -tp barcelona
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-phenomII-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp shanghai
CFLAGSOPTF90  := -tp shanghai
CFLAGSOPTC    := -tp shanghai
CFLAGSOPTCXX  := -tp shanghai
LDFLAGSOPT    := -tp shanghai
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(opteron4100|opteron6100)-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp istanbul
CFLAGSOPTF90  := -tp istanbul
CFLAGSOPTC    := -tp istanbul
CFLAGSOPTCXX  := -tp istanbul
LDFLAGSOPT    := -tp istanbul
endif

# AMD Family 15: Bulldozer core
# CPUs based on AMD Family 15h cores with x86-64 instruction set
# support. (This supersets FMA4, AVX, XOP, LWP, AES, PCL_MUL, CX16,
# MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2, ABM and 64-bit
# instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-bulldozer-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp bulldozer
CFLAGSOPTF90  := -tp bulldozer
CFLAGSOPTC    := -tp bulldozer
CFLAGSOPTCXX  := -tp bulldozer
LDFLAGSOPT    := -tp bulldozer
endif

# AMD Family 15h: Piledriver core
# AMD Family 15h core based CPUs with x86-64 instruction set
# support. (This supersets BMI, TBM, F16C, FMA, AVX, XOP, LWP, AES,
# PCL_MUL, CX16, MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2,
# ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-piledriver-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp bulldozer
CFLAGSOPTF90  := -tp bulldozer
CFLAGSOPTC    := -tp bulldozer
CFLAGSOPTCXX  := -tp bulldozer
LDFLAGSOPT    := -tp bulldozer
endif

# AMD Family 15h: Steamroller core
# AMD Family 15h core based CPUs with x86-64 instruction set
# support. (This supersets BMI, TBM, F16C, FMA, AVX, XOP, LWP, AES,
# PCL_MUL, CX16, MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2,
# ABM and 64-bit instruction set extensions.
ifeq ($(call match,$(ID),(pc|pc64)-steamroller-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp bulldozer
CFLAGSOPTF90  := -tp bulldozer
CFLAGSOPTC    := -tp bulldozer
CFLAGSOPTCXX  := -tp bulldozer
LDFLAGSOPT    := -tp bulldozer
endif

# AMD Family 14h: Bobcat core
# CPUs based on AMD Family 14h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSSE3, SSE4A, CX16,
# ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-bobcat-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp bulldozer
CFLAGSOPTF90  := -tp bulldozer
CFLAGSOPTC    := -tp bulldozer
CFLAGSOPTCXX  := -tp bulldozer
LDFLAGSOPT    := -tp bulldozer
endif

# AMD Family 16h: Jaguar core
# CPUs based on AMD Family 16h cores with x86-64 instruction set
# support. This includes MOVBE, F16C, BMI, AVX, PCL_MUL, AES, SSE4.2,
# SSE4.1, CX16, ABM, SSE4A, SSSE3, SSE3, SSE2, SSE, MMX and 64-bit
# instruction set extensions.
ifeq ($(call match,$(ID),(pc|pc64)-jaguar-.*-pgi-.*),yes)
CFLAGSOPTF77  := -tp bulldozer
CFLAGSOPTF90  := -tp bulldozer
CFLAGSOPTC    := -tp bulldozer
CFLAGSOPTCXX  := -tp bulldozer
LDFLAGSOPT    := -tp bulldozer
endif

##############################################################################
# Automatic CPU type detection
#
##############################################################################

ifneq (,$(findstring NATIVE,$(OPT)))
MESSAGE  := $(MESSAGE) \
	    echo '*** Message: Used auto-tuning by PGI compiler!';
CFLAGSOPTF77 := 
CFLAGSOPTF90 := 
CFLAGSOPTC   := 
CFLAGSOPTCXX := 
LDFLAGSOPT   := 
endif
