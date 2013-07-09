# -*- mode: makefile -*-

##############################################################################
# SunStudio Compiler suite with individual settings for x86/x86_64 CPU IDs
#
##############################################################################

# Option -xtarget=cpu-type specifies that the code is compiled with
# optimization for the specified machine type "cpu-type".

##############################################################################
# Intel x86/x86_64 CPUs
#
##############################################################################

# List of supported cpu-types and minimally required compiler version
# CPU-TYPE    min. SunStudio version
# generic     2009
# native      2009
# pentium     2009
# pentium_pro 2009
# pentium3    2009
# pentium4    2009
# woodcrest   2009
# penryn      2009
# nehalem     2009
# sandybridge 12.3
# westmere    2009

# Original Intel i386 CPU.
ifeq ($(call match,$(ID),pc-i386-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=generic
CFLAGSOPTF90 := -xtarget=generic
CFLAGSOPTC   := -xtarget=generic
CFLAGSOPTCXX := -xtarget=generic
LDFLAGSOPT   := -xtarget=generic
endif

# Intel i486 CPU. (No scheduling is implemented for this chip.)
ifeq ($(call match,$(ID),pc-i486-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=generic
CFLAGSOPTF90 := -xtarget=generic
CFLAGSOPTC   := -xtarget=generic
CFLAGSOPTCXX := -xtarget=generic
LDFLAGSOPT   := -xtarget=generic
endif

# Intel Pentium CPU with or without MMX support.
ifeq ($(call match,$(ID),pc-(i568|pentium)-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=pentium
CFLAGSOPTF90 := -xtarget=pentium
CFLAGSOPTC   := -xtarget=pentium
CFLAGSOPTCXX := -xtarget=pentium
LDFLAGSOPT   := -xtarget=pentium
endif

# Intel Pentium Pro CPU.
ifeq ($(call match,$(ID),pc-pentiumpro-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=pentium_pro
CFLAGSOPTF90 := -xtarget=pentium_pro
CFLAGSOPTC   := -xtarget=pentium_pro
CFLAGSOPTCXX := -xtarget=pentium_pro
LDFLAGSOPT   := -xtarget=pentium_pro
endif

# Intel Pentium II CPU, based on Pentium Pro core with MMX instruction
# set support.
ifeq ($(call match,$(ID),pc-pentium2-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=pentium_pro
CFLAGSOPTF90 := -xtarget=pentium_pro
CFLAGSOPTC   := -xtarget=pentium_pro
CFLAGSOPTCXX := -xtarget=pentium_pro
LDFLAGSOPT   := -xtarget=pentium_pro
endif

# Intel Pentium III CPU, based on Pentium Pro core with MMX and SSE
# instruction set support.
ifeq ($(call match,$(ID),pc-pentium3-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=pentium3
CFLAGSOPTF90  := -xtarget=pentium3
CFLAGSOPTC    := -xtarget=pentium3
CFLAGSOPTCXX  := -xtarget=pentium3
LDFLAGSOPT    := -xtarget=pentium3
endif

# Intel Pentium M; low-power version of Intel Pentium III CPU with MMX, 
# SSE and SSE2 instruction set support. Used by Centrino notebooks.
ifeq ($(call match,$(ID),pc-pentiumm-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=pentium3
CFLAGSOPTF90  := -xtarget=pentium3
CFLAGSOPTC    := -xtarget=pentium3
CFLAGSOPTCXX  := -xtarget=pentium3
LDFLAGSOPT    := -xtarget=pentium3
endif

# Intel Pentium 4 CPU with MMX, SSE and SSE2 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-pentium4-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=pentium4
CFLAGSOPTF90  := -xtarget=pentium4
CFLAGSOPTC    := -xtarget=pentium4
CFLAGSOPTCXX  := -xtarget=pentium4
LDFLAGSOPT    := -xtarget=pentium4
endif

# Intel CoreSolo/Duo. Improved version of Intel Pentium 4 CPU with MMX,
# SSE, SSE2 and SSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-coresolo-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=pentium4
CFLAGSOPTF90  := -xtarget=pentium4
CFLAGSOPTC    := -xtarget=pentium4
CFLAGSOPTCXX  := -xtarget=pentium4
LDFLAGSOPT    := -xtarget=pentium4
endif

# Intel Core 2 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-coreduo-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=woodcrest
CFLAGSOPTF90  := -xtarget=woodcrest
CFLAGSOPTC    := -xtarget=woodcrest
CFLAGSOPTCXX  := -xtarget=woodcrest
LDFLAGSOPT    := -xtarget=woodcrest
endif

# Intel Core 2 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-penryn-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=penryn
CFLAGSOPTF90  := -xtarget=penryn
CFLAGSOPTC    := -xtarget=penryn
CFLAGSOPTCXX  := -xtarget=penryn
LDFLAGSOPT    := -xtarget=penryn
endif

# Intel Nehalem CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-nehalem-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=nehalem
CFLAGSOPTF90  := -xtarget=nehalem
CFLAGSOPTC    := -xtarget=nehalem
CFLAGSOPTCXX  := -xtarget=nehalem
LDFLAGSOPT    := -xtarget=nehalem
endif

# Intel Core i7 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3,
# SSSE3, SSE4.1, SSE4.2, AVX, AES and PCLMUL instruction set
# support.
ifeq ($(call match,$(ID),(pc|pc64)-sandybridge-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=sandybridge
CFLAGSOPTF90 := -xtarget=sandybridge
CFLAGSOPTC   := -xtarget=sandybridge
CFLAGSOPTCXX := -xtarget=sandybridge
LDFLAGSOPT   := -xtarget=sandybridge
endif

# Intel Core CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3, SSSE3,
# SSE4.1, SSE4.2, AVX, AES, PCLMUL, FSGSBASE, RDRND and F16C
# instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-ivybridge-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=sandybridge
CFLAGSOPTF90 := -xtarget=sandybridge
CFLAGSOPTC   := -xtarget=sandybridge
CFLAGSOPTCXX := -xtarget=sandybridge
LDFLAGSOPT   := -xtarget=sandybridge
endif

# Intel Core CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3, SSSE3,
# SSE4.1, SSE4.2, AVX, AES, PCLMUL, FSGSBASE, RDRND and F16C
# instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-haswell-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=sandybridge
CFLAGSOPTF90 := -xtarget=sandybridge
CFLAGSOPTC   := -xtarget=sandybridge
CFLAGSOPTCXX := -xtarget=sandybridge
LDFLAGSOPT   := -xtarget=sandybridge
endif

# Intel Atom CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-atom-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=sandybridge
CFLAGSOPTF90 := -xtarget=sandybridge
CFLAGSOPTC   := -xtarget=sandybridge
CFLAGSOPTCXX := -xtarget=sandybridge
LDFLAGSOPT   := -xtarget=sandybridge
endif

##############################################################################
# AMD x86/x86_64 CPUs
#
##############################################################################

# List of supported cpu-types and minimally required compiler version
# CPU-TYPE    min. SunStudio version
# generic     2009
# native      2009
# barcelone   2009
# opteron     2009

# AMD K6 CPU with MMX instruction set support and improved versions of
# AMD K6 CPU with MMX and 3DNow! instruction set support,
# respectively.
ifeq ($(call match,$(ID),pc-k6-.*-sunstudio-.*),yes)
CFLAGSOPTF77 := -xtarget=generic
CFLAGSOPTF90 := -xtarget=generic
CFLAGSOPTC   := -xtarget=generic
CFLAGSOPTCXX := -xtarget=generic
LDFLAGSOPT   := -xtarget=generic
endif

# AMD Athlon CPU with MMX, 3dNOW!, enhanced 3DNow! and SSE prefetch
# instructions support.
ifeq ($(call match,$(ID),pc-athlon-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=generic
CFLAGSOPTF90  := -xtarget=generic
CFLAGSOPTC    := -xtarget=generic
CFLAGSOPTCXX  := -xtarget=generic
LDFLAGSOPT    := -xtarget=generic
endif

# Improved AMD Athlon CPU with MMX, 3DNow!, enhanced 3DNow! and full
# SSE instruction set support.
ifeq ($(call match,$(ID),pc-athlonxp-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=generic
CFLAGSOPTF90  := -xtarget=generic
CFLAGSOPTC    := -xtarget=generic
CFLAGSOPTCXX  := -xtarget=generic
LDFLAGSOPT    := -xtarget=generic
endif

# Processors based on the AMD K8 core with x86-64 instruction set
# support, including the AMD Opteron, Athlon 64, and Athlon 64 FX
# processors. (This supersets MMX, SSE, SSE2, 3DNow!, enhanced 3DNow!
# and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(athlon64|athlon64x2|turion64|turion64x2)-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=opteron
CFLAGSOPTF90  := -xtarget=opteron
CFLAGSOPTC    := -xtarget=opteron
CFLAGSOPTCXX  := -xtarget=opteron
LDFLAGSOPT    := -xtarget=opteron
endif

# Processors based on the AMD K8 core with x86-64 instruction set
# support, including the AMD Opteron, Athlon 64, and Athlon 64 FX
# processors. (This supersets MMX, SSE, SSE2, 3DNow!, enhanced 3DNow!
# and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(opteron|opteronx2)-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=opteron
CFLAGSOPTF90  := -xtarget=opteron
CFLAGSOPTC    := -xtarget=opteron
CFLAGSOPTCXX  := -xtarget=opteron
LDFLAGSOPT    := -xtarget=opteron
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-phenom-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=barcelona
CFLAGSOPTF90  := -xtarget=barcelona
CFLAGSOPTC    := -xtarget=barcelona
CFLAGSOPTCXX  := -xtarget=barcelona
LDFLAGSOPT    := -xtarget=barcelona
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-phenomII-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=barcelona
CFLAGSOPTF90  := -xtarget=barcelona
CFLAGSOPTC    := -xtarget=barcelona
CFLAGSOPTCXX  := -xtarget=barcelona
LDFLAGSOPT    := -xtarget=barcelona
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(opteron4100|opteron6100)-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=barcelona
CFLAGSOPTF90  := -xtarget=barcelona
CFLAGSOPTC    := -xtarget=barcelona
CFLAGSOPTCXX  := -xtarget=barcelona
LDFLAGSOPT    := -xtarget=barcelona
endif

# AMD Family 15: Bulldozer core
# CPUs based on AMD Family 15h cores with x86-64 instruction set
# support. (This supersets FMA4, AVX, XOP, LWP, AES, PCL_MUL, CX16,
# MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2, ABM and 64-bit
# instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-bulldozer-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=barcelona
CFLAGSOPTF90  := -xtarget=barcelona
CFLAGSOPTC    := -xtarget=barcelona
CFLAGSOPTCXX  := -xtarget=barcelona
LDFLAGSOPT    := -xtarget=barcelona
endif

# AMD Family 15h: Piledriver core
# AMD Family 15h core based CPUs with x86-64 instruction set
# support. (This supersets BMI, TBM, F16C, FMA, AVX, XOP, LWP, AES,
# PCL_MUL, CX16, MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2,
# ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-piledriver-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=barcelona
CFLAGSOPTF90  := -xtarget=barcelona
CFLAGSOPTC    := -xtarget=barcelona
CFLAGSOPTCXX  := -xtarget=barcelona
LDFLAGSOPT    := -xtarget=barcelona
endif

# AMD Family 15h: Steamroller core
# AMD Family 15h core based CPUs with x86-64 instruction set
# support. (This supersets BMI, TBM, F16C, FMA, AVX, XOP, LWP, AES,
# PCL_MUL, CX16, MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2,
# ABM and 64-bit instruction set extensions.
ifeq ($(call match,$(ID),(pc|pc64)-steamroller-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=barcelona
CFLAGSOPTF90  := -xtarget=barcelona
CFLAGSOPTC    := -xtarget=barcelona
CFLAGSOPTCXX  := -xtarget=barcelona
LDFLAGSOPT    := -xtarget=barcelona
endif

# AMD Family 14h: Bobcat core
# CPUs based on AMD Family 14h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSSE3, SSE4A, CX16,
# ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-bobcat-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=barcelona
CFLAGSOPTF90  := -xtarget=barcelona
CFLAGSOPTC    := -xtarget=barcelona
CFLAGSOPTCXX  := -xtarget=barcelona
LDFLAGSOPT    := -xtarget=barcelona
endif

# AMD Family 16h: Jaguar core
# CPUs based on AMD Family 16h cores with x86-64 instruction set
# support. This includes MOVBE, F16C, BMI, AVX, PCL_MUL, AES, SSE4.2,
# SSE4.1, CX16, ABM, SSE4A, SSSE3, SSE3, SSE2, SSE, MMX and 64-bit
# instruction set extensions.
ifeq ($(call match,$(ID),(pc|pc64)-jaguar-.*-sunstudio-.*),yes)
CFLAGSOPTF77  := -xtarget=barcelona
CFLAGSOPTF90  := -xtarget=barcelona
CFLAGSOPTC    := -xtarget=barcelona
CFLAGSOPTCXX  := -xtarget=barcelona
LDFLAGSOPT    := -xtarget=barcelona
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
