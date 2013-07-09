# -*- mode: makefile -*-

##############################################################################
# G95 (www.g95.org) with individual settings for x86/x86_64 CPU IDs
#
##############################################################################

# Option -march=cpu-type specifies that the code is compiled with
# optimization for the specified machine type "cpu-type". Since G95 is
# based on GCC 4.1.2, the handy option -march=native is not
# available. Moreover, optimization for newer CPUs is not possible.

##############################################################################
# Intel x86/x86_64 CPUs
#
##############################################################################

# List of supported cpu-types and minimally required compiler version
# CPU-TYPE     min. GCC version   min. G95 version
# i386         2.95.x             0.9
# i486         2.95.x             0.9
# i586/pentium 2.95.x             0.9
# pentiumpro   2.95.x             0.9
# pentium2     3.1.x              0.9
# pentium3     3.1.x              0.9
# pentiumm     3.4.1              0.9
# pentium4     3.1.x              0.9
# prescott     3.4.1              0.9
# nocona       3.4.1              0.9
# core2        4.3.x               -
# corei7       4.6.x               -
# corei7-avx   4.6.x               -
# core-avx-i   4.6.x               -
# core-avx2    4.7.x               -
# atom         4.5.x               -

# Original Intel i386 CPU.
ifeq ($(call match,$(ID),pc-i386-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=i386
CFLAGSOPTF90 := -march=i386
CFLAGSOPTC   := -march=i386
CFLAGSOPTCXX := -march=i386
LDFLAGSOPT   := -march=i386
endif

# Intel i486 CPU. (No scheduling is implemented for this chip.)
ifeq ($(call match,$(ID),pc-i486-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=i486
CFLAGSOPTF90 := -march=i486
CFLAGSOPTC   := -march=i486
CFLAGSOPTCXX := -march=i486
LDFLAGSOPT   := -march=i486
endif

# Intel Pentium CPU with or without MMX support.
ifeq ($(call match,$(ID),pc-(i568|pentium)-.*-g95-.*),yes)
ifeq ($(call match,$(core),"pentium-mmx"),yes)
CFLAGSOPTF77 := -march=pentium-mmx
CFLAGSOPTF90 := -march=pentium-mmx
CFLAGSOPTC   := -march=pentium-mmx
CFLAGSOPTCXX := -march=pentium-mmx
LDFLAGSOPT   := -march=pentium-mmx
else
CFLAGSOPTF77 := -march=pentium
CFLAGSOPTF90 := -march=pentium
CFLAGSOPTC   := -march=pentium
CFLAGSOPTCXX := -march=pentium
LDFLAGSOPT   := -march=pentium
endif
endif

# Intel Pentium Pro CPU.
ifeq ($(call match,$(ID),pc-pentiumpro-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=pentiumpro
CFLAGSOPTF90 := -march=pentiumpro
CFLAGSOPTC   := -march=pentiumpro
CFLAGSOPTCXX := -march=pentiumpro
LDFLAGSOPT   := -march=pentiumpro
endif

# Intel Pentium II CPU, based on Pentium Pro core with MMX instruction
# set support.
ifeq ($(call match,$(ID),pc-pentium2-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=pentium2
CFLAGSOPTF90 := -march=pentium2
CFLAGSOPTC   := -march=pentium2
CFLAGSOPTCXX := -march=pentium2
LDFLAGSOPT   := -march=pentium2
endif

# Intel Pentium III CPU, based on Pentium Pro core with MMX and SSE
# instruction set support.
ifeq ($(call match,$(ID),pc-pentium3-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=pentium3
CFLAGSOPTF90 := -march=pentium3
CFLAGSOPTC   := -march=pentium3
CFLAGSOPTCXX := -march=pentium3
LDFLAGSOPT   := -march=pentium3
endif

# Intel Pentium M; low-power version of Intel Pentium III CPU with MMX, 
# SSE and SSE2 instruction set support. Used by Centrino notebooks.
ifeq ($(call match,$(ID),pc-pentiumm-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=pentium-m
CFLAGSOPTF90 := -march=pentium-m
CFLAGSOPTC   := -march=pentium-m
CFLAGSOPTCXX := -march=pentium-m
LDFLAGSOPT   := -march=pentium-m
endif

# Intel Pentium 4 CPU with MMX, SSE and SSE2 instruction set support.
ifeq ($(call match,$(ID),pc-pentium4-.*-g95-.*),yes)
ifeq ($(call match,$(core),"prescott"),yes)
CFLAGSOPTF77 := -march=prescott
CFLAGSOPTF90 := -march=prescott
CFLAGSOPTC   := -march=prescott
CFLAGSOPTCXX := -march=prescott
LDFLAGSOPT   := -march=prescott
else
ifeq ($(call match,$(core),"nocona"),yes)
CFLAGSOPTF77 := -march=nocona
CFLAGSOPTF90 := -march=nocona
CFLAGSOPTC   := -march=nocona
CFLAGSOPTCXX := -march=nocona
LDFLAGSOPT   := -march=nocona
else
CFLAGSOPTF77 := -march=pentium4
CFLAGSOPTF90 := -march=pentium4
CFLAGSOPTC   := -march=pentium4
CFLAGSOPTCXX := -march=pentium4
LDFLAGSOPT   := -march=pentium4
endif
endif
endif

# Intel CoreSolo/Duo. Improved version of Intel Pentium 4 CPU with MMX,
# SSE, SSE2 and SSE3 instruction set support.
ifeq ($(call match,$(ID),pc-coresolo-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=prescott -mtune=prescott -msse
CFLAGSOPTF90 := -march=prescott -mtune=prescott -msse
CFLAGSOPTC   := -march=prescott -mtune=prescott -msse
CFLAGSOPTCXX := -march=prescott -mtune=prescott -msse
LDFLAGSOPT   := -march=prescott -mtune=prescott -msse
endif

# Intel Core 2 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support. Compiler versions below 4.3 do not
# support -march=core2 so that -march=prescott has to be used.
ifeq ($(call match,$(ID),pc-(coreduo|penryn)-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=prescott
CFLAGSOPTF90 := -march=prescott
LDFLAGSOPT   := -march=prescott
ifeq ($(call gccminversion,4,3),yes)
CFLAGSOPTC   := -march=core2
CFLAGSOPTCXX := -march=core2
else
CFLAGSOPTC   := -march=prescott
CFLAGSOPTCXX := -march=prescott
endif
endif

# Intel Core 2 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support. Compiler versions below 4.3 do not
# support -march=core2 so that -march=nocona has to be used.
ifeq ($(call match,$(ID),pc64-(coreduo|penryn)-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=nocona
CFLAGSOPTF90 := -march=nocona
LDFLAGSOPT   := -march=nocona
ifeq ($(call gccminversion,4,3),yes)
CFLAGSOPTC   := -march=core2
CFLAGSOPTCXX := -march=core2
else
CFLAGSOPTC   := -march=prescott
CFLAGSOPTCXX := -march=prescott
endif
endif

# Intel Nehalem CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support. Compiler versions below 4.3 do not
# support -march=core2 so that -march=prescott has to be used.
ifeq ($(call match,$(ID),pc-nehalem-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=prescott
CFLAGSOPTF90 := -march=prescott
LDFLAGSOPT   := -march=prescott
ifeq ($(call gccminversion,4,3),yes)
CFLAGSOPTC   := -march=core2
CFLAGSOPTCXX := -march=core2
else
CFLAGSOPTC   := -march=prescott
CFLAGSOPTCXX := -march=prescott
endif
endif

# Intel Nehalem CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support. Compiler versions below 4.3 do not
# support -march=core2 so that -march=nocona has to be used.
ifeq ($(call match,$(ID),pc64-nehalem-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=nocona
CFLAGSOPTF90 := -march=nocona
LDFLAGSOPT   := -march=nocona
ifeq ($(call gccminversion,4,3),yes)
CFLAGSOPTC   := -march=core2
CFLAGSOPTCXX := -march=core2
else
CFLAGSOPTC   := -march=nocona
CFLAGSOPTCXX := -march=nocona
endif
endif

# Intel Core i7 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3,
# SSSE3, SSE4.1, SSE4.2, AVX, AES and PCLMUL instruction set
# support. Compiler versions below 4.6 do not support
# -march=corei7-avx so that -march=prescott has to be used.
ifeq ($(call match,$(ID),pc-sandybridge-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=prescott
CFLAGSOPTF90 := -march=prescott
LDFLAGSOPT   := -march=prescott
ifeq ($(call gccminversion,4,6),yes)
CFLAGSOPTC   := -march=corei7-avx
CFLAGSOPTCXX := -march=corei7-avx
else
CFLAGSOPTC   := -march=prescott
CFLAGSOPTCXX := -march=prescott
endif
endif

# Intel Core i7 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3,
# SSSE3, SSE4.1, SSE4.2, AVX, AES and PCLMUL instruction set
# support. Compiler versions below 4.6 do not support
# -march=corei7-avx so that -march=prescott has to be used.
ifeq ($(call match,$(ID),pc64-sandybridge-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=nocona
CFLAGSOPTF90 := -march=nocona
LDFLAGSOPT   := -march=nocona
ifeq ($(call gccminversion,4,6),yes)
CFLAGSOPTC   := -march=corei7-avx
CFLAGSOPTCXX := -march=corei7-avx
else
CFLAGSOPTC   := -march=nocona
CFLAGSOPTCXX := -march=nocona
endif
endif

# Intel Core CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3, SSSE3,
# SSE4.1, SSE4.2, AVX, AES, PCLMUL, FSGSBASE, RDRND and F16C
# instruction set support. Compiler versions below 4.6 do not support
# -march=core-avx-i so that -march=prescott has to be used.
ifeq ($(call match,$(ID),pc-ivybridge-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=prescott
CFLAGSOPTF90 := -march=prescott
LDFLAGSOPT   := -march=prescott
ifeq ($(call gccminversion,4,6),yes)
CFLAGSOPTC   := -march=core-avx-i
CFLAGSOPTCXX := -march=core-avx-i
else
CFLAGSOPTC   := -march=prescott
CFLAGSOPTCX := -march=prescott
endif
endif

# Intel Core CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3, SSSE3,
# SSE4.1, SSE4.2, AVX, AES, PCLMUL, FSGSBASE, RDRND and F16C
# instruction set support. Compiler versions below 4.6 do not support
# -march=core-avx-i so that -march=nocona has to be used.
ifeq ($(call match,$(ID),pc64-ivybridge-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=nocona
CFLAGSOPTF90 := -march=nocona
LDFLAGSOPT   := -march=nocona
ifeq ($(call gccminversion,4,6),yes)
CFLAGSOPTC   := -march=core-avx-i
CFLAGSOPTCXX := -march=core-avx-i
else
CFLAGSOPTC   := -march=nocona
CFLAGSOPTCXX := -march=nocona
endif
endif

# Intel Core CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3, SSSE3,
# SSE4.1, SSE4.2, AVX, AES, PCLMUL, FSGSBASE, RDRND and F16C
# instruction set support. Compiler versions below 4.7 do not support
# -march=core-avx2 so that -march=prescott has to be used.
ifeq ($(call match,$(ID),pc-haswell-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=prescott
CFLAGSOPTF90 := -march=prescott
LDFLAGSOPT   := -march=prescott
ifeq ($(call gccminversion,4,7),yes)
CFLAGSOPTC   := -march=core-avx2
CFLAGSOPTCXX := -march=core-avx2
else
CFLAGSOPTC   := -march=prescott
CFLAGSOPTCXX := -march=prescott
endif
endif

# Intel Core CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3, SSSE3,
# SSE4.1, SSE4.2, AVX, AES, PCLMUL, FSGSBASE, RDRND and F16C
# instruction set support. Compiler versions below 4.7 do not support
# -march=core-avx2 so that -march=nocona has to be used.
ifeq ($(call match,$(ID),pc64-haswell-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=nocona
CFLAGSOPTF90 := -march=nocona
LDFLAGSOPT   := -march=nocona
ifeq ($(call gccminversion,4,7),yes)
CFLAGSOPTC   := -march=core-avx2
CFLAGSOPTCXX := -march=core-avx2
else
CFLAGSOPTC   := -march=nocona
CFLAGSOPTCXX := -march=nocona
endif
endif

# Intel Atom CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support. Compiler versions below 4.3 do not
# support -march=atom so that -march=prescott has to be used.
ifeq ($(call match,$(ID),pc-atom-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=prescott
CFLAGSOPTF90 := -march=prescott
LDFLAGSOPT   := -march=prescott
ifeq ($(call gccminversion,4,5),yes)
CFLAGSOPTC   := -march=atom
CFLAGSOPTCXX := -march=atom
else
CFLAGSOPTC   := -march=prescott
CFLAGSOPTCXX := -march=prescott
endif
endif

# Intel Atom CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support. Compiler versions below 4.3 do not
# support -march=atom so that -march=nocona has to be used.
ifeq ($(call match,$(ID),pc64-atom-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=nocona
CFLAGSOPTF90 := -march=nocona
LDFLAGSOPT   := -march=nocona
ifeq ($(call gccminversion,4,5),yes)
CFLAGSOPTC   := -march=atom
CFLAGSOPTCXX := -march=atom
else
CFLAGSOPTC   := -march=nocona
CFLAGSOPTCXX := -march=nocona
endif
endif

##############################################################################
# AMD x86/x86_64 CPUs
#
##############################################################################

# List of supported cpu-types and minimally required compiler version
# CPU-TYPE      min. GCC version   min. G95 version
# k6            2.95.x             0.9
# k6-2          3.1.x              0.9
# k6-3          3.1.x              0.9
# athlon        3.0.x              0.9
# athlon-tbird  3.1.x              0.9
# athlon-4      3.1.x              0.9
# athlon-xp     3.1.x              0.9
# athlon-mp     3.1.x              0.9
# k8            3.4.x              0.9
# opteron       3.4.x              0.9
# athlon64      3.4.0              0.9
# athlon-fx     3.4.x              0.9
# k8-sse3       4.3.x               -
# opteron-sse3  4.3.x               -
# athlon64-sse3 4.3.x               -
# amdfam10      4.3.x               -
# barcelona     4.3.x               -
# bdver1        4.6.x               -
# bdver2        4.7.x               -
# bdver3        4.8.x               -
# btver1        4.6.x               -
# btver2        4.8.x               -

# AMD K6 CPU with MMX instruction set support and improved versions of
# AMD K6 CPU with MMX and 3DNow! instruction set support,
# respectively.
ifeq ($(call match,$(ID),pc-k6-.*-g95-.*),yes)
ifeq ($(call match,$(core),"k6-2"),yes)
CFLAGSOPTF77 := -march=k6-2
CFLAGSOPTF90 := -march=k6-2
CFLAGSOPTC   := -march=k6-2
CFLAGSOPTCXX := -march=k6-2
LDFLAGSOPT   := -march=k6-2
else
ifeq ($(call match,$(core),"k6-3"),yes)
CFLAGSOPTF77 := -march=k6-3
CFLAGSOPTF90 := -march=k6-3
CFLAGSOPTC   := -march=k6-3
CFLAGSOPTCXX := -march=k6-3
LDFLAGSOPT   := -march=k6-3
else
CFLAGSOPTF77 := -march=k6
CFLAGSOPTF90 := -march=k6
CFLAGSOPTC   := -march=k6
CFLAGSOPTCXX := -march=k6
LDFLAGSOPT   := -march=k6
endif
endif
endif

# AMD Athlon CPU with MMX, 3dNOW!, enhanced 3DNow! and SSE prefetch
# instructions support.
ifeq ($(call match,$(ID),pc-athlon-.*-g95-.*),yes)
ifeq ($(call match,$(core),"thunderbird"),yes)
CFLAGSOPTF77 := -march=athlon-tbird
CFLAGSOPTF90 := -march=athlon-tbird
CFLAGSOPTC   := -march=athlon-tbird
CFLAGSOPTCXX := -march=athlon-tbird
LDFLAGSOPT   := -march=athlon-tbird
else
CFLAGSOPTF77 := -march=athlon
CFLAGSOPTF90 := -march=athlon
CFLAGSOPTC   := -march=athlon
CFLAGSOPTCXX := -march=athlon
LDFLAGSOPT   := -march=athlon
endif
endif

# Improved AMD Athlon CPU with MMX, 3DNow!, enhanced 3DNow! and full
# SSE instruction set support.
ifeq ($(call match,$(ID),pc-athlonxp-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=athlon-xp
CFLAGSOPTF90 := -march=athlon-xp
CFLAGSOPTC   := -march=athlon-xp
CFLAGSOPTCXX := -march=athlon-xp
LDFLAGSOPT   := -march=athlon-xp
endif

# Processors based on the AMD K8 core with x86-64 instruction set
# support, including the AMD Opteron, Athlon 64, and Athlon 64 FX
# processors. (This supersets MMX, SSE, SSE2, 3DNow!, enhanced 3DNow!
# and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(athlon64|athlon64x2|turion64|turion64x2)-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=athlon64
CFLAGSOPTF90 := -march=athlon64
CFLAGSOPTC   := -march=athlon64
CFLAGSOPTCXX := -march=athlon64
LDFLAGSOPT   := -march=athlon64
endif

# Processors based on the AMD K8 core with x86-64 instruction set
# support, including the AMD Opteron, Athlon 64, and Athlon 64 FX
# processors. (This supersets MMX, SSE, SSE2, 3DNow!, enhanced 3DNow!
# and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(opteron|opteronx2)-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=opteron
CFLAGSOPTF90 := -march=opteron
CFLAGSOPTC   := -march=opteron
CFLAGSOPTCXX := -march=opteron
LDFLAGSOPT   := -march=opteron
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-phenom-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=opteron
CFLAGSOPTF90 := -march=opteron
LDFLAGSOPT   := -march=opteron
ifeq ($(call gccminversion,4,3),yes)
CFLAGSOPTC   := -march=amdfam10
CFLAGSOPTCXX := -march=amdfam10
else
CFLAGSOPTC   := -march=opteron
CFLAGSOPTCXX := -march=opteron
endif
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-phenomII-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=opteron
CFLAGSOPTF90 := -march=opteron
LDFLAGSOPT   := -march=opteron
ifeq ($(call gccminversion,4,3),yes)
CFLAGSOPTC   := -march=amdfam10
CFLAGSOPTCXX := -march=amdfam10
else
CFLAGSOPTC   := -march=opteron
CFLAGSOPTCXX := -march=opteron
endif
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(opteron4100|opteron6100)-.*-g95-.*),yes)
CFLAGSOPTF77 := -march=opteron
CFLAGSOPTF90 := -march=opteron
LDFLAGSOPT   := -march=opteron
ifeq ($(call gccminversion,4,3),yes)
CFLAGSOPTC   := -march=barcelona
CFLAGSOPTCXX := -march=barcelona
else
CFLAGSOPTC   := -march=opteron
CFLAGSOPTCXX := -march=opteron
endif
endif

# AMD Family 15: Bulldozer core
# CPUs based on AMD Family 15h cores with x86-64 instruction set
# support. (This supersets FMA4, AVX, XOP, LWP, AES, PCL_MUL, CX16,
# MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2, ABM and 64-bit
# instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-bulldozer-.*-gcc-.*),yes)
CFLAGSOPTF77 := -march=opteron
CFLAGSOPTF90 := -march=opteron
LDFLAGSOPT   := -march=opteron
ifeq ($(call gccminversion,4,6),yes)
CFLAGSOPTC   := -march=bdver1
CFLAGSOPTCXX := -march=bdver1
else
ifeq ($(call gccminversion,4,3),yes)
CFLAGSOPTC   := -march=amdfam10
CFLAGSOPTCXX := -march=amdfam10
else
CFLAGSOPTC   := -march=opteron
CFLAGSOPTCXX := -march=opteron
endif
endif
endif

# AMD Family 15h: Piledriver core
# AMD Family 15h core based CPUs with x86-64 instruction set
# support. (This supersets BMI, TBM, F16C, FMA, AVX, XOP, LWP, AES,
# PCL_MUL, CX16, MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2,
# ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-piledriver-.*-gcc-.*),yes)
CFLAGSOPTF77 := -march=opteron
CFLAGSOPTF90 := -march=opteron
LDFLAGSOPT   := -march=opteron
ifeq ($(call gccminversion,4,7),yes)
CFLAGSOPTC   := -march=bdver2
CFLAGSOPTCXX := -march=bdver2
else
ifeq ($(call gccminversion,4,3),yes)
CFLAGSOPTC   := -march=amdfam10
CFLAGSOPTCXX := -march=amdfam10
else
CFLAGSOPTC   := -march=opteron
CFLAGSOPTCXX := -march=opteron
endif
endif
endif

# AMD Family 15h: Steamroller core
# AMD Family 15h core based CPUs with x86-64 instruction set
# support. (This supersets BMI, TBM, F16C, FMA, AVX, XOP, LWP, AES,
# PCL_MUL, CX16, MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2,
# ABM and 64-bit instruction set extensions.
ifeq ($(call match,$(ID),(pc|pc64)-steamroller-.*-gcc-.*),yes)
CFLAGSOPTF77 := -march=opteron
CFLAGSOPTF90 := -march=opteron
LDFLAGSOPT   := -march=opteron
ifeq ($(call gccminversion,4,8),yes)
CFLAGSOPTC   := -march=bdver3
CFLAGSOPTCXX := -march=bdver3
else
ifeq ($(call gccminversion,4,3),yes)
CFLAGSOPTC   := -march=amdfam10
CFLAGSOPTCXX := -march=amdfam10
else
CFLAGSOPTC   := -march=opteron
CFLAGSOPTCXX := -march=opteron
endif
endif
endif

# AMD Family 14h: Bobcat core
# CPUs based on AMD Family 14h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSSE3, SSE4A, CX16,
# ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-bobcat-.*-gcc-.*),yes)
CFLAGSOPTF77 := -march=opteron
CFLAGSOPTF90 := -march=opteron
LDFLAGSOPT   := -march=opteron
ifeq ($(call gccminversion,4,6),yes)
CFLAGSOPTC   := -march=btver1
CFLAGSOPTCXX := -march=btver1
else
ifeq ($(call gccminversion,4,3),yes)
CFLAGSOPTC   := -march=amdfam10
CFLAGSOPTCXX := -march=amdfam10
else
CFLAGSOPTC   := -march=opteron
CFLAGSOPTCXX := -march=opteron
endif
endif
endif


# AMD Family 16h: Jaguar core
# CPUs based on AMD Family 16h cores with x86-64 instruction set
# support. This includes MOVBE, F16C, BMI, AVX, PCL_MUL, AES, SSE4.2,
# SSE4.1, CX16, ABM, SSE4A, SSSE3, SSE3, SSE2, SSE, MMX and 64-bit
# instruction set extensions.
ifeq ($(call match,$(ID),(pc|pc64)-jaguar-.*-gcc-.*),yes)
CFLAGSOPTF77 := -march=opteron
CFLAGSOPTF90 := -march=opteron
LDFLAGSOPT   := -march=opteron
ifeq ($(call gccminversion,4,8),yes)
CFLAGSOPTC   := -march=btver2
CFLAGSOPTCXX := -march=btver2
else
ifeq ($(call gccminversion,4,3),yes)
CFLAGSOPTC   := -march=amdfam10
CFLAGSOPTCXX := -march=amdfam10
else
CFLAGSOPTC   := -march=opteron
CFLAGSOPTCXX := -march=opteron
endif
endif
endif

##############################################################################
# Automatic CPU type detection
#
##############################################################################

ifneq (,$(findstring NATIVE,$(OPT)))
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: No auto-tuning for g95 compiler available!';
endif
