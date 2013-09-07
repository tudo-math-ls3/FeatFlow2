# -*- mode: makefile -*-

##############################################################################
# Intel Compiler suite with individual settings for x86/x86_64 CPU IDs
#
##############################################################################

# The intel compiler suite frequently changed the naming convention
# for optimization flags so that it is necessary to carefully check
# for the compiler versions.
#
# Intel Compiler suite prior to and including 8.1
#
# -tpp5  Intel Pentium CPU
# -tpp6  Intel Pentium Pro, II, II CPU
# -tpp7  Intel Pentium 4 CPU
#
#
# Intel Compiler suite prior to and including 9.1
#
# -mtune=<cpu>  optimize for a specific cpu
#        pentium    - optimize for Pentium(R) processor
#        pentiumpro - optimize for Pentium(R) Pro, Pentium(R) II and Pentium(R)
#                     III processors
#        pentium4   - optimize for Pentium(R) 4 processor (DEFAULT)
# -march=<cpu>  generate code excusively for a given <cpu>
#        pentiumpro - Pentium(R) Pro and Pentium(R) II processor instructions
#        pentiumii  - MMX(TM)instructions
#        pentiumiii - streaming SIMD extensions
#        pentium4   - Pentium(R) 4 New Instructions
# -ax<codes> generate code specialized for processors specified by <codes>
#            while also generating generic IA-32 code.  <codes> includes
#            one or more of the following characters:
#     W  Intel Pentium 4 and compatible Intel processors
#     P  Intel Core(TM) Duo processors, Intel Core(TM) Solo processors, Intel
#        Pentium 4 and compatible Intel processors with Streaming SIMD
#        Extensions 3 (SSE3) instruction support
# -x<codes>  generate specialized code to run exclusively on processors
#            indicated by <codes> as described above.
# -tune <keyword> 
#        pn1 - optimize for Pentium(R) processor
#        pn2 - optimize for Pentium(R) Pro, Pentium(R) II, and 
#              Pentium(R) III processors
#        pn3 - same as pn2
#        pn4 - optimize for Pentium(R) 4 processor (DEFAULT)
# -arch <keyword>
#        pn1 - optimize for Pentium(R) processor
#        pn2 - optimize for Pentium(R) Pro, Pentium(R) II, and 
#              Pentium(R) III processors
#        pn3 - same as pn2
#        pn4 - optimize for Pentium(R) 4 processor (DEFAULT)
# -msse3 generate code for Intel Core(TM) Duo processors, Intel Core(TM) Solo
#        processors, Intel Pentium 4 and compatible Intel processors with
#        Streaming SIMD Extensions 3 (SSE3) instruction support
#
#
# Intel Compiler suite prior to and including 10.1
#
# -mtune=<cpu>  optimize for a specific cpu
#        pentium3   - optimize for Intel Pentium(R) III processors
#        pentium4   - optimize for Intel Pentium(R) 4 processor (DEFAULT)
#        core2      - optimize for Intel Core(TM)2 processor family
# -march=<cpu>  generate code excusively for a given <cpu>
#        pentium3   - streaming SIMD extensions
#        pentium4   - Intel Pentium(R) 4 New Instructions
#        core2      - Intel Core(TM)2 processor family
# -ax<codes> generate code specialized for processors specified by <codes>
#            while also generating generic IA-32 instructions.  <codes>
#            includes one or more of the following characters:
#     W  Intel Pentium 4 and compatible Intel processors
#     P  Intel(R) Core(TM) processor family with Streaming SIMD
#        Extensions 3 (SSE3) instruction support
#     T  Intel(R) Core(TM)2 processor family with SSSE3
#     S  Future Intel processors supporting SSE4 Vectorizing Compiler and 
#        Media Accelerator instructions
# -x<codes>  generate specialized code to run exclusively on processors
#            indicated by <codes> as described below
#     W  Intel Pentium 4 and compatible Intel processors
#     P  Intel(R) Core(TM) processor family with Streaming SIMD
#        Extensions 3 (SSE3) instruction support
#     T  Intel(R) Core(TM)2 processor family with SSSE3
#     O  Intel(R) Core(TM) processor family.  Code is expected to run properly
#        on any processor that supports SSE3, SSE2 and SSE instruction sets
#     S  Future Intel processors supporting SSE4 Vectorizing Compiler and 
#        Media Accelerator instructions
# -tune <keyword> 
#        pn1 - optimize for Pentium(R) processor
#        pn2 - optimize for Pentium(R) Pro, Pentium(R) II, and 
#              Pentium(R) III processors
#        pn3 - same as pn2
#        pn4 - optimize for Pentium(R) 4 processor (DEFAULT)
# -arch <keyword>
#        pn1 - optimize for Pentium(R) processor
#        pn2 - optimize for Pentium(R) Pro, Pentium(R) II, and 
#              Pentium(R) III processors
#        pn3 - same as pn2
#        pn4 - optimize for Pentium(R) 4 processor (DEFAULT)
# -msse3 generate code for Intel(R) Core(TM) Duo processors, Intel(R) Core(TM)
#        Solo processors, Intel Pentium 4 and compatible Intel processors with
#        Streaming SIMD Extensions 3 (SSE3) instruction support
#
#
# Intel Compiler suite 11.1, 12.0, 12.1 and 13.0
#
# -x<code>  generate specialized code to run exclusively on processors
#           indicated by <code> as described below
#             SSE2    May generate Intel(R) SSE2 and SSE instructions for Intel
#                     processors.  Optimizes for the Intel NetBurst(R)
#                     microarchitecture.
#             SSE3    May generate Intel(R) SSE3, SSE2, and SSE instructions for
#                     Intel processors.  Optimizes for the enhanced Pentium(R) M 
#                     processor microarchitecture and Intel NetBurst(R)
#                     microarchitecture. 
#             SSSE3   May generate Intel(R) SSSE3, SSE3, SSE2, and SSE
#                     instructions for Intel processors.  Optimizes for the
#                     Intel(R) Core(TM) microarchitecture.
#             SSE4.1  May generate Intel(R) SSE4 Vectorizing Compiler and Media
#                     Accelerator instructions for Intel processors.  May 
#                     generate Intel(R) SSSE3, SSE3, SSE2, and SSE instructions
#                     and it may optimize for Intel(R) 45nm Hi-k next generation
#                     Intel Core(TM) microarchitecture.
#             SSE4.2  May generate Intel(R) SSE4 Efficient Accelerated String
#                     and Text Processing instructions supported by Intel(R)
#                     Core(TM) i7 processors.  May generate Intel(R) SSE4 
#                     Vectorizing Compiler and Media Accelerator, Intel(R) SSSE3,
#                     SSE3, SSE2, and SSE instructions and it may optimize for
#                     the Intel(R) Core(TM) processor family.
#             AVX     May generate Intel(R) Advanced Vector Extensions (Intel(R)
#                     AVX), Intel(R) SSE4.2, SSE4.1, SSSE3, SSE3,
#                     SSE2, and SSE instructions for Intel(R) processors.
#                     Optimizes for a future Intel processor.
#             CORE-AVX2
#                     May generate Intel(R) Advanced Vector Extensions 2
#                     (Intel(R) AVX2), Intel(R) AVX, SSE4.2, SSE4.1, SSSE3, SSE3,
#                     SSE2, and SSE instructions for Intel(R) processors.
#                     Optimizes for a future Intel processor.
#             CORE-AVX-I
#                     May generate Intel(R) Advanced Vector Extensions (Intel(R)
#                     AVX), including instructions in Intel(R) Core 2(TM)
#                     processors in process technology smaller than 32nm,
#                     Intel(R) SSE4.2, SSE4.1, SSSE3, SSE3, SSE2, and SSE
#                     instructions for Intel(R) processors. Optimizes for a
#                     future Intel processor. 
#             SSSE3_ATOM
#                     May generate MOVBE instructions for Intel processors,
#                     depending on the setting of option -minstruction.
#                     May also generate Intel(R) SSSE3, SSE3, SSE2, and SSE
#                     instructions for Intel processors. Optimizes for the
#                     Intel(R) Atom(TM) processor and Intel(R) Centrino(R)
#                     Atom(TM) Processor Technology.
# -xHost    generate instructions for the highest instruction set and processor
#           available on the compilation host machine
# -arch <code>
#           generate specialized code to optimize for processors indicated by
#           <code> as described below
#             SSE2    May generate Intel(R) SSE2 and SSE instructions
#             SSE3    May generate Intel(R) SSE3, SSE2 and SSE instructions
#             SSSE3   May generate Intel(R) SSSE3, SSE3, SSE2 and SSE
#                     instructions
#             SSE4.1  May generate Intel(R) SSE4.1, SSSE3, SSE3, SSE2 and SSE
#                     instructions
#             SSE4.2  May generate Intel(R) SSE4.2, SSE4.1, SSSE3, SSE3, SSE2 and
#                     SSE instructions
#             AVX     May generate Intel(R) AVX, SSE4.2, SSE4.1, SSSE3, SSE3,
#                     SSE2 and SSE instructions
# -mcpu=<cpu>
#           same as -mtune=<cpu>
# -mtune=<cpu>
#           optimize for a specific <cpu>
#             pentium3  - optimize for Pentium(R) III processors
#             pentium4  - optimize for Pentium(R) 4 processor (DEFAULT)
# -march=<cpu>
#           generate code exclusively for a given <cpu>
#             core-avx2  - processors that support Intel(R) Advanced Vector
#                          Extensions 2 (Intel(R) AVX2)
#             core-avx-i - processors that support Intel(R) Advanced Vector
#                          Extensions (Intel(R) AVX), including instructions in
#                          Intel(R) Core 2(TM) processors in process technology
#                          smaller than 32nm
#             corei7-avx - processors that support Intel(R) Advanced Vector
#                          Extensions (Intel(R) AVX)
#             corei7     - processors that support Intel(R) SSE4 Efficient
#                          Accelerated String and Text Processing instructions
#             atom       - processors that support MOVBE instructions
#             core2      - Intel(R) Core 2(TM) processor family
#             pentium-m  - Intel(R) Pentium(R) M processors
#             pentium4   - Intel(R) Pentium(R) 4 processors
#             pentium3   - Intel(R) Pentium(R) III processors (Linux only)
# -msse3    May generate Intel(R) SSE3, SSE2, and SSE instructions
# -mssse3   May generate Intel(R) SSSE3, SSE3, SSE2, and SSE instructions
# -msse4    Enable -msse4.2
# -msse4.1  May generate Intel(R) SSE4.1, SSSE3, SSE3, SSE2, and SSE instructions
# -msse4.2  May generate Intel(R) SSE4.2, SSE4.1, SSSE3, SSE3, SSE2, and SSE
#           instructions
# -mavx     May generate Intel(R) AVX, SSE4.2, SSE4.1, SSSE3, SSE3, SSE2, and SSE
#           instructions

##############################################################################
# Intel x86/x86_64 CPUs
#
##############################################################################

# Original Intel i386 CPU.
ifeq ($(call match,$(ID),pc-i386-.*-intel-.*),yes)
ifeq ($(call intelminversion,10,1),yes)
CFLAGSOPTF77 := -mia32
CFLAGSOPTF90 := -mia32
CFLAGSOPTC   := -mia32
CFLAGSOPTCXX := -mia32
LDFLAGSOPT   := -mia32
else
CFLAGSOPTF77 := 
CFLAGSOPTF90 := 
CFLAGSOPTC   := 
CFLAGSOPTCXX := 
LDFLAGSOPT   := 
endif
endif

# Intel i486 CPU. (No scheduling is implemented for this chip.)
ifeq ($(call match,$(ID),pc-i486-.*-intel-.*),yes)
ifeq ($(call intelminversion,10,1),yes)
CFLAGSOPTF77 := -mia32
CFLAGSOPTF90 := -mia32
CFLAGSOPTC   := -mia32
CFLAGSOPTCXX := -mia32
LDFLAGSOPT   := -mia32
else
CFLAGSOPTF77 := 
CFLAGSOPTF90 := 
CFLAGSOPTC   := 
CFLAGSOPTCXX := 
LDFLAGSOPT   := 
endif
endif

# Intel Pentium CPU with or without MMX support.
ifeq ($(call match,$(ID),pc-(i568|pentium)-.*-intel-.*),yes)
ifeq ($(call intelminversion,10,1),yes)
CFLAGSOPTF77 := -mia32
CFLAGSOPTF90 := -mia32
CFLAGSOPTC   := -mia32
CFLAGSOPTCXX := -mia32
LDFLAGSOPT   := -mia32
else
CFLAGSOPTF77 := -tpp5
CFLAGSOPTF90 := -tpp5
CFLAGSOPTC   := -tpp5
CFLAGSOPTCXX := -tpp5
LDFLAGSOPT   := -tpp5
endif
endif

# Intel Pentium Pro CPU.
ifeq ($(call match,$(ID),pc-pentiumpro-.*-intel-.*),yes)
ifeq ($(call intelminversion,10,1),yes)
CFLAGSOPTF77 := -mia32
CFLAGSOPTF90 := -mia32
CFLAGSOPTC   := -mia32
CFLAGSOPTCXX := -mia32
LDFLAGSOPT   := -mia32
else
CFLAGSOPTF77 := -tpp6
CFLAGSOPTF90 := -tpp6
CFLAGSOPTC   := -tpp6
CFLAGSOPTCXX := -tpp6
LDFLAGSOPT   := -tpp6
endif
endif

# Intel Pentium II CPU, based on Pentium Pro core with MMX instruction
# set support.
ifeq ($(call match,$(ID),pc-pentium2-.*-intel-.*),yes)
ifeq ($(call intelminversion,10,1),yes)
CFLAGSOPTF77 := -mia32
CFLAGSOPTF90 := -mia32
CFLAGSOPTC   := -mia32
CFLAGSOPTCXX := -mia32
LDFLAGSOPT   := -mia32
else
CFLAGSOPTF77 := -tpp6
CFLAGSOPTF90 := -tpp6
CFLAGSOPTC   := -tpp6
CFLAGSOPTCXX := -tpp6
LDFLAGSOPT   := -tpp6
endif
endif

# Intel Pentium III CPU, based on Pentium Pro core with MMX and SSE
# instruction set support.
ifeq ($(call match,$(ID),pc-pentium3-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -xSSE
CFLAGSOPTF90  := -xSSE
CFLAGSOPTC    := -xSSE
CFLAGSOPTCXX  := -xSSE
LDFLAGSOPT    := -xSSE
else
CFLAGSOPTF77  := -xK
CFLAGSOPTF90  := -xK
CFLAGSOPTC    := -xK
CFLAGSOPTCXX  := -xK
LDFLAGSOPT    := -xK
endif
endif

# Intel Pentium M; low-power version of Intel Pentium III CPU with MMX, 
# SSE and SSE2 instruction set support. Used by Centrino notebooks.
ifeq ($(call match,$(ID),pc-pentiumm-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -xSSE2
CFLAGSOPTF90  := -xSSE2
CFLAGSOPTC    := -xSSE2
CFLAGSOPTCXX  := -xSSE2
LDFLAGSOPT    := -xSSE2
else
CFLAGSOPTF77  := -xB
CFLAGSOPTF90  := -xB
CFLAGSOPTC    := -xB
CFLAGSOPTCXX  := -xB
LDFLAGSOPT    := -xB
endif
endif

# Intel Pentium 4 CPU with MMX, SSE and SSE2 instruction set support.
ifeq ($(call match,$(ID),pc-pentium4-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -xSSE2
CFLAGSOPTF90  := -xSSE2
CFLAGSOPTC    := -xSSE2
CFLAGSOPTCXX  := -xSSE2
LDFLAGSOPT    := -xSSE2
else
CFLAGSOPTF77  := -xN
CFLAGSOPTF90  := -xN
CFLAGSOPTC    := -xN
CFLAGSOPTCXX  := -xN
LDFLAGSOPT    := -xN
endif
endif

# Intel Pentium 4 CPU with MMX, SSE and SSE2 instruction set support.
ifeq ($(call match,$(ID),pc64-pentium4-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -xSSE2
CFLAGSOPTF90  := -xSSE2
CFLAGSOPTC    := -xSSE2
CFLAGSOPTCXX  := -xSSE2
LDFLAGSOPT    := -xSSE2
else
CFLAGSOPTF77  := -xW
CFLAGSOPTF90  := -xW
CFLAGSOPTC    := -xW
CFLAGSOPTCXX  := -xW
LDFLAGSOPT    := -xW
endif
endif

# Intel CoreSolo/Duo. Improved version of Intel Pentium 4 CPU with MMX,
# SSE, SSE2 and SSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-coresolo-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -msse3
CFLAGSOPTF90  := -msse3
CFLAGSOPTC    := -msse3
CFLAGSOPTCXX  := -msse3
LDFLAGSOPT    := -msse3
else
CFLAGSOPTF77  := -xO
CFLAGSOPTF90  := -xO
CFLAGSOPTC    := -xO
CFLAGSOPTCXX  := -xO
LDFLAGSOPT    := -xO
endif
endif

# Intel Core 2 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),pc-(coreduo|penryn)-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -xSSE3
CFLAGSOPTF90  := -xSSE3
CFLAGSOPTC    := -xSSE3
CFLAGSOPTCXX  := -xSSE3
LDFLAGSOPT    := -xSSE3
else
CFLAGSOPTF77  := -xP
CFLAGSOPTF90  := -xP
CFLAGSOPTC    := -xP
CFLAGSOPTCXX  := -xP
LDFLAGSOPT    := -xP
endif
endif

# Intel Core 2 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),pc64-(coreduo|penryn)-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -xSSE3
CFLAGSOPTF90  := -xSSE3
CFLAGSOPTC    := -xSSE3
CFLAGSOPTCXX  := -xSSE3
LDFLAGSOPT    := -xSSE3
else
CFLAGSOPTF77  := -xT
CFLAGSOPTF90  := -xT
CFLAGSOPTC    := -xT
CFLAGSOPTCXX  := -xT
LDFLAGSOPT    := -xT
endif
endif

# There settings are reported to work
ifeq ($(call match,$(ID),(pc|pc64)-penryn-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -xSSE4.1
CFLAGSOPTF90  := -xSSE4.1
CFLAGSOPTC    := -xSSE4.1
CFLAGSOPTCXX  := -xSSE4.1
LDFLAGSOPT    := -xSSE4.1
else
CFLAGSOPTF77  := -xT
CFLAGSOPTF90  := -xT
CFLAGSOPTC    := -xT
CFLAGSOPTCXX  := -xT
LDFLAGSOPT    := -xT
endif
endif

# Intel Nehalem CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-nehalem-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -xSSE4.2
CFLAGSOPTF90  := -xSSE4.2
CFLAGSOPTC    := -xSSE4.2
CFLAGSOPTCXX  := -xSSE4.2
LDFLAGSOPT    := -xSSE4.2
else
CFLAGSOPTF77  := -xT
CFLAGSOPTF90  := -xT
CFLAGSOPTC    := -xT
CFLAGSOPTCXX  := -xT
LDFLAGSOPT    := -xT
endif
endif

# Intel Core i7 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3,
# SSSE3, SSE4.1, SSE4.2, AVX, AES and PCLMUL instruction set
# support.
ifeq ($(call match,$(ID),(pc|pc64)-sandybridge-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77 := -xAVX
CFLAGSOPTF90 := -xAVX
CFLAGSOPTC   := -xAVX
CFLAGSOPTCXX := -xAVX
LDFLAGSOPT   := -xAVX
else
CFLAGSOPTF77 := -xS
CFLAGSOPTF90 := -xS
CFLAGSOPTC   := -xS
CFLAGSOPTCXX := -xS
LDFLAGSOPT   := -xS
endif
endif

ifeq ($(call match,$(ID),(pc|pc64)-sandybridge_ep-.*-intel-.*),yes)
ifeq ($(shell grep "\bavx\b" /proc/cpuinfo),)
# No AVX support in kernel, only available from 2.6.30 onwards
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: CPU has AVX features, but the (current) kernel not yet. Disabling AVX optimisations.'; \
	    echo '***          In case of cross-compilation, this auto-disabling might be wrong if the target kernel'; \
	    echo '***          is newer than the kernel of the current machine, namely >= 2.6.30'; \
CFLAGSOPTF77 := -D__pentium4 -D__pentium4__ -D__tune_pentium4__ -D__SSE2__ -D__SSE3__ -D__SSSE3__ -D__SSE4_1__ -D__SSE4_2__ -D__SSE__ -D__MMX__
CFLAGSOPTF90 := -D__pentium4 -D__pentium4__ -D__tune_pentium4__ -D__SSE2__ -D__SSE3__ -D__SSSE3__ -D__SSE4_1__ -D__SSE4_2__ -D__SSE__ -D__MMX__
CFLAGSOPTC   := -D__pentium4 -D__pentium4__ -D__tune_pentium4__ -D__SSE2__ -D__SSE3__ -D__SSSE3__ -D__SSE4_1__ -D__SSE4_2__ -D__SSE__ -D__MMX__
CFLAGSOPTCXX := -D__pentium4 -D__pentium4__ -D__tune_pentium4__ -D__SSE2__ -D__SSE3__ -D__SSSE3__ -D__SSE4_1__ -D__SSE4_2__ -D__SSE__ -D__MMX__
LDFLAGSOPT   := -D__pentium4 -D__pentium4__ -D__tune_pentium4__ -D__SSE2__ -D__SSE3__ -D__SSSE3__ -D__SSE4_1__ -D__SSE4_2__ -D__SSE__ -D__MMX__
else
ifeq ($(call intelminversion,11,1),yes)
# AVX support in kernel and compiler version
CFLAGSOPTF77 := -xAVX
CFLAGSOPTF90 := -xAVX
CFLAGSOPTC   := -xAVX
CFLAGSOPTCXX := -xAVX
LDFLAGSOPT   := -xAVX
else
CFLAGSOPTF77 := -xS
CFLAGSOPTF90 := -xS
CFLAGSOPTC   := -xS
CFLAGSOPTCXX := -xS
LDFLAGSOPT   := -xS
endif
endif
endif

# Intel Core CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3, SSSE3,
# SSE4.1, SSE4.2, AVX, AES, PCLMUL, FSGSBASE, RDRND and F16C
# instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-ivybridge-.*-intel-.*),yes)
ifeq ($(call intelminversion,12,1),yes)
CFLAGSOPTF77 := -xCORE-AVX-I
CFLAGSOPTF90 := -xCORE-AVX-I
CFLAGSOPTC   := -xCORE-AVX-I
CFLAGSOPTCXX := -xCORE-AVX-I
LDFLAGSOPT   := -xCORE-AVX-I
else
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77 := -xAVX
CFLAGSOPTF90 := -xAVX
CFLAGSOPTC   := -xAVX
CFLAGSOPTCXX := -xAVX
LDFLAGSOPT   := -xAVX
else
CFLAGSOPTF77 := -xS
CFLAGSOPTF90 := -xS
CFLAGSOPTC   := -xS
CFLAGSOPTCXX := -xS
LDFLAGSOPT   := -xS
endif
endif
endif

# Intel Core CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3, SSSE3,
# SSE4.1, SSE4.2, AVX, AES, PCLMUL, FSGSBASE, RDRND and F16C
# instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-haswell-.*-intel-.*),yes)
ifeq ($(call intelminversion,13,0),yes)
CFLAGSOPTF77 := -xCORE-AVX2
CFLAGSOPTF90 := -xCORE-AVX2
CFLAGSOPTC   := -xCORE-AVX2
CFLAGSOPTCXX := -xCORE-AVX2
LDFLAGSOPT   := -xCORE-AVX2
else
ifeq ($(call intelminversion,12,1),yes)
CFLAGSOPTF77 := -xCORE-AVX-I
CFLAGSOPTF90 := -xCORE-AVX-I
CFLAGSOPTC   := -xCORE-AVX-I
CFLAGSOPTCXX := -xCORE-AVX-I
LDFLAGSOPT   := -xCORE-AVX-I
else
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77 := -xAVX
CFLAGSOPTF90 := -xAVX
CFLAGSOPTC   := -xAVX
CFLAGSOPTCXX := -xAVX
LDFLAGSOPT   := -xAVX
else
CFLAGSOPTF77 := -xS
CFLAGSOPTF90 := -xS
CFLAGSOPTC   := -xS
CFLAGSOPTCXX := -xS
LDFLAGSOPT   := -xS
endif
endif
endif
endif

# Intel Atom CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and
# SSSE3 instruction set support.
ifeq ($(call match,$(ID),(pc|pc64)-atom-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77 := -xSSE3_ATOM
CFLAGSOPTF90 := -xSSE3_ATOM
CFLAGSOPTC   := -xSSE3_ATOM
CFLAGSOPTCXX := -xSSE3_ATOM
LDFLAGSOPT   := -xSSE3_ATOM
else
CFLAGSOPTF77 := -xSSE3_ATOM
CFLAGSOPTF90 := -xSSE3_ATOM
CFLAGSOPTC   := -xSSE3_ATOM
CFLAGSOPTCXX := -xSSE3_ATOM
LDFLAGSOPT   := -xSSE3_ATOM
endif
endif

##############################################################################
# AMD x86/x86_64 CPUs
#
##############################################################################

# AMD K6 CPU with MMX instruction set support and improved versions of
# AMD K6 CPU with MMX and 3DNow! instruction set support,
# respectively.
ifeq ($(call match,$(ID),pc-k6-.*-intel-.*),yes)
CFLAGSOPTF77 := -march=anyx86 -m3dnow
CFLAGSOPTF90 := -march=anyx86 -m3dnow
CFLAGSOPTC   := -march=anyx86 -m3dnow
CFLAGSOPTCXX := -march=anyx86 -m3dnow
LDFLAGSOPT   := -march=anyx86 -m3dnow
endif

# AMD Athlon CPU with MMX, 3dNOW!, enhanced 3DNow! and SSE prefetch
# instructions support.
ifeq ($(call match,$(ID),pc-athlon-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -xSSE
CFLAGSOPTF90  := -xSSE
CFLAGSOPTC    := -xSSE
CFLAGSOPTCXX  := -xSSE
LDFLAGSOPT    := -xSSE
else
CFLAGSOPTF77  := -axK
CFLAGSOPTF90  := -axK
CFLAGSOPTC    := -axK
CFLAGSOPTCXX  := -axK
LDFLAGSOPT    := -axK
endif
endif

# Improved AMD Athlon CPU with MMX, 3DNow!, enhanced 3DNow! and full
# SSE instruction set support.
ifeq ($(call match,$(ID),pc-athlonxp-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -xSSE
CFLAGSOPTF90  := -xSSE
CFLAGSOPTC    := -xSSE
CFLAGSOPTCXX  := -xSSE
LDFLAGSOPT    := -xSSE
else
CFLAGSOPTF77  := -axK
CFLAGSOPTF90  := -axK
CFLAGSOPTC    := -axK
CFLAGSOPTCXX  := -axK
LDFLAGSOPT    := -axK
endif
endif

# Processors based on the AMD K8 core with x86-64 instruction set
# support, including the AMD Opteron, Athlon 64, and Athlon 64 FX
# processors. (This supersets MMX, SSE, SSE2, 3DNow!, enhanced 3DNow!
# and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(athlon64|athlon64x2|turion64|turion64x2)-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -xSSE2
CFLAGSOPTF90  := -xSSE2
CFLAGSOPTC    := -xSSE2
CFLAGSOPTCXX  := -xSSE2
LDFLAGSOPT    := -xSSE2
else
CFLAGSOPTF77  := -axP
CFLAGSOPTF90  := -axP
CFLAGSOPTC    := -axP
CFLAGSOPTCXX  := -axP
LDFLAGSOPT    := -axP
endif
endif

# Processors based on the AMD K8 core with x86-64 instruction set
# support, including the AMD Opteron, Athlon 64, and Athlon 64 FX
# processors. (This supersets MMX, SSE, SSE2, 3DNow!, enhanced 3DNow!
# and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(opteron|opteronx2)-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -xSSE2
CFLAGSOPTF90  := -xSSE2
CFLAGSOPTC    := -xSSE2
CFLAGSOPTCXX  := -xSSE2
LDFLAGSOPT    := -xSSE2
else
CFLAGSOPTF77  := -axW
CFLAGSOPTF90  := -axW
CFLAGSOPTC    := -axW
CFLAGSOPTCXX  := -axW
LDFLAGSOPT    := -axW
endif
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-phenom-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -xSSE4.1
CFLAGSOPTF90  := -xSSE4.1
CFLAGSOPTC    := -xSSE4.1
CFLAGSOPTCXX  := -xSSE4.1
LDFLAGSOPT    := -xSSE4.1
else
CFLAGSOPTF77  := -xT
CFLAGSOPTF90  := -xT
CFLAGSOPTC    := -xT
CFLAGSOPTCXX  := -xT
LDFLAGSOPT    := -xT
endif
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-phenomII-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -xSSE4.1
CFLAGSOPTF90  := -xSSE4.1
CFLAGSOPTC    := -xSSE4.1
CFLAGSOPTCXX  := -xSSE4.1
LDFLAGSOPT    := -xSSE4.1
else
CFLAGSOPTF77  := -xT
CFLAGSOPTF90  := -xT
CFLAGSOPTC    := -xT
CFLAGSOPTCXX  := -xT
LDFLAGS       := -xT
endif
endif

# CPUs based on AMD Family 10h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSE4A, 3DNow!,
# enhanced 3DNow!, ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-(opteron4100|opteron6100)-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -msse3
CFLAGSOPTF90  := -msse3
cflagsoptc    := -msse3
CFLAGSOPTCXX  := -msse3
LDFLAGSOPT    := -mSSE3
else
CFLAGSOPTF77  := -xT
CFLAGSOPTF90  := -xT
CFLAGSOPTC    := -xT
CFLAGSOPTCXX  := -xT
LDFLAGSOPT    := -xT
endif
endif

# AMD Family 15: Bulldozer core
# CPUs based on AMD Family 15h cores with x86-64 instruction set
# support. (This supersets FMA4, AVX, XOP, LWP, AES, PCL_MUL, CX16,
# MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2, ABM and 64-bit
# instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-bulldozer-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -msse4.2
CFLAGSOPTF90  := -msse4.2
CFLAGSOPTC    := -msse4.2
CFLAGSOPTCXX  := -msse4.2
LDFLAGSOPT    := -msse4.2
else
CFLAGSOPTF77  := -xT
CFLAGSOPTF90  := -xT
CFLAGSOPTC    := -xT
CFLAGSOPTCXX  := -xT
LDFLAGSOPT    := -xT
endif
endif

# AMD Family 15h: Piledriver core
# AMD Family 15h core based CPUs with x86-64 instruction set
# support. (This supersets BMI, TBM, F16C, FMA, AVX, XOP, LWP, AES,
# PCL_MUL, CX16, MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2,
# ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-piledriver-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -msse4.2
CFLAGSOPTF90  := -msse4.2
CFLAGSOPTC    := -msse4.2
CFLAGSOPTCXX  := -msse4.2
LDFLAGSOPT    := -msse4.2
else
CFLAGSOPTF77  := -xT
CFLAGSOPTF90  := -xT
CFLAGSOPTC    := -xT
CFLAGSOPTCXX  := -xT
LDFLAGSOPT    := -xT
endif
endif

# AMD Family 15h: Steamroller core
# AMD Family 15h core based CPUs with x86-64 instruction set
# support. (This supersets BMI, TBM, F16C, FMA, AVX, XOP, LWP, AES,
# PCL_MUL, CX16, MMX, SSE, SSE2, SSE3, SSE4A, SSSE3, SSE4.1, SSE4.2,
# ABM and 64-bit instruction set extensions.
ifeq ($(call match,$(ID),(pc|pc64)-steamroller-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -msse4.2
CFLAGSOPTF90  := -msse4.2
CFLAGSOPTC    := -msse4.2
CFLAGSOPTCXX  := -msse4.2
LDFLAGSOPT    := -msse4.2
else
CFLAGSOPTF77  := -xT
CFLAGSOPTF90  := -xT
CFLAGSOPTC    := -xT
CFLAGSOPTCXX  := -xT
LDFLAGSOPT    := -xT
endif
endif

# AMD Family 14h: Bobcat core
# CPUs based on AMD Family 14h cores with x86-64 instruction set
# support. (This supersets MMX, SSE, SSE2, SSE3, SSSE3, SSE4A, CX16,
# ABM and 64-bit instruction set extensions.)
ifeq ($(call match,$(ID),(pc|pc64)-bobcat-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -msse4.1
CFLAGSOPTF90  := -msse4.1
CFLAGSOPTC    := -msse4.1
CFLAGSOPTCXX  := -msse4.1
LDFLAGSOPT    := -msse4.1
else
CFLAGSOPTF77  := -xT
CFLAGSOPTF90  := -xT
CFLAGSOPTC    := -xT
CFLAGSOPTCXX  := -xT
LDFLAGSOPT    := -xT
endif
endif


# AMD Family 16h: Jaguar core
# CPUs based on AMD Family 16h cores with x86-64 instruction set
# support. This includes MOVBE, F16C, BMI, AVX, PCL_MUL, AES, SSE4.2,
# SSE4.1, CX16, ABM, SSE4A, SSSE3, SSE3, SSE2, SSE, MMX and 64-bit
# instruction set extensions.
ifeq ($(call match,$(ID),(pc|pc64)-jaguar-.*-intel-.*),yes)
ifeq ($(call intelminversion,11,1),yes)
CFLAGSOPTF77  := -msse4.2
CFLAGSOPTF90  := -msse4.2
CFLAGSOPTC    := -msse4.2
CFLAGSOPTCXX  := -msse4.2
LDFLAGSOPT    := -msse4.2
else
CFLAGSOPTF77  := -xT
CFLAGSOPTF90  := -xT
CFLAGSOPTC    := -xT
CFLAGSOPTCXX  := -xT
LDFLAGSOPT    := -xT
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
