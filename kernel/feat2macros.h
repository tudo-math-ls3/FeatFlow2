#ifndef _FEAT2MACROS_H_
#define _FEAT2MACROS_H_

#if 0
!##############################################################################
!# <name> feat2macros </name>
!# ****************************************************************************
!#
!# <purpose>
!# This header file defines a couple of preprocessor
!# macros used throught the complete Featflow2 kernel.
!# </purpose>
!##############################################################################
#endif

#include "kernel/feat2constants.h"

#if 0
!###############################################################################
! External module files cannot be included by the standard use
! statement since the configure script would generate a rule for
! compiling the file mymod.mod from mymod.f90 which, does not exist.
!###############################################################################
#endif

#define FEAT2_PP_EXTERNAL_USE(module) use module


#if 0
!###############################################################################
! The following macros can be used to identify the operation system.
! The macros are based on the reference:
! http://nadeausoftware.com/articles/2012/01/c_c_tip_how_use_compiler_predefined_macros_detect_operating_system
!###############################################################################
#endif

#if (defined(_AIX))
#define FEAT2_PP_OS_IS_AIX()     1
#else
#define FEAT2_PP_OS_IS_AIX()     0
#endif

#if (defined(__CYGWIN__) && !defined(_WIN32))
#define FEAT2_PP_OS_IS_CYGWIN()  1
#else
#define FEAT2_PP_OS_IS_CYGWIN()  0
#endif

#if (defined(__hpux))
#define FEAT2_PP_OS_IS_HPUX()    1
#else
#define FEAT2_PP_OS_IS_HPUX()    0
#endif

#if (defined(__linux__))
#define FEAT2_PP_OS_IS_LINUX()   1
#else
#define FEAT2_PP_OS_IS_LINUX()   0
#endif

#if (defined(__APPLE__) && defined(__MACH__))
#define FEAT2_PP_OS_IS_OSX()     1
#else
#define FEAT2_PP_OS_IS_OSX()     0
#endif

#if (defined(__sun) && defined(__SVR4))
#define FEAT2_PP_OS_IS_SOLARIS() 1
#else
#define FEAT2_PP_OS_IS_SOLARIS() 0
#endif

#if (defined(_WIN32))
#define FEAT2_PP_OS_IS_WIN()     1
#else
#define FEAT2_PP_OS_IS_WIN()     0
#endif

#if (defined(_WIN32) && !defined(_WIN64))
#define FEAT2_PP_OS_IS_WIN32()   1
#else
#define FEAT2_PP_OS_IS_WIN32()   0
#endif

#if (defined(_WIN64))
#define FEAT2_PP_OS_IS_WIN64()   1
#else
#define FEAT2_PP_OS_IS_WIN64()   0
#endif

#if 0
!###############################################################################
! All current UNIX-style OSes, including BSD, Linux, OSX, and Solaris
!###############################################################################
#endif

#if (!defined(_WIN32) && (defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__))))
#define FEAT2_PP_OS_IS_UNIX() 1
#else
#define FEAT2_PP_OS_IS_UNIX() 0
#endif

#if 0
!###############################################################################
! The following macros can be used to check if the operating system is 32bit or 64bit
!###############################################################################
#endif

#if 0
!###############################################################################
! LP64 machine, OS X or Linux or Unix or LLP64 machine, Windows
!###############################################################################
#endif

#if (defined(_WIN64) || defined(_LP64) || defined(__LP64__))
#define FEAT2_PP_OS_IS_64BIT()   1
#else
#define FEAT2_PP_OS_IS_64BIT()   0
#endif

#if 0
!###############################################################################
! 32-bit machine, Windows or Linux or OS X or Unix
!###############################################################################
#endif

#define FEAT2_PP_OS_IS_32BIT() 0
#if !FEAT2_PP_OS_IS_64BIT()
#if (defined(_WIN32) || defined(_ILD32))
#undef FEAT2_PP_OS_IS_32BIT
#define FEAT2_PP_OS_IS_32BIT() 1
#endif
#endif


#if 0
!###############################################################################
! The following macro tries to detect the programming language automatically
!
! Example: FEAT2_PP_AUTO_LANGUAGE()
!          expands to the content of the global variable LANGUAGE if defined;
!          if variable LANGUAGE is not set then we try to identify the
!          programming language by checking various preprocessor macros which
!          are typically predefined by Fortran compilers.
!          It is also possible to overwrite the global LANGUAGE variable by
!          hand by passing the optional argument language.
!###############################################################################
#endif

#define FEAT2_PP_AUTO_LANGUAGE(language) FEAT2_PP_GET_LANGUAGE_OVERRIDE(language)

#if 0
!-------------------------------------------------------------------------------
! Auxiliary macros that should not be called directly
!-------------------------------------------------------------------------------
#endif

#define FEAT2_PP_GET_LANGUAGE_OVERRIDE(language) FEAT2_PP_GET_LANGUAGE_##language()
#define FEAT2_PP_GET_LANGUAGE_1()                LANGUAGE_C
#define FEAT2_PP_GET_LANGUAGE_2()                LANGUAGE_F
#define FEAT2_PP_GET_LANGUAGE_()                 FEAT2_PP_GET_LANGUAGE(LANGUAGE)
#define FEAT2_PP_GET_LANGUAGE(language)          FEAT2_PP_GET_LANGUAGE_I(language)
#define FEAT2_PP_GET_LANGUAGE_I(language)        FEAT2_PP_GET_LANGUAGE_##language()

#if defined(USE_PREPROC_F90CPP) ||                     \
    defined(__GFORTRAN__) || defined(__G95__) ||       \
    defined(_LANGUAGE_FORTRAN) ||                      \
    defined(__SUNPRO_F90) || defined (__SUNPRO_F95) || \
    defined(__INTEL_COMPILER) && !defined(__ICC)
#define FEAT2_PP_GET_LANGUAGE_LANGUAGE()  LANGUAGE_F
#else
#define FEAT2_PP_GET_LANGUAGE_LANGUAGE()  LANGUAGE_C
#endif


#if 0
!###############################################################################
! The following macro expands a given constant with the correct data type
!
! Example: FEAT2_PP_CONST(1.0, DOUBLE_PREC)
!          expands to 1.0_DP in Fortran-mode and 1.0E in C-mode
!
!          FEAT2_PP_CONST(2.0, SINGLE_PREC)
!          expands to 2.0_SP in Fortran-mode and 2.0F in C-mode
!###############################################################################
#endif

#define FEAT2_PP_CONST(value,precision)   FEAT2_PP_CONST_LANG(value,precision,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_CONST_C(value,precision) FEAT2_PP_CONST_LANG(value,precision,LANGUAGE_C)
#define FEAT2_PP_CONST_F(value,precision) FEAT2_PP_CONST_LANG(value,precision,LANGUAGE_F)

#if 0
!-------------------------------------------------------------------------------
! Auxiliary macros that should not be called directly
!-------------------------------------------------------------------------------
#endif

#define FEAT2_PP_CONST_LANG(value,precision,language)   FEAT2_PP_CONST_LANG_I(value,precision,language)
#define FEAT2_PP_CONST_LANG_I(value,precision,language) FEAT2_PP_CONST_##precision##_##language(value)

#define FEAT2_PP_CONST_1_1(value) value##l
#define FEAT2_PP_CONST_2_1(value) value
#define FEAT2_PP_CONST_3_1(value) value##f
#ifdef ENABLE_LARGEINT
#define FEAT2_PP_CONST_4_1(value) value##l
#else
#define FEAT2_PP_CONST_4_1(value) value
#endif
#define FEAT2_PP_CONST_5_1(value) value
#define FEAT2_PP_CONST_6_1(value) value
#define FEAT2_PP_CONST_7_1(value) value
#define FEAT2_PP_CONST_8_1(value) value##l

#define FEAT2_PP_CONST_CONCAT_F_I(value,prec_f) value##_##prec_f
#define FEAT2_PP_CONST_CONCAT_F(value,prec_f) FEAT2_PP_CONST_CONCAT_F_I(value,prec_f)

#define FEAT2_PP_CONST_1_2(value) FEAT2_PP_CONST_CONCAT_F(value,QUAD_PREC_F)
#define FEAT2_PP_CONST_2_2(value) FEAT2_PP_CONST_CONCAT_F(value,DOUBLE_PREC_F)
#define FEAT2_PP_CONST_3_2(value) FEAT2_PP_CONST_CONCAT_F(value,SINGLE_PREC_F)
#define FEAT2_PP_CONST_4_2(value) FEAT2_PP_CONST_CONCAT_F(value,INTEGER_PREC_F)
#define FEAT2_PP_CONST_5_2(value) FEAT2_PP_CONST_CONCAT_F(value,INT8_PREC_F)
#define FEAT2_PP_CONST_6_2(value) FEAT2_PP_CONST_CONCAT_F(value,INT16_PREC_F)
#define FEAT2_PP_CONST_7_2(value) FEAT2_PP_CONST_CONCAT_F(value,INT32_PREC_F)
#define FEAT2_PP_CONST_8_2(value) FEAT2_PP_CONST_CONCAT_F(value,INT64_PREC_F)


#if 0
!###############################################################################
! The following macros decodes one or multiple integers into one identifier
!
! Example: FEAT2_PP_ID3(2,5,7,10)  ->  2 + 10*5 + 7*10*10 = 752
!
! Note that it is of course possible to use bases other than 10.
!###############################################################################
#endif

#define FEAT2_PP_ID_I(id,base) (id*base)

#define FEAT2_PP_ID1(id1,base)                 (id1)
#define FEAT2_PP_ID2(id1,id2,base)             ((id1)+FEAT2_PP_ID_I(id2,base))
#define FEAT2_PP_ID3(id1,id2,id3,base)         FEAT2_PP_ID2(id1,FEAT2_PP_ID2(id2,id3,base),base)
#define FEAT2_PP_ID4(id1,id2,id3,id4,base)     FEAT2_PP_ID2(id1,FEAT2_PP_ID3(id2,id3,id4,base),base)
#define FEAT2_PP_ID5(id1,id2,id3,id4,id5,base) FEAT2_PP_ID2(id1,FEAT2_PP_ID4(id2,id3,id4,id5,base),base)


#if 0
!###############################################################################
! The following macro concatenates two tokens
!
! Example: FEAT2_PP_CONCAT(prefix_,suffix)
!          expands to prefix_suffix
!
! Note that use of this macro is recommended since '##' is not supported by all
! compiler suites, e.g., gfortran/g95 if invoked with -cpp flag internally call
! cpp in traditional mode which does not support concatenation using '##'.
!###############################################################################
#endif

#define FEAT2_PP_CONCAT(prefix,suffix) FEAT2_PP_CONCAT_I(prefix,suffix)

#if (defined(__GFORTRAN__) || defined(__G95__) || defined(_LANGUAGE_FORTRAN)) && !defined(__STDC__)
#define FEAT2_PP_CONCAT_I(prefix,suffix) prefix/**/suffix
#else
#define FEAT2_PP_CONCAT_I(prefix,suffix) prefix##suffix
#endif


#if 0
!###############################################################################
! The following macro concatenates two tokens
!
! Example: FEAT2_PP_STRING(string)
!          expands to "string"
!
! Note that use of this macro is recommended since '#' is not supported by all
! compiler suites, e.g., gfortran/g95 if invoked with -cpp flag internally call
! cpp in traditional mode which does not support stringification using '#'.
!###############################################################################
#endif

#define FEAT2_PP_STRING(string) FEAT2_PP_STRING_I(string)

#define FEAT2_PP_STRING_I(string) FEAT2_PP_STRING_II(string)
#if (defined(__GFORTRAN__) || defined(__G95__) || defined(_LANGUAGE_FORTRAN)) && !defined(__STDC__)
#define FEAT2_PP_STRING_II(string) "string"
#else
#define FEAT2_PP_STRING_II(string) #string
#endif

#endif
