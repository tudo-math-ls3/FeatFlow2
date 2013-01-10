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

#if 0
! External module files cannot be included by the standard use
! statement since the configure script would generate a rule for
! compiling the file mymod.mod from mymod.f90 which, does not exist.
#endif

#define FEAT2_PP_EXTERNAL_USE(module) use module


#if 0
!-------------------------------------------------------------------------------
! The following macros can be used to identify the operation system.
! The macros are based on the reference:
! http://nadeausoftware.com/articles/2012/01/c_c_tip_how_use_compiler_predefined_macros_detect_operating_system
!-------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------
! All current UNIX-style OSes, including BSD, Linux, OSX, and Solaris
!-------------------------------------------------------------------------------
#endif

#if (!defined(_WIN32) && (defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__))))
#define FEAT2_PP_OS_IS_UNIX() 1
#else
#define FEAT2_PP_OS_IS_UNIX() 0
#endif

#if 0
!-------------------------------------------------------------------------------
! The following macros can be used to check if the operating system is 32bit or 64bit
!-------------------------------------------------------------------------------
#endif

#if 0
!-------------------------------------------------------------------------------
! LP64 machine, OS X or Linux or Unix or LLP64 machine, Windows
!-------------------------------------------------------------------------------
#endif

#if (defined(_WIN64) || defined(_LP64) || defined(__LP64__))
#define FEAT2_PP_OS_IS_64BIT()   1
#else
#define FEAT2_PP_OS_IS_64BIT()   0
#endif

#if 0
!-------------------------------------------------------------------------------
! 32-bit machine, Windows or Linux or OS X or Unix
!-------------------------------------------------------------------------------
#endif

#define FEAT2_PP_OS_IS_32BIT() 0
#if !FEAT2_PP_OS_IS_64BIT()
#if (defined(_WIN32) || defined(_ILD32))
#undef FEAT2_PP_OS_IS_32BIT()
#define FEAT2_PP_OS_IS_32BIT() 1
#endif
#endif

#endif
