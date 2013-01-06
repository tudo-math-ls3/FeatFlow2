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

#define FEAT2_PP_OS_IS_AIX()     (defined(_AIX))
#define FEAT2_PP_OS_IS_CYGWIN()  (defined(__CYGWIN__) && !defined(_WIN32))
#define FEAT2_PP_OS_IS_HPUX()    (defined(__hpux))
#define FEAT2_PP_OS_IS_LINUX()   (defined(__linux__))
#define FEAT2_PP_OS_IS_OSX()     (defined(__APPLE__) && defined(__MACH__))
#define FEAT2_PP_OS_IS_SOLARIS() (defined(__sun) && defined(__SVR4))

#define FEAT2_PP_OS_IS_WIN()     (defined(_WIN32))
#define FEAT2_PP_OS_IS_WIN32()   (defined(_WIN32) && !defined(_WIN64))
#define FEAT2_PP_OS_IS_WIN64()   (defined(_WIN64))

#if 0
!-------------------------------------------------------------------------------
! All current UNIX-style OSes, including BSD, Linux, OSX, and Solaris
!-------------------------------------------------------------------------------
#endif

#define FEAT2_PP_OS_IS_UNIX()    (!defined(_WIN32) && (defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__))))


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
#define FEAT2_PP_OS_IS_64BIT()   (defined(_WIN64) || defined(_LP64) || defined(__LP64__))

#if 0
!-------------------------------------------------------------------------------
! 32-bit machine, Windows or Linux or OS X or Unix
!-------------------------------------------------------------------------------
#endif
#define FEAT2_PP_OS_IS_32BIT()   (!FEAT2_PP_OS_IS_64BIT() && (defined(_WIN32) || defined(_ILD32)))
#endif
