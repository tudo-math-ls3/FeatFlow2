!-*- mode: f90; -*-

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

#endif
