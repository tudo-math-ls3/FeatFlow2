#ifndef _FEAT2CONSTANS_H_
#define _FEAT2CONSTANS_H_

#if 0
!-*- mode: f90; -*-
!##############################################################################
!# <name> feat2constants </name>
!# ****************************************************************************
!#
!# <purpose>
!# This header file provides global preprocessor constants.
!#
!# </purpose>
!##############################################################################
#endif


#if 0
!##############################################################################
! Supported programming languages: C or Fortran
!##############################################################################
#endif

#define LANGUAGE_C 1
#define LANGUAGE_F 2


#if 0
!##############################################################################
! Supported array index addressing: C- or Fortran-style
!##############################################################################
#endif

#define IDXADDR_C LANGUAGE_C
#define IDXADDR_F LANGUAGE_F


#if 0
!##############################################################################
! Supported memory layouts: Column/Row Major Ordering
!##############################################################################
#endif

#define COLUMN_MAJOR_ORDER 1
#define ROW_MAJOR_ORDER    2


#if 0
!##############################################################################
! Supported precisions: Quad, Double, Single, I8, I16, I32, I64
!##############################################################################
#endif

#define QUAD_PREC   1
#define DOUBLE_PREC 2
#define SINGLE_PREC 3
#define INT8_PREC   4
#define INT16_PREC  5
#define INT32_PREC  6
#define INT64_PREC  7

#endif
