#ifndef _FEAT2CONSTANTS_H_
#define _FEAT2CONSTANTS_H_

#if 0
!##############################################################################
!# <name> feat2constants </name>
!# ****************************************************************************
!#
!# <purpose>
!# This header file defines a couple of preprocessor
!# constants used throught the complete Featflow2 kernel.
!# </purpose>
!##############################################################################
#endif


#if 0
!##############################################################################
! Undefine the keyword 'MYNEWLINE' (i.e. keep it) if this file is preprocessed 
! by f90cpp and simply remove it if another preprocessor is used instead.
!##############################################################################
#endif

#ifdef USE_PREPROC_F90CPP
#undef MYNEWLINE
#else
#define MYNEWLINE
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
! Supported precisions: Quad, Double, Single, Integer, I8, I16, I32, I64
!##############################################################################
#endif

#define QUAD_PREC    1
#define DOUBLE_PREC  2
#define SINGLE_PREC  3
#define INTEGER_PREC 4
#define INT8_PREC    5
#define INT16_PREC   6
#define INT32_PREC   7
#define INT64_PREC   8

      
#define QUAD_PREC_F    QP
#define DOUBLE_PREC_F  DP
#define SINGLE_PREC_F  SP

#define INT8_PREC_F    I8
#define INT16_PREC_F   I16
#define INT32_PREC_F   I32
#define INT64_PREC_F   I64
#ifdef ENABLE_LARGEINT
#define INTEGER_PREC_F I64
#else
#define INTEGER_PREC_F I32
#endif  

      
#define QUAD_PREC_C    long double
#define DOUBLE_PREC_C  double
#define SINGLE_PREC_C  float
#define INT8_PREC_C    signed char
#define INT16_PREC_C   short
#define INT32_PREC_C   int
#define INT64_PREC_C   long long int
#ifdef ENABLE_LARGEINT
#define INTEGER_PREC_C long long int
#else
#define INTEGER_PREC_C int
#endif     

#endif
