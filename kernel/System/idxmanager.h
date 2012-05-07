#ifndef _IDXMANAGER_H_
#define _IDXMANAGER_H_

#include "../feat2constants.h"

#if 0
!##############################################################################
!# <name> idxmanager </name>
!# ****************************************************************************
!#
!# <purpose>
!# This header file provides a transparent index manager which allows to
!# address multi-dimensional arrays in different programming languages in
!# a uniform manner. This is particularly useful if ...
!#
!# a) ... different programming languages are mixed in one application
!# b) ... different memory layouts should be tested without rewriting the code
!# c) ... many more nice features to come
!#
!# </purpose>
!##############################################################################
#endif


#if 0
!##############################################################################
! Check for prerequisites
!##############################################################################
#endif

#ifndef LANGUAGE
#eror "Constant LANGUAGE must be defined!"
#endif

#ifndef IDXADDR
#error "Constant IDXADDR must be defined!"
#endif


#if 0
!##############################################################################
! Fortran array indices typically start at position 1, whereas the first
! index of an array in C is 0. This mismatch can be compensated by setting
! __IDFOFFSET__ to +/-1 so that arrays can be addressed uniformly from C and Fortran.
! If both LANGUAGE and IDXADDR coincide, then __IDFOFFSET__ is set to 0, otherwise
! it is set to +1 or -1 depending on the combination LANGUAGE/IDXADDR.
!##############################################################################
#endif

#if LANGUAGE == IDXADDR
#define __IDFOFFSET__ 0
#else
#if (LANGUAGE == LANGUAGE_C) && (IDXADDR == LANGUAGE_F)
#define __IDFOFFSET__ -1
#elif (LANGUAGE == LANGUAGE_F) && (IDXADDR == LANGUAGE_C)
#define __IDFOFFSET__  1
#else
#error "Unsupported combination LANGUAGE/IDXADDR!"
#endif
#endif


#if 0
!##############################################################################
! Supported memory layouts: Column/Row Major Ordering
!
! Default memory layout for Fortran programming language is Column Major
! Ordering and Row Major Ordering for C.
!##############################################################################
#endif

#ifndef MEMORY_LAYOUT

#if LANGUAGE == LANGUAGE_F
#define MEMORY_LAYOUT COLUMN_MAJOR_ORDER
#elif LANGUAGE == LANGUAGE_C
#define MEMORY_LAYOUT ROW_MAJOR_ORDER
#else
#error "Unsupported programming language"
#endif

#endif


#if 0
!##############################################################################
! Low-level index wrapper macros for addressing multi-dimensional arrays
! from Fortran and C-code in a unified manner as array(IDXn(1,2,3)).
!##############################################################################
#endif

#if __IDFOFFSET__ == 0
#define __IDX1_C__(i1) i1
#define __IDX1_F__(i1) i1
#else
#define __IDX1_C__(i1) (i1+(__IDFOFFSET__))
#define __IDX1_F__(i1) (i1+(__IDFOFFSET__))
#endif

#if __IDFOFFSET__ == 0
#define __IDX2_C__(i1,i2,n1,n2) ((n2)*(i1)+(i2))
#define __IDX2_F__(i1,i2,n1,n2) i1,i2
#else
#define __IDX2_C__(i1,i2,n1,n2) ((n2)*(i1+(__IDFOFFSET__))+(i2)+(__IDFOFFSET__))
#define __IDX2_F__(i1,i2,n1,n2) i1+(__IDFOFFSET__),i2+(__IDFOFFSET__)
#endif

#if __IDFOFFSET__ == 0
#define __IDX3_C__(i1,i2,i3,n1,n2,n3) ((n3)*(n2)*(i1)+(n3)*(i2)+i3)
#define __IDX3_F__(i1,i2,i3,n1,n2,n3) i1,i2,i3
#else
#define __IDX3_C__(i1,i2,i3,n1,n2,n3) ((n3)*(n2)*(i1+(__IDFOFFSET__))+(n3)*(i2+(__IDFOFFSET__))+(i3)+(__IDFOFFSET__))
#define __IDX3_F__(i1,i2,i3,n1,n2,n3) i1+(__IDFOFFSET__),i2+(__IDFOFFSET__),i3+(__IDFOFFSET__)
#endif


#if 0
!##############################################################################
! Generic index wrapper macros for addressing multi-dimensional arrays
! from Fortran and C-code in a unified manner as array(IDXn(1,2,3)).
!
! These macros take into account:
! - the programming language: C or Fortran
! - the desired memory layout: row/column major order
!
! With these macros (multi-dimensional) arrays can be addressed in a unified
! manner both from C and Fortran. Consider the generic statement
!
!   array(IDX2(i,j,n,m))
!
! which is transformed into:
!
!   array(i,j)  if Fortran is used with column major order
!   array(j,i)  if Fortran is used with row major order
!   array[n*j+i-1]    if C is used with column major order
!   array[m*i+j-1]    if C is used with row major order
!
!##############################################################################
#endif

#if LANGUAGE == LANGUAGE_F

#define IDX1_FORWARD(a,i1)                a(__IDX1_F__(i1))
#define IDX2_FORWARD(a,i1,i2,n1,n2)       a(__IDX2_F__(i1,i2,n1,n2))
#define IDX3_FORWARD(a,i1,i2,i3,n1,n2,n3) a(__IDX3_F__(i1,i2,i3,n1,n2,n3))

#define IDX1_REVERSE(a,i1)                a(__IDX1_F__(i1))
#define IDX2_REVERSE(a,i1,i2,n1,n2)       a(__IDX2_F__(i2,i1,n2,n1))
#define IDX3_REVERSE(a,i1,i2,i3,n1,n2,n3) a(__IDX3_F__(i3,i2,i1,n3,n2,n1))

#if MEMORY_LAYOUT == COLUMN_MAJOR_ORDER

#define IDX1(a,i1)                IDX1_FORWARD(a,i1)
#define IDX2(a,i1,i2,n1,n2)       IDX2_FORWARD(a,i1,i2,n1,n2)
#define IDX3(a,i1,i2,i3,n1,n2,n3) IDX3_FORWARD(a,i1,i2,i3,n1,n2,n3)

#define IDX1T(a,i1)                IDX1_REVERSE(a,i1)
#define IDX2T(a,i1,i2,n1,n2)       IDX2_REVERSE(a,i1,i2,n1,n2)
#define IDX3T(a,i1,i2,i3,n1,n2,n3) IDX3_REVERSE(a,i1,i2,i3,n1,n2,n3)

#elif MEMORY_LAYOUT == ROW_MAJOR_ORDER

#define IDX1(a,i1)                IDX1_REVERSE(a,i1)
#define IDX2(a,i1,i2,n1,n2)       IDX2_REVERSE(a,i1,i2,n1,n2)
#define IDX3(a,i1,i2,i3,n1,n2,n3) IDX3_REVERSE(a,i1,i2,i3,n1,n2,n3)

#define IDX1T(a,i1)                IDX1_FORWARD(a,i1)
#define IDX2T(a,i1,i2,n1,n2)       IDX2_FORWARD(a,i1,i2,n1,n2)
#define IDX3T(a,i1,i2,i3,n1,n2,n3) IDX3_FORWARD(a,i1,i2,i3,n1,n2,n3)

#else
#error "Unsupported memory layout"
#endif

#elif LANGUAGE == LANGUAGE_C

#define IDX1_FORWARD(a,i1)                a[__IDX1_C__(i1)]
#define IDX2_FORWARD(a,i1,i2,n1,n2)       a[__IDX2_C__(i1,i2,n1,n2)]
#define IDX3_FORWARD(a,i1,i2,i3,n1,n2,n3) a[__IDX3_C__(i1,i2,i3,n1,n2,n3)]

#define IDX1_REVERSE(a,i1)                a[__IDX1_C__(i1)]
#define IDX2_REVERSE(a,i1,i2,n1,n2)       a[__IDX2_C__(i2,i1,n2,n1)]
#define IDX3_REVERSE(a,i1,i2,i3,n1,n2,n3) a[__IDX3_C__(i3,i2,i1,n3,n2,n1)]

#if MEMORY_LAYOUT == COLUMN_MAJOR_ORDER

#define IDX1(a,i1)                IDX1_REVERSE(a,i1)
#define IDX2(a,i1,i2,n1,n2)       IDX2_REVERSE(a,i1,i2,n1,n2)
#define IDX3(a,i1,i2,i3,n1,n2,n3) IDX3_REVERSE(a,i1,i2,i3,n1,n2,n3)

#define IDX1T(a,i1)                IDX1_FORWARD(a,i1)
#define IDX2T(a,i1,i2,n1,n2)       IDX2_FORWARD(a,i1,i2,n1,n2)
#define IDX3T(a,i1,i2,i3,n1,n2,n3) IDX3_FORWARD(a,i1,i2,i3,n1,n2,n3)

#elif MEMORY_LAYOUT == ROW_MAJOR_ORDER

#define IDX1(a,i1)                IDX1_FORWARD(a,i1)
#define IDX2(a,i1,i2,n1,n2)       IDX2_FORWARD(a,i1,i2,n1,n2)
#define IDX3(a,i1,i2,i3,n1,n2,n3) IDX3_FORWARD(a,i1,i2,i3,n1,n2,n3)

#define IDX1T(a,i1)                IDX1_REVERSE(a,i1)
#define IDX2T(a,i1,i2,n1,n2)       IDX2_REVERSE(a,i1,i2,n1,n2)
#define IDX3T(a,i1,i2,i3,n1,n2,n3) IDX3_REVERSE(a,i1,i2,i3,n1,n2,n3)

#else
#error "Unsupported memory layout"
#endif

#else
#error "Unsupported programming language"
#endif


#define IDX1_ALLOCATE(a,n1)       IDX1(a,n1)
#define IDX2_ALLOCATE(a,n1,n2)    IDX2(a,n1,n2,n1,n2)
#define IDX3_ALLOCATE(a,n1,n2,n3) IDX3(a,n1,n2,n3,n1,n2,n3)

#define IDX1_MALLOC(a,n1)         IDX1(a,n1)
#define IDX2_MALLOC(a,n1,n2)      IDX2(a,n1,n2,n1,n2)
#define IDX3_MALLOC(a,n1,n2,n3)   IDX3(a,n1,n2,n3,n1,n2,n3)

#endif
