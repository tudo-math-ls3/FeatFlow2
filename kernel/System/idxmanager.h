#ifndef _IDXMANAGER_H_
#define _IDXMANAGER_H_

#if 0
!-*- mode: f90; -*-
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
! Supported programming languages: C or Fortran
!
! Default programming language is Fortran
!##############################################################################
#endif

#define LANGUAGE_C 1
#define LANGUAGE_F 2

#ifndef LANGUAGE
#define LANGUAGE LANGUAGE_F
#endif

#if 0
! Fortran array indices typically start at position 1, whereas the first
! index of an array in C is 0. This mismatch can be compensated by setting
! IOFFSET to -1 so that arrays can be addressed uniformly starting at 1.
! If C-programmers wish to adress arrays starting at position 0, then
! IOFFSET can be set to 0 before this header file is included.
#endif

#if LANGUAGE == LANGUAGE_C
#ifndef IOFFSET
#define IOFFSET (-1)
#endif
#endif


#if 0
!##############################################################################
! Supported memory layouts: Column/Row Major Ordering
!
! Default memory layout for Fortran programming language is Column Major
! Ordering and Row Major Ordering for C. If no programming language is defined
! then Row Major Ordering is used as default memory layout.
!##############################################################################
#endif

#define COLUMN_MAJOR_ORDER 1
#define ROW_MAJOR_ORDER    2

#ifndef MEMORY_LAYOUT
#ifdef LANGUAGE
#if LANGUAGE == LANGUAGE_C
#define MEMORY_LAYOUT ROW_MAJOR_ORDER
#elif LANGUAGE == LANGUAGE_F
#define MEMORY_LAYOUT COLUMN_MAJOR_ORDER
#else
#error "Unsupported programming language"
#endif
#else
#define MEMORY_LAYOUT COLUMN_MAJOR_ORDER
#endif
#endif


#if 0
!##############################################################################
! Low-level index wrapper macros for addressing multi-dimensional arrays
! from Fortran and C-code in a unified manner as array(IDXn(1,2,3)).
!##############################################################################
#endif

#define IDX1_FORWARD_C(i1) (i1+IOFFSET)
#define IDX1_REVERSE_C(i1) (i1+IOFFSET)
#define IDX1_FORWARD_F(i1) i1
#define IDX1_REVERSE_F(i1) i1

#define IDX2_FORWARD_C(i1,i2,n1,n2) ((n2)*(i1+IOFFSET)+i2+IOFFSET)
#define IDX2_REVERSE_C(i1,i2,n1,n2) ((n1)*(i2+IOFFSET)+i1+IOFFSET)
#define IDX2_FORWARD_F(i1,i2,n1,n2) i1,i2
#define IDX2_REVERSE_F(i1,i2,n1,n2) i2,i1

#define IDX3_FORWARD_C(i1,i2,i3,n1,n2,n3) ((n3)*(n2)*(i1+IOFFSET)+(n3)*(i2+IOFFSET)+i3+IOFFSET)
#define IDX3_REVERSE_C(i1,i2,i3,n1,n2,n3) ((n1)*(n2)*(i3+IOFFSET)+(n1)*(i2+IOFFSET)+i1+IOFFSET)
#define IDX3_FORWARD_F(i1,i2,i3,n1,n2,n3) i1,i2,i3
#define IDX3_REVERSE_F(i1,i2,i3,n1,n2,n3) i3,i2,i1


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

#if MEMORY_LAYOUT == COLUMN_MAJOR_ORDER

#if LANGUAGE == LANGUAGE_F

#define IDX1(i1)                IDX1_FORWARD_F(i1)
#define IDX2(i1,i2,n1,n2)       IDX2_FORWARD_F(i1,i2,n1,n2)
#define IDX3(i1,i2,i3,n1,n2,n3) IDX3_FORWARD_F(i1,i2,i3,n1,n2,n3)

#elif LANGUAGE == LANGUAGE_C

#define IDX1(i1)                IDX1_REVERSE_C(i1)
#define IDX2(i1,i2,n1,n2)       IDX2_REVERSE_C(i1,i2,n1,n2)
#define IDX3(i1,i2,i3,n1,n2,n3) IDX3_REVERSE_C(i1,i2,i3,n1,n2,n3)

#else
#error "Unsupported programming language"
#endif

#elif MEMORY_LAYOUT == ROW_MAJOR_ORDER

#if LANGUAGE == LANGUAGE_F

#define IDX1(i1)                IDX1_REVERSE_F(i1)
#define IDX2(i1,i2,n1,n2)       IDX2_REVERSE_F(i1,i2,n1,n2)
#define IDX3(i1,i2,i3,n1,n2,n3) IDX3_REVERSE_F(i1,i2,i3,n1,n2,n3)

#elif LANGUAGE == LANGUAGE_C

#define IDX1(i1)                IDX1_FORWARD_C(i1)
#define IDX2(i1,i2,n1,n2)       IDX2_FORWARD_C(i1,i2,n1,n2)
#define IDX3(i1,i2,i3,n1,n2,n3) IDX3_FORWARD_C(i1,i2,i3,n1,n2,n3)

#else
#error "Unsupported programming language"
#endif

#else
#error "Unsupported memory layout"
#endif

#define IDX1_ALLOCATE(n1)       IDX1(n1)
#define IDX2_ALLOCATE(n1,n2)    IDX2(n1,n2,n1,n2)
#define IDX3_ALLOCATE(n1,n2,n3) IDX3(n1,n2,n3,n1,n2,n3)

#define IDX1_MALLOC(n1)       IDX1(n1)+1
#define IDX2_MALLOC(n1,n2)    IDX2(n1,n2,n1,n2)+1
#define IDX3_MALLOC(n1,n2,n3) IDX3(n1,n2,n3,n1,n2,n3)+1

#endif
