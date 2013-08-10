#ifndef _IDXMANAGER_H_
#define _IDXMANAGER_H_

#include "kernel/feat2constants.h"
#include "kernel/feat2macros.h"

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
! Supported array index addressing: C- or Fortran-style
!##############################################################################
#endif

#define IDXADDR_C LANGUAGE_C
#define IDXADDR_F LANGUAGE_F


#if 0
!###############################################################################
! The following macro tries to detect the array index addressing automatically
!
! Example: AUTO_IDXADDR()
!          expands to the content of variable IDXADDR if defined;
!          if variable IDXADDR is not set then the default value depends on the
!          programming language which is determined via FEAT2_PP_AUTO_LANGUAGE.
!          It is possible to overwrite the default behaviour of the language
!          auto-detection macro by passing the optioanl argument language.
!###############################################################################
#endif

#define AUTO_IDXADDR(language) GET_IDXADDR(IDXADDR,language)

#if 0
!-------------------------------------------------------------------------------
! Auxiliary macros that should not be called directly
!-------------------------------------------------------------------------------
#endif

#define GET_IDXADDR(idxaddr,language)   GET_IDXADDR_I(idxaddr,language)
#define GET_IDXADDR_I(idxaddr,language) GET_IDXADDR_##idxaddr(language)

#define GET_IDXADDR_1(language)        IDXADDR_C
#define GET_IDXADDR_2(language)        IDXADDR_F
#define GET_IDXADDR_IDXADDR(language)  FEAT2_PP_AUTO_LANGUAGE(language)


#if 0
!##############################################################################
! Supported memory layouts: Column/Row Major Ordering
!##############################################################################
#endif

#define ROW_MAJOR_ORDER    LANGUAGE_C
#define COLUMN_MAJOR_ORDER LANGUAGE_F


#if 0
!###############################################################################
! The following macro tries to detect the memory layout automatically
!
! Example: AUTO_MEMORY_LAYOUT()
!          expands to the content of variable MEMORY_LAYOUT if defined;
!          if variable MEMORY_LAYOUT is not set then the default value depends
!          on the programming language which is determined via
!          FEAT2_PP_AUTO_LANGUAGE.
!          It is possible to overwrite the default behaviour of the language
!          auto-detection macro by passing the optioanl argument language.
!###############################################################################
#endif

#define AUTO_MEMORY_LAYOUT(language) GET_MEMORY_LAYOUT(MEMORY_LAYOUT,language)

#if 0
!-------------------------------------------------------------------------------
! Auxiliary macros that should not be called directly
!-------------------------------------------------------------------------------
#endif

#define GET_MEMORY_LAYOUT(memorylayout,language)   GET_MEMORY_LAYOUT_I(memorylayout,language)
#define GET_MEMORY_LAYOUT_I(memorylayout,language) GET_MEMORY_LAYOUT_##memorylayout(language)

#define GET_MEMORY_LAYOUT_1(language)             ROW_MAJOR_ORDER
#define GET_MEMORY_LAYOUT_2(language)             COLUMN_MAJOR_ORDER
#define GET_MEMORY_LAYOUT_MEMORY_LAYOUT(language) FEAT2_PP_AUTO_LANGUAGE(language)


#if 0
!##############################################################################
! Fortran array indices typically start at position 1, whereas the
! first index of an array in C is 0. This mismatch can be compensated
! by using the AUTO_IDXOFFSET macro which computes the index (-1,0,1)
! based on the programming language and the index addressing scheme.
!
! Example: AUTO_IDXOFFSET()
!          expands to MINUS_ONE, PLUS_ONE, or ZERO, whereby the programming
!          language and the index addressing scheme are detected automatically
!          using the macros FEAT2_PP_AUTO_LANGUAGE and AUTO_IDXADDR, resp.
!
!          GET_IDXOFFSET(language,idxaddr)
!          expands to MINUS_ONE, PLUS_ONE, or ZERO, whereby the programming
!          language and the index addressing scheme are passed as arguments.
!###############################################################################
#endif

#define AUTO_IDXOFFSET()                GET_IDXOFFSET(FEAT2_PP_AUTO_LANGUAGE(),AUTO_IDXADDR())
#define GET_IDXOFFSET(language,idxaddr) GET_IDXOFFSET_I(language,idxaddr)

#if 0
!-------------------------------------------------------------------------------
! Auxiliary macros that should not be called directly
!-------------------------------------------------------------------------------
#endif

#define GET_IDXOFFSET_I(language,idxaddr) GET_IDXOFFSET_##language##_##idxaddr

#define GET_IDXOFFSET_1_1 ZERO
#define GET_IDXOFFSET_1_2 MINUS_ONE
#define GET_IDXOFFSET_2_1 PLUS_ONE
#define GET_IDXOFFSET_2_2 ZERO


#if 0
!##############################################################################
! Low-level index wrapper macros for addressing multi-dimensional arrays
! from Fortran and C-code in a unified manner as array(IDXn(1,2,3)).
!
! NOTE: Do not use these macros directly unless you really know what you are
!       doing! Instead, make use of the IDX{1,2,3}[T] macros defined below.
!##############################################################################
#endif

#define IDX1_C(i1) GET_IDX1_C(i1,GET_IDXOFFSET(LANGUAGE_C,AUTO_IDXADDR(LANGUAGE_C)))
#define IDX1_F(i1) GET_IDX1_F(i1,GET_IDXOFFSET(LANGUAGE_F,AUTO_IDXADDR(LANGUAGE_F)))

#if 0
!-------------------------------------------------------------------------------
#endif

#define GET_IDX1_C(i1,idxoffset)   GET_IDX1_C_I(i1,idxoffset)
#define GET_IDX1_C_I(i1,idxoffset) GET_IDX1_C_##idxoffset(i1)

#define GET_IDX1_C_ZERO(i1)      (i1)
#define GET_IDX1_C_PLUS_ONE(i1)  (i1+1)
#define GET_IDX1_C_MINUS_ONE(i1) (i1-1)


#define GET_IDX1_F(i1,idxoffset)   GET_IDX1_F_I(i1,idxoffset)
#define GET_IDX1_F_I(i1,idxoffset) GET_IDX1_F_##idxoffset(i1)

#define GET_IDX1_F_ZERO(i1)      i1
#define GET_IDX1_F_PLUS_ONE(i1)  i1+1
#define GET_IDX1_F_MINUS_ONE(i1) i1-1

#if 0
!-------------------------------------------------------------------------------
#endif

#define IDX2_C(i1,i2,n1,n2) GET_IDX2_C(i1,i2,n1,n2,GET_IDXOFFSET(LANGUAGE_C,AUTO_IDXADDR(LANGUAGE_C)))
#define IDX2_F(i1,i2,n1,n2) GET_IDX2_F(i1,i2,n1,n2,GET_IDXOFFSET(LANGUAGE_F,AUTO_IDXADDR(LANGUAGE_F)))

#if 0
!-------------------------------------------------------------------------------
#endif

#define GET_IDX2_C(i1,i2,n1,n2,idxoffset)   GET_IDX2_C_I(i1,i2,n1,n2,idxoffset)
#define GET_IDX2_C_I(i1,i2,n1,n2,idxoffset) GET_IDX2_C_##idxoffset(i1,i2,n1,n2)

#define GET_IDX2_C_ZERO(i1,i2,n1,n2)        ((n2)*(i1)+(i2))
#define GET_IDX2_C_PLUS_ONE(i1,i2,n1,n2)    ((n2)*(i1+1)+(i2+1))
#define GET_IDX2_C_MINUS_ONE(i1,i2,n1,n2)   ((n2)*(i1-1)+(i2-1))


#define GET_IDX2_F(i1,i2,n1,n2,idxoffset)   GET_IDX2_F_I(i1,i2,n1,n2,idxoffset)
#define GET_IDX2_F_I(i1,i2,n1,n2,idxoffset) GET_IDX2_F_##idxoffset(i1,i2,n1,n2)

#define GET_IDX2_F_ZERO(i1,i2,n1,n2)        i1,i2
#define GET_IDX2_F_PLUS_ONE(i1,i2,n1,n2)    i1+1,i2+1
#define GET_IDX2_F_MINUS_ONE(i1,i2,n1,n2)   i1-1,i2-1

#if 0
!-------------------------------------------------------------------------------
#endif

#define IDX3_C(i1,i2,i3,n1,n2,n3) GET_IDX3_C(i1,i2,i3,n1,n2,n3,GET_IDXOFFSET(LANGUAGE_C,AUTO_IDXADDR(LANGUAGE_C)))
#define IDX3_F(i1,i2,i3,n1,n2,n3) GET_IDX3_F(i1,i2,i3,n1,n2,n3,GET_IDXOFFSET(LANGUAGE_F,AUTO_IDXADDR(LANGUAGE_F)))

#if 0
!-------------------------------------------------------------------------------
#endif

#define GET_IDX3_C(i1,i2,i3,n1,n2,n3,idxoffset)   GET_IDX3_C_I(i1,i2,i3,n1,n2,n3,idxoffset)
#define GET_IDX3_C_I(i1,i2,i3,n1,n2,n3,idxoffset) GET_IDX3_C_##idxoffset(i1,i2,i3,n1,n2,n3)

#define GET_IDX3_C_ZERO(i1,i2,i3,n1,n2,n3)        ((n3)*(n2)*(i1)+(n3)*(i2)+(i3))
#define GET_IDX3_C_PLUS_ONE(i1,i2,i3,n1,n2,n3)    ((n3)*(n2)*(i1+1)+(n3)*(i2+1)+(i3+1))
#define GET_IDX3_C_MINUS_ONE(i1,i2,i3,n1,n2,n3)   ((n3)*(n2)*(i1-1)+(n3)*(i2-1)+(i3-1))


#define GET_IDX3_F(i1,i2,i3,n1,n2,n3,idxoffset)   GET_IDX3_F_I(i1,i2,i3,n1,n2,n3,idxoffset)
#define GET_IDX3_F_I(i1,i2,i3,n1,n2,n3,idxoffset) GET_IDX3_F_##idxoffset(i1,i2,i3,n1,n2,n3)

#define GET_IDX3_F_ZERO(i1,i2,i3,n1,n2,n3)        i1,i2,i3
#define GET_IDX3_F_PLUS_ONE(i1,i2,i3,n1,n2,n3)    i1+1,i2+1,i3+1
#define GET_IDX3_F_MINUS_ONE(i1,i2,i3,n1,n2,n3)   i1-1,i2-1,i3-1


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

#if 0
!##############################################################################
! 1D-arrays: IDX1, IDX1T, IDX1_FORWARD, IDX1_REVERSE
!
! All other macros are auxiliary ones and should not be called directly
!##############################################################################
#endif

#define IDX1_FORWARD(a,i1)             IDX1_FORWARD_I(a,i1,FEAT2_PP_AUTO_LANGUAGE())
#define IDX1_FORWARD_I(a,i1,language)  IDX1_FORWARD_II(a,i1,language)
#define IDX1_FORWARD_II(a,i1,language) IDX1_FORWARD_##language(a,i1)
#define IDX1_FORWARD_1(a,i1)           a[IDX1_C(i1)]
#define IDX1_FORWARD_2(a,i1)           a(IDX1_F(i1))

#define IDX1_REVERSE(a,i1)             IDX1_REVERSE_I(a,i1,FEAT2_PP_AUTO_LANGUAGE())
#define IDX1_REVERSE_I(a,i1,language)  IDX1_REVERSE_II(a,i1,language)
#define IDX1_REVERSE_II(a,i1,language) IDX1_REVERSE_##language(a,i1)
#define IDX1_REVERSE_1(a,i1)           a[IDX1_C(i1)]
#define IDX1_REVERSE_2(a,i1)           a(IDX1_F(i1))


#define IDX1(a,i1)                             GET_IDX1(a,i1,FEAT2_PP_AUTO_LANGUAGE(),AUTO_MEMORY_LAYOUT())
#define GET_IDX1(a,i1,language,memorylayout)   GET_IDX1_I(a,i1,language,memorylayout)
#define GET_IDX1_I(a,i1,language,memorylayout) GET_IDX1_##language##_##memorylayout(a,i1)

#define IDX1T(a,i1)                             GET_IDX1T(a,i1,FEAT2_PP_AUTO_LANGUAGE(),AUTO_MEMORY_LAYOUT())
#define GET_IDX1T(a,i1,language,memorylayout)   GET_IDX1T_I(a,i1,language,memorylayout)
#define GET_IDX1T_I(a,i1,language,memorylayout) GET_IDX1T_##language##_##memorylayout(a,i1)

#define GET_IDX1_1_1(a,i1) IDX1_FORWARD(a,i1)
#define GET_IDX1_1_2(a,i1) IDX1_REVERSE(a,i1)
#define GET_IDX1_2_1(a,i1) IDX1_REVERSE(a,i1)
#define GET_IDX1_2_2(a,i1) IDX1_FORWARD(a,i1)

#define GET_IDX1T_1_1(a,i1) IDX1_REVERSE(a,i1)
#define GET_IDX1T_1_2(a,i1) IDX1_FORWARD(a,i1)
#define GET_IDX1T_2_1(a,i1) IDX1_FORWARD(a,i1)
#define GET_IDX1T_2_2(a,i1) IDX1_REVERSE(a,i1)

#if 0
!##############################################################################
! 2D-arrays: IDX2, IDX2T, IDX2_FORWARD, IDX2_REVERSE
!
! All other macros are auxiliary ones and should not be called directly
!##############################################################################
#endif

#define IDX2_FORWARD(a,i1,i2,n1,n2)             IDX2_FORWARD_I(a,i1,i2,n1,n2,FEAT2_PP_AUTO_LANGUAGE())
#define IDX2_FORWARD_I(a,i1,i2,n1,n2,language)  IDX2_FORWARD_II(a,i1,i2,n1,n2,language)
#define IDX2_FORWARD_II(a,i1,i2,n1,n2,language) IDX2_FORWARD_##language(a,i1,i2,n1,n2)
#define IDX2_FORWARD_1(a,i1,i2,n1,n2)           a[IDX2_C(i1,i2,n1,n2)]
#define IDX2_FORWARD_2(a,i1,i2,n1,n2)           a(IDX2_F(i1,i2,n1,n2))

#define IDX2_REVERSE(a,i1,i2,n1,n2)             IDX2_REVERSE_I(a,i1,i2,n1,n2,FEAT2_PP_AUTO_LANGUAGE())
#define IDX2_REVERSE_I(a,i1,i2,n1,n2,language)  IDX2_REVERSE_II(a,i1,i2,n1,n2,language)
#define IDX2_REVERSE_II(a,i1,i2,n1,n2,language) IDX2_REVERSE_##language(a,i1,i2,n1,n2)
#define IDX2_REVERSE_1(a,i1,i2,n1,n2)           a[IDX2_C(i2,i1,n2,n1)]
#define IDX2_REVERSE_2(a,i1,i2,n1,n2)           a(IDX2_F(i2,i1,n2,n1))


#define IDX2(a,i1,i2,n1,n2)                             GET_IDX2(a,i1,i2,n1,n2,FEAT2_PP_AUTO_LANGUAGE(),AUTO_MEMORY_LAYOUT())
#define GET_IDX2(a,i1,i2,n1,n2,language,memorylayout)   GET_IDX2_I(a,i1,i2,n1,n2,language,memorylayout)
#define GET_IDX2_I(a,i1,i2,n1,n2,language,memorylayout) GET_IDX2_##language##_##memorylayout(a,i1,i2,n1,n2)

#define IDX2T(a,i1,i2,n1,n2)                             GET_IDX2T(a,i1,i2,n1,n2,FEAT2_PP_AUTO_LANGUAGE(),AUTO_MEMORY_LAYOUT())
#define GET_IDX2T(a,i1,i2,n1,n2,language,memorylayout)   GET_IDX2T_I(a,i1,i2,n1,n2,language,memorylayout)
#define GET_IDX2T_I(a,i1,i2,n1,n2,language,memorylayout) GET_IDX2T_##language##_##memorylayout(a,i1,i2,n1,n2)

#define GET_IDX2_1_1(a,i1,i2,n1,n2) IDX2_FORWARD(a,i1,i2,n1,n2)
#define GET_IDX2_1_2(a,i1,i2,n1,n2) IDX2_REVERSE(a,i1,i2,n1,n2)
#define GET_IDX2_2_1(a,i1,i2,n1,n2) IDX2_REVERSE(a,i1,i2,n1,n2)
#define GET_IDX2_2_2(a,i1,i2,n1,n2) IDX2_FORWARD(a,i1,i2,n1,n2)

#define GET_IDX2T_1_1(a,i1,i2,n1,n2) IDX2T_REVERSE(a,i1,i2,n1,n2)
#define GET_IDX2T_1_2(a,i1,i2,n1,n2) IDX2T_FORWARD(a,i1,i2,n1,n2)
#define GET_IDX2T_2_1(a,i1,i2,n1,n2) IDX2T_FORWARD(a,i1,i2,n1,n2)
#define GET_IDX2T_2_2(a,i1,i2,n1,n2) IDX2T_REVERSE(a,i1,i2,n1,n2)

#if 0
!##############################################################################
! 3D-arrays: IDX3, IDX3T, IDX3_FORWARD, IDX3_REVERSE
!
! All other macros are auxiliary ones and should not be called directly
!##############################################################################
#endif

#define IDX3_FORWARD(a,i1,i2,i3,n1,n2,n3)             IDX3_FORWARD_I(a,i1,i2,i3,n1,n2,n3,FEAT2_PP_AUTO_LANGUAGE())
#define IDX3_FORWARD_I(a,i1,i2,i3,n1,n2,n3,language)  IDX3_FORWARD_II(a,i1,i2,i3,n1,n2,n3,language)
#define IDX3_FORWARD_II(a,i1,i2,i3,n1,n2,n3,language) IDX3_FORWARD_##language(a,i1,i2,i3,n1,n2,n3)
#define IDX3_FORWARD_1(a,i1,i2,i3,n1,n2,n3)           a[IDX3_C(i1,i2,i3,n1,n2,n3)]
#define IDX3_FORWARD_2(a,i1,i2,i3,n1,n2,n3)           a(IDX3_F(i1,i2,i3,n1,n2,n3))

#define IDX3_REVERSE(a,i1,i2,i3,n1,n2,n3)             IDX3_REVERSE_I(a,i1,i2,i3,n1,n2,n3,FEAT2_PP_AUTO_LANGUAGE())
#define IDX3_REVERSE_I(a,i1,i2,i3,n1,n2,n3,language)  IDX3_REVERSE_II(a,i1,i2,i3,n1,n2,n3,language)
#define IDX3_REVERSE_II(a,i1,i2,i3,n1,n2,n3,language) IDX3_REVERSE_##language(a,i1,i2,i3,n1,n2,n3)
#define IDX3_REVERSE_1(a,i1,i2,i3,n1,n2,n3)           a[IDX3_C(i3,i2,i1,n3,n2,n1)]
#define IDX3_REVERSE_2(a,i1,i2,i3,n1,n2,n3)           a(IDX3_F(i3,i2,i1,n3,n2,n1))


#define IDX3(a,i1,i2,i3,n1,n2,n3)                             GET_IDX3(a,i1,i2,i3,n1,n2,n3,FEAT2_PP_AUTO_LANGUAGE(),AUTO_MEMORY_LAYOUT())
#define GET_IDX3(a,i1,i2,i3,n1,n2,n3,language,memorylayout)   GET_IDX3_I(a,i1,i2,i3,n1,n2,n3,language,memorylayout)
#define GET_IDX3_I(a,i1,i2,i3,n1,n2,n3,language,memorylayout) GET_IDX3_##language##_##memorylayout(a,i1,i2,i3,n1,n2,n3)

#define IDX3T(a,i1,i2,i3,n1,n2,n3)                             GET_IDX3T(a,i1,i2,i3,n1,n2,n3,FEAT2_PP_AUTO_LANGUAGE(),AUTO_MEMORY_LAYOUT())
#define GET_IDX3T(a,i1,i2,i3,n1,n2,n3,language,memorylayout)   GET_IDX3T_I(a,i1,i2,i3,n1,n2,n3,language,memorylayout)
#define GET_IDX3T_I(a,i1,i2,i3,n1,n2,n3,language,memorylayout) GET_IDX3T_##language##_##memorylayout(a,i1,i2,i3,n1,n2,n3)

#define GET_IDX3_1_1(a,i1,i2,i3,n1,n2,n3) IDX3_FORWARD(a,i1,i2,i3,n1,n2,n3)
#define GET_IDX3_1_2(a,i1,i2,i3,n1,n2,n3) IDX3_REVERSE(a,i1,i2,i3,n1,n2,n3)
#define GET_IDX3_2_1(a,i1,i2,i3,n1,n2,n3) IDX3_REVERSE(a,i1,i2,i3,n1,n2,n3)
#define GET_IDX3_2_2(a,i1,i2,i3,n1,n2,n3) IDX3_FORWARD(a,i1,i2,i3,n1,n2,n3)

#define GET_IDX3T_1_1(a,i1,i2,i3,n1,n2,n3) IDX3T_REVERSE(a,i1,i2,i3,n1,n2,n3)
#define GET_IDX3T_1_2(a,i1,i2,i3,n1,n2,n3) IDX3T_FORWARD(a,i1,i2,i3,n1,n2,n3)
#define GET_IDX3T_2_1(a,i1,i2,i3,n1,n2,n3) IDX3T_FORWARD(a,i1,i2,i3,n1,n2,n3)
#define GET_IDX3T_2_2(a,i1,i2,i3,n1,n2,n3) IDX3T_REVERSE(a,i1,i2,i3,n1,n2,n3)


#if 0
!##############################################################################
! Generic index wrapper macros for allocating multi-dimensional arrays
! from Fortran and C-code in a unified manner as follows:
!
! allocate(IDX3_ALLOCATE(a,2,3,4))
! malloc(IDX3_ALLOCATE(a,2,3,4))
!
!##############################################################################
#endif

#define IDX1_ALLOCATE(a,n1)       IDX1(a,n1)
#define IDX2_ALLOCATE(a,n1,n2)    IDX2(a,n1,n2,n1,n2)
#define IDX3_ALLOCATE(a,n1,n2,n3) IDX3(a,n1,n2,n3,n1,n2,n3)

#define IDX1_MALLOC(a,n1)         IDX1(a,n1)
#define IDX2_MALLOC(a,n1,n2)      IDX2(a,n1,n2,n1,n2)
#define IDX3_MALLOC(a,n1,n2,n3)   IDX3(a,n1,n2,n3,n1,n2,n3)

#endif
