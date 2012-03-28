#ifndef _FMATH_H_
#define _FMATH_H_

#include "../feat2constants.h"

#if 0
!##############################################################################
!# <name> fmath </name>
!# ****************************************************************************
!#
!# <purpose>
!# This header file provides preprocessor constants and macros for common
!# mathematical functions which have different function names in Fortran an C.
!#
!# </purpose>
!##############################################################################
#endif


#if 0
!##############################################################################
! Mathematical functions
!##############################################################################
#endif

#if LANGUAGE == LANGUAGE_F

#define POW(a,b)   (a)**(b)
#define MOD(a,b)   mod(a,b)

#define EQ(a,b)    (a).eq.(b)
#define NE(a,b)    (a).ne.(b)
#define LT(a,b)    (a).lt.(b)
#define GT(a,b)    (a).gt.(b)
#define LE(a,b)    (a).le.(b)
#define GE(a,b)    (a).ge.(b)
#define NOT(a)     .not.(a)
#define AND(a,b)   (a).and.(b)
#define OR(a,b)    (a).or.(b)

#define IAND(a,b)  iand(a,b)
#define IOR(a,b)   ior(a,b)
#define IEOR(a,b)  ieor(a,b)
#define INOT(a)    inot(a)
#define ISHFT(a,b) ishft(a,b)


#elif LANGUAGE == LANGUAGE_C

#define POW(a,b)   pow(a,b)
#define MOD(a,b)   a%b

#define EQ(a,b)    (a)==(b)
#define NE(a,b)    (a)!=(b)
#define LT(a,b)    (a)<(b)
#define GT(a,b)    (a)>(b)
#define LE(a,b)    (a)<=(b)
#define GE(a,b)    (a)>=(b)
#define NOT(a)     !(a)
#define AND(a,b)   (a)&&(b)
#define OR(a,b)    (a)||(b)

#define IAND(a,b)  (a)&(b)
#define IOR(a,b)   (a)|(b)
#define IEOR(a,b)  (a)^(b)
#define INOT(a)    ~(a)
#define ISHFT(a,b) (b>0 ? a<<b : a>>b)


#else
#error "Unsupported programming language!"
#endif


#if 0
!##############################################################################
! Mathematical constants
!##############################################################################
#endif

#if LANGUAGE == LANGUAGE_F

#define CONST(a,prec) a##_##prec

#elif LANGUAGE == LANGUAGE_C

#define CONST(a,prec) a

#else
#error "Unsupported programming language!"
#endif


#ifdef REAL_PREC

#if REAL_PREC == QUAD_PREC
#define RCONST(a) CONST(a,QP)
#elif REAL_PREC == DOUBLE_PREC
#define RCONST(a) CONST(a,DP)
#elif REAL_PREC == SINGLE_PREC
#define RCONST(a) CONST(a,SP)
#else
#error "Unsupported real kind precision!"
#endif

#else
#define RCONST(a) a
#endif

#ifdef INT_PREC

#if INT_PREC == INT8_PREC
#define ICONST(a) CONST(a,I8)
#elif INT_PREC == INT16_PREC
#define ICONST(a) CONST(a,I16)
#elif INT_PREC == INT32_PREC
#define ICONST(a) CONST(a,I32)
#elif INT_PREC == INT64_PREC
#define ICONST(a) CONST(a,I64)
#else
#error "Unsupported integer kind precision!"
#endif

#else
#define ICONST(a) a
#endif

#endif
