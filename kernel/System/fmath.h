#ifndef _FMATH_H_
#define _FMATH_H_

#include "../feat2constants.h"
#include "../feat2macros.h"

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

#define POW(a,b)              FEAT2_PP_POW(a,b,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_POW(a,b,l)   FEAT2_PP_POW_I(a,b,l)
#define FEAT2_PP_POW_I(a,b,l) FEAT2_PP_POW##_##l(a,b)
#define FEAT2_PP_POW_1(a,b)   pow(a,b)
#define FEAT2_PP_POW_2(a,b)   (a)**(b)

#define MOD(a,b)              FEAT2_PP_MOD(a,b,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_MOD(a,b,l)   FEAT2_PP_MOD_I(a,b,l)
#define FEAT2_PP_MOD_I(a,b,l) FEAT2_PP_MOD##_##l(a,b)
#define FEAT2_PP_MOD_1(a,b)   a%b
#define FEAT2_PP_MOD_2(a,b)   mod(a,b)

#define EQ(a,b)              FEAT2_PP_EQ(a,b,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_EQ(a,b,l)   FEAT2_PP_EQ_I(a,b,l)
#define FEAT2_PP_EQ_I(a,b,l) FEAT2_PP_EQ##_##l(a,b)
#define FEAT2_PP_EQ_1(a,b)   (a)==(b)
#define FEAT2_PP_EQ_2(a,b)   (a).eq.(b)

#define NE(a,b)              FEAT2_PP_NE(a,b,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_NE(a,b,l)   FEAT2_PP_NE_I(a,b,l)
#define FEAT2_PP_NE_I(a,b,l) FEAT2_PP_NE##_##l(a,b)
#define FEAT2_PP_NE_1(a,b)   (a)!=(b)
#define FEAT2_PP_NE_2(a,b)   (a).ne.(b)

#define LT(a,b)              FEAT2_PP_LT(a,b,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_LT(a,b,l)   FEAT2_PP_LT_I(a,b,l)
#define FEAT2_PP_LT_I(a,b,l) FEAT2_PP_LT##_##l(a,b)
#define FEAT2_PP_LT_1(a,b)   (a)<(b)
#define FEAT2_PP_LT_2(a,b)   (a).lt.(b)

#define GT(a,b)              FEAT2_PP_GT(a,b,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_GT(a,b,l)   FEAT2_PP_GT_I(a,b,l)
#define FEAT2_PP_GT_I(a,b,l) FEAT2_PP_GT##_##l(a,b)
#define FEAT2_PP_GT_1(a,b)   (a)>(b)
#define FEAT2_PP_GT_2(a,b)   (a).gt.(b)

#define LE(a,b)              FEAT2_PP_LE(a,b,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_LE(a,b,l)   FEAT2_PP_LE_I(a,b,l)
#define FEAT2_PP_LE_I(a,b,l) FEAT2_PP_LE##_##l(a,b)
#define FEAT2_PP_LE_1(a,b)   (a)<=(b)
#define FEAT2_PP_LE_2(a,b)   (a).le.(b)

#define GE(a,b)              FEAT2_PP_GE(a,b,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_GE(a,b,l)   FEAT2_PP_GE_I(a,b,l)
#define FEAT2_PP_GE_I(a,b,l) FEAT2_PP_GE##_##l(a,b)
#define FEAT2_PP_GE_1(a,b)   (a)>=(b)
#define FEAT2_PP_GE_2(a,b)   (a).ge.(b)

#define NOT(a)               FEAT2_PP_NOT(a,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_NOT(a,l)    FEAT2_PP_NOT_I(a,l)
#define FEAT2_PP_NOT_I(a,l)  FEAT2_PP_NOT##_##l(a)
#define FEAT2_PP_NOT_1(a)    !(a)
#define FEAT2_PP_NOT_2(a)    .not.(a)

#define AND(a,b)              FEAT2_PP_AND(a,b,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_AND(a,b,l)   FEAT2_PP_AND_I(a,b,l)
#define FEAT2_PP_AND_I(a,b,l) FEAT2_PP_AND##_##l(a,b)
#define FEAT2_PP_AND_1(a,b)   (a)&&(b)
#define FEAT2_PP_AND_2(a,b)   (a).and.(b)

#define OR(a,b)              FEAT2_PP_OR(a,b,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_OR(a,b,l)   FEAT2_PP_OR_I(a,b,l)
#define FEAT2_PP_OR_I(a,b,l) FEAT2_PP_OR##_##l(a,b)
#define FEAT2_PP_OR_1(a,b)   (a)||(b)
#define FEAT2_PP_OR_2(a,b)   (a).or.(b)

#define IAND(a,b)              FEAT2_PP_IAND(a,b,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_IAND(a,b,l)   FEAT2_PP_IAND_I(a,b,l)
#define FEAT2_PP_IAND_I(a,b,l) FEAT2_PP_IAND##_##l(a,b)
#define FEAT2_PP_IAND_1(a,b)   (a)&(b)
#define FEAT2_PP_IAND_2(a,b)   iand(a,b)

#define IOR(a,b)              FEAT2_PP_IOR(a,b,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_IOR(a,b,l)   FEAT2_PP_IOR_I(a,b,l)
#define FEAT2_PP_IOR_I(a,b,l) FEAT2_PP_IOR##_##l(a,b)
#define FEAT2_PP_IOR_1(a,b)   (a)|(b)
#define FEAT2_PP_IOR_2(a,b)   ior(a,b)

#define IEOR(a,b)              FEAT2_PP_IEOR(a,b,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_IEOR(a,b,l)   FEAT2_PP_IEOR_I(a,b,l)
#define FEAT2_PP_IEOR_I(a,b,l) FEAT2_PP_IEOR##_##l(a,b)
#define FEAT2_PP_IEOR_1(a,b)   (a)^(b)
#define FEAT2_PP_IEOR_2(a,b)   ieor(a,b)

#define INOT(a)               FEAT2_PP_INOT(a,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_INOT(a,l)    FEAT2_PP_INOT_I(a,l)
#define FEAT2_PP_INOT_I(a,l)  FEAT2_PP_INOT##_##l(a)
#define FEAT2_PP_INOT_1(a)    ~(a)
#define FEAT2_PP_INOT_2(a)    inot(a)

#define ISHFT(a,b)              FEAT2_PP_ISHFT(a,b,FEAT2_PP_AUTO_LANGUAGE())
#define FEAT2_PP_ISHFT(a,b,l)   FEAT2_PP_ISHFT_I(a,b,l)
#define FEAT2_PP_ISHFT_I(a,b,l) FEAT2_PP_ISHFT##_##l(a,b)
#define FEAT2_PP_ISHFT_1(a,b)   (b>0 ? a<<b : a>>b)
#define FEAT2_PP_ISHFT_2(a,b)   ishft(a,b)



#if 0
!##############################################################################
! Mathematical constants
!##############################################################################
#endif

#define CONST(a,prec) FEAT2_PP_CONST(a,prec)
#define QCONST(a)     FEAT2_PP_CONST(a,QUAD_PREC)
#define DCONST(a)     FEAT2_PP_CONST(a,DOUBLE_PREC)
#define SCONST(a)     FEAT2_PP_CONST(a,SINGLE_PREC)
#define ICONST(a)     FEAT2_PP_CONST(a,INTEGER_PREC)

#endif
