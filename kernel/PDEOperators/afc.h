#include "kernel/feat2constants.h" 

#define LANGUAGE LANGUAGE_F

#if 0
! ------------------------------------------------------------------------------
! Convention: all internal identifiers have two leading and two
! trailing underscores and are undefined at the end of this file
! ------------------------------------------------------------------------------
#endif

#ifdef __AFCName__ 
#undef __AFCName__
#endif

#ifdef __AFCData__
#undef __AFCData__
#endif

#ifdef __AFCType__
#undef __AFCType__
#endif

#if TemplateType_AFC == SINGLE_PREC
#define __AFCName__ SP
#define __AFCData__ real(SP)
#define __AFCType__ SP
#elif TemplateType_AFC == DOUBLE_PREC
#define __AFCName__ DP
#define __AFCData__ real(DP)
#define __AFCType__ DP
#elif TemplateType_AFC == QUAD_PREC
#define __AFCName__ QP
#define __AFCData__ real(QP)
#define __AFCType__ QP
#endif


#if 0
! ------------------------------------------------------------------------------
#endif

#ifdef __VectorName__ 
#undef __VectorName__
#endif

#ifdef __VectorData__
#undef __VectorData__
#endif

#ifdef __VectorType__
#undef __VectorType__
#endif

#if TemplateType_Vector == SINGLE_PREC
#define __VectorName__ SP
#define __VectorData__ real(SP)
#define __VectorType__ SP
#elif TemplateType_Vector == DOUBLE_PREC
#define __VectorName__ DP
#define __VectorData__ real(DP)
#define __VectorType__ DP
#elif TemplateType_Vector == QUAD_PREC
#define __VectorName__ QP
#define __VectorData__ real(QP)
#define __VectorType__ QP
#endif


#if 0
! ------------------------------------------------------------------------------
#endif

#ifdef __MatrixName__ 
#undef __MatrixName__
#endif

#ifdef __MatrixData__
#undef __MatrixData__
#endif

#ifdef __MatrixType__
#undef __MatrixType__
#endif

#if TemplateType_Matrix == SINGLE_PREC
#define __MatrixName__ SP
#define __MatrixData__ real(SP)
#define __MatrixType__ SP
#elif TemplateType_Matrix == DOUBLE_PREC
#define __MatrixName__ DP
#define __MatrixData__ real(DP)
#define __MatrixType__ DP
#elif TemplateType_Matrix == QUAD_PREC
#define __MatrixName__ QP
#define __MatrixData__ real(QP)
#define __MatrixType__ QP
#endif
