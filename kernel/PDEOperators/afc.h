#include "../feat2constants.h" 

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
