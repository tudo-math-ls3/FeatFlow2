#ifndef _TEMPLATE_H_
#define _TEMPLATE_H_

#if 0
!##############################################################################
!# T         : name of the template parameter
!# T_STORAGE : name of the storage identifier
!#             if this constant is undefined, then the implementation is
!#             generated for the derived type T_TYPE defined in a different
!#             module file; for many dynamic data structures the derived type
!#             has to provide routines to compares two items of the same type,
!#             e.g. .eq., .ne., .gt., .lt., .ge., .le. and so on.
!# T_TYPE    : name of the derived type
!#
!#
!# Example : intrinsic type - integer -
!#
!# #define T          Int
!# #define T_STORAGE  ST_INT
!# #define T_TYPE     integer
!#
!# This will instantiate a pseudo-templated data structure for integer data.
!#
!# Example : non-intrinsic type - t_scalarMatrix -
!#
!# #define T          ScMat
!# #undef  T_STORAGE
!# #define T_TYPE     t_scalarMatrix
!#
!# This will instantiate a pseudo-templated data structure for the derived type
!# 't_scalarMatrix' which is originally defined in module 'linearsystemscalar'.
!#
!# Some data structures provide support for auxiliary data, e.g. a linked
!# list is addressed via its key value but each list item has some double
!# data. Thus, each parameter described above has a D_xxx counterpart which
!# specifies the data part in more detail, that is
!#
!# D         : name of the template parameter
!# D_STORAGE : name of the storage identifier
!#             if this constant is undefined, then the implementation is
!#             generated for the derived type D_TYPE defined in a different
!#             module file; for many dynamic data structures the derived type
!#             has to provide routines to compares two items of the same type,
!#             e.g. .eq., .ne., .gt., .lt., .ge., .le. and so on.
!# D_TYPE    : name of the derived type
!#
!##############################################################################
#endif

#ifndef T
#error "Definition of T is missing!"
#endif

#ifndef T_TYPE
#error "Definition of T_TYPE is missing!"
#endif

#if 0
!###############################################################################
!# Macro definition for templating without auxiliary data type
!###############################################################################
#endif

#define FEAT2_PP_TEMPLATE_T_(name,type) name##type
#define FEAT2_PP_TEMPLATE_T(name,type) FEAT2_PP_TEMPLATE_T_(name,type)

#ifdef T_STORAGE
#define FEAT2_PP_TTYPE(t_type) t_type
#else
#define FEAT2_PP_TTYPE(t_type) type(t_type)
#endif


#if 0
!###############################################################################
!# Macro definition for templating with auxiliary data type
!###############################################################################
#endif

#ifdef D
#define FEAT2_PP_TEMPLATE_TD_(name,type,data) name##type##_##data
#define FEAT2_PP_TEMPLATE_TD(name,type,data) FEAT2_PP_TEMPLATE_TD_(name,type,data)
#else
#define FEAT2_PP_TEMPLATE_TD(name,type,data) FEAT2_PP_TEMPLATE_T_(name,type)
#endif

#ifdef D_STORAGE
#define FEAT2_PP_DTYPE(d_type) d_type
#else
#define FEAT2_PP_DTYPE(d_type) type(d_type)
#endif

#endif
