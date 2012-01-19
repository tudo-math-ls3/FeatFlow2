#ifndef _TEMPLATE_H_
#define _TEMPLATE_H_

#if 0
!##############################################################################
!# T         : name of the template parameter
!# T_STORAGE : name of the storage identifier
!#             if this constant is undefined, then the implementation
!#             is generated for the derived type T_NAME which must be
!#             defined in the module T_MODULE
!# T_TYPE    : name of the derived type
!# T_MODULE  : name of the module where the derived type is defined
!#
!#
!# Example : intrinsic type - integer -
!#
!# #define T          Int
!# #define T_STORAGE  ST_INT
!# #define T_TYPE     integer
!# #undef  T_MODULE
!#
!# This will instantiate a pseudo-templated data structure for integer data.
!#
!# Example : non-intrinsic type - t_scalarMatrix -
!#
!# #define T          ScMat
!# #undef  T_STORAGE
!# #define T_TYPE     t_scalarMatrix
!# #define T_MODULE   linearsystemscalar
!#
!# This will instantiate a pseudo-templated data structure for the derived type
!# 't_scalarMatrix' which is originally defined in module 'linearsystemscalar'.
!# In order to make sure that the configure scripts is able to create the
!# dependency list correctly, the following lines have to be added:
!#
!# #if 0
!# use linearsystemscalar
!# #endif
!#
!# This ensures that the configure script considers module 'linearsystemscalar'
!# as prerequisite but the use-statement is neglected by the compiler.
!#
!# Some data structures provide support for auxiliary data, e.g. a linked
!# list is addressed via its key value but each list item has some double
!# data. Thus, each parameter described above has a D_xxx counterpart which
!# specifies the data part in more detail.
!#
!##############################################################################
#endif

#ifndef T
#error "Definition of T is missing!"
#endif

#ifndef T_TYPE
#error "Definition of T_TYPE is missing!"
#endif

#ifndef T_STORAGE
#ifndef T_MODULE
#error "Definition of T_STORAGE and T_MODULE is missing!"
#endif
#endif

#if 0
!# Macro definition for templating without data type
#endif

#define TEMPLATE_(name,type) name##type
#define TEMPLATE(name,type) TEMPLATE_(name,type)

#ifdef T_STORAGE
#define TTYPE(t_type) t_type
#else
#define TTYPE(t_type) type(t_type)
#endif


#if 0
!# Macro definition for templating with data type
#endif

#ifdef D
#define TEMPLATE_D_(name,type,data) name##type##_##data
#define TEMPLATE_D(name,type,data) TEMPLATE_D_(name,type,data)
#else
#define TEMPLATE_D(name,type,data) TEMPLATE_(name,type)
#endif

#ifdef D_STORAGE
#define DTYPE(d_type) d_type
#else
#define DTYPE(d_type) type(d_type)
#endif


#if 0
!# Macro definition for hiding module files from configure script
#endif

#define __external_use__(module) use module


#endif
