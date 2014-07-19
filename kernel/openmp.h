#ifndef _OPENMP_H_
#define _OPENMP_H_

#if 0
!##############################################################################
!# <name> openmp </name>
!# ****************************************************************************
!#
!# <purpose>
!# This header file defines a couple of preprocessor macros to simplify
!# enabling/disabling special OpenMP features which are only available
!# in newer versions of the OpenMP standard
!# </purpose>
!##############################################################################
#endif

#ifdef HAS_OPENMP20
#define omp20(yes,no) yes
#else
#define omp20(yes,no) no
#endif

#ifdef HAS_OPENMP25
#define omp25(yes,no) yes
#else
#define omp25(yes,no) no
#endif

#ifdef HAS_OPENMP30
#define omp30(yes,no) yes
#else
#define omp30(yes,no) no
#endif

#ifdef HAS_OPENMP31
#define omp31(yes,no) yes
#else
#define omp31(yes,no) no
#endif

#ifdef HAS_OPENMP40
#define omp40(yes,no) yes
#else
#define omp40(yes,no) no
#endif

#define omp(spec,yes,no) omp ## spec(yes,no)

#endif
