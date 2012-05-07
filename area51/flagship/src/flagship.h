#if 0
!##############################################################################
!# ****************************************************************************
!# <name> flagship </name>
!# ****************************************************************************
!#
!# <purpose>
!# This is the main header file for the application
!#
!# </purpose>
!##############################################################################
#endif

#ifndef _FLAGSHIP_H_
#define _FLAGSHIP_H_

#if 0
!##############################################################################
! Default configuration: Fortran programming language
!##############################################################################
#endif

#ifndef LANGUAGE
#define LANGUAGE LANGUAGE_F
#endif

#ifndef IDXADDR
#define IDXADDR LANGUAGE_F
#endif

#include "../../kernel/feat2constants.h"
#include "../../kernel/feat2macros.h"
#include "../../kernel/System/fmath.h"
#include "../../kernel/System/idxmanager.h"

#if 0
!##############################################################################
! Preserve manual line breaks if preprocessing is done by F90CPP script
!##############################################################################
#endif

#if LANGUAGE == LANGUAGE_F
#ifdef USE_PREPROC_F90CPP
#undef MYNEWLINE
#else
#define MYNEWLINE
#endif
#else
#define MYNEWLINE
#endif

#if 0
!##############################################################################
! Use double precision and default integer
!##############################################################################
#endif

#define REAL_PREC DOUBLE_PREC
#undef  INT_PREC

#if 0
!##############################################################################
!# Configuration of global data arrays in host and device memory
!##############################################################################
#endif

#define AOS  0
#define SOA  1

#if 0
! DO NOT EDIT THESE SETTINGS SINCE ALL KERNEL ROUTINES ARE HARD-CODED TO AOS !!!
#endif
#define IDIAGLIST_HOST    AOS
#define IEDGELIST_HOST    AOS
#define COEFFSATDIAG_HOST AOS
#define COEFFSATEDGE_HOST AOS

#if 0
! THESE SETTINGS AFFECT THE DEVICE IMPLEMENTATION ONLY
#endif
#define IDIAGLIST_DEVICE    AOS
#define IEDGELIST_DEVICE    AOS
#define COEFFSATDIAG_DEVICE AOS
#define COEFFSATEDGE_DEVICE AOS

#if 0
!##############################################################################
! Constants for entropy fixes
!##############################################################################
#endif

#define HARTEN_HYMAN_ENTROPYFIX 1001
#define HARTEN_ENTROPYFIX       1002

#if 0
!###############################################################################
! Transport model
!###############################################################################
#endif

#define TRANSP_USE_IBP
#define TRANSP_USE_GFEM_AT_BOUNDARY

#if 0
!###############################################################################
! Hydrodynamic model
!###############################################################################
#endif

#define HYDRO_USE_IBP
#define HYDRO_USE_GFEM_AT_BOUNDARY
#define HYDRO_GAMMA 1.4
#undef  HYDRO_USE_ENTROPYFIX
#define HYDRO_HARTEN_ENTROPYFIX 0.0

#if 0
!###############################################################################
! Magnetohydrodynamic model
!###############################################################################
#endif

#define MHD_USE_IBP
#define MHD_USE_GFEM_AT_BOUNDARY
#define MHD_GAMMA 1.4
#define MHD_VACUUM_PERM 1.0
#define MHD_XMAGFIELD_CONST 0.75
#undef  MHD_USE_ENTROPYFIX
#define MHD_HARTEN_ENTROPYFIX

#endif
