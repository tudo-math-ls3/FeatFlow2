#if 0
!-*- mode: f90; -*-
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

#include "../../../kernel/System/feat2constants.h"
#include "../../../kernel/System/idxmanager.h"

#if defined(LANGUAGE) && LANGUAGE == LANGUAGE_F

#if 0
!##############################################################################
! Preserve manual line breaks if preprocessing is done by F90CPP script
!##############################################################################
#endif

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

#endif
