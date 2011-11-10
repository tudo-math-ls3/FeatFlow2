#ifndef _HYDRO_CALLBACK_H_
#define _HYDRO_CALLBACK_H_

#if 0
!-*- mode: f90; -*-
!##############################################################################
!# ****************************************************************************
!# <name> hydro_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This is the header file for callback routines
!#
!# </purpose>
!##############################################################################
#endif

#if 0
!##############################################################################
! Macro: fcb_calcVectorEdgeSys
!
! User-defined callback function is supplied to routine gfsys_buildVectorEdge
! if coprocessor support is enabled.
!##############################################################################
#endif

#ifdef ENABLE_COPROCESSOR_SUPPORT
#define __fcb_calcVectorEdgeSys__(callbackRoutine)	\
  ,fcb_calcVectorEdgeSys=callbackRoutine
#else
#define __fcb_calcVectorEdgeSys__(callbackRoutine)
#endif

#endif
