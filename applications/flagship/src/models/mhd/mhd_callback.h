

#ifndef _MHD_CALLBACK_H_
#define _MHD_CALLBACK_H_

#if 0
!##############################################################################
!# ****************************************************************************
!# <name> mhd_callback </name>
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
! Macros: fcb_calcVectorEdgeSys
!         fcb_calcOperatorEdgeSys
!
! User-defined callback functions are supplied to the routines
! gfsys_buildVectorEdge and gfsys_buildOperatorEdge
! if coprocessor support is enabled.
!##############################################################################
#endif

#ifdef ENABLE_COPROCESSOR_SUPPORT
#define COPROC_FCB_CALCVECTOREDGESYS(callbackRoutine)	\
  ,fcb_calcVectorEdgeSys=callbackRoutine
#define COPROC_FCB_CALCOPERATOREDGESYS(callbackRoutine)	\
  ,fcb_calcOperatorEdgeSys=callbackRoutine
#else
#define COPROC_FCB_CALCVECTOREDGESYS(callbackRoutine)
#define COPROC_FCB_CALCOPERATOREDGESYS(callbackRoutine)
#endif

#endif
