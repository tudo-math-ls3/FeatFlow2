#ifndef _THERMODYNAMICS_H_
#define _THERMODYNAMICS_H_

#if 0
!-*- mode: f90; -*-
!##############################################################################
!# ****************************************************************************
!# <name> thermodynamic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This is the thermodynamic header file
!#
!# </purpose>
!##############################################################################
#endif


#include "../flagship.h"


#if 0
!##############################################################################
! Compile-time constants
!##############################################################################
#endif

#if 0
! Speficy ratio of specific heats (if required)
#endif
#ifndef THERMODYN_GAMMA
#define THERMODYN_GAMMA 1.4_DP
#endif


#if 0
!##############################################################################
! Velocity field
!##############################################################################
#endif

#if 0
! Compute x-velocity from conservative variables in 1D, 2D, and 3D
#endif
#define XVELOCITY1(U,IDX)              (XMOMENTUM1(U,IDX)/DENSITY1(U,IDX))
#define XVELOCITY2(U,IDX,i,n1,n2)      (XMOMENTUM2(U,IDX,i,n1,n2)/DENSITY2(U,IDX,i,n1,n2))
#define XVELOCITY3(U,IDX,i,j,n1,n2,n3) (XMOMENTUM3(U,IDX,i,j,n1,n2,n3)/DENSITY3(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute y-velocity from conservative variables in 2D, and 3D
#endif
#define YVELOCITY1(U,IDX)              (YMOMENTUM1(U,IDX)/DENSITY1(U,IDX))
#define YVELOCITY2(U,IDX,i,n1,n2)      (YMOMENTUM2(U,IDX,i,n1,n2)/DENSITY2(U,IDX,i,n1,n2))
#define YVELOCITY3(U,IDX,i,j,n1,n2,n3) (YMOMENTUM3(U,IDX,i,j,n1,n2,n3)/DENSITY3(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute z-velocity from conservative variables in 3D
#endif
#define ZVELOCITY1(U,IDX)              (ZMOMENTUM1(U,IDX)/DENSITY1(U,IDX))
#define ZVELOCITY2(U,IDX,i,n1,n2)      (ZMOMENTUM2(U,IDX,i,n1,n2)/DENSITY2(U,IDX,i,n1,n2))
#define ZVELOCITY3(U,IDX,i,j,n1,n2,n3) (ZMOMENTUM3(U,IDX,i,j,n1,n2,n3)/DENSITY3(U,IDX,i,j,n1,n2,n3))


#if 0
!##############################################################################
! Magnitude of velocity field
!##############################################################################
#endif

#if 0
! Compute magnitude of velocity field in 1D
#endif
#define VELMAGNITUDE1_1D(U,IDX)              (abs(XVELOCITY1(U,IDX)))
#define VELMAGNITUDE2_1D(U,IDX,i,n1,n2)      (abs(XVELOCITY2(U,IDX,i,n1,n2)))
#define VELMAGNITUDE3_1D(U,IDX,i,j,n1,n2,n3) (abs(XVELOCITY3(U,IDX,i,j,n1,n2,n3)))

#if 0
! Compute magnitude of velocity field in 2D
#endif
#define VELMAGNITUDE1_2D(U,IDX)              (sqrt(XVELOCITY1(U,IDX)**2+MYNEWLINE YVELOCITY1(U,IDX)**2))
#define VELMAGNITUDE2_2D(U,IDX,i,n1,n2)      (sqrt(XVELOCITY2(U,IDX,i,n1,n2)**2+MYNEWLINE YVELOCITY2(U,IDX,i,n1,n2)**2))
#define VELMAGNITUDE3_2D(U,IDX,i,j,n1,n2,n3) (sqrt(XVELOCITY3(U,IDX,i,j,n1,n2,n3)**2+MYNEWLINE YVELOCITY3(U,IDX,i,j,n1,n2,n3)**2))

#if 0
! Compute magnitude of velocity field in 3D
#endif
#define VELMAGNITUDE1_3D(U,IDX)              (sqrt(XVELOCITY1(U,IDX)**2+MYNEWLINE YVELOCITY1(U,IDX)**2+MYNEWLINE ZVELOCITY1(U,IDX)**2))
#define VELMAGNITUDE2_3D(U,IDX,i,n1,n2)      (sqrt(XVELOCITY2(U,IDX,i,n1,n2)**2+MYNEWLINE YVELOCITY2(U,IDX,i,n1,n2)**2+MYNEWLINE ZVELOCITY2(U,IDX,i,n1,n2)**2))
#define VELMAGNITUDE3_3D(U,IDX,i,j,n1,n2,n3) (sqrt(XVELOCITY3(U,IDX,i,j,n1,n2,n3)**2+MYNEWLINE YVELOCITY3(U,IDX,i,j,n1,n2,n3)**2+MYNEWLINE ZVELOCITY3(U,IDX,i,j,n1,n2,n3)**2))

#if THERMODYN_NDIM == 1
#define VELMAGNITUDE1 VELMAGNITUDE1_1D
#define VELMAGNITUDE2 VELMAGNITUDE2_1D
#define VELMAGNITUDE3 VELMAGNITUDE3_1D
#elif THERMODYN_NDIM == 2
#define VELMAGNITUDE1 VELMAGNITUDE1_2D
#define VELMAGNITUDE2 VELMAGNITUDE2_2D
#define VELMAGNITUDE3 VELMAGNITUDE3_2D
#elif THERMODYN_NDIM == 3
#define VELMAGNITUDE1 VELMAGNITUDE1_3D
#define VELMAGNITUDE2 VELMAGNITUDE2_3D
#define VELMAGNITUDE3 VELMAGNITUDE3_3D
#else
#ifdef THERMODYN_NDIM
#error "Invalid spatial dimension specified!"
#endif
#endif


#if 0
!##############################################################################
! Pressure
!##############################################################################
#endif

#if 0
! Compute pressure from conservative variables in 1D
#endif
#define PRESSURE1_1D(U,IDX)              (((THERMODYN_GAMMA)-1.0_DP)*INTERNALENERGY1_1D(U,IDX))
#define PRESSURE2_1D(U,IDX,i,n1,n2)      (((THERMODYN_GAMMA)-1.0_DP)*INTERNALENERGY2_1D(U,IDX,i,n1,n2))
#define PRESSURE3_1D(U,IDX,i,j,n1,n2,n3) (((THERMODYN_GAMMA)-1.0_DP)*INTERNALENERGY3_1D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute pressure from conservative variables in 2D
#endif
#define PRESSURE1_2D(U,IDX)              (((THERMODYN_GAMMA)-1.0_DP)*INTERNALENERGY1_2D(U,IDX))
#define PRESSURE2_2D(U,IDX,i,n1,n2)      (((THERMODYN_GAMMA)-1.0_DP)*INTERNALENERGY2_2D(U,IDX,i,n1,n2))
#define PRESSURE3_2D(U,IDX,i,j,n1,n2,n3) (((THERMODYN_GAMMA)-1.0_DP)*INTERNALENERGY3_2D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute pressure from conservative variables in 3D
#endif
#define PRESSURE1_3D(U,IDX)              (((THERMODYN_GAMMA)-1.0_DP)*INTERNALENERGY1_3D(U,IDX))
#define PRESSURE2_3D(U,IDX,i,n1,n2)      (((THERMODYN_GAMMA)-1.0_DP)*INTERNALENERGY2_3D(U,IDX,i,n1,n2))
#define PRESSURE3_3D(U,IDX,i,j,n1,n2,n3) (((THERMODYN_GAMMA)-1.0_DP)*INTERNALENERGY3_3D(U,IDX,i,j,n1,n2,n3))


#if THERMODYN_NDIM == 1
#define PRESSURE1 PRESSURE1_1D
#define PRESSURE2 PRESSURE2_1D
#define PRESSURE3 PRESSURE3_1D
#elif THERMODYN_NDIM == 2
#define PRESSURE1 PRESSURE1_2D
#define PRESSURE2 PRESSURE2_2D
#define PRESSURE3 PRESSURE3_2D
#elif THERMODYN_NDIM == 3
#define PRESSURE1 PRESSURE1_3D
#define PRESSURE2 PRESSURE2_3D
#define PRESSURE3 PRESSURE3_3D
#else
#ifdef THERMODYN_NDIM
#error "Invalid spatial dimension specified!"
#endif
#endif


#if 0
!##############################################################################
! Internal Energy
!##############################################################################
#endif

#if 0
! Compute internal energy from conservative variables in 1D
#endif
#define INTERNALENERGY1_1D(U,IDX)              (TOTALENERGY1(U,IDX)-MYNEWLINE KINETICENERGY1_1D(U,IDX))
#define INTERNALENERGY2_1D(U,IDX,i,n1,n2)      (TOTALENERGY2(U,IDX,i,n1,n2)-MYNEWLINE KINETICENERGY2_1D(U,IDX,i,n1,n2))
#define INTERNALENERGY3_1D(U,IDX,i,j,n1,n2,n3) (TOTALENERGY3(U,IDX,i,j,n1,n2,n3)-MYNEWLINE KINETICENERGY3_1D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute internal energy from conservative variables in 2D
#endif
#define INTERNALENERGY1_2D(U,IDX)              (TOTALENERGY1(U,IDX)-MYNEWLINE KINETICENERGY1_2D(U,IDX))
#define INTERNALENERGY2_2D(U,IDX,i,n1,n2)      (TOTALENERGY2(U,IDX,i,n1,n2)-MYNEWLINE KINETICENERGY2_2D(U,IDX,i,n1,n2))
#define INTERNALENERGY3_2D(U,IDX,i,j,n1,n2,n3) (TOTALENERGY3(U,IDX,i,j,n1,n2,n3)-MYNEWLINE KINETICENERGY3_2D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute internal energy from conservative variables in 3D
#endif
#define INTERNALENERGY1_3D(U,IDX)              (TOTALENERGY1(U,IDX)-MYNEWLINE KINETICENERGY1_3D(U,IDX))
#define INTERNALENERGY2_3D(U,IDX,i,n1,n2)      (TOTALENERGY2(U,IDX,i,n1,n2)-MYNEWLINE KINETICENERGY2_3D(U,IDX,i,n1,n2))
#define INTERNALENERGY3_3D(U,IDX,i,j,n1,n2,n3) (TOTALENERGY3(U,IDX,i,j,n1,n2,n3)-MYNEWLINE KINETICENERGY3_3D(U,IDX,i,j,n1,n2,n3))

#if THERMODYN_NDIM == 1
#define INTERNALENERGY1 INTERNALENERGY1_1D
#define INTERNALENERGY2 INTERNALENERGY2_1D
#define INTERNALENERGY3 INTERNALENERGY3_1D
#elif THERMODYN_NDIM == 2
#define INTERNALENERGY1 INTERNALENERGY1_2D
#define INTERNALENERGY2 INTERNALENERGY2_2D
#define INTERNALENERGY3 INTERNALENERGY3_2D
#elif THERMODYN_NDIM == 3
#define INTERNALENERGY1 INTERNALENERGY1_3D
#define INTERNALENERGY2 INTERNALENERGY2_3D
#define INTERNALENERGY3 INTERNALENERGY3_3D
#else
#ifdef THERMODYN_NDIM
#error "Invalid spatial dimension specified!"
#endif
#endif


#if 0
!##############################################################################
! Kinetic Energy
!##############################################################################
#endif

#if 0
! Compute kinetic energy from conservative variables in 1D
#endif
#define KINETICENERGY1_1D(U,IDX)              (0.5_DP*(XMOMENTUM1(U,IDX)**2)/MYNEWLINE DENSITY1(U,IDX))
#define KINETICENERGY2_1D(U,IDX,i,n1,n2)      (0.5_DP*(XMOMENTUM2(U,IDX,i,n1,n2)**2)/MYNEWLINE DENSITY2(U,IDX,i,n1,n2))
#define KINETICENERGY3_1D(U,IDX,i,j,n1,n2,n3) (0.5_DP*(XMOMENTUM3(U,IDX,i,j,n1,n2,n3)**2)/MYNEWLINE DENSITY3(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute kinetic energy from conservative variables in 2D
#endif
#define KINETICENERGY1_2D(U,IDX)              (0.5_DP*(XMOMENTUM1(U,IDX)**2+MYNEWLINE YMOMENTUM1(U,IDX)**2)/MYNEWLINE DENSITY1(U,IDX))
#define KINETICENERGY2_2D(U,IDX,i,n1,n2)      (0.5_DP*(XMOMENTUM2(U,IDX,i,n1,n2)**2+MYNEWLINE YMOMENTUM2(U,IDX,i,n1,n2)**2)/MYNEWLINE DENSITY2(U,IDX,i,n1,n2))
#define KINETICENERGY3_2D(U,IDX,i,j,n1,n2,n3) (0.5_DP*(XMOMENTUM3(U,IDX,i,j,n1,n2,n3)**2+MYNEWLINE YMOMENTUM3(U,IDX,i,j,n1,n2,n3)**2)/MYNEWLINE DENSITY3(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute kinetic energy from conservative variables in 3D
#endif
#define KINETICENERGY1_3D(U,IDX)              (0.5_DP*(XMOMENTUM1(U,IDX)**2+MYNEWLINE YMOMENTUM1(U,IDX)**2+MYNEWLINE ZMOMENTUM1(U,IDX)**2)/MYNEWLINE DENSITY1(U,IDX))
#define KINETICENERGY2_3D(U,IDX,i,n1,n2)      (0.5_DP*(XMOMENTUM2(U,IDX,i,n1,n2)**2+MYNEWLINE YMOMENTUM2(U,IDX,i,n1,n2)**2+MYNEWLINE ZMOMENTUM2(U,IDX,i,n1,n2)**2)/MYNEWLINE DENSITY2(U,IDX,i,n1,n2))
#define KINETICENERGY3_3D(U,IDX,i,j,n1,n2,n3) (0.5_DP*(XMOMENTUM3(U,IDX,i,j,n1,n2,n3)**2+MYNEWLINE YMOMENTUM3(U,IDX,i,j,n1,n2,n3)**2+MYNEWLINE ZMOMENTUM3(U,IDX,i,j,n1,n2,n3)**2)/MYNEWLINE DENSITY3(U,IDX,i,j,n1,n2,n3))

#if THERMODYN_NDIM == 1
#define KINETICENERGY1 KINETICENERGY1_1D
#define KINETICENERGY2 KINETICENERGY2_1D
#define KINETICENERGY3 KINETICENERGY3_1D
#elif THERMODYN_NDIM == 2
#define KINETICENERGY1 KINETICENERGY1_2D
#define KINETICENERGY2 KINETICENERGY2_2D
#define KINETICENERGY3 KINETICENERGY3_2D
#elif THERMODYN_NDIM == 3
#define KINETICENERGY1 KINETICENERGY1_3D
#define KINETICENERGY2 KINETICENERGY2_3D
#define KINETICENERGY3 KINETICENERGY3_3D
#else
#ifdef THERMODYN_NDIM
#error "Invalid spatial dimension specified!"
#endif
#endif


#if 0
!##############################################################################
! Specific total energy, aka, total energy per unit mass
!##############################################################################
#endif

#define SPECIFICTOTALENERGY1(U,IDX)              (TOTALENERGY1(U,IDX)/DENSITY1(U,IDX))
#define SPECIFICTOTALENERGY2(U,IDX,i,n1,n2)      (TOTALENERGY2(U,IDX,i,n1,n2)/DENSITY2(U,IDX,i,n1,n2))
#define SPECIFICTOTALENERGY3(U,IDX,i,j,n1,n2,n3) (TOTALENERGY3(U,IDX,i,j,n1,n2,n3)/DENSITY3(U,IDX,i,j,n1,n2,n3))


#if 0
!##############################################################################
! Speed of sound
!##############################################################################
#endif

#if 0
! Compute speed of sound for a thermally perfect gas in 1D
#endif
#define SOUNDSPEED1_1D(U,IDX)              (sqrt((THERMODYN_GAMMA)*MYNEWLINE PRESSURE1_1D(U,IDX)/MYNEWLINE DENSITY1(U,IDX)))
#define SOUNDSPEED2_1D(U,IDX,i,n1,n2)      (sqrt((THERMODYN_GAMMA)*MYNEWLINE PRESSURE2_1D(U,IDX,i,n1,n2)/MYNEWLINE DENSITY2(U,IDX,i,n1,n2)))
#define SOUNDSPEED3_1D(U,IDX,i,j,n1,n2,n3) (sqrt((THERMODYN_GAMMA)*MYNEWLINE PRESSURE3_1D(U,IDX,i,j,n1,n2,n3)/MYNEWLINE DENSITY3(U,IDX,i,j,n1,n2,n3)))

#if 0
! Compute speed of sound for a thermally perfect gas in 2D
#endif
#define SOUNDSPEED1_2D(U,IDX)              (sqrt((THERMODYN_GAMMA)*MYNEWLINE PRESSURE1_2D(U,IDX)/MYNEWLINE DENSITY1(U,IDX)))
#define SOUNDSPEED2_2D(U,IDX,i,n1,n2)      (sqrt((THERMODYN_GAMMA)*MYNEWLINE PRESSURE2_2D(U,IDX,i,n1,n2)/MYNEWLINE DENSITY2(U,IDX,i,n1,n2)))
#define SOUNDSPEED3_2D(U,IDX,i,j,n1,n2,n3) (sqrt((THERMODYN_GAMMA)*MYNEWLINE PRESSURE3_2D(U,IDX,i,j,n1,n2,n3)/MYNEWLINE DENSITY3(U,IDX,i,j,n1,n2,n3)))

#if 0
! Compute speed of sound for a thermally perfect gas in 3D
#endif
#define SOUNDSPEED1_3D(U,IDX)              (sqrt((THERMODYN_GAMMA)*MYNEWLINE PRESSURE1_3D(U,IDX)/MYNEWLINE DENSITY1(U,IDX)))
#define SOUNDSPEED2_3D(U,IDX,i,n1,n2)      (sqrt((THERMODYN_GAMMA)*MYNEWLINE PRESSURE2_3D(U,IDX,i,n1,n2)/MYNEWLINE DENSITY2(U,IDX,i,n1,n2)))
#define SOUNDSPEED3_3D(U,IDX,i,j,n1,n2,n3) (sqrt((THERMODYN_GAMMA)*MYNEWLINE PRESSURE3_3D(U,IDX,i,j,n1,n2,n3)/MYNEWLINE DENSITY3(U,IDX,i,j,n1,n2,n3)))

#if THERMODYN_NDIM == 1
#define SOUNDSPEED1 SOUNDSPEED1_1D
#define SOUNDSPEED2 SOUNDSPEED2_1D
#define SOUNDSPEED3 SOUNDSPEED3_1D
#elif THERMODYN_NDIM == 2
#define SOUNDSPEED1 SOUNDSPEED1_2D
#define SOUNDSPEED2 SOUNDSPEED2_2D
#define SOUNDSPEED3 SOUNDSPEED3_2D
#elif THERMODYN_NDIM == 3
#define SOUNDSPEED1 SOUNDSPEED1_3D
#define SOUNDSPEED2 SOUNDSPEED2_3D
#define SOUNDSPEED3 SOUNDSPEED3_3D
#else
#ifdef THERMODYN_NDIM
#error "Invalid spatial dimension specified!"
#endif
#endif


#if 0
!##############################################################################
! Mach number
!##############################################################################
#endif

#if 0
! Compute the Mach number in 1D
#endif
#define MACHNUMBER1_1D(U,IDX)              (abs(XVELOCITY1(U,IDX))/MYNEWLINE SOUNDSPEED1_1D(U,IDX))
#define MACHNUMBER2_1D(U,IDX,i,n1,n2)      (abs(XVELOCITY2(U,IDX,i,n1,n2))/MYNEWLINE SOUNDSPEED2_1D(U,IDX,i,n1,n2))
#define MACHNUMBER3_1D(U,IDX,i,j,n1,n2,n3) (abs(XVELOCITY3(U,IDX,i,j,n1,n2,n3))/MYNEWLINE SOUNDSPEED3_1D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute the Mach number in 2D
#endif
#define MACHNUMBER1_2D(U,IDX)              (sqrt(XVELOCITY1(U,IDX)**2+MYNEWLINE YVELOCITY1(U,IDX)**2)/MYNEWLINE	SOUNDSPEED1_2D(U,IDX))
#define MACHNUMBER2_2D(U,IDX,i,n1,n2)      (sqrt(XVELOCITY2(U,IDX,i,n1,n2)**2+MYNEWLINE YVELOCITY2(U,IDX,i,n1,n2)**2)/MYNEWLINE	SOUNDSPEED2_2D(U,IDX,i,n1,n2))
#define MACHNUMBER3_2D(U,IDX,i,j,n1,n2,n3) (sqrt(XVELOCITY3(U,IDX,i,j,n1,n2,n3)+MYNEWLINE YVELOCITY3(U,IDX,i,j,n1,n2,n3)**2)/MYNEWLINE SOUNDSPEED3_2D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute the Mach number in 3D
#endif
#define MACHNUMBER1_3D(U,IDX)              (sqrt(XVELOCITY1(U,IDX)**2+MYNEWLINE YVELOCITY1(U,IDX)**2+MYNEWLINE ZVELOCITY1(U,IDX)**2)/MYNEWLINE	SOUNDSPEED1_3D(U,IDX))
#define MACHNUMBER2_3D(U,IDX,i,n1,n2)      (sqrt(XVELOCITY2(U,IDX,i,n1,n2)**2+MYNEWLINE YVELOCITY2(U,IDX,i,n1,n2)**2+MYNEWLINE ZVELOCITY2(U,IDX,i,n1,n2)**2)/MYNEWLINE SOUNDSPEED2_3D(U,IDX,i,n1,n2))
#define MACHNUMBER3_3D(U,IDX,i,j,n1,n2,n3) (sqrt(XVELOCITY3(U,IDX,i,j,n1,n2,n3)+MYNEWLINE YVELOCITY3(U,IDX,i,j,n1,n2,n3)**2+MYNEWLINE ZVELOCITY3(U,IDX,i,j,n1,n2,n3)**2)/MYNEWLINE SOUNDSPEED3_3D(U,IDX,i,j,n1,n2,n3))

#if THERMODYN_NDIM == 1
#define MACHNUMBER1 MACHNUMBER1_1D
#define MACHNUMBER2 MACHNUMBER2_1D
#define MACHNUMBER3 MACHNUMBER3_1D
#elif THERMODYN_NDIM == 2
#define MACHNUMBER1 MACHNUMBER1_2D
#define MACHNUMBER2 MACHNUMBER2_2D
#define MACHNUMBER3 MACHNUMBER3_2D
#elif THERMODYN_NDIM == 3
#define MACHNUMBER1 MACHNUMBER1_3D
#define MACHNUMBER2 MACHNUMBER2_3D
#define MACHNUMBER3 MACHNUMBER3_3D
#else
#ifdef THERMODYN_NDIM
#error "Invalid spatial dimension specified!"
#endif
#endif

#endif


