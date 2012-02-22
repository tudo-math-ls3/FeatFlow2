!-*- mode: f90; -*-

#ifndef _MAGNETOHYDRODYNAMICS_OLD_H_
#define _MAGNETOHYDRODYNAMICS_OLD_H_

#if 0
!##############################################################################
!# ****************************************************************************
!# <name> magnetohydrodynamics </name>
!# ****************************************************************************
!#
!# <purpose>
!# This is the magnetohydrodynamic header file
!#
!# </purpose>
!##############################################################################
#endif

#include "../../../../kernel/feat2constants.h"
#include "../../../../kernel/feat2macros.h"
#include "../../../../kernel/System/fmath.h"
#include "../flagship.h"


#if 0
!##############################################################################
! Compile-time constants
!##############################################################################
#endif

#if 0
! Speficy ratio of specific heats (if required)
#endif
#ifndef MAGNETOHYDRODYN_GAMMA
#define MAGNETOHYDRODYN_GAMMA 1.4
#endif

#if 0
! Vacuum permittivity (if required)
#endif
#ifndef MAGNETOHYDRODYN_VACUUM_PERM
#define MAGNETOHYDRODYN_VACUUM_PERM 1.0
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
! Compute y-velocity from conservative variables in 1D, 2D, and 3D
#endif
#define YVELOCITY1(U,IDX)              (YMOMENTUM1(U,IDX)/DENSITY1(U,IDX))
#define YVELOCITY2(U,IDX,i,n1,n2)      (YMOMENTUM2(U,IDX,i,n1,n2)/DENSITY2(U,IDX,i,n1,n2))
#define YVELOCITY3(U,IDX,i,j,n1,n2,n3) (YMOMENTUM3(U,IDX,i,j,n1,n2,n3)/DENSITY3(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute z-velocity from conservative variables in 1D, 2D, and 3D
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
! Compute magnitude of velocity field in 1D, 2D, and 3D
#endif
#define VELMAGNITUDE1(U,IDX)              (sqrt(POW(XVELOCITY1(U,IDX),2)+MYNEWLINE POW(YVELOCITY1(U,IDX),2)+MYNEWLINE POW(ZVELOCITY1(U,IDX),2)))
#define VELMAGNITUDE2(U,IDX,i,n1,n2)      (sqrt(POW(XVELOCITY2(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YVELOCITY2(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZVELOCITY2(U,IDX,i,n1,n2),2)))
#define VELMAGNITUDE3(U,IDX,i,j,n1,n2,n3) (sqrt(POW(XVELOCITY3(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YVELOCITY3(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZVELOCITY3(U,IDX,i,j,n1,n2,n3),2)))


#if 0
!##############################################################################
! Magnitude of magnetic field
!##############################################################################
#endif

#if 0
! Compute magnitude of magnetic field in 1D, 2D, and 3D
#endif
#define MAGFIELDMAGNITUDE1(U,IDX)              (sqrt(POW(XMAGFIELD1(U,IDX),2)+MYNEWLINE POW(YMAGFIELD1(U,IDX),2)+MYNEWLINE POW(ZMAGFIELD1(U,IDX),2)))
#define MAGFIELDMAGNITUDE2(U,IDX,i,n1,n2)      (sqrt(POW(XMAGFIELD2(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YMAGFIELD2(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZMAGFIELD2(U,IDX,i,n1,n2),2)))
#define MAGFIELDMAGNITUDE3(U,IDX,i,j,n1,n2,n3) (sqrt(POW(XMAGFIELD3(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YMAGFIELD3(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZMAGFIELD3(U,IDX,i,j,n1,n2,n3),2)))


#if 0
!##############################################################################
! Pressure
!##############################################################################
#endif

#if 0
! Compute pressure from conservative variables in 1D, 2D, and 3D
#endif
#define PRESSURE1(U,IDX)              (((RCONST(MAGNETOHYDRODYN_GAMMA))-RCONST(1.0))*INTERNALENERGY1(U,IDX))
#define PRESSURE2(U,IDX,i,n1,n2)      (((RCONST(MAGNETOHYDRODYN_GAMMA))-RCONST(1.0))*INTERNALENERGY2(U,IDX,i,n1,n2))
#define PRESSURE3(U,IDX,i,j,n1,n2,n3) (((RCONST(MAGNETOHYDRODYN_GAMMA))-RCONST(1.0))*INTERNALENERGY3(U,IDX,i,j,n1,n2,n3))


#if 0
!##############################################################################
! Magnetic pressure
!##############################################################################
#endif

#define MAGNETICPRESSURE1    MAGNETICENERGY1
#define MAGNETICPRESSURE2    MAGNETICENERGY2
#define MAGNETICPRESSURE3    MAGNETICENERGY3


#if 0
!##############################################################################
! Total pressure
!##############################################################################
#endif

#if 0
! Compute total pressure from conservative variables in 1D, 2D, and 3D
#endif
#define TOTALPRESSURE1(U,IDX)              PRESSURE1(U,IDX)+MYNEWLINE MAGNETICPRESSURE1(U,IDX)
#define TOTALPRESSURE2(U,IDX,i,n1,n2)      PRESSURE2(U,IDX,i,n1,n2)+MYNEWLINE MAGNETICPRESSURE2(U,IDX,i,n1,n2)
#define TOTALPRESSURE3(U,IDX,i,j,n1,n2,n3) PRESSURE3(U,IDX,i,j,n1,n2,n3)+MYNEWLINE MAGNETICPRESSURE3(U,IDX,i,j,n1,n2,n3)


#if 0
!##############################################################################
! Kinetic Energy
!##############################################################################
#endif

#if 0
! Compute kinetic energy from conservative variables in 1D, 2D, and 3D
#endif
#define KINETICENERGY1(U,IDX)              (RCONST(0.5)*(POW(XMOMENTUM1(U,IDX),2)+MYNEWLINE POW(YMOMENTUM1(U,IDX),2)+MYNEWLINE POW(ZMOMENTUM1(U,IDX),2))/MYNEWLINE DENSITY1(U,IDX))
#define KINETICENERGY2(U,IDX,i,n1,n2)      (RCONST(0.5)*(POW(XMOMENTUM2(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YMOMENTUM2(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZMOMENTUM2(U,IDX,i,n1,n2),2))/MYNEWLINE DENSITY2(U,IDX,i,n1,n2))
#define KINETICENERGY3(U,IDX,i,j,n1,n2,n3) (RCONST(0.5)*(POW(XMOMENTUM3(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YMOMENTUM3(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZMOMENTUM3(U,IDX,i,j,n1,n2,n3),2))/MYNEWLINE DENSITY3(U,IDX,i,j,n1,n2,n3))


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
! Internal Energy
!##############################################################################
#endif

#if 0
! Compute internal energy from conservative variables in 1D, 2D, and 3D
#endif
#define INTERNALENERGY1(U,IDX)              (TOTALENERGY1(U,IDX)-MYNEWLINE KINETICENERGY1(U,IDX)-MYNEWLINE MAGNETICENERGY1(U,IDX))
#define INTERNALENERGY2(U,IDX,i,n1,n2)      (TOTALENERGY2(U,IDX,i,n1,n2)-MYNEWLINE KINETICENERGY2(U,IDX,i,n1,n2)-MYNEWLINE MAGNETICENERGY2(U,IDX,i,n1,n2))
#define INTERNALENERGY3(U,IDX,i,j,n1,n2,n3) (TOTALENERGY3(U,IDX,i,j,n1,n2,n3)-MYNEWLINE KINETICENERGY3(U,IDX,i,j,n1,n2,n3)-MYNEWLINE MAGNETICENERGY3(U,IDX,i,j,n1,n2,n3))


#if 0
!##############################################################################
! Magnetic energy
!##############################################################################
#endif

#if 0
  ! Compute magnitude energy in 1D, 2D, and 3D
#endif
#define MAGNETICENERGY1(U,IDX)              (POW(XMAGFIELD1(U,IDX),2)+MYNEWLINE POW(YMAGFIELD1(U,IDX),2)+MYNEWLINE POW(ZMAGFIELD1(U,IDX),2))/(RCONST(2.0)*MAGNETOHYDRODYN_VACUUM_PERM)
#define MAGNETICENERGY2(U,IDX,i,n1,n2)      (POW(XMAGFIELD2(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YMAGFIELD2(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZMAGFIELD2(U,IDX,i,n1,n2),2))/(RCONST(2.0)*MAGNETOHYDRODYN_VACUUM_PERM)
#define MAGNETICENERGY3(U,IDX,i,j,n1,n2,n3) (POW(XMAGFIELD3(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YMAGFIELD3(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZMAGFIELD3(U,IDX,i,j,n1,n2,n3),2))/(RCONST(2.0)*MAGNETOHYDRODYN_VACUUM_PERM)


#if 0
!##############################################################################
! Speed of sound
!
! There are two components in 1D and three components in 2D and 3D
!##############################################################################
#endif

#if 0
  ! Compute speed of sound for a thermally perfect gas in 1D, 2D, and 3D
#endif
#define SOUNDSPEED1(U,IDX)              (sqrt((RCONST(MAGNETOHYDRODYN_GAMMA))*MYNEWLINE PRESSURE1(U,IDX)/MYNEWLINE DENSITY1(U,IDX)))
#define SOUNDSPEED2(U,IDX,i,n1,n2)      (sqrt((RCONST(MAGNETOHYDRODYN_GAMMA))*MYNEWLINE PRESSURE2(U,IDX,i,n1,n2)/MYNEWLINE DENSITY2(U,IDX,i,n1,n2)))
#define SOUNDSPEED3(U,IDX,i,j,n1,n2,n3) (sqrt((RCONST(MAGNETOHYDRODYN_GAMMA))*MYNEWLINE PRESSURE3(U,IDX,i,j,n1,n2,n3)/MYNEWLINE DENSITY3(U,IDX,i,j,n1,n2,n3)))


#if 0
!##############################################################################
! Mach number
!
! There are two components in 1D and three components in 2D and 3D
!##############################################################################
#endif

#if 0
  ! Compute the Mach number in 1D, 2D, and 3D
#endif
#define MACHNUMBER1(U,IDX)              (sqrt(POW(XVELOCITY1(U,IDX),2)+MYNEWLINE POW(YVELOCITY1(U,IDX),2)+MYNEWLINE POW(ZVELOCITY1(U,IDX),2))/MYNEWLINE	SOUNDSPEED1(U,IDX))
#define MACHNUMBER2(U,IDX,i,n1,n2)      (sqrt(POW(XVELOCITY2(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YVELOCITY2(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZVELOCITY2(U,IDX,i,n1,n2),2))/MYNEWLINE SOUNDSPEED2(U,IDX,i,n1,n2))
#define MACHNUMBER3(U,IDX,i,j,n1,n2,n3) (sqrt(POW(XVELOCITY3(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YVELOCITY3(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZVELOCITY3(U,IDX,i,j,n1,n2,n3),2))/MYNEWLINE SOUNDSPEED3(U,IDX,i,j,n1,n2,n3))

#endif
