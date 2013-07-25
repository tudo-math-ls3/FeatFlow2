

#ifndef _MAGNETOHYDRODYNAMICS_H_
#define _MAGNETOHYDRODYNAMICS_H_

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
#define XVELOCITY1_1D(U,IDX)              (XMOMENTUM1_1D(U,IDX)/DENSITY1_1D(U,IDX))
#define XVELOCITY2_1D(U,IDX,i,n1,n2)      (XMOMENTUM2_1D(U,IDX,i,n1,n2)/DENSITY2_1D(U,IDX,i,n1,n2))
#define XVELOCITY3_1D(U,IDX,i,j,n1,n2,n3) (XMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3)/DENSITY3_1D(U,IDX,i,j,n1,n2,n3))

#define XVELOCITY1_2D(U,IDX)              (XMOMENTUM1_2D(U,IDX)/DENSITY1_2D(U,IDX))
#define XVELOCITY2_2D(U,IDX,i,n1,n2)      (XMOMENTUM2_2D(U,IDX,i,n1,n2)/DENSITY2_2D(U,IDX,i,n1,n2))
#define XVELOCITY3_2D(U,IDX,i,j,n1,n2,n3) (XMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)/DENSITY3_2D(U,IDX,i,j,n1,n2,n3))

#define XVELOCITY1_3D(U,IDX)              (XMOMENTUM1_3D(U,IDX)/DENSITY1_3D(U,IDX))
#define XVELOCITY2_3D(U,IDX,i,n1,n2)      (XMOMENTUM2_3D(U,IDX,i,n1,n2)/DENSITY2_3D(U,IDX,i,n1,n2))
#define XVELOCITY3_3D(U,IDX,i,j,n1,n2,n3) (XMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)/DENSITY3_3D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute y-velocity from conservative variables in 1D, 2D, and 3D
#endif
#define YVELOCITY1_1D(U,IDX)              (YMOMENTUM1_1D(U,IDX)/DENSITY1_1D(U,IDX))
#define YVELOCITY2_1D(U,IDX,i,n1,n2)      (YMOMENTUM2_1D(U,IDX,i,n1,n2)/DENSITY2_1D(U,IDX,i,n1,n2))
#define YVELOCITY3_1D(U,IDX,i,j,n1,n2,n3) (YMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3)/DENSITY3_1D(U,IDX,i,j,n1,n2,n3))

#define YVELOCITY1_2D(U,IDX)              (YMOMENTUM1_2D(U,IDX)/DENSITY1_2D(U,IDX))
#define YVELOCITY2_2D(U,IDX,i,n1,n2)      (YMOMENTUM2_2D(U,IDX,i,n1,n2)/DENSITY2_2D(U,IDX,i,n1,n2))
#define YVELOCITY3_2D(U,IDX,i,j,n1,n2,n3) (YMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)/DENSITY3_2D(U,IDX,i,j,n1,n2,n3))

#define YVELOCITY1_3D(U,IDX)              (YMOMENTUM1_3D(U,IDX)/DENSITY1_3D(U,IDX))
#define YVELOCITY2_3D(U,IDX,i,n1,n2)      (YMOMENTUM2_3D(U,IDX,i,n1,n2)/DENSITY2_3D(U,IDX,i,n1,n2))
#define YVELOCITY3_3D(U,IDX,i,j,n1,n2,n3) (YMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)/DENSITY3_3D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute z-velocity from conservative variables in 1D, 2D, and 3D
#endif
#define ZVELOCITY1_1D(U,IDX)              (ZMOMENTUM1_1D(U,IDX)/DENSITY1_1D(U,IDX))
#define ZVELOCITY2_1D(U,IDX,i,n1,n2)      (ZMOMENTUM2_1D(U,IDX,i,n1,n2)/DENSITY2_1D(U,IDX,i,n1,n2))
#define ZVELOCITY3_1D(U,IDX,i,j,n1,n2,n3) (ZMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3)/DENSITY3_1D(U,IDX,i,j,n1,n2,n3))

#define ZVELOCITY1_2D(U,IDX)              (ZMOMENTUM1_2D(U,IDX)/DENSITY1_2D(U,IDX))
#define ZVELOCITY2_2D(U,IDX,i,n1,n2)      (ZMOMENTUM2_2D(U,IDX,i,n1,n2)/DENSITY2_2D(U,IDX,i,n1,n2))
#define ZVELOCITY3_2D(U,IDX,i,j,n1,n2,n3) (ZMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)/DENSITY3_2D(U,IDX,i,j,n1,n2,n3))

#define ZVELOCITY1_3D(U,IDX)              (ZMOMENTUM1_3D(U,IDX)/DENSITY1_3D(U,IDX))
#define ZVELOCITY2_3D(U,IDX,i,n1,n2)      (ZMOMENTUM2_3D(U,IDX,i,n1,n2)/DENSITY2_3D(U,IDX,i,n1,n2))
#define ZVELOCITY3_3D(U,IDX,i,j,n1,n2,n3) (ZMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)/DENSITY3_3D(U,IDX,i,j,n1,n2,n3))


#if 0
!##############################################################################
! Magnitude of velocity field
!##############################################################################
#endif

#if 0
! Compute magnitude of velocity field in 1D, 2D, and 3D
#endif
#define VELMAGNITUDE1_1D(U,IDX)              (sqrt(POW(XVELOCITY1_1D(U,IDX),2)+MYNEWLINE POW(YVELOCITY1_1D(U,IDX),2)+MYNEWLINE POW(ZVELOCITY1_1D(U,IDX),2)))
#define VELMAGNITUDE2_1D(U,IDX,i,n1,n2)      (sqrt(POW(XVELOCITY2_1D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YVELOCITY2_1D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZVELOCITY2_1D(U,IDX,i,n1,n2),2)))
#define VELMAGNITUDE3_1D(U,IDX,i,j,n1,n2,n3) (sqrt(POW(XVELOCITY3_1D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YVELOCITY3_1D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZVELOCITY3_1D(U,IDX,i,j,n1,n2,n3),2)))

#define VELMAGNITUDE1_2D(U,IDX)              (sqrt(POW(XVELOCITY1_2D(U,IDX),2)+MYNEWLINE POW(YVELOCITY1_2D(U,IDX),2)+MYNEWLINE POW(ZVELOCITY1_2D(U,IDX),2)))
#define VELMAGNITUDE2_2D(U,IDX,i,n1,n2)      (sqrt(POW(XVELOCITY2_2D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YVELOCITY2_2D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZVELOCITY2_2D(U,IDX,i,n1,n2),2)))
#define VELMAGNITUDE3_2D(U,IDX,i,j,n1,n2,n3) (sqrt(POW(XVELOCITY3_2D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YVELOCITY3_2D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZVELOCITY3_2D(U,IDX,i,j,n1,n2,n3),2)))
      
#define VELMAGNITUDE1_3D(U,IDX)              (sqrt(POW(XVELOCITY1_3D(U,IDX),2)+MYNEWLINE POW(YVELOCITY1_3D(U,IDX),2)+MYNEWLINE POW(ZVELOCITY1_3D(U,IDX),2)))
#define VELMAGNITUDE2_3D(U,IDX,i,n1,n2)      (sqrt(POW(XVELOCITY2_3D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YVELOCITY2_3D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZVELOCITY2_3D(U,IDX,i,n1,n2),2)))
#define VELMAGNITUDE3_3D(U,IDX,i,j,n1,n2,n3) (sqrt(POW(XVELOCITY3_3D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YVELOCITY3_3D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZVELOCITY3_3D(U,IDX,i,j,n1,n2,n3),2)))


#if 0
!##############################################################################
! Magnitude of magnetic field
!##############################################################################
#endif

#if 0
! Compute magnitude of magnetic field in 1D, 2D, and 3D
#endif
#define MAGFIELDMAGNITUDE1_1D(U,IDX)              (sqrt(POW(XMAGFIELD1_1D(U,IDX),2)+MYNEWLINE POW(YMAGFIELD1_1D(U,IDX),2)+MYNEWLINE POW(ZMAGFIELD1_1D(U,IDX),2)))
#define MAGFIELDMAGNITUDE2_1D(U,IDX,i,n1,n2)      (sqrt(POW(XMAGFIELD2_1D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YMAGFIELD2_1D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZMAGFIELD2_1D(U,IDX,i,n1,n2),2)))
#define MAGFIELDMAGNITUDE3_1D(U,IDX,i,j,n1,n2,n3) (sqrt(POW(XMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3),2)))

#define MAGFIELDMAGNITUDE1_2D(U,IDX)              (sqrt(POW(XMAGFIELD1_2D(U,IDX),2)+MYNEWLINE POW(YMAGFIELD1_2D(U,IDX),2)+MYNEWLINE POW(ZMAGFIELD1_2D(U,IDX),2)))
#define MAGFIELDMAGNITUDE2_2D(U,IDX,i,n1,n2)      (sqrt(POW(XMAGFIELD2_2D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YMAGFIELD2_2D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZMAGFIELD2_2D(U,IDX,i,n1,n2),2)))
#define MAGFIELDMAGNITUDE3_2D(U,IDX,i,j,n1,n2,n3) (sqrt(POW(XMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3),2)))

#define MAGFIELDMAGNITUDE1_3D(U,IDX)              (sqrt(POW(XMAGFIELD1_3D(U,IDX),2)+MYNEWLINE POW(YMAGFIELD1_3D(U,IDX),2)+MYNEWLINE POW(ZMAGFIELD1_3D(U,IDX),2)))
#define MAGFIELDMAGNITUDE2_3D(U,IDX,i,n1,n2)      (sqrt(POW(XMAGFIELD2_3D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YMAGFIELD2_3D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZMAGFIELD2_3D(U,IDX,i,n1,n2),2)))
#define MAGFIELDMAGNITUDE3_3D(U,IDX,i,j,n1,n2,n3) (sqrt(POW(XMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3),2)))


#if 0
!##############################################################################
! Pressure
!##############################################################################
#endif

#if 0
! Compute pressure from conservative variables in 1D, 2D, and 3D
#endif
#define PRESSURE1_1D(U,IDX)              (((DCONST(MAGNETOHYDRODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY1_1D(U,IDX))
#define PRESSURE2_1D(U,IDX,i,n1,n2)      (((DCONST(MAGNETOHYDRODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY2_1D(U,IDX,i,n1,n2))
#define PRESSURE3_1D(U,IDX,i,j,n1,n2,n3) (((DCONST(MAGNETOHYDRODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY3_1D(U,IDX,i,j,n1,n2,n3))

#define PRESSURE1_2D(U,IDX)              (((DCONST(MAGNETOHYDRODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY1_2D(U,IDX))
#define PRESSURE2_2D(U,IDX,i,n1,n2)      (((DCONST(MAGNETOHYDRODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY2_2D(U,IDX,i,n1,n2))
#define PRESSURE3_2D(U,IDX,i,j,n1,n2,n3) (((DCONST(MAGNETOHYDRODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY3_2D(U,IDX,i,j,n1,n2,n3))

#define PRESSURE1_3D(U,IDX)              (((DCONST(MAGNETOHYDRODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY1_3D(U,IDX))
#define PRESSURE2_3D(U,IDX,i,n1,n2)      (((DCONST(MAGNETOHYDRODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY2_3D(U,IDX,i,n1,n2))
#define PRESSURE3_3D(U,IDX,i,j,n1,n2,n3) (((DCONST(MAGNETOHYDRODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY3_3D(U,IDX,i,j,n1,n2,n3))


#if 0
!##############################################################################
! Magnetic pressure
!##############################################################################
#endif

#define MAGNETICPRESSURE1_1D    MAGNETICENERGY1_1D
#define MAGNETICPRESSURE2_1D    MAGNETICENERGY2_1D
#define MAGNETICPRESSURE3_1D    MAGNETICENERGY3_1D

#define MAGNETICPRESSURE1_2D    MAGNETICENERGY1_2D
#define MAGNETICPRESSURE2_2D    MAGNETICENERGY2_2D
#define MAGNETICPRESSURE3_2D    MAGNETICENERGY3_2D

#define MAGNETICPRESSURE1_3D    MAGNETICENERGY1_3D
#define MAGNETICPRESSURE2_3D    MAGNETICENERGY2_3D
#define MAGNETICPRESSURE3_3D    MAGNETICENERGY3_3D


#if 0
!##############################################################################
! Total pressure
!##############################################################################
#endif

#if 0
! Compute total pressure from conservative variables in 1D, 2D, and 3D
#endif
#define TOTALPRESSURE1_1D(U,IDX)              PRESSURE1_1D(U,IDX)+MYNEWLINE MAGNETICPRESSURE1_1D(U,IDX)
#define TOTALPRESSURE2_1D(U,IDX,i,n1,n2)      PRESSURE2_1D(U,IDX,i,n1,n2)+MYNEWLINE MAGNETICPRESSURE2_1D(U,IDX,i,n1,n2)
#define TOTALPRESSURE3_1D(U,IDX,i,j,n1,n2,n3) PRESSURE3_1D(U,IDX,i,j,n1,n2,n3)+MYNEWLINE MAGNETICPRESSURE3_1D(U,IDX,i,j,n1,n2,n3)

#define TOTALPRESSURE1_2D(U,IDX)              PRESSURE1_2D(U,IDX)+MYNEWLINE MAGNETICPRESSURE1_2D(U,IDX)
#define TOTALPRESSURE2_2D(U,IDX,i,n1,n2)      PRESSURE2_2D(U,IDX,i,n1,n2)+MYNEWLINE MAGNETICPRESSURE2_2D(U,IDX,i,n1,n2)
#define TOTALPRESSURE3_2D(U,IDX,i,j,n1,n2,n3) PRESSURE3_2D(U,IDX,i,j,n1,n2,n3)+MYNEWLINE MAGNETICPRESSURE3_2D(U,IDX,i,j,n1,n2,n3)

#define TOTALPRESSURE1_3D(U,IDX)              PRESSURE1_3D(U,IDX)+MYNEWLINE MAGNETICPRESSURE1_3D(U,IDX)
#define TOTALPRESSURE2_3D(U,IDX,i,n1,n2)      PRESSURE2_3D(U,IDX,i,n1,n2)+MYNEWLINE MAGNETICPRESSURE2_3D(U,IDX,i,n1,n2)
#define TOTALPRESSURE3_3D(U,IDX,i,j,n1,n2,n3) PRESSURE3_3D(U,IDX,i,j,n1,n2,n3)+MYNEWLINE MAGNETICPRESSURE3_3D(U,IDX,i,j,n1,n2,n3)


#if 0
!##############################################################################
! Kinetic Energy
!##############################################################################
#endif

#if 0
! Compute kinetic energy from conservative variables in 1D, 2D, and 3D
#endif
#define KINETICENERGY1_1D(U,IDX)              (DCONST(0.5)*(POW(XMOMENTUM1_1D(U,IDX),2)+MYNEWLINE POW(YMOMENTUM1_1D(U,IDX),2)+MYNEWLINE POW(ZMOMENTUM1_1D(U,IDX),2))/MYNEWLINE DENSITY1_1D(U,IDX))
#define KINETICENERGY2_1D(U,IDX,i,n1,n2)      (DCONST(0.5)*(POW(XMOMENTUM2_1D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YMOMENTUM2_1D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZMOMENTUM2_1D(U,IDX,i,n1,n2),2))/MYNEWLINE DENSITY2_1D(U,IDX,i,n1,n2))
#define KINETICENERGY3_1D(U,IDX,i,j,n1,n2,n3) (DCONST(0.5)*(POW(XMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3),2))/MYNEWLINE DENSITY3_1D(U,IDX,i,j,n1,n2,n3))

#define KINETICENERGY1_2D(U,IDX)              (DCONST(0.5)*(POW(XMOMENTUM1_2D(U,IDX),2)+MYNEWLINE POW(YMOMENTUM1_2D(U,IDX),2)+MYNEWLINE POW(ZMOMENTUM1_2D(U,IDX),2))/MYNEWLINE DENSITY1_2D(U,IDX))
#define KINETICENERGY2_2D(U,IDX,i,n1,n2)      (DCONST(0.5)*(POW(XMOMENTUM2_2D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YMOMENTUM2_2D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZMOMENTUM2_2D(U,IDX,i,n1,n2),2))/MYNEWLINE DENSITY2_2D(U,IDX,i,n1,n2))
#define KINETICENERGY3_2D(U,IDX,i,j,n1,n2,n3) (DCONST(0.5)*(POW(XMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3),2))/MYNEWLINE DENSITY3_2D(U,IDX,i,j,n1,n2,n3))
      
#define KINETICENERGY1_3D(U,IDX)              (DCONST(0.5)*(POW(XMOMENTUM1_3D(U,IDX),2)+MYNEWLINE POW(YMOMENTUM1_3D(U,IDX),2)+MYNEWLINE POW(ZMOMENTUM1_3D(U,IDX),2))/MYNEWLINE DENSITY1_3D(U,IDX))
#define KINETICENERGY2_3D(U,IDX,i,n1,n2)      (DCONST(0.5)*(POW(XMOMENTUM2_3D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YMOMENTUM2_3D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZMOMENTUM2_3D(U,IDX,i,n1,n2),2))/MYNEWLINE DENSITY2_3D(U,IDX,i,n1,n2))
#define KINETICENERGY3_3D(U,IDX,i,j,n1,n2,n3) (DCONST(0.5)*(POW(XMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3),2))/MYNEWLINE DENSITY3_3D(U,IDX,i,j,n1,n2,n3))


#if 0
!##############################################################################
! Specific total energy, aka, total energy per unit mass
!##############################################################################
#endif

#define SPECIFICTOTALENERGY1_1D(U,IDX)              (TOTALENERGY1_1D(U,IDX)/DENSITY1_1D(U,IDX))
#define SPECIFICTOTALENERGY2_1D(U,IDX,i,n1,n2)      (TOTALENERGY2_1D(U,IDX,i,n1,n2)/DENSITY2_1D(U,IDX,i,n1,n2))
#define SPECIFICTOTALENERGY3_1D(U,IDX,i,j,n1,n2,n3) (TOTALENERGY3_1D(U,IDX,i,j,n1,n2,n3)/DENSITY3_1D(U,IDX,i,j,n1,n2,n3))

#define SPECIFICTOTALENERGY1_2D(U,IDX)              (TOTALENERGY1_2D(U,IDX)/DENSITY1_2D(U,IDX))
#define SPECIFICTOTALENERGY2_2D(U,IDX,i,n1,n2)      (TOTALENERGY2_2D(U,IDX,i,n1,n2)/DENSITY2_2D(U,IDX,i,n1,n2))
#define SPECIFICTOTALENERGY3_2D(U,IDX,i,j,n1,n2,n3) (TOTALENERGY3_2D(U,IDX,i,j,n1,n2,n3)/DENSITY3_2D(U,IDX,i,j,n1,n2,n3))

#define SPECIFICTOTALENERGY1_3D(U,IDX)              (TOTALENERGY1_3D(U,IDX)/DENSITY1_3D(U,IDX))
#define SPECIFICTOTALENERGY2_3D(U,IDX,i,n1,n2)      (TOTALENERGY2_3D(U,IDX,i,n1,n2)/DENSITY2_3D(U,IDX,i,n1,n2))
#define SPECIFICTOTALENERGY3_3D(U,IDX,i,j,n1,n2,n3) (TOTALENERGY3_3D(U,IDX,i,j,n1,n2,n3)/DENSITY3_3D(U,IDX,i,j,n1,n2,n3))


#if 0
!##############################################################################
! Internal Energy
!##############################################################################
#endif

#if 0
! Compute internal energy from conservative variables in 1D, 2D, and 3D
#endif
#define INTERNALENERGY1_1D(U,IDX)              (TOTALENERGY1_1D(U,IDX)-MYNEWLINE KINETICENERGY1_1D(U,IDX)-MYNEWLINE MAGNETICENERGY1_1D(U,IDX))
#define INTERNALENERGY2_1D(U,IDX,i,n1,n2)      (TOTALENERGY2_1D(U,IDX,i,n1,n2)-MYNEWLINE KINETICENERGY2_1D(U,IDX,i,n1,n2)-MYNEWLINE MAGNETICENERGY2_1D(U,IDX,i,n1,n2))
#define INTERNALENERGY3_1D(U,IDX,i,j,n1,n2,n3) (TOTALENERGY3_1D(U,IDX,i,j,n1,n2,n3)-MYNEWLINE KINETICENERGY3_1D(U,IDX,i,j,n1,n2,n3)-MYNEWLINE MAGNETICENERGY3_1D(U,IDX,i,j,n1,n2,n3))

#define INTERNALENERGY1_2D(U,IDX)              (TOTALENERGY1_2D(U,IDX)-MYNEWLINE KINETICENERGY1_2D(U,IDX)-MYNEWLINE MAGNETICENERGY1_2D(U,IDX))
#define INTERNALENERGY2_2D(U,IDX,i,n1,n2)      (TOTALENERGY2_2D(U,IDX,i,n1,n2)-MYNEWLINE KINETICENERGY2_2D(U,IDX,i,n1,n2)-MYNEWLINE MAGNETICENERGY2_2D(U,IDX,i,n1,n2))
#define INTERNALENERGY3_2D(U,IDX,i,j,n1,n2,n3) (TOTALENERGY3_2D(U,IDX,i,j,n1,n2,n3)-MYNEWLINE KINETICENERGY3_2D(U,IDX,i,j,n1,n2,n3)-MYNEWLINE MAGNETICENERGY3_2D(U,IDX,i,j,n1,n2,n3))
      
#define INTERNALENERGY1_3D(U,IDX)              (TOTALENERGY1_3D(U,IDX)-MYNEWLINE KINETICENERGY1_3D(U,IDX)-MYNEWLINE MAGNETICENERGY1_3D(U,IDX))
#define INTERNALENERGY2_3D(U,IDX,i,n1,n2)      (TOTALENERGY2_3D(U,IDX,i,n1,n2)-MYNEWLINE KINETICENERGY2_3D(U,IDX,i,n1,n2)-MYNEWLINE MAGNETICENERGY2_3D(U,IDX,i,n1,n2))
#define INTERNALENERGY3_3D(U,IDX,i,j,n1,n2,n3) (TOTALENERGY3_3D(U,IDX,i,j,n1,n2,n3)-MYNEWLINE KINETICENERGY3_3D(U,IDX,i,j,n1,n2,n3)-MYNEWLINE MAGNETICENERGY3_3D(U,IDX,i,j,n1,n2,n3))


#if 0
!##############################################################################
! Magnetic energy
!##############################################################################
#endif

#if 0
  ! Compute magnitude energy in 1D, 2D, and 3D
#endif
#define MAGNETICENERGY1_1D(U,IDX)              (POW(XMAGFIELD1_1D(U,IDX),2)+MYNEWLINE POW(YMAGFIELD1_1D(U,IDX),2)+MYNEWLINE POW(ZMAGFIELD1_1D(U,IDX),2))/(DCONST(2.0)*MAGNETOHYDRODYN_VACUUM_PERM)
#define MAGNETICENERGY2_1D(U,IDX,i,n1,n2)      (POW(XMAGFIELD2_1D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YMAGFIELD2_1D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZMAGFIELD2_1D(U,IDX,i,n1,n2),2))/(DCONST(2.0)*MAGNETOHYDRODYN_VACUUM_PERM)
#define MAGNETICENERGY3_1D(U,IDX,i,j,n1,n2,n3) (POW(XMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZMAGFIELD3_1D(U,IDX,i,j,n1,n2,n3),2))/(DCONST(2.0)*MAGNETOHYDRODYN_VACUUM_PERM)

#define MAGNETICENERGY1_2D(U,IDX)              (POW(XMAGFIELD1_2D(U,IDX),2)+MYNEWLINE POW(YMAGFIELD1_2D(U,IDX),2)+MYNEWLINE POW(ZMAGFIELD1_2D(U,IDX),2))/(DCONST(2.0)*MAGNETOHYDRODYN_VACUUM_PERM)
#define MAGNETICENERGY2_2D(U,IDX,i,n1,n2)      (POW(XMAGFIELD2_2D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YMAGFIELD2_2D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZMAGFIELD2_2D(U,IDX,i,n1,n2),2))/(DCONST(2.0)*MAGNETOHYDRODYN_VACUUM_PERM)
#define MAGNETICENERGY3_2D(U,IDX,i,j,n1,n2,n3) (POW(XMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZMAGFIELD3_2D(U,IDX,i,j,n1,n2,n3),2))/(DCONST(2.0)*MAGNETOHYDRODYN_VACUUM_PERM)

#define MAGNETICENERGY1_3D(U,IDX)              (POW(XMAGFIELD1_3D(U,IDX),2)+MYNEWLINE POW(YMAGFIELD1_3D(U,IDX),2)+MYNEWLINE POW(ZMAGFIELD1_3D(U,IDX),2))/(DCONST(2.0)*MAGNETOHYDRODYN_VACUUM_PERM)
#define MAGNETICENERGY2_3D(U,IDX,i,n1,n2)      (POW(XMAGFIELD2_3D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YMAGFIELD2_3D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZMAGFIELD2_3D(U,IDX,i,n1,n2),2))/(DCONST(2.0)*MAGNETOHYDRODYN_VACUUM_PERM)
#define MAGNETICENERGY3_3D(U,IDX,i,j,n1,n2,n3) (POW(XMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZMAGFIELD3_3D(U,IDX,i,j,n1,n2,n3),2))/(DCONST(2.0)*MAGNETOHYDRODYN_VACUUM_PERM)


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
#define SOUNDSPEED1_1D(U,IDX)              (sqrt((DCONST(MAGNETOHYDRODYN_GAMMA))*MYNEWLINE PRESSURE1_1D(U,IDX)/MYNEWLINE DENSITY1_1D(U,IDX)))
#define SOUNDSPEED2_1D(U,IDX,i,n1,n2)      (sqrt((DCONST(MAGNETOHYDRODYN_GAMMA))*MYNEWLINE PRESSURE2_1D(U,IDX,i,n1,n2)/MYNEWLINE DENSITY2_1D(U,IDX,i,n1,n2)))
#define SOUNDSPEED3_1D(U,IDX,i,j,n1,n2,n3) (sqrt((DCONST(MAGNETOHYDRODYN_GAMMA))*MYNEWLINE PRESSURE3_1D(U,IDX,i,j,n1,n2,n3)/MYNEWLINE DENSITY3_1D(U,IDX,i,j,n1,n2,n3)))

#define SOUNDSPEED1_2D(U,IDX)              (sqrt((DCONST(MAGNETOHYDRODYN_GAMMA))*MYNEWLINE PRESSURE1_2D(U,IDX)/MYNEWLINE DENSITY1_2D(U,IDX)))
#define SOUNDSPEED2_2D(U,IDX,i,n1,n2)      (sqrt((DCONST(MAGNETOHYDRODYN_GAMMA))*MYNEWLINE PRESSURE2_2D(U,IDX,i,n1,n2)/MYNEWLINE DENSITY2_2D(U,IDX,i,n1,n2)))
#define SOUNDSPEED3_2D(U,IDX,i,j,n1,n2,n3) (sqrt((DCONST(MAGNETOHYDRODYN_GAMMA))*MYNEWLINE PRESSURE3_2D(U,IDX,i,j,n1,n2,n3)/MYNEWLINE DENSITY3_2D(U,IDX,i,j,n1,n2,n3)))

#define SOUNDSPEED1_3D(U,IDX)              (sqrt((DCONST(MAGNETOHYDRODYN_GAMMA))*MYNEWLINE PRESSURE1_3D(U,IDX)/MYNEWLINE DENSITY1_3D(U,IDX)))
#define SOUNDSPEED2_3D(U,IDX,i,n1,n2)      (sqrt((DCONST(MAGNETOHYDRODYN_GAMMA))*MYNEWLINE PRESSURE2_3D(U,IDX,i,n1,n2)/MYNEWLINE DENSITY2_3D(U,IDX,i,n1,n2)))
#define SOUNDSPEED3_3D(U,IDX,i,j,n1,n2,n3) (sqrt((DCONST(MAGNETOHYDRODYN_GAMMA))*MYNEWLINE PRESSURE3_3D(U,IDX,i,j,n1,n2,n3)/MYNEWLINE DENSITY3_3D(U,IDX,i,j,n1,n2,n3)))


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
#define MACHNUMBER1_1D(U,IDX)              (sqrt(POW(XVELOCITY1_1D(U,IDX),2)+MYNEWLINE POW(YVELOCITY1_1D(U,IDX),2)+MYNEWLINE POW(ZVELOCITY1_1D(U,IDX),2))/MYNEWLINE	SOUNDSPEED1_1D(U,IDX))
#define MACHNUMBER2_1D(U,IDX,i,n1,n2)      (sqrt(POW(XVELOCITY2_1D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YVELOCITY2_1D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZVELOCITY2_1D(U,IDX,i,n1,n2),2))/MYNEWLINE SOUNDSPEED2_1D(U,IDX,i,n1,n2))
#define MACHNUMBER3_1D(U,IDX,i,j,n1,n2,n3) (sqrt(POW(XVELOCITY3_1D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YVELOCITY3_1D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZVELOCITY3_1D(U,IDX,i,j,n1,n2,n3),2))/MYNEWLINE SOUNDSPEED3_1D(U,IDX,i,j,n1,n2,n3))

#define MACHNUMBER1_2D(U,IDX)              (sqrt(POW(XVELOCITY1_2D(U,IDX),2)+MYNEWLINE POW(YVELOCITY1_2D(U,IDX),2)+MYNEWLINE POW(ZVELOCITY1_2D(U,IDX),2))/MYNEWLINE	SOUNDSPEED1_2D(U,IDX))
#define MACHNUMBER2_2D(U,IDX,i,n1,n2)      (sqrt(POW(XVELOCITY2_2D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YVELOCITY2_2D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZVELOCITY2_2D(U,IDX,i,n1,n2),2))/MYNEWLINE SOUNDSPEED2_2D(U,IDX,i,n1,n2))
#define MACHNUMBER3_2D(U,IDX,i,j,n1,n2,n3) (sqrt(POW(XVELOCITY3_2D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YVELOCITY3_2D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZVELOCITY3_2D(U,IDX,i,j,n1,n2,n3),2))/MYNEWLINE SOUNDSPEED3_2D(U,IDX,i,j,n1,n2,n3))

#define MACHNUMBER1_3D(U,IDX)              (sqrt(POW(XVELOCITY1_3D(U,IDX),2)+MYNEWLINE POW(YVELOCITY1_3D(U,IDX),2)+MYNEWLINE POW(ZVELOCITY1_3D(U,IDX),2))/MYNEWLINE	SOUNDSPEED1_3D(U,IDX))
#define MACHNUMBER2_3D(U,IDX,i,n1,n2)      (sqrt(POW(XVELOCITY2_3D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(YVELOCITY2_3D(U,IDX,i,n1,n2),2)+MYNEWLINE POW(ZVELOCITY2_3D(U,IDX,i,n1,n2),2))/MYNEWLINE SOUNDSPEED2_3D(U,IDX,i,n1,n2))
#define MACHNUMBER3_3D(U,IDX,i,j,n1,n2,n3) (sqrt(POW(XVELOCITY3_3D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(YVELOCITY3_3D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE POW(ZVELOCITY3_3D(U,IDX,i,j,n1,n2,n3),2))/MYNEWLINE SOUNDSPEED3_3D(U,IDX,i,j,n1,n2,n3))

#endif
