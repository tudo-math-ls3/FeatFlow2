

#ifndef _THERMODYNAMICS_H_
#define _THERMODYNAMICS_H_

#if 0
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

#include "flagship.h"


#if 0
!##############################################################################
! Compile-time constants
!##############################################################################
#endif

#if 0
! Speficy ratio of specific heats (if required)
#endif
#ifndef THERMODYN_GAMMA
#define THERMODYN_GAMMA 1.4
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
! Compute y-velocity from conservative variables in 2D, and 3D
#endif
#define YVELOCITY1_2D(U,IDX)              (YMOMENTUM1_2D(U,IDX)/DENSITY1_2D(U,IDX))
#define YVELOCITY2_2D(U,IDX,i,n1,n2)      (YMOMENTUM2_2D(U,IDX,i,n1,n2)/DENSITY2_2D(U,IDX,i,n1,n2))
#define YVELOCITY3_2D(U,IDX,i,j,n1,n2,n3) (YMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3)/DENSITY3_2D(U,IDX,i,j,n1,n2,n3))

#define YVELOCITY1_3D(U,IDX)              (YMOMENTUM1_3D(U,IDX)/DENSITY1_3D(U,IDX))
#define YVELOCITY2_3D(U,IDX,i,n1,n2)      (YMOMENTUM2_3D(U,IDX,i,n1,n2)/DENSITY2_3D(U,IDX,i,n1,n2))
#define YVELOCITY3_3D(U,IDX,i,j,n1,n2,n3) (YMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)/DENSITY3_3D(U,IDX,i,j,n1,n2,n3))
      
#if 0
! Compute z-velocity from conservative variables in 3D
#endif
#define ZVELOCITY1_3D(U,IDX)              (ZMOMENTUM1_3D(U,IDX)/DENSITY1_3D(U,IDX))
#define ZVELOCITY2_3D(U,IDX,i,n1,n2)      (ZMOMENTUM2_3D(U,IDX,i,n1,n2)/DENSITY2_3D(U,IDX,i,n1,n2))
#define ZVELOCITY3_3D(U,IDX,i,j,n1,n2,n3) (ZMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3)/DENSITY3_3D(U,IDX,i,j,n1,n2,n3))


#if 0
!##############################################################################
! Magnitude of velocity field
!##############################################################################
#endif

#if 0
! Compute magnitude of velocity field in 1D
#endif
#define VELMAGNITUDE1_1D(U,IDX)              (abs(XVELOCITY1_1D(U,IDX)))
#define VELMAGNITUDE2_1D(U,IDX,i,n1,n2)      (abs(XVELOCITY2_1D(U,IDX,i,n1,n2)))
#define VELMAGNITUDE3_1D(U,IDX,i,j,n1,n2,n3) (abs(XVELOCITY3_1D(U,IDX,i,j,n1,n2,n3)))

#if 0
! Compute magnitude of velocity field in 2D
#endif
#define VELMAGNITUDE1_2D(U,IDX)              (sqrt(POW(XVELOCITY1_2D(U,IDX),2)+MYNEWLINE \
                                                   POW(YVELOCITY1_2D(U,IDX),2)))
#define VELMAGNITUDE2_2D(U,IDX,i,n1,n2)      (sqrt(POW(XVELOCITY2_2D(U,IDX,i,n1,n2),2)+MYNEWLINE \
                                                   POW(YVELOCITY2_2D(U,IDX,i,n1,n2),2)))
#define VELMAGNITUDE3_2D(U,IDX,i,j,n1,n2,n3) (sqrt(POW(XVELOCITY3_2D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE \
                                                   POW(YVELOCITY3_2D(U,IDX,i,j,n1,n2,n3),2)))

#if 0
! Compute magnitude of velocity field in 3D
#endif
#define VELMAGNITUDE1_3D(U,IDX)              (sqrt(POW(XVELOCITY1_3D(U,IDX),2)+MYNEWLINE \
                                                   POW(YVELOCITY1_3D(U,IDX),2)+MYNEWLINE \
                                                   POW(ZVELOCITY1_3D(U,IDX),2)))
#define VELMAGNITUDE2_3D(U,IDX,i,n1,n2)      (sqrt(POW(XVELOCITY2_3D(U,IDX,i,n1,n2),2)+MYNEWLINE \
                                                   POW(YVELOCITY2_3D(U,IDX,i,n1,n2),2)+MYNEWLINE \
                                                   POW(ZVELOCITY2_3D(U,IDX,i,n1,n2),2)))
#define VELMAGNITUDE3_3D(U,IDX,i,j,n1,n2,n3) (sqrt(POW(XVELOCITY3_3D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE \
                                                   POW(YVELOCITY3_3D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE \
                                                   POW(ZVELOCITY3_3D(U,IDX,i,j,n1,n2,n3),2)))


#if 0
!##############################################################################
! Pressure
!##############################################################################
#endif

#if 0
! Compute pressure from conservative variables in 1D
#endif
#define PRESSURE1_1D(U,IDX)              (((DCONST(THERMODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY1_1D(U,IDX))
#define PRESSURE2_1D(U,IDX,i,n1,n2)      (((DCONST(THERMODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY2_1D(U,IDX,i,n1,n2))
#define PRESSURE3_1D(U,IDX,i,j,n1,n2,n3) (((DCONST(THERMODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY3_1D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute pressure from conservative variables in 2D
#endif
#define PRESSURE1_2D(U,IDX)              (((DCONST(THERMODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY1_2D(U,IDX))
#define PRESSURE2_2D(U,IDX,i,n1,n2)      (((DCONST(THERMODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY2_2D(U,IDX,i,n1,n2))
#define PRESSURE3_2D(U,IDX,i,j,n1,n2,n3) (((DCONST(THERMODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY3_2D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute pressure from conservative variables in 3D
#endif
#define PRESSURE1_3D(U,IDX)              (((DCONST(THERMODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY1_3D(U,IDX))
#define PRESSURE2_3D(U,IDX,i,n1,n2)      (((DCONST(THERMODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY2_3D(U,IDX,i,n1,n2))
#define PRESSURE3_3D(U,IDX,i,j,n1,n2,n3) (((DCONST(THERMODYN_GAMMA))-DCONST(1.0))*INTERNALENERGY3_3D(U,IDX,i,j,n1,n2,n3))


#if 0
!##############################################################################
! Internal Energy
!##############################################################################
#endif

#if 0
! Compute internal energy from conservative variables in 1D
#endif
#define INTERNALENERGY1_1D(U,IDX)              (TOTALENERGY1_1D(U,IDX)-MYNEWLINE \
                                                KINETICENERGY1_1D(U,IDX))
#define INTERNALENERGY2_1D(U,IDX,i,n1,n2)      (TOTALENERGY2_1D(U,IDX,i,n1,n2)-MYNEWLINE \
                                                KINETICENERGY2_1D(U,IDX,i,n1,n2))
#define INTERNALENERGY3_1D(U,IDX,i,j,n1,n2,n3) (TOTALENERGY3_1D(U,IDX,i,j,n1,n2,n3)-MYNEWLINE \
                                                KINETICENERGY3_1D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute internal energy from conservative variables in 2D
#endif
#define INTERNALENERGY1_2D(U,IDX)              (TOTALENERGY1_2D(U,IDX)-MYNEWLINE \
                                                KINETICENERGY1_2D(U,IDX))
#define INTERNALENERGY2_2D(U,IDX,i,n1,n2)      (TOTALENERGY2_2D(U,IDX,i,n1,n2)-MYNEWLINE \
                                                KINETICENERGY2_2D(U,IDX,i,n1,n2))
#define INTERNALENERGY3_2D(U,IDX,i,j,n1,n2,n3) (TOTALENERGY3_2D(U,IDX,i,j,n1,n2,n3)-MYNEWLINE \
                                                KINETICENERGY3_2D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute internal energy from conservative variables in 3D
#endif
#define INTERNALENERGY1_3D(U,IDX)              (TOTALENERGY1_3D(U,IDX)-MYNEWLINE \
                                                KINETICENERGY1_3D(U,IDX))
#define INTERNALENERGY2_3D(U,IDX,i,n1,n2)      (TOTALENERGY2_3D(U,IDX,i,n1,n2)-MYNEWLINE \
                                                KINETICENERGY2_3D(U,IDX,i,n1,n2))
#define INTERNALENERGY3_3D(U,IDX,i,j,n1,n2,n3) (TOTALENERGY3_3D(U,IDX,i,j,n1,n2,n3)-MYNEWLINE \
                                                KINETICENERGY3_3D(U,IDX,i,j,n1,n2,n3))


#if 0
!##############################################################################
! Kinetic Energy
!##############################################################################
#endif

#if 0
! Compute kinetic energy from conservative variables in 1D
#endif
#define KINETICENERGY1_1D(U,IDX)              (DCONST(0.5)*(POW(XMOMENTUM1_1D(U,IDX),2))/MYNEWLINE \
                                               DENSITY1_1D(U,IDX))
#define KINETICENERGY2_1D(U,IDX,i,n1,n2)      (DCONST(0.5)*(POW(XMOMENTUM2_1D(U,IDX,i,n1,n2),2))/MYNEWLINE \
                                               DENSITY2_1D(U,IDX,i,n1,n2))
#define KINETICENERGY3_1D(U,IDX,i,j,n1,n2,n3) (DCONST(0.5)*(POW(XMOMENTUM3_1D(U,IDX,i,j,n1,n2,n3),2))/MYNEWLINE \
                                               DENSITY3_1D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute kinetic energy from conservative variables in 2D
#endif
#define KINETICENERGY1_2D(U,IDX)              (DCONST(0.5)*(POW(XMOMENTUM1_2D(U,IDX),2)+MYNEWLINE \
                                                            POW(YMOMENTUM1_2D(U,IDX),2))/MYNEWLINE \
                                               DENSITY1_2D(U,IDX))
#define KINETICENERGY2_2D(U,IDX,i,n1,n2)      (DCONST(0.5)*(POW(XMOMENTUM2_2D(U,IDX,i,n1,n2),2)+MYNEWLINE \
                                                            POW(YMOMENTUM2_2D(U,IDX,i,n1,n2),2))/MYNEWLINE \
                                               DENSITY2_2D(U,IDX,i,n1,n2))
#define KINETICENERGY3_2D(U,IDX,i,j,n1,n2,n3) (DCONST(0.5)*(POW(XMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE \
                                                            POW(YMOMENTUM3_2D(U,IDX,i,j,n1,n2,n3),2))/MYNEWLINE \
                                               DENSITY3_2D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute kinetic energy from conservative variables in 3D
#endif
#define KINETICENERGY1_3D(U,IDX)              (DCONST(0.5)*(POW(XMOMENTUM1_3D(U,IDX),2)+MYNEWLINE \
                                                            POW(YMOMENTUM1_3D(U,IDX),2)+MYNEWLINE \
                                                            POW(ZMOMENTUM1_3D(U,IDX),2))/MYNEWLINE \
                                               DENSITY1_3D(U,IDX))
#define KINETICENERGY2_3D(U,IDX,i,n1,n2)      (DCONST(0.5)*(POW(XMOMENTUM2_3D(U,IDX,i,n1,n2),2)+MYNEWLINE \
                                                            POW(YMOMENTUM2_3D(U,IDX,i,n1,n2),2)+MYNEWLINE \
                                                            POW(ZMOMENTUM2_3D(U,IDX,i,n1,n2),2))/MYNEWLINE \
                                               DENSITY2_3D(U,IDX,i,n1,n2))
#define KINETICENERGY3_3D(U,IDX,i,j,n1,n2,n3) (DCONST(0.5)*(POW(XMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE \
                                                            POW(YMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3),2)+MYNEWLINE \
                                                            POW(ZMOMENTUM3_3D(U,IDX,i,j,n1,n2,n3),2))/MYNEWLINE \
                                               DENSITY3_3D(U,IDX,i,j,n1,n2,n3))


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
! Speed of sound
!##############################################################################
#endif

#if 0
! Compute speed of sound for a thermally perfect gas in 1D
#endif
#define SOUNDSPEED1_1D(U,IDX)              (sqrt((DCONST(THERMODYN_GAMMA))*MYNEWLINE \
                                                 PRESSURE1_1D(U,IDX)/MYNEWLINE \
                                                 DENSITY1_1D(U,IDX)))
#define SOUNDSPEED2_1D(U,IDX,i,n1,n2)      (sqrt((DCONST(THERMODYN_GAMMA))*MYNEWLINE \
                                                 PRESSURE2_1D(U,IDX,i,n1,n2)/MYNEWLINE \
                                                 DENSITY2_1D(U,IDX,i,n1,n2)))
#define SOUNDSPEED3_1D(U,IDX,i,j,n1,n2,n3) (sqrt((DCONST(THERMODYN_GAMMA))*MYNEWLINE \
                                                 PRESSURE3_1D(U,IDX,i,j,n1,n2,n3)/MYNEWLINE \
                                                 DENSITY3_1D(U,IDX,i,j,n1,n2,n3)))

#if 0
! Compute speed of sound for a thermally perfect gas in 2D
#endif
#define SOUNDSPEED1_2D(U,IDX)              (sqrt((DCONST(THERMODYN_GAMMA))*MYNEWLINE \
                                                 PRESSURE1_2D(U,IDX)/MYNEWLINE \
                                                 DENSITY1_2D(U,IDX)))
#define SOUNDSPEED2_2D(U,IDX,i,n1,n2)      (sqrt((DCONST(THERMODYN_GAMMA))*MYNEWLINE \
                                                 PRESSURE2_2D(U,IDX,i,n1,n2)/MYNEWLINE \
                                                 DENSITY2_2D(U,IDX,i,n1,n2)))
#define SOUNDSPEED3_2D(U,IDX,i,j,n1,n2,n3) (sqrt((DCONST(THERMODYN_GAMMA))*MYNEWLINE \
                                                 PRESSURE3_2D(U,IDX,i,j,n1,n2,n3)/MYNEWLINE \
                                                 DENSITY3_2D(U,IDX,i,j,n1,n2,n3)))

#if 0
! Compute speed of sound for a thermally perfect gas in 3D
#endif
#define SOUNDSPEED1_3D(U,IDX)              (sqrt((DCONST(THERMODYN_GAMMA))*MYNEWLINE \
                                                 PRESSURE1_3D(U,IDX)/MYNEWLINE \
                                                 DENSITY1_3D(U,IDX)))
#define SOUNDSPEED2_3D(U,IDX,i,n1,n2)      (sqrt((DCONST(THERMODYN_GAMMA))*MYNEWLINE \
                                                 PRESSURE2_3D(U,IDX,i,n1,n2)/MYNEWLINE \
                                                 DENSITY2_3D(U,IDX,i,n1,n2)))
#define SOUNDSPEED3_3D(U,IDX,i,j,n1,n2,n3) (sqrt((DCONST(THERMODYN_GAMMA))*MYNEWLINE \
                                                 PRESSURE3_3D(U,IDX,i,j,n1,n2,n3)/MYNEWLINE \
                                                 DENSITY3_3D(U,IDX,i,j,n1,n2,n3)))


#if 0
!##############################################################################
! Mach number
!##############################################################################
#endif

#if 0
! Compute the Mach number in 1D
#endif
#define MACHNUMBER1_1D(U,IDX)              (VELMAGNITUDE1_1D(U,IDX)/MYNEWLINE \
                                            SOUNDSPEED1_1D(U,IDX))
#define MACHNUMBER2_1D(U,IDX,i,n1,n2)      (VELMAGNITUDE2_1D(U,IDX,i,n1,n2)/MYNEWLINE \
                                            SOUNDSPEED2_1D(U,IDX,i,n1,n2))
#define MACHNUMBER3_1D(U,IDX,i,j,n1,n2,n3) (VELMAGNITUDE3_1D(U,IDX,i,j,n1,n2,n3)/MYNEWLINE \
                                            SOUNDSPEED3_1D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute the Mach number in 2D
#endif
#define MACHNUMBER1_2D(U,IDX)              (VELMAGNITUDE1_2D(U,IDX)/MYNEWLINE \
                                            SOUNDSPEED1_2D(U,IDX))
#define MACHNUMBER2_2D(U,IDX,i,n1,n2)      (VELMAGNITUDE2_2D(U,IDX,i,n1,n2)/MYNEWLINE \
                                            SOUNDSPEED2_2D(U,IDX,i,n1,n2))
#define MACHNUMBER3_2D(U,IDX,i,j,n1,n2,n3) (VELMAGNITUDE3_2D(U,IDX,i,j,n1,n2,n3)/MYNEWLINE \
                                            SOUNDSPEED3_2D(U,IDX,i,j,n1,n2,n3))

#if 0
! Compute the Mach number in 3D
#endif
#define MACHNUMBER1_3D(U,IDX)              (VELMAGNITUDE1_3D(U,IDX)/MYNEWLINE \
                                            SOUNDSPEED1_3D(U,IDX))
#define MACHNUMBER2_3D(U,IDX,i,n1,n2)      (VELMAGNITUDE2_3D(U,IDX,i,n1,n2)/MYNEWLINE \
                                            SOUNDSPEED2_3D(U,IDX,i,n1,n2))
#define MACHNUMBER3_3D(U,IDX,i,j,n1,n2,n3) (VELMAGNITUDE3_3D(U,IDX,i,j,n1,n2,n3)/MYNEWLINE \
                                            SOUNDSPEED3_3D(U,IDX,i,j,n1,n2,n3))
      
#endif


