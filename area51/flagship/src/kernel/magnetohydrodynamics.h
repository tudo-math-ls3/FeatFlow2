#if 0
! -*- mode: f90; -*-
!##############################################################################
!# ****************************************************************************
!# <name> magnetohydrodynamics </name>
!# ****************************************************************************
!#
!# <purpose>
!# This header file contains magnetohydrodynamics constants and pre-processor macros
!#
!# </purpose>
!##############################################################################
#endif

#ifndef _MAGNETOHYDRODYNAMICS_H_
#define _MAGNETOHYDRODYNAMICS_H_

#include "../flagship.h"


#if 0
!##############################################################################
! Compile-time constants
!##############################################################################
#endif

#if 0
! Universal gas constant
#endif
#define UNIVERSAL_GAS_CONSTANT   8.3145_DP


#if 0
!##############################################################################
! Equations of state
!##############################################################################
#endif

#if 0
! Constants for equations of state
#endif
#define EOS_PERFECT_GAS   1
#define EOS_TAMMANN       2
#define EOS_VANDERWAALS   3

#if 0
! Constants for types of gases
#endif
#ifdef PERFECT_GAS
#define THERMALLY_IDEAL_GAS
#define CALORICALLY_IDEAL_GAS
#define POLYTROPIC_GAS
#define HAS_HOMOGENEITY_PROPERTY
#endif

#if 0
! Set equation of state (if required)
#endif
#ifndef EOS
#ifdef PERFECT_GAS
#define EOS  EOS_PERFECT_GAS
#else
#error "Equation of state must be specified!"
#endif
#endif


#if 0
!##############################################################################
! Velocity field
!##############################################################################
#endif

#if 0
! Compute x-velocity from conservative variables in 1D, 2D, and 3D
#endif
#define X_VELOCITY_FROM_CONSVAR(u,nvar) \
  (X_MOMENTUM_FROM_CONSVAR(u,nvar)/ MYNEWLINE\
   DENSITY_FROM_CONSVAR(u,nvar))
#define X_VELOCITY_1T_FROM_CONSVAR(u,nvar,i1) \
  (X_MOMENTUM_1T_FROM_CONSVAR(u,nvar,i1)/ MYNEWLINE\
   DENSITY_1T_FROM_CONSVAR(u,nvar,i1))
#define X_VELOCITY_1L_FROM_CONSVAR(u,nvar,i1) \
  (X_MOMENTUM_1L_FROM_CONSVAR(u,nvar,i1)/ MYNEWLINE\
   DENSITY_1L_FROM_CONSVAR(u,nvar,i1))
#define X_VELOCITY_2T_FROM_CONSVAR(u,nvar,i1,i2) \
  (X_MOMENTUM_2T_FROM_CONSVAR(u,nvar,i1,i2)/ MYNEWLINE\
   DENSITY_2T_FROM_CONSVAR(u,nvar,i1,i2))
#define X_VELOCITY_2L_FROM_CONSVAR(u,nvar,i1,i2) \
  (X_MOMENTUM_2L_FROM_CONSVAR(u,nvar,i1,i2)/ MYNEWLINE\
   DENSITY_2L_FROM_CONSVAR(u,nvar,i1,i2))

#if 0
! Compute y-velocity from conservative variables in 2D, and 3D
#endif
#define Y_VELOCITY_FROM_CONSVAR(u,nvar) \
  (Y_MOMENTUM_FROM_CONSVAR(u,nvar)/ MYNEWLINE\
   DENSITY_FROM_CONSVAR(u,nvar))
#define Y_VELOCITY_1T_FROM_CONSVAR(u,nvar,i1) \
  (Y_MOMENTUM_1T_FROM_CONSVAR(u,nvar,i1)/ MYNEWLINE\
   DENSITY_1T_FROM_CONSVAR(u,nvar,i1))
#define Y_VELOCITY_1L_FROM_CONSVAR(u,nvar,i1) \
  (Y_MOMENTUM_1L_FROM_CONSVAR(u,nvar,i1)/ MYNEWLINE\
   DENSITY_1L_FROM_CONSVAR(u,nvar,i1))
#define Y_VELOCITY_2T_FROM_CONSVAR(u,nvar,i1,i2) \
  (Y_MOMENTUM_2T_FROM_CONSVAR(u,nvar,i1,i2)/ MYNEWLINE\
   DENSITY_2T_FROM_CONSVAR(u,nvar,i1,i2))
#define Y_VELOCITY_2L_FROM_CONSVAR(u,nvar,i1,i2) \
  (Y_MOMENTUM_2L_FROM_CONSVAR(u,nvar,i1,i2)/ MYNEWLINE\
   DENSITY_2L_FROM_CONSVAR(u,nvar,i1,i2))

#if 0
! Compute z-velocity from conservative variables in 3D
#endif
#define Z_VELOCITY_FROM_CONSVAR(u,nvar) \
  (Z_MOMENTUM_FROM_CONSVAR(u,nvar)/ MYNEWLINE\
   DENSITY_FROM_CONSVAR(u,nvar))
#define Z_VELOCITY_1T_FROM_CONSVAR(u,nvar,i1) \
  (Z_MOMENTUM_1T_FROM_CONSVAR(u,nvar,i1)/ MYNEWLINE\
   DENSITY_1T_FROM_CONSVAR(u,nvar,i1))
#define Z_VELOCITY_1L_FROM_CONSVAR(u,nvar,i1) \
  (Z_MOMENTUM_1L_FROM_CONSVAR(u,nvar,i1)/ MYNEWLINE\
   DENSITY_1L_FROM_CONSVAR(u,nvar,i1))
#define Z_VELOCITY_2T_FROM_CONSVAR(u,nvar,i1,i2) \
  (Z_MOMENTUM_2T_FROM_CONSVAR(u,nvar,i1,i2)/ MYNEWLINE\
   DENSITY_2T_FROM_CONSVAR(u,nvar,i1,i2))
#define Z_VELOCITY_2L_FROM_CONSVAR(u,nvar,i1,i2) \
  (Z_MOMENTUM_2L_FROM_CONSVAR(u,nvar,i1,i2)/ MYNEWLINE\
   DENSITY_2L_FROM_CONSVAR(u,nvar,i1,i2))


#if 0
!##############################################################################
! Magnitude of velocity field
!
! There are always three momentum/velocity components.
!##############################################################################
#endif

#if 0
  ! Compute magnitude of velocity field in 1D, 2D, and 3D
#endif
#define VELOCITY_MAGNITUDE_FROM_CONSVAR(u,nvar) \
  (sqrt(X_VELOCITY_FROM_CONSVAR(u,nvar)**2+ MYNEWLINE\
        Y_VELOCITY_FROM_CONSVAR(u,nvar)**2+ MYNEWLINE\
        Z_VELOCITY_FROM_CONSVAR(u,nvar)**2))
#define VELOCITY_MAGNITUDE_1T_FROM_CONSVAR(u,nvar,i1) \
  (sqrt(X_VELOCITY_1T_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Y_VELOCITY_1T_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Z_VELOCITY_1T_FROM_CONSVAR(u,nvar,i1)**2))
#define VELOCITY_MAGNITUDE_1L_FROM_CONSVAR(u,nvar,i1) \
  (sqrt(X_VELOCITY_1L_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Y_VELOCITY_1L_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Z_VELOCITY_1L_FROM_CONSVAR(u,nvar,i1)**2))
#define VELOCITY_MAGNITUDE_2T_FROM_CONSVAR(u,nvar,i1,i2) \
  (sqrt(X_VELOCITY_2T_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Y_VELOCITY_2T_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Z_VELOCITY_2T_FROM_CONSVAR(u,nvar,i1,i2)**2))
#define VELOCITY_MAGNITUDE_2L_FROM_CONSVAR(u,nvar,i1,i2) \
  (sqrt(X_VELOCITY_2L_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Y_VELOCITY_2L_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Z_VELOCITY_2L_FROM_CONSVAR(u,nvar,i1,i2)**2))


#if 0
!##############################################################################
! Magnitude of magnetic field
!
! There are two components in 1D and three components in 2D and 3D
!##############################################################################
#endif

#if 0
! Compute magnitude of magnetic field in 1D
#endif
#define MAGNETICFIELD_MAGNITUDE_FROM_CONSVAR_1D(u,nvar) \
  (sqrt(X_MAGNETICFIELD_FROM_CONSVAR_1D(u,nvar)**2+ MYNEWLINE\
        Y_MAGNETICFIELD_FROM_CONSVAR_1D(u,nvar)**2+ MYNEWLINE\
        Z_MAGNETICFIELD_FROM_CONSVAR_1D(u,nvar)**2))
#define MAGNETICFIELD_MAGNITUDE_1T_FROM_CONSVAR_1D(u,nvar,i1) \
  (sqrt(X_MAGNETICFIELD_1T_FROM_CONSVAR_1D(u,nvar,i1)**2+ MYNEWLINE\
        Y_MAGNETICFIELD_1T_FROM_CONSVAR_1D(u,nvar,i1)**2+ MYNEWLINE\
        Z_MAGNETICFIELD_1T_FROM_CONSVAR_1D(u,nvar,i1)**2))
#define MAGNETICFIELD_MAGNITUDE_1L_FROM_CONSVAR_1D(u,nvar,i1) \
  (sqrt(X_MAGNETICFIELD_1L_FROM_CONSVAR_1D(u,nvar,i1)**2+ MYNEWLINE\
        Y_MAGNETICFIELD_1L_FROM_CONSVAR_1D(u,nvar,i1)**2+ MYNEWLINE\
        Z_MAGNETICFIELD_1L_FROM_CONSVAR_1D(u,nvar,i1)**2))
#define MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR_1D(u,nvar,i1,i2) \
  (sqrt(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(u,nvar,i1,i2)**2+ MYNEWLINE\
        Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(u,nvar,i1,i2)**2+ MYNEWLINE\
        Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(u,nvar,i1,i2)**2))
#define MAGNETICFIELD_MAGNITUDE_2L_FROM_CONSVAR_1D(u,nvar,i1,i2) \
  (sqrt(X_MAGNETICFIELD_2L_FROM_CONSVAR_1D(u,nvar,i1,i2)**2+ MYNEWLINE\
        Y_MAGNETICFIELD_2L_FROM_CONSVAR_1D(u,nvar,i1,i2)**2+ MYNEWLINE\
        Z_MAGNETICFIELD_2L_FROM_CONSVAR_1D(u,nvar,i1,i2)**2))

#if 0
! Compute magnitude of magnetic field in 2D and 3D
#endif
#define MAGNETICFIELD_MAGNITUDE_FROM_CONSVAR(u,nvar) \
  (sqrt(X_MAGNETICFIELD_FROM_CONSVAR(u,nvar)**2+ MYNEWLINE\
        Y_MAGNETICFIELD_FROM_CONSVAR(u,nvar)**2+ MYNEWLINE\
        Z_MAGNETICFIELD_FROM_CONSVAR(u,nvar)**2))
#define MAGNETICFIELD_MAGNITUDE_1T_FROM_CONSVAR(u,nvar,i1) \
  (sqrt(X_MAGNETICFIELD_1T_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Y_MAGNETICFIELD_1T_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Z_MAGNETICFIELD_1T_FROM_CONSVAR(u,nvar,i1)**2))
#define MAGNETICFIELD_MAGNITUDE_1L_FROM_CONSVAR(u,nvar,i1) \
  (sqrt(X_MAGNETICFIELD_1L_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Y_MAGNETICFIELD_1L_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Z_MAGNETICFIELD_1L_FROM_CONSVAR(u,nvar,i1)**2))
#define MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(u,nvar,i1,i2) \
  (sqrt(X_MAGNETICFIELD_2T_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Y_MAGNETICFIELD_2T_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Z_MAGNETICFIELD_2T_FROM_CONSVAR(u,nvar,i1,i2)**2))
#define MAGNETICFIELD_MAGNITUDE_2L_FROM_CONSVAR(u,nvar,i1,i2) \
  (sqrt(X_MAGNETICFIELD_2L_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Y_MAGNETICFIELD_2L_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Z_MAGNETICFIELD_2L_FROM_CONSVAR(u,nvar,i1,i2)**2))


#if 0
!##############################################################################
! Pressure
!
! The total pressure involves the magnitic field so that we
! have to distinguish between 1D and multiple dimensions
!##############################################################################
#endif

#if (EOS == EOS_PERFECT_GAS)

#if 0
! Use equation of state for perfect gases
#endif

#if 0
  ! Compute pressure from conservative variables in 1D
#endif
#define PRESSURE_FROM_CONSVAR_1D(u,nvar) \
  ((GAMMA-1.0)*INTERNAL_ENERGY_FROM_CONSVAR_1D(u,nvar))
#define PRESSURE_1T_FROM_CONSVAR_1D(u,nvar,i1) \
  ((GAMMA-1.0)*INTERNAL_ENERGY_1T_FROM_CONSVAR_1D(u,nvar,i1))
#define PRESSURE_1L_FROM_CONSVAR_1D(u,nvar,i1) \
  ((GAMMA-1.0)*INTERNAL_ENERGY_1L_FROM_CONSVAR_1D(u,nvar,i1))
#define PRESSURE_2T_FROM_CONSVAR_1D(u,nvar,i1,i2) \
  ((GAMMA-1.0)*INTERNAL_ENERGY_2T_FROM_CONSVAR_1D(u,nvar,i1,i2))
#define PRESSURE_2L_FROM_CONSVAR_1D(u,nvar,i1,i2) \
  ((GAMMA-1.0)*INTERNAL_ENERGY_2L_FROM_CONSVAR_1D(u,nvar,i1,i2))

#if 0
  ! Compute pressure from conservative variables in 2D and 3D
#endif
#define PRESSURE_FROM_CONSVAR(u,nvar) \
  ((GAMMA-1.0)*INTERNAL_ENERGY_FROM_CONSVAR(u,nvar))
#define PRESSURE_1T_FROM_CONSVAR(u,nvar,i1) \
  ((GAMMA-1.0)*INTERNAL_ENERGY_1T_FROM_CONSVAR(u,nvar,i1))
#define PRESSURE_1L_FROM_CONSVAR(u,nvar,i1) \
  ((GAMMA-1.0)*INTERNAL_ENERGY_1L_FROM_CONSVAR(u,nvar,i1))
#define PRESSURE_2T_FROM_CONSVAR(u,nvar,i1,i2) \
  ((GAMMA-1.0)*INTERNAL_ENERGY_2T_FROM_CONSVAR(u,nvar,i1,i2))
#define PRESSURE_2L_FROM_CONSVAR(u,nvar,i1,i2) \
  ((GAMMA-1.0)*INTERNAL_ENERGY_2L_FROM_CONSVAR(u,nvar,i1,i2))

#else
#error "Invalid equation of state"
#endif


#if 0
!##############################################################################
! Magnetic pressure
!
! The magnetic pressure involves the magnitic field so that we
! have to distinguish between 1D and multiple dimensions
!##############################################################################
#endif

#define MAGNETIC_PRESSURE_FROM_CONSVAR_1D    MAGNETIC_ENERGY_FROM_CONSVAR_1D
#define MAGNETIC_PRESSURE_1T_FROM_CONSVAR_1D MAGNETIC_ENERGY_1T_FROM_CONSVAR_1D
#define MAGNETIC_PRESSURE_1L_FROM_CONSVAR_1D MAGNETIC_ENERGY_1L_FROM_CONSVAR_1D
#define MAGNETIC_PRESSURE_2T_FROM_CONSVAR_1D MAGNETIC_ENERGY_2T_FROM_CONSVAR_1D
#define MAGNETIC_PRESSURE_2L_FROM_CONSVAR_1D MAGNETIC_ENERGY_2L_FROM_CONSVAR_1D

#define MAGNETIC_PRESSURE_FROM_CONSVAR    MAGNETIC_ENERGY_FROM_CONSVAR
#define MAGNETIC_PRESSURE_1T_FROM_CONSVAR MAGNETIC_ENERGY_1T_FROM_CONSVAR
#define MAGNETIC_PRESSURE_1L_FROM_CONSVAR MAGNETIC_ENERGY_1L_FROM_CONSVAR
#define MAGNETIC_PRESSURE_2T_FROM_CONSVAR MAGNETIC_ENERGY_2T_FROM_CONSVAR
#define MAGNETIC_PRESSURE_2L_FROM_CONSVAR MAGNETIC_ENERGY_2L_FROM_CONSVAR


#if 0
!##############################################################################
! Total pressure
!
! The total pressure involves the magnitic field so that we
! have to distinguish between 1D and multiple dimensions
!##############################################################################
#endif

#if 0
! Compute total pressure from conservative variables in 1D
#endif
#define TOTAL_PRESSURE_FROM_CONSVAR_1D(u,nvar) \
  (PRESSURE_FROM_CONSVAR_1D(u,nvar)+MYNEWLINE\
   MAGNETIC_PRESSURE_FROM_CONSVAR_1D(u,nvar))
#define TOTAL_PRESSURE_1T_FROM_CONSVAR_1D(u,nvar,i1) \
  (PRESSURE_1T_FROM_CONSVAR_1D(u,nvar,i1)+MYNEWLINE\
   MAGNETIC_PRESSURE_1T_FROM_CONSVAR_1D(u,nvar,i1))
#define TOTAL_PRESSURE_1L_FROM_CONSVAR_1D(u,nvar,i1) \
  (PRESSURE_1L_FROM_CONSVAR_1D(u,nvar,i1)+MYNEWLINE\
   MAGNETIC_PRESSURE_1L_FROM_CONSVAR_1D(u,nvar,i1))
#define TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(u,nvar,i1,i2) \
  (PRESSURE_2T_FROM_CONSVAR_1D(u,nvar,i1,i2)+MYNEWLINE\
   MAGNETIC_PRESSURE_2T_FROM_CONSVAR_1D(u,nvar,i1,i2))
#define TOTAL_PRESSURE_2L_FROM_CONSVAR_1D(u,nvar,i1,i2) \
  (PRESSURE_2L_FROM_CONSVAR_1D(u,nvar,i1,i2)+MYNEWLINE\
   MAGNETIC_PRESSURE_2L_FROM_CONSVAR_1D(u,nvar,i1,i2))

#if 0
! Compute total pressure from conservative variables in 2D and 3D
#endif
#define TOTAL_PRESSURE_FROM_CONSVAR(u,nvar) \
  (PRESSURE_FROM_CONSVAR(u,nvar)+MYNEWLINE\
   MAGNETIC_PRESSURE_FROM_CONSVAR(u,nvar))
#define TOTAL_PRESSURE_1T_FROM_CONSVAR(u,nvar,i1) \
  (PRESSURE_1T_FROM_CONSVAR(u,nvar,i1)+MYNEWLINE\
   MAGNETIC_PRESSURE_1T_FROM_CONSVAR(u,nvar,i1))
#define TOTAL_PRESSURE_1L_FROM_CONSVAR(u,nvar,i1) \
  (PRESSURE_1L_FROM_CONSVAR(u,nvar,i1)+MYNEWLINE\
   MAGNETIC_PRESSURE_1L_FROM_CONSVAR(u,nvar,i1))
#define TOTAL_PRESSURE_2T_FROM_CONSVAR(u,nvar,i1,i2) \
  (PRESSURE_2T_FROM_CONSVAR(u,nvar,i1,i2)+MYNEWLINE\
   MAGNETIC_PRESSURE_2T_FROM_CONSVAR(u,nvar,i1,i2))
#define TOTAL_PRESSURE_2L_FROM_CONSVAR(u,nvar,i1,i2) \
  (PRESSURE_2L_FROM_CONSVAR(u,nvar,i1,i2)+MYNEWLINE\
   MAGNETIC_PRESSURE_2L_FROM_CONSVAR(u,nvar,i1,i2))


#if 0
!##############################################################################
! Kinetic Energy
!
! There are always three momentum/velocity components.
!##############################################################################
#endif

#if 0
  ! Compute kinetic energy from conservative variables in 1D, 2D, and 3D
#endif
#define KINETIC_ENERGY_FROM_CONSVAR(u,nvar) \
  (0.5* (X_MOMENTUM_FROM_CONSVAR(u,nvar)**2+ MYNEWLINE\
         Y_MOMENTUM_FROM_CONSVAR(u,nvar)**2+ MYNEWLINE\
         Z_MOMENTUM_FROM_CONSVAR(u,nvar)**2)/ MYNEWLINE\
         DENSITY_FROM_CONSVAR(u,nvar))
#define KINETIC_ENERGY_1T_FROM_CONSVAR(u,nvar,i1) \
  (0.5* (X_MOMENTUM_1T_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
         Y_MOMENTUM_1T_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
         Z_MOMENTUM_1T_FROM_CONSVAR(u,nvar,i1)**2)/ MYNEWLINE\
         DENSITY_1T_FROM_CONSVAR(u,nvar,i1))
#define KINETIC_ENERGY_1L_FROM_CONSVAR(u,nvar,i1) \
  (0.5* (X_MOMENTUM_1L_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
         Y_MOMENTUM_1L_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
         Z_MOMENTUM_1L_FROM_CONSVAR(u,nvar,i1)**2)/ MYNEWLINE\
         DENSITY_1L_FROM_CONSVAR(u,nvar,i1))
#define KINETIC_ENERGY_2T_FROM_CONSVAR(u,nvar,i1,i2) \
  (0.5* (X_MOMENTUM_2T_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
         Y_MOMENTUM_2T_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
         Z_MOMENTUM_2T_FROM_CONSVAR(u,nvar,i1,2)**2)/ MYNEWLINE\
         DENSITY_2T_FROM_CONSVAR(u,nvar,i1,i2))
#define KINETIC_ENERGY_2L_FROM_CONSVAR(u,nvar,i1,i2) \
  (0.5* (X_MOMENTUM_2L_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
         Y_MOMENTUM_2L_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
         Z_MOMENTUM_2L_FROM_CONSVAR(u,nvar,i1,i2)**2)/ MYNEWLINE\
         DENSITY_2L_FROM_CONSVAR(u,nvar,i1,i2))


#if 0
!##############################################################################
! Internal Energy
!
! The internal energy involves the magnitic field so that we
! have to distinguish between 1D and multiple dimensions
!##############################################################################
#endif

#if 0
! Compute internal energy from conservative variables in 1D
#endif
#define INTERNAL_ENERGY_FROM_CONSVAR_1D(u,nvar) \
  (TOTAL_ENERGY_FROM_CONSVAR(u,nvar)- MYNEWLINE\
   KINETIC_ENERGY_FROM_CONSVAR(u,nvar)- MYNEWLINE\
   MAGNETIC_ENERGY_FROM_CONSVAR_1D(u,nvar))
#define INTERNAL_ENERGY_1T_FROM_CONSVAR_1D(u,nvar,i1) \
  (TOTAL_ENERGY_1T_FROM_CONSVAR(u,nvar,i1)- MYNEWLINE\
   KINETIC_ENERGY_1T_FROM_CONSVAR(u,nvar,i1)- MYNEWLINE\
   MAGNETIC_ENERGY_1T_FROM_CONSVAR_1D(u,nvar,i1))
#define INTERNAL_ENERGY_1L_FROM_CONSVAR_1D(u,nvar,i1) \
  (TOTAL_ENERGY_1L_FROM_CONSVAR(u,nvar,i1)- MYNEWLINE\
   KINETIC_ENERGY_1L_FROM_CONSVAR(u,nvar,i1)- MYNEWLINE\
   MAGNETIC_ENERGY_1L_FROM_CONSVAR_1D(u,nvar,i1))
#define INTERNAL_ENERGY_2T_FROM_CONSVAR_1D(u,nvar,i1,i2) \
  (TOTAL_ENERGY_2T_FROM_CONSVAR(u,nvar,i1,i2)- MYNEWLINE\
   KINETIC_ENERGY_2T_FROM_CONSVAR(u,nvar,i1,i2)- MYNEWLINE\
   MAGNETIC_ENERGY_2T_FROM_CONSVAR_1D(u,nvar,i1,i2))
#define INTERNAL_ENERGY_2L_FROM_CONSVAR_1D(u,nvar,i1,i2) \
  (TOTAL_ENERGY_2L_FROM_CONSVAR(u,nvar,i1,i2)- MYNEWLINE\
   KINETIC_ENERGY_2L_FROM_CONSVAR(u,nvar,i1,i2)- MYNEWLINE\
   MAGNETIC_ENERGY_2L_FROM_CONSVAR_1D(u,nvar,i1,i2))

#if 0
! Compute internal energy from conservative variables in 2D and 3D
#endif
#define INTERNAL_ENERGY_FROM_CONSVAR(u,nvar) \
  (TOTAL_ENERGY_FROM_CONSVAR(u,nvar)- MYNEWLINE\
   KINETIC_ENERGY_FROM_CONSVAR(u,nvar)- MYNEWLINE\
   MAGNETIC_ENERGY_FROM_CONSVAR(u,nvar))
#define INTERNAL_ENERGY_1T_FROM_CONSVAR(u,nvar,i1) \
  (TOTAL_ENERGY_1T_FROM_CONSVAR(u,nvar,i1)- MYNEWLINE\
   KINETIC_ENERGY_1T_FROM_CONSVAR(u,nvar,i1)- MYNEWLINE\
   MAGNETIC_ENERGY_1T_FROM_CONSVAR(u,nvar,i1))
#define INTERNAL_ENERGY_1L_FROM_CONSVAR(u,nvar,i1) \
  (TOTAL_ENERGY_1L_FROM_CONSVAR(u,nvar,i1)- MYNEWLINE\
   KINETIC_ENERGY_1L_FROM_CONSVAR(u,nvar,i1)- MYNEWLINE\
   MAGNETIC_ENERGY_1L_FROM_CONSVAR(u,nvar,i1))
#define INTERNAL_ENERGY_2T_FROM_CONSVAR(u,nvar,i1,i2) \
  (TOTAL_ENERGY_2T_FROM_CONSVAR(u,nvar,i1,i2)- MYNEWLINE\
   KINETIC_ENERGY_2T_FROM_CONSVAR(u,nvar,i1,i2)- MYNEWLINE\
   MAGNETIC_ENERGY_2T_FROM_CONSVAR(u,nvar,i1,i2))
#define INTERNAL_ENERGY_2L_FROM_CONSVAR(u,nvar,i1,i2) \
  (TOTAL_ENERGY_2L_FROM_CONSVAR(u,nvar,i1,i2)- MYNEWLINE\
   KINETIC_ENERGY_2L_FROM_CONSVAR(u,nvar,i1,i2)- MYNEWLINE\
   MAGNETIC_ENERGY_2L_FROM_CONSVAR(u,nvar,i1,i2))


#if 0
!##############################################################################
! Magnetic energy
!
! There are two components in 1D and three components in 2D and 3D
!##############################################################################
#endif

#if 0
! Compute magnetic energy in 1D
#endif
#define MAGNETIC_ENERGY_FROM_CONSVAR_1D(u,nvar) \
  ((X_MAGNETICFIELD_FROM_CONSVAR_1D(u,nvar)**2+ MYNEWLINE\
    Y_MAGNETICFIELD_FROM_CONSVAR_1D(u,nvar)**2+ MYNEWLINE\
    Z_MAGNETICFIELD_FROM_CONSVAR_1D(u,nvar)**2)/(2.0*VACUUM_PERM))
#define MAGNETIC_ENERGY_1T_FROM_CONSVAR_1D(u,nvar,i1) \
  ((X_MAGNETICFIELD_1T_FROM_CONSVAR_1D(u,nvar,i1)**2+ MYNEWLINE\
    Y_MAGNETICFIELD_1T_FROM_CONSVAR_1D(u,nvar,i1)**2+ MYNEWLINE\
    Z_MAGNETICFIELD_1T_FROM_CONSVAR_1D(u,nvar,i1)**2)/(2.0*VACUUM_PERM))
#define MAGNETIC_ENERGY_1L_FROM_CONSVAR_1D(u,nvar,i1) \
  ((X_MAGNETICFIELD_1L_FROM_CONSVAR_1D(u,nvar,i1)**2+ MYNEWLINE\
    Y_MAGNETICFIELD_1L_FROM_CONSVAR_1D(u,nvar,i1)**2+ MYNEWLINE\
    Z_MAGNETICFIELD_1L_FROM_CONSVAR_1D(u,nvar,i1)**2)/(2.0*VACUUM_PERM))
#define MAGNETIC_ENERGY_2T_FROM_CONSVAR_1D(u,nvar,i1,i2) \
  ((X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(u,nvar,i1,i2)**2+ MYNEWLINE\
    Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(u,nvar,i1,i2)**2+ MYNEWLINE\
    Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(u,nvar,i1,i2)**2)/(2.0*VACUUM_PERM))
#define MAGNETIC_ENERGY_2L_FROM_CONSVAR_1D(u,nvar,i1,i2) \
  ((X_MAGNETICFIELD_2L_FROM_CONSVAR_1D(u,nvar,i1,i2)**2+ MYNEWLINE\
    Y_MAGNETICFIELD_2L_FROM_CONSVAR_1D(u,nvar,i1,i2)**2+ MYNEWLINE\
    Z_MAGNETICFIELD_2L_FROM_CONSVAR_1D(u,nvar,i1,i2)**2)/(2.0*VACUUM_PERM))

#if 0
! Compute magnitude energy in 2D and 3D
#endif
#define MAGNETIC_ENERGY_FROM_CONSVAR(u,nvar) \
  ((X_MAGNETICFIELD_FROM_CONSVAR(u,nvar)**2+ MYNEWLINE\
    Y_MAGNETICFIELD_FROM_CONSVAR(u,nvar)**2+ MYNEWLINE\
    Z_MAGNETICFIELD_FROM_CONSVAR(u,nvar)**2)/(2.0*VACUUM_PERM))
#define MAGNETIC_ENERGY_1T_FROM_CONSVAR(u,nvar,i1) \
  ((X_MAGNETICFIELD_1T_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
    Y_MAGNETICFIELD_1T_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
    Z_MAGNETICFIELD_1T_FROM_CONSVAR(u,nvar,i1)**2)/(2.0*VACUUM_PERM))
#define MAGNETIC_ENERGY_1L_FROM_CONSVAR(u,nvar,i1) \
  ((X_MAGNETICFIELD_1L_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
    Y_MAGNETICFIELD_1L_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
    Z_MAGNETICFIELD_1L_FROM_CONSVAR(u,nvar,i1)**2)/(2.0*VACUUM_PERM))
#define MAGNETIC_ENERGY_2T_FROM_CONSVAR(u,nvar,i1,i2) \
  ((X_MAGNETICFIELD_2T_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
    Y_MAGNETICFIELD_2T_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
    Z_MAGNETICFIELD_2T_FROM_CONSVAR(u,nvar,i1,i2)**2)/(2.0*VACUUM_PERM))
#define MAGNETIC_ENERGY_2L_FROM_CONSVAR(u,nvar,i1,i2) \
  ((X_MAGNETICFIELD_2L_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
    Y_MAGNETICFIELD_2L_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
    Z_MAGNETICFIELD_2L_FROM_CONSVAR(u,nvar,i1,i2)**2)/(2.0*VACUUM_PERM))


#if 0
!##############################################################################
! Speed of sound
!
! There are two components in 1D and three components in 2D and 3D
!##############################################################################
#endif

#ifdef THERMALLY_IDEAL_GAS

#if 0
! Compute speed of sound for a thermally perfect gas in 1D
#endif
#define SPEED_OF_SOUND_FROM_CONSVAR_1D(u,nvar) \
  (sqrt(GAMMA* MYNEWLINE\
   PRESSURE_FROM_CONSVAR_1D(u,nvar)/ MYNEWLINE\
   DENSITY_FROM_CONSVAR(u,nvar)))
#define SPEED_OF_SOUND_1T_FROM_CONSVAR_1D(u,nvar,i1) \
  (sqrt(GAMMA* MYNEWLINE\
   PRESSURE_1T_FROM_CONSVAR_1D(u,nvar,i1)/ MYNEWLINE\
   DENSITY_1T_FROM_CONSVAR(u,nvar,i1)))
#define SPEED_OF_SOUND_1L_FROM_CONSVAR_1D(u,nvar,i1) \
  (sqrt(GAMMA* MYNEWLINE\
   PRESSURE_1L_FROM_CONSVAR_1D(u,nvar,i1)/ MYNEWLINE\
   DENSITY_1L_FROM_CONSVAR(u,nvar,i1)))
#define SPEED_OF_SOUND_2T_FROM_CONSVAR_1D(u,nvar,i1,i2) \
  (sqrt(GAMMA* MYNEWLINE\
   PRESSURE_2T_FROM_CONSVAR_1D(u,nvar,i1,i2)/ MYNEWLINE\
   DENSITY_2T_FROM_CONSVAR(u,nvar,i1,i2)))
#define SPEED_OF_SOUND_2L_FROM_CONSVAR_1D(u,nvar,i1,i2) \
  (sqrt(GAMMA* MYNEWLINE\
   PRESSURE_2L_FROM_CONSVAR_1D(u,nvar,i1,i2)/ MYNEWLINE\
   DENSITY_2L_FROM_CONSVAR(u,nvar,i1,i2)))

#if 0
! Compute speed of sound for a thermally perfect gas in 2D and 3D
#endif
#define SPEED_OF_SOUND_FROM_CONSVAR(u,nvar) \
  (sqrt(GAMMA* MYNEWLINE\
   PRESSURE_FROM_CONSVAR(u,nvar)/ MYNEWLINE\
   DENSITY_FROM_CONSVAR(u,nvar)))
#define SPEED_OF_SOUND_1T_FROM_CONSVAR(u,nvar,i1) \
  (sqrt(GAMMA* MYNEWLINE\
   PRESSURE_1T_FROM_CONSVAR(u,nvar,i1)/ MYNEWLINE\
   DENSITY_1T_FROM_CONSVAR(u,nvar,i1)))
#define SPEED_OF_SOUND_1L_FROM_CONSVAR(u,nvar,i1) \
  (sqrt(GAMMA* MYNEWLINE\
   PRESSURE_1L_FROM_CONSVAR(u,nvar,i1)/ MYNEWLINE\
   DENSITY_1L_FROM_CONSVAR(u,nvar,i1)))
#define SPEED_OF_SOUND_2T_FROM_CONSVAR(u,nvar,i1,i2) \
  (sqrt(GAMMA* MYNEWLINE\
   PRESSURE_2T_FROM_CONSVAR(u,nvar,i1,i2)/ MYNEWLINE\
   DENSITY_2T_FROM_CONSVAR(u,nvar,i1,i2)))
#define SPEED_OF_SOUND_2L_FROM_CONSVAR(u,nvar,i1,i2) \
  (sqrt(GAMMA* MYNEWLINE\
   PRESSURE_2L_FROM_CONSVAR(u,nvar,i1,i2)/ MYNEWLINE\
   DENSITY_2L_FROM_CONSVAR(u,nvar,i1,i2)))

#endif

#if 0
!##############################################################################
! Mach number
!
! There are two components in 1D and three components in 2D and 3D
!##############################################################################
#endif

#if 0
! Compute the Mach number in 1D
#endif
#define MACH_NUMBER_FROM_CONSVAR_1D(u,nvar) \
  (sqrt(X_VELOCITY_FROM_CONSVAR(u,nvar)**2+ MYNEWLINE\
        Y_VELOCITY_FROM_CONSVAR(u,nvar)**2+ MYNEWLINE\
        Z_VELOCITY_FROM_CONSVAR(u,nvar)**2)/ MYNEWLINE\
        SPEED_OF_SOUND_FROM_CONSVAR_1D(u,nvar))
#define MACH_NUMBER_1T_FROM_CONSVAR_1D(u,nvar,i1) \
  (sqrt(X_VELOCITY_1T_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Y_VELOCITY_1T_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Z_VELOCITY_1T_FROM_CONSVAR(u,nvar,i1)**2)/ MYNEWLINE\
        SPEED_OF_SOUND_1T_FROM_CONSVAR_1D(u,nvar,i1))
#define MACH_NUMBER_1L_FROM_CONSVAR_1D(u,nvar,i1) \
  (sqrt(X_VELOCITY_1L_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Y_VELOCITY_1L_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Z_VELOCITY_1L_FROM_CONSVAR(u,nvar,i1)**2)/ MYNEWLINE\
        SPEED_OF_SOUND_1L_FROM_CONSVAR_1D(u,nvar,i1))
#define MACH_NUMBER_2T_FROM_CONSVAR_1D(u,nvar,i1,i2) \
  (sqrt(X_VELOCITY_2T_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Y_VELOCITY_2T_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Z_VELOCITY_2T_FROM_CONSVAR(u,nvar,i1,i2)**2)/ MYNEWLINE\
        SPEED_OF_SOUND_2T_FROM_CONSVAR_1D(u,nvar,i1,i2))
#define MACH_NUMBER_2L_FROM_CONSVAR_1D(u,nvar,i1,i2) \
  (sqrt(X_VELOCITY_2L_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Y_VELOCITY_2L_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Z_VELOCITY_2L_FROM_CONSVAR(u,nvar,i1,i2)**2)/ MYNEWLINE\
        SPEED_OF_SOUND_2L_FROM_CONSVAR_1D(u,nvar,i1,i2))

#if 0
! Compute the Mach number in 2D and 3D
#endif
#define MACH_NUMBER_FROM_CONSVAR(u,nvar) \
  (sqrt(X_VELOCITY_FROM_CONSVAR(u,nvar)**2+ MYNEWLINE\
        Y_VELOCITY_FROM_CONSVAR(u,nvar)**2+ MYNEWLINE\
        Z_VELOCITY_FROM_CONSVAR(u,nvar)**2)/ MYNEWLINE\
        SPEED_OF_SOUND_FROM_CONSVAR(u,nvar))
#define MACH_NUMBER_1T_FROM_CONSVAR(u,nvar,i1) \
  (sqrt(X_VELOCITY_1T_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Y_VELOCITY_1T_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Z_VELOCITY_1T_FROM_CONSVAR(u,nvar,i1)**2)/ MYNEWLINE\
        SPEED_OF_SOUND_1T_FROM_CONSVAR(u,nvar,i1))
#define MACH_NUMBER_1L_FROM_CONSVAR(u,nvar,i1) \
  (sqrt(X_VELOCITY_1L_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Y_VELOCITY_1L_FROM_CONSVAR(u,nvar,i1)**2+ MYNEWLINE\
        Z_VELOCITY_1L_FROM_CONSVAR(u,nvar,i1)**2)/ MYNEWLINE\
        SPEED_OF_SOUND_1L_FROM_CONSVAR(u,nvar,i1))
#define MACH_NUMBER_2T_FROM_CONSVAR(u,nvar,i1,i2) \
  (sqrt(X_VELOCITY_2T_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Y_VELOCITY_2T_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Z_VELOCITY_2T_FROM_CONSVAR(u,nvar,i1,i2)**2)/ MYNEWLINE\
        SPEED_OF_SOUND_2T_FROM_CONSVAR(u,nvar,i1,i2))
#define MACH_NUMBER_2L_FROM_CONSVAR(u,nvar,i1,i2) \
  (sqrt(X_VELOCITY_2L_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Y_VELOCITY_2L_FROM_CONSVAR(u,nvar,i1,i2)**2+ MYNEWLINE\
        Z_VELOCITY_2L_FROM_CONSVAR(u,nvar,i1,i2)**2)/ MYNEWLINE\
        SPEED_OF_SOUND_2L_FROM_CONSVAR(u,nvar,i1,i2))

#endif
