#if 0
! -*- mode: f90; -*-
!##############################################################################
!# ****************************************************************************
!# <name> mhd </name>
!# ****************************************************************************
!#
!# <purpose>
!# This is the main header file for the compressible MHD model
!#
!# </purpose>
!##############################################################################
#endif

#ifndef _MHD_H_
#define _MHD_H_

#include "../../flagship.h"

#if 0
! Equation of state for perfect gas
#endif
#define PERFECT_GAS

#if 0
! Vacuum permittivity
#endif
#define VACUUM_PERM 1.0

#if 0
! Ratio of specific heats for air at sea-level
#endif
#define GAMMA 1.4_DP

#if 0
! In one dimension, the x-component of the magnetic field is constant.
#endif
#define X_MAGNETICFIELD_CONSTANT_1D 3.0/4.0

#if 0
!##############################################################################
! Conservative variables
!##############################################################################
#endif

#if 0
! Compute density from conservative variables in 1D,2D, and 3D
#endif
#define DENSITY_FROM_CONSVAR(u,nvar) (u(1))
#define DENSITY_1T_FROM_CONSVAR(u,nvar,i1) (u(1,i1))
#define DENSITY_1L_FROM_CONSVAR(u,nvar,i1) (u(i1,1))
#define DENSITY_2T_FROM_CONSVAR(u,nvar,i1,i2) (u(1,i1,i2))
#define DENSITY_2L_FROM_CONSVAR(u,nvar,i1,i2) (u(i1,i2,1))

#if 0
! Compute x-momentum from conservative variables in 1D, 2D, and 3D
#endif
#define X_MOMENTUM_FROM_CONSVAR(u,nvar) (u(2))
#define X_MOMENTUM_1T_FROM_CONSVAR(u,nvar,i1) (u(2,i1))
#define X_MOMENTUM_1L_FROM_CONSVAR(u,nvar,i1) (u(i1,2))
#define X_MOMENTUM_2T_FROM_CONSVAR(u,nvar,i1,i2) (u(2,i1,i2))
#define X_MOMENTUM_2L_FROM_CONSVAR(u,nvar,i1,i2) (u(i1,i2,2))

#if 0
! Compute y-momentum from conservative variables in 2D, and 3D
#endif
#define Y_MOMENTUM_FROM_CONSVAR(u,nvar) (u(3))
#define Y_MOMENTUM_1T_FROM_CONSVAR(u,nvar,i1) (u(3,i1))
#define Y_MOMENTUM_1L_FROM_CONSVAR(u,nvar,i1) (u(i1,3))
#define Y_MOMENTUM_2T_FROM_CONSVAR(u,nvar,i1,i2) (u(3,i1,i2))
#define Y_MOMENTUM_2L_FROM_CONSVAR(u,nvar,i1,i2) (u(i1,i2,3))

#if 0
! Compute z-momentum from conservative variables in 3D
#endif
#define Z_MOMENTUM_FROM_CONSVAR(u,nvar) (u(4))
#define Z_MOMENTUM_1T_FROM_CONSVAR(u,nvar,i1) (u(4,i1))
#define Z_MOMENTUM_1L_FROM_CONSVAR(u,nvar,i1) (u(i1,4))
#define Z_MOMENTUM_2T_FROM_CONSVAR(u,nvar,i1,i2) (u(4,i1,i2))
#define Z_MOMENTUM_2L_FROM_CONSVAR(u,nvar,i1,i2) (u(i1,i2,4))

#if 0
! Compute x-component of the magnetic fiels from conservative variables in 1D
#endif
#define X_MAGNETICFIELD_FROM_CONSVAR_1D(u,nvar) (X_MAGNETICFIELD_CONSTANT_1D)
#define X_MAGNETICFIELD_1T_FROM_CONSVAR_1D(u,nvar,i1) (X_MAGNETICFIELD_CONSTANT_1D)
#define X_MAGNETICFIELD_1L_FROM_CONSVAR_1D(u,nvar,i1) (X_MAGNETICFIELD_CONSTANT_1D)
#define X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(u,nvar,i1,i2) (X_MAGNETICFIELD_CONSTANT_1D)
#define X_MAGNETICFIELD_2L_FROM_CONSVAR_1D(u,nvar,i1,i2) (X_MAGNETICFIELD_CONSTANT_1D)

#if 0
! Compute x-component of the magnetic fiels from conservative variables in 2D, and 3D
#endif
#define X_MAGNETICFIELD_FROM_CONSVAR(u,nvar) (u(nvar/2+1))
#define X_MAGNETICFIELD_1T_FROM_CONSVAR(u,nvar,i1) (u(nvar/2+1,i1))
#define X_MAGNETICFIELD_1L_FROM_CONSVAR(u,nvar,i1) (u(i1,nvar/2+1))
#define X_MAGNETICFIELD_2T_FROM_CONSVAR(u,nvar,i1,i2) (u(nvar/2+1,i1,i2))
#define X_MAGNETICFIELD_2L_FROM_CONSVAR(u,nvar,i1,i2) (u(i1,i2,nvar/2+1))

#if 0
! Compute y-component of the magnetic fiels from conservative variables in 1D
#endif
#define Y_MAGNETICFIELD_FROM_CONSVAR_1D(u,nvar) (u((nvar-1)/2+2))
#define Y_MAGNETICFIELD_1T_FROM_CONSVAR_1D(u,nvar,i1) (u((nvar-1)/2+2,i1))
#define Y_MAGNETICFIELD_1L_FROM_CONSVAR_1D(u,nvar,i1) (u(i1,(nvar-1)/2+2))
#define Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(u,nvar,i1,i2) (u((nvar-1)/2+2,i1,i2))
#define Y_MAGNETICFIELD_2L_FROM_CONSVAR_1D(u,nvar,i1,i2) (u(i1,i2,(nvar-1)/2+2))

#if 0
! Compute y-component of the magnetic fiels from conservative variables in 2D, and 3D
#endif
#define Y_MAGNETICFIELD_FROM_CONSVAR(u,nvar) (u(nvar/2+2))
#define Y_MAGNETICFIELD_1T_FROM_CONSVAR(u,nvar,i1) (u(nvar/2+2,i1))
#define Y_MAGNETICFIELD_1L_FROM_CONSVAR(u,nvar,i1) (u(i1,nvar/2+2))
#define Y_MAGNETICFIELD_2T_FROM_CONSVAR(u,nvar,i1,i2) (u(nvar/2+2,i1,i2))
#define Y_MAGNETICFIELD_2L_FROM_CONSVAR(u,nvar,i1,i2) (u(i1,i2,nvar/2+2))

#if 0
! Compute z-component of the magnetic fiels from conservative variables in 1D
#endif
#define Z_MAGNETICFIELD_FROM_CONSVAR_1D(u,nvar) (u((nvar-1)/2+3))
#define Z_MAGNETICFIELD_1T_FROM_CONSVAR_1D(u,nvar,i1) (u((nvar-1)/2+3,i1))
#define Z_MAGNETICFIELD_1L_FROM_CONSVAR_1D(u,nvar,i1) (u(i1,(nvar-1)/2+3))
#define Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(u,nvar,i1,i2) (u((nvar-1)/2+3,i1,i2))
#define Z_MAGNETICFIELD_2L_FROM_CONSVAR_1D(u,nvar,i1,i2) (u(i1,i2,(nvar-1)/2+3))

#if 0
! Compute z-component of the magnetic fiels from conservative variables in 2D, and 3D
#endif
#define Z_MAGNETICFIELD_FROM_CONSVAR(u,nvar) (u(nvar/2+3))
#define Z_MAGNETICFIELD_1T_FROM_CONSVAR(u,nvar,i1) (u(nvar/2+3,i1))
#define Z_MAGNETICFIELD_1L_FROM_CONSVAR(u,nvar,i1) (u(i1,nvar/2+3))
#define Z_MAGNETICFIELD_2T_FROM_CONSVAR(u,nvar,i1,i2) (u(nvar/2+3,i1,i2))
#define Z_MAGNETICFIELD_2L_FROM_CONSVAR(u,nvar,i1,i2) (u(i1,i2,nvar/2+3))

#if 0
! Compute total energy from conservative variables in 1D, 2D, and 3D
#endif
#define TOTAL_ENERGY_FROM_CONSVAR(u,nvar) (u(nvar))
#define TOTAL_ENERGY_1T_FROM_CONSVAR(u,nvar,i1) (u(nvar,i1))
#define TOTAL_ENERGY_1L_FROM_CONSVAR(u,nvar,i1) (u(i1,nvar))
#define TOTAL_ENERGY_2T_FROM_CONSVAR(u,nvar,i1,i2) (u(nvar,i1,i2))
#define TOTAL_ENERGY_2L_FROM_CONSVAR(u,nvar,i1,i2) (u(i1,i2,nvar))

#if 0
!##############################################################################
! Roe average values
!##############################################################################
#endif

#define ROE_MEAN_RATIO(ul,ur) (ul/ur)
#define ROE_MEAN_VALUE(ul,ur,ratio) ((ratio*ul+ur)/(ratio+1))

#if 0
! Include magnetohydrodynamic header file
#endif

#include "../../kernel/magnetohydrodynamics.h"

#endif
