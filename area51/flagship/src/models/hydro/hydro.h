#if 0
!-*- mode: f90; -*-
!##############################################################################
!# ****************************************************************************
!# <name> hydro </name>
!# ****************************************************************************
!#
!# <purpose>
!# This is the main header file for the hydrodynamic model
!#
!# </purpose>
!##############################################################################
#endif

#ifndef _HYDRO_H_
#define _HYDRO_H_

#include "../../flagship.h"

#if 0
! Equation of state for perfect gas
#endif
#define PERFECT_GAS

#if 0
! Ratio of specific heats for air at sea-level
#endif
#define GAMMA 1.4_DP

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
! Compute total energy from conservative variables in 1D, 2d, and 3D
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
! Include thermodynamic header file
#endif

#include "thermodynamics.h"

#endif
