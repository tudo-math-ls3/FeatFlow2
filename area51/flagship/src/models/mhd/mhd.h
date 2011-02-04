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
#ifdef MHD_VACUUM_PERM
#define VACUUM_PERM MHD_VACUUM_PERM
#else
#define VACUUM_PERM 1.0_DP
#endif

#if 0
! Ratio of specific heats for air at sea-level
#endif
#ifdef MHD_GAMMA
#define GAMMA MHD_GAMMA
#else
#define GAMMA 1.4_DP
#endif

#if 0
! In one dimension, the x-component of the magnetic field is constant.
#endif
#ifdef MHD_X_MAGNETICFIELD_CONSTANT_1D
#define X_MAGNETICFIELD_CONSTANT_1D MHD_X_MAGNETICFIELD_CONSTANT_1D
#else
#define X_MAGNETICFIELD_CONSTANT_1D 1.0_DP
#endif


#if 0
!##############################################################################
! Conservative variables
!##############################################################################
#endif

#if 0
! Compute density from conservative variables in 1D,2D, and 3D
#endif
#define DENSITY_FROM_CONSVAR(U,nvar) (U(1))
#define DENSITY_1T_FROM_CONSVAR(U,nvar,i1) (U(1,i1))
#define DENSITY_1L_FROM_CONSVAR(U,nvar,i1) (U(i1,1))
#define DENSITY_2T_FROM_CONSVAR(U,nvar,i1,i2) (U(1,i1,i2))
#define DENSITY_2L_FROM_CONSVAR(U,nvar,i1,i2) (U(i1,i2,1))

#if 0
! Compute x-momentum from conservative variables in 1D, 2D, and 3D
#endif
#define X_MOMENTUM_FROM_CONSVAR(U,nvar) (U(2))
#define X_MOMENTUM_1T_FROM_CONSVAR(U,nvar,i1) (U(2,i1))
#define X_MOMENTUM_1L_FROM_CONSVAR(U,nvar,i1) (U(i1,2))
#define X_MOMENTUM_2T_FROM_CONSVAR(U,nvar,i1,i2) (U(2,i1,i2))
#define X_MOMENTUM_2L_FROM_CONSVAR(U,nvar,i1,i2) (U(i1,i2,2))

#if 0
! Compute y-momentum from conservative variables in 2D, and 3D
#endif
#define Y_MOMENTUM_FROM_CONSVAR(U,nvar) (U(3))
#define Y_MOMENTUM_1T_FROM_CONSVAR(U,nvar,i1) (U(3,i1))
#define Y_MOMENTUM_1L_FROM_CONSVAR(U,nvar,i1) (U(i1,3))
#define Y_MOMENTUM_2T_FROM_CONSVAR(U,nvar,i1,i2) (U(3,i1,i2))
#define Y_MOMENTUM_2L_FROM_CONSVAR(U,nvar,i1,i2) (U(i1,i2,3))

#if 0
! Compute z-momentum from conservative variables in 3D
#endif
#define Z_MOMENTUM_FROM_CONSVAR(U,nvar) (U(4))
#define Z_MOMENTUM_1T_FROM_CONSVAR(U,nvar,i1) (U(4,i1))
#define Z_MOMENTUM_1L_FROM_CONSVAR(U,nvar,i1) (U(i1,4))
#define Z_MOMENTUM_2T_FROM_CONSVAR(U,nvar,i1,i2) (U(4,i1,i2))
#define Z_MOMENTUM_2L_FROM_CONSVAR(U,nvar,i1,i2) (U(i1,i2,4))

#if 0
! Compute x-component of the magnetic fiels from conservative variables in 1D
#endif
#define X_MAGNETICFIELD_FROM_CONSVAR_1D(U,nvar) (X_MAGNETICFIELD_CONSTANT_1D)
#define X_MAGNETICFIELD_1T_FROM_CONSVAR_1D(U,nvar,i1) (X_MAGNETICFIELD_CONSTANT_1D)
#define X_MAGNETICFIELD_1L_FROM_CONSVAR_1D(U,nvar,i1) (X_MAGNETICFIELD_CONSTANT_1D)
#define X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,nvar,i1,i2) (X_MAGNETICFIELD_CONSTANT_1D)
#define X_MAGNETICFIELD_2L_FROM_CONSVAR_1D(U,nvar,i1,i2) (X_MAGNETICFIELD_CONSTANT_1D)

#if 0
! Compute x-component of the magnetic fiels from conservative variables in 2D, and 3D
#endif
#define X_MAGNETICFIELD_FROM_CONSVAR(U,nvar) (U(5))
#define X_MAGNETICFIELD_1T_FROM_CONSVAR(U,nvar,i1) (U(5,i1))
#define X_MAGNETICFIELD_1L_FROM_CONSVAR(U,nvar,i1) (U(i1,5))
#define X_MAGNETICFIELD_2T_FROM_CONSVAR(U,nvar,i1,i2) (U(5,i1,i2))
#define X_MAGNETICFIELD_2L_FROM_CONSVAR(U,nvar,i1,i2) (U(i1,i2,5))

#if 0
! Compute y-component of the magnetic fiels from conservative variables in 1D
#endif
#define Y_MAGNETICFIELD_FROM_CONSVAR_1D(U,nvar) (U(5))
#define Y_MAGNETICFIELD_1T_FROM_CONSVAR_1D(U,nvar,i1) (U(5,i1))
#define Y_MAGNETICFIELD_1L_FROM_CONSVAR_1D(U,nvar,i1) (U(i1,5))
#define Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,nvar,i1,i2) (U(5,i1,i2))
#define Y_MAGNETICFIELD_2L_FROM_CONSVAR_1D(U,nvar,i1,i2) (U(i1,i2,5))

#if 0
! Compute y-component of the magnetic fiels from conservative variables in 2D, and 3D
#endif
#define Y_MAGNETICFIELD_FROM_CONSVAR(U,nvar) (U(6))
#define Y_MAGNETICFIELD_1T_FROM_CONSVAR(U,nvar,i1) (U(6,i1))
#define Y_MAGNETICFIELD_1L_FROM_CONSVAR(U,nvar,i1) (U(i1,6))
#define Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,nvar,i1,i2) (U(6,i1,i2))
#define Y_MAGNETICFIELD_2L_FROM_CONSVAR(U,nvar,i1,i2) (U(i1,i2,6))

#if 0
! Compute z-component of the magnetic fiels from conservative variables in 1D
#endif
#define Z_MAGNETICFIELD_FROM_CONSVAR_1D(U,nvar) (U(6))
#define Z_MAGNETICFIELD_1T_FROM_CONSVAR_1D(U,nvar,i1) (U(6,i1))
#define Z_MAGNETICFIELD_1L_FROM_CONSVAR_1D(U,nvar,i1) (U(i1,6))
#define Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,nvar,i1,i2) (U(6,i1,i2))
#define Z_MAGNETICFIELD_2L_FROM_CONSVAR_1D(U,nvar,i1,i2) (U(i1,i2,6))

#if 0
! Compute z-component of the magnetic fiels from conservative variables in 2D, and 3D
#endif
#define Z_MAGNETICFIELD_FROM_CONSVAR(U,nvar) (U(7))
#define Z_MAGNETICFIELD_1T_FROM_CONSVAR(U,nvar,i1) (U(7,i1))
#define Z_MAGNETICFIELD_1L_FROM_CONSVAR(U,nvar,i1) (U(i1,7))
#define Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,nvar,i1,i2) (U(7,i1,i2))
#define Z_MAGNETICFIELD_2L_FROM_CONSVAR(U,nvar,i1,i2) (U(i1,i2,7))

#if 0
! Compute total energy from conservative variables in 1D, 2D, and 3D
#endif
#define TOTAL_ENERGY_FROM_CONSVAR(U,nvar) (U(nvar))
#define TOTAL_ENERGY_1T_FROM_CONSVAR(U,nvar,i1) (U(nvar,i1))
#define TOTAL_ENERGY_1L_FROM_CONSVAR(U,nvar,i1) (U(i1,nvar))
#define TOTAL_ENERGY_2T_FROM_CONSVAR(U,nvar,i1,i2) (U(nvar,i1,i2))
#define TOTAL_ENERGY_2L_FROM_CONSVAR(U,nvar,i1,i2) (U(i1,i2,nvar))

#if 0
!##############################################################################
! Roe average values
!##############################################################################
#endif

#define ROE_MEAN_RATIO(ul,ur) (sqrt(ul/ur))
#define ROE_MEAN_VALUE(ul,ur,ratio) ((ratio*ul+MYNEWLINE ur)/MYNEWLINE (ratio+1))

#if 0
!##############################################################################
! Include magnetohydrodynamic header file
!##############################################################################
#endif

#include "../../kernel/magnetohydrodynamics.h"

#if 0
!##############################################################################
! Fluxes and matrices in 1D
!##############################################################################
#endif

#if 0
! Flux in x-direction for ideal MHD in 1D
#endif
#define FLUX_MHD_2T_XDIR_1D(F,U,i,idx,ui,vi,wi,pi,qi)\
  F(1) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,i,idx);MYNEWLINE\
  F(2) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,i,idx)*ui+pi;MYNEWLINE\
  F(3) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,i,idx)*ui-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx)*MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx);MYNEWLINE\
  F(4) = Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,i,idx)*ui-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx)*MYNEWLINE\
         Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx);MYNEWLINE\
  F(5) = Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx)*ui-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx)*vi;MYNEWLINE\
  F(6) = Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx)*ui-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx)*wi;MYNEWLINE\
  F(7) = (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR1D,i,idx)+pi)*ui-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx)*qi

#if 0
! Flux difference in x-direction for ideal MHD in 1D
#endif
#define FLUXDIFF_MHD_2T_XDIR_1D(F,U,i,j,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)\
  F(1) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,i,idx)-MYNEWLINE\
         X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,j,idx);MYNEWLINE\
  F(2) = (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,i,idx)*ui+pi)-MYNEWLINE\
         (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,j,idx)*uj+pj);MYNEWLINE\
  F(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,i,idx)*ui-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx)*MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx))-MYNEWLINE\
         (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,j,idx)*uj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,j,idx)*MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,j,idx));MYNEWLINE\
  F(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,i,idx)*ui-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx)*MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx))-MYNEWLINE\
         (Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,j,idx)*uj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,j,idx)*MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,j,idx));MYNEWLINE\
  F(5) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx)*ui-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx)*vi)-MYNEWLINE\
         (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,j,idx)*uj-MYNEWLINE\
	  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,j,idx)*vj);MYNEWLINE\
  F(6) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx)*ui-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx)*wi)-MYNEWLINE\
         (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,j,idx)*uj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,j,idx)*wj);MYNEWLINE\
  F(7) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR1D,i,idx)+pi)*ui-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,i,idx)*qi)-MYNEWLINE\
         ((TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR1D,j,idx)+pj)*uj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(U,NVAR1D,j,idx)*qj)

#if 0
!##############################################################################
! Fluxes and matrices in 2D
!##############################################################################
#endif

#if 0
! Flux in x-direction for ideal MHD in 2D
#endif
#define FLUX_MHD_2T_XDIR_2D(F,U,i,idx,ui,vi,wi,pi,qi)\
  F(1) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx);MYNEWLINE\
  F(2) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui+pi-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)**2;MYNEWLINE\
  F(3) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx);MYNEWLINE\
  F(4) = Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*MYNEWLINE\
         Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx);MYNEWLINE\
  F(5) = 0.0_DP;MYNEWLINE\
  F(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi;MYNEWLINE\
  F(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*wi;MYNEWLINE\
  F(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR2D,i,idx)+pi)*ui-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*qi

#if 0
! Flux in y-direction for ideal MHD in 2D
#endif
#define FLUX_MHD_2T_YDIR_2D(F,U,i,idx,ui,vi,wi,pi,qi)\
  F(1) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx);MYNEWLINE\
  F(2) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi-MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx);MYNEWLINE\
  F(3) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi+pi-MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)**2;MYNEWLINE\
  F(4) = Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi-MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*MYNEWLINE\
         Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx);MYNEWLINE\
  F(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi-MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui;MYNEWLINE\
  F(6) = 0.0_DP;MYNEWLINE\
  F(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi-MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*wi;MYNEWLINE\
  F(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR2D,i,idx)+pi)*vi-MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*qi

#if 0
! Flux difference in x-direction for ideal MHD in 2D
#endif
#define FLUXDIFF_MHD_2T_XDIR_2D(F,U,i,j,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)\
  F(1) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)-MYNEWLINE\
         X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,j,idx);MYNEWLINE\
  F(2) = (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui+pi-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)**2)-MYNEWLINE\
         (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*uj+pj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)**2);MYNEWLINE\
  F(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx))-MYNEWLINE\
         (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*uj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx));MYNEWLINE\
  F(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx))-MYNEWLINE\
         (Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*uj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx));MYNEWLINE\
  F(5) = 0.0_DP;MYNEWLINE\
  F(6) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi)-MYNEWLINE\
         (Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*uj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*vj);MYNEWLINE\
  F(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*wi)-MYNEWLINE\
         (Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*uj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*wj);MYNEWLINE\
  F(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR2D,i,idx)+pi)*ui-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*qi)-MYNEWLINE\
         ((TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR2D,j,idx)+pj)*uj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*qj)

#if 0
! Flux difference in y-direction for ideal MHD in 2D
#endif
#define FLUXDIFF_MHD_2T_YDIR_2D(F,U,i,j,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)\
  F(1) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)-MYNEWLINE\
          Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,j,idx);MYNEWLINE\
  F(2) = (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx))-MYNEWLINE\
         (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*vj-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx));MYNEWLINE\
  F(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi+pi-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)**2)-MYNEWLINE\
         (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*vj+pj-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)**2);MYNEWLINE\
  F(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx))-MYNEWLINE\
         (Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*vj-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx));MYNEWLINE\
  F(5) = (X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui)-MYNEWLINE\
         (X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*vj-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*uj);MYNEWLINE\
  F(6) = 0.0_DP;MYNEWLINE\
  F(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*wi)-MYNEWLINE\
         (Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*vj-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*wj);MYNEWLINE\
  F(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR2D,i,idx)+pi)*vi-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*qi)-MYNEWLINE\
         ((TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR2D,j,idx)+pj)*vj-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*qj)

#if 0
!##############################################################################
! Fluxes and matrices in 3D
!##############################################################################
#endif

#if 0
! Flux in x-direction for ideal MHD in 3D
#endif
#define FLUX_MHD_2T_XDIR_3D(F,U,i,idx,ui,vi,wi,pi,qi)\
  F(1) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx);MYNEWLINE\
  F(2) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui+pi-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)**2;MYNEWLINE\
  F(3) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx);MYNEWLINE\
  F(4) = Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*MYNEWLINE\
         Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx);MYNEWLINE\
  F(5) = 0.0_DP;MYNEWLINE\
  F(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi;MYNEWLINE\
  F(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi;MYNEWLINE\
  F(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,i,idx)+pi)*ui-MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*qi

#if 0
! Flux in y-direction for ideal MHD in 3D
#endif
#define FLUX_MHD_2T_YDIR_3D(F,U,i,idx,ui,vi,wi,pi,qi)\
  F(1) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx);MYNEWLINE\
  F(2) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi-MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx);MYNEWLINE\
  F(3) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi+pi-MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)**2;MYNEWLINE\
  F(4) = Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi-MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*MYNEWLINE\
         Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx);MYNEWLINE\
  F(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi-MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui;MYNEWLINE\
  F(6) = 0.0_DP;MYNEWLINE\
  F(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi-MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi;MYNEWLINE\
  F(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,i,idx)+pi)*vi-MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*qi

#if 0
! Flux in y-direction for ideal MHD in 3D
#endif
#define FLUX_MHD_2T_ZDIR_3D(F,U,i,idx,ui,vi,wi,pi,qi)\
  F(1) = Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx);\
  F(2) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi-MYNEWLINE\
         Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*MYNEWLINE\
         X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx);MYNEWLINE\
  F(3) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi-MYNEWLINE\
         Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*MYNEWLINE\
         Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx);MYNEWLINE\
  F(4) = Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi+pi-MYNEWLINE\
         Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)**2;MYNEWLINE\
  F(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi-MYNEWLINE\
         Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui;MYNEWLINE\
  F(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi-MYNEWLINE\
         Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi;MYNEWLINE\
  F(7) = 0.0_DP;MYNEWLINE\
  F(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,i,idx)+pi)*wi-MYNEWLINE\
         Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*qi

#if 0
! Flux difference in x-direction for ideal MHD in 3D
#endif
#define FLUXDIFF_MHD_2T_XDIR_3D(F,U,i,j,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)\
  F(1) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)-MYNEWLINE\
         X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx);MYNEWLINE\
  F(2) = (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui+pi-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)**2)-MYNEWLINE\
         (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*uj+pj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)**2);MYNEWLINE\
  F(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx))-MYNEWLINE\
         (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*uj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx));MYNEWLINE\
  F(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx))-MYNEWLINE\
         (Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*uj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx));MYNEWLINE\
  F(5) = 0.0_DP;MYNEWLINE\
  F(6) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi)-MYNEWLINE\
         (Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*uj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*vj);MYNEWLINE\
  F(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi)-MYNEWLINE\
         (Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*uj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*wj);MYNEWLINE\
  F(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,i,idx)+pi)*ui-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*qi)-MYNEWLINE\
         ((TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,j,idx)+pj)*uj-MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*qj)

#if 0
! Flux difference in y-direction for ideal MHD in 3D
#endif
#define FLUXDIFF_MHD_2T_YDIR_3D(F,U,i,j,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)\
  F(1) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)-MYNEWLINE\
          Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx);MYNEWLINE\
  F(2) = (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx))-MYNEWLINE\
         (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*vj-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx));MYNEWLINE\
  F(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi+pi-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)**2)-MYNEWLINE\
         (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*vj+pj-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)**2);MYNEWLINE\
  F(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx))-MYNEWLINE\
         (Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*vj-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx));MYNEWLINE\
  F(5) = (X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui)-MYNEWLINE\
         (X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*vj-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*uj);MYNEWLINE\
  F(6) = 0.0_DP;MYNEWLINE\
  F(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi)-MYNEWLINE\
         (Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*vj-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*wj);MYNEWLINE\
  F(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,i,idx)+pi)*vi-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*qi)-MYNEWLINE\
         ((TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,j,idx)+pj)*vj-MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*qj)

#if 0
! Flux difference in z-direction for ideal MHD in 3D
#endif
#define FLUXDIFF_MHD_2T_ZDIR_3D(F,U,i,j,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)\
  F(1) = Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)-MYNEWLINE\
         Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx);MYNEWLINE\
  F(2) = (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi-MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx))-MYNEWLINE\
          (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*wj-MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*MYNEWLINE\
          X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx));MYNEWLINE\
  F(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi-MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx))-MYNEWLINE\
          (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*wj-MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*MYNEWLINE\
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx));MYNEWLINE\
  F(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi+pi-MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)**2)-MYNEWLINE\
          (Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*wj+pj-MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)**2);MYNEWLINE\
  F(5) = (X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi-MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui)-MYNEWLINE\
          (X_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*wj-MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*uj);MYNEWLINE\
  F(6) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi-MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi)-MYNEWLINE\
          (Y_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*wj-MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*vj);MYNEWLINE\
  F(7) = 0.0_DP;MYNEWLINE\
  F(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,i,idx)+pi)*wi-MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*qi)-MYNEWLINE\
          ((TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,j,idx)+pj)*wj-MYNEWLINE\
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*qj)

#endif
