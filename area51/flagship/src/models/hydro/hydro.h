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
#ifdef HYDRO_GAMMA
#define GAMMA HYDRO_GAMMA
#else
#define GAMMA 1.4_DP
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
! Compute total energy from conservative variables in 1D, 2d, and 3D
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

#define ROE_MEAN_RATIO(ul,ur) (ul/ur)
#define ROE_MEAN_VALUE(ul,ur,ratio) ((ratio*ul+MYNEWLINE ur)/MYNEWLINE (ratio+1))

#if 0
!##############################################################################
! Include thermodynamic header file
!##############################################################################
#endif

#include "../../kernel/thermodynamics.h"

#if 0
!##############################################################################
! Fluxes and matrices in 1D
!##############################################################################
#endif

#if 0
! Flux in x-direction for 
#endif
#define FLUX_HYDRO_2T_XDIR_1D(F,U,i,idx,ui,pi)\
  F(1) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,i,idx);MYNEWLINE\
  F(2) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,i,idx)*ui+pi;MYNEWLINE\
  F(3) = TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR1D,i,idx)*ui+pi*ui

#if 0
! Flux difference in x-direction for inviscid hydrodynamics in 1D
#endif
#define FLUXDIFF_HYDRO_2T_XDIR_1D(F,U,i,j,idx,ui,uj,pi,pj)\
  F(1) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,i,idx)-MYNEWLINE\
         X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,j,idx);MYNEWLINE\
  F(2) = (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,i,idx)*ui+pi)-MYNEWLINE\
         (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR1D,j,idx)*uj+pj);MYNEWLINE\
  F(3) = (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR1D,i,idx)*ui+pi*ui)-MYNEWLINE\
         (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR1D,j,idx)*uj+pj*uj)

#if 0
! Diagonal matrix for inviscid hydrodynamics in 1D
#endif
#define MATRIXDIAG_HYDRO_2T_1D(K,i,idx,dscale,CX,ui)\
  K(1,i,idx) = 0.0;MYNEWLINE\
  K(2,i,idx) = dscale*(3.0-GAMMA)*ui*CX;MYNEWLINE\
  K(3,i,idx) = dscale*GAMMA*ui*CX

#if 0
! Full matrix for inviscid hydrodynamics in 1D
#endif
#define MATRIX_HYDRO_2T_1D(K,i,idx,dscale,CX,ui,Ei)\
  K(1,i,idx) = 0.0;MYNEWLINE\
  K(2,i,idx) = dscale*(GAMMA-3.0)/2.0*ui*ui*CX;MYNEWLINE\
  K(3,i,idx) = dscale*((GAMMA-1.0)*ui*ui-GAMMA*Ei)*ui*CX;MYNEWLINE\
  K(4,i,idx) = dscale*CX;MYNEWLINE\
  K(5,i,idx) = dscale*(3.0-GAMMA)*ui*CX;MYNEWLINE\
  K(6,i,idx) = dscale*(GAMMA*Ei-3.0*(GAMMA-1.0)/2.0*ui*ui)*CX;MYNEWLINE\
  K(7,i,idx) = 0.0;MYNEWLINE\
  K(8,i,idx) = dscale*(GAMMA-1.0)*CX;MYNEWLINE\
  K(9,i,idx) = dscale*GAMMA*ui*CX


#if 0
!##############################################################################
! Fluxes and matrices in 2D
!##############################################################################
#endif

#if 0
! Flux in x-direction for inviscid hydrodynamics in 2D
#endif
#define FLUX_HYDRO_2T_XDIR_2D(F,U,i,idx,ui,pi)\
  F(1) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx);MYNEWLINE\
  F(2) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui+pi;MYNEWLINE\
  F(3) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui;MYNEWLINE\
  F(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui+pi*ui
  
#if 0
! Flux in y-direction for inviscid hydrodynamics in 2D
#endif
#define FLUX_HYDRO_2T_YDIR_2D(F,U,i,idx,vi,pi)\
  F(1) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx);MYNEWLINE\
  F(2) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi;MYNEWLINE\
  F(3) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi+pi;MYNEWLINE\
  F(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi+pi*vi

#if 0
! Flux difference in x-direction for inviscid hydrodynamics in 2D
#endif
#define FLUXDIFF_HYDRO_2T_XDIR_2D(F,U,i,j,idx,ui,uj,pi,pj)\
  F(1) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)-MYNEWLINE\
         X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,j,idx);MYNEWLINE\
  F(2) = (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui+pi)-MYNEWLINE\
         (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*uj+pj);MYNEWLINE\
  F(3) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui-MYNEWLINE\
	 Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*uj;MYNEWLINE\
  F(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*ui+pi*ui)-MYNEWLINE\
         (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*uj+pj*uj)

#if 0
! Flux difference in y-direction for inviscid hydrodynamics in 2D
#endif
#define FLUXDIFF_HYDRO_2T_YDIR_2D(F,U,i,j,idx,vi,vj,pi,pj)\
  F(1) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)-MYNEWLINE\
         Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,j,idx);MYNEWLINE\
  F(2) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi-MYNEWLINE\
         X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*vj;MYNEWLINE\
  F(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi+pi)-MYNEWLINE\
	 (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*vj+pj);MYNEWLINE\
  F(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR2D,i,idx)*vi+pi*vi)-MYNEWLINE\
         (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR2D,j,idx)*vj+pj*vj)

#if 0
! Diagonal matrix for inviscid hydrodynamics in 2D
#endif
#define MATRIXDIAG_HYDRO_2T_2D(K,i,idx,dscale,CX,CY,ui,vi)\
  K(1,i,idx) = 0.0;MYNEWLINE\
  K(2,i,idx) = dscale*((3.0-GAMMA)*ui*CX+MYNEWLINE vi*CY);MYNEWLINE\
  K(3,i,idx) = dscale*(ui*CX+MYNEWLINE (3.0-GAMMA)*vi*CY);MYNEWLINE\
  K(4,i,idx) = dscale*(GAMMA*(ui*CX+MYNEWLINE vi*CY))

#if 0
! Full matrix for inviscid hydrodynamics in 2D
#endif
#define MATRIX_HYDRO_2T_2D(K,i,idx,dscale,CX,CY,ui,vi,Ei)\
  K( 1,i,idx) = 0.0;MYNEWLINE\
  K( 2,i,idx) = dscale*(((GAMMA-1.0)/2.0*(ui*ui+vi*vi)-ui*ui)*CX-MYNEWLINE ui*vi*CY);MYNEWLINE\
  K( 3,i,idx) = dscale*(((GAMMA-1.0)/2.0*(ui*ui+vi*vi)-vi*vi)*CY-MYNEWLINE ui*vi*CX);MYNEWLINE\
  K( 4,i,idx) = dscale*((GAMMA-1.0)*(ui*ui+vi*vi)-GAMMA*Ei)*(ui*CX+MYNEWLINE vi*CY);MYNEWLINE\
  K( 5,i,idx) = dscale*CX;MYNEWLINE\
  K( 6,i,idx) = dscale*((3.0-GAMMA)*ui*CX+MYNEWLINE vi*CY);MYNEWLINE\
  K( 7,i,idx) = dscale*(vi*CX-MYNEWLINE (GAMMA-1.0)*ui*CY);MYNEWLINE\
  K( 8,i,idx) = dscale*((GAMMA*Ei-(GAMMA-1.0)/2.0*(ui*ui+vi*vi))*CX-MYNEWLINE (GAMMA-1.0)*ui*(ui*CX+MYNEWLINE vi*CY));MYNEWLINE\
  K( 9,i,idx) = dscale*CY;MYNEWLINE\
  K(10,i,idx) = dscale*(ui*CY-MYNEWLINE (GAMMA-1.0)*vi*CX);MYNEWLINE\
  K(11,i,idx) = dscale*(ui*CX+MYNEWLINE (3.0-GAMMA)*vi*CY);MYNEWLINE\
  K(12,i,idx) = dscale*((GAMMA*Ei-(GAMMA-1.0)/2.0*(ui*ui+vi*vi))*CY-MYNEWLINE (GAMMA-1.0)*vi*(ui*CX+MYNEWLINE vi*CY));MYNEWLINE\
  K(13,i,idx) = 0.0;MYNEWLINE\
  K(14,i,idx) = dscale*(GAMMA-1.0)*CX;MYNEWLINE\
  K(15,i,idx) = dscale*(GAMMA-1.0)*CY;MYNEWLINE\
  K(16,i,idx) = dscale*(GAMMA*(ui*CX+MYNEWLINE vi*CY))

#if 0
!##############################################################################
! Fluxes and matrices in 3D
!##############################################################################
#endif

#if 0
! Flux in x-direction for inviscid hydrodynamics in 3D
#endif
#define FLUX_HYDRO_2T_XDIR_3D(F,U,i,idx,ui,pi)\
  F(1) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx);MYNEWLINE\
  F(2) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui+pi;MYNEWLINE\
  F(3) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui;MYNEWLINE\
  F(4) = Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui;MYNEWLINE\
  F(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui+pi*ui

#if 0
! Flux in y-direction for inviscid hydrodynamics in 3D
#endif
#define FLUX_HYDRO_2T_YDIR_3D(F,U,i,idx,vi,pi)\
  F(1) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx);MYNEWLINE\
  F(2) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi;MYNEWLINE\
  F(3) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi+pi;MYNEWLINE\
  F(4) = Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi;MYNEWLINE\
  F(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi+pi*vi

#if 0
! Flux in z-direction for inviscid hydrodynamics in 3D
#endif
#define FLUX_HYDRO_2T_ZDIR_3D(F,U,i,idx,wi,pi)\
  F(1) = Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx);MYNEWLINE\
  F(2) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi;MYNEWLINE\
  F(3) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi;MYNEWLINE\
  F(4) = Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi+pi;MYNEWLINE\
  F(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi+pi*wi

#if 0
! Flux difference in x-direction for inviscid hydrodynamics in 3D
#endif
#define FLUXDIFF_HYDRO_2T_XDIR_3D(F,U,i,j,idx,ui,uj,pi,pj)\
  F(1) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)-MYNEWLINE\
         X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx);MYNEWLINE\
  F(2) = (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui+pi)-MYNEWLINE\
         (X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*uj+pj);MYNEWLINE\
  F(3) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui-MYNEWLINE\
         Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*uj;MYNEWLINE\
  F(4) = Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui-MYNEWLINE\
         Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*uj;MYNEWLINE\
  F(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*ui+pi*ui)-MYNEWLINE\
         (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*uj+pj*uj)

#if 0
! Flux difference in y-direction for inviscid hydrodynamics in 3D
#endif
#define FLUXDIFF_HYDRO_2T_YDIR_3D(F,U,i,j,idx,vi,vj,pi,pj)\
  F(1) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)-MYNEWLINE\
         Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx);MYNEWLINE\
  F(2) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi-MYNEWLINE\
	 X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*vj;MYNEWLINE\
  F(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi+pi)-MYNEWLINE\
         (Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*vj+pj);MYNEWLINE\
  F(4) = Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi-MYNEWLINE\
	 Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*vj;MYNEWLINE\
  F(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*vi+pi*vi)-MYNEWLINE\
         (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*vj+pj*vj)

#if 0
! Flux difference in z-direction for inviscid hydrodynamics in 3D
#endif
#define FLUXDIFF_HYDRO_2T_ZDIR_3D(F,U,i,j,idx,wi,wj,pi,pj)\
  F(1) = Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)-MYNEWLINE\
         Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx);MYNEWLINE\
  F(2) = X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi-MYNEWLINE\
	 X_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*wj;MYNEWLINE\
  F(3) = Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi-MYNEWLINE\
	 Y_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*wj;MYNEWLINE\
  F(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi+pi)-MYNEWLINE\
         (Z_MOMENTUM_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*wj+pj);MYNEWLINE\
  F(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,i,idx)*wi+pi*wi)-MYNEWLINE\
         (TOTAL_ENERGY_2T_FROM_CONSVAR(U,NVAR3D,j,idx)*wj+pj*wj)

#if 0
! Diagonal matrix for inviscid hydrodynamics in 3D
#endif
#define MATRIXDIAG_HYDRO_2T_3D(K,i,idx,dscale,CX,CY,CZ,ui,vi,wi)\
  K(1,i,idx) = 0.0;MYNEWLINE\
  K(2,i,idx) = dscale*((3.0-GAMMA)*ui*CX+MYNEWLINE vi*CY+MYNEWLINE wi*CZ);MYNEWLINE\
  K(3,i,idx) = dscale*(ui*CX+MYNEWLINE (3.0-GAMMA)*vi*CY+MYNEWLINE wi*CZ);MYNEWLINE\
  K(4,i,idx) = dscale*(ui*CX+MYNEWLINE vi*CY+MYNEWLINE (3.0-GAMMA)*wi*CZ);MYNEWLINE\
  K(5,i,idx) = dscale*(GAMMA*(ui*CX+MYNEWLINE vi*CY+MYNEWLINE wi*CZ))

#if 0
! Diagonal matrix for inviscid hydrodynamics in 3D
#endif
#define MATRIX_HYDRO_2T_3D(K,i,idx,dscale,CX,CY,CZ,ui,vi,wi,Ei)	\
  K( 1,i,idx) = 0.0;MYNEWLINE\
  K( 2,i,idx) = dscale*(((GAMMA-1.0)/2.0*(ui*ui+vi*vi+wi*wi)-ui*ui)*CX-MYNEWLINE ui*vi*CY-MYNEWLINE ui*wi*CZ);MYNEWLINE\
  K( 3,i,idx) = dscale*(((GAMMA-1.0)/2.0*(ui*ui+vi*vi+wi*wi)-vi*vi)*CY-MYNEWLINE ui*vi*CX-MYNEWLINE vi*wi*CZ);MYNEWLINE\
  K( 4,i,idx) = dscale*(((GAMMA-1.0)/2.0*(ui*ui+vi*vi+wi*wi)-wi*wi)*CZ-MYNEWLINE ui*wi*CX-MYNEWLINE vi*wi*CY);MYNEWLINE\
  K( 5,i,idx) = dscale*((GAMMA-1.0)*(ui*ui+vi*vi+wi*wi)-MYNEWLINE GAMMA*Ei)*(ui*CX+MYNEWLINE vi*CY+MYNEWLINE wi*CZ);MYNEWLINE\
  K( 6,i,idx) = dscale*CX;MYNEWLINE\
  K( 7,i,idx) = dscale*((3.0-GAMMA)*ui*CX+MYNEWLINE vi*CY+MYNEWLINE wi*CZ);MYNEWLINE\
  K( 8,i,idx) = dscale*(vi*CX-MYNEWLINE (GAMMA-1.0)*ui*CY);MYNEWLINE\
  K( 9,i,idx) = dscale*(wi*CX-MYNEWLINE (GAMMA-1.0)*ui*CZ);MYNEWLINE\
  K(10,i,idx) = dscale*((GAMMA*Ei-(GAMMA-1.0)/2.0*(ui*ui+vi*vi+wi*wi))*CX-MYNEWLINE(GAMMA-1.0)*ui*(ui*CX+MYNEWLINE vi*CY+MYNEWLINE wi*CZ));MYNEWLINE\
  K(11,i,idx) = dscale*CY;MYNEWLINE\
  K(12,i,idx) = dscale*(ui*CY-MYNEWLINE (GAMMA-1.0)*vi*CX);MYNEWLINE\
  K(13,i,idx) = dscale*((3.0-GAMMA)*vi*CY+MYNEWLINE ui*CX+MYNEWLINE wi*CZ);MYNEWLINE\
  K(14,i,idx) = dscale*(wi*CY-MYNEWLINE (GAMMA-1.0)*vi*CZ);MYNEWLINE\
  K(15,i,idx) = dscale*((GAMMA*Ei-(GAMMA-1.0)/2.0*(ui*ui+vi*vi+wi*wi))*CY-MYNEWLINE (GAMMA-1.0)*vi*(ui*CX+MYNEWLINE vi*CY+MYNEWLINE wi*CZ));MYNEWLINE\
  K(16,i,idx) = dscale*CZ;MYNEWLINE\
  K(17,i,idx) = dscale*(ui*CZ-MYNEWLINE (GAMMA-1.0)*wi*CX);MYNEWLINE\
  K(18,i,idx) = dscale*(vi*CZ-MYNEWLINE (GAMMA-1.0)*wi*CY);MYNEWLINE\
  K(19,i,idx) = dscale*((3.0-GAMMA)*wi*CZ+MYNEWLINE ui*CX+MYNEWLINE vi*CY);MYNEWLINE\
  K(20,i,idx) = dscale*((GAMMA*Ei-(GAMMA-1.0)/2.0*(ui*ui+vi*vi+wi*wi))*CZ-MYNEWLINE (GAMMA-1.0)*wi*(ui*CX+MYNEWLINE vi*CY+MYNEWLINE wi*CZ));MYNEWLINE\
  K(21,i,idx) = 0.0;MYNEWLINE\
  K(22,i,idx) = dscale*(GAMMA-1.0)*CX;MYNEWLINE\
  K(23,i,idx) = dscale*(GAMMA-1.0)*CY;MYNEWLINE\
  K(24,i,idx) = dscale*(GAMMA-1.0)*CZ;MYNEWLINE\
  K(25,i,idx) = dscale*(GAMMA*(ui*CX+MYNEWLINE vi*CY+MYNEWLINE wi*CZ))
#endif
