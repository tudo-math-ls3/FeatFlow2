!##############################################################################
!# ****************************************************************************
!# <name> hydro_callback3d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible Euler/Navier-Stokes equations in 3D.
!#
!# The following callback functions are available:
!#
!# 1.) hydro_calcFluxGal3d_sim
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#
!# 2.) hydro_calcFluxGalNoBdr3d_sim
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#        without assembling the symmetric boundary contribution
!#
!# 3.) hydro_calcFluxScDiss3d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting scalar artificial viscosities
!#
!# 4.) hydro_calcFluxScDissDiSp3d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting scalar artificial viscosities based on
!#        dimensional splitting approach
!#
!# 5.) hydro_calcFluxRoeDiss3d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities
!#
!# 6.) hydro_calcFluxRoeDissDiSp3d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities based on
!#        dimensional splitting approach
!#
!# 7.) hydro_calcFluxRusDiss3d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion
!#
!# 8.) hydro_calcFluxRusDissDiSp3d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion
!#
!# 9.) hydro_calcMatDiagMatD3d_sim
!#     -> Computes local matrix for diagonal entry
!#
!# 10.) hydro_calcMatDiag3d_sim
!#      -> Computes local matrix for diagonal entry
!#
!# 11.) hydro_calcMatGalMatD3d_sim
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 12.) hydro_calcMatGal3d_sim
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 13.) hydro_calcMatScDissMatD3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 14.) hydro_calcMatScDiss3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 15.) hydro_calcMatRoeDissMatD3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 16.) hydro_calcMatRoeDiss3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 17.) hydro_calcMatRusDissMatD3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov artificial viscosities
!#
!# 18.) hydro_calcMatRusDiss3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov flux artificial viscosities
!#
!# 19.) hydro_calcCharacteristics3d_sim
!#      -> Computes characteristic variables
!#
!# 20.) hydro_calcFluxFCTScalarDiss3d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting scalar artificial viscosities
!#
!# 21.) hydro_calcFluxFCTRoeDiss3d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting tensorial artificial viscosities
!#
!# 22.) hydro_calcFluxFCTRusanov3d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting the Rusanov artificial viscosities
!#
!# 23.) hydro_trafoFluxDensity3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 24.) hydro_trafoDiffDensity3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 25.) hydro_trafoFluxEnergy3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 26.) hydro_trafoDiffEnergy3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 27.) hydro_trafoFluxPressure3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 28.) hydro_trafoFluxVelocity3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 29.) hydro_trafoDiffVelocity3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 30.) hydro_trafoFluxMomentum3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 31.) hydro_trafoDiffMomentum3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 32.) hydro_trafoFluxDenEng3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 33.) hydro_trafoDiffDenEng3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 34.) hydro_trafoFluxDenPre3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 35.) hydro_trafoDiffDenPre3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 36.) hydro_trafoFluxDenPreVel3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 37.) hydro_trafoDiffDenPreVel3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure 
!#         and the velocity
!#
!# 38.) hydro_calcBoundaryvalues3d
!#      -> Computes the boundary values for a given node
!#
!# 39.) hydro_hadaptCallbackScalar3d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 3D, whereby the vector is stored in interleave format
!#
!# 40.) hydro_hadaptCallbackBlock3d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 3D, whereby the vector is stored in block format
!#
!# </purpose>
!##############################################################################

module hydro_callback3d

#include "hydro.h"

  use collection
  use hydro_basic
  use flagship_callback
  use fsystem
  use genoutput
  use graph
  use groupfemsystem
  use hadaptaux
  use linearsystemblock
  use linearsystemscalar
  use problem
  use solveraux
  use storage

  implicit none

  private
  public :: hydro_calcFluxGal3d_sim
  public :: hydro_calcFluxGalNoBdr3d_sim
  public :: hydro_calcFluxScDiss3d_sim
  public :: hydro_calcFluxScDissDiSp3d_sim
  public :: hydro_calcFluxRoeDiss3d_sim
  public :: hydro_calcFluxRoeDissDiSp3d_sim
  public :: hydro_calcFluxRusDiss3d_sim
  public :: hydro_calcFluxRusDissDiSp3d_sim
  public :: hydro_calcMatDiagMatD3d_sim
  public :: hydro_calcMatDiag3d_sim
  public :: hydro_calcMatGalMatD3d_sim
  public :: hydro_calcMatGal3d_sim
  public :: hydro_calcMatScDissMatD3d_sim
  public :: hydro_calcMatScDiss3d_sim
  public :: hydro_calcMatRoeDissMatD3d_sim
  public :: hydro_calcMatRoeDiss3d_sim
  public :: hydro_calcMatRusDissMatD3d_sim
  public :: hydro_calcMatRusDiss3d_sim
  public :: hydro_calcCharacteristics3d_sim
  public :: hydro_calcFluxFCTScalarDiss3d
  public :: hydro_calcFluxFCTRoeDiss3d
  public :: hydro_calcFluxFCTRusanov3d
  public :: hydro_trafoFluxDensity3d_sim
  public :: hydro_trafoFluxEnergy3d_sim
  public :: hydro_trafoFluxPressure3d_sim
  public :: hydro_trafoFluxVelocity3d_sim
  public :: hydro_trafoFluxMomentum3d_sim
  public :: hydro_trafoFluxDenEng3d_sim
  public :: hydro_trafoFluxDenPre3d_sim
  public :: hydro_trafoFluxDenPreVel3d_sim
  public :: hydro_trafoDiffDensity3d_sim
  public :: hydro_trafoDiffEnergy3d_sim
  public :: hydro_trafoDiffPressure3d_sim
  public :: hydro_trafoDiffVelocity3d_sim
  public :: hydro_trafoDiffMomentum3d_sim
  public :: hydro_trafoDiffDenEng3d_sim
  public :: hydro_trafoDiffDenPre3d_sim
  public :: hydro_trafoDiffDenPreVel3d_sim
  public :: hydro_calcBoundaryvalues3d
  public :: hydro_hadaptCallbackScalar3d
  public :: hydro_hadaptCallbackBlock3d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGal3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the standard
    ! Galerkin discretisation in 3D.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

    ! local variables
#ifdef HYDRO_USE_IBP
    real(DP), dimension(NVAR3D) :: dF1_i,dF2_i,dF1_j,dF2_j,dF3_i,dF3_j
#else
    real(DP), dimension(NVAR3D) :: dF1_ij,dF2_ij,dF3_ij
#endif
    real(DP) :: pi,pj,ui,vi,uj,vj,wi,wj
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \        / rho*v         \         / rho*w         \
      !      | rho*u*u + p   |        | rho*u*v       |         | rho*u*w       |
      ! F1 = | rho*v*u       |   F2 = | rho*v*v + p   |   F3 =  | rho*v*w       |
      !      | rho*w*u       |        | rho*w*v       |         | rho*w*w + p   |
      !      \ rho*E*u + p*u /        \ rho*E*v + p*v /         \ rho*E*w + p*w /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui
      dF1_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui
      dF1_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi*ui

      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj*uj

      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi
      dF2_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi*vi

      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj*vj

      ! Compute fluxes for z-direction
      dF2_i(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi
      dF2_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi*wi

      dF2_j(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj
      dF2_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj*wj

      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*dF2_j+&
                                         DmatrixCoeffsAtEdge(3,2,idx)*dF3_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*dF2_i-&
                                         DmatrixCoeffsAtEdge(3,1,idx)*dF3_i )
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj*uj)

      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj)
      dF2_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj*vj)

      ! Compute flux difference for z-direction
      dF2_ij(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi)-&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj)
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi*wi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj*wj)

      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij+&
                                          DmatrixCoeffsAtEdge(3,1,idx)*dF3_ij)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij+&
                                          DmatrixCoeffsAtEdge(3,2,idx)*dF3_ij)
#endif
    end do

  end subroutine hydro_calcFluxGal3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGalNoBdr3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the TVD
    ! discretisation in 3D. The symmetric boundary contributions
    ! are neglected and incorporated in the antidiffusive flux.
    ! Hence, this is simply the standard Galerkin flux for the
    ! skew-symmetric internal contributions.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

  ! local variables
  real(DP), dimension(NVAR3D) :: dF1_ij,dF2_ij,dF3_ij
  real(DP) :: pi,pj,ui,vi,uj,vj,wi,wj
  integer :: idx
  
    
    do idx = 1, size(DfluxesAtEdge,3)

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \        / rho*v         \         / rho*w         \
      !      | rho*u*u + p   |        | rho*u*v       |         | rho*u*w       |
      ! F1 = | rho*v*u       |   F2 = | rho*v*v + p   |   F3 =  | rho*v*w       |
      !      | rho*w*u       |        | rho*w*v       |         | rho*w*w + p   |
      !      \ rho*E*u + p*u /        \ rho*E*v + p*v /         \ rho*E*w + p*w /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj*uj)

      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj)
      dF2_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj*vj)

      ! Compute flux difference for z-direction
      dF2_ij(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi)-&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj)
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi*wi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj*wj)

      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale *&
          (0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))*dF1_ij+&
           0.5*(DmatrixCoeffsAtEdge(2,1,idx)-DmatrixCoeffsAtEdge(2,2,idx))*dF2_ij+&
           0.5*(DmatrixCoeffsAtEdge(3,1,idx)-DmatrixCoeffsAtEdge(3,2,idx))*dF3_ij)
      DfluxesAtEdge(:,2,idx) = DfluxesAtEdge(:,1,idx)
    end do

  end subroutine hydro_calcFluxGalNoBdr3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxScDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using scalar dissipation.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

    ! local variables
#ifdef HYDRO_USE_IBP
    real(DP), dimension(NVAR3D) :: dF1_i,dF2_i,dF1_j,dF2_j,dF3_i,dF3_j
#else
    real(DP), dimension(NVAR3D) :: dF1_ij,dF2_ij,dF3_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: pi,pj,ui,vi,uj,vj,wi,wj,d_ij,H_ij,q_ij,u_ij,v_ij,w_ij,anorm,vel_ij,c_ij,aux
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \        / rho*v         \         / rho*w         \
      !      | rho*u*u + p   |        | rho*u*v       |         | rho*u*w       |
      ! F1 = | rho*v*u       |   F2 = | rho*v*v + p   |   F3 =  | rho*v*w       |
      !      | rho*w*u       |        | rho*w*v       |         | rho*w*w + p   |
      !      \ rho*E*u + p*u /        \ rho*E*v + p*v /         \ rho*E*w + p*w /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui
      dF1_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui
      dF1_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi*ui

      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj*uj

      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi
      dF2_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi*vi

      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj*vj

      ! Compute fluxes for z-direction
      dF2_i(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi
      dF2_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi*wi

      dF2_j(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj
      dF2_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj*wj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj*uj)

      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj)
      dF2_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj*vj)

      ! Compute flux difference for z-direction
      dF2_ij(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi)-&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj)
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi*wi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj*wj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral radius
      !-------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      w_ij = ROE_MEAN_VALUE(wi,wj,aux)
      H_ij = ROE_MEAN_VALUE(MYNEWLINE
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)+pi)/MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),MYNEWLINE
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)+pj)/MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx), aux)

      ! Compute auxiliary variables
      vel_ij = u_ij*a(1) + v_ij*a(2) + w_ij*a(3)
      q_ij   = 0.5*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      c_ij   = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif

      ! Scalar dissipation
      d_ij = abs(vel_ij) + anorm*c_ij

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*dF2_j+&
                                         DmatrixCoeffsAtEdge(3,2,idx)*dF3_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*dF2_i-&
                                         DmatrixCoeffsAtEdge(3,1,idx)*dF3_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij+&
                                          DmatrixCoeffsAtEdge(3,1,idx)*dF3_ij+ Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij+&
                                          DmatrixCoeffsAtEdge(3,2,idx)*dF3_ij+ Diff)
#endif
    end do

  end subroutine hydro_calcFluxScDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxScDissDiSp3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using scalar dissipation,
    ! whereby dimensional splitting is employed.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

    ! local variables
#ifdef HYDRO_USE_IBP
    real(DP), dimension(NVAR3D) :: dF1_i,dF2_i,dF1_j,dF2_j,dF3_i,dF3_j
#else
    real(DP), dimension(NVAR3D) :: dF1_ij,dF2_ij,dF3_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: pi,pj,ui,vi,uj,vj,wi,wj,d_ij,H_ij,q_ij,u_ij,v_ij,w_ij,aux,c_ij
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \        / rho*v         \         / rho*w         \
      !      | rho*u*u + p   |        | rho*u*v       |         | rho*u*w       |
      ! F1 = | rho*v*u       |   F2 = | rho*v*v + p   |   F3 =  | rho*v*w       |
      !      | rho*w*u       |        | rho*w*v       |         | rho*w*w + p   |
      !      \ rho*E*u + p*u /        \ rho*E*v + p*v /         \ rho*E*w + p*w /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui
      dF1_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui
      dF1_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi*ui

      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj*uj

      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi
      dF2_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi*vi

      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj*vj

      ! Compute fluxes for z-direction
      dF2_i(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi
      dF2_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi*wi

      dF2_j(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj
      dF2_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj*wj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj*uj)

      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj)
      dF2_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj*vj)

      ! Compute flux difference for z-direction
      dF2_ij(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi)-&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj)
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi*wi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj*wj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral radius
      !-------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      
      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      w_ij = ROE_MEAN_VALUE(vi,vj,aux)
      H_ij = ROE_MEAN_VALUE(MYNEWLINE
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)+pi)/MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),MYNEWLINE
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)+pj)/MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx), aux)

      ! Compute auxiliary variable
      q_ij = 0.5*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Scalar dissipation for x- and y-direction
      d_ij = ( abs(a(1)*u_ij) + abs(a(1))*c_ij +&
               abs(a(2)*v_ij) + abs(a(2))*c_ij +&
               abs(a(3)*w_ij) + abs(a(3))*c_ij )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*dF2_j+&
                                         DmatrixCoeffsAtEdge(3,2,idx)*dF3_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*dF2_i-&
                                         DmatrixCoeffsAtEdge(3,1,idx)*dF3_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij+&
                                          DmatrixCoeffsAtEdge(3,1,idx)*dF3_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij+&
                                          DmatrixCoeffsAtEdge(3,2,idx)*dF3_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxScDissDiSp3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRoeDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using tensorial dissipation.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

    ! local variables
#ifdef HYDRO_USE_IBP
    real(DP), dimension(NVAR3D) :: dF1_i,dF2_i,dF1_j,dF2_j,dF3_i,dF3_j
#else
    real(DP), dimension(NVAR3D) :: dF1_ij,dF2_ij,dF3_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: pi,pj,ui,vi,uj,vj,wi,wj,u_ij,v_ij,w_ij,H_ij,q_ij,c_ij,c2_ij,vel_ij
    real(DP) :: aux,aux1,aux2,anorm
    real(DP) :: l1,l2,l3,l4,l5,w1,w2,w3,w4,w5
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \        / rho*v         \         / rho*w         \
      !      | rho*u*u + p   |        | rho*u*v       |         | rho*u*w       |
      ! F1 = | rho*v*u       |   F2 = | rho*v*v + p   |   F3 =  | rho*v*w       |
      !      | rho*w*u       |        | rho*w*v       |         | rho*w*w + p   |
      !      \ rho*E*u + p*u /        \ rho*E*v + p*v /         \ rho*E*w + p*w /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui
      dF1_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui
      dF1_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi*ui

      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj*uj

      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi
      dF2_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi*vi

      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj*vj

      ! Compute fluxes for z-direction
      dF2_i(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi
      dF2_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi*wi

      dF2_j(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj
      dF2_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj*wj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj*uj)

      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj)
      dF2_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj*vj)

      ! Compute flux difference for z-direction
      dF2_ij(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi)-&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj)
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi*wi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj*wj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor by Roe
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      if (anorm .gt. SYS_EPSREAL) then
        
        ! Normalize the skew-symmetric coefficient
        a = a/anorm
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)+pi)/MYNEWLINE
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)+pj)/MYNEWLINE
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx), aux)

        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2) + w_ij*a(3)
        q_ij   = 0.5*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

        ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
        c2_ij  = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij   = sqrt(c2_ij)

        ! Compute eigenvalues
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)
        l5 = abs(vel_ij)

        ! Compute solution difference U_j-U_i
        Diff = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = (GAMMA-1.0)/2.0/c2_ij*(q_ij*Diff(1)-u_ij*Diff(2)-v_ij*Diff(3)-w_ij*Diff(4)+Diff(5))
        aux2 = 0.5*(vel_ij*Diff(1)-a(1)*Diff(2)-a(2)*Diff(3)-a(3)*Diff(4))/c_ij

        ! Get the dimension with largest coefficient
        select case(maxloc(a,1))
        case(1)
          ! Compute characteristic variables multiplied by the corresponding eigenvalue
          w1 = l1 * (aux1 + aux2)
          w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/c2_ij)*Diff(1)+&
                     (GAMMA-1.0)*(u_ij*Diff(2)+v_ij*Diff(3)+w_ij*Diff(4)-Diff(5))/c2_ij)
          w3 = l3 * (aux1 - aux2)
          w4 = l4 * ( (v_ij-vel_ij*a(2))*Diff(1)+a(1)*a(2)*Diff(2)+&
                      (a(2)*a(2)-1.0)*Diff(3)+a(2)*a(3)*Diff(4) )/a(1)
          w5 = l5 * ( (vel_ij*a(3)-w_ij)*Diff(1)-a(1)*a(3)*Diff(2)-&
                      a(2)*a(3)*Diff(3)+(1.0-a(3)*a(3))*Diff(4) )/a(1)

          ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
          Diff(1) = anorm * ( w1 + w2 + w3 + w4 )
          Diff(2) = anorm * ( (u_ij-c_ij*a(1))*w1 + u_ij*w2 + (u_ij+c_ij*a(1))*w3 + a(2)*w4 - a(3)*w5 )
          Diff(3) = anorm * ( (v_ij-c_ij*a(1))*w1 + v_ij*w2 + (v_ij+c_ij*a(1))*w3 - a(1)*w4 )
          Diff(4) = anorm * ( (w_ij-c_ij*a(1))*w1 + w_ij*w2 + (w_ij+c_ij*a(1))*w3 + a(2)*w5)
          Diff(5) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3&
                              + (u_ij*a(2)-v_ij*a(1))*w4 + (w_ij*a(1)-u_ij*a(3))*w5 )
        case(2)
          ! Compute characteristic variables multiplied by the corresponding eigenvalue
          w1 = l1 * (aux1 + aux2)
          w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/c2_ij)*Diff(1)+&
                     (GAMMA-1.0)*(u_ij*Diff(2)+v_ij*Diff(3)+w_ij*Diff(4)-Diff(5))/c2_ij)
          w3 = l3 * (aux1 - aux2)
          w4 = l4 * ( (vel_ij*a(1)-u_ij)*Diff(1)+(1.0-a(1)*a(1))*Diff(2)-&
                      a(1)*a(2)*Diff(3)-a(1)*a(3)*Diff(4) )/a(2)
          w5 = l5 * ( (w_ij-vel_ij*a(3))*Diff(1)+a(1)*a(3)*Diff(2)+&
                      a(2)*a(3)*Diff(3)+(a(3)*a(3)-1.0)*Diff(4) )/a(2)

          ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
          Diff(1) = anorm * ( w1 + w2 + w3 + w4 )
          Diff(2) = anorm * ( (u_ij-c_ij*a(1))*w1 + u_ij*w2 + (u_ij+c_ij*a(1))*w3 + a(2)*w4 )
          Diff(3) = anorm * ( (v_ij-c_ij*a(1))*w1 + v_ij*w2 + (v_ij+c_ij*a(1))*w3 - a(1)*w4 + a(3)*w5 )
          Diff(4) = anorm * ( (w_ij-c_ij*a(1))*w1 + w_ij*w2 + (w_ij+c_ij*a(1))*w3 - a(2)*w5)
          Diff(5) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3&
                              + (u_ij*a(2)-v_ij*a(1))*w4 + (v_ij*a(3)-w_ij*a(2))*w5 )
        case(3)
          ! Compute characteristic variables multiplied by the corresponding eigenvalue
          w1 = l1 * (aux1 + aux2)
          w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/c2_ij)*Diff(1)+&
                     (GAMMA-1.0)*(u_ij*Diff(2)+v_ij*Diff(3)+w_ij*Diff(4)-Diff(5))/c2_ij)
          w3 = l3 * (aux1 - aux2)
          w4 = l4 * ( (u_ij-vel_ij*a(1))*Diff(1)+(a(1)*a(1)-11.0)*Diff(2)+&
                      a(1)*a(2)*Diff(3)+a(1)*a(3)*Diff(4) )/a(3)
          w5 = l5 * ( (vel_ij*a(2)-v_ij)*Diff(1)-a(1)*a(2)*Diff(2)+&
                      (1.0-a(2)*a(2))*Diff(3)-a(2)*a(3)*Diff(4) )/a(3)

          ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
          Diff(1) = anorm * ( w1 + w2 + w3 + w4 )
          Diff(2) = anorm * ( (u_ij-c_ij*a(1))*w1 + u_ij*w2 + (u_ij+c_ij*a(1))*w3 - a(3)*w4 )
          Diff(3) = anorm * ( (v_ij-c_ij*a(1))*w1 + v_ij*w2 + (v_ij+c_ij*a(1))*w3 + a(3)*w5 )
          Diff(4) = anorm * ( (w_ij-c_ij*a(1))*w1 + w_ij*w2 + (w_ij+c_ij*a(1))*w3 + a(1)*w4 - a(2)*w5)
          Diff(5) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3&
                              + (w_ij*a(1)-u_ij*a(3))*w4 + (v_ij*a(3)-w_ij*a(2))*w5 )
        end select

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*dF2_j+&
                                           DmatrixCoeffsAtEdge(3,2,idx)*dF3_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*dF2_i-&
                                           DmatrixCoeffsAtEdge(3,1,idx)*dF3_i + Diff)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij+&
                                            DmatrixCoeffsAtEdge(3,1,idx)*dF3_ij + Diff)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij+&
                                            DmatrixCoeffsAtEdge(3,2,idx)*dF3_ij + Diff)
#endif
      else
        
#ifdef HYDRO_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*dF2_j+&
                                           DmatrixCoeffsAtEdge(3,2,idx)*dF3_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*dF2_i-&
                                           DmatrixCoeffsAtEdge(3,1,idx)*dF3_i)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij+&
                                            DmatrixCoeffsAtEdge(3,1,idx)*dF3_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij+&
                                            DmatrixCoeffsAtEdge(3,2,idx)*dF3_ij)
#endif
      end if
    end do

  end subroutine hydro_calcFluxRoeDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRoeDissDiSp3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using tensorial dissipation,
    ! whereby dimensional splitting is employed.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

    ! local variables
#ifdef HYDRO_USE_IBP
    real(DP), dimension(NVAR3D) :: dF1_i,dF2_i,dF1_j,dF2_j,dF3_i,dF3_j
#else
    real(DP), dimension(NVAR3D) :: dF1_ij,dF2_ij,dF3_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff1,Diff2,Diff3
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: pi,pj,ui,vi,uj,vj,wi,wj,u_ij,v_ij,w_ij,H_ij,q_ij,c_ij,c2_ij,vel_ij
    real(DP) :: aux,aux1,aux2,anorm
    real(DP) :: l1,l2,l3,l4,l5,w1,w2,w3,w4,w5
    integer :: idx
    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \        / rho*v         \         / rho*w         \
      !      | rho*u*u + p   |        | rho*u*v       |         | rho*u*w       |
      ! F1 = | rho*v*u       |   F2 = | rho*v*v + p   |   F3 =  | rho*v*w       |
      !      | rho*w*u       |        | rho*w*v       |         | rho*w*w + p   |
      !      \ rho*E*u + p*u /        \ rho*E*v + p*v /         \ rho*E*w + p*w /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui
      dF1_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui
      dF1_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi*ui

      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj*uj

      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi
      dF2_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi*vi

      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj*vj

      ! Compute fluxes for z-direction
      dF2_i(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi
      dF2_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi*wi

      dF2_j(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj
      dF2_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj*wj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj*uj)

      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj)
      dF2_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj*vj)

      ! Compute flux difference for z-direction
      dF2_ij(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi)-&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj)
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi*wi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj*wj)
#endif
      
      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor by Roe
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      if (anorm .gt. SYS_EPSREAL) then
        
        ! Normalize the skew-symmetric coefficient
        a = abs(a)
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)+pi)/MYNEWLINE
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)+pj)/MYNEWLINE
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx), aux)

        ! Compute auxiliary variables
        q_ij = 0.5*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

        ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
        c2_ij  = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij   = sqrt(c2_ij)

        ! ***      ***
        ! *** TODO ***
        ! ***      ***

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*dF2_j+&
                                           DmatrixCoeffsAtEdge(3,2,idx)*dF3_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*dF2_i-&
                                           DmatrixCoeffsAtEdge(3,1,idx)*dF3_i + Diff1+Diff2+Diff3)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij+&
                                            DmatrixCoeffsAtEdge(3,1,idx)*dF3_ij + Diff1+Diff2+Diff3)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij+&
                                            DmatrixCoeffsAtEdge(3,2,idx)*dF3_ij + Diff1+Diff2+Diff3)
#endif
      else
        
#ifdef HYDRO_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*dF2_j+&
                                           DmatrixCoeffsAtEdge(3,2,idx)*dF3_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*dF2_i-&
                                           DmatrixCoeffsAtEdge(3,1,idx)*dF3_i)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij+&
                                            DmatrixCoeffsAtEdge(3,1,idx)*dF3_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij+&
                                            DmatrixCoeffsAtEdge(3,2,idx)*dF3_ij)
#endif
      end if
    end do

  end subroutine hydro_calcFluxRoeDissDiSp3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRusDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using the Rusanov dissipation.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

    ! local variables
#ifdef HYDRO_USE_IBP
    real(DP), dimension(NVAR3D) :: dF1_i,dF2_i,dF1_j,dF2_j,dF3_i,dF3_j
#else
    real(DP), dimension(NVAR3D) :: dF1_ij,dF2_ij,dF3_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP) :: pi,pj,ui,vi,uj,vj,wi,wj,d_ij,ci,cj,Ei,Ej
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \        / rho*v         \         / rho*w         \
      !      | rho*u*u + p   |        | rho*u*v       |         | rho*u*w       |
      ! F1 = | rho*v*u       |   F2 = | rho*v*v + p   |   F3 =  | rho*v*w       |
      !      | rho*w*u       |        | rho*w*v       |         | rho*w*w + p   |
      !      \ rho*E*u + p*u /        \ rho*E*v + p*v /         \ rho*E*w + p*w /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui
      dF1_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui
      dF1_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi*ui

      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj*uj

      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi
      dF2_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi*vi

      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj*vj

      ! Compute fluxes for z-direction
      dF2_i(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi
      dF2_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi*wi

      dF2_j(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj
      dF2_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj*wj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj*uj)

      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj)
      dF2_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj*vj)

      ! Compute flux difference for z-direction
      dF2_ij(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi)-&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj)
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi*wi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj*wj)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov type
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*(ui*ui+vi*vi+wi*wi)), SYS_EPSREAL))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*(uj*uj+vj*vj+wj*wj)), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Scalar dissipation for the Rusanov flux
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj+&
                      DmatrixCoeffsAtEdge(2,1,idx)*vj+&
                      DmatrixCoeffsAtEdge(3,1,idx)*wj)+&
                 sqrt(DmatrixCoeffsAtEdge(1,1,idx)**2+&
                      DmatrixCoeffsAtEdge(2,1,idx)**2+&
                      DmatrixCoeffsAtEdge(3,1,idx)**2)*cj,&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui+&
                      DmatrixCoeffsAtEdge(2,2,idx)*vi+&
                      DmatrixCoeffsAtEdge(3,2,idx)*wi)+&
                 sqrt(DmatrixCoeffsAtEdge(1,2,idx)**2+&
                      DmatrixCoeffsAtEdge(2,2,idx)**2+&
                      DmatrixCoeffsAtEdge(3,2,idx)**2)*ci )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------
      
#ifdef HYDRO_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*dF2_j+&
                                         DmatrixCoeffsAtEdge(3,2,idx)*dF3_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*dF2_i-&
                                         DmatrixCoeffsAtEdge(3,1,idx)*dF3_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij+&
                                          DmatrixCoeffsAtEdge(3,1,idx)*dF3_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij+&
                                          DmatrixCoeffsAtEdge(3,2,idx)*dF3_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxRusDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRusDissDiSp3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using the Rusanov dissipation,
    ! whereby dimensional splitting is employed.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

    ! local variables
#ifdef HYDRO_USE_IBP
    real(DP), dimension(NVAR3D) :: dF1_i,dF2_i,dF1_j,dF2_j,dF3_i,dF3_j
#else
    real(DP), dimension(NVAR3D) :: dF1_ij,dF2_ij,dF3_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP) :: pi,pj,ui,vi,uj,vj,wi,wj,d_ij,ci,cj,Ei,Ej
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \        / rho*v         \         / rho*w         \
      !      | rho*u*u + p   |        | rho*u*v       |         | rho*u*w       |
      ! F1 = | rho*v*u       |   F2 = | rho*v*v + p   |   F3 =  | rho*v*w       |
      !      | rho*w*u       |        | rho*w*v       |         | rho*w*w + p   |
      !      \ rho*E*u + p*u /        \ rho*E*v + p*v /         \ rho*E*w + p*w /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui
      dF1_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui
      dF1_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi*ui

      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj*uj

      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi
      dF2_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi*vi

      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj*vj

      ! Compute fluxes for z-direction
      dF2_i(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi
      dF2_i(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi*wi

      dF2_j(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj
      dF2_j(5) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj*wj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj
      dF1_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj + pj*uj)

      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj)
      dF2_ij(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj + pj*vj)

      ! Compute flux difference for z-direction
      dF2_ij(1) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)-&
                  Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi)-&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj)
      dF2_ij(5) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi + pi*wi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj + pj*wj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov type
      !-------------------------------------------------------------------------

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*(ui*ui+vi*vi+wi*wi)), SYS_EPSREAL))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*(uj*uj+vj*vj+wj*wj)), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Scalar dissipation
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj)+&
                  abs(DmatrixCoeffsAtEdge(1,1,idx))*cj,&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui)+&
                  abs(DmatrixCoeffsAtEdge(1,2,idx))*ci )&
           + max( abs(DmatrixCoeffsAtEdge(2,1,idx)*vj)+&
                  abs(DmatrixCoeffsAtEdge(2,1,idx))*cj,&
                  abs(DmatrixCoeffsAtEdge(2,2,idx)*vi)+&
                  abs(DmatrixCoeffsAtEdge(2,2,idx))*ci )&
           + max( abs(DmatrixCoeffsAtEdge(3,1,idx)*wj)+&
                  abs(DmatrixCoeffsAtEdge(3,1,idx))*cj,&
                  abs(DmatrixCoeffsAtEdge(3,2,idx)*wi)+&
                  abs(DmatrixCoeffsAtEdge(3,2,idx))*ci )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))
      
      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*dF2_j+&
                                         DmatrixCoeffsAtEdge(3,2,idx)*dF3_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*dF2_i-&
                                         DmatrixCoeffsAtEdge(3,1,idx)*dF3_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij+&
                                          DmatrixCoeffsAtEdge(3,1,idx)*dF3_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij+&
                                          DmatrixCoeffsAtEdge(3,2,idx)*dF3_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxRusDissDiSp3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiagMatD3d_sim(DdataAtNode, DmatrixCoeffsAtNode,&
      IverticesAtNode, dscale, DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 3D
!</description>

!<input>
  ! Nodal solution values for all nodes under consideration
  real(DP), dimension(:,:), intent(in) :: DdataAtNode

  ! Entries of the coefficient matrices for all nodes under consideration
  real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode

  ! Numbers of vertices and matrix entries for all nodes under consideration
  integer, dimension(:,:), intent(in) :: IverticesAtNode

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all nodes under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtNode
!</output>
!</subroutine>

  end subroutine hydro_calcMatDiagMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiag3d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 3D
!</description>

!<input>
  ! Nodal solution values for all nodes under consideration
  real(DP), dimension(:,:), intent(in) :: DdataAtNode

  ! Entries of the coefficient matrices for all nodes under consideration
  real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode

  ! Numbers of vertices and matrix entries for all nodes under consideration
  integer, dimension(:,:), intent(in) :: IverticesAtNode

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all nodes under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtNode
!</output>
!</subroutine>

  end subroutine hydro_calcMatDiag3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGalMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine hydro_calcMatGalMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGal3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine hydro_calcMatGal3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDissMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies scalar artificial viscosities in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine hydro_calcMatScDissMatD3d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDiss3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies scalar artificial viscosities in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine hydro_calcMatScDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDissMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine hydro_calcMatRoeDissMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDiss3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine hydro_calcMatRoeDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDissMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine hydro_calcMatRusDissMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDiss3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine hydro_calcMatRusDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcCharacteristics3d_sim(Dweight, DdataAtEdge,&
      DcharVariablesAtEdge, DeigenvaluesAtEdge,&
      DrightEigenvectorsAtEdge, DleftEigenvectorsAtEdge, rcollection)

!<description>
    ! This subroutine computes the characteristic variables in 3D
!</description>

!<input>
    ! Weighting coefficient for wave-decomposition
    real(DP), dimension(:), intent(in)  :: Dweight
    
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! OPTIONAL: Characteristic variables for all edges under consideration
    !   DIMENSION(nvar,nedge)
    ! with nvar the number of variables at each edge
    real(DP), dimension(:,:), intent(out), optional :: DcharVariablesAtEdge
    
    ! OPTIONAL: Eigenvalues for all edges under consideration
    !   DIMENSION(nvar,nedge)
    ! with nvar the number of variables at each edge
    real(DP), dimension(:,:), intent(out), optional :: DeigenvaluesAtEdge
    
    ! OPTIONAL: Matrices of left eigenvectors for all edges under consideration
    !   DIMENSION(nvar*nvar,nedge)
    ! with nvar the number of variables at each edge
    real(DP), dimension(:,:), intent(out), optional :: DleftEigenvectorsAtEdge
    
    ! OPTIONAL: Matrices of right eigenvectors for all edges under consideration
    !   DIMENSION(nvar*nvar,nedge)
    ! with nvar the number of variables at each edge
    real(DP), dimension(:,:), intent(out), optional :: DrightEigenvectorsAtEdge
!</output>
!</subroutine>

  end subroutine hydro_calcCharacteristics3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTScalarDiss3d(&
      U1_i, U1_j, U2_i, U2_j, C_ij, C_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 3D using scalar dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficients
    real(DP), intent(in) :: dscale1,dscale2

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! raw antidiffusive flux
    real(DP), dimension(:), intent(out) :: F_ij
!</output>
!</subroutine>

  end subroutine hydro_calcFluxFCTScalarDiss3d

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRoeDiss3d(&
      U1_i, U1_j, U2_i, U2_j, C_ij, C_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 3D using tensorial dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficients
    real(DP), intent(in) :: dscale1,dscale2

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! raw antidiffusive flux
    real(DP), dimension(:), intent(out) :: F_ij
!</output>
!</subroutine>

  end subroutine hydro_calcFluxFCTRoeDiss3d

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRusanov3d(&
      U1_i, U1_j, U2_i, U2_j, C_ij, C_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 3D using the Rusanov dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficients
    real(DP), intent(in) :: dscale1,dscale2

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! raw antidiffusive flux
    real(DP), dimension(:), intent(out) :: F_ij
!</output>
!</subroutine>

  end subroutine hydro_calcFluxFCTRusanov3d

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDensity3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
    end do

  end subroutine hydro_trafoFluxDensity3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDensity3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed density difference
      DtransformedDataAtEdge(1,idx) =&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
    end do

  end subroutine hydro_trafoDiffDensity3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxEnergy3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the energy in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed total energy fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
    end do

  end subroutine hydro_trafoFluxEnergy3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffEnergy3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the energy in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed total density difference
      DtransformedDataAtEdge(1,idx) =&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
    end do

  end subroutine hydro_trafoDiffEnergy3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxPressure3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the pressure in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(1,1,idx) = (GAMMA-1.0)*&
          (0.5*(ui*ui+vi*vi)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vi*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wi*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))
      DtransformedFluxesAtEdge(1,2,idx) =-(GAMMA-1.0)*&
          (0.5*(uj*uj+vj*vj)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vj*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wj*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do

  end subroutine hydro_trafoFluxPressure3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffPressure3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the pressure in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(1,idx) =&
          PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
    end do

  end subroutine hydro_trafoDiffPressure3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxVelocity3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the x-, y- and z-velocity
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      
      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(1,1,idx) =&
          (X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -(X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Transformed velocity fluxes in y-direction
      DtransformedFluxesAtEdge(2,1,idx) =&
          (Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -(Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Transformed velocity fluxes in z-direction
      DtransformedFluxesAtEdge(3,1,idx) =&
          (Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      DtransformedFluxesAtEdge(3,2,idx) =&
         -(Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
    end do
    
  end subroutine hydro_trafoFluxVelocity3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffVelocity3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the x-, y- and z-velocity
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)

      ! Transformed velocity difference in x-direction
      DtransformedDataAtEdge(1,idx) =&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed velocity difference in y-direction
      DtransformedDataAtEdge(2,idx) =&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed velocity difference in z-direction
      DtransformedDataAtEdge(3,idx) =&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
    end do
    
  end subroutine hydro_trafoDiffVelocity3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxMomentum3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the x-, y- and z-momentum
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed momentum fluxes in x-direction
      DtransformedFluxesAtEdge(1,1,idx) =&
          X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)

      ! Transformed momentum fluxes in y-direction
      DtransformedFluxesAtEdge(2,1,idx) =&
          Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)

      ! Transformed momentum fluxes in z-direction
      DtransformedFluxesAtEdge(3,1,idx) =&
          Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(3,2,idx) =&
         -Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
    end do
    
  end subroutine hydro_trafoFluxMomentum3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffMomentum3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the x-, y- and z-momentum
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

     ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed momentum difference in x-direction
      DtransformedDataAtEdge(1,idx) =&
          X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed momentum difference in y-direction
      DtransformedDataAtEdge(2,idx) =&
          Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed momentum difference in z-direction
      DtransformedDataAtEdge(3,idx) =&
          Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
    end do
    
  end subroutine hydro_trafoDiffMomentum3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenEng3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)

      ! Transformed total energy fluxes
      DtransformedFluxesAtEdge(2,1,idx) =&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
    end do

  end subroutine hydro_trafoFluxDenEng3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenEng3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)

      ! Transformed density difference
      DtransformedDataAtEdge(1,idx) =&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed total energy difference
      DtransformedDataAtEdge(2,idx) =&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
    end do

  end subroutine hydro_trafoDiffDenEng3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenPre3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      
      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(2,1,idx) = (GAMMA-1.0)*&
          (0.5*(ui*ui+vi*vi)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vi*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wi*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))
      DtransformedFluxesAtEdge(2,2,idx) =-(GAMMA-1.0)*&
          (0.5*(uj*uj+vj*vj)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vj*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wj*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do

  end subroutine hydro_trafoFluxDenPre3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenPre3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)

      ! Transformed density difference
      DtransformedDataAtEdge(1,idx) =&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(2,idx) =&
          PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
    end do
    
  end subroutine hydro_trafoDiffDenPre3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenPreVel3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation
    ! of the given flux into primitive variables in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)

      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(2,1,idx) =&
          (X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -(X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Transformed velocity fluxes in y-direction
      DtransformedFluxesAtEdge(3,1,idx) =&
          (Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      DtransformedFluxesAtEdge(3,2,idx) =&
         -(Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Transformed velocity fluxes in z-direction
      DtransformedFluxesAtEdge(4,1,idx) =&
          (Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      DtransformedFluxesAtEdge(4,2,idx) =&
         -(Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(5,1,idx) =(GAMMA-1.0)*&
          (0.5*(ui*ui+vi*vi)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vi*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wi*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))
      DtransformedFluxesAtEdge(5,2,idx) =-(GAMMA-1.0)*&
          (0.5*(uj*uj+vj*vj)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vj*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wj*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do

  end subroutine hydro_trafoFluxDenPreVel3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenPreVel3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density, pressure and velocity in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed density difference
      DtransformedDataAtEdge(2,idx) =&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed velocity difference in x-direction
      DtransformedDataAtEdge(2,idx) =&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      
      ! Transformed velocity difference in y-direction
      DtransformedDataAtEdge(3,idx) =&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed velocity difference in z-direction
      DtransformedDataAtEdge(4,idx) =&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed pressure difference
      DtransformedDataAtEdge(5,idx) =&
          PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
    end do

  end subroutine hydro_trafoDiffDenPreVel3d_sim

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcBoundaryvalues3d(DbdrNormal, DpointNormal,&
      DbdrValue, ibdrCondType, Du, Du0, istatus)

!<description>
    ! This subroutine computes the boundary values for a given node in 3D
!</description>

!<input>
    ! normal vector at the boundary
    real(DP), dimension(:), intent(in) :: DbdrNormal

    ! normal vector at the point on the boundary
    real(DP), dimension(:), intent(in) :: DpointNormal

    ! evaluated boundary values
    real(DP), dimension(:), intent(in) :: DbdrValue

    ! initial solution from the previous time step
    real(DP), dimension(:), intent(in) :: Du0

    ! type of boundary condition
    integer, intent(in) :: ibdrCondType
!</input>

!<inputoutput>
    ! computed boundary values
    real(DP), dimension(:), intent(inout) :: Du

    ! OPTIONAL: status of the callback function
    integer, intent(inout), optional :: istatus
!</inputoutput>
!</subroutine>

  end subroutine hydro_calcBoundaryvalues3d

  !*****************************************************************************

!<subroutine>

  subroutine hydro_hadaptCallbackScalar3d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 3D. The solution vector is assumed
    ! to be store in scalar interleave format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolution
    real(DP), dimension(:), pointer, save :: p_Dsolution
    integer :: ivar


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve solution vector from colletion
      rsolution => rcollection%p_rvectorQuickAccess1

      ! Check if solution is stored in interleave format
      if (rsolution%nblocks .ne. 1) then
        call output_line('Vector is not in interleave format!',&
            OU_CLASS_WARNING,OU_MODE_STD,'hydro_hadaptCallbackScalar3d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      do ivar = 1, NVAR3D
        p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR3D+ivar) = &
            0.5_DP*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR3D+ivar)+&
                    p_Dsolution((rcollection%IquickAccess(3)-1)*NVAR3D+ivar))
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      do ivar = 1, NVAR3D
        p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR3D+ivar) = &
            0.25_DP*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR3D+ivar)+&
                     p_Dsolution((rcollection%IquickAccess(3)-1)*NVAR3D+ivar)+&
                     p_Dsolution((rcollection%IquickAccess(4)-1)*NVAR3D+ivar)+&
                     p_Dsolution((rcollection%IquickAccess(5)-1)*NVAR3D+ivar))
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        do ivar = 1, NVAR3D
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR3D+ivar) = &
              p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR3D+ivar)
        end do
      else
        do ivar = 1, NVAR3D
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR3D+ivar) = 0.0_DP
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)

    end select

  end subroutine hydro_hadaptCallbackScalar3d

  !*****************************************************************************

!<subroutine>

  subroutine hydro_hadaptCallbackBlock3d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 3D. The solution vector is assumed
    ! to be store in block format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolution
    real(DP), dimension(:), pointer, save :: p_Dsolution
    integer :: ivar,neq


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve solution vector from colletion
      rsolution => rcollection%p_rvectorQuickAccess1

      ! Check if solution is stored in interleave format
      if (rsolution%nblocks .ne. NVAR3D) then
        call output_line('Vector is not in block format!',&
            OU_CLASS_WARNING,OU_MODE_STD,'hydro_hadaptCallbackBlock3d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      neq = rsolution%NEQ/NVAR3D
      do ivar = 1, NVAR3D
        p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = &
            0.5_DP*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
                    p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(3)) )
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      neq = rsolution%NEQ/NVAR3D
      do ivar = 1, NVAR3D
        p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) =&
            0.25_DP*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
                     p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(3))+&
                     p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(4))+&
                     p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(5)) )
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        neq = rsolution%NEQ/NVAR3D
        do ivar = 1, NVAR3D
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = &
              p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))
        end do
      else
        neq = rsolution%NEQ/NVAR3D
        do ivar = 1, NVAR3D
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = 0.0_DP
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)

    end select

  end subroutine hydro_hadaptCallbackBlock3d

end module hydro_callback3d
