!##############################################################################
!# ****************************************************************************
!# <name> mhd_callback1d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible MHDequations in 1D.
!#
!# The following callback functions are available:
!#
!# 1.) mhd_calcFluxGal1d_sim
!#     -> Computes fluxes for standard Galerkin scheme
!#
!# 2.) mhd_calcFluxGalNoBdr1d_sim
!#     -> Computes fluxes for standard Galerkin scheme without
!#        assembling the symmetric boundary contribution
!#
!# 3.) mhd_calcFluxScDiss1d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting scalar artificial viscosities
!#
!# 4.) mhd_calcFluxRoeDiss1d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities
!#
!# 5.) mhd_calcFluxRusDiss1d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion
!#
!# 6.) mhd_calcMatDiagMatD1d_sim
!#     -> Computes local matrix for diagonal entry
!#
!# 7.) mhd_calcMatDiag1d_sim
!#     -> Computes local matrix for diagonal entry
!#
!# 8.) mhd_calcMatGalMatD1d_sim
!#     -> Computes local matrices for standard Galerkin scheme
!#
!# 9.) mhd_calcMatGal1d_sim
!#     -> Computes local matrices for standard Galerkin scheme
!#
!# 10.) mhd_calcMatScDissMatD1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 11.) mhd_calcMatScDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 12.) mhd_calcMatRoeDissMatD1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 13.) mhd_calcMatRoeDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 14.) mhd_calcMatRusDissMatD1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov artificial viscosities
!#
!# 15.) mhd_calcMatRusDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov flux artificial viscosities
!#
!# 16.) mhd_calcCharacteristics1d_sim
!#      -> Computes characteristic variables
!#
!# 17.) mhd_calcFluxFCTScalarDiss1d
!#      -> Computes fluxes for FCT algorithm
!#         adopting scalar artificial viscosities
!#
!# 18.) mhd_calcFluxFCTTensorDiss1d
!#      -> Computes fluxes for FCT algorithm
!#         adopting tensorial artificial viscosities
!#
!# 19.) mhd_calcFluxFCTRusanov1d
!#      -> Computes fluxes for FCT algorithm
!#         adopting the Rusanov artificial viscosities
!#
!# 20.) mhd_trafoFluxDensity1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 21.) mhd_trafoDiffDensity1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 22.) mhd_trafoFluxEnergy1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 23.) mhd_trafoDiffEnergy1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 24.) mhd_trafoFluxPressure1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 25.) mhd_trafoDiffPressure1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the pressure
!#
!# 26.) mhd_trafoFluxVelocity1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 27.) mhd_trafoDiffVelocity1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 28.) mhd_trafoFluxMomentum1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 29.) mhd_trafoDiffMomentum1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 30.) mhd_trafoFluxDenEng1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 31.) mhd_trafoDiffDenEng1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 32.) mhd_trafoFluxDenPre1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 33.) mhd_trafoDiffDenPre1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 34.) mhd_trafoFluxDenPreVel1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 35.) mhd_trafoDiffDenPreVel1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure 
!#         and the velocity
!#
!# 36.) mhd_calcBoundaryvalues1d
!#      -> Computes the boundary values for a given node
!#
!# 37.) mhd_hadaptCallbackScalar1d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 1D, whereby the vector is stored in interleave format
!#
!# 38.) mhd_hadaptCallbackBlock1d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 1D, whereby the vector is stored in block format
!#
!# </purpose>
!##############################################################################

module mhd_callback1d

#include "mhd.h"

  use boundarycondaux
  use collection
  use mhd_basic
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
  public :: mhd_calcFluxGal1d_sim
  public :: mhd_calcFluxGalNoBdr1d_sim
  public :: mhd_calcFluxScDiss1d_sim
  public :: mhd_calcFluxRoeDiss1d_sim
  public :: mhd_calcFluxRusDiss1d_sim
  public :: mhd_calcMatDiagMatD1d_sim
  public :: mhd_calcMatDiag1d_sim
  public :: mhd_calcMatGalMatD1d_sim
  public :: mhd_calcMatGal1d_sim
  public :: mhd_calcMatScDissMatD1d_sim
  public :: mhd_calcMatScDiss1d_sim
  public :: mhd_calcMatRoeDissMatD1d_sim
  public :: mhd_calcMatRoeDiss1d_sim
  public :: mhd_calcMatRusDissMatD1d_sim
  public :: mhd_calcMatRusDiss1d_sim
  public :: mhd_calcCharacteristics1d_sim
  public :: mhd_calcFluxFCTScalarDiss1d
  public :: mhd_calcFluxFCTTensorDiss1d
  public :: mhd_calcFluxFCTRusanov1d
  public :: mhd_trafoFluxDensity1d_sim
  public :: mhd_trafoFluxEnergy1d_sim
  public :: mhd_trafoFluxPressure1d_sim
  public :: mhd_trafoFluxVelocity1d_sim
  public :: mhd_trafoFluxMomentum1d_sim
  public :: mhd_trafoFluxDenEng1d_sim
  public :: mhd_trafoFluxDenPre1d_sim
  public :: mhd_trafoFluxDenPreVel1d_sim
  public :: mhd_trafoDiffDensity1d_sim
  public :: mhd_trafoDiffEnergy1d_sim
  public :: mhd_trafoDiffPressure1d_sim
  public :: mhd_trafoDiffVelocity1d_sim
  public :: mhd_trafoDiffMomentum1d_sim
  public :: mhd_trafoDiffDenEng1d_sim
  public :: mhd_trafoDiffDenPre1d_sim
  public :: mhd_trafoDiffDenPreVel1d_sim
  public :: mhd_calcBoundaryvalues1d
  public :: mhd_hadaptCallbackScalar1d
  public :: mhd_hadaptCallbackBlock1d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGal1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the standard Galerkin
    ! discretisation in 1D.
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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR1D) :: dF_i,dF_j
#else
    real(DP), dimension(NVAR1D) :: dF_ij
#endif
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx


    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !
      !     /   rho*u            \
      !     | rho*u*u + pT       |
      !     | rho*v*u - Bx*By    |
      ! F = | rho*w*u - Bx*Bz    |
      !     |    By*u - Bx*v     |
      !     |    Bz*u - Bx*w     |
      !     \ (rho*E+pT)*u -Bx*q /
      !
      ! pT = p + q, q = 1/(2*nu)*|B|^2
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj
      
#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      dF_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      dF_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui + pi
      dF_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      dF_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      dF_i(5) = Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi
      dF_i(6) = Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*wi
      dF_i(7) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx) + pi)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*qi
      
      dF_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      dF_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj + pj
      dF_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)
      dF_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)
      dF_j(5) = Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj
      dF_j(6) = Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj
      dF_j(7) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx) + pj)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*qj
      
      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF_i )
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      ! Compute flux difference for x-direction
      dF_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)-&
                 X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      dF_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui + pi) -&
                 (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj + pj)
      dF_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                  Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)) -&
                 (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                  Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))
      dF_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                  Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)) -&
                 (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                  Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))
      dF_ij(5) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi) -&
                 (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj)
      dF_ij(6) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj) -&
                 (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj)
      dF_ij(7) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx) + pi)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*qi) -&
                 ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx) + pj)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*qj)
                 
      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * DmatrixCoeffsAtEdge(1,1,idx)*dF_ij
      DfluxesAtEdge(:,2,idx) = -dscale * DmatrixCoeffsAtEdge(1,2,idx)*dF_ij      
#endif
    end do

  end subroutine mhd_calcFluxGal1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGalNoBdr1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the TVD discretisation
    ! in 1D. The symmetric boundary contributions are neglected and
    ! incorporated in the antidiffusive flux.  Hence, this is simply
    ! the standard Galerkin flux for the skew-symmetric internal
    ! contributions.
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
    real(DP), dimension(NVAR1D) :: dF_ij
    real(DP) ::  pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx


    do idx = 1, size(DfluxesAtEdge,3)

      !-------------------------------------------------------------------------
      !
      !     /   rho*u            \
      !     | rho*u*u + pT       |
      !     | rho*v*u - Bx*By    |
      ! F = | rho*w*u - Bx*Bz    |
      !     |    By*u - Bx*v     |
      !     |    Bz*u - Bx*w     |
      !     \ (rho*E+pT)*u -Bx*q /
      !
      ! pT = p + q, q = 1/(2*nu)*|B|^2
      !-------------------------------------------------------------------------
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj
      
      ! Compute flux difference for x-direction
      dF_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)-&
                 X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      dF_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui + pi) -&
                 (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj + pj)
      dF_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                  Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)) -&
                 (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                  Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))
      dF_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                  Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)) -&
                 (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                  Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))
      dF_ij(5) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi) -&
                 (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj)
      dF_ij(6) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj) -&
                 (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj)
      dF_ij(7) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx) + pi)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*qi) -&
                 ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx) + pj)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*qj)

      ! Assemble symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale *&
          0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))*dF_ij
      DfluxesAtEdge(:,2,idx) = DfluxesAtEdge(:,1,idx)
    end do

  end subroutine mhd_calcFluxGalNoBdr1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxScDiss1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 1D using scalar dissipation.
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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR1D) :: dF_i,dF_j
#else
    real(DP), dimension(NVAR1D) :: dF_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj,cf_ij
    real(DP) :: aux,u_ij,q_ij,H_ij,X_ij
    real(DP) :: aPow2_ij,bPow2_ij,bxPow2_ij,astPow2_ij
    real(DP) :: anorm,d_ij
    integer :: idx

    do idx = 1, size(DfluxesAtEdge,3)
          
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !
      !     /   rho*u            \
      !     | rho*u*u + pT       |
      !     | rho*v*u - Bx*By    |
      ! F = | rho*w*u - Bx*Bz    |
      !     |    By*u - Bx*v     |
      !     |    Bz*u - Bx*w     |
      !     \ (rho*E+pT)*u -Bx*q /
      !
      ! pT = p + q, q = 1/(2*nu)*|B|^2
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      dF_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      dF_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui + pi
      dF_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      dF_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      dF_i(5) = Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi
      dF_i(6) = Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*wi
      dF_i(7) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx) + pi)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*qi
      
      dF_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      dF_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj + pj
      dF_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)
      dF_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)
      dF_j(5) = Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj
      dF_j(6) = Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj
      dF_j(7) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx) + pj)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*qj
#else
      ! Compute flux difference for x-direction
      dF_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)-&
                 X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      dF_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui + pi) -&
                 (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj + pj)
      dF_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                  Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)) -&
                 (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                  Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))
      dF_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                  Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)) -&
                 (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                  Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))
      dF_ij(5) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi) -&
                 (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj)
      dF_ij(6) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj) -&
                 (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj)
      dF_ij(7) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx) + pi)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*qi) -&
                 ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx) + pj)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral radius
      !
      ! There are seven eigenvalues
      !
      !   u-cf, u-ca, u-cs, u, u+cs, u+ca, u+cf,
      !
      ! where u us the x-velocity component and ca, cs and cf are the
      ! velocities of th Alfveen waves, the slow and fast waves.
      !
      ! The largest in magnitude eigenvalue is |u|+cf
      !-------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))
      anorm = abs(a(1)) ! = sqrt(a(1)*a(1))

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(MYNEWLINE
              DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
              DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      H_ij = ROE_MEAN_VALUE(MYNEWLINE
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+MYNEWLINE
              TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))/MYNEWLINE
              DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+MYNEWLINE
              TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))/MYNEWLINE
              DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx), aux)
        
      ! Compute the square of the Roe-averaged speed of the Alfven waves.
      ! Note that left and right states are interchanged!
      bxPow2_ij = (ROE_MEAN_VALUE(MYNEWLINE
                    X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),MYNEWLINE
                    X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux))**2/MYNEWLINE
                  ROE_MEAN_VALUE(MYNEWLINE
                   DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),MYNEWLINE
                   DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),aux)
    
      ! Compute the density-averaged magnetic field
      X_ij = ((X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-MYNEWLINE
               X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2+MYNEWLINE
              (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-MYNEWLINE
               Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2+MYNEWLINE
              (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-MYNEWLINE
               Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2)/MYNEWLINE
              (2.0*(sqrt(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))+MYNEWLINE
                    sqrt(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))))

      ! Compute the square of the Roe-averaged magnetic field.
      ! Note that left and right states are interchanged!
      bPow2_ij = (ROE_MEAN_VALUE(MYNEWLINE
                   X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),MYNEWLINE
                   X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2+MYNEWLINE
                  ROE_MEAN_VALUE(MYNEWLINE
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),MYNEWLINE
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2+MYNEWLINE
                  ROE_MEAN_VALUE(MYNEWLINE
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),MYNEWLINE
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2)/MYNEWLINE
                  ROE_MEAN_VALUE(MYNEWLINE
                   DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),MYNEWLINE
                   DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),aux)

      ! Compute the magnitude of the Roe-averaged velocity
      q_ij = ROE_MEAN_VALUE(MYNEWLINE
              X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
              X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2+MYNEWLINE
             ROE_MEAN_VALUE(MYNEWLINE
              Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
              Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2+MYNEWLINE
             ROE_MEAN_VALUE(MYNEWLINE
              Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
              Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2

      ! Compute the Roe-averaged speed of sound
#ifdef THERMALLY_IDEAL_GAS
      aPow2_ij = (2.0-GAMMA)*X_ij + (GAMMA-1.0)*(H_ij-0.5*q_ij-bPow2_ij)
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Compute auxiliary variables
      astPow2_ij = aPow2_ij+bPow2_ij
      aux = sqrt(astPow2_ij**2-4.0*aPow2_ij*bxPow2_ij)
            
      ! Compute the Roe-averagred speed of the fast waves
      cf_ij = sqrt(0.5*(astPow2_ij+aux))

      ! Scalar dissipation
      d_ij = abs(u_ij*a(1)) + anorm*cf_ij
      
      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------
      
#ifdef MHD_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_ij + Diff)
#endif
    end do

  end subroutine mhd_calcFluxScDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRoeDiss1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 1D using tensorial dissipation.
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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR1D) :: dF_i,dF_j
#else
    real(DP), dimension(NVAR1D) :: dF_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj,ca_ij,cs_ij,cf_ij
    real(DP) :: aux,u_ij,q_ij,H_ij,X_ij
    real(DP) :: aPow2_ij,bPow2_ij,bxPow2_ij,astPow2_ij
    real(DP) :: l1,l2,l3,l4,l5,l6,l7,w1,w2,w3,w4,w5,w6,w7
    real(DP) :: anorm
    integer :: idx
    
    
    do idx = 1, size(DfluxesAtEdge,3)
          
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !
      !     /   rho*u            \
      !     | rho*u*u + pT       |
      !     | rho*v*u - Bx*By    |
      ! F = | rho*w*u - Bx*Bz    |
      !     |    By*u - Bx*v     |
      !     |    Bz*u - Bx*w     |
      !     \ (rho*E+pT)*u -Bx*q /
      !
      ! pT = p + q, q = 1/(2*nu)*|B|^2
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      dF_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      dF_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui + pi
      dF_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      dF_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      dF_i(5) = Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi
      dF_i(6) = Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*wi
      dF_i(7) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx) + pi)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*qi
      
      dF_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      dF_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj + pj
      dF_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)
      dF_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)
      dF_j(5) = Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj
      dF_j(6) = Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj
      dF_j(7) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx) + pj)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*qj
#else
      ! Compute flux difference for x-direction
      dF_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)-&
                 X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      dF_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui + pi) -&
                 (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj + pj)
      dF_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                  Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)) -&
                 (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                  Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))
      dF_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                  Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)) -&
                 (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                  Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))
      dF_ij(5) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi) -&
                 (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj)
      dF_ij(6) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj) -&
                 (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj)
      dF_ij(7) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx) + pi)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*qi) -&
                 ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx) + pj)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor by Roe
      !
      ! There are seven eigenvalues
      !
      !   u-cf, u-ca, u-cs, u, u+cs, u+ca, u+cf,
      !
      ! where u us the x-velocity component and ca, cs and cf are the
      ! velocities of th Alfveen waves, the slow and fast waves.
      !-------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))
      anorm = abs(a(1)) ! = sqrt(a(1)*a(1))

      if (anorm .gt. SYS_EPSREAL) then

        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+MYNEWLINE
                TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))/MYNEWLINE
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+MYNEWLINE
                TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))/MYNEWLINE
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx), aux)
        
        ! Compute the square of the Roe-averaged speed of the Alfven waves.
        ! Note that left and right states are interchanged!
        bxPow2_ij = (ROE_MEAN_VALUE(MYNEWLINE
                      X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),MYNEWLINE
                      X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux))**2/MYNEWLINE
                    ROE_MEAN_VALUE(MYNEWLINE
                     DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),MYNEWLINE
                     DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),aux)
        ca_ij = sqrt(bxPow2_ij)
    
        ! Compute the density-averaged magnetic field
        X_ij = ((X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-MYNEWLINE
                 X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2+MYNEWLINE
                (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-MYNEWLINE
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2+MYNEWLINE
                (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-MYNEWLINE
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2)/MYNEWLINE
                (2.0*(sqrt(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))+MYNEWLINE
                      sqrt(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))))

        ! Compute the square of the Roe-averaged magnetic field.
        ! Note that left and right states are interchanged!
        bPow2_ij = (ROE_MEAN_VALUE(MYNEWLINE
                     X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),MYNEWLINE
                     X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2+MYNEWLINE
                    ROE_MEAN_VALUE(MYNEWLINE
                     Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),MYNEWLINE
                     Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2+MYNEWLINE
                    ROE_MEAN_VALUE(MYNEWLINE
                     Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),MYNEWLINE
                     Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2)/MYNEWLINE
                    ROE_MEAN_VALUE(MYNEWLINE
                     DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),MYNEWLINE
                     DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),aux)

        ! Compute the magnitude of the Roe-averaged velocity
        q_ij = ROE_MEAN_VALUE(MYNEWLINE
                X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
                X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2+MYNEWLINE
               ROE_MEAN_VALUE(MYNEWLINE
                Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
                Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2+MYNEWLINE
               ROE_MEAN_VALUE(MYNEWLINE
                Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
                Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2

        ! Compute the Roe-averaged speed of sound
#ifdef THERMALLY_IDEAL_GAS
        aPow2_ij = (2.0-GAMMA)*X_ij+(GAMMA-1.0)*(H_ij-0.5*q_ij-bPow2_ij)
#else
#error "Speed of sound must be implemented!"
#endif

        ! Compute auxiliary quantities
        astPow2_ij = aPow2_ij+bPow2_ij
        aux = sqrt(astPow2_ij**2-4.0*aPow2_ij*bxPow2_ij)

        ! Compute the Roe-averagred speed of the slow and fast waves
        cf_ij = sqrt(0.5*(astPow2_ij+aux))
        cs_ij = sqrt(0.5*(astPow2_ij-aux))
        
        ! Compute eigenvalues
        l1 = abs(u_ij-cf_ij)
        l2 = abs(u_ij-ca_ij)
        l3 = abs(u_ij-cs_ij)
        l4 = abs(u_ij)
        l5 = abs(u_ij+cs_ij)
        l6 = abs(u_ij+ca_ij)
        l7 = abs(u_ij+cf_ij)
        
        ! Compute solution difference U_j-U_i
        Diff = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)

        ! Compute characteristic variables multiplied by the corresponding eigenvalue

        ! TODO

        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"

        ! TODO

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF_i + Diff)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else        
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF_ij + Diff)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_ij + Diff)
#endif
      else

#ifdef MHD_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF_i )
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_ij)
#endif
      end if
    end do

  end subroutine mhd_calcFluxRoeDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRusDiss1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)


!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 1D using the Rusanov dissipation.
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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR1D) :: dF_i,dF_j
#else
    real(DP), dimension(NVAR1D) :: dF_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: cai,caj,cfi,cfj,d_ij
    real(DP) :: aPow2i,aPow2j,astPow2i,astPow2j
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !
      !     /   rho*u            \
      !     | rho*u*u + pT       |
      !     | rho*v*u - Bx*By    |
      ! F = | rho*w*u - Bx*Bz    |
      !     |    By*u - Bx*v     |
      !     |    Bz*u - Bx*w     |
      !     \ (rho*E+pT)*u -Bx*q /
      !
      ! pT = p + q, q = 1/(2*nu)*|B|^2
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      dF_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      dF_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui + pi
      dF_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      dF_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      dF_i(5) = Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi
      dF_i(6) = Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*wi
      dF_i(7) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx) + pi)*ui -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*qi
      
      dF_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      dF_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj + pj
      dF_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)
      dF_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)
      dF_j(5) = Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj
      dF_j(6) = Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj
      dF_j(7) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx) + pj)*uj -&
                X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*qj
#else
      ! Compute flux difference for x-direction
      dF_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)-&
                 X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      dF_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui + pi) -&
                 (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj + pj)
      dF_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                  Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)) -&
                 (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                  Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))
      dF_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
                  Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)) -&
                 (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
                  Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))
      dF_ij(5) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi) -&
                 (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj)
      dF_ij(6) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj) -&
                 (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj)
      dF_ij(7) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx) + pi)*ui -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*qi) -&
                 ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx) + pj)*uj -&
                  X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !
      ! There are seven eigenvalues
      !
      !   u-cf, u-ca, u-cs, u, u+cs, u+ca, u+cf,
      !
      ! where u us the x-velocity component and ca, cs and cf are the
      ! velocities of th Alfveen waves, the slow and fast waves. Since
      !
      !   cf >= ca >= cs >= 0
      !
      ! it suffices to consider only the two eigenvalues
      !
      !   u-cf and u+cf
      !
      ! to construct the Rusanov fluxes
      ! -------------------------------------------------------------------------
      
      ! Compute the speed of the Alfven waves
      cai = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))
      caj = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      aPow2i = GAMMA*PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      aPow2j = GAMMA*PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
#else
#error "Speed of sound must be implemented!"
#endif

      ! Compute auxiliary quantities
      astPow2i = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)/MYNEWLINE
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx) + aPow2i
      astPow2j = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)/MYNEWLINE
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx) + aPow2j

      ! Compute the speed of the fast waves
      cfi = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*cai**2)))
      cfj = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*caj**2)))
            
      ! Scalar dissipation for the Rusanov flux
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj)+&
                  abs(DmatrixCoeffsAtEdge(1,1,idx))*cfj,&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui)+&
                  abs(DmatrixCoeffsAtEdge(1,2,idx))*cfi )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_ij + Diff)
#endif
    end do

  end subroutine mhd_calcFluxRusDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiagMatD1d_sim(DdataAtNode, DmatrixCoeffsAtNode,&
      IverticesAtNode, dscale, DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 1D
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

    ! local variable
    real(DP) :: ui
    integer :: inode

    do inode = 1, size(DcoefficientsAtNode,3)
      
      ! Compute velocity
      ui = X_VELOCITY_1T_FROM_CONSVAR(DdataAtNode,NVAR1D,inode)
      
#ifdef MHD_USE_IBP      
      ! Compute Galerkin coefficient $K_ii = diag(A_i)*C_{ii}$
      DcoefficientsAtNode(1,1,inode) = 0.0
      DcoefficientsAtNode(2,1,inode) = dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(3,1,inode) = dscale * GAMMA*ui*DmatrixCoeffsAtNode(1,inode)
#else
      ! Compute Galerkin coefficient $K_ii = -diag(A_i)*C_{ii}$
      DcoefficientsAtNode(1,1,inode) = 0.0
      DcoefficientsAtNode(2,1,inode) = -dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(3,1,inode) = -dscale * GAMMA*ui*DmatrixCoeffsAtNode(1,inode)
#endif
    end do

  end subroutine mhd_calcMatDiagMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiag1d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 1D
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

    ! local variable
    real(DP) :: ui,Ei,uPow2i
    integer :: inode

    do inode = 1, size(DcoefficientsAtNode,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_1T_FROM_CONSVAR(DdataAtNode,NVAR1D,inode)
      Ei = TOTAL_ENERGY_1T_FROM_CONSVAR(DdataAtNode,NVAR1D,inode)
      uPow2i = ui*ui

#ifdef MHD_USE_IBP      
      ! Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      DcoefficientsAtNode(1,1,inode) = 0.0
      DcoefficientsAtNode(2,1,inode) = dscale * (GAMMA-3.0)/2.0*uPow2i*DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(3,1,inode) = dscale * ((GAMMA-1.0)*uPow2i-GAMMA*Ei)*&
                                                ui*DmatrixCoeffsAtNode(1,inode)

      DcoefficientsAtNode(4,1,inode) = dscale * DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(5,1,inode) = dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(6,1,inode) = dscale * (GAMMA*Ei-3.0*(GAMMA-1.0)/2.0*&
                                                uPow2i)*DmatrixCoeffsAtNode(1,inode)

      DcoefficientsAtNode(7,1,inode) = 0.0
      DcoefficientsAtNode(8,1,inode) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(9,1,inode) = dscale * GAMMA*ui*DmatrixCoeffsAtNode(1,inode)
#else
      ! Compute Galerkin coefficient $K_ii = -A_i*C_{ii}$
      DcoefficientsAtNode(1,1,inode) = 0.0
      DcoefficientsAtNode(2,1,inode) = -dscale * (GAMMA-3.0)/2.0*uPow2i*DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(3,1,inode) = -dscale * ((GAMMA-1.0)*uPow2i-GAMMA*Ei)*&
                                                 ui*DmatrixCoeffsAtNode(1,inode)
      
      DcoefficientsAtNode(4,1,inode) = -dscale * DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(5,1,inode) = -dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(6,1,inode) = -dscale * (GAMMA*Ei-3.0*(GAMMA-1.0)/2.0*&
                                                 uPow2i)*DmatrixCoeffsAtNode(1,inode)
      
      DcoefficientsAtNode(7,1,inode) = 0.0
      DcoefficientsAtNode(8,1,inode) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(9,1,inode) = -dscale * GAMMA*ui*DmatrixCoeffsAtNode(1,inode)
#endif
    end do
    
  end subroutine mhd_calcMatDiag1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGalMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 1D
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

    ! local variable
    real(DP) :: ui,uj
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef MHD_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,2,idx)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,1,idx)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = -dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = -dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = -dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = -dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)
#endif
    end do

  end subroutine mhd_calcMatGalMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGal1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 1D
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

    ! local variable
    real(DP) :: ui,uj,Ei,Ej,uPow2i,uPow2j
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      Ei = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uPow2i = ui*ui
      
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      Ej = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      uPow2j = uj*uj
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef MHD_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = dscale * (GAMMA-3.0)/2.0*uPow2j*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * ((GAMMA-1.0)*uPow2j-GAMMA*Ej)*&
                                              uj*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(4,2,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(5,2,idx) = dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(6,2,idx) = dscale * (GAMMA*Ej-3.0*(GAMMA-1.0)/2.0*&
                                              uPow2j)*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(7,2,idx) = 0.0
      DcoefficientsAtEdge(8,2,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(9,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,2,idx)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = dscale * (GAMMA-3.0)/2.0*uPow2i*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * ((GAMMA-1.0)*uPow2i-GAMMA*Ei)*&
                                              ui*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(4,3,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(5,3,idx) = dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(6,3,idx) = dscale * (GAMMA*Ei-3.0*(GAMMA-1.0)/2.0*&
                                              uPow2i)*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(7,3,idx) = 0.0
      DcoefficientsAtEdge(8,3,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(9,3,idx) = dscale * GAMMA*uj*ui*DmatrixCoeffsAtEdge(1,1,idx)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = -dscale * (GAMMA-3.0)/2.0*uPow2j*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = -dscale * ((GAMMA-1.0)*uPow2j-GAMMA*Ej)*&
                                               uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(4,2,idx) = -dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(5,2,idx) = -dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(6,2,idx) = -dscale * (GAMMA*Ej-3.0*(GAMMA-1.0)/2.0*&
                                               uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(7,2,idx) = 0.0
      DcoefficientsAtEdge(8,2,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(9,2,idx) = -dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = -dscale * (GAMMA-3.0)/2.0*uPow2i*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = -dscale * ((GAMMA-1.0)*uPow2i-GAMMA*Ei)*&
                                               ui*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(4,3,idx) = -dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(5,3,idx) = -dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(6,3,idx) = -dscale * (GAMMA*Ei-3.0*(GAMMA-1.0)/2.0*&
                                               uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(7,3,idx) = 0.0
      DcoefficientsAtEdge(8,3,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(9,3,idx) = -dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)
#endif
    end do

  end subroutine mhd_calcMatGal1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDissMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies scalar artificial viscosities in 1D
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

    ! local variable
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: anorm,aux,H_ij,q_ij,ui,uj,u_ij
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef MHD_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,2,idx)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,1,idx)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = -dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = -dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = -dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = -dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)
#endif
      
      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ji}-C_{ij})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))
      anorm = abs(a(1))
      
      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx), aux)
        
        ! Compute auxiliary values
        q_ij = 0.5*(u_ij*u_ij)
        
        ! Compute scalar dissipation
#ifdef THERMALLY_IDEAL_GAS
        DcoefficientsAtEdge(:,1,idx) = dscale * (abs(a(1)*u_ij) +&
            anorm*sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)))
#else
#error "Speed of sound must be implemented!"
#endif
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0
        
      end if
    end do

  end subroutine mhd_calcMatScDissMatD1d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDiss1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies scalar artificial viscosities in 1D
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

    ! local variable
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: anorm,aux,Ei,Ej,H_ij,q_ij,ui,uj,u_ij,uPow2i,uPow2j
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
    
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      Ei = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uPow2i = ui*ui
      
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      Ej = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      uPow2j = uj*uj

#ifdef MHD_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = dscale * (GAMMA-3.0)/2.0*uPow2j*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * ((GAMMA-1.0)*uPow2j-GAMMA*Ej)*&
                                              uj*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(4,2,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(5,2,idx) = dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(6,2,idx) = dscale * (GAMMA*Ej-3.0*(GAMMA-1.0)/2.0*&
                                              uPow2j)*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(7,2,idx) = 0.0
      DcoefficientsAtEdge(8,2,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(9,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,2,idx)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = dscale * (GAMMA-3.0)/2.0*uPow2i*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * ((GAMMA-1.0)*uPow2i-GAMMA*Ei)*&
                                              ui*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(4,3,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(5,3,idx) = dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(6,3,idx) = dscale * (GAMMA*Ei-3.0*(GAMMA-1.0)/2.0*&
                                              uPow2i)*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(7,3,idx) = 0.0
      DcoefficientsAtEdge(8,3,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(9,3,idx) = dscale * GAMMA*uj*ui*DmatrixCoeffsAtEdge(1,1,idx)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = -dscale * (GAMMA-3.0)/2.0*uPow2j*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = -dscale * ((GAMMA-1.0)*uPow2j-GAMMA*Ej)*&
                                               uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(4,2,idx) = -dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(5,2,idx) = -dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(6,2,idx) = -dscale * (GAMMA*Ej-3.0*(GAMMA-1.0)/2.0*&
                                               uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(7,2,idx) = 0.0
      DcoefficientsAtEdge(8,2,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(9,2,idx) = -dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = -dscale * (GAMMA-3.0)/2.0*uPow2i*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = -dscale * ((GAMMA-1.0)*uPow2i-GAMMA*Ei)*&
                                               ui*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(4,3,idx) = -dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(5,3,idx) = -dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(6,3,idx) = -dscale * (GAMMA*Ei-3.0*(GAMMA-1.0)/2.0*&
                                               uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(7,3,idx) = 0.0
      DcoefficientsAtEdge(8,3,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(9,3,idx) = -dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))
      anorm = abs(a(1))
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx), aux)
        
        ! Compute auxiliary values
        q_ij = 0.5*(u_ij*u_ij)
        
        ! Compute scalar dissipation
#ifdef THERMALLY_IDEAL_GAS
        aux = dscale * (abs(a(1)*u_ij) +&
            anorm*sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)))
#else
#error "Speed of sound must be implemented!"
#endif
        
        DcoefficientsAtEdge(1,1,idx) = aux
        DcoefficientsAtEdge(5,1,idx) = aux
        DcoefficientsAtEdge(9,1,idx) = aux
      end if
    end do
    
  end subroutine mhd_calcMatScDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDissMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 1D
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

    ! local variable
    real(DP), dimension(NVAR1D,NVAR1D) :: R_ij,L_ij
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: aux,H_ij,q_ij,ui,uj,u_ij
    real(DP) :: l1,l2,l3,anorm,c_ij,cPow2_ij,b1,b2
    integer :: idx
    
    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef MHD_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,2,idx)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,1,idx)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = -dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = -dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = -dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = -dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(1,1,idx)-&
                  DmatrixCoeffsAtEdge(1,2,idx))
      anorm = abs(a(1))
      
      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx), aux)
                
        ! Compute auxiliary values
        q_ij  = 0.5*(u_ij*u_ij)
#ifdef THERMALLY_IDEAL_GAS
        cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(cPow2_ij)
        
        b2   = (GAMMA-1.0)/cPow2_ij
        b1    = b2*q_ij

        ! Diagonal matrix of eigenvalues
        l1 = abs(u_ij-c_ij)
        l2 = abs(u_ij)
        l3 = abs(u_ij+c_ij)
        
        ! Matrix of right eigenvectors
        R_ij(1,1) =  l1
        R_ij(2,1) =  l1*(u_ij-c_ij)
        R_ij(3,1) =  l1*(H_ij-c_ij*u_ij)
        
        R_ij(1,2) =  l2
        R_ij(2,2) =  l2*u_ij
        R_ij(3,2) =  l2*q_ij
        
        R_ij(1,3) =  l3
        R_ij(2,3) =  l3*(u_ij+c_ij)
        R_ij(3,3) =  l3*(H_ij+c_ij*u_ij)
        
        ! Matrix of left eigenvectors
        L_ij(1,1) = 0.5*(b1+u_ij/c_ij)
        L_ij(2,1) = 1.0-b1
        L_ij(3,1) = 0.5*(b1-u_ij/c_ij)
        
        L_ij(1,2) = -0.5*(b2*u_ij+1/c_ij)
        L_ij(2,2) =  b2*u_ij
        L_ij(3,2) = -0.5*(b2*u_ij-1/c_ij)
        
        L_ij(1,3) = 0.5*b2
        L_ij(2,3) = -b2
        L_ij(3,3) = 0.5*b2

        ! Include scaling parameter
        anorm = dscale*anorm

        ! Compute tensorial dissipation D_ij = diag(R_ij*|Lbd_ij|*L_ij)*I
        DcoefficientsAtEdge(1,1,idx) = anorm*( R_ij(1,1)*L_ij(1,1)+&
            R_ij(1,2)*L_ij(2,1)+R_ij(1,3)*L_ij(3,1)  )
        DcoefficientsAtEdge(2,1,idx) = anorm*( R_ij(2,1)*L_ij(1,2)+&
            R_ij(2,2)*L_ij(2,2)+R_ij(2,3)*L_ij(3,2)  )
        DcoefficientsAtEdge(3,1,idx) = anorm*( R_ij(3,1)*L_ij(1,3)+&
            R_ij(3,2)*L_ij(2,3)+R_ij(3,3)*L_ij(3,3)  )
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0

      end if
    end do

  end subroutine mhd_calcMatRoeDissMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDiss1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 1D
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

    ! local variable
    real(DP), dimension(NVAR1D,NVAR1D) :: R_ij,L_ij
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: aux,Ei,Ej,H_ij,q_ij,ui,uj,u_ij
    real(DP) :: l1,l2,l3,anorm,c_ij,cPow2_ij,b1,b2,uPow2i,uPow2j
    integer :: idx,i,j,k

    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      Ei = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uPow2i = ui*ui
      
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      Ej = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      uPow2j = uj*uj

#ifdef MHD_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = dscale * (GAMMA-3.0)/2.0*uPow2j*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * ((GAMMA-1.0)*uPow2j-GAMMA*Ej)*&
                                              uj*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(4,2,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(5,2,idx) = dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(6,2,idx) = dscale * (GAMMA*Ej-3.0*(GAMMA-1.0)/2.0*&
                                              uPow2j)*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(7,2,idx) = 0.0
      DcoefficientsAtEdge(8,2,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(9,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,2,idx)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = dscale * (GAMMA-3.0)/2.0*uPow2i*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * ((GAMMA-1.0)*uPow2i-GAMMA*Ei)*&
                                              ui*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(4,3,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(5,3,idx) = dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(6,3,idx) = dscale * (GAMMA*Ei-3.0*(GAMMA-1.0)/2.0*&
                                              uPow2i)*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(7,3,idx) = 0.0
      DcoefficientsAtEdge(8,3,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(9,3,idx) = dscale * GAMMA*uj*ui*DmatrixCoeffsAtEdge(1,1,idx)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = -dscale * (GAMMA-3.0)/2.0*uPow2j*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = -dscale * ((GAMMA-1.0)*uPow2j-GAMMA*Ej)*&
                                               uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(4,2,idx) = -dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(5,2,idx) = -dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(6,2,idx) = -dscale * (GAMMA*Ej-3.0*(GAMMA-1.0)/2.0*&
                                               uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(7,2,idx) = 0.0
      DcoefficientsAtEdge(8,2,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(9,2,idx) = -dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = -dscale * (GAMMA-3.0)/2.0*uPow2i*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = -dscale * ((GAMMA-1.0)*uPow2i-GAMMA*Ei)*&
                                               ui*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(4,3,idx) = -dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(5,3,idx) = -dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(6,3,idx) = -dscale * (GAMMA*Ei-3.0*(GAMMA-1.0)/2.0*&
                                               uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(7,3,idx) = 0.0
      DcoefficientsAtEdge(8,3,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(9,3,idx) = -dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))
      anorm = abs(a(1))

      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx), aux)
        
        ! Compute auxiliary values
        q_ij  = 0.5*u_ij*u_ij
#ifdef THERMALLY_IDEAL_GAS
        cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(cPow2_ij)

        b2   = (GAMMA-1.0)/cPow2_ij
        b1    = b2*q_ij
        
        ! Diagonal matrix of eigenvalues
        l1 = abs(u_ij-c_ij)
        l2 = abs(u_ij)
        l3 = abs(u_ij+c_ij)
        
        ! Matrix of right eigenvectors
        R_ij(1,1) =  l1
        R_ij(2,1) =  l1*(u_ij-c_ij)
        R_ij(3,1) =  l1*(H_ij-c_ij*u_ij)
        
        R_ij(1,2) =  l2
        R_ij(2,2) =  l2*u_ij
        R_ij(3,2) =  l2*q_ij
        
        R_ij(1,3) =  l3
        R_ij(2,3) =  l3*(u_ij+c_ij)
        R_ij(3,3) =  l3*(H_ij+c_ij*u_ij)
        
        ! Matrix of left eigenvectors
        L_ij(1,1) = 0.5*(b1+u_ij/c_ij)
        L_ij(2,1) = 1.0-b1
        L_ij(3,1) = 0.5*(b1-u_ij/c_ij)
        
        L_ij(1,2) = -0.5*(b2*u_ij+1/c_ij)
        L_ij(2,2) =  b2*u_ij
        L_ij(3,2) = -0.5*(b2*u_ij-1/c_ij)
        
        L_ij(1,3) = 0.5*b2
        L_ij(2,3) = -b2
        L_ij(3,3) = 0.5*b2
        
        ! Include scaling parameter
        anorm = dscale*anorm

        ! Compute tensorial dissipation D_ij = R_ij*|Lbd_ij|*L_ij
        do i = 1, NVAR1D
          do j = 1, NVAR1D
            aux = 0.0
            do k = 1, NVAR1D
              aux = aux + R_ij(i,k)*L_ij(k,j)
            end do
            DcoefficientsAtEdge(NVAR1D*(j-1)+i,1,idx) = anorm*aux
          end do
        end do
        
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0
        
      end if
    end do

  end subroutine mhd_calcMatRoeDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDissMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 1D
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

    ! local variable
    real(DP) :: ui,uj,ci,cj,Ei,Ej
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      Ei = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      Ej = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

#ifdef MHD_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,2,idx)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,1,idx)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = -dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = -dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = -dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = -dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*ui*ui), SYS_EPSREAL))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*uj*uj), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Compute dissipation tensor D_ij
      DcoefficientsAtEdge(:,1,idx) = dscale *&
          max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj)+&
               abs(DmatrixCoeffsAtEdge(1,1,idx))*cj,&
               abs(DmatrixCoeffsAtEdge(1,2,idx)*ui)+&
               abs(DmatrixCoeffsAtEdge(1,2,idx))*ci )
    end do
    
  end subroutine mhd_calcMatRusDissMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDiss1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 1D
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

    ! local variable
    real(DP) :: ui,uj,ci,cj,Ei,Ej,uPow2i,uPow2j,aux
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      Ei = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uPow2i = ui*ui

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      Ej = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      uPow2j = uj*uj

#ifdef MHD_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = dscale * (GAMMA-3.0)/2.0*uPow2j*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * ((GAMMA-1.0)*uPow2j-GAMMA*Ej)*&
                                              uj*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(4,2,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(5,2,idx) = dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(6,2,idx) = dscale * (GAMMA*Ej-3.0*(GAMMA-1.0)/2.0*&
                                              uPow2j)*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(7,2,idx) = 0.0
      DcoefficientsAtEdge(8,2,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(9,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,2,idx)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = dscale * (GAMMA-3.0)/2.0*uPow2i*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * ((GAMMA-1.0)*uPow2i-GAMMA*Ei)*&
                                              ui*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(4,3,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(5,3,idx) = dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(6,3,idx) = dscale * (GAMMA*Ei-3.0*(GAMMA-1.0)/2.0*&
                                              uPow2i)*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(7,3,idx) = 0.0
      DcoefficientsAtEdge(8,3,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(9,3,idx) = dscale * GAMMA*uj*ui*DmatrixCoeffsAtEdge(1,1,idx)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = -dscale * (GAMMA-3.0)/2.0*uPow2j*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = -dscale * ((GAMMA-1.0)*uPow2j-GAMMA*Ej)*&
                                               uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(4,2,idx) = -dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(5,2,idx) = -dscale * (3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(6,2,idx) = -dscale * (GAMMA*Ej-3.0*(GAMMA-1.0)/2.0*&
                                               uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(7,2,idx) = 0.0
      DcoefficientsAtEdge(8,2,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(9,2,idx) = -dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = -dscale * (GAMMA-3.0)/2.0*uPow2i*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = -dscale * ((GAMMA-1.0)*uPow2i-GAMMA*Ei)*&
                                               ui*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(4,3,idx) = -dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(5,3,idx) = -dscale * (3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(6,3,idx) = -dscale * (GAMMA*Ei-3.0*(GAMMA-1.0)/2.0*&
                                               uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(7,3,idx) = 0.0
      DcoefficientsAtEdge(8,3,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(9,3,idx) = -dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*ui*ui), SYS_EPSREAL))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*uj*uj), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif

      ! Compute dissipation tensor D_ij
      aux = dscale * max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj)+&
                          abs(DmatrixCoeffsAtEdge(1,1,idx))*cj,&
                          abs(DmatrixCoeffsAtEdge(1,2,idx)*ui)+&
                          abs(DmatrixCoeffsAtEdge(1,2,idx))*ci )

      DcoefficientsAtEdge(:,1,idx) = 0.0
      DcoefficientsAtEdge(1,1,idx) = aux
      DcoefficientsAtEdge(5,1,idx) = aux
      DcoefficientsAtEdge(9,1,idx) = aux
    end do
    
  end subroutine mhd_calcMatRusDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcCharacteristics1d_sim(Dweight, DdataAtEdge,&
      DcharVariablesAtEdge, DeigenvaluesAtEdge,&
      DrightEigenvectorsAtEdge, DleftEigenvectorsAtEdge, rcollection)

!<description>
    ! This subroutine computes the characteristic variables in 1D
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

    ! local variables
    real(DP), dimension(NVAR1D) :: Diff
    real(DP) :: u_ij,H_ij,q_ij,c_ij,aux,ui,uj,cPow2_ij,anorm,b1_ij,b2_ij
    integer :: idx

    ! Compute norm of weighting coefficient
    anorm = abs(Dweight(1))

    ! Check if weighting coefficient is zero
    if (anorm .le. SYS_EPSREAL) then
      if (present(DcharVariablesAtEdge))     DcharVariablesAtEdge     = 0.0
      if (present(DeigenvaluesAtEdge))       DeigenvaluesAtEdge       = 0.0
      if (present(DrightEigenvectorsAtEdge)) DrightEigenvectorsAtEdge = 0.0
      if (present(DleftEigenvectorsAtEdge))  DleftEigenvectorsAtEdge  = 0.0

      ! That's it
      return
    end if
    
    
    ! Do we have to compute characteristic variables
    if (present(DcharVariablesAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)
        
        ! Compute velocities
        ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
        uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx), aux)
               
        ! Compute auxiliary variables
        q_ij  = 0.5*(u_ij*u_ij)
#ifdef THERMALLY_IDEAL_GAS
        cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(cPow2_ij)

        b2_ij = (GAMMA-1.0)/cPow2_ij
        b1_ij  = b2_ij*q_ij
        
        ! Compute solution difference U_j-U_i
        Diff = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
          
        ! Compute characteristic variables
        DcharVariablesAtEdge(1,idx) = anorm * 0.5 *&
            (       (b1_ij+u_ij/c_ij)*Diff(1)-&
                (b2_ij*u_ij+1.0/c_ij)*Diff(2)+&
                                b2_ij*Diff(3))
        DcharVariablesAtEdge(2,idx) = anorm *&
            (             (1.0-b1_ij)*Diff(1)+&
                           b2_ij*u_ij*Diff(2)-&
                                b2_ij*Diff(3) )
        DcharVariablesAtEdge(3,idx) = anorm * 0.5 *&
            (       (b1_ij-u_ij/c_ij)*Diff(1)-&
                (b2_ij*u_ij-1.0/c_ij)*Diff(2)+&
                                b2_ij*Diff(3) )
      end do
    end if


    ! Do we have to compute eigenvalues
    if (present(DeigenvaluesAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)
        
        ! Compute velocities
        ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
        uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx), aux)
        
        ! Compute auxiliary variable
#ifdef THERMALLY_IDEAL_GAS
        c_ij = sqrt(max((GAMMA-1.0)*(H_ij-0.5*(u_ij*u_ij)), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif

        ! Compute eigenvalues
        DeigenvaluesAtEdge(1,idx) = u_ij-c_ij
        DeigenvaluesAtEdge(2,idx) = u_ij
        DeigenvaluesAtEdge(3,idx) = u_ij+c_ij
      end do
    end if


    ! Do we have to compute right eigenvectors
    if (present(DrightEigenvectorsAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)
        
        ! Compute velocities
        ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
        uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx), aux)

        ! Compute auxiliary variables
        q_ij  = 0.5*(u_ij*u_ij)
#ifdef THERMALLY_IDEAL_GAS
        c_ij  = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif

        ! Compute right eigenvectors
        DrightEigenvectorsAtEdge(1,idx) =  1.0
        DrightEigenvectorsAtEdge(2,idx) =  u_ij-c_ij
        DrightEigenvectorsAtEdge(3,idx) =  H_ij-u_ij*c_ij

        DrightEigenvectorsAtEdge(4,idx) =  1.0
        DrightEigenvectorsAtEdge(5,idx) =  u_ij
        DrightEigenvectorsAtEdge(6,idx) =  q_ij

        DrightEigenvectorsAtEdge(7,idx) =  1.0
        DrightEigenvectorsAtEdge(8,idx) =  u_ij+c_ij
        DrightEigenvectorsAtEdge(9,idx) =  H_ij+u_ij*c_ij
      end do
    end if


    ! Do we have to compute left eigenvectors
    if (present(DleftEigenvectorsAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)
        
        ! Compute velocities
        ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
        uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx), aux)

        ! Compute auxiliary variables
        q_ij  = 0.5*(u_ij*u_ij)
#ifdef THERMALLY_IDEAL_GAS
        cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(cPow2_ij)

        b2_ij = (GAMMA-1.0)/cPow2_ij
        b1_ij = b2_ij*q_ij

        ! Compute left eigenvectors
        DleftEigenvectorsAtEdge(1,idx) = 0.5 * (b1_ij+u_ij/c_ij)
        DleftEigenvectorsAtEdge(2,idx) =       (1.0-b1_ij)
        DleftEigenvectorsAtEdge(3,idx) = 0.5 * (b1_ij-u_ij/c_ij)

        DleftEigenvectorsAtEdge(4,idx) =-0.5 * (b2_ij*u_ij+1.0/c_ij)
        DleftEigenvectorsAtEdge(5,idx) =       (b2_ij*u_ij)
        DleftEigenvectorsAtEdge(6,idx) =-0.5 * (b2_ij*u_ij-1.0/c_ij)

        DleftEigenvectorsAtEdge(7,idx) = 0.5*b2_ij
        DleftEigenvectorsAtEdge(8,idx) =    -b2_ij
        DleftEigenvectorsAtEdge(9,idx) = 0.5*b2_ij
      end do
    end if
    
  end subroutine mhd_calcCharacteristics1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTScalarDiss1d(&
      U1_i, U1_j, U2_i, U2_j, Coeff_ij, Coeff_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 1D using scalar dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: Coeff_ij,Coeff_ji

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

    ! local variables
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: ui,uj,vi,vj,wi,wj,cf_ij
    real(DP) :: aux,u_ij,q_ij,H_ij,X_ij
    real(DP) :: aPow2_ij,bPow2_ij,bxPow2_ij,astPow2_ij
    real(DP) :: anorm,d_ij

    ! Compute velocities
    ui = X_VELOCITY_FROM_CONSVAR(U2_i,NVAR1D)
    vi = Y_VELOCITY_FROM_CONSVAR(U2_j,NVAR1D)
    wi = Z_VELOCITY_FROM_CONSVAR(U2_i,NVAR1D)
    uj = X_VELOCITY_FROM_CONSVAR(U2_j,NVAR1D)
    vj = Y_VELOCITY_FROM_CONSVAR(U2_i,NVAR1D)
    wj = Z_VELOCITY_FROM_CONSVAR(U2_j,NVAR1D)

    ! Compute skew-symmetric coefficient
    a = 0.5*(Coeff_ij-Coeff_ji)

    ! Compute Roe mean values
    aux  = ROE_MEAN_RATIO(MYNEWLINE
            DENSITY_FROM_CONSVAR(U2_i,NVAR1D),MYNEWLINE
            DENSITY_FROM_CONSVAR(U2_j,NVAR1D))
    u_ij = ROE_MEAN_VALUE(ui,uj,aux)
    H_ij = ROE_MEAN_VALUE(MYNEWLINE
           (TOTAL_ENERGY_FROM_CONSVAR(U2_i,NVAR1D)+MYNEWLINE
            TOTAL_PRESSURE_FROM_CONSVAR_1D(U2_i,NVAR1D))/MYNEWLINE
            DENSITY_FROM_CONSVAR(U2_i,NVAR1D),MYNEWLINE
           (TOTAL_ENERGY_FROM_CONSVAR(U2_j,NVAR1D)+MYNEWLINE
            TOTAL_PRESSURE_FROM_CONSVAR_1D(U2_j,NVAR1D))/MYNEWLINE
            DENSITY_FROM_CONSVAR(U2_j,NVAR1D), aux)

    ! Compute the square of the Roe-averaged speed of the Alfven waves.
    ! Note that left and right states are interchanged!
    bxPow2_ij = (ROE_MEAN_VALUE(MYNEWLINE
                  X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D),MYNEWLINE
                  X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D),aux))**2/MYNEWLINE
                ROE_MEAN_VALUE(MYNEWLINE
                 DENSITY_FROM_CONSVAR(U2_j,NVAR1D),MYNEWLINE
                 DENSITY_FROM_CONSVAR(U2_i,NVAR1D),aux)

    ! Compute the density-averaged magnetic field
    X_ij = ((X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D)-MYNEWLINE
             X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D))**2 +MYNEWLINE
            (Y_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D)-MYNEWLINE
             Y_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D))**2 +MYNEWLINE
            (Z_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D)-MYNEWLINE
             Z_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D))**2)/MYNEWLINE
            (2.0*(sqrt(X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D))+MYNEWLINE
                  sqrt(X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D))))

    ! Compute the square of the Roe-averaged magnetic field.
    ! Note that left and right states are interchanged!
    bPow2_ij = (ROE_MEAN_VALUE(MYNEWLINE
                 X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D),MYNEWLINE
                 X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D),aux)**2+MYNEWLINE
                ROE_MEAN_VALUE(MYNEWLINE
                 Y_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D),MYNEWLINE
                 Y_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D),aux)**2+MYNEWLINE
                ROE_MEAN_VALUE(MYNEWLINE
                 Z_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D),MYNEWLINE
                 Z_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D),aux)**2)/MYNEWLINE
                ROE_MEAN_VALUE(MYNEWLINE
                 DENSITY_FROM_CONSVAR(U2_j,NVAR1D),MYNEWLINE
                 DENSITY_FROM_CONSVAR(U2_i,NVAR1D),aux)

    ! Compute the magnitude of the Roe-averaged velocity
    q_ij = ROE_MEAN_VALUE(MYNEWLINE
            X_VELOCITY_FROM_CONSVAR(U2_i,NVAR1D),MYNEWLINE
            X_VELOCITY_FROM_CONSVAR(U2_j,NVAR1D),aux)**2+MYNEWLINE
           ROE_MEAN_VALUE(MYNEWLINE
            Y_VELOCITY_FROM_CONSVAR(U2_i,NVAR1D),MYNEWLINE
            Y_VELOCITY_FROM_CONSVAR(U2_j,NVAR1D),aux)**2+MYNEWLINE
           ROE_MEAN_VALUE(MYNEWLINE
            Z_VELOCITY_FROM_CONSVAR(U2_i,NVAR1D),MYNEWLINE
            Z_VELOCITY_FROM_CONSVAR(U2_j,NVAR1D),aux)**2

    ! Compute the Roe-averaged speed of sound
#ifdef THERMALLY_IDEAL_GAS
    aPow2_ij = (2.0-GAMMA)*X_ij + (GAMMA-1.0)*(H_ij-0.5*q_ij-bPow2_ij)
#else
#error "Speed of sound must be implemented!"
#endif

    ! Compute auxiliary variables
    astPow2_ij = aPow2_ij+bPow2_ij
    aux = sqrt(astPow2_ij**2-4.0*aPow2_ij*bxPow2_ij)
            
    ! Compute the Roe-averagred speed of the fast waves
    cf_ij = sqrt(0.5*(astPow2_ij+aux))

    ! Scalar dissipation
    d_ij = abs(u_ij*a(1)) + anorm*cf_ij

    ! Compute conservative fluxes
    F_ij = dscale1*(U1_i-U1_j) + dscale2*d_ij*(U2_i-U2_j)

  end subroutine mhd_calcFluxFCTScalarDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTTensorDiss1d(&
      U1_i, U1_j, U2_i, U2_j, Coeff_ij, Coeff_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 1D using tensorial dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: Coeff_ij,Coeff_ji

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

    ! local variables
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: aux,u_ij,q_ij,H_ij,X_ij,ca_ij,cs_ij,cf_ij
    real(DP) :: aPow2_ij,bPow2_ij,bxPow2_ij,astPow2_ij
    real(DP) :: l1,l2,l3,l4,l5,l6,l7,w1,w2,w3,w4,w5,w6,w7
    real(DP) :: anorm

    ! Compute skew-symmetric coefficient
    a = 0.5*(Coeff_ij-Coeff_ji); anorm = abs(a(1))

    ! Compute Roe mean values
    aux  = ROE_MEAN_RATIO(MYNEWLINE
            DENSITY_FROM_CONSVAR(U2_i,NVAR1D),MYNEWLINE
            DENSITY_FROM_CONSVAR(U2_j,NVAR1D))
    u_ij = ROE_MEAN_VALUE(MYNEWLINE
            X_VELOCITY_FROM_CONSVAR(U2_i,NVAR1D),MYNEWLINE
            X_VELOCITY_FROM_CONSVAR(U2_j,NVAR1D),aux)
    H_ij = ROE_MEAN_VALUE(MYNEWLINE
           (TOTAL_ENERGY_FROM_CONSVAR(U2_i,NVAR1D)+MYNEWLINE
            TOTAL_PRESSURE_FROM_CONSVAR_1D(U2_i,NVAR1D))/MYNEWLINE
            DENSITY_FROM_CONSVAR(U2_i,NVAR1D),MYNEWLINE
           (TOTAL_ENERGY_FROM_CONSVAR(U2_j,NVAR1D)+MYNEWLINE
            TOTAL_PRESSURE_FROM_CONSVAR_1D(U2_j,NVAR1D))/MYNEWLINE
            DENSITY_FROM_CONSVAR(U2_j,NVAR1D), aux)

    ! Compute the square of the Roe-averaged speed of the Alfven waves.
    ! Note that left and right states are interchanged!
    bxPow2_ij = (ROE_MEAN_VALUE(MYNEWLINE
                  X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D),MYNEWLINE
                  X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D),aux))**2/MYNEWLINE
                ROE_MEAN_VALUE(MYNEWLINE
                 DENSITY_FROM_CONSVAR(U2_j,NVAR1D),MYNEWLINE
                 DENSITY_FROM_CONSVAR(U2_i,NVAR1D),aux)
    ca_ij = sqrt(bxPow2_ij)
    
    ! Compute the density-averaged magnetic field
    X_ij = ((X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D)-MYNEWLINE
             X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D))**2+MYNEWLINE
            (Y_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D)-MYNEWLINE
             Y_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D))**2+MYNEWLINE
            (Z_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D)-MYNEWLINE
             Z_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D))**2)/MYNEWLINE
            (2.0*(sqrt(X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D))+MYNEWLINE
                  sqrt(X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D))))

    ! Compute the square of the Roe-averaged magnetic field.
    ! Note that left and right states are interchanged!
    bPow2_ij = (ROE_MEAN_VALUE(MYNEWLINE
                 X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D),MYNEWLINE
                 X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D),aux)**2+MYNEWLINE
                ROE_MEAN_VALUE(MYNEWLINE
                 Y_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D),MYNEWLINE
                 Y_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D),aux)**2+MYNEWLINE
                ROE_MEAN_VALUE(MYNEWLINE
                 Z_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D),MYNEWLINE
                 Z_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D),aux)**2)/MYNEWLINE
                ROE_MEAN_VALUE(MYNEWLINE
                 DENSITY_FROM_CONSVAR(U2_j,NVAR1D),MYNEWLINE
                 DENSITY_FROM_CONSVAR(U2_i,NVAR1D),aux)

    ! Compute the magnitude of the Roe-averaged velocity
    q_ij = ROE_MEAN_VALUE(MYNEWLINE
            X_VELOCITY_FROM_CONSVAR(U2_i,NVAR1D),MYNEWLINE
            X_VELOCITY_FROM_CONSVAR(U2_j,NVAR1D),aux)**2+MYNEWLINE
           ROE_MEAN_VALUE(MYNEWLINE
            Y_VELOCITY_FROM_CONSVAR(U2_i,NVAR1D),MYNEWLINE
            Y_VELOCITY_FROM_CONSVAR(U2_j,NVAR1D),aux)**2+MYNEWLINE
           ROE_MEAN_VALUE(MYNEWLINE
            Z_VELOCITY_FROM_CONSVAR(U2_i,NVAR1D),MYNEWLINE
            Z_VELOCITY_FROM_CONSVAR(U2_j,NVAR1D),aux)**2

        ! Compute the Roe-averaged speed of sound
#ifdef THERMALLY_IDEAL_GAS
    aPow2_ij = (2.0-GAMMA)*X_ij+(GAMMA-1.0)*(H_ij-0.5*q_ij-bPow2_ij)
#else
#error "Speed of sound must be implemented!"
#endif

    ! Compute auxiliary quantities
    astPow2_ij = aPow2_ij+bPow2_ij
    aux = sqrt(astPow2_ij**2-4.0*aPow2_ij*bxPow2_ij)

    ! Compute the Roe-averagred speed of the slow and fast waves
    cf_ij = sqrt(0.5*(astPow2_ij+aux))
    cs_ij = sqrt(0.5*(astPow2_ij-aux))
    
    ! Compute eigenvalues
    l1 = abs(u_ij-cf_ij)
    l2 = abs(u_ij-ca_ij)
    l3 = abs(u_ij-cs_ij)
    l4 = abs(u_ij)
    l5 = abs(u_ij+cs_ij)
    l6 = abs(u_ij+ca_ij)
    l7 = abs(u_ij+cf_ij)
    
    ! Compute solution difference U2_i-U2_j
    Diff = U2_i-U2_j

    ! Compute characteristic variables multiplied by the corresponding eigenvalue

    ! TODO

    ! Compute "R_ij * |Lbd_ij| * L_ij * dU"

    ! TODO

    ! Compute conservative fluxes
    F_ij = dscale1*(U1_i-U1_j) + dscale2*anorm*Diff

  end subroutine mhd_calcFluxFCTTensorDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTRusanov1d(&
      U1_i, U1_j, U2_i, U2_j, Coeff_ij, Coeff_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 1D using the Rusanov dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: Coeff_ij,Coeff_ji

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

    ! local variables
    real(DP) :: ui,uj,cai,caj,aPow2i,aPow2j,astPow2i,astPow2j,cfi,cfj,d_ij
    
    ! Compute velocities
    ui = X_VELOCITY_FROM_CONSVAR(U2_i,NVAR1D)
    uj = X_VELOCITY_FROM_CONSVAR(U2_j,NVAR1D)
  
    ! Compute the speed of the Alfven waves
    cai = abs(X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_i,NVAR1D))
    caj = abs(X_MAGNETICFIELD_FROM_CONSVAR_1D(U2_j,NVAR1D))

    ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
    aPow2i = GAMMA*PRESSURE_FROM_CONSVAR_1D(U2_i,NVAR1D)/MYNEWLINE
             DENSITY_FROM_CONSVAR(U2_i,NVAR1D)
    aPow2j = GAMMA*PRESSURE_FROM_CONSVAR_1D(U2_j,NVAR1D)/MYNEWLINE
             DENSITY_FROM_CONSVAR(U2_j,NVAR1D)
#else
#error "Speed of sound must be implemented!"
#endif

    ! Compute auxiliary quantities
    astPow2i = MAGNETICFIELD_MAGNITUDE_FROM_CONSVAR_1D(U2_i,NVAR1D)/MYNEWLINE
               DENSITY_FROM_CONSVAR(U2_i,NVAR1D) + aPow2i
    astPow2j = MAGNETICFIELD_MAGNITUDE_FROM_CONSVAR_1D(U2_j,NVAR1D)/MYNEWLINE
               DENSITY_FROM_CONSVAR(U2_j,NVAR1D) + aPow2j

    ! Compute the speed of the fast waves
    cfi = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*cai**2)))
    cfj = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*caj**2)))

    ! Scalar dissipation for the Rusanov flux
    d_ij = max( abs(Coeff_ij(1)*uj) + abs(Coeff_ij(1))*cfj,&
                abs(Coeff_ji(1)*ui) + abs(Coeff_ji(1))*cfi )

    ! Compute conservative fluxes
    F_ij = dscale1*(U1_i-U1_j) + dscale2*d_ij*(U2_i-U2_j)

  end subroutine mhd_calcFluxFCTRusanov1d

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDensity1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density in 1D
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
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
    end do

  end subroutine mhd_trafoFluxDensity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDensity1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density in 1D
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
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine mhd_trafoDiffDensity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxEnergy1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the energy in 1D
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
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
    end do

  end subroutine mhd_trafoFluxEnergy1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffEnergy1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the energy in 1D
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
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine mhd_trafoDiffEnergy1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxPressure1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the pressure in 1D
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
    real(DP) :: ui,uj
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(1,1,idx) = (GAMMA-1.0)*&
          (0.5*ui*ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
      DtransformedFluxesAtEdge(1,2,idx) =-(GAMMA-1.0)*&
          (0.5*uj*uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do
    
  end subroutine mhd_trafoFluxPressure1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffPressure1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the pressure in 1D
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
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine mhd_trafoDiffPressure1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxVelocity1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the x-velocity
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
    real(DP) :: ui,uj

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(1,1,idx) =&
          (X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -(X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
    end do
    
  end subroutine mhd_trafoFluxVelocity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffVelocity1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the x-velocity
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
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine mhd_trafoDiffVelocity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxMomentum1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the x-momentum
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
          X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
    end do
    
  end subroutine mhd_trafoFluxMomentum1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffMomentum1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the x-momentum
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
          X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do
    
  end subroutine mhd_trafoDiffMomentum1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDenEng1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 1D
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
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)

      ! Transformed total energy fluxes
      DtransformedFluxesAtEdge(2,1,idx) =&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
    end do

  end subroutine mhd_trafoFluxDenEng1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDenEng1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 1D
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
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)

      ! Transformed total energy difference
      DtransformedDataAtEdge(2,idx) =&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine mhd_trafoDiffDenEng1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDenPre1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 1D
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
    real(DP) :: ui,uj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
           
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)

      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(2,1,idx) = (GAMMA-1.0)*&
          (0.5*ui*ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
      DtransformedFluxesAtEdge(2,2,idx) =-(GAMMA-1.0)*&
          (0.5*uj*uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do

  end subroutine mhd_trafoFluxDenPre1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDenPre1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 1D
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
      DtransformedDataAtEdge(1,idx) = &
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(2,idx) =&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine mhd_trafoDiffDenPre1d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDenPreVel1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density, pressure and velocity in 1D
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
    real(DP) :: ui,uj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)

      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(2,1,idx) =&
          (X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -(X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(3,1,idx) = (GAMMA-1.0)*&
          (0.5*ui*ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
      DtransformedFluxesAtEdge(3,2,idx) =-(GAMMA-1.0)*&
          (0.5*uj*uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do
      
  end subroutine mhd_trafoFluxDenPreVel1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDenPreVel1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density, pressure and velocity in 1D
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
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      
      ! Transformed velocity difference in x-direction
      DtransformedDataAtEdge(2,idx) =&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(3,idx) =&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine mhd_trafoDiffDenPreVel1d_sim

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcBoundaryvalues1d(DbdrNormal, DpointNormal,&
      DbdrValue, ibdrCondType, Du, Du0, istatus)

!<description>
    ! This subroutine computes the boundary values for a given node in 1D
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

  end subroutine mhd_calcBoundaryvalues1d

  !*****************************************************************************

!<subroutine>

  subroutine mhd_hadaptCallbackScalar1d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 1D. The solution vector is assumed
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
            OU_CLASS_WARNING,OU_MODE_STD,'mhd_hadaptCallbackScalar1d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR1D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR1D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR1D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR1D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      do ivar = 1, NVAR1D
        p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR1D+ivar) = &
            0.5*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR1D+ivar)+&
                    p_Dsolution((rcollection%IquickAccess(3)-1)*NVAR1D+ivar))
      end do

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        do ivar = 1, NVAR1D
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR1D+ivar) = &
              p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR1D+ivar)
        end do
      else
        do ivar = 1, NVAR1D
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR1D+ivar) = 0.0
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)

    end select

  end subroutine mhd_hadaptCallbackScalar1d

  !*****************************************************************************

!<subroutine>

  subroutine mhd_hadaptCallbackBlock1d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 1D. The solution vector is assumed
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
      if (rsolution%nblocks .ne. NVAR1D) then
        call output_line('Vector is not in block format!',&
            OU_CLASS_WARNING,OU_MODE_STD,'mhd_hadaptCallbackBlock1d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR1D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR1D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR1D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR1D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      neq = rsolution%NEQ/NVAR1D
      do ivar = 1, NVAR1D
        p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = &
            0.5*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
                    p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(3)) )
      end do

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        neq = rsolution%NEQ/NVAR1D
        do ivar = 1, NVAR1D
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = &
              p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))
        end do
      else
        neq = rsolution%NEQ/NVAR1D
        do ivar = 1, NVAR1D
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = 0.0
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)

    end select

  end subroutine mhd_hadaptCallbackBlock1d

end module mhd_callback1d
