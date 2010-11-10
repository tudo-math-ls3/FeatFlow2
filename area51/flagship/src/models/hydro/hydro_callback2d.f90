!##############################################################################
!# ****************************************************************************
!# <name> hydro_callback2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible Euler/Navier-Stokes equations in 2D.
!#
!# The following callback functions are available:
!#
!# 1.) hydro_calcFluxGal2d_sim
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#
!# 2.) hydro_calcFluxGalNoBdr2d_sim
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#        without assembling the symmetric boundary contribution
!#
!# 3.) hydro_calcFluxScDiss2d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting scalar artificial viscosities
!#
!# 4.) hydro_calcFluxScDissDiSp2d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting scalar artificial viscosities based on
!#        dimensional splitting approach
!#
!# 5.) hydro_calcFluxRoeDiss2d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities
!#
!# 6.) hydro_calcFluxRoeDissDiSp2d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities based on
!#        dimensional splitting approach
!#
!# 7.) hydro_calcFluxRusDiss2d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion based on
!#        dimensional splitting approach
!#
!# 8.) hydro_calcFluxRusDissDiSp2d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion
!#
!# 9.) hydro_calcMatDiagMatD2d_sim
!#     -> Computes local matrix for diagonal entry
!#
!# 10.) hydro_calcMatDiag2d_sim
!#      -> Computes local matrix for diagonal entry
!#
!# 11.) hydro_calcMatGalMatD2d_sim
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 12.) hydro_calcMatGal2d_sim
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 13.) hydro_calcMatScDissMatD2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 14.) hydro_calcMatScDiss2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 15.) hydro_calcMatRoeDissMatD2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 16.) hydro_calcMatRoeDiss2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 17.) hydro_calcMatRusDissMatD2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov artificial viscosities
!#
!# 18.) hydro_calcMatRusDiss2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov flux artificial viscosities
!#
!# 19.) hydro_calcCharacteristics2d_sim
!#      -> Computes characteristic variables
!#
!# 20.) hydro_calcFluxFCTScalarDiss2d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting scalar artificial viscosities
!#
!# 21.) hydro_calcFluxFCTRoeDiss2d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting tensorial artificial viscosities
!#
!# 22.) hydro_calcFluxFCTRusanov2d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting the Rusanov artificial viscosities
!#
!# 23.) hydro_trafoFluxDensity2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 24.) hydro_trafoDiffDensity2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 25.) hydro_trafoFluxEnergy2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 26.) hydro_trafoDiffEnergy2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 27.) hydro_trafoFluxPressure2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 28.) hydro_trafoFluxVelocity2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 29.) hydro_trafoDiffVelocity2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 30.) hydro_trafoFluxMomentum2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 31.) hydro_trafoDiffMomentum2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 32.) hydro_trafoFluxDenEng2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 33.) hydro_trafoDiffDenEng2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 34.) hydro_trafoFluxDenPre2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 35.) hydro_trafoDiffDenPre2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 36.) hydro_trafoFluxDenPreVel2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 37.) hydro_trafoDiffDenPreVel2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure 
!#         and the velocity
!#
!# 38.) hydro_calcBoundaryvalues2d
!#      -> Computes the boundary values for a given node
!#
!# 39.) hydro_hadaptCallbackScalar2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in interleave format
!#
!# 40.) hydro_hadaptCallbackBlock2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in block format
!#
!# 41.) hydro_coeffVectorBdr2d_sim
!#      -> Calculates the coefficients for the linear form in 2D
!#
!# 42.) hydro_coeffMatrixBdr2d_sim
!#      -> Calculates the coefficients for the bilinear form in 2D
!# </purpose>
!##############################################################################

module hydro_callback2d

#include "hydro.h"

  use boundarycondaux
  use collection
  use derivatives
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
  public :: hydro_calcFluxGal2d_sim
  public :: hydro_calcFluxGalNoBdr2d_sim
  public :: hydro_calcFluxScDiss2d_sim
  public :: hydro_calcFluxScDissDiSp2d_sim
  public :: hydro_calcFluxRoeDiss2d_sim
  public :: hydro_calcFluxRoeDissDiSp2d_sim
  public :: hydro_calcFluxRusDiss2d_sim
  public :: hydro_calcFluxRusDissDiSp2d_sim
  public :: hydro_calcMatDiagMatD2d_sim
  public :: hydro_calcMatDiag2d_sim
  public :: hydro_calcMatGalMatD2d_sim
  public :: hydro_calcMatGal2d_sim
  public :: hydro_calcMatScDissMatD2d_sim
  public :: hydro_calcMatScDiss2d_sim
  public :: hydro_calcMatRoeDissMatD2d_sim
  public :: hydro_calcMatRoeDiss2d_sim
  public :: hydro_calcMatRusDissMatD2d_sim
  public :: hydro_calcMatRusDiss2d_sim
  public :: hydro_calcCharacteristics2d_sim
  public :: hydro_calcFluxFCTScalarDiss2d
  public :: hydro_calcFluxFCTRoeDiss2d
  public :: hydro_calcFluxFCTRusanov2d
  public :: hydro_trafoFluxDensity2d_sim
  public :: hydro_trafoFluxEnergy2d_sim
  public :: hydro_trafoFluxPressure2d_sim
  public :: hydro_trafoFluxVelocity2d_sim
  public :: hydro_trafoFluxMomentum2d_sim
  public :: hydro_trafoFluxDenEng2d_sim
  public :: hydro_trafoFluxDenPre2d_sim
  public :: hydro_trafoFluxDenPreVel2d_sim
  public :: hydro_trafoDiffDensity2d_sim
  public :: hydro_trafoDiffEnergy2d_sim
  public :: hydro_trafoDiffPressure2d_sim
  public :: hydro_trafoDiffVelocity2d_sim
  public :: hydro_trafoDiffMomentum2d_sim
  public :: hydro_trafoDiffDenEng2d_sim
  public :: hydro_trafoDiffDenPre2d_sim
  public :: hydro_trafoDiffDenPreVel2d_sim
  public :: hydro_calcBoundaryvalues2d
  public :: hydro_coeffVectorBdr2d_sim
  public :: hydro_hadaptCallbackScalar2d
  public :: hydro_hadaptCallbackBlock2d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGal2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the standard
    ! Galerkin discretisation in 2D.
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
    real(DP), dimension(NVAR2D) :: dF1_i,dF2_i,dF1_j,dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij
#endif
    real(DP) :: pi,pj,ui,vi,uj,vj
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \              / rho*v         \
      ! F1 = | rho*u*u + p   |   and   F2 = | rho*u*v       |
      !      | rho*v*u       |              | rho*v*v + p   |
      !      \ rho*E*u + p*u /              \ rho*E*v + p*v /
      !-------------------------------------------------------------------------
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx)
      
#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui
      dF1_i(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi*ui
      
      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF1_j(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj*uj
      
      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi
      dF2_i(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi*vi
      
      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj
      dF2_j(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj*vj

      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*dF2_i )
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF1_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj*uj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj)
      dF2_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj*vj)
      
      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij)
#endif
    end do

  end subroutine hydro_calcFluxGal2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGalNoBdr2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the TVD
    ! discretisation in 2D. The symmetric boundary contributions
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
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij
    real(DP) :: pi,pj,ui,vi,uj,vj
    integer :: idx

    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \              / rho*v         \
      ! F1 = | rho*u*u + p   |   and   F2 = | rho*u*v       |
      !      | rho*v*u       |              | rho*v*v + p   |
      !      \ rho*E*u + p*u /              \ rho*E*v + p*v /
      !-------------------------------------------------------------------------
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx)

      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF1_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj*uj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj)
      dF2_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj*vj)
      
      ! Assemble symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale *&
          (0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))*dF1_ij+&
           0.5*(DmatrixCoeffsAtEdge(2,1,idx)-DmatrixCoeffsAtEdge(2,2,idx))*dF2_ij)
      DfluxesAtEdge(:,2,idx) = DfluxesAtEdge(:,1,idx)
    end do

  end subroutine hydro_calcFluxGalNoBdr2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxScDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using scalar dissipation.
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
    real(DP), dimension(NVAR2D) :: dF1_i,dF2_i,dF1_j,dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: pi,pj,ui,vi,uj,vj,d_ij,H_ij,q_ij,u_ij,v_ij,anorm,vel_ij,c_ij,aux
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \              / rho*v         \
      ! F1 = | rho*u*u + p   |   and   F2 = | rho*u*v       |
      !      | rho*v*u       |              | rho*v*v + p   |
      !      \ rho*E*u + p*u /              \ rho*E*v + p*v /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui
      dF1_i(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi*ui
      
      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF1_j(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj*uj
      
      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi
      dF2_i(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi*vi
      
      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj
      dF2_j(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj*vj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF1_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj*uj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj)
      dF2_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj*vj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral radius
      !-------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      H_ij = ROE_MEAN_VALUE(MYNEWLINE
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)+pi)/MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+pj)/MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx), aux)
      
      ! Compute auxiliary variables
      vel_ij = u_ij*a(1) + v_ij*a(2)
      q_ij   = 0.5*(u_ij*u_ij+v_ij*v_ij)

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL))
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
                                         DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*dF2_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxScDissDiSp2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using scalar dissipation,
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
    real(DP), dimension(NVAR2D) :: dF1_i,dF2_i,dF1_j,dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: pi,pj,ui,vi,uj,vj,d_ij,H_ij,q_ij,u_ij,v_ij,aux,c_ij
    integer :: idx
    
    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \              / rho*v         \
      ! F1 = | rho*u*u + p   |   and   F2 = | rho*u*v       |
      !      | rho*v*u       |              | rho*v*v + p   |
      !      \ rho*E*u + p*u /              \ rho*E*v + p*v /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui
      dF1_i(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi*ui
      
      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF1_j(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj*uj
      
      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi
      dF2_i(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi*vi
      
      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj
      dF2_j(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj*vj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF1_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj*uj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj)
      dF2_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj*vj)
#endif
      
      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral radius
      !-------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      
      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      H_ij = ROE_MEAN_VALUE(MYNEWLINE
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)+pi)/MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+pj)/MYNEWLINE
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx), aux)
      
      ! Compute auxiliary variable
      q_ij = 0.5*(u_ij*u_ij+v_ij*v_ij)

! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Scalar dissipation for x- and y-direction
      d_ij = ( abs(a(1)*u_ij) + abs(a(1))*c_ij +&
               abs(a(2)*v_ij) + abs(a(2))*c_ij )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*dF2_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxScDissDiSp2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRoeDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using tensorial dissipation.
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
    real(DP), dimension(NVAR2D) :: dF1_i,dF2_i,dF1_j,dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: pi,pj,ui,vi,uj,vj,u_ij,v_ij,H_ij,q_ij,c_ij,c2_ij,vel_ij
    real(DP) :: aux,aux1,aux2,anorm
    real(DP) :: l1,l2,l3,l4,w1,w2,w3,w4
    integer :: idx


    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \              / rho*v         \
      ! F1 = | rho*u*u + p   |   and   F2 = | rho*u*v       |
      !      | rho*v*u       |              | rho*v*v + p   |
      !      \ rho*E*u + p*u /              \ rho*E*v + p*v /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui
      dF1_i(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi*ui
      
      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF1_j(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj*uj
      
      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi
      dF2_i(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi*vi
      
      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj
      dF2_j(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj*vj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF1_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj*uj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj)
      dF2_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj*vj)
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
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)+pi)/MYNEWLINE
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+pj)/MYNEWLINE
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx), aux)

        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2)
        q_ij   = 0.5*(u_ij*u_ij+v_ij*v_ij)

        ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
        c2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(c2_ij)
        
        ! Compute eigenvalues
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)
        
        ! Compute solution difference U_j-U_i
        Diff = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = (GAMMA-1.0)/2.0/c2_ij*(q_ij*Diff(1)-u_ij*Diff(2)-v_ij*Diff(3)+Diff(4))
        aux2 = 0.5*(vel_ij*Diff(1)-a(1)*Diff(2)-a(2)*Diff(3))/c_ij
        
        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/c2_ij)*Diff(1)+(GAMMA-1.0)*(u_ij*Diff(2)+v_ij*Diff(3)-Diff(4))/c2_ij)
        w3 = l3 * (aux1 - aux2)
        w4 = l4 * ((a(1)*v_ij-a(2)*u_ij)*Diff(1)+a(2)*Diff(2)-a(1)*Diff(3))
        
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        Diff(1) = anorm * ( w1 + w2 + w3 )
        Diff(2) = anorm * ( (u_ij-c_ij*a(1))*w1 + u_ij*w2 + (u_ij+c_ij*a(1))*w3 + a(2)*w4 )
        Diff(3) = anorm * ( (v_ij-c_ij*a(2))*w1 + v_ij*w2 + (v_ij+c_ij*a(2))*w3 - a(1)*w4 )
        Diff(4) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3&
                            + (u_ij*a(2)-v_ij*a(1))*w4 )

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*dF2_i + Diff)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij + Diff)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij + Diff)
#endif
      else
        
#ifdef HYDRO_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*dF2_i)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij)
#endif
      end if
    end do

  end subroutine hydro_calcFluxRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRoeDissDiSp2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using tensorial dissipation,
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
    real(DP), dimension(NVAR2D) :: dF1_i,dF2_i,dF1_j,dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff1,Diff2
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: pi,pj,ui,vi,uj,vj,u_ij,v_ij,H_ij,q_ij,c_ij,c2_ij
    real(DP) :: aux,aux1,aux2,anorm
    real(DP) :: l1,l2,l3,l4,w1,w2,w3,w4
    integer :: idx

    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \              / rho*v         \
      ! F1 = | rho*u*u + p   |   and   F2 = | rho*u*v       |
      !      | rho*v*u       |              | rho*v*v + p   |
      !      \ rho*E*u + p*u /              \ rho*E*v + p*v /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui
      dF1_i(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi*ui
      
      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF1_j(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj*uj
      
      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi
      dF2_i(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi*vi
      
      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj
      dF2_j(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj*vj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF1_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj*uj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj)
      dF2_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj*vj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor by Roe
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      if (anorm .gt. SYS_EPSREAL) then

        ! Compute the absolute value
        a = abs(a)
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)+pi)/MYNEWLINE
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+pj)/MYNEWLINE
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx), aux)

        ! Compute auxiliary variables
        q_ij = 0.5*(u_ij*u_ij+v_ij*v_ij)

        ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
        c2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(c2_ij)

        !-----------------------------------------------------------------------
        ! Dimensional splitting: x-direction
        !-----------------------------------------------------------------------
        
        ! Compute eigenvalues
        l1 = abs(u_ij-c_ij)
        l2 = abs(u_ij)
        l3 = abs(u_ij+c_ij)
        l4 = abs(u_ij)
        
        ! Compute solution difference U_j-U_i
        Diff1 = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = (GAMMA-1.0)/2.0/c2_ij*(q_ij*Diff1(1)-u_ij*Diff1(2)-v_ij*Diff1(3)+Diff1(4))
        aux2 = 0.5*(u_ij*Diff1(1)-Diff1(2))/c_ij
        
        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/c2_ij)*Diff1(1)+&
                   (GAMMA-1.0)*(u_ij*Diff1(2)+v_ij*Diff1(3)-Diff1(4))/c2_ij)
        w3 = l3 * (aux1 - aux2)
        w4 = l4 * (v_ij*Diff1(1)-Diff1(3))
        
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        Diff1(1) = a(1) * ( w1 + w2 + w3 )
        Diff1(2) = a(1) * ( (u_ij-c_ij)*w1 + u_ij*w2 + (u_ij+c_ij)*w3 )
        Diff1(3) = a(1) * ( v_ij*w1 + v_ij*w2 + v_ij*w3 - w4 )
        Diff1(4) = a(1) * ( (H_ij-c_ij*u_ij)*w1 + q_ij*w2 + (H_ij+c_ij*u_ij)*w3 - v_ij*w4 )
        
        !-----------------------------------------------------------------------
        ! Dimensional splitting: y-direction
        !-----------------------------------------------------------------------

        ! Compute eigenvalues
        l1 = abs(v_ij-c_ij)
        l2 = abs(v_ij)
        l3 = abs(v_ij+c_ij)
        l4 = abs(v_ij)
        
        ! Compute solution difference U_j-U_i
        Diff2 = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = (GAMMA-1.0)/2.0/c2_ij*(q_ij*Diff2(1)-u_ij*Diff2(2)-v_ij*Diff2(3)+Diff2(4))
        aux2 = 0.5*(v_ij*Diff2(1)-Diff2(3))/c_ij

        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/c2_ij)*Diff2(1)+&
                   (GAMMA-1.0)*(u_ij*Diff2(2)+v_ij*Diff2(3)-Diff2(4))/c2_ij)
        w3 = l3 * (aux1 - aux2)
        w4 = l4 * (-u_ij*Diff2(1)+Diff2(2))
        
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        Diff2(1) = a(2) * ( w1 + w2 + w3 )
        Diff2(2) = a(2) * ( u_ij*w1 + u_ij*w2 + u_ij*w3 + w4 )
        Diff2(3) = a(2) * ( (v_ij-c_ij)*w1 + v_ij*w2 + (v_ij+c_ij)*w3 )
        Diff2(4) = a(2) * ( (H_ij-c_ij*v_ij)*w1 + q_ij*w2 + (H_ij+c_ij*v_ij)*w3 + u_ij*w4 )

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*dF2_i+&
                                           Diff1+Diff2)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij+&
                                            Diff1+Diff2)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij+&
                                            Diff1+Diff2)
#endif
      else

#ifdef HYDRO_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*dF2_i)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij)
#endif
      end if
    end do

  end subroutine hydro_calcFluxRoeDissDiSp2d_sim
 
  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRusDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using the Rusanov dissipation.
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
    real(DP), dimension(NVAR2D) :: dF1_i,dF2_i,dF1_j,dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: pi,pj,ui,vi,uj,vj,d_ij,ci,cj,Ei,Ej
    integer :: idx

    do idx = 1, size(DfluxesAtEdge,3)

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \              / rho*v         \
      ! F1 = | rho*u*u + p   |   and   F2 = | rho*u*v       |
      !      | rho*v*u       |              | rho*v*v + p   |
      !      \ rho*E*u + p*u /              \ rho*E*v + p*v /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui
      dF1_i(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi*ui
      
      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF1_j(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj*uj
      
      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi
      dF2_i(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi*vi
      
      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj
      dF2_j(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj*vj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF1_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj*uj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj)
      dF2_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj*vj)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov type
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*(ui*ui+vi*vi)), SYS_EPSREAL))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*(uj*uj+vj*vj)), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Scalar dissipation for the Rusanov flux
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj+&
                      DmatrixCoeffsAtEdge(2,1,idx)*vj)+&
                 sqrt(DmatrixCoeffsAtEdge(1,1,idx)**2+&
                      DmatrixCoeffsAtEdge(2,1,idx)**2)*cj,&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui+&
                      DmatrixCoeffsAtEdge(2,2,idx)*vi)+&
                 sqrt(DmatrixCoeffsAtEdge(1,2,idx)**2+&
                      DmatrixCoeffsAtEdge(2,2,idx)**2)*ci )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------
      
#ifdef HYDRO_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*dF2_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxRusDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRusDissDiSp2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using the Rusanov dissipation,
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
    real(DP), dimension(NVAR2D) :: dF1_i,dF2_i,dF1_j,dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: pi,pj,ui,vi,uj,vj,d_ij,ci,cj,Ei,Ej
    integer :: idx
    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u         \              / rho*v         \
      ! F1 = | rho*u*u + p   |   and   F2 = | rho*u*v       |
      !      | rho*v*u       |              | rho*v*v + p   |
      !      \ rho*E*u + p*u /              \ rho*E*v + p*v /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui
      dF1_i(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi*ui
      
      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF1_j(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj*uj
      
      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi
      dF2_i(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi*vi
      
      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj
      dF2_j(4) = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj*vj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi)-&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj)
      dF1_ij(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF1_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi*ui)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj*uj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)-&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi-&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi)-&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj)
      dF2_ij(4) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi*vi)-&
                  (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj*vj)
#endif
      
      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov type
      !-------------------------------------------------------------------------

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*(ui*ui+vi*vi)), SYS_EPSREAL))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*(uj*uj+vj*vj)), SYS_EPSREAL))
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
                  abs(DmatrixCoeffsAtEdge(2,2,idx))*ci )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))
      
      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*dF2_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxRusDissDiSp2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiagMatD2d_sim(DdataAtNode, DmatrixCoeffsAtNode,&
      IverticesAtNode, dscale, DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 2D
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
    real(DP) :: ui,vi
    integer :: inode

    do inode = 1, size(DcoefficientsAtNode,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_1T_FROM_CONSVAR(DdataAtNode,NVAR2D,inode)
      vi = Y_VELOCITY_1T_FROM_CONSVAR(DdataAtNode,NVAR2D,inode)
      
#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ii = diag(A_i)*C_{ii}$
      DcoefficientsAtNode(1,1,inode) = 0.0
      DcoefficientsAtNode(2,1,inode) = dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtNode(1,inode)+&
                                                             vi*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode(3,1,inode) = dscale * (            ui*DmatrixCoeffsAtNode(1,inode)+&
                                                 (3.0-GAMMA)*vi*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode(4,1,inode) = dscale * (GAMMA*(ui*DmatrixCoeffsAtNode(1,inode)+&
                                                        vi*DmatrixCoeffsAtNode(2,inode)))
#else
      ! Compute Galerkin coefficient $K_ii = -diag(A_i)*C_{ii}$
      DcoefficientsAtNode(1,1,inode) = 0.0
      DcoefficientsAtNode(2,1,inode) = -dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtNode(1,inode)+&
                                                              vi*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode(3,1,inode) = -dscale * (            ui*DmatrixCoeffsAtNode(1,inode)+&
                                                  (3.0-GAMMA)*vi*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode(4,1,inode) = -dscale * (GAMMA*(ui*DmatrixCoeffsAtNode(1,inode)+&
                                                         vi*DmatrixCoeffsAtNode(2,inode)))
#endif
    end do

  end subroutine hydro_calcMatDiagMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiag2d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 2D
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
    real(DP) :: ui,vi,qi,Ei,uvi,uPow2i,vPow2i,aux
    integer :: inode

    do inode = 1, size(DcoefficientsAtNode,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_1T_FROM_CONSVAR(DdataAtNode,NVAR2D,inode)
      vi = Y_VELOCITY_1T_FROM_CONSVAR(DdataAtNode,NVAR2D,inode)
      Ei = TOTAL_ENERGY_1T_FROM_CONSVAR(DdataAtNode,NVAR2D,inode)
      uvi = ui*vi; qi = ui*ui+vi*vi; uPow2i = ui*ui; vPow2i = vi*vi
      aux = ui*DmatrixCoeffsAtNode(1,inode)+vi*DmatrixCoeffsAtNode(2,inode)

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      DcoefficientsAtNode( 1,1,inode) = 0.0
      DcoefficientsAtNode( 2,1,inode) = dscale * (((GAMMA-1.0)/2.0*qi-uPow2i)*DmatrixCoeffsAtNode(1,inode)-&
                                                                          uvi*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode( 3,1,inode) = dscale * (((GAMMA-1.0)/2.0*qi-vPow2i)*DmatrixCoeffsAtNode(2,inode)-&
                                                                          uvi*DmatrixCoeffsAtNode(1,inode))
      DcoefficientsAtNode( 4,1,inode) = dscale * ((GAMMA-1.0)*qi-GAMMA*Ei)*aux
      
      DcoefficientsAtNode( 5,1,inode) = dscale * DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode( 6,1,inode) = dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtNode(1,inode)+&
                                                              vi*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode( 7,1,inode) = dscale * (            vi*DmatrixCoeffsAtNode(1,inode)-&
                                                  (GAMMA-1.0)*ui*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode( 8,1,inode) = dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtNode(1,inode)-&
                                                  (GAMMA-1.0)*ui*aux)
      
      DcoefficientsAtNode( 9,1,inode) = dscale * DmatrixCoeffsAtNode(2,inode)
      DcoefficientsAtNode(10,1,inode) = dscale * (            ui*DmatrixCoeffsAtNode(2,inode)-&
                                                  (GAMMA-1.0)*vi*DmatrixCoeffsAtNode(1,inode))
      DcoefficientsAtNode(11,1,inode) = dscale * (            ui*DmatrixCoeffsAtNode(1,inode)+&
                                                  (3.0-GAMMA)*vi*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode(12,1,inode) = dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtNode(2,inode)-&
                                                  (GAMMA-1.0)*vi*aux)
      
      DcoefficientsAtNode(13,1,inode) = 0.0
      DcoefficientsAtNode(14,1,inode) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(15,1,inode) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtNode(2,inode)
      DcoefficientsAtNode(16,1,inode) = dscale * (GAMMA*(ui*DmatrixCoeffsAtNode(1,inode)+&
                                                         vi*DmatrixCoeffsAtNode(2,inode)))
#else
      ! Compute Galerkin coefficient $K_ii = -A_i*C_{ii}$
      DcoefficientsAtNode( 1,1,inode) = 0.0
      DcoefficientsAtNode( 2,1,inode) =-dscale * (((GAMMA-1.0)/2.0*qi-uPow2i)*DmatrixCoeffsAtNode(1,inode)-&
                                                                          uvi*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode( 3,1,inode) =-dscale * (((GAMMA-1.0)/2.0*qi-vPow2i)*DmatrixCoeffsAtNode(2,inode)-&
                                                                          uvi*DmatrixCoeffsAtNode(1,inode))
      DcoefficientsAtNode( 4,1,inode) =-dscale * ((GAMMA-1.0)*qi-GAMMA*Ei)*aux
      
      DcoefficientsAtNode( 5,1,inode) =-dscale * DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode( 6,1,inode) =-dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtNode(1,inode)+&
                                                              vi*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode( 7,1,inode) =-dscale * (            vi*DmatrixCoeffsAtNode(1,inode)-&
                                                  (GAMMA-1.0)*ui*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode( 8,1,inode) =-dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtNode(1,inode)-&
                                                  (GAMMA-1.0)*ui*aux)
      
      DcoefficientsAtNode( 9,1,inode) =-dscale * DmatrixCoeffsAtNode(2,inode)
      DcoefficientsAtNode(10,1,inode) =-dscale * (            ui*DmatrixCoeffsAtNode(2,inode)-&
                                                  (GAMMA-1.0)*vi*DmatrixCoeffsAtNode(1,inode))
      DcoefficientsAtNode(11,1,inode) =-dscale * (            ui*DmatrixCoeffsAtNode(1,inode)+&
                                                  (3.0-GAMMA)*vi*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode(12,1,inode) =-dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtNode(2,inode)-&
                                                  (GAMMA-1.0)*vi*aux)
      
      DcoefficientsAtNode(13,1,inode) = 0.0
      DcoefficientsAtNode(14,1,inode) =-dscale * (GAMMA-1.0)*DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(15,1,inode) =-dscale * (GAMMA-1.0)*DmatrixCoeffsAtNode(2,inode)
      DcoefficientsAtNode(16,1,inode) =-dscale * (GAMMA*(ui*DmatrixCoeffsAtNode(1,inode)+&
                                                         vi*DmatrixCoeffsAtNode(2,inode)))
#endif
    end do

  end subroutine hydro_calcMatDiag2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGalMatD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 2D
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
    real(DP) :: ui,uj,vi,vj
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                           vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(3,2,idx) = dscale * (            uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                               (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(4,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                      vj*DmatrixCoeffsAtEdge(2,2,idx)))
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                           vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(3,3,idx) = dscale * (            ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                               (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(4,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                      vi*DmatrixCoeffsAtEdge(2,1,idx)))
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = -dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                            vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(3,2,idx) = -dscale * (            uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(4,2,idx) = -dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                       vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = -dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                            vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(3,3,idx) = -dscale * (            ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(4,3,idx) = -dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                       vi*DmatrixCoeffsAtEdge(2,2,idx)))
#endif
    end do

  end subroutine hydro_calcMatGalMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGal2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 2D
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
    real(DP) :: Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj,uPow2i,uPow2j,vPow2i,vPow2j,aux1,aux2
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      Ei = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uvi = ui*vi; qi = ui*ui+vi*vi; uPow2i = ui*ui; vPow2i = vi*vi

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      Ej = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      uvj = uj*vj; qj = uj*uj+vj*vj; uPow2j = uj*uj; vPow2j = vj*vj
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP
      aux1 = uj*DmatrixCoeffsAtEdge(1,2,idx)+vj*DmatrixCoeffsAtEdge(2,2,idx)
      aux2 = ui*DmatrixCoeffsAtEdge(1,1,idx)+vi*DmatrixCoeffsAtEdge(2,1,idx)

      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      DcoefficientsAtEdge( 1,2,idx) = 0.0
      DcoefficientsAtEdge( 2,2,idx) = dscale * (((GAMMA-1.0)/2.0*qj-uPow2j)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                                        uvj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 3,2,idx) = dscale * (((GAMMA-1.0)/2.0*qj-vPow2j)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                                        uvj*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge( 4,2,idx) = dscale * ((GAMMA-1.0)*qj-GAMMA*Ej)*aux1
      
      DcoefficientsAtEdge( 5,2,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge( 6,2,idx) = dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                            vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 7,2,idx) = dscale * (            vj*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                (GAMMA-1.0)*uj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 8,2,idx) = dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                (GAMMA-1.0)*uj*aux1)
      
      DcoefficientsAtEdge( 9,2,idx) = dscale * DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(10,2,idx) = dscale * (            uj*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                (GAMMA-1.0)*vj*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge(11,2,idx) = dscale * (            uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(12,2,idx) = dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                (GAMMA-1.0)*vj*aux1)
      
      DcoefficientsAtEdge(13,2,idx) = 0.0
      DcoefficientsAtEdge(14,2,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(15,2,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(16,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                       vj*DmatrixCoeffsAtEdge(2,2,idx)))
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      DcoefficientsAtEdge( 1,3,idx) = 0.0
      DcoefficientsAtEdge( 2,3,idx) = dscale * (((GAMMA-1.0)*qi-uPow2i)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                                    uvi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 3,3,idx) = dscale * (((GAMMA-1.0)*qi-vPow2i)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                                    uvi*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge( 4,3,idx) = dscale * ((GAMMA-1.0)*qi-GAMMA*Ei)*aux2
      
      DcoefficientsAtEdge( 5,3,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge( 6,3,idx) = dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                            vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 7,3,idx) = dscale * (            vi*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                (GAMMA-1.0)*ui*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 8,3,idx) = dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                (GAMMA-1.0)*ui*aux2)
      
      DcoefficientsAtEdge( 9,3,idx) = dscale * DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(10,3,idx) = dscale * (            ui*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                (GAMMA-1.0)*vi*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge(11,3,idx) = dscale * (            ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(12,3,idx) = dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                (GAMMA-1.0)*vi*aux2)
      
      DcoefficientsAtEdge(13,3,idx) = 0.0
      DcoefficientsAtEdge(14,3,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(15,3,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(16,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                       vi*DmatrixCoeffsAtEdge(2,1,idx)))
#else
      aux1 = uj*DmatrixCoeffsAtEdge(1,1,idx)+vj*DmatrixCoeffsAtEdge(2,1,idx)
      aux2 = ui*DmatrixCoeffsAtEdge(1,2,idx)+vi*DmatrixCoeffsAtEdge(2,2,idx)

      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      DcoefficientsAtEdge( 1,2,idx) = 0.0
      DcoefficientsAtEdge( 2,2,idx) = -dscale * (((GAMMA-1.0)/2.0*qj-uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                                         uvj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 3,2,idx) = -dscale * (((GAMMA-1.0)/2.0*qj-vPow2j)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                                         uvj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge( 4,2,idx) = -dscale * ((GAMMA-1.0)*qj-GAMMA*Ej)*aux1
      
      DcoefficientsAtEdge( 5,2,idx) = -dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge( 6,2,idx) = -dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                             vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 7,2,idx) = -dscale * (            vj*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                 (GAMMA-1.0)*uj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 8,2,idx) = -dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                 (GAMMA-1.0)*uj*aux1)
      
      DcoefficientsAtEdge( 9,2,idx) = -dscale * DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(10,2,idx) = -dscale * (            uj*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                 (GAMMA-1.0)*vj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge(11,2,idx) = -dscale * (            uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                 (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(12,2,idx) = -dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                 (GAMMA-1.0)*vj*aux1)
      
      DcoefficientsAtEdge(13,2,idx) = 0.0
      DcoefficientsAtEdge(14,2,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(15,2,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(16,2,idx) = -dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                        vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      DcoefficientsAtEdge( 1,3,idx) = 0.0
      DcoefficientsAtEdge( 2,3,idx) = -dscale * (((GAMMA-1.0)*qi-uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                                     uvi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 3,3,idx) = -dscale * (((GAMMA-1.0)*qi-vPow2i)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                                     uvi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge( 4,3,idx) = -dscale * ((GAMMA-1.0)*qi-GAMMA*Ei)*aux2
      
      DcoefficientsAtEdge( 5,3,idx) = -dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge( 6,3,idx) = -dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                             vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 7,3,idx) = -dscale * (            vi*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                 (GAMMA-1.0)*ui*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 8,3,idx) = -dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                 (GAMMA-1.0)*ui*aux2)
      
      DcoefficientsAtEdge( 9,3,idx) = -dscale * DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(10,3,idx) = -dscale * (            ui*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                 (GAMMA-1.0)*vi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge(11,3,idx) = -dscale * (            ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                 (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(12,3,idx) = -dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                 (GAMMA-1.0)*vi*aux2)
      
      DcoefficientsAtEdge(13,3,idx) = 0.0
      DcoefficientsAtEdge(14,3,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(15,3,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(16,3,idx) = -dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                        vi*DmatrixCoeffsAtEdge(2,2,idx)))
#endif
    end do
      
  end subroutine hydro_calcMatGal2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDissMatD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies scalar artificial viscosities in 2D
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
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: anorm,aux,H_ij,q_ij,ui,uj,u_ij,vi,vj,v_ij,c_ij,vel_ij
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                           vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(3,2,idx) = dscale * (           uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                               (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(4,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                      vj*DmatrixCoeffsAtEdge(2,2,idx)))
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                           vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(3,3,idx) = dscale * (            ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                               (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(4,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                      vi*DmatrixCoeffsAtEdge(2,1,idx)))
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = -dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                            vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(3,2,idx) = -dscale * (            uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(4,2,idx) = -dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                       vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = -dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                            vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(3,3,idx) = -dscale * (            ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(4,3,idx) = -dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                       vi*DmatrixCoeffsAtEdge(2,2,idx)))
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx), aux)
               
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2)
        q_ij   = 0.5*(u_ij*u_ij+v_ij*v_ij)

        ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
        c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif

        ! Compute scalar dissipation
        DcoefficientsAtEdge(:,1,idx) = dscale * (abs(vel_ij) + anorm*c_ij)

      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0

      end if
    end do

  end subroutine hydro_calcMatScDissMatD2d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDiss2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies scalar artificial viscosities in 2D
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
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: anorm,aux,aux1,aux2,H_ij,q_ij,u_ij,v_ij,vel_ij,c_ij
    real(DP) :: Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj,uPow2i,uPow2j,vPow2i,vPow2j
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      Ei = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uvi = ui*vi; qi = ui*ui+vi*vi; uPow2i = ui*ui; vPow2i = vi*vi

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      Ej = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      uvj = uj*vj; qj = uj*uj+vj*vj; uPow2j = uj*uj; vPow2j = vj*vj
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP
      aux1 = uj*DmatrixCoeffsAtEdge(1,2,idx)+vj*DmatrixCoeffsAtEdge(2,2,idx)
      aux2 = ui*DmatrixCoeffsAtEdge(1,1,idx)+vi*DmatrixCoeffsAtEdge(2,1,idx)
      
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      DcoefficientsAtEdge( 1,2,idx) = 0.0
      DcoefficientsAtEdge( 2,2,idx) = dscale * (((GAMMA-1.0)/2.0*qj-uPow2j)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                                        uvj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 3,2,idx) = dscale * (((GAMMA-1.0)/2.0*qj-vPow2j)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                                        uvj*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge( 4,2,idx) = dscale * ((GAMMA-1.0)*qj-GAMMA*Ej)*aux1
      
      DcoefficientsAtEdge( 5,2,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge( 6,2,idx) = dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                            vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 7,2,idx) = dscale * (            vj*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                (GAMMA-1.0)*uj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 8,2,idx) = dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                (GAMMA-1.0)*uj*aux1)
      
      DcoefficientsAtEdge( 9,2,idx) = dscale * DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(10,2,idx) = dscale * (            uj*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                (GAMMA-1.0)*vj*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge(11,2,idx) = dscale * (            uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(12,2,idx) = dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                (GAMMA-1.0)*vj*aux1)
      
      DcoefficientsAtEdge(13,2,idx) = 0.0
      DcoefficientsAtEdge(14,2,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(15,2,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(16,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                       vj*DmatrixCoeffsAtEdge(2,2,idx)))
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      DcoefficientsAtEdge( 1,3,idx) = 0.0
      DcoefficientsAtEdge( 2,3,idx) = dscale * (((GAMMA-1.0)*qi-uPow2i)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                                    uvi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 3,3,idx) = dscale * (((GAMMA-1.0)*qi-vPow2i)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                                    uvi*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge( 4,3,idx) = dscale * ((GAMMA-1.0)*qi-GAMMA*Ei)*aux2
      
      DcoefficientsAtEdge( 5,3,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge( 6,3,idx) = dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                            vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 7,3,idx) = dscale * (            vi*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                (GAMMA-1.0)*ui*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 8,3,idx) = dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                (GAMMA-1.0)*ui*aux2)
      
      DcoefficientsAtEdge( 9,3,idx) = dscale * DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(10,3,idx) = dscale * (            ui*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                (GAMMA-1.0)*vi*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge(11,3,idx) = dscale * (            ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(12,3,idx) = dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                (GAMMA-1.0)*vi*aux2)
      
      DcoefficientsAtEdge(13,3,idx) = 0.0
      DcoefficientsAtEdge(14,3,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(15,3,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(16,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                       vi*DmatrixCoeffsAtEdge(2,1,idx)))
#else
      aux1 = uj*DmatrixCoeffsAtEdge(1,1,idx)+vj*DmatrixCoeffsAtEdge(2,1,idx)
      aux2 = ui*DmatrixCoeffsAtEdge(1,2,idx)+vi*DmatrixCoeffsAtEdge(2,2,idx)
      
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      DcoefficientsAtEdge( 1,2,idx) = 0.0
      DcoefficientsAtEdge( 2,2,idx) = -dscale * (((GAMMA-1.0)/2.0*qj-uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                                         uvj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 3,2,idx) = -dscale * (((GAMMA-1.0)/2.0*qj-vPow2j)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                                         uvj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge( 4,2,idx) = -dscale * ((GAMMA-1.0)*qj-GAMMA*Ej)*aux1
      
      DcoefficientsAtEdge( 5,2,idx) = -dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge( 6,2,idx) = -dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                             vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 7,2,idx) = -dscale * (            vj*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                 (GAMMA-1.0)*uj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 8,2,idx) = -dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                 (GAMMA-1.0)*uj*aux1)
      
      DcoefficientsAtEdge( 9,2,idx) = -dscale * DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(10,2,idx) = -dscale * (            uj*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                 (GAMMA-1.0)*vj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge(11,2,idx) = -dscale * (            uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                 (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(12,2,idx) = -dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                 (GAMMA-1.0)*vj*aux1)
      
      DcoefficientsAtEdge(13,2,idx) = 0.0
      DcoefficientsAtEdge(14,2,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(15,2,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(16,2,idx) = -dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                        vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      DcoefficientsAtEdge( 1,3,idx) = 0.0
      DcoefficientsAtEdge( 2,3,idx) = -dscale * (((GAMMA-1.0)*qi-uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                                     uvi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 3,3,idx) = -dscale * (((GAMMA-1.0)*qi-vPow2i)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                                     uvi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge( 4,3,idx) = -dscale * ((GAMMA-1.0)*qi-GAMMA*Ei)*aux2
      
      DcoefficientsAtEdge( 5,3,idx) = -dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge( 6,3,idx) = -dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                             vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 7,3,idx) = -dscale * (            vi*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                 (GAMMA-1.0)*ui*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 8,3,idx) = -dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                 (GAMMA-1.0)*ui*aux2)
      
      DcoefficientsAtEdge( 9,3,idx) = -dscale * DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(10,3,idx) = -dscale * (            ui*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                 (GAMMA-1.0)*vi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge(11,3,idx) = -dscale * (            ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                 (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(12,3,idx) = -dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                 (GAMMA-1.0)*vi*aux2)
      
      DcoefficientsAtEdge(13,3,idx) = 0.0
      DcoefficientsAtEdge(14,3,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(15,3,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(16,3,idx) = -dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                        vi*DmatrixCoeffsAtEdge(2,2,idx)))
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx), aux)
        
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2)
        q_ij   = 0.5*(u_ij*u_ij+v_ij*v_ij)


! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
        c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif

        ! Compute scalar dissipation
        aux = dscale * (abs(vel_ij) + anorm*c_ij)
              
        DcoefficientsAtEdge( 1,1,idx) = aux
        DcoefficientsAtEdge( 6,1,idx) = aux
        DcoefficientsAtEdge(11,1,idx) = aux
        DcoefficientsAtEdge(16,1,idx) = aux

      end if
    end do

  end subroutine hydro_calcMatScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDissMatD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 2D
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
    real(DP), dimension(NVAR2D,NVAR2D) :: R_ij,L_ij
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: aux,H_ij,q_ij,ui,uj,u_ij,vi,vj,v_ij
    real(DP) :: l1,l2,l3,l4,anorm,c1,c2,c_ij,cPow2_ij,vel
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                           vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(3,2,idx) = dscale * (            uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                               (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(4,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                      vj*DmatrixCoeffsAtEdge(2,2,idx)))
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                           vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(3,3,idx) = dscale * (            ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                               (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(4,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                      vi*DmatrixCoeffsAtEdge(2,1,idx)))
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = -dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                            vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(3,2,idx) = -dscale * (            uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(4,2,idx) = -dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                       vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = -dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                            vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(3,3,idx) = -dscale * (            ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(4,3,idx) = -dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                       vi*DmatrixCoeffsAtEdge(2,2,idx)))
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL) then

        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx), aux)
        
        ! Compute auxiliary values
        c1    = a(1)/anorm
        c2    = a(2)/anorm
        vel   = c1*u_ij+c2*v_ij
        q_ij  = 0.5*(u_ij*u_ij+v_ij*v_ij)
#ifdef THERMALLY_IDEAL_GAS
        cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij= sqrt(cPow2_ij)
        
        ! Diagonal matrix of eigenvalues
        l1 = abs(vel-c_ij)
        l2 = abs(vel)
        l3 = abs(vel+c_ij)
        l4 = abs(vel)

        ! Matrix of right eigenvectors
        R_ij(1,1) =  l1
        R_ij(2,1) =  l1*(u_ij-c_ij*c1)
        R_ij(3,1) =  l1*(v_ij-c_ij*c2)
        R_ij(4,1) =  l1*(H_ij-c_ij*vel)
        
        R_ij(1,2) =  l2
        R_ij(2,2) =  l2*u_ij
        R_ij(3,2) =  l2*v_ij
        R_ij(4,2) =  l2*q_ij
        
        R_ij(1,3) =  l3
        R_ij(2,3) =  l3*(u_ij+c_ij*c1)
        R_ij(3,3) =  l3*(v_ij+c_ij*c2)
        R_ij(4,3) =  l3*(H_ij+c_ij*vel)
        
        R_ij(1,4) =  0.0
        R_ij(2,4) =  l4*c2
        R_ij(3,4) = -l4*c1
        R_ij(4,4) =  l4*(u_ij*c2-v_ij*c1)
        
        ! Matrix of left eigenvectors
        L_ij(1,1) = 0.5*((GAMMA-1.0)*q_ij+c_ij*vel)/cPow2_ij
        L_ij(2,1) = (cPow2_ij-(GAMMA-1.0)*q_ij)/cPow2_ij
        L_ij(3,1) = 0.5*((GAMMA-1.0)*q_ij-c_ij*vel)/cPow2_ij
        L_ij(4,1) = v_ij*c1-u_ij*c2
        
        L_ij(1,2) = 0.5*(-(GAMMA-1.0)*u_ij-c_ij*c1)/cPow2_ij
        L_ij(2,2) = (GAMMA-1.0)*u_ij/cPow2_ij
        L_ij(3,2) = 0.5*(-(GAMMA-1.0)*u_ij+c_ij*c1)/cPow2_ij
        L_ij(4,2) = c2

        L_ij(1,3) = 0.5*(-(GAMMA-1.0)*v_ij-c_ij*c2)/cPow2_ij
        L_ij(2,3) = (GAMMA-1.0)*v_ij/cPow2_ij
        L_ij(3,3) = 0.5*(-(GAMMA-1.0)*v_ij+c_ij*c2)/cPow2_ij
        L_ij(4,3) = -c1
        
        L_ij(1,4) =  (GAMMA-1.0)/2.0/cPow2_ij
        L_ij(2,4) = -(GAMMA-1.0)/cPow2_ij
        L_ij(3,4) =  (GAMMA-1.0)/2.0/cPow2_ij
        L_ij(4,4) =  0.0
        
        ! Include scaling parameter
        anorm = dscale*anorm
        
        ! Compute tensorial dissipation D_ij = diag(R_ij*|Lbd_ij|*L_ij)*I
        DcoefficientsAtEdge(1,1,idx) = anorm*( R_ij(1,1)*L_ij(1,1)+&
            R_ij(1,2)*L_ij(2,1)+R_ij(1,3)*L_ij(3,1)+R_ij(1,4)*L_ij(4,1)  )
        DcoefficientsAtEdge(2,1,idx) = anorm*( R_ij(2,1)*L_ij(1,2)+&
            R_ij(2,2)*L_ij(2,2)+R_ij(2,3)*L_ij(3,2)+R_ij(2,4)*L_ij(4,2)  )
        DcoefficientsAtEdge(3,1,idx) = anorm*( R_ij(3,1)*L_ij(1,3)+&
            R_ij(3,2)*L_ij(2,3)+R_ij(3,3)*L_ij(3,3)+R_ij(3,4)*L_ij(4,3)  )
        DcoefficientsAtEdge(4,1,idx) = anorm*( R_ij(4,1)*L_ij(1,4)+&
            R_ij(4,2)*L_ij(2,4)+R_ij(4,3)*L_ij(3,4)+R_ij(4,4)*L_ij(4,4)  )
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0
        
      end if
    end do

  end subroutine hydro_calcMatRoeDissMatD2d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDiss2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 2D
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
    real(DP), dimension(NVAR2D,NVAR2D) :: R_ij,L_ij
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: anorm,aux,H_ij,q_ij,u_ij,v_ij,vel,c1,c2,cPow2_ij,c_ij,l1,l2,l3,l4
    real(DP) :: Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj,uPow2i,uPow2j,vPow2i,vPow2j,aux1,aux2
    integer :: idx,i,j,k

    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      Ei = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uvi = ui*vi; qi = ui*ui+vi*vi; uPow2i = ui*ui; vPow2i = vi*vi

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      Ej = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      uvj = uj*vj; qj = uj*uj+vj*vj; uPow2j = uj*uj; vPow2j = vj*vj
      
#ifdef HYDRO_USE_IBP
      aux1 = uj*DmatrixCoeffsAtEdge(1,2,idx)+vj*DmatrixCoeffsAtEdge(2,2,idx)
      aux2 = ui*DmatrixCoeffsAtEdge(1,1,idx)+vi*DmatrixCoeffsAtEdge(2,1,idx)

      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      DcoefficientsAtEdge( 1,2,idx) = 0.0
      DcoefficientsAtEdge( 2,2,idx) = dscale * (((GAMMA-1.0)/2.0*qj-uPow2j)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                                        uvj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 3,2,idx) = dscale * (((GAMMA-1.0)/2.0*qj-vPow2j)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                                        uvj*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge( 4,2,idx) = dscale * ((GAMMA-1.0)*qj-GAMMA*Ej)*aux1
      
      DcoefficientsAtEdge( 5,2,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge( 6,2,idx) = dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                            vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 7,2,idx) = dscale * (            vj*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                (GAMMA-1.0)*uj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 8,2,idx) = dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                (GAMMA-1.0)*uj*aux1)
      
      DcoefficientsAtEdge( 9,2,idx) = dscale * DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(10,2,idx) = dscale * (            uj*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                (GAMMA-1.0)*vj*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge(11,2,idx) = dscale * (            uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(12,2,idx) = dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                (GAMMA-1.0)*vj*aux1)
      
      DcoefficientsAtEdge(13,2,idx) = 0.0
      DcoefficientsAtEdge(14,2,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(15,2,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(16,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                       vj*DmatrixCoeffsAtEdge(2,2,idx)))
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      DcoefficientsAtEdge( 1,3,idx) = 0.0
      DcoefficientsAtEdge( 2,3,idx) = dscale * (((GAMMA-1.0)*qi-uPow2i)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                                    uvi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 3,3,idx) = dscale * (((GAMMA-1.0)*qi-vPow2i)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                                    uvi*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge( 4,3,idx) = dscale * ((GAMMA-1.0)*qi-GAMMA*Ei)*aux2
      
      DcoefficientsAtEdge( 5,3,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge( 6,3,idx) = dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                            vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 7,3,idx) = dscale * (            vi*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                (GAMMA-1.0)*ui*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 8,3,idx) = dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                (GAMMA-1.0)*ui*aux2)
      
      DcoefficientsAtEdge( 9,3,idx) = dscale * DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(10,3,idx) = dscale * (            ui*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                (GAMMA-1.0)*vi*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge(11,3,idx) = dscale * (            ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(12,3,idx) = dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                (GAMMA-1.0)*vi*aux2)
      
      DcoefficientsAtEdge(13,3,idx) = 0.0
      DcoefficientsAtEdge(14,3,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(15,3,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(16,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                       vi*DmatrixCoeffsAtEdge(2,1,idx)))
#else
      aux1 = uj*DmatrixCoeffsAtEdge(1,1,idx)+vj*DmatrixCoeffsAtEdge(2,1,idx)
      aux2 = ui*DmatrixCoeffsAtEdge(1,2,idx)+vi*DmatrixCoeffsAtEdge(2,2,idx)

      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      DcoefficientsAtEdge( 1,2,idx) = 0.0
      DcoefficientsAtEdge( 2,2,idx) = -dscale * (((GAMMA-1.0)/2.0*qj-uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                                         uvj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 3,2,idx) = -dscale * (((GAMMA-1.0)/2.0*qj-vPow2j)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                                         uvj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge( 4,2,idx) = -dscale * ((GAMMA-1.0)*qj-GAMMA*Ej)*aux1
      
      DcoefficientsAtEdge( 5,2,idx) = -dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge( 6,2,idx) = -dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                             vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 7,2,idx) = -dscale * (            vj*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                 (GAMMA-1.0)*uj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 8,2,idx) = -dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                 (GAMMA-1.0)*uj*aux1)
      
      DcoefficientsAtEdge( 9,2,idx) = -dscale * DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(10,2,idx) = -dscale * (            uj*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                 (GAMMA-1.0)*vj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge(11,2,idx) = -dscale * (            uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                 (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(12,2,idx) = -dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                 (GAMMA-1.0)*vj*aux1)
      
      DcoefficientsAtEdge(13,2,idx) = 0.0
      DcoefficientsAtEdge(14,2,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(15,2,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(16,2,idx) = -dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                        vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      DcoefficientsAtEdge( 1,3,idx) = 0.0
      DcoefficientsAtEdge( 2,3,idx) = -dscale * (((GAMMA-1.0)*qi-uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                                     uvi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 3,3,idx) = -dscale * (((GAMMA-1.0)*qi-vPow2i)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                                     uvi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge( 4,3,idx) = -dscale * ((GAMMA-1.0)*qi-GAMMA*Ei)*aux2
      
      DcoefficientsAtEdge( 5,3,idx) = -dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge( 6,3,idx) = -dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                             vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 7,3,idx) = -dscale * (            vi*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                 (GAMMA-1.0)*ui*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 8,3,idx) = -dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                 (GAMMA-1.0)*ui*aux2)
      
      DcoefficientsAtEdge( 9,3,idx) = -dscale * DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(10,3,idx) = -dscale * (            ui*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                 (GAMMA-1.0)*vi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge(11,3,idx) = -dscale * (            ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                 (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(12,3,idx) = -dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                 (GAMMA-1.0)*vi*aux2)
      
      DcoefficientsAtEdge(13,3,idx) = 0.0
      DcoefficientsAtEdge(14,3,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(15,3,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(16,3,idx) = -dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                        vi*DmatrixCoeffsAtEdge(2,2,idx)))
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx), aux)
        
        ! Compute auxiliary variables
        c1    = a(1)/anorm
        c2    = a(2)/anorm
        vel   = c1*u_ij+c2*v_ij
        q_ij  = 0.5*(u_ij*u_ij+v_ij*v_ij)
#ifdef THERMALLY_IDEAL_GAS
        cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(cPow2_ij)
        
        ! Diagonal matrix of eigenvalues
        l1 = abs(vel-c_ij)
        l2 = abs(vel)
        l3 = abs(vel+c_ij)
        l4 = abs(vel)
        
        ! Matrix of right eigenvectors
        R_ij(1,1) =  l1
        R_ij(2,1) =  l1*(u_ij-c_ij*c1)
        R_ij(3,1) =  l1*(v_ij-c_ij*c2)
        R_ij(4,1) =  l1*(H_ij-c_ij*vel)
        
        R_ij(1,2) =  l2
        R_ij(2,2) =  l2*u_ij
        R_ij(3,2) =  l2*v_ij
        R_ij(4,2) =  l2*q_ij
        
        R_ij(1,3) =  l3
        R_ij(2,3) =  l3*(u_ij+c_ij*c1)
        R_ij(3,3) =  l3*(v_ij+c_ij*c2)
        R_ij(4,3) =  l3*(H_ij+c_ij*vel)
        
        R_ij(1,4) =  0.0
        R_ij(2,4) =  l4*c2
        R_ij(3,4) = -l4*c1
        R_ij(4,4) =  l4*(u_ij*c2-v_ij*c1)
        
        ! Matrix of left eigenvectors
        L_ij(1,1) = 0.5*((GAMMA-1.0)*q_ij+c_ij*vel)/cPow2_ij
        L_ij(2,1) = (cPow2_ij-(GAMMA-1.0)*q_ij)/cPow2_ij
        L_ij(3,1) = 0.5*((GAMMA-1.0)*q_ij-c_ij*vel)/cPow2_ij
        L_ij(4,1) = v_ij*c1-u_ij*c2
        
        L_ij(1,2) = 0.5*(-(GAMMA-1.0)*u_ij-c_ij*c1)/cPow2_ij
        L_ij(2,2) = (GAMMA-1.0)*u_ij/cPow2_ij
        L_ij(3,2) = 0.5*(-(GAMMA-1.0)*u_ij+c_ij*c1)/cPow2_ij
        L_ij(4,2) = c2
        
        L_ij(1,3) = 0.5*(-(GAMMA-1.0)*v_ij-c_ij*c2)/cPow2_ij
        L_ij(2,3) = (GAMMA-1.0)*v_ij/cPow2_ij
        L_ij(3,3) = 0.5*(-(GAMMA-1.0)*v_ij+c_ij*c2)/cPow2_ij
        L_ij(4,3) = -c1
        
        L_ij(1,4) =  (GAMMA-1.0)/2.0/cPow2_ij
        L_ij(2,4) = -(GAMMA-1.0)/cPow2_ij
        L_ij(3,4) =  (GAMMA-1.0)/2.0/cPow2_ij
        L_ij(4,4) =  0.0
        
        ! Include scaling parameter
        anorm = dscale*anorm

        ! Compute tensorial dissipation D_ij = R_ij*|Lbd_ij|*L_ij
        do i = 1, NVAR2D
          do j = 1, NVAR2D
            aux = 0.0
            do k = 1, NVAR2D
              aux = aux + R_ij(i,k)*L_ij(k,j)
            end do
            DcoefficientsAtEdge(NVAR2D*(j-1)+i,1,idx) = anorm*aux
          end do
        end do
        
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0
        
      end if
    end do

  end subroutine hydro_calcMatRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDissMatD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 2D
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
    real(DP) :: ui,uj,vi,vj,ci,cj,Ei,Ej
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      Ei = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      Ej = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                           vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(3,2,idx) = dscale * (            uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                               (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(4,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                      vj*DmatrixCoeffsAtEdge(2,2,idx)))
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                           vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(3,3,idx) = dscale * (            ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                               (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(4,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                      vi*DmatrixCoeffsAtEdge(2,1,idx)))
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      DcoefficientsAtEdge(1,2,idx) = 0.0
      DcoefficientsAtEdge(2,2,idx) = -dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                            vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(3,2,idx) = -dscale * (            uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(4,2,idx) = -dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                       vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      DcoefficientsAtEdge(1,3,idx) = 0.0
      DcoefficientsAtEdge(2,3,idx) = -dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                            vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(3,3,idx) = -dscale * (            ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(4,3,idx) = -dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                       vi*DmatrixCoeffsAtEdge(2,2,idx)))
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*(ui*ui+vi*vi)), SYS_EPSREAL))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*(uj*uj+vj*vj)), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Compute dissipation tensor D_ij
      DcoefficientsAtEdge(:,1,idx) = dscale *&
          max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj+&
                   DmatrixCoeffsAtEdge(2,1,idx)*vj) +&
                   sqrt(DmatrixCoeffsAtEdge(1,1,idx)**2+&
                        DmatrixCoeffsAtEdge(2,1,idx)**2)*cj,&
               abs(DmatrixCoeffsAtEdge(1,2,idx)*ui+&
                   DmatrixCoeffsAtEdge(2,2,idx)*vi) +&
                   sqrt(DmatrixCoeffsAtEdge(1,2,idx)**2+&
                        DmatrixCoeffsAtEdge(2,2,idx)**2)*ci )
    end do

  end subroutine hydro_calcMatRusDissMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDiss2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 2D
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
    real(DP) :: ci,cj,Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj
    real(DP) :: uPow2i,uPow2j,vPow2i,vPow2j,aux1,aux2
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      Ei = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uvi = ui*vi; qi = ui*ui+vi*vi; uPow2i = ui*ui; vPow2i = vi*vi

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      Ej = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      uvj = uj*vj; qj = uj*uj+vj*vj; uPow2j = uj*uj; vPow2j = vj*vj

#ifdef HYDRO_USE_IBP
      aux1 = uj*DmatrixCoeffsAtEdge(1,2,idx)+vj*DmatrixCoeffsAtEdge(2,2,idx)
      aux2 = ui*DmatrixCoeffsAtEdge(1,1,idx)+vi*DmatrixCoeffsAtEdge(2,1,idx)

      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      DcoefficientsAtEdge( 1,2,idx) = 0.0
      DcoefficientsAtEdge( 2,2,idx) = dscale * (((GAMMA-1.0)/2.0*qj-uPow2j)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                                        uvj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 3,2,idx) = dscale * (((GAMMA-1.0)/2.0*qj-vPow2j)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                                        uvj*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge( 4,2,idx) = dscale * ((GAMMA-1.0)*qj-GAMMA*Ej)*aux1
      
      DcoefficientsAtEdge( 5,2,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge( 6,2,idx) = dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                            vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 7,2,idx) = dscale * (            vj*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                (GAMMA-1.0)*uj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 8,2,idx) = dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                (GAMMA-1.0)*uj*aux1)
      
      DcoefficientsAtEdge( 9,2,idx) = dscale * DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(10,2,idx) = dscale * (            uj*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                (GAMMA-1.0)*vj*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge(11,2,idx) = dscale * (            uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(12,2,idx) = dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                (GAMMA-1.0)*vj*aux1)
      
      DcoefficientsAtEdge(13,2,idx) = 0.0
      DcoefficientsAtEdge(14,2,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(15,2,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(16,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                       vj*DmatrixCoeffsAtEdge(2,2,idx)))
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      DcoefficientsAtEdge( 1,3,idx) = 0.0
      DcoefficientsAtEdge( 2,3,idx) = dscale * (((GAMMA-1.0)*qi-uPow2i)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                                    uvi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 3,3,idx) = dscale * (((GAMMA-1.0)*qi-vPow2i)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                                    uvi*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge( 4,3,idx) = dscale * ((GAMMA-1.0)*qi-GAMMA*Ei)*aux2
      
      DcoefficientsAtEdge( 5,3,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge( 6,3,idx) = dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                            vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 7,3,idx) = dscale * (            vi*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                (GAMMA-1.0)*ui*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 8,3,idx) = dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                (GAMMA-1.0)*ui*aux2)
      
      DcoefficientsAtEdge( 9,3,idx) = dscale * DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(10,3,idx) = dscale * (            ui*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                (GAMMA-1.0)*vi*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge(11,3,idx) = dscale * (            ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(12,3,idx) = dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                (GAMMA-1.0)*vi*aux2)
      
      DcoefficientsAtEdge(13,3,idx) = 0.0
      DcoefficientsAtEdge(14,3,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(15,3,idx) = dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(16,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                       vi*DmatrixCoeffsAtEdge(2,1,idx)))
#else
      aux1 = uj*DmatrixCoeffsAtEdge(1,1,idx)+vj*DmatrixCoeffsAtEdge(2,1,idx)
      aux2 = ui*DmatrixCoeffsAtEdge(1,2,idx)+vi*DmatrixCoeffsAtEdge(2,2,idx)

      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      DcoefficientsAtEdge( 1,2,idx) = 0.0
      DcoefficientsAtEdge( 2,2,idx) = -dscale * (((GAMMA-1.0)/2.0*qj-uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                                         uvj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 3,2,idx) = -dscale * (((GAMMA-1.0)/2.0*qj-vPow2j)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                                         uvj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge( 4,2,idx) = -dscale * ((GAMMA-1.0)*qj-GAMMA*Ej)*aux1
      
      DcoefficientsAtEdge( 5,2,idx) = -dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge( 6,2,idx) = -dscale * ((3.0-GAMMA)*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                             vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 7,2,idx) = -dscale * (            vj*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                 (GAMMA-1.0)*uj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 8,2,idx) = -dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                 (GAMMA-1.0)*uj*aux1)
      
      DcoefficientsAtEdge( 9,2,idx) = -dscale * DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(10,2,idx) = -dscale * (            uj*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                 (GAMMA-1.0)*vj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge(11,2,idx) = -dscale * (            uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                 (3.0-GAMMA)*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(12,2,idx) = -dscale * ((GAMMA*Ej-(GAMMA-1.0)/2.0*qj)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                 (GAMMA-1.0)*vj*aux1)
      
      DcoefficientsAtEdge(13,2,idx) = 0.0
      DcoefficientsAtEdge(14,2,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(15,2,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(16,2,idx) = -dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                        vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      DcoefficientsAtEdge( 1,3,idx) = 0.0
      DcoefficientsAtEdge( 2,3,idx) = -dscale * (((GAMMA-1.0)*qi-uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                                     uvi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 3,3,idx) = -dscale * (((GAMMA-1.0)*qi-vPow2i)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                                     uvi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge( 4,3,idx) = -dscale * ((GAMMA-1.0)*qi-GAMMA*Ei)*aux2
      
      DcoefficientsAtEdge( 5,3,idx) = -dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge( 6,3,idx) = -dscale * ((3.0-GAMMA)*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                             vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 7,3,idx) = -dscale * (            vi*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                 (GAMMA-1.0)*ui*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 8,3,idx) = -dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                 (GAMMA-1.0)*ui*aux2)
      
      DcoefficientsAtEdge( 9,3,idx) = -dscale * DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(10,3,idx) = -dscale * (            ui*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                 (GAMMA-1.0)*vi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge(11,3,idx) = -dscale * (            ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                 (3.0-GAMMA)*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(12,3,idx) = -dscale * ((GAMMA*Ei-(GAMMA-1.0)/2.0*qi)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                 (GAMMA-1.0)*vi*aux2)
      
      DcoefficientsAtEdge(13,3,idx) = 0.0
      DcoefficientsAtEdge(14,3,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(15,3,idx) = -dscale * (GAMMA-1.0)*DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(16,3,idx) = -dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                        vi*DmatrixCoeffsAtEdge(2,2,idx)))
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*(ui*ui+vi*vi)), SYS_EPSREAL))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*(uj*uj+vj*vj)), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Compute dissipation tensor D_ij
      aux1 = dscale *&
          max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj+&
                   DmatrixCoeffsAtEdge(2,1,idx)*vj) +&
                   sqrt(DmatrixCoeffsAtEdge(1,1,idx)**2+&
                        DmatrixCoeffsAtEdge(2,1,idx)**2)*cj,&
               abs(DmatrixCoeffsAtEdge(1,2,idx)*ui+&
                   DmatrixCoeffsAtEdge(2,2,idx)*vi) +&
                   sqrt(DmatrixCoeffsAtEdge(1,2,idx)**2+&
                        DmatrixCoeffsAtEdge(2,2,idx)**2)*ci )

      DcoefficientsAtEdge( :,1,idx) = 0.0
      DcoefficientsAtEdge( 1,1,idx) = aux1
      DcoefficientsAtEdge( 6,1,idx) = aux1
      DcoefficientsAtEdge(11,1,idx) = aux1
      DcoefficientsAtEdge(16,1,idx) = aux1
    end do

  end subroutine hydro_calcMatRusDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcCharacteristics2d_sim(Dweight, DdataAtEdge,&
      DcharVariablesAtEdge, DeigenvaluesAtEdge,&
      DrightEigenvectorsAtEdge, DleftEigenvectorsAtEdge, rcollection)

!<description>
    ! This subroutine computes the characteristic variables in 2D
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
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: u_ij,v_ij,H_ij,q_ij,c_ij,aux,aux1,aux2,ui,uj,vi,vj,cPow2_ij,a1,a2,anorm
    integer :: idx

    ! Compute norm of weighting coefficient
    anorm = sqrt(Dweight(1)*Dweight(1)+Dweight(2)*Dweight(2))

    ! Check if weighting coefficient is zero
    if (anorm .le. SYS_EPSREAL) then
      if (present(DcharVariablesAtEdge))     DcharVariablesAtEdge     = 0.0
      if (present(DeigenvaluesAtEdge))       DeigenvaluesAtEdge       = 0.0
      if (present(DrightEigenvectorsAtEdge)) DrightEigenvectorsAtEdge = 0.0
      if (present(DleftEigenvectorsAtEdge))  DleftEigenvectorsAtEdge  = 0.0

      ! That`s it
      return
    end if

    ! Compute normalised weighting coefficient
    a1  = Dweight(1)/anorm
    a2  = Dweight(2)/anorm

    ! Do we have to compute characteristic variables
    if (present(DcharVariablesAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)
        
        ! Compute velocities
        ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
        vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
        uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
        vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx), aux)
        
        ! Compute auxiliary variables
        q_ij  = 0.5*(u_ij*u_ij+v_ij*v_ij)
#ifdef THERMALLY_IDEAL_GAS
        cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(cPow2_ij)
        aux  = a1*u_ij+a2*v_ij

        ! Compute solution difference U_j-U_i
        Diff = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = (GAMMA-1.0)/2.0/cPow2_ij*(q_ij*Diff(1)-u_ij*Diff(2)-v_ij*Diff(3)+Diff(4) )
        aux2 = 0.5*(aux*Diff(1)-a1*Diff(2)-a2*Diff(3) )/c_ij

        ! Compute characteristic variables
        DcharVariablesAtEdge(1,idx) = anorm * (aux1 + aux2)
        DcharVariablesAtEdge(2,idx) = anorm * ((1.0-(GAMMA-1.0)*q_ij/cPow2_ij)*Diff(1)+&
                                                             (GAMMA-1.0)*(u_ij*Diff(2)+&
                                                                          v_ij*Diff(3)-&
                                                                      Diff(4))/cPow2_ij)
        DcharVariablesAtEdge(3,idx) = anorm * (aux1 - aux2)
        DcharVariablesAtEdge(4,idx) = anorm * ((a1*v_ij-a2*u_ij)*Diff(1)+&
                                                              a2*Diff(2)-&
                                                              a1*Diff(3))
      end do
    end if

    ! Do we have to compute eigenvalues
    if (present(DeigenvaluesAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)

        ! Compute velocities
        ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
        vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
        uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
        vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx), aux)
        
        ! Compute auxiliary variables
        q_ij  = 0.5*(u_ij*u_ij+v_ij*v_ij)
#ifdef THERMALLY_IDEAL_GAS
        c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
        aux  = a1*u_ij+a2*v_ij

        ! Compute eigenvalues
        DeigenvaluesAtEdge(1,idx) = aux-c_ij
        DeigenvaluesAtEdge(2,idx) = aux
        DeigenvaluesAtEdge(3,idx) = aux+c_ij
        DeigenvaluesAtEdge(4,idx) = aux
      end do
    end if

    ! Do we have to compute right eigenvectors
    if (present(DrightEigenvectorsAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)

        ! Compute velocities
        ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
        vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
        uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
        vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx), aux)
        
        ! Compute auxiliary variables
        q_ij  = 0.5*(u_ij*u_ij+v_ij*v_ij)
#ifdef THERMALLY_IDEAL_GAS
        c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
        aux  = a1*u_ij+a2*v_ij

        ! Compute right eigenvectors
        DrightEigenvectorsAtEdge( 1,idx) =  1.0
        DrightEigenvectorsAtEdge( 2,idx) =  u_ij-c_ij*a1
        DrightEigenvectorsAtEdge( 3,idx) =  v_ij-c_ij*a2
        DrightEigenvectorsAtEdge( 4,idx) =  H_ij-c_ij*aux

        DrightEigenvectorsAtEdge( 5,idx) =  1.0
        DrightEigenvectorsAtEdge( 6,idx) =  u_ij
        DrightEigenvectorsAtEdge( 7,idx) =  v_ij
        DrightEigenvectorsAtEdge( 8,idx) =  q_ij

        DrightEigenvectorsAtEdge( 9,idx) =  1.0
        DrightEigenvectorsAtEdge(10,idx) =  u_ij+c_ij*a1
        DrightEigenvectorsAtEdge(11,idx) =  v_ij+c_ij*a2
        DrightEigenvectorsAtEdge(12,idx) =  H_ij+c_ij*aux

        DrightEigenvectorsAtEdge(13,idx) =  0.0
        DrightEigenvectorsAtEdge(14,idx) =  a2
        DrightEigenvectorsAtEdge(15,idx) = -a1
        DrightEigenvectorsAtEdge(16,idx) =  u_ij*a2-v_ij*a1
      end do
    end if

    ! Do we have to compute left eigenvectors
    if (present(DleftEigenvectorsAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)

        ! Compute velocities
        ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
        vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
        uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
        vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx), aux)
        
        ! Compute auxiliary variables
        q_ij  = 0.5*(u_ij*u_ij+v_ij*v_ij)
#ifdef THERMALLY_IDEAL_GAS
        cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(cPow2_ij)
        aux  = a1*u_ij+a2*v_ij

        ! Compute left eigenvectors
        DleftEigenvectorsAtEdge( 1,idx) =  0.5*((GAMMA-1.0)*q_ij+c_ij*aux)/cPow2_ij
        DleftEigenvectorsAtEdge( 2,idx) = (cPow2_ij-(GAMMA-1.0)*q_ij)/cPow2_ij
        DleftEigenvectorsAtEdge( 3,idx) =  0.5*((GAMMA-1.0)*q_ij-c_ij*aux)/cPow2_ij
        DleftEigenvectorsAtEdge( 4,idx) =  v_ij*a1-u_ij*a2

        DleftEigenvectorsAtEdge( 5,idx) =  0.5*(-(GAMMA-1.0)*u_ij-c_ij*a1)/cPow2_ij
        DleftEigenvectorsAtEdge( 6,idx) =  (GAMMA-1.0)*u_ij/cPow2_ij
        DleftEigenvectorsAtEdge( 7,idx) =  0.5*(-(GAMMA-1.0)*u_ij+c_ij*a1)/cPow2_ij
        DleftEigenvectorsAtEdge( 8,idx) =  a2

        DleftEigenvectorsAtEdge( 9,idx) =  0.5*(-(GAMMA-1.0)*v_ij-c_ij*a2)/cPow2_ij
        DleftEigenvectorsAtEdge(10,idx) =  (GAMMA-1.0)*v_ij/cPow2_ij
        DleftEigenvectorsAtEdge(11,idx) =  0.5*(-(GAMMA-1.0)*v_ij+c_ij*a2)/cPow2_ij
        DleftEigenvectorsAtEdge(12,idx) = -a1

        DleftEigenvectorsAtEdge(13,idx) =  (GAMMA-1.0)/2.0/cPow2_ij
        DleftEigenvectorsAtEdge(14,idx) = -(GAMMA-1.0)/cPow2_ij
        DleftEigenvectorsAtEdge(15,idx) =  (GAMMA-1.0)/2.0/cPow2_ij
        DleftEigenvectorsAtEdge(16,idx) =  0.0
      end do
    end if

  end subroutine hydro_calcCharacteristics2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTScalarDiss2d(&
      U1_i, U1_j, U2_i, U2_j, Coeff_ij, Coeff_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 2D using scalar dissipation.
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
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: ui,vi,uj,vj
    real(DP) :: d_ij,H_ij,q_ij,u_ij,v_ij,aux,vel_ij,c_ij,anorm

    ! Compute velocities
    ui = X_VELOCITY_FROM_CONSVAR(U2_i,NVAR2D)
    vi = Y_VELOCITY_FROM_CONSVAR(U2_i,NVAR2D)

    uj = X_VELOCITY_FROM_CONSVAR(U2_j,NVAR2D)
    vj = Y_VELOCITY_FROM_CONSVAR(U2_j,NVAR2D)

    ! Compute skew-symmetric coefficient
    a = 0.5*(Coeff_ij-Coeff_ji); anorm = sqrt(a(1)*a(1)+a(2)*a(2))

    ! Compute Roe mean values
    aux  = ROE_MEAN_RATIO(MYNEWLINE
           DENSITY_FROM_CONSVAR(U2_i,NVAR2D), DENSITY_FROM_CONSVAR(U2_j,NVAR2D))
    u_ij = ROE_MEAN_VALUE(ui,uj,aux)
    v_ij = ROE_MEAN_VALUE(vi,vj,aux)
    H_ij = ROE_MEAN_VALUE(MYNEWLINE
           (TOTAL_ENERGY_FROM_CONSVAR(U2_i,NVAR2D)+MYNEWLINE
            PRESSURE_FROM_CONSVAR_2D(U2_i,NVAR2D))/DENSITY_FROM_CONSVAR(U2_i,NVAR2D),MYNEWLINE
           (TOTAL_ENERGY_FROM_CONSVAR(U2_j,NVAR2D)+MYNEWLINE
            PRESSURE_FROM_CONSVAR_2D(U2_j,NVAR2D))/DENSITY_FROM_CONSVAR(U2_j,NVAR2D), aux)
            
    ! Compute auxiliary variables
    vel_ij = u_ij*a(1) + v_ij*a(2)
    q_ij   = 0.5*(u_ij*u_ij+v_ij*v_ij)
#ifdef THERMALLY_IDEAL_GAS
    c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif

    ! Scalar dissipation
    d_ij = abs(vel_ij) + anorm*c_ij

    ! Compute conservative fluxes
    F_ij = dscale1*(U1_i-U1_j) + dscale2*d_ij*(U2_i-U2_j)

  end subroutine hydro_calcFluxFCTScalarDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRoeDiss2d(&
      U1_i, U1_j, U2_i, U2_j, Coeff_ij, Coeff_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 2D using tensorial dissipation.
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
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: ui,vi,uj,vj
    real(DP) :: aux,aux1,aux2,uPow2,vPow2,H_ij,q_ij,u_ij,v_ij
    real(DP) :: anorm,l1,l2,l3,l4,w1,w2,w3,w4,cPow2_ij,c_ij,vel_ij

    ! Compute velocities
    ui = X_VELOCITY_FROM_CONSVAR(U2_i,NVAR2D)
    vi = Y_VELOCITY_FROM_CONSVAR(U2_i,NVAR2D)

    uj = X_VELOCITY_FROM_CONSVAR(U2_j,NVAR2D)
    vj = Y_VELOCITY_FROM_CONSVAR(U2_j,NVAR2D)
    
    ! Compute the skew-symmetric coefficient
    a = 0.5*(Coeff_ij-Coeff_ji); anorm = sqrt(a(1)*a(1)+a(2)*a(2))

    if (anorm .gt. SYS_EPSREAL) then

      ! Normalize the skew-symmetric coefficient
      a = a/anorm

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(MYNEWLINE
             DENSITY_FROM_CONSVAR(U2_i,NVAR2D), DENSITY_FROM_CONSVAR(U2_j,NVAR2D))
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      H_ij = ROE_MEAN_VALUE(MYNEWLINE
             (TOTAL_ENERGY_FROM_CONSVAR(U2_i,NVAR2D)+MYNEWLINE
              PRESSURE_FROM_CONSVAR_2D(U2_i,NVAR2D))/DENSITY_FROM_CONSVAR(U2_i,NVAR2D),MYNEWLINE
             (TOTAL_ENERGY_FROM_CONSVAR(U2_j,NVAR2D)+MYNEWLINE
              PRESSURE_FROM_CONSVAR_2D(U2_j,NVAR2D))/DENSITY_FROM_CONSVAR(U2_j,NVAR2D), aux)

      ! Compute auxiliary variables
      vel_ij = u_ij*a(1) + v_ij*a(2)
      uPow2  = u_ij*u_ij
      vPow2  = v_ij*v_ij
      q_ij   = 0.5*(uPow2+vPow2)
#ifdef THERMALLY_IDEAL_GAS
      cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
      c_ij  = sqrt(cPow2_ij)

      ! Compute eigenvalues
      l1 = abs(vel_ij-c_ij)
      l2 = abs(vel_ij)
      l3 = abs(vel_ij+c_ij)
      l4 = abs(vel_ij)

      ! Compute solution difference U2_i-U2_j
      Diff = U2_i-U2_j

      ! Compute auxiliary quantities for characteristic variables
      aux1 = (GAMMA-1.0)/2.0/cPow2_ij*(q_ij*Diff(1)-u_ij*Diff(2)-v_ij*Diff(3)+Diff(4))
      aux2 = 0.5*(vel_ij*Diff(1)-a(1)*Diff(2)-a(2)*Diff(3))/c_ij

      ! Compute characteristic variables multiplied by the corresponding eigenvalue
      w1 = l1 * (aux1 + aux2)
      w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/cPow2_ij)*Diff(1)+(GAMMA-1.0)*(u_ij*Diff(2)+v_ij*Diff(3)-Diff(4))/cPow2_ij)
      w3 = l3 * (aux1 - aux2)
      w4 = l4 * ((a(1)*v_ij-a(2)*u_ij)*Diff(1)+a(2)*Diff(2)-a(1)*Diff(3))

      ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
      Diff(1) = w1 + w2 + w3
      Diff(2) = (u_ij-c_ij*a(1))*w1 + u_ij*w2 + (u_ij+c_ij*a(1))*w3 + a(2)*w4
      Diff(3) = (v_ij-c_ij*a(2))*w1 + v_ij*w2 + (v_ij+c_ij*a(2))*w3 - a(1)*w4
      Diff(4) = (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3 + (u_ij*a(2)-v_ij*a(1))*w4

      ! Compute conservative flux
      F_ij = dscale1*(U1_i-U1_j) + dscale2*anorm*Diff

    else

      ! Compute conservative flux without spatial contribution
      F_ij = dscale1*(U1_i-U1_j)

    end if
  end subroutine hydro_calcFluxFCTRoeDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRusanov2d(&
      U1_i, U1_j, U2_i, U2_j, Coeff_ij, Coeff_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 2D using the Rusanov dissipation.
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
    real(DP) :: ui,vi,uj,vj
    real(DP) :: d_ij,ci,cj,Ei,Ej

    ! Compute velocities and energy
    ui = X_VELOCITY_FROM_CONSVAR(U2_i,NVAR2D)
    vi = Y_VELOCITY_FROM_CONSVAR(U2_i,NVAR2D)
    Ei = TOTAL_ENERGY_FROM_CONSVAR(U2_i,NVAR2D)

    uj = X_VELOCITY_FROM_CONSVAR(U2_j,NVAR2D)
    vj = Y_VELOCITY_FROM_CONSVAR(U2_j,NVAR2D)
    Ej = TOTAL_ENERGY_FROM_CONSVAR(U2_j,NVAR2D)

    ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
    ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*(ui*ui+vi*vi)), SYS_EPSREAL))
    cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*(uj*uj+vj*vj)), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif

    ! Scalar dissipation for the Rusanov flux
    d_ij = max( abs(Coeff_ij(1)*uj+Coeff_ij(2)*vj) +&
                sqrt(Coeff_ij(1)*Coeff_ij(1)+Coeff_ij(2)*Coeff_ij(2))*cj,&
                abs(Coeff_ji(1)*ui+Coeff_ji(2)*vi) +&
                sqrt(Coeff_ji(1)*Coeff_ji(1)+Coeff_ji(2)*Coeff_ji(2))*ci )

    ! Compute conservative fluxes
    F_ij = dscale1*(U1_i-U1_j) + dscale2*d_ij*(U2_i-U2_j)

  end subroutine hydro_calcFluxFCTRusanov2d

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDensity2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density in 2D
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
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
    end do

  end subroutine hydro_trafoFluxDensity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDensity2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density in 2D
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
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
    end do

  end subroutine hydro_trafoDiffDensity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxEnergy2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the energy in 2D
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
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
    end do

  end subroutine hydro_trafoFluxEnergy2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffEnergy2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the energy in 2D
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
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
    end do

  end subroutine hydro_trafoDiffEnergy2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxPressure2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the pressure in 2D
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
    real(DP) :: ui,uj,vi,vj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      
      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(1,1,idx) = (GAMMA-1.0)*&
          (0.5*(ui*ui+vi*vi)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vi*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))
      DtransformedFluxesAtEdge(1,2,idx) =-(GAMMA-1.0)*&
          (0.5*(uj*uj+vj*vj)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vj*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do

  end subroutine hydro_trafoFluxPressure2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffPressure2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the pressure in 2D
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
          PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx)
    end do

  end subroutine hydro_trafoDiffPressure2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxVelocity2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the x- and y-velocity
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
    real(DP) :: ui,uj,vi,vj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(1,1,idx) =&
          (X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -(X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Transformed velocity fluxes in y-direction
      DtransformedFluxesAtEdge(2,1,idx) =&
          (Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -(Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
    end do
    
  end subroutine hydro_trafoFluxVelocity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffVelocity2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the x- and y-velocity
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
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)

      ! Transformed velocity difference in y-direction
      DtransformedDataAtEdge(2,idx) =&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
    end do
    
  end subroutine hydro_trafoDiffVelocity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxMomentum2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the x- and y-momentum
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
          X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)

      ! Transformed momentum fluxes in y-direction
      DtransformedFluxesAtEdge(2,1,idx) =&
          Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
    end do
    
  end subroutine hydro_trafoFluxMomentum2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffMomentum2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the x- and y-momentum
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
          X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)

      ! Transformed momentum difference in y-direction
      DtransformedDataAtEdge(2,idx) =&
          Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
    end do
    
  end subroutine hydro_trafoDiffMomentum2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenEng2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 2D
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
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)

      ! Transformed total energy fluxes
      DtransformedFluxesAtEdge(2,1,idx) =&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
    end do

  end subroutine hydro_trafoFluxDenEng2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenEng2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 2D
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
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)

      ! Transformed total energy difference
      DtransformedDataAtEdge(2,idx) =&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
    end do

  end subroutine hydro_trafoDiffDenEng2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenPre2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 2D
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
    real(DP) :: ui,uj,vi,vj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
      
      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(2,1,idx) = (GAMMA-1.0)*&
          (0.5*(ui*ui+vi*vi)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vi*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))
      DtransformedFluxesAtEdge(2,2,idx) =-(GAMMA-1.0)*&
          (0.5*(uj*uj+vj*vj)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vj*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do

  end subroutine hydro_trafoFluxDenPre2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenPre2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 2D
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
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(2,idx) =&
          PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx)
    end do
    
  end subroutine hydro_trafoDiffDenPre2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenPreVel2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation
    ! of the given flux into primitive variables in 2D
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
    real(DP) :: ui,uj,vi,vj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)

      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(2,1,idx) =&
          (X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -(X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Transformed velocity fluxes in y-direction
      DtransformedFluxesAtEdge(3,1,idx) =&
          (Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      DtransformedFluxesAtEdge(3,2,idx) =&
         -(Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(4,1,idx) =(GAMMA-1.0)*&
          (0.5*(ui*ui+vi*vi)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vi*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))
      DtransformedFluxesAtEdge(4,2,idx) =-(GAMMA-1.0)*&
          (0.5*(uj*uj+vj*vj)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vj*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do

  end subroutine hydro_trafoFluxDenPreVel2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenPreVel2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density, pressure and velocity in 2D
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
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)

      ! Transformed velocity difference in x-direction
      DtransformedDataAtEdge(2,idx) =&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      
      ! Transformed velocity difference in y-direction
      DtransformedDataAtEdge(3,idx) =&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)

      ! Transformed pressure difference
      DtransformedDataAtEdge(4,idx) =&
          PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR_2D(DdataAtEdge,NVAR2D,1,idx)
    end do

  end subroutine hydro_trafoDiffDenPreVel2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcBoundaryvalues2d(DbdrNormal, DpointNormal,&
      DbdrValue, ibdrCondType, Du, Du0, istatus)

!<description>
    ! This subroutine computes the boundary values for a given node in 2D
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

    ! local variables
    real(DP), dimension(NVAR2D) :: W,Wu,Winf    ! Riemann invariants, eigenvalues, etc.
    real(DP) :: rho,v1,v2,p,E,c,v1_0,v2_0       ! primitive variables
    real(DP) :: v1_b,v2_b,vn_b,vn,vt,pstar,ps   ! velocities and boundary values
    real(DP) :: cup,f,fd,ge,qrt                 ! auxiliary variables ...
    real(DP) :: pold,ppv,prat,ptl,ptr,vdiff,vm  ! ... for the Riemann solver
    real(DP) :: auxA,auxB,aux,dnx2,dny2,dnxy
    integer:: ite

    ! What type of boundary condition is given?
    select case(iand(ibdrCondType, BDRC_TYPEMASK))
    case(BDRC_FREESLIP, BDRC_RLXFREESLIP)
      !-------------------------------------------------------------------------

      ! The wall boundary conditions follow algorithm II from the paper
      !
      !    `High-order accurate implementation of solid wall
      !     boundary conditions in curved geometries`
      !     L. Krivodonova and M. Berger, J. Comput. Physics 211, (2006) 492-512
      !
      ! From the computed primitive values U=[rho, v1, v2, p] the boundary
      ! values are determined as follows:
      !
      !     $$\rho_b = \rho$
      !     $$v1_b   = v_1*(n_y^2-n_x^2)-2*n_x*n_y*v_2$$
      !     $$v2_b   = v_2*(n_x^2-n_y^2)-2*n_x*n_y*v_1$$
      !     $$p_b    = p$
      !
      ! where $n=[n_x,n_y]$ denotes the physical normal vector which is given
      ! analytically, i.e. it is more accurate than the finite element normal.
      ! The Riemann problem Riem(U, U_b, N) in the direction of the numerical
      ! normal vector $N=[N_x,N_y]$ is solved exactly. Due to the identical
      ! pressure and density in the two states, the exact solution consists if
      ! either two shocks or two rarefaction waves.
      !
      ! Note that the relaxed wall boundary conditions is intended to prevent
      ! impulsive start, that is, the fluid is allowed to seep through the wall
      ! at startup and the normal velocity is gradually driven to zero as the
      ! flow evolves. This technique is presented and analyzed by Lyra:
      !
      !    `Unstructured Grid Adaptive Algorithms for
      !     Fluid Dynamics and Heat Conduction`
      !     P.R.M. Lyra, PhD thesis, University of Wales, Swansea, 1994.
      !
      ! In the framework of algorithm II by Krivodonova and Berger the boundary
      ! values are determined as follows:
      !
      ! rho_b = rho
      ! v1_b  = v_1*(ny^2-nx^2)-2*nx*ny*v_2+2*c*(v_1^n*n_x^2+v_2^n*n_x*n_y)
      ! v2_b  = v_2*(nx^2-ny^2)-2*nx*ny*v_1+2*c*(v_2^n*n_y^2+v_1^n*n_x*n_y)
      ! p_b   = p

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      v2  = Du(3)/rho
      E   = Du(4)/rho
      p   = (GAMMA-1.0)*rho*(E-0.5*(v1*v1+v2*v2))
      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))

      ! Precompute auxiliary data
      dnxy = DbdrNormal(1)*DbdrNormal(2)
      dnx2 = DbdrNormal(1)*DbdrNormal(1)
      dny2 = DbdrNormal(2)*DbdrNormal(2)

      if (ibdrCondType .eq. BDRC_FREESLIP) then
        ! Compute reflected velocities at the boundary
        v1_b = (dny2-dnx2)*v1 - 2.0*dnxy*v2
        v2_b = (dnx2-dny2)*v2 - 2.0*dnxy*v1
      else
        ! Compute initial velocity from previous time step
        v1_0 = Du0(2)/Du0(1)
        v2_0 = Du0(3)/Du0(1)

        ! Compute semi-reflected velocities at the boundary
        v1_b = (dny2-dnx2)*v1-2.0*dnxy*v2 + 2*DbdrValue(1)*(v1_0*dnx2+v2_0*dnxy)
        v2_b = (dnx2-dny2)*v2-2.0*dnxy*v1 + 2*DbdrValue(1)*(v2_0*dny2+v1_0*dnxy)
      end if

      ! Compute normal velocities at the boundary and the ghost state
      ! w.r.t. the numerical/approximate  outward unit normal vector
      vn   = DpointNormal(1)*v1   + DpointNormal(2)*v2
      vn_b = DpointNormal(1)*v1_b + DpointNormal(2)*v2_b

      ! Compute the tangential velocity depending on the sign of N*v
      if (vn .gt. 0.0) then
        vt = DpointNormal(2)*v1   - DpointNormal(1)*v2
      else
        vt = DpointNormal(2)*v1_b - DpointNormal(1)*v2_b
      end if


      !-------------------------------------------------------------------------
      ! Calculate the pressure in the star region
      !
      ! Note that the pressure equation can only be solved if the pressure
      ! positivity condition is satisfied, that is
      !
      !     $$\frac{2}{\gamma-1}(c+c_b)>v_b-v$$
      !
      ! Otherwise, the Riemann problem gives rise to vacuum so that the
      ! "star region" does no longer exist and the standard procedure fails.
      !
      ! Here and below, the left state corresponds to the interior value
      ! and the right state corresponds to the ghost values since the unit
      ! normal vector is directed outward to the boundary.
      !-------------------------------------------------------------------------

      ! Check the pressure positivity condition
      if (2.0*(2.0/(GAMMA-1.0))*c .le. vn_b-vn) then
        if (present(istatus)) then
          istatus = -ibdrCondType
          return
        else
          call output_line('Riemann solver failed due to vacuum',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcBoundaryvalues2d')
          call sys_halt()
        end if
      end if

      ! Provide a guess value for pressure in the "star region"
      ! by using the PVRS Riemann solver as suggested by Toro

      cup  = rho*c
      ppv  = p+0.5*(vn-vn_b)*cup
      ppv  = max(0.0, ppv)

      if (ppv .eq. p) then

        ! Select guessed pressure from PVRS Riemann solver
        pstar = ppv
      else
        if (ppv .lt. p) then

          ! Guess pressure from the Two-Rarefaction Riemann solver
          vm    = 0.5*(vn+vn_b)
          ptl   = 1.0 + (GAMMA-1.0)/2.0*(vn-vm)/c
          ptr   = 1.0 + (GAMMA-1.0)/2.0*(vm-vn_b)/c
          pstar = 0.5*(p*ptl + p*ptr)**(2.0*GAMMA/(GAMMA-1.0))
        else

          ! Guess pressure from the Two-Shock Riemann solver
          ! with PVRS as estimated pressure value
          ge    = sqrt((2.0/(GAMMA+1.0)/rho)/((GAMMA-1.0)/(GAMMA+1.0)*p+ppv))
          pstar = p - 0.5*(vn_b-vn)/ge
        end if
      end if

      ! Initialize solution difference and pressure
      vdiff = (vn_b-vn)/2.0
      pold  = pstar

      newton: do ite = 1, 100

        ! Compute pressure function f(pold) and its derivative f1(pold)
        if (pold .le. p) then

          ! Rarefaction wave
          prat = pold/p

          f  = (2.0/(GAMMA-1.0))*c*(prat**((GAMMA-1.0)/(2.0*GAMMA)) - 1.0)
          fd = (1.0/(rho*c))*prat**(-(GAMMA+1.0)/(2.0*GAMMA))
        else

          ! Shock wave
          auxA = 2.0/(GAMMA+1.0)/rho
          auxB = (GAMMA-1.0)/(GAMMA+1.0)*p
          qrt  = sqrt(auxA/(auxB + pold))

          f  = (pold-p)*qrt
          fd = (1.0 - 0.5*(pold - p)/(auxB + pold))*qrt
        end if

        pstar = pold - (f+vdiff)/fd
        if (pstar .lt. 0.0) then
          pold = 1.0E-6
          cycle newton
        end if

        aux = 2.0*abs((pstar-pold)/(pstar+pold))
        if (aux .le. 1.0E-6)  exit newton

        pold = pstar

      end do newton

      ! Check if Newton`s method converged
      if (ite .ge. 100) then
        if (present(istatus)) then
          istatus = -ibdrCondType
          return
        else
          call output_line('Riemann solver failed due to divergence in' // &
              ' Newton-Raphson iteration',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcBoundaryvalues2d')
          call sys_halt()
        end if
      end if

      !-------------------------------------------------------------------------
      ! Calculate the velocity in the star region
      !-------------------------------------------------------------------------

      ! Note that the contribution fR-fL vanishes due to constant states
      vn = 0.5*(vn+vn_b)


      !-------------------------------------------------------------------------
      ! Calculate the density in the star region
      !-------------------------------------------------------------------------

      if (pstar .le. p) then

        ! Rarefaction wave
        rho = rho*(pstar/p)**(1.0/GAMMA)
      else

        ! Shock wave
        rho = rho*(pstar/p+(GAMMA-1.0)/(GAMMA+1.0))/((GAMMA-1.0)/(GAMMA+1.0)*(pstar/p)+1.0)
      end if

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DpointNormal(2)*vt+DpointNormal(1)*vn)
      Du(3) = rho*(-DpointNormal(1)*vt+DpointNormal(2)*vn)
      Du(4) = pstar/(GAMMA-1.0)+0.5*rho*(vn*vn+vt*vt)


    case(BDRC_VISCOUSWALL)
      !-------------------------------------------------------------------------

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      v2  = Du(3)/rho
      E   = Du(4)/rho
      p   = (GAMMA-1.0)*rho*(E-0.5*(v1*v1+v2*v2))

      ! Update the solution vector and let vn:=0 and vt:=0
      Du(2) = 0.0
      Du(3) = 0.0
      Du(4) = p/(GAMMA-1.0)


    case(BDRC_SUPERINLET)
      !-------------------------------------------------------------------------

      ! The free stream primitive variables are Deval=[rho,v1,v2,p]
      rho = DbdrValue(1)
      p   = DbdrValue(4)
      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*DbdrValue(2)+DbdrNormal(2)*DbdrValue(3)
      vt  = DbdrNormal(2)*DbdrValue(2)-DbdrNormal(1)*DbdrValue(3)

      ! Compute Riemann invariants based on the free stream values
      W(1) = vn-2*c/(GAMMA-1.0)
      W(2) = vn+2*c/(GAMMA-1.0)
      W(3) = p/(rho**GAMMA)
      W(4) = vt

      ! Transform back into conservative variables
      vn   = 0.5*(W(1)+W(2))
      c    = 0.25*(GAMMA-1.0)*(W(2)-W(1))
      rho  = (c*c/GAMMA/W(3))**(1.0/(GAMMA-1.0))
      p    = rho*c*c/GAMMA
      vt   = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/(GAMMA-1.0)+0.5*rho*(vn*vn+vt*vt)


    case(BDRC_FREESTREAM)
      !-------------------------------------------------------------------------

      ! The free stream primitive variables are Deval=[rho,v1,v2,p]
      rho = DbdrValue(1)
      p   = DbdrValue(4)
      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*DbdrValue(2)+DbdrNormal(2)*DbdrValue(3)
      vt  = DbdrNormal(2)*DbdrValue(2)-DbdrNormal(1)*DbdrValue(3)

      ! Compute Riemann invariants based on the free stream values
      Winf(1) = vn-2.0*c/(GAMMA-1.0)
      Winf(2) = vn+2.0*c/(GAMMA-1.0)
      Winf(3) = p/(rho**GAMMA)
      Winf(4) = vt

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      v2  = Du(3)/rho
      E   = Du(4)/rho
      p   = (GAMMA-1.0)*rho*(E-0.5*(v1*v1+v2*v2))

      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*v1+DbdrNormal(2)*v2
      vt  = DbdrNormal(2)*v1-DbdrNormal(1)*v2

      ! Compute Riemann invariants based on the solution values
      Wu(1) = vn-2.0*c/(GAMMA-1.0)
      Wu(2) = vn+2.0*c/(GAMMA-1.0)
      Wu(3) = p/(rho**GAMMA)
      Wu(4) = vt

      ! Adopt free stream/computed values depending on the sign of the eigenvalue
      W(1) = merge(Winf(1), Wu(1), vn <  c)
      W(2) = merge(Winf(2), Wu(2), vn < -c)
      W(3) = merge(Winf(3), Wu(3), vn <  SYS_EPSREAL)
      W(4) = merge(Winf(4), Wu(4), vn <  SYS_EPSREAL)

      ! Transform back into conservative variables
      vn   = 0.5*(W(1)+W(2))
      c    = 0.25*(GAMMA-1.0)*(W(2)-W(1))
      rho  = (c*c/GAMMA/W(3))**(1.0/(GAMMA-1.0))
      p    = rho*c*c/GAMMA
      vt   = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/(GAMMA-1.0)+0.5*rho*(vn*vn+vt*vt)


    case(BDRC_SUBINLET)
      !-------------------------------------------------------------------------

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      v2  = Du(3)/rho
      E   = Du(4)/rho
      p   = (GAMMA-1.0)*rho*(E-0.5*(v1*v1+v2*v2))

      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*v1+DbdrNormal(2)*v2
      vt  = DbdrNormal(2)*v1-DbdrNormal(1)*v2

      ! The specified density and pressure is Deval=[rho,p]
      rho = DbdrValue(1)
      p   = DbdrValue(2)

      ! Compute Riemann invariants
      W(1) = vn-2.0*c/(GAMMA-1.0)
      W(2) = vn+2.0*c/(GAMMA-1.0)
      W(3) = p/(rho**GAMMA)
      W(4) = vt

      ! Transform back into conservative variables
      vn   = 0.5*(W(1)+W(2))
      c    = 0.25*(GAMMA-1.0)*(W(2)-W(1))
      rho  = (c*c/GAMMA/W(3))**(1.0/(GAMMA-1.0))
      p    = rho*c*c/GAMMA
      vt   = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/(GAMMA-1.0)+0.5*rho*(vn*vn+vt*vt)


    case(BDRC_SUBOUTLET)
      !-------------------------------------------------------------------------

      ! The subsonic outlet conditions follow the thesis
      !
      ! `Adaptive Finite Element Solution Algorithm
      !  for the Euler Equations`, R.A. Shapiro

      ! The specified exit static/pressure is Deval=[ps]
      ps = DbdrValue(1)

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      v2  = Du(3)/rho
      E   = Du(4)/rho
      p   = (GAMMA-1.0)*rho*(E-0.5*(v1*v1+v2*v2))

      vn  = DbdrNormal(1)*v1+DbdrNormal(2)*v2
      vt  = DbdrNormal(2)*v1-DbdrNormal(1)*v2
      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))

      ! Compute Riemann invariants based on the solution values and prescribed exit pressure
      W(2) = 2*c/(GAMMA-1.0)-vn
      W(3) = p/(rho**GAMMA)
      W(4) = vt
      W(1) = 4/(GAMMA-1.0)*sqrt(max(GAMMA*ps/rho*(p/ps)**(1.0/GAMMA), SYS_EPSREAL))-W(2)

      ! Transform back into conservative variables
      vn  = 0.5*(W(1)-W(2))
      c   = 0.25*(GAMMA-1.0)*(W(1)+W(2))
      rho = (c*c/GAMMA/W(3))**(1.0/(GAMMA-1.0))
      p   = rho*c*c/GAMMA
      vt  = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/(GAMMA-1.0)+0.5*rho*(vn*vn+vt*vt)


    case DEFAULT
      call output_line('Unsupported type of boundary condition!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcBoundaryvalues2d')
      call sys_halt()
    end select

  end subroutine hydro_calcBoundaryvalues2d

  ! ***************************************************************************

!<subroutine>

  subroutine hydro_coeffVectorBdr2d_sim(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

    use basicgeometry
    use boundary
    use boundarycondaux
    use collection
    use domainintegration
    use feevaluation
    use fparser
    use scalarpde
    use spatialdiscretisation
    use triangulation

!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! This routine handles the constant velocities in the primal problem.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(nblocks,itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:,:), intent(out) :: Dcoefficients
!</output>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution
    real(DP), dimension(:,:), pointer :: Daux1
    real(DP), dimension(:,:,:), pointer :: Daux2
    real(DP), dimension(NVAR2D) :: DstateI,DstateM,Dflux,Diff
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dnx,dny,dtime,dscale,pI,cI,rM,pM,cM,dvnI,dvtI,dvnM,dvtM,w1,w4
    integer :: ibdrtype,isegment,iel,ipoint,ndim,ivar,nvar,iexpr,nmaxExpr

#ifndef HYDRO_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DHYDRO_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'hydro_coeffVectorBdr2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first quick access string
    ! value holds the name of the function parser in the collection.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first quick access vector
    ! points to the solution vector
    p_rsolution => rcollection%p_rvectorQuickAccess1

    ! Check if the solution is given in block or interleaved format
    if (p_rsolution%nblocks .eq. 1) then
      nvar = p_rsolution%RvectorBlock(1)%NVAR
    else
      nvar = p_rsolution%nblocks
    end if
    
    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first three quick access integer values hold:
    ! - the type of boundary condition
    ! - the segment number
    ! - the maximum number of expressions
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)
    nmaxExpr = rcollection%IquickAccess(3)

    if (p_rsolution%nblocks .eq. 1) then

      !-------------------------------------------------------------------------
      ! Solution is stored in interleaved format
      !-------------------------------------------------------------------------
      
      ! Allocate temporal memory
      allocate(Daux1(ubound(Dpoints,2)*nvar, ubound(Dpoints,3)))
      
      ! Evaluate the solution in the cubature points on the boundary
      call fevl_evaluate_sim(DER_FUNC, Daux1, p_rsolution%RvectorBlock(1),&
          Dpoints, rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! What type of boundary conditions are we?
      select case(iand(ibdrtype, BDRC_TYPEMASK))
        
      case (BDRC_FREESTREAM)
        !-----------------------------------------------------------------------
        ! Free-stream boundary conditions:
        !
        ! Compute the Riemann invariants based on the computed (internal)
        ! state vector and the given freestream state vector and select
        ! the Riemman invariant for each characteristic fields based on 
        ! the sign of the corresponding eigenvalue.
        
        ! Initialize values for function parser
        Dvalue = 0.0
        Dvalue(NDIM3D+1) = dtime

        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)

        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)

            ! Get the normal vector in the point from the boundary
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)
            
            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute free stream values from function parser given in
            ! term of the primitive variables [rho,v1,v2,p]
            do iexpr = 1, 4
              call fparser_evalFunction(p_rfparser,&
                  nmaxExpr*(isegment-1)+iexpr, Dvalue, DstateM(iexpr))
            end do

            ! Compute auxiliary quantities based on free stream state vector
            rM = DstateM(1)
            pM = DstateM(4)
            cM = sqrt(max(GAMMA*pM/rM, SYS_EPSREAL))
            dvnM =  dnx*DstateM(2)+dny*DstateM(3)
            dvtM = -dny*DstateM(2)+dnx*DstateM(3)
            
            ! Compute auxiliary quantities based on internal state vector
            pI = (GAMMA-1.0)*(Daux1((ipoint-1)*NVAR2D+4,iel)-0.5*&
                     (Daux1((ipoint-1)*NVAR2D+2,iel)**2+&
                      Daux1((ipoint-1)*NVAR2D+3,iel)**2))/&
                 Daux1((ipoint-1)*NVAR2D+1,iel)
            cI = sqrt(max(GAMMA*pI/Daux1((ipoint-1)*NVAR2D+1,iel), SYS_EPSREAL))

            ! Compute the normal and tangential velocities based
            ! on internal state vector
            dvnI = ( dnx*Daux1((ipoint-1)*NVAR2D+2,iel)+&
                     dny*Daux1((ipoint-1)*NVAR2D+3,iel) )/&
                     Daux1((ipoint-1)*NVAR2D+1,iel)
            dvtI = (-dny*Daux1((ipoint-1)*NVAR2D+2,iel)+&
                     dnx*Daux1((ipoint-1)*NVAR2D+3,iel) )/&
                     Daux1((ipoint-1)*NVAR2D+1,iel)
            
            ! Select free stream or computed Riemann invariant depending
            ! on the sign of the corresponding eigenvalue
            if (dvnI .lt. cI) then
              DstateM(1) = dvnM-(2.0/(GAMMA-1.0))*cM
            else
              DstateM(1) = dvnI-(2.0/(GAMMA-1.0))*cI
            end if

            if (dvnI .lt. SYS_EPSREAL) then
              DstateM(2) = pM/(rM**GAMMA)
              DstateM(3) = dvtM
            else
              DstateM(2) = pI/(Daux1((ipoint-1)*NVAR2D+1,iel)**GAMMA)
              DstateM(3) = dvtI
            end if

            if (dvnI .lt. -cI) then
              DstateM(4) = dvnM+(2.0/(GAMMA-1.0))*cM
            else
              DstateM(4) = dvnI+(2.0/(GAMMA-1.0))*cI
            end if

            ! Convert Riemann invariants into conservative state variables
            cM = 0.25*(GAMMA-1.0)*(DstateM(4)-DstateM(1))
            rM = (cM*cM/(GAMMA*DstateM(2)))**(1.0/(GAMMA-1.0))
            pM = rM*cM*cM/GAMMA
            dvnM = 0.5*(DstateM(1)+DstateM(4))
            dvtM = DstateM(3)

            ! Setup the state vector based on Riemann invariants
            DstateM(1) = rM
            DstateM(2) = rM*( dnx*dvnM+dny*dvtM)
            DstateM(3) = rM*(-dny*dvnM+dnx*dvtM)
            DstateM(4) = pM/(GAMMA-1.0)+0.5*(dvnM*dvnM+dvtM*dvtM)
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux1((ipoint-1)*NVAR2D+1,iel)
            DstateI(2) = Daux1((ipoint-1)*NVAR2D+2,iel)
            DstateI(3) = Daux1((ipoint-1)*NVAR2D+3,iel)
            DstateI(4) = Daux1((ipoint-1)*NVAR2D+4,iel)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*0.5*(Dflux-Diff)
          end do
        end do

        
      case (BDRC_FREESLIP)
        !-----------------------------------------------------------------------
        ! Free-slip boundary condition:
        !
        ! Compute the mirrored state vector based on the values of the
        ! computed state vector and use an approximate Riemann solver
        
        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)
            
            ! Get the normal vector in the point from the boundary
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux1((ipoint-1)*NVAR2D+1,iel)
            DstateI(2) = Daux1((ipoint-1)*NVAR2D+2,iel)
            DstateI(3) = Daux1((ipoint-1)*NVAR2D+3,iel)
            DstateI(4) = Daux1((ipoint-1)*NVAR2D+4,iel)
            
            ! Compute the normal and tangential velocities based
            ! on the internal state vector
            dvnI = ( dnx*Daux1((ipoint-1)*NVAR2D+2,iel)+&
                     dny*Daux1((ipoint-1)*NVAR2D+3,iel) )/&
                     Daux1((ipoint-1)*NVAR2D+1,iel)
            dvtI = (-dny*Daux1((ipoint-1)*NVAR2D+2,iel)+&
                     dnx*Daux1((ipoint-1)*NVAR2D+3,iel) )/&
                     Daux1((ipoint-1)*NVAR2D+1,iel)

            ! Compute the mirrored state vector
            DstateM(1) = DstateI(1)
            DstateM(2) = DstateM(1)*(-dvnI*dnx - dvtI*dny)
            DstateM(3) = DstateM(1)*(-dvnI*dny + dvtI*dnx)
            DstateM(4) = DstateI(4)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)

            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*0.5*(Dflux-Diff)
          end do
        end do


      case (BDRC_SUPERINLET)
        !-----------------------------------------------------------------------
        ! Supersonic inlet boundary conditions:
        !
        ! Prescribe the state vector in conservative variables
        
        ! Initialize values for function parser
        Dvalue = 0.0
        Dvalue(NDIM3D+1) = dtime

        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)

        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)

            ! Get the normal vector in the point from the boundary
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)

            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute boundary values from function parser given in
            ! term of the primitive variables [rho,v1,v2,p]
            do iexpr = 1, 4
              call fparser_evalFunction(p_rfparser,&
                  nmaxExpr*(isegment-1)+iexpr, Dvalue, DstateM(iexpr))
            end do

            ! Compute convervative variables
            DstateM(4) = DstateM(4)*(1.0/(GAMMA-1.0))&
                       + DstateM(1)*0.5*(DstateM(2)**2+DstateM(3)**2)
            DstateM(2) = DstateM(1)*DstateM(2)
            DstateM(3) = DstateM(1)*DstateM(3)

            ! Setup the computed internal state vector
            DstateI(1) = Daux1((ipoint-1)*NVAR2D+1,iel)
            DstateI(2) = Daux1((ipoint-1)*NVAR2D+2,iel)
            DstateI(3) = Daux1((ipoint-1)*NVAR2D+3,iel)
            DstateI(4) = Daux1((ipoint-1)*NVAR2D+4,iel)

            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*0.5*(Dflux-Diff)
          end do
        end do

        
      case (BDRC_SUPEROUTLET)
        !-----------------------------------------------------------------------
        ! Supersonic outlet boundary conditions:
        !
        ! Evaluate the boundary fluxes based on the computed state vector

        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)

            ! Get the normal vector in the point from the boundary
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)
        
            ! Setup the computed internal state vector
            DstateI(1) = Daux1((ipoint-1)*NVAR2D+1,iel)
            DstateI(2) = Daux1((ipoint-1)*NVAR2D+2,iel)
            DstateI(3) = Daux1((ipoint-1)*NVAR2D+3,iel)
            DstateI(4) = Daux1((ipoint-1)*NVAR2D+4,iel)

            ! Assemble Galerkin fluxes at the boundary
            call doGalerkinFlux(DstateI, dnx, dny, Dflux)

            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*Dflux
          end do
        end do

        
      case (BDRC_SUBINLET)
        !-----------------------------------------------------------------------
        ! Subsonic pressure-density inlet boundary conditions:
        !
        ! Prescribe the density, pressure and tangential velocity at the inlet
        
        ! Initialize values for function parser
        Dvalue = 0.0
        Dvalue(NDIM3D+1) = dtime

        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)

        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)

            ! Get the normal vector in the point from the boundary
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)

            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute boundary values from function parser given in
            ! term of the density, pressure and tangential velocity
            do iexpr = 1, 3
              call fparser_evalFunction(p_rfparser,&
                  nmaxExpr*(isegment-1)+iexpr, Dvalue, DstateM(iexpr))
            end do

            ! Compute auxiliary quantities based on prescribed boundary values
            rM = DstateM(1)
            pM = DstateM(2)
            cM = sqrt(max(GAMMA*pM/rM, SYS_EPSREAL))
            dvtM = DstateM(3)

            ! Compute the normal velocity based on the internal state vector
            dvnI = ( dnx*Daux1((ipoint-1)*NVAR2D+2,iel)+&
                     dny*Daux1((ipoint-1)*NVAR2D+3,iel) )/&
                     Daux1((ipoint-1)*NVAR2D+1,iel)

            ! Compute the speed of sound based on the internal state vector
            cI = sqrt(max(GAMMA*pI/Daux1((ipoint-1)*NVAR2D+1,iel), SYS_EPSREAL))

            ! Compute fourth Riemann invariant based on the internal state vector
            w4 = dvnI+(2.0/(GAMMA-1.0))*cI

            ! Compute the first Riemann invariant based on the first Riemann
            ! invariant and the prescribed boundary values
            w1 = w4-2*(2.0/(GAMMA-1.0))*cM

            ! Setup the state vector based on Rimann invariants
            DstateM(1) = rM
            DstateM(2) = rM*(dnx*0.5*(w1+w4)-dny*dvtM)
            DstateM(3) = rM*(dny*0.5*(w1+w4)+dnx*dvtM)
            DstateM(4) = (1.0/(GAMMA-1.0))*pM+0.5*rM*((dnx*0.5*(w1+w4)-dny*dvtM)**2+&
                                          (dny*0.5*(w1+w4)+dnx*dvtM)**2)

            ! Setup the computed internal state vector
            DstateI(1) = Daux1((ipoint-1)*NVAR2D+1,iel)
            DstateI(2) = Daux1((ipoint-1)*NVAR2D+2,iel)
            DstateI(3) = Daux1((ipoint-1)*NVAR2D+3,iel)
            DstateI(4) = Daux1((ipoint-1)*NVAR2D+4,iel)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*0.5*(Dflux-Diff)
          end do
        end do


      case (BDRC_SUBOUTLET)
        !-----------------------------------------------------------------------
        ! Subsonic pressure outlet boundary condition:
        !
        ! Prescribe the pressure at the outlet

        ! Initialize values for function parser
        Dvalue = 0.0
        Dvalue(NDIM3D+1) = dtime

        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)

        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)

            ! Get the normal vector in the point from the boundary
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)

            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute pressure value from function parser
            call fparser_evalFunction(p_rfparser,&
                nmaxExpr*(isegment-1)+1, Dvalue, pM)

            ! Compute auxiliary quantities based on internal state vector
            pI = (GAMMA-1.0)*(Daux1((ipoint-1)*NVAR2D+4,iel)-0.5*&
                     (Daux1((ipoint-1)*NVAR2D+2,iel)**2+&
                      Daux1((ipoint-1)*NVAR2D+3,iel)**2))/&
                 Daux1((ipoint-1)*NVAR2D+1,iel)
            cI = sqrt(max(GAMMA*pI/Daux1((ipoint-1)*NVAR2D+1,iel), SYS_EPSREAL))

            ! Compute the normal and tangential velocities based
            ! on internal state vector
            dvnI = ( dnx*Daux1((ipoint-1)*NVAR2D+2,iel)+&
                     dny*Daux1((ipoint-1)*NVAR2D+3,iel) )/&
                     Daux1((ipoint-1)*NVAR2D+1,iel)
            dvtI = (-dny*Daux1((ipoint-1)*NVAR2D+2,iel)+&
                     dnx*Daux1((ipoint-1)*NVAR2D+3,iel) )/&
                     Daux1((ipoint-1)*NVAR2D+1,iel)
            
            ! Compute three Riemann invariants based on internal state vector
            DstateM(2) = pI/(Daux1((ipoint-1)*NVAR2D+1,iel)**GAMMA)
            DstateM(3) = dvtI
            DstateM(4) = dvnI+(2.0/(GAMMA-1.0))*cI
            
            ! Compute first Riemann invariant based on fourth Riemann invariant,
            ! the computed density and pressure and the prescribed exit pressure
            DstateM(1) = DstateM(4)-2*(2.0/(GAMMA-1.0))*sqrt(max(SYS_EPSREAL,&
                GAMMA*pM/Daux1((ipoint-1)*NVAR2D+1,iel)*(pI/pM)**(1.0/GAMMA)))

            ! Convert Riemann invariants into conservative state variables
            cM = 0.25*(GAMMA-1.0)*(DstateM(4)-DstateM(1))
            rM = (cM*cM/(GAMMA*DstateM(2)))**(1.0/(GAMMA-1.0))
            pM = rM*cM*cM/GAMMA
            dvnM = 0.5*(DstateM(1)+DstateM(4))
            dvtM = DstateM(3)

            ! Setup the state vector based on Riemann invariants
            DstateM(1) = rM
            DstateM(2) = rM*( dnx*dvnM+dny*dvtM)
            DstateM(3) = rM*(-dny*dvnM+dnx*dvtM)
            DstateM(4) = pM/(GAMMA-1.0)+0.5*(dvnM*dvnM+dvtM*dvtM)
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux1((ipoint-1)*NVAR2D+1,iel)
            DstateI(2) = Daux1((ipoint-1)*NVAR2D+2,iel)
            DstateI(3) = Daux1((ipoint-1)*NVAR2D+3,iel)
            DstateI(4) = Daux1((ipoint-1)*NVAR2D+4,iel)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*0.5*(Dflux-Diff)
          end do
        end do
                
      case default
        call output_line('Invalid type of boundary conditions!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_coeffVectorBdr2d_sim')
        call sys_halt()
        
      end select

      ! Deallocate temporal memory
      deallocate(Daux1)
        
    else

      !-------------------------------------------------------------------------
      ! Solution is stored in block format
      !-------------------------------------------------------------------------
      
      ! Allocate temporal memory
      allocate(Daux2(ubound(Dpoints,2), ubound(Dpoints,3), nvar))
      
      ! Evaluate the solution in the cubature points on the boundary
      do ivar = 1, nvar
        call fevl_evaluate_sim(DER_FUNC, Daux2(:,:,ivar),&
            p_rsolution%RvectorBlock(ivar), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      end do

      ! What type of boundary conditions are we?
      select case(iand(ibdrtype, BDRC_TYPEMASK))
        
      case (BDRC_FREESTREAM)
        !-----------------------------------------------------------------------
        ! Free-stream boundary conditions:
        !
        ! Compute the Riemann invariants based on the computed (internal)
        ! state vector and the given freestream state vector and select
        ! the Riemman invariant for each characteristic fields based on 
        ! the sign of the corresponding eigenvalue.

        ! Initialize values for function parser
        Dvalue = 0.0
        Dvalue(NDIM3D+1) = dtime

        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)

        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)

            ! Get the normal vector in the point from the boundary
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)
            
            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute free stream values from function parser given in
            ! term of the primitive variables [rho,v1,v2,p]
            do iexpr = 1, 4
              call fparser_evalFunction(p_rfparser,&
                  nmaxExpr*(isegment-1)+iexpr, Dvalue, DstateM(iexpr))
            end do

            ! Compute auxiliary quantities based on free stream state vector
            rM = DstateM(1)
            pM = DstateM(4)
            cM = sqrt(max(GAMMA*pM/rM, SYS_EPSREAL))
            dvnM =  dnx*DstateM(2)+dny*DstateM(3)
            dvtM = -dny*DstateM(2)+dnx*DstateM(3)

            ! Compute auxiliary quantities based on internal state vector
            pI = (GAMMA-1.0)*(Daux2(ipoint,iel,4)-0.5*&
                     (Daux2(ipoint,iel,2)**2+&
                      Daux2(ipoint,iel,3)**2))/Daux2(ipoint,iel,1)
            cI = sqrt(max(GAMMA*pI/Daux2(ipoint,iel,1), SYS_EPSREAL))

            ! Compute the normal and tangential velocities based
            ! on internal state vector
            dvnI = ( dnx*Daux2(ipoint,iel,2)+&
                     dny*Daux2(ipoint,iel,3))/Daux2(ipoint,iel,1)
            dvtI = (-dny*Daux2(ipoint,iel,2)+&
                     dnx*Daux2(ipoint,iel,3))/Daux2(ipoint,iel,1)

            ! Select free stream or computed Riemann invariant depending
            ! on the sign of the corresponding eigenvalue
            if (dvnI .lt. cI) then
              DstateM(1) = dvnM-(2.0/(GAMMA-1.0))*cM
            else
              DstateM(1) = dvnI-(2.0/(GAMMA-1.0))*cI
            end if

            if (dvnI .lt. SYS_EPSREAL) then
              DstateM(2) = pM/(rM**GAMMA)
              DstateM(3) = dvtM
            else
              DstateM(2) = pI/(Daux2(ipoint,iel,1)**GAMMA)
              DstateM(3) = dvtI
            end if

            if (dvnI .lt. -cI) then
              DstateM(4) = dvnM+(2.0/(GAMMA-1.0))*cM
            else
              DstateM(4) = dvnI+(2.0/(GAMMA-1.0))*cI
            end if

            ! Convert Riemann invariants into conservative state variables
            cM = 0.25*(GAMMA-1.0)*(DstateM(4)-DstateM(1))
            rM = (cM*cM/(GAMMA*DstateM(2)))**(1.0/(GAMMA-1.0))
            pM = rM*cM*cM/GAMMA
            dvnM = 0.5*(DstateM(1)+DstateM(4))
            dvtM = DstateM(3)

            ! Setup the state vector based on Riemann invariants
            DstateM(1) = rM
            DstateM(2) = rM*( dnx*dvnM+dny*dvtM)
            DstateM(3) = rM*(-dny*dvnM+dnx*dvtM)
            DstateM(4) = pM/(GAMMA-1.0)+0.5*(dvnM*dvnM+dvtM*dvtM)
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux2(ipoint,iel,1)
            DstateI(2) = Daux2(ipoint,iel,2)
            DstateI(3) = Daux2(ipoint,iel,3)
            DstateI(4) = Daux2(ipoint,iel,4)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*0.5*(Dflux-Diff)
          end do
        end do


      case (BDRC_FREESLIP)
       
        !-----------------------------------------------------------------------
        ! Free-slip boundary condition:
        !
        ! Compute the mirrored state vector based on the values of the
        ! computed state vector and use an approximate Riemann solver
        
        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)
            
            ! Get the normal vector in the point from the boundary
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux2(ipoint,iel,1)
            DstateI(2) = Daux2(ipoint,iel,2)
            DstateI(3) = Daux2(ipoint,iel,3)
            DstateI(4) = Daux2(ipoint,iel,4)

            ! Compute the normal and tangential velocities based
            ! on the internal state vector
            dvnI = ( dnx*Daux2(ipoint,iel,2)+&
                     dny*Daux2(ipoint,iel,3) )/Daux2(ipoint,iel,1)
            dvtI = (-dny*Daux2(ipoint,iel,2)+&
                      dnx*Daux2(ipoint,iel,3) )/Daux2(ipoint,iel,1)

            ! Compute the mirrored state vector
            DstateM(1) = DstateI(1)
            DstateM(2) = DstateM(1)*(-dvnI*dnx - dvtI*dny)
            DstateM(3) = DstateM(1)*(-dvnI*dny + dvtI*dnx)
            DstateM(4) = DstateI(4)

            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*0.5*(Dflux-Diff)
          end do
        end do


        case (BDRC_SUPERINLET)
        !-----------------------------------------------------------------------
        ! Supersonic inlet boundary conditions:
        !
        ! Prescribe the state vector in conservative variables
        
        ! Initialize values for function parser
        Dvalue = 0.0
        Dvalue(NDIM3D+1) = dtime

        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)

        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)

            ! Get the normal vector in the point from the boundary
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)

            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute boundary values from function parser given in
            ! term of the primitive variables [rho,v1,v2,p]
            do iexpr = 1, 4
              call fparser_evalFunction(p_rfparser,&
                  nmaxExpr*(isegment-1)+iexpr, Dvalue, DstateM(iexpr))
            end do

            ! Compute convervative variables
            DstateM(4) = DstateM(4)*(1.0/(GAMMA-1.0))&
                       + DstateM(1)*0.5*(DstateM(2)**2+DstateM(3)**2)
            DstateM(2) = DstateM(1)*DstateM(2)
            DstateM(3) = DstateM(1)*DstateM(3)
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux2(ipoint,iel,1)
            DstateI(2) = Daux2(ipoint,iel,2)
            DstateI(3) = Daux2(ipoint,iel,3)
            DstateI(4) = Daux2(ipoint,iel,4)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*0.5*(Dflux-Diff)
          end do
        end do


      case (BDRC_SUPEROUTLET)
        !-----------------------------------------------------------------------
        ! Supersonic outlet boundary conditions:
        !
        ! Evaluate the boundary fluxes based on the computed state vector

        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)

            ! Get the normal vector in the point from the boundary
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)
        
            ! Setup the computed internal state vector
            DstateI(1) = Daux2(ipoint,iel,1)
            DstateI(2) = Daux2(ipoint,iel,2)
            DstateI(3) = Daux2(ipoint,iel,3)
            DstateI(4) = Daux2(ipoint,iel,4)

            ! Assemble Galerkin fluxes at the boundary
            call doGalerkinFlux(DstateI, dnx, dny, Dflux)

            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*Dflux
          end do
        end do


      case (BDRC_SUBINLET)
        !-----------------------------------------------------------------------
        ! Subsonic pressure-density inlet boundary conditions:
        !
        ! Prescribe the density, pressure and tangential velocity at the inlet
        
        ! Initialize values for function parser
        Dvalue = 0.0
        Dvalue(NDIM3D+1) = dtime

        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)

        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)

            ! Get the normal vector in the point from the boundary
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)

            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute boundary values from function parser given in
            ! term of the density, pressure and tangential velocity
            do iexpr = 1, 3
              call fparser_evalFunction(p_rfparser,&
                  nmaxExpr*(isegment-1)+iexpr, Dvalue, DstateM(iexpr))
            end do

            ! Compute auxiliary quantities based on prescribed boundary values
            rM = DstateM(1)
            pM = DstateM(2)
            cM = sqrt(max(GAMMA*pM/rM, SYS_EPSREAL))
            dvtM = DstateM(3)

            ! Compute the normal velocity based on the internal state vector
            dvnI = ( dnx*Daux2(ipoint,iel,2)+&
                     dny*Daux2(ipoint,iel,3) )/Daux2(ipoint,iel,1)

            ! Compute the speed of sound based on the internal state vector
            cI = sqrt(max(GAMMA*pI/Daux2(ipoint,iel,1), SYS_EPSREAL))

            ! Compute fourth Riemann invariant based on the internal state vector
            w4 = dvnI+(2.0/(GAMMA-1.0))*cI

            ! Compute the first Riemann invariant based on the first Riemann
            ! invariant and the prescribed boundary values
            w1 = w4-2*(2.0/(GAMMA-1.0))*cM

            ! Setup the state vector based on Rimann invariants
            DstateM(1) = rM
            DstateM(2) = rM*(dnx*0.5*(w1+w4)-dny*dvtM)
            DstateM(3) = rM*(dny*0.5*(w1+w4)+dnx*dvtM)
            DstateM(4) = (1.0/(GAMMA-1.0))*pM+0.5*rM*((dnx*0.5*(w1+w4)-dny*dvtM)**2+&
                                          (dny*0.5*(w1+w4)+dnx*dvtM)**2)

            ! Setup the computed internal state vector
            DstateI(1) = Daux2(ipoint,iel,1)
            DstateI(2) = Daux2(ipoint,iel,2)
            DstateI(3) = Daux2(ipoint,iel,3)
            DstateI(4) = Daux2(ipoint,iel,4)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*0.5*(Dflux-Diff)
          end do
        end do


      case (BDRC_SUBOUTLET)
        !-----------------------------------------------------------------------
        ! Subsonic pressure outlet boundary condition:
        !
        ! Prescribe the pressure at the outlet

        ! Initialize values for function parser
        Dvalue = 0.0
        Dvalue(NDIM3D+1) = dtime

        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)

        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)

            ! Get the normal vector in the point from the boundary
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)

            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute pressure value from function parser
            call fparser_evalFunction(p_rfparser,&
                nmaxExpr*(isegment-1)+1, Dvalue, pM)

            ! Compute auxiliary quantities based on internal state vector
            pI = (GAMMA-1.0)*(Daux2(ipoint,iel,4)-0.5*&
                     (Daux2(ipoint,iel,2)**2+&
                      Daux2(ipoint,iel,3)**2))/Daux2(ipoint,iel,1)
            cI = sqrt(max(GAMMA*pI/Daux2(ipoint,iel,1), SYS_EPSREAL))

            ! Compute the normal and tangential velocities based
            ! on internal state vector
            dvnI = ( dnx*Daux2(ipoint,iel,2)+&
                     dny*Daux2(ipoint,iel,3) )/Daux2(ipoint,iel,1)
            dvtI = (-dny*Daux2(ipoint,iel,2)+&
                     dnx*Daux2(ipoint,iel,3) )/Daux2(ipoint,iel,1)
            
            ! Compute three Riemann invariants based on internal state vector
            DstateM(2) = pI/(Daux2(ipoint,iel,1)**GAMMA)
            DstateM(3) = dvtI
            DstateM(4) = dvnI+(2.0/(GAMMA-1.0))*cI
            
            ! Compute first Riemann invariant based on fourth Riemann invariant,
            ! the computed density and pressure and the prescribed exit pressure
            DstateM(1) = DstateM(4)-2*(2.0/(GAMMA-1.0))*sqrt(max(SYS_EPSREAL,&
                GAMMA*pM/Daux2(ipoint,iel,1)*(pI/pM)**(1.0/GAMMA)))

            ! Convert Riemann invariants into conservative state variables
            cM = 0.25*(GAMMA-1.0)*(DstateM(4)-DstateM(1))
            rM = (cM*cM/(GAMMA*DstateM(2)))**(1.0/(GAMMA-1.0))
            pM = rM*cM*cM/GAMMA
            dvnM = 0.5*(DstateM(1)+DstateM(4))
            dvtM = DstateM(3)

            ! Setup the state vector based on Riemann invariants
            DstateM(1) = rM
            DstateM(2) = rM*( dnx*dvnM+dny*dvtM)
            DstateM(3) = rM*(-dny*dvnM+dnx*dvtM)
            DstateM(4) = pM/(GAMMA-1.0)+0.5*(dvnM*dvnM+dvtM*dvtM)
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux2(ipoint,iel,1)
            DstateI(2) = Daux2(ipoint,iel,2)
            DstateI(3) = Daux2(ipoint,iel,3)
            DstateI(4) = Daux2(ipoint,iel,4)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*0.5*(Dflux-Diff)
          end do
        end do

      case default
        call output_line('Invalid type of boundary conditions!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_coeffVectorBdr2d_sim')
        call sys_halt()
        
      end select

      ! Deallocate temporal memory
      deallocate(Daux2)
      
    end if

  contains

    ! Here come the working routines

    !***************************************************************************
    ! Approximate Riemann solver along the outward unit normal
    !***************************************************************************

    subroutine doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)

      ! input parameters
      real(DP), dimension(NVAR2D), intent(in) :: DstateI, DstateM
      real(DP), intent(in) :: dnx, dny

      ! output parameters
      real(DP), dimension(NVAR2D), intent(out) :: Dflux, Diff

      ! local variables
      real(DP) :: uI,vI,ru2I,rv2I,uM,vM,ru2M,rv2M,hI,hM
      real(DP) :: cPow2,uPow2,vPow2,c_IM,h_IM,q_IM,u_IM,v_IM
      real(DP) :: l1,l2,l3,l4,w1,w2,w3,w4,aux,aux1,aux2,dveln
      
      ! Compute auxiliary quantities
      uI = DstateI(2)/DstateI(1)
      vI = DstateI(3)/DstateI(1)
      ru2I = uI*DstateI(2)
      rv2I = vI*DstateI(3)
      
      ! Compute auxiliary quantities
      uM = DstateM(2)/DstateM(1)
      vM = DstateM(3)/DstateM(1)
      ru2M = uM*DstateM(2)
      rv2M = vM*DstateM(3)

      ! Calculate $\frac12{\bf n}\cdot[{\bf F}(U_I)+{\bf F}(U_M)]$
      Dflux(1) = dnx*(DstateI(2)+DstateM(2))+&
                 dny*(DstateI(3)+DstateM(3))
      Dflux(2) = dnx*((GAMMA-1.0)*(DstateI(4)+DstateM(4))-&
                      (GAMMA-3.0)/2.0*(ru2I+ru2M)-&
                      (GAMMA-1.0)/2.0*(rv2I+rv2M))+&
                 dny*(DstateI(2)*vI+DstateM(2)*vM)
      Dflux(3) = dnx*(DstateI(3)*uI+DstateM(3)*uM)+&
                 dny*((GAMMA-1.0)*(DstateI(4)+DstateM(4))-&
                      (GAMMA-3.0)/2.0*(rv2I+rv2M)-&
                      (GAMMA-1.0)/2.0*(ru2I+ru2M))
      Dflux(4) = dnx*((GAMMA*DstateI(4)-(GAMMA-1.0)/2.0*(ru2I+rv2I))*uI+&
                      (GAMMA*DstateM(4)-(GAMMA-1.0)/2.0*(ru2M+rv2M))*uM)+&
                 dny*((GAMMA*DstateI(4)-(GAMMA-1.0)/2.0*(ru2I+rv2I))*vI+&
                      (GAMMA*DstateM(4)-(GAMMA-1.0)/2.0*(ru2M+rv2M))*vM)
      
      ! Compute Roe mean values
      aux  = sqrt(max(DstateI(1)/DstateM(1), SYS_EPSREAL))
      u_IM = (aux*uI+uM)/(aux+1.0)
      v_IM = (aux*vI+vM)/(aux+1.0)
      hI   = GAMMA*DstateI(4)/DstateI(1)-(GAMMA-1.0)/2.0*(uI**2+vI**2)
      hM   = GAMMA*DstateM(4)/DstateM(1)-(GAMMA-1.0)/2.0*(uM**2+vM**2)
      h_IM =(aux*hI+hM)/(aux+1.0)
      
      ! Compute auxiliary variables
      uPow2 = u_IM*u_IM
      vPow2 = v_IM*v_IM
      q_IM  = 0.5*(uPow2+vPow2)
      cPow2 = max((GAMMA-1.0)*(H_IM-q_IM), SYS_EPSREAL)
      c_IM  = sqrt(cPow2)

      ! Compute normal velocity
      dveln =  dnx*uI+dny*vI

      ! Compute eigenvalues
      l1 = abs(dveln-c_IM)
      l2 = abs(dveln)
      l3 = abs(dveln+c_IM)
      l4 = abs(dveln)
      
      ! Compute solution difference U_M-U_I
      Diff = DstateM-DstateI
      
      ! Compute auxiliary quantities for characteristic variables
      aux1 = (GAMMA-1.0)/2.0/cPow2*(q_IM*Diff(1)-u_IM*Diff(2)-v_IM*Diff(3)+Diff(4))
      aux2 = 0.5*(aux*Diff(1)-dnx*Diff(2)-dny*Diff(3))/c_IM
      
      ! Compute characteristic variables multiplied by the corresponding eigenvalue
      w1 = l1 * (aux1 + aux2)
      w2 = l2 * ((1.0-(GAMMA-1.0)*q_IM/cPow2)*Diff(1)+(GAMMA-1.0)*(u_IM*Diff(2)+v_IM*Diff(3)-Diff(4))/cPow2)
      w3 = l3 * (aux1 - aux2)
      w4 = l4 * ((dnx*v_IM-dny*u_IM)*Diff(1)+dny*Diff(2)-dnx*Diff(3))
      
      ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
      Diff(1) = w1 + w2 + w3
      Diff(2) = (u_IM-c_IM*dnx)*w1 + u_IM*w2 + (u_IM+c_IM*dnx)*w3 + dny*w4
      Diff(3) = (v_IM-c_IM*dny)*w1 + v_IM*w2 + (v_IM+c_IM*dny)*w3 - dnx*w4
      Diff(4) = (H_IM-c_IM*dveln)*w1  + q_IM*w2 +&
                (H_IM+c_IM*dveln)*w3  + (u_IM*dny-v_IM*dnx)*w4

    end subroutine doRiemannSolver

    !***************************************************************************
    ! Compute the Galerkin flux (ued for supersonic outflow)
    !***************************************************************************

    subroutine doGalerkinFlux(Dstate, dnx, dny, Dflux)

      ! input parameters
      real(DP), dimension(NVAR2D), intent(in) :: Dstate
      real(DP), intent(in) :: dnx, dny

      ! output parameters
      real(DP), dimension(NVAR2D), intent(out) :: Dflux

      ! local variables
      real(DP) :: u,v,ru2,rv2
      
      
      ! Compute auxiliary quantities
      u = Dstate(2)/Dstate(1)
      v = Dstate(3)/Dstate(1)
      ru2 = u*Dstate(2)
      rv2 = v*Dstate(3)
      
      ! Calculate ${\bf n}\cdot{\bf F}(U)$
      Dflux(1) = dnx*Dstate(2)+dny*Dstate(3)
      Dflux(2) = dnx*((GAMMA-1.0)*Dstate(4)-(GAMMA-3.0)/2.0*ru2-&
                      (GAMMA-1.0)/2.0*rv2)+dny*Dstate(2)*v
      Dflux(3) = dnx*Dstate(3)*u+dny*((GAMMA-1.0)*Dstate(4)-&
                     (GAMMA-3.0)/2.0*rv2-(GAMMA-1.0)/2.0*ru2)
      Dflux(4) = dnx*(GAMMA*Dstate(4)-(GAMMA-1.0)/2.0*(ru2+rv2))*u+&
                 dny*(GAMMA*Dstate(4)-(GAMMA-1.0)/2.0*(ru2+rv2))*v

    end subroutine doGalerkinFlux
    
  end subroutine hydro_coeffVectorBdr2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine hydro_hadaptCallbackScalar2d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D. The solution vector is assumed
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
            OU_CLASS_WARNING,OU_MODE_STD,'hydro_hadaptCallbackScalar2d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR2D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR2D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      do ivar = 1, NVAR2D
        p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR2D+ivar) = &
            0.5*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR2D+ivar)+&
                    p_Dsolution((rcollection%IquickAccess(3)-1)*NVAR2D+ivar))
      end do

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR2D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      do ivar = 1, NVAR2D
        p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR2D+ivar) = &
            0.25*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR2D+ivar)+&
                     p_Dsolution((rcollection%IquickAccess(3)-1)*NVAR2D+ivar)+&
                     p_Dsolution((rcollection%IquickAccess(4)-1)*NVAR2D+ivar)+&
                     p_Dsolution((rcollection%IquickAccess(5)-1)*NVAR2D+ivar))
      end do

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        do ivar = 1, NVAR2D
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR2D+ivar) = &
              p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR2D+ivar)
        end do
      else
        do ivar = 1, NVAR2D
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR2D+ivar) = 0.0
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)

    end select

  end subroutine hydro_hadaptCallbackScalar2d

  !*****************************************************************************

!<subroutine>

  subroutine hydro_hadaptCallbackBlock2d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D. The solution vector is assumed
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
      if (rsolution%nblocks .ne. NVAR2D) then
        call output_line('Vector is not in block format!',&
            OU_CLASS_WARNING,OU_MODE_STD,'hydro_hadaptCallbackBlock2d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR2D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution, NVAR2D*rcollection%IquickAccess(1),&
            .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      neq = rsolution%NEQ/NVAR2D
      do ivar = 1, NVAR2D
        p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = &
            0.5*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
                    p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(3)) )
      end do

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR2D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      neq = rsolution%NEQ/NVAR2D
      do ivar = 1, NVAR2D
        p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) =&
            0.25*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
                     p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(3))+&
                     p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(4))+&
                     p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(5)) )
      end do

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        neq = rsolution%NEQ/NVAR2D
        do ivar = 1, NVAR2D
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = &
              p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))
        end do
      else
        neq = rsolution%NEQ/NVAR2D
        do ivar = 1, NVAR2D
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = 0.0
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)

    end select

  end subroutine hydro_hadaptCallbackBlock2d

end module hydro_callback2d
