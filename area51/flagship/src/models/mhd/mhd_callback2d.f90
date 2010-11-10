!##############################################################################
!# ****************************************************************************
!# <name> mhd_callback2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible MHDequations in 2D.
!#
!# The following callback functions are available:
!#
!# 1.) mhd_calcFluxGal2d_sim
!#     -> Computes fluxes for standard Galerkin scheme
!#
!# 2.) mhd_calcFluxGalNoBdr2d_sim
!#     -> Computes fluxes for standard Galerkin scheme without
!#        assembling the symmetric boundary contribution
!#
!# 3.) mhd_calcFluxScDiss2d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting scalar artificial viscosities
!#
!# 4.) mhd_calcFluxScDissDiSp2d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting scalar artificial viscosities based on
!#        dimensional splitting approach
!#
!# 5.) mhd_calcFluxRoeDiss2d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities
!#
!# 6.) mhd_calcFluxRoeDissDiSp2d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities based on
!#        dimensional splitting approach
!#
!# 7.) mhd_calcFluxRusDiss2d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion based on
!#        dimensional splitting approach
!#
!# 8.) mhd_calcFluxRusDissDiSp2d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion
!#
!# 9.) mhd_calcMatDiagMatD2d_sim
!#     -> Computes local matrix for diagonal entry
!#
!# 10.) mhd_calcMatDiag2d_sim
!#      -> Computes local matrix for diagonal entry
!#
!# 11.) mhd_calcMatGalMatD2d_sim
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 12.) mhd_calcMatGal2d_sim
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 13.) mhd_calcMatScDissMatD2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 14.) mhd_calcMatScDiss2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 15.) mhd_calcMatRoeDissMatD2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 16.) mhd_calcMatRoeDiss2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 17.) mhd_calcMatRusDissMatD2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov artificial viscosities
!#
!# 18.) mhd_calcMatRusDiss2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov flux artificial viscosities
!#
!# 19.) mhd_calcCharacteristics2d_sim
!#      -> Computes characteristic variables
!#
!# 20.) mhd_calcFluxFCTScalarDiss2d
!#      -> Computes fluxes for FCT algorithm
!#         adopting scalar artificial viscosities
!#
!# 21.) mhd_calcFluxFCTRoeDiss2d
!#      -> Computes fluxes for FCT algorithm
!#         adopting tensorial artificial viscosities
!#
!# 22.) mhd_calcFluxFCTRusanov2d
!#      -> Computes fluxes for FCT algorithm
!#         adopting the Rusanov artificial viscosities
!#
!# 23.) mhd_trafoFluxDensity2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 24.) mhd_trafoDiffDensity2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 25.) mhd_trafoFluxEnergy2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 26.) mhd_trafoDiffEnergy2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 27.) mhd_trafoFluxPressure2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 28.) mhd_trafoFluxVelocity2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 29.) mhd_trafoDiffVelocity2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 30.) mhd_trafoFluxMomentum2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 31.) mhd_trafoDiffMomentum2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 32.) mhd_trafoFluxDenEng2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 33.) mhd_trafoDiffDenEng2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 34.) mhd_trafoFluxDenPre2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 35.) mhd_trafoDiffDenPre2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 36.) mhd_trafoFluxDenPreVel2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 37.) mhd_trafoDiffDenPreVel2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure 
!#         and the velocity
!#
!# 38.) mhd_calcBoundaryvalues2d
!#      -> Computes the boundary values for a given node
!#
!# 39.) mhd_hadaptCallbackScalar2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in interleave format
!#
!# 40.) mhd_hadaptCallbackBlock2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in block format
!#
!# 41.) mhd_coeffVectorBdr2d_sim
!#      -> Calculates the coefficients for the linear form in 2D
!#
!# 42.) mhd_coeffMatrixBdr2d_sim
!#      -> Calculates the coefficients for the bilinear form in 2D
!# </purpose>
!##############################################################################

module mhd_callback2d

#include "mhd.h"

  use boundarycondaux
  use collection
  use derivatives
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
  public :: mhd_calcFluxGal2d_sim
  public :: mhd_calcFluxGalNoBdr2d_sim
  public :: mhd_calcFluxScDiss2d_sim
  public :: mhd_calcFluxScDissDiSp2d_sim
  public :: mhd_calcFluxRoeDiss2d_sim
  public :: mhd_calcFluxRoeDissDiSp2d_sim
  public :: mhd_calcFluxRusDiss2d_sim
  public :: mhd_calcFluxRusDissDiSp2d_sim
  public :: mhd_calcMatDiagMatD2d_sim
  public :: mhd_calcMatDiag2d_sim
  public :: mhd_calcMatGalMatD2d_sim
  public :: mhd_calcMatGal2d_sim
  public :: mhd_calcMatScDissMatD2d_sim
  public :: mhd_calcMatScDiss2d_sim
  public :: mhd_calcMatRoeDissMatD2d_sim
  public :: mhd_calcMatRoeDiss2d_sim
  public :: mhd_calcMatRusDissMatD2d_sim
  public :: mhd_calcMatRusDiss2d_sim
  public :: mhd_calcCharacteristics2d_sim
  public :: mhd_calcFluxFCTScalarDiss2d
  public :: mhd_calcFluxFCTRoeDiss2d
  public :: mhd_calcFluxFCTRusanov2d
  public :: mhd_trafoFluxDensity2d_sim
  public :: mhd_trafoFluxEnergy2d_sim
  public :: mhd_trafoFluxPressure2d_sim
  public :: mhd_trafoFluxVelocity2d_sim
  public :: mhd_trafoFluxMomentum2d_sim
  public :: mhd_trafoFluxDenEng2d_sim
  public :: mhd_trafoFluxDenPre2d_sim
  public :: mhd_trafoFluxDenPreVel2d_sim
  public :: mhd_trafoDiffDensity2d_sim
  public :: mhd_trafoDiffEnergy2d_sim
  public :: mhd_trafoDiffPressure2d_sim
  public :: mhd_trafoDiffVelocity2d_sim
  public :: mhd_trafoDiffMomentum2d_sim
  public :: mhd_trafoDiffDenEng2d_sim
  public :: mhd_trafoDiffDenPre2d_sim
  public :: mhd_trafoDiffDenPreVel2d_sim
  public :: mhd_calcBoundaryvalues2d
  public :: mhd_coeffVectorBdr2d_sim
  public :: mhd_hadaptCallbackScalar2d
  public :: mhd_hadaptCallbackBlock2d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGal2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the standard
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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR2D) :: dF1_i,dF2_i,dF1_j,dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij
#endif
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      /   rho*u             \
      !      | rho*u*u + pT - Bx^2 |
      !      | rho*v*u - Bx*By     |
      ! F1 = | rho*w*u - Bx*Bz     |
      !      |         0           |
      !      |    By*u - Bx*v      |
      !      |    Bz*u - Bx*w      |
      !      \ (rho*E+pT)*u -Bx*q  /
      !
      !      /   rho*v             \
      !      | rho*u*v - By*Bx     |
      !      | rho*v*v + pT - By^2 |
      ! F2 = | rho*w*v - By*Bz     |
      !      |    Bx*v - By*u      |
      !      |         0           |
      !      |    Bz*v - By*w      |
      !      \ (rho*E+pT)*u -By*q  /
      !
      ! pT = p + q, q = 1/(2*nu)*|B|^2
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj
      
#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(5) = 0.0
      dF1_i(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi
      dF1_i(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      dF1_i(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi
      
      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(5) = 0.0
      dF1_j(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF1_j(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj
      dF1_j(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj

      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui
      dF2_i(6) = 0.0
      dF2_i(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      dF2_i(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi

      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF2_j(6) = 0.0
      dF2_j(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj
      dF2_j(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj
      
      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*dF2_i )
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF1_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(5) = 0.0
      dF1_ij(6) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi) -&
                  (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj)
      dF1_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF1_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(5) = (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui) -&
                  (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj)
      dF2_ij(6) = 0.0
      dF2_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF2_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
           
      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij)
#endif
    end do

  end subroutine mhd_calcFluxGal2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGalNoBdr2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the TVD
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
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      /   rho*u             \
      !      | rho*u*u + pT - Bx^2 |
      !      | rho*v*u - Bx*By     |
      ! F1 = | rho*w*u - Bx*Bz     |
      !      |         0           |
      !      |    By*u - Bx*v      |
      !      |    Bz*u - Bx*w      |
      !      \ (rho*E+pT)*u -Bx*q  /
      !
      !      /   rho*v             \
      !      | rho*u*v - By*Bx     |
      !      | rho*v*v + pT - By^2 |
      ! F2 = | rho*w*v - By*Bz     |
      !      |    Bx*v - By*u      |
      !      |         0           |
      !      |    Bz*v - By*w      |
      !      \ (rho*E+pT)*u -By*q  /
      !
      ! pT = p + q, q = 1/(2*nu)*|B|^2
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj

      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF1_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(5) = 0.0
      dF1_ij(6) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi) -&
                  (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj)
      dF1_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF1_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(5) = (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui) -&
                  (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj)
      dF2_ij(6) = 0.0
      dF2_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF2_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
      
      ! Assemble symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale *&
          (0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))*dF1_ij+&
           0.5*(DmatrixCoeffsAtEdge(2,1,idx)-DmatrixCoeffsAtEdge(2,2,idx))*dF2_ij)
      DfluxesAtEdge(:,2,idx) = DfluxesAtEdge(:,1,idx)
    end do

  end subroutine mhd_calcFluxGalNoBdr2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxScDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the fluxes for the
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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR2D) :: dF1_i,dF2_i,dF1_j,dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: d_ij,H_ij,q_ij,u_ij,v_ij,anorm,vel_ij,c_ij,aux
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)

     !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      /   rho*u             \
      !      | rho*u*u + pT - Bx^2 |
      !      | rho*v*u - Bx*By     |
      ! F1 = | rho*w*u - Bx*Bz     |
      !      |         0           |
      !      |    By*u - Bx*v      |
      !      |    Bz*u - Bx*w      |
      !      \ (rho*E+pT)*u -Bx*q  /
      !
      !      /   rho*v             \
      !      | rho*u*v - By*Bx     |
      !      | rho*v*v + pT - By^2 |
      ! F2 = | rho*w*v - By*Bz     |
      !      |    Bx*v - By*u      |
      !      |         0           |
      !      |    Bz*v - By*w      |
      !      \ (rho*E+pT)*u -By*q  /
      !
      ! pT = p + q, q = 1/(2*nu)*|B|^2
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(5) = 0.0
      dF1_i(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi
      dF1_i(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      dF1_i(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi
      
      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(5) = 0.0
      dF1_j(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF1_j(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj
      dF1_j(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj

      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui
      dF2_i(6) = 0.0
      dF2_i(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      dF2_i(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi

      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF2_j(6) = 0.0
      dF2_j(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj
      dF2_j(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF1_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(5) = 0.0
      dF1_ij(6) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi) -&
                  (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj)
      dF1_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF1_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(5) = (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui) -&
                  (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj)
      dF2_ij(6) = 0.0
      dF2_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF2_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
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

#ifdef MHD_USE_IBP
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

  end subroutine mhd_calcFluxScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxScDissDiSp2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    

!<description>
    ! This subroutine computes the fluxes for the
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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR2D) :: dF1_i,dF2_i,dF1_j,dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: d_ij,H_ij,q_ij,u_ij,v_ij,aux,c_ij
    integer :: idx
    
    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      /   rho*u             \
      !      | rho*u*u + pT - Bx^2 |
      !      | rho*v*u - Bx*By     |
      ! F1 = | rho*w*u - Bx*Bz     |
      !      |         0           |
      !      |    By*u - Bx*v      |
      !      |    Bz*u - Bx*w      |
      !      \ (rho*E+pT)*u -Bx*q  /
      !
      !      /   rho*v             \
      !      | rho*u*v - By*Bx     |
      !      | rho*v*v + pT - By^2 |
      ! F2 = | rho*w*v - By*Bz     |
      !      |    Bx*v - By*u      |
      !      |         0           |
      !      |    Bz*v - By*w      |
      !      \ (rho*E+pT)*u -By*q  /
      !
      ! pT = p + q, q = 1/(2*nu)*|B|^2
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(5) = 0.0
      dF1_i(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi
      dF1_i(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      dF1_i(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi
      
      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(5) = 0.0
      dF1_j(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF1_j(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj
      dF1_j(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj

      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui
      dF2_i(6) = 0.0
      dF2_i(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      dF2_i(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi

      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF2_j(6) = 0.0
      dF2_j(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj
      dF2_j(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF1_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(5) = 0.0
      dF1_ij(6) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi) -&
                  (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj)
      dF1_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF1_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(5) = (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui) -&
                  (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj)
      dF2_ij(6) = 0.0
      dF2_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF2_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
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

#ifdef MHD_USE_IBP
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

  end subroutine mhd_calcFluxScDissDiSp2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRoeDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the
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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR2D) :: dF1_i,dF2_i,dF1_j,dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: u_ij,v_ij,H_ij,q_ij,c_ij,c2_ij,vel_ij
    real(DP) :: aux,aux1,aux2,anorm
    real(DP) :: l1,l2,l3,l4,w1,w2,w3,w4
    integer :: idx


    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      /   rho*u             \
      !      | rho*u*u + pT - Bx^2 |
      !      | rho*v*u - Bx*By     |
      ! F1 = | rho*w*u - Bx*Bz     |
      !      |         0           |
      !      |    By*u - Bx*v      |
      !      |    Bz*u - Bx*w      |
      !      \ (rho*E+pT)*u -Bx*q  /
      !
      !      /   rho*v             \
      !      | rho*u*v - By*Bx     |
      !      | rho*v*v + pT - By^2 |
      ! F2 = | rho*w*v - By*Bz     |
      !      |    Bx*v - By*u      |
      !      |         0           |
      !      |    Bz*v - By*w      |
      !      \ (rho*E+pT)*u -By*q  /
      !
      ! pT = p + q, q = 1/(2*nu)*|B|^2
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(5) = 0.0
      dF1_i(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi
      dF1_i(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      dF1_i(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi
      
      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(5) = 0.0
      dF1_j(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF1_j(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj
      dF1_j(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj

      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui
      dF2_i(6) = 0.0
      dF2_i(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      dF2_i(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi

      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF2_j(6) = 0.0
      dF2_j(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj
      dF2_j(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF1_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(5) = 0.0
      dF1_ij(6) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi) -&
                  (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj)
      dF1_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF1_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(5) = (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui) -&
                  (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj)
      dF2_ij(6) = 0.0
      dF2_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF2_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
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

#ifdef MHD_USE_IBP
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
        
#ifdef MHD_USE_IBP
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

  end subroutine mhd_calcFluxRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRoeDissDiSp2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    

!<description>
    ! This subroutine computes the fluxes for the
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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR2D) :: dF1_i,dF2_i,dF1_j,dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff1,Diff2
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: u_ij,v_ij,H_ij,q_ij,c_ij,c2_ij
    real(DP) :: aux,aux1,aux2,anorm
    real(DP) :: l1,l2,l3,l4,w1,w2,w3,w4
    integer :: idx

    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      /   rho*u             \
      !      | rho*u*u + pT - Bx^2 |
      !      | rho*v*u - Bx*By     |
      ! F1 = | rho*w*u - Bx*Bz     |
      !      |         0           |
      !      |    By*u - Bx*v      |
      !      |    Bz*u - Bx*w      |
      !      \ (rho*E+pT)*u -Bx*q  /
      !
      !      /   rho*v             \
      !      | rho*u*v - By*Bx     |
      !      | rho*v*v + pT - By^2 |
      ! F2 = | rho*w*v - By*Bz     |
      !      |    Bx*v - By*u      |
      !      |         0           |
      !      |    Bz*v - By*w      |
      !      \ (rho*E+pT)*u -By*q  /
      !
      ! pT = p + q, q = 1/(2*nu)*|B|^2
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(5) = 0.0
      dF1_i(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi
      dF1_i(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      dF1_i(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi
      
      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(5) = 0.0
      dF1_j(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF1_j(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj
      dF1_j(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj

      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui
      dF2_i(6) = 0.0
      dF2_i(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      dF2_i(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi

      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF2_j(6) = 0.0
      dF2_j(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj
      dF2_j(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF1_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(5) = 0.0
      dF1_ij(6) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi) -&
                  (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj)
      dF1_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF1_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(5) = (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui) -&
                  (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj)
      dF2_ij(6) = 0.0
      dF2_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF2_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
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

#ifdef MHD_USE_IBP
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

#ifdef MHD_USE_IBP
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

  end subroutine mhd_calcFluxRoeDissDiSp2d_sim
 
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRusDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the fluxes for the
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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR2D) :: dF1_i,dF2_i,dF1_j,dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: ca1i,ca1j,ca2i,ca2j,cf1i,cf1j,cf2i,cf2j,d_ij
    real(DP) :: aPow2i,aPow2j,astPow2i,astPow2j
    integer :: idx

    do idx = 1, size(DfluxesAtEdge,3)

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      /   rho*u             \
      !      | rho*u*u + pT - Bx^2 |
      !      | rho*v*u - Bx*By     |
      ! F1 = | rho*w*u - Bx*Bz     |
      !      |         0           |
      !      |    By*u - Bx*v      |
      !      |    Bz*u - Bx*w      |
      !      \ (rho*E+pT)*u -Bx*q  /
      !
      !      /   rho*v             \
      !      | rho*u*v - By*Bx     |
      !      | rho*v*v + pT - By^2 |
      ! F2 = | rho*w*v - By*Bz     |
      !      |    Bx*v - By*u      |
      !      |         0           |
      !      |    Bz*v - By*w      |
      !      \ (rho*E+pT)*u -By*q  /
      !
      ! pT = p + q, q = 1/(2*nu)*|B|^2
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(5) = 0.0
      dF1_i(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi
      dF1_i(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      dF1_i(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi
      
      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(5) = 0.0
      dF1_j(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF1_j(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj
      dF1_j(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj

      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui
      dF2_i(6) = 0.0
      dF2_i(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      dF2_i(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi

      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF2_j(6) = 0.0
      dF2_j(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj
      dF2_j(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF1_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(5) = 0.0
      dF1_ij(6) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi) -&
                  (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj)
      dF1_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF1_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(5) = (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui) -&
                  (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj)
      dF2_ij(6) = 0.0
      dF2_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF2_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
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
      
      ! Compute the speed of the Alfven waves in x-direction
      ca1i = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx))
      ca1j = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))

      ! Compute the speed of the Alfven waves in y-direction
      ca2i = abs(Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx))
      ca2j = abs(Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))

! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      aPow2i = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      aPow2j = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
#else
#error "Speed of sound must be implemented!"
#endif

      ! Compute auxiliary quantities
      astPow2i = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)/MYNEWLINE
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + aPow2i
      astPow2j = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)/MYNEWLINE
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + aPow2j

      ! Compute the speed of the fast waves in x-direction
      cf1i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca1i**2)))
      cf1j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca1j**2)))

      ! Compute the speed of the fast waves in y-direction
      cf2i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca2i**2)))
      cf2j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca2j**2)))
      
      ! Scalar dissipation for the Rusanov flux
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj+&
                      DmatrixCoeffsAtEdge(2,1,idx)*vj)+&
                 sqrt((DmatrixCoeffsAtEdge(1,1,idx)**2)*cf1j+&
                      (DmatrixCoeffsAtEdge(2,1,idx)**2)*cf2j),&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui+&
                      DmatrixCoeffsAtEdge(2,2,idx)*vi)+&
                 sqrt((DmatrixCoeffsAtEdge(1,2,idx)**2)*cf1i+&
                      (DmatrixCoeffsAtEdge(2,2,idx)**2)*cf2i) )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------
      
#ifdef MHD_USE_IBP
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

  end subroutine mhd_calcFluxRusDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRusDissDiSp2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the viscid fluxes for the
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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR2D) :: dF1_i,dF2_i,dF1_j,dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: ca1i,ca1j,ca2i,ca2j,cf1i,cf1j,cf2i,cf2j,d_ij
    real(DP) :: aPow2i,aPow2j,astPow2i,astPow2j
    integer :: idx
    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      /   rho*u             \
      !      | rho*u*u + pT - Bx^2 |
      !      | rho*v*u - Bx*By     |
      ! F1 = | rho*w*u - Bx*Bz     |
      !      |         0           |
      !      |    By*u - Bx*v      |
      !      |    Bz*u - Bx*w      |
      !      \ (rho*E+pT)*u -Bx*q  /
      !
      !      /   rho*v             \
      !      | rho*u*v - By*Bx     |
      !      | rho*v*v + pT - By^2 |
      ! F2 = | rho*w*v - By*Bz     |
      !      |    Bx*v - By*u      |
      !      |         0           |
      !      |    Bz*v - By*w      |
      !      \ (rho*E+pT)*u -By*q  /
      !
      ! pT = p + q, q = 1/(2*nu)*|B|^2
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      dF1_i(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2
      dF1_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF1_i(5) = 0.0
      dF1_i(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi
      dF1_i(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      dF1_i(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*ui -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi
      
      dF1_j(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2
      dF1_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_j(5) = 0.0
      dF1_j(6) = Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj
      dF1_j(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj
      dF1_j(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*uj -&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj

      ! Compute fluxes for y-direction
      dF2_i(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2
      dF2_i(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      dF2_i(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui
      dF2_i(6) = 0.0
      dF2_i(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi
      dF2_i(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*vi -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi

      dF2_j(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(2) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(3) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2
      dF2_j(4) = Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_j(5) = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj
      dF2_j(6) = 0.0
      dF2_j(7) = Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj
      dF2_j(8) = (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*vj -&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF1_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui + pi -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj + pj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF1_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF1_ij(5) = 0.0
      dF1_ij(6) = (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi) -&
                  (Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj)
      dF1_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF1_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*ui -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*uj -&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
                        
      ! Compute flux difference for y-direction
      dF2_ij(1) = Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) -&
                  Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      dF2_ij(2) = (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(3) = (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi + pi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)**2) -&
                  (Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj + pj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)**2)
      dF2_ij(4) = (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)) -&
                  (Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))
      dF2_ij(5) = (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*ui) -&
                  (X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*uj)
      dF2_ij(6) = 0.0
      dF2_ij(7) = (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*wi) -&
                  (Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*wj)
      dF2_ij(8) = ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + pi)*vi -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*qi) -&
                  ((TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + pj)*vj -&
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*qj)
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

      ! Compute the speed of the Alfven waves in x-direction
      ca1i = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx))
      ca1j = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))

      ! Compute the speed of the Alfven waves in y-direction
      ca2i = abs(Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx))
      ca2j = abs(Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))

! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      aPow2i = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      aPow2j = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
#else
#error "Speed of sound must be implemented!"
#endif

      ! Compute auxiliary quantities
      astPow2i = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)/MYNEWLINE
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + aPow2i
      astPow2j = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)/MYNEWLINE
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + aPow2j

      ! Compute the speed of the fast waves in x-direction
      cf1i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca1i**2)))
      cf1j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca1j**2)))

      ! Compute the speed of the fast waves in y-direction
      cf2i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca2i**2)))
      cf2j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca2j**2)))

      ! Scalar dissipation
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj)+&
                  abs(DmatrixCoeffsAtEdge(1,1,idx))*cf1j,&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui)+&
                  abs(DmatrixCoeffsAtEdge(1,2,idx))*cf1i )&
           + max( abs(DmatrixCoeffsAtEdge(2,1,idx)*vj)+&
                  abs(DmatrixCoeffsAtEdge(2,1,idx))*cf2j,&
                  abs(DmatrixCoeffsAtEdge(2,2,idx)*vi)+&
                  abs(DmatrixCoeffsAtEdge(2,2,idx))*cf2i )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))
      
      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
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

  end subroutine mhd_calcFluxRusDissDiSp2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiagMatD2d_sim(DdataAtNode, DmatrixCoeffsAtNode,&
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
      vi = X_VELOCITY_1T_FROM_CONSVAR(DdataAtNode,NVAR2D,inode)
      
#ifdef MHD_USE_IBP
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

  end subroutine mhd_calcMatDiagMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiag2d_sim(DdataAtNode,&
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

#ifdef MHD_USE_IBP      
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

  end subroutine mhd_calcMatDiag2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGalMatD2d_sim(DdataAtEdge,&
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

#ifdef MHD_USE_IBP      
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

  end subroutine mhd_calcMatGalMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGal2d_sim(DdataAtEdge,&
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

#ifdef MHD_USE_IBP
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
      
  end subroutine mhd_calcMatGal2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDissMatD2d_sim(DdataAtEdge,&
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

#ifdef MHD_USE_IBP      
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
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
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

  end subroutine mhd_calcMatScDissMatD2d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDiss2d_sim(DdataAtEdge,&
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

#ifdef MHD_USE_IBP
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
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
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

  end subroutine mhd_calcMatScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDissMatD2d_sim(DdataAtEdge,&
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

#ifdef MHD_USE_IBP      
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
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
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

  end subroutine mhd_calcMatRoeDissMatD2d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDiss2d_sim(DdataAtEdge,&
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
      
#ifdef MHD_USE_IBP
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
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
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

  end subroutine mhd_calcMatRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDissMatD2d_sim(DdataAtEdge,&
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

#ifdef MHD_USE_IBP
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

  end subroutine mhd_calcMatRusDissMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDiss2d_sim(DdataAtEdge,&
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

#ifdef MHD_USE_IBP      
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

      DcoefficientsAtEdge(:,1,idx) = 0.0
      DcoefficientsAtEdge( 1,1,idx) = aux1
      DcoefficientsAtEdge( 6,1,idx) = aux1
      DcoefficientsAtEdge(11,1,idx) = aux1
      DcoefficientsAtEdge(16,1,idx) = aux1
    end do

  end subroutine mhd_calcMatRusDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcCharacteristics2d_sim(Dweight, DdataAtEdge,&
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
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
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
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
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
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
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
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx))/MYNEWLINE
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx),MYNEWLINE
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)+MYNEWLINE
               PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))/MYNEWLINE
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

  end subroutine mhd_calcCharacteristics2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTScalarDiss2d(&
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
            PRESSURE_FROM_CONSVAR(U2_i,NVAR2D))/DENSITY_FROM_CONSVAR(U2_i,NVAR2D),MYNEWLINE
           (TOTAL_ENERGY_FROM_CONSVAR(U2_j,NVAR2D)+MYNEWLINE
            PRESSURE_FROM_CONSVAR(U2_j,NVAR2D))/DENSITY_FROM_CONSVAR(U2_j,NVAR2D), aux)
            
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

  end subroutine mhd_calcFluxFCTScalarDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTRoeDiss2d(&
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
              PRESSURE_FROM_CONSVAR(U2_i,NVAR2D))/DENSITY_FROM_CONSVAR(U2_i,NVAR2D),MYNEWLINE
             (TOTAL_ENERGY_FROM_CONSVAR(U2_j,NVAR2D)+MYNEWLINE
              PRESSURE_FROM_CONSVAR(U2_j,NVAR2D))/DENSITY_FROM_CONSVAR(U2_j,NVAR2D), aux)

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
  end subroutine mhd_calcFluxFCTRoeDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTRusanov2d(&
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

  end subroutine mhd_calcFluxFCTRusanov2d

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDensity2d_sim(DdataAtEdge,&
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

  end subroutine mhd_trafoFluxDensity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDensity2d_sim(DdataAtEdge,&
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

  end subroutine mhd_trafoDiffDensity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxEnergy2d_sim(DdataAtEdge,&
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

  end subroutine mhd_trafoFluxEnergy2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffEnergy2d_sim(DdataAtEdge,&
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

  end subroutine mhd_trafoDiffEnergy2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxPressure2d_sim(DdataAtEdge,&
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

  end subroutine mhd_trafoFluxPressure2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffPressure2d_sim(DdataAtEdge,&
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
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
    end do

  end subroutine mhd_trafoDiffPressure2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxVelocity2d_sim(DdataAtEdge,&
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
    
  end subroutine mhd_trafoFluxVelocity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffVelocity2d_sim(DdataAtEdge,&
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
    
  end subroutine mhd_trafoDiffVelocity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxMomentum2d_sim(DdataAtEdge,&
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
    
  end subroutine mhd_trafoFluxMomentum2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffMomentum2d_sim(DdataAtEdge,&
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
    
  end subroutine mhd_trafoDiffMomentum2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDenEng2d_sim(DdataAtEdge,&
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

  end subroutine mhd_trafoFluxDenEng2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDenEng2d_sim(DdataAtEdge,&
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

  end subroutine mhd_trafoDiffDenEng2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDenPre2d_sim(DdataAtEdge,&
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

  end subroutine mhd_trafoFluxDenPre2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDenPre2d_sim(DdataAtEdge,&
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
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
    end do
    
  end subroutine mhd_trafoDiffDenPre2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDenPreVel2d_sim(DdataAtEdge,&
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

  end subroutine mhd_trafoFluxDenPreVel2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDenPreVel2d_sim(DdataAtEdge,&
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
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
    end do

  end subroutine mhd_trafoDiffDenPreVel2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcBoundaryvalues2d(DbdrNormal, DpointNormal,&
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

  

  end subroutine mhd_calcBoundaryvalues2d

  ! ***************************************************************************

!<subroutine>

  subroutine mhd_coeffVectorBdr2d_sim(rdiscretisation, rform,&
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
    
  end subroutine mhd_coeffVectorBdr2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine mhd_hadaptCallbackScalar2d(iOperation, rcollection)

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
            OU_CLASS_WARNING,OU_MODE_STD,'mhd_hadaptCallbackScalar2d')
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

  end subroutine mhd_hadaptCallbackScalar2d

  !*****************************************************************************

!<subroutine>

  subroutine mhd_hadaptCallbackBlock2d(iOperation, rcollection)

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
            OU_CLASS_WARNING,OU_MODE_STD,'mhd_hadaptCallbackBlock2d')
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

  end subroutine mhd_hadaptCallbackBlock2d

end module mhd_callback2d
