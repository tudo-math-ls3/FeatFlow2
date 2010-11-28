!##############################################################################
!# ****************************************************************************
!# <name> mhd_callback3d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible MHDequations in 3D.
!#
!# The following callback functions are available:
!#
!# 1.) mhd_calcFluxGal3d_sim
!#     -> Computes fluxes for standard Galerkin scheme
!#
!# 2.) mhd_calcFluxGalNoBdr3d_sim
!#     -> Computes fluxes for standard Galerkin scheme without
!#        assembling the symmetric boundary contribution
!#
!# 3.) mhd_calcFluxScDiss3d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting scalar artificial viscosities
!#
!# 4.) mhd_calcFluxScDissDiSp3d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting scalar artificial viscosities based on
!#        dimensional splitting approach
!#
!# 5.) mhd_calcFluxRoeDiss3d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities
!#
!# 6.) mhd_calcFluxRoeDissDiSp3d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities based on
!#        dimensional splitting approach
!#
!# 7.) mhd_calcFluxRusDiss3d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion
!#
!# 8.) mhd_calcFluxRusDissDiSp3d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion
!#
!# 9.) mhd_calcMatDiagMatD3d_sim
!#     -> Computes local matrix for diagonal entry
!#
!# 10.) mhd_calcMatDiag3d_sim
!#      -> Computes local matrix for diagonal entry
!#
!# 11.) mhd_calcMatGalMatD3d_sim
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 12.) mhd_calcMatGal3d_sim
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 13.) mhd_calcMatScDissMatD3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 14.) mhd_calcMatScDiss3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 15.) mhd_calcMatRoeDissMatD3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 16.) mhd_calcMatRoeDiss3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 17.) mhd_calcMatRusDissMatD3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov artificial viscosities
!#
!# 18.) mhd_calcMatRusDiss3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov flux artificial viscosities
!#
!# 19.) mhd_calcCharacteristics3d_sim
!#      -> Computes characteristic variables
!#
!# 20.) mhd_calcFluxFCTScDiss3d_sim
!#      -> Computes fluxes for FCT algorithm
!#         adopting scalar artificial viscosities
!#
!# 21.) mhd_calcFluxFCTRoeDiss3d_sim
!#      -> Computes fluxes for FCT algorithm
!#         adopting tensorial artificial viscosities
!#
!# 22.) mhd_calcFluxFCTRusDiss3d_sim
!#      -> Computes fluxes for FCT algorithm
!#         adopting the Rusanov artificial viscosities
!#
!# 23.) mhd_trafoFluxDensity3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 24.) mhd_trafoDiffDensity3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 25.) mhd_trafoFluxEnergy3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 26.) mhd_trafoDiffEnergy3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 27.) mhd_trafoFluxPressure3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 28.) mhd_trafoDiffPressure3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the pressure
!#
!# 29.) mhd_trafoFluxVelocity3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 30.) mhd_trafoDiffVelocity3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 31.) mhd_trafoFluxMomentum3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 32.) mhd_trafoDiffMomentum3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 33.) mhd_trafoFluxDenEng3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 34.) mhd_trafoDiffDenEng3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 35.) mhd_trafoFluxDenPre3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 36.) mhd_trafoDiffDenPre3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 37.) mhd_trafoFluxDenPreVel3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 38.) mhd_trafoDiffDenPreVel3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure 
!#         and the velocity
!#
!# 39.) mhd_trafoDiffMagfield3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the magnetic field
!#
!# 40.) mhd_trafoFluxMagfield3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the magnetic field
!#
!# 41.) mhd_calcBoundaryvalues3d
!#      -> Computes the boundary values for a given node
!#
!# 42.) mhd_hadaptCallbackScalar3d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 3D, whereby the vector is stored in interleave format
!#
!# 43.) mhd_hadaptCallbackBlock3d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 3D, whereby the vector is stored in block format
!#
!# </purpose>
!##############################################################################

module mhd_callback3d

#include "mhd.h"

  use collection
  use flagship_callback
  use fsystem
  use genoutput
  use graph
  use groupfemsystem
  use hadaptaux
  use linearsystemblock
  use linearsystemscalar
  use mhd_basic
  use problem
  use solveraux
  use storage

  ! The following 3D-routine coincide with their 2D-versions
  use mhd_callback2d, only : &
      mhd_trafoFluxDensity3d_sim => mhd_trafoFluxDensity2d_sim
  use mhd_callback2d, only : &
      mhd_trafoDiffDensity3d_sim => mhd_trafoDiffDensity2d_sim
  use mhd_callback2d, only : &
      mhd_trafoFluxEnergy3d_sim => mhd_trafoFluxEnergy2d_sim
  use mhd_callback2d, only : &
      mhd_trafoDiffEnergy3d_sim => mhd_trafoDiffEnergy2d_sim
  use mhd_callback2d, only : &
      mhd_trafoFluxPressure3d_sim => mhd_trafoFluxPressure2d_sim
  use mhd_callback2d, only : &
      mhd_trafoDiffPressure3d_sim => mhd_trafoDiffPressure2d_sim
  use mhd_callback2d, only : &
      mhd_trafoFluxVelocity3d_sim => mhd_trafoFluxVelocity2d_sim
  use mhd_callback2d, only : &
      mhd_trafoDiffVelocity3d_sim => mhd_trafoDiffVelocity2d_sim
  use mhd_callback2d, only : &
      mhd_trafoFluxMomentum3d_sim => mhd_trafoFluxMomentum2d_sim
  use mhd_callback2d, only : &
      mhd_trafoDiffMomentum3d_sim => mhd_trafoDiffMomentum2d_sim
  use mhd_callback2d, only : &
      mhd_trafoFluxDenEng3d_sim => mhd_trafoFluxDenEng2d_sim
  use mhd_callback2d, only : &
      mhd_trafoDiffDenEng3d_sim => mhd_trafoDiffDenEng2d_sim
  use mhd_callback2d, only : &
      mhd_trafoFluxDenPre3d_sim => mhd_trafoFluxDenPre2d_sim
  use mhd_callback2d, only : &
      mhd_trafoDiffDenPre3d_sim => mhd_trafoDiffDenPre2d_sim
  use mhd_callback2d, only : &
      mhd_trafoFluxDenPreVel3d_sim => mhd_trafoFluxDenPreVel2d_sim
  use mhd_callback2d, only : &
      mhd_trafoDiffDenPreVel3d_sim => mhd_trafoDiffDenPreVel2d_sim
  use mhd_callback2d, only : &
      mhd_trafoDiffMagfield3d_sim => mhd_trafoDiffMagfield2d_sim
  use mhd_callback2d, only : &
      mhd_trafoFluxMagfield3d_sim => mhd_trafoFluxMagfield2d_sim

  implicit none

  private
  public :: mhd_calcFluxGal3d_sim
  public :: mhd_calcFluxGalNoBdr3d_sim
  public :: mhd_calcFluxScDiss3d_sim
  public :: mhd_calcFluxScDissDiSp3d_sim
  public :: mhd_calcFluxRoeDiss3d_sim
  public :: mhd_calcFluxRoeDissDiSp3d_sim
  public :: mhd_calcFluxRusDiss3d_sim
  public :: mhd_calcFluxRusDissDiSp3d_sim
  public :: mhd_calcMatDiagMatD3d_sim
  public :: mhd_calcMatDiag3d_sim
  public :: mhd_calcMatGalMatD3d_sim
  public :: mhd_calcMatGal3d_sim
  public :: mhd_calcMatScDissMatD3d_sim
  public :: mhd_calcMatScDiss3d_sim
  public :: mhd_calcMatRoeDissMatD3d_sim
  public :: mhd_calcMatRoeDiss3d_sim
  public :: mhd_calcMatRusDissMatD3d_sim
  public :: mhd_calcMatRusDiss3d_sim
  public :: mhd_calcCharacteristics3d_sim
  public :: mhd_calcFluxFCTScDiss3d_sim
  public :: mhd_calcFluxFCTRoeDiss3d_sim
  public :: mhd_calcFluxFCTRusDiss3d_sim
  public :: mhd_trafoFluxDensity3d_sim
  public :: mhd_trafoFluxEnergy3d_sim
  public :: mhd_trafoFluxPressure3d_sim
  public :: mhd_trafoFluxVelocity3d_sim
  public :: mhd_trafoFluxMomentum3d_sim
  public :: mhd_trafoFluxDenEng3d_sim
  public :: mhd_trafoFluxDenPre3d_sim
  public :: mhd_trafoFluxDenPreVel3d_sim
  public :: mhd_trafoFluxMagfield3d_sim
  public :: mhd_trafoDiffDensity3d_sim
  public :: mhd_trafoDiffEnergy3d_sim
  public :: mhd_trafoDiffPressure3d_sim
  public :: mhd_trafoDiffVelocity3d_sim
  public :: mhd_trafoDiffMomentum3d_sim
  public :: mhd_trafoDiffDenEng3d_sim
  public :: mhd_trafoDiffDenPre3d_sim
  public :: mhd_trafoDiffDenPreVel3d_sim
  public :: mhd_trafoDiffMagfield3d_sim
  public :: mhd_calcBoundaryvalues3d
  public :: mhd_hadaptCallbackScalar3d
  public :: mhd_hadaptCallbackBlock3d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGal3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the standard Galerkin
    ! discretisation in 3D.
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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_MHD_2T_XDIR_3D(Fxi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_3D(Fxj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      FLUX_MHD_2T_YDIR_3D(Fyi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_YDIR_3D(Fyj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      FLUX_MHD_2T_ZDIR_3D(Fzi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_ZDIR_3D(Fzj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                         DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                         DmatrixCoeffsAtEdge(3,1,idx)*Fzi )
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)

      ! Compute flux difference for z-direction
      FLUXDIFF_MHD_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)

      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij)
#endif
    end do

  end subroutine mhd_calcFluxGal3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGalNoBdr3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the TVD discretisation
    ! in 3D. The symmetric boundary contributions are neglected and
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
  real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
  real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
  integer :: idx
  
    
    do idx = 1, size(DfluxesAtEdge,3)

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj

      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)

      ! Compute flux difference for z-direction
      FLUXDIFF_MHD_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)

      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale *&
          (0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))*Fx_ij+&
           0.5*(DmatrixCoeffsAtEdge(2,1,idx)-DmatrixCoeffsAtEdge(2,2,idx))*Fy_ij+&
           0.5*(DmatrixCoeffsAtEdge(3,1,idx)-DmatrixCoeffsAtEdge(3,2,idx))*Fz_ij)
      DfluxesAtEdge(:,2,idx) = DfluxesAtEdge(:,1,idx)
    end do

  end subroutine mhd_calcFluxGalNoBdr3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxScDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 3D using scalar artificial viscosities proportional to the
    ! spectral radius (largest eigenvalue) of the Roe-matrix.
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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_MHD_2T_XDIR_3D(Fxi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_3D(Fxj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      FLUX_MHD_2T_YDIR_3D(Fyi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_YDIR_3D(Fyj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      FLUX_MHD_2T_ZDIR_3D(Fzi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_ZDIR_3D(Fzj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)

      ! Compute flux difference for z-direction
      FLUXDIFF_MHD_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar artificial dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------

      Diff = 0

      !!! TODO !!!

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                         DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                         DmatrixCoeffsAtEdge(3,1,idx)*Fzi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij+ Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij+ Diff)
#endif
    end do

  end subroutine mhd_calcFluxScDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxScDissDiSp3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 3D using scalar artificial viscosities proportional to the
    ! spectral radius (largest eigenvalue) of the Roe-matrix, whereby
    ! dimensional splitting is employed.
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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: anorm
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_MHD_2T_XDIR_3D(Fxi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_3D(Fxj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      FLUX_MHD_2T_YDIR_3D(Fyi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_YDIR_3D(Fyj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      FLUX_MHD_2T_ZDIR_3D(Fzi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_ZDIR_3D(Fzj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)

      ! Compute flux difference for z-direction
      FLUXDIFF_MHD_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar artificial dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      Diff = 0

      !!! TODO !!!

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                         DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                         DmatrixCoeffsAtEdge(3,1,idx)*Fzi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij + Diff)
#endif
    end do

  end subroutine mhd_calcFluxScDissDiSp3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRoeDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 3D using tensorial artificial viscosities of Roe-type.
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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: anorm
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_MHD_2T_XDIR_3D(Fxi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_3D(Fxj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      FLUX_MHD_2T_YDIR_3D(Fyi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_YDIR_3D(Fyj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      FLUX_MHD_2T_ZDIR_3D(Fzi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_ZDIR_3D(Fzj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)

      ! Compute flux difference for z-direction
      FLUXDIFF_MHD_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      if (anorm .gt. SYS_EPSREAL) then
        
        ! Normalize the skew-symmetric coefficient
        a = a/anorm
        
        Diff = 0

        !! TODO !!

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef MHD_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                           DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                           DmatrixCoeffsAtEdge(3,1,idx)*Fzi + Diff)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij + Diff)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij + Diff)
#endif
      else
        
#ifdef MHD_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                           DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                           DmatrixCoeffsAtEdge(3,1,idx)*Fzi)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij)
#endif
      end if
    end do

  end subroutine mhd_calcFluxRoeDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRoeDissDiSp3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 3D using tensorial artificial viscosities of Roe-type, whereby
    ! dimensional splitting is employed.
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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: DiffX,DiffY,DiffZ
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: anorm
    integer :: idx
    

    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_MHD_2T_XDIR_3D(Fxi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_3D(Fxj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      FLUX_MHD_2T_YDIR_3D(Fyi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_YDIR_3D(Fyj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      FLUX_MHD_2T_ZDIR_3D(Fzi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_ZDIR_3D(Fzj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)

      ! Compute flux difference for z-direction
      FLUXDIFF_MHD_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
#endif
      
      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      if (anorm .gt. SYS_EPSREAL) then
        
        ! Normalize the skew-symmetric coefficient
        a = a/anorm
        
        DiffX = 0; DiffY = 0; DiffZ = 0

        !!! TODO !!!

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef MHD_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                           DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                           DmatrixCoeffsAtEdge(3,1,idx)*Fzi+&
                                           DiffX+DiffY+DiffZ)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij+&
                                            DiffX+DiffY+DiffZ)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij+&
                                            DiffX+DiffY+DiffZ)
#endif
      else
        
#ifdef MHD_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                           DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                           DmatrixCoeffsAtEdge(3,1,idx)*Fzi)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij)
#endif
      end if
    end do

  end subroutine mhd_calcFluxRoeDissDiSp3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRusDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 3D using scalar artificial viscosities of Rusanov-type.
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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: ca1i,ca1j,ca2i,ca2j,ca3i,ca3j,cf1i,cf1j,cf2i,cf2j,cf3i,cf3j,d_ij
    real(DP) :: aPow2i,aPow2j,astPow2i,astPow2j
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes      
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_MHD_2T_XDIR_3D(Fxi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_3D(Fxj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      FLUX_MHD_2T_YDIR_3D(Fyi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_YDIR_3D(Fyj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      FLUX_MHD_2T_ZDIR_3D(Fzi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_ZDIR_3D(Fzj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)

      ! Compute flux difference for z-direction
      FLUXDIFF_MHD_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
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
      ca1i = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx))
      ca1j = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))

      ! Compute the speed of the Alfven waves in y-direction
      ca2i = abs(Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx))
      ca2j = abs(Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))

      ! Compute the speed of the Alfven waves in z-direction
      ca3i = abs(Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx))
      ca3j = abs(Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
      
! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      aPow2i = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      aPow2j = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
#else
#error "Speed of sound must be implemented!"
#endif

      ! Compute auxiliary quantities
      astPow2i = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)/\
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx) + aPow2i
      astPow2j = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)/\
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx) + aPow2j

      ! Compute the speed of the fast waves in x-direction
      cf1i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca1i**2)))
      cf1j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca1j**2)))

      ! Compute the speed of the fast waves in y-direction
      cf2i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca2i**2)))
      cf2j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca2j**2)))

      ! Compute the speed of the fast waves in z-direction
      cf3i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca3i**2)))
      cf3j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca3j**2)))

      ! Compute scalar dissipation for the Rusanov flux
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj+&
                      DmatrixCoeffsAtEdge(2,1,idx)*vj+&
                      DmatrixCoeffsAtEdge(3,1,idx)*wj)+&
                 sqrt((DmatrixCoeffsAtEdge(1,1,idx)**2)*cf1j+&
                      (DmatrixCoeffsAtEdge(2,1,idx)**2)*cf2j+&
                      (DmatrixCoeffsAtEdge(3,1,idx)**2)*cf3j),&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui+&
                      DmatrixCoeffsAtEdge(2,2,idx)*vi+&
                      DmatrixCoeffsAtEdge(3,2,idx)*wi)+&
                 sqrt((DmatrixCoeffsAtEdge(1,2,idx)**2)*cf1i+&
                      (DmatrixCoeffsAtEdge(2,2,idx)**2)*cf2i+&
                      (DmatrixCoeffsAtEdge(3,2,idx)**2)*cf3i) )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------
      
#ifdef MHD_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                         DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                         DmatrixCoeffsAtEdge(3,1,idx)*Fzi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij + Diff)
#endif
    end do

  end subroutine mhd_calcFluxRusDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRusDissDiSp3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 3D using scalar artificial viscosities of Rusanov-type, whereby
    ! dimensional splitting is employed.
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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: ca1i,ca1j,ca2i,ca2j,ca3i,ca3j,cf1i,cf1j,cf2i,cf2j,cf3i,cf3j,d_ij
    real(DP) :: aPow2i,aPow2j,astPow2i,astPow2j
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_MHD_2T_XDIR_3D(Fxi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_3D(Fxj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      FLUX_MHD_2T_YDIR_3D(Fyi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_YDIR_3D(Fyj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      FLUX_MHD_2T_ZDIR_3D(Fzi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_ZDIR_3D(Fzj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)

      ! Compute flux difference for z-direction
      FLUXDIFF_MHD_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
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
      ca1i = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx))
      ca1j = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))

      ! Compute the speed of the Alfven waves in y-direction
      ca2i = abs(Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx))
      ca2j = abs(Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))

      ! Compute the speed of the Alfven waves in z-direction
      ca3i = abs(Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx))
      ca3j = abs(Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))

! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      aPow2i = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      aPow2j = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
#else
#error "Speed of sound must be implemented!"
#endif

      ! Compute auxiliary quantities
      astPow2i = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)/\
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx) + aPow2i
      astPow2j = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)/\
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx) + aPow2j

      ! Compute the speed of the fast waves in x-direction
      cf1i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca1i**2)))
      cf1j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca1j**2)))

      ! Compute the speed of the fast waves in y-direction
      cf2i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca2i**2)))
      cf2j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca2j**2)))

      ! Compute the speed of the fast waves in z-direction
      cf3i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca3i**2)))
      cf3j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca3j**2)))
      
      ! Compute calar dissipation with dimensional splitting
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj)+&
                  abs(DmatrixCoeffsAtEdge(1,1,idx))*cf1j,&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui)+&
                  abs(DmatrixCoeffsAtEdge(1,2,idx))*cf1i )&
           + max( abs(DmatrixCoeffsAtEdge(2,1,idx)*vj)+&
                  abs(DmatrixCoeffsAtEdge(2,1,idx))*cf2j,&
                  abs(DmatrixCoeffsAtEdge(2,2,idx)*vi)+&
                  abs(DmatrixCoeffsAtEdge(2,2,idx))*cf2i )&
           + max( abs(DmatrixCoeffsAtEdge(3,1,idx)*wj)+&
                  abs(DmatrixCoeffsAtEdge(3,1,idx))*cf3j,&
                  abs(DmatrixCoeffsAtEdge(3,2,idx)*wi)+&
                  abs(DmatrixCoeffsAtEdge(3,2,idx))*cf3i )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))
      
      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                         DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                         DmatrixCoeffsAtEdge(3,1,idx)*Fzi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij + Diff)
#endif
    end do

  end subroutine mhd_calcFluxRusDissDiSp3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiagMatD3d_sim(DdataAtNode, DmatrixCoeffsAtNode,&
      IverticesAtNode, dscale, DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 3D.
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


  end subroutine mhd_calcMatDiagMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiag3d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 3D.
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


  end subroutine mhd_calcMatDiag3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGalMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 3D.
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


  end subroutine mhd_calcMatGalMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGal3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D.
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

  end subroutine mhd_calcMatGal3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDissMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 3D and applies scalar artificial viscosities proportional to
    ! the spectral radius (largest eigenvalue) of the Roe-matrix.
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


  end subroutine mhd_calcMatScDissMatD3d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDiss3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D and applies
    ! scalar artificial viscosities proportional to the spectral
    ! radius (largest eigenvalue) of the Roe-matrix.
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


  end subroutine mhd_calcMatScDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDissMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D and applies
    ! tensorial artificial viscosities of Roe-type.
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


  end subroutine mhd_calcMatRoeDissMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDiss3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D and applies
    ! tensorial artificial viscosities of Roe-type.
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


  end subroutine mhd_calcMatRoeDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDissMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 3D and applies the scalar artificial viscosities of Rusanov-type.
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


  end subroutine mhd_calcMatRusDissMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDiss3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D applies
    ! scalar artificial viscosities of Rusanov-type.
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


  end subroutine mhd_calcMatRusDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcCharacteristics3d_sim(Dweight, DdataAtEdge,&
      DcharVariablesAtEdge, DeigenvaluesAtEdge,&
      DrightEigenvectorsAtEdge, DleftEigenvectorsAtEdge, rcollection)

!<description>
    ! This subroutine computes the characteristic variables in 3D.
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


  end subroutine mhd_calcCharacteristics3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTScDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for FCT
    ! algorithms in 3D using scalar dissipation proportional to the
    ! spectral radius of the Roe-matrix.
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
  !   DIMENSION(nvar,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>


  end subroutine mhd_calcFluxFCTScDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTRoeDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes forFCT
    ! algorithms in 3D using tensorial dissipation of Roe-type.
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
  !   DIMENSION(nvar,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

  end subroutine mhd_calcFluxFCTRoeDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTRusDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for FCT
    ! algorithms in 3D using scalar dissipation of Rusanov-type.
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
  !   DIMENSION(nvar,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

    ! local variables
    real(DP) :: ui,uj,vi,vj,wi,wj
    real(DP) :: ca1i,ca1j,ca2i,ca2j,ca3i,ca3j,cf1i,cf1j,cf2i,cf2j,cf3i,cf3j,d_ij
    real(DP) :: aPow2i,aPow2j,astPow2i,astPow2j
    integer :: idx
    

    do idx = 1, size(DfluxesAtEdge,2)

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute the speed of the Alfven waves in x-direction
      ca1i = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx))
      ca1j = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))

      ! Compute the speed of the Alfven waves in y-direction
      ca2i = abs(Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx))
      ca2j = abs(Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))

      ! Compute the speed of the Alfven waves in z-direction
      ca3i = abs(Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx))
      ca3j = abs(Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
      
! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      aPow2i = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      aPow2j = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
#else
#error "Speed of sound must be implemented!"
#endif

      ! Compute auxiliary quantities
      astPow2i = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)/\
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx) + aPow2i
      astPow2j = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)/\
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx) + aPow2j

      ! Compute the speed of the fast waves in x-direction
      cf1i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca1i**2)))
      cf1j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca1j**2)))

      ! Compute the speed of the fast waves in y-direction
      cf2i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca2i**2)))
      cf2j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca2j**2)))

      ! Compute the speed of the fast waves in z-direction
      cf3i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca3i**2)))
      cf3j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca3j**2)))

      ! Compute scalar dissipation for the Rusanov flux
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj+&
                      DmatrixCoeffsAtEdge(2,1,idx)*vj+&
                      DmatrixCoeffsAtEdge(3,1,idx)*wj)+&
                 sqrt((DmatrixCoeffsAtEdge(1,1,idx)**2)*cf1j+&
                      (DmatrixCoeffsAtEdge(2,1,idx)**2)*cf2j+&
                      (DmatrixCoeffsAtEdge(3,1,idx)**2)*cf3j),&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui+&
                      DmatrixCoeffsAtEdge(2,2,idx)*vi+&
                      DmatrixCoeffsAtEdge(3,2,idx)*wi)+&
                 sqrt((DmatrixCoeffsAtEdge(1,2,idx)**2)*cf1i+&
                      (DmatrixCoeffsAtEdge(2,2,idx)**2)*cf2i+&
                      (DmatrixCoeffsAtEdge(3,2,idx)**2)*cf3i) )

      ! Compute conservative fluxes
      DfluxesAtEdge(:,idx) = dscale*d_ij*(DdataAtEdge(:,1,idx)-DdataAtEdge(:,2,idx))
    end do

  end subroutine mhd_calcFluxFCTRusDiss3d_sim
 
  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcBoundaryvalues3d(DbdrNormal, DpointNormal,&
      DbdrValue, ibdrCondType, Du, Du0, istatus)

!<description>
    ! This subroutine computes the boundary values for a given node in 3D.
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

  end subroutine mhd_calcBoundaryvalues3d

  !*****************************************************************************

!<subroutine>

  subroutine mhd_hadaptCallbackScalar3d(iOperation, rcollection)

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
            OU_CLASS_WARNING,OU_MODE_STD,'mhd_hadaptCallbackScalar3d')
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

  end subroutine mhd_hadaptCallbackScalar3d

  !*****************************************************************************

!<subroutine>

  subroutine mhd_hadaptCallbackBlock3d(iOperation, rcollection)

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
            OU_CLASS_WARNING,OU_MODE_STD,'mhd_hadaptCallbackBlock3d')
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

  end subroutine mhd_hadaptCallbackBlock3d

end module mhd_callback3d
