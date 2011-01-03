
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
!# 20.) mhd_calcFluxFCTScDiss2d_sim
!#      -> Computes fluxes for FCT algorithm
!#         adopting scalar artificial viscosities
!#
!# 21.) mhd_calcFluxFCTRoeDiss2d_sim
!#      -> Computes fluxes for FCT algorithm
!#         adopting tensorial artificial viscosities
!#
!# 22.) mhd_calcFluxFCTRusDiss2d_sim
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
!# 28.) mhd_trafoDiffPressure2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the presure
!#
!# 29.) mhd_trafoFluxVelocity2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 30.) mhd_trafoDiffVelocity2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 31.) mhd_trafoFluxMomentum2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 32.) mhd_trafoDiffMomentum2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 33.) mhd_trafoFluxDenEng2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 34.) mhd_trafoDiffDenEng2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 35.) mhd_trafoFluxDenPre2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 36.) mhd_trafoDiffDenPre2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 37.) mhd_trafoFluxDenPreVel2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 38.) mhd_trafoDiffDenPreVel2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure 
!#         and the velocity
!#
!# 39.) mhd_trafoDiffMagfield2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the magnetic field
!#
!# 40.) mhd_trafoFluxMagfield2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the magnetic field
!#
!# 41.) mhd_calcBoundaryvalues2d
!#      -> Computes the boundary values for a given node
!#
!# 42.) mhd_hadaptCallbackScalar2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in interleave format
!#
!# 43.) mhd_hadaptCallbackBlock2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in block format
!#
!# 44.) mhd_coeffVectorBdr2d_sim
!#      -> Calculates the coefficients for the linear form in 2D
!#
!# 45.) mhd_coeffMatrixBdr2d_sim
!#      -> Calculates the coefficients for the bilinear form in 2D
!# </purpose>
!##############################################################################

module mhd_callback2d

#include "mhd.h"

  use boundarycondaux
  use collection
  use derivatives
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
  public :: mhd_calcFluxFCTScDiss2d_sim
  public :: mhd_calcFluxFCTRoeDiss2d_sim
  public :: mhd_calcFluxFCTRusDiss2d_sim
  public :: mhd_trafoFluxDensity2d_sim
  public :: mhd_trafoFluxEnergy2d_sim
  public :: mhd_trafoFluxPressure2d_sim
  public :: mhd_trafoFluxVelocity2d_sim
  public :: mhd_trafoFluxMomentum2d_sim
  public :: mhd_trafoFluxDenEng2d_sim
  public :: mhd_trafoFluxDenPre2d_sim
  public :: mhd_trafoFluxDenPreVel2d_sim
  public :: mhd_trafoFluxMagfield2d_sim
  public :: mhd_trafoDiffDensity2d_sim
  public :: mhd_trafoDiffEnergy2d_sim
  public :: mhd_trafoDiffPressure2d_sim
  public :: mhd_trafoDiffVelocity2d_sim
  public :: mhd_trafoDiffMomentum2d_sim
  public :: mhd_trafoDiffDenEng2d_sim
  public :: mhd_trafoDiffDenPre2d_sim
  public :: mhd_trafoDiffDenPreVel2d_sim
  public :: mhd_trafoDiffMagfield2d_sim
  public :: mhd_calcBoundaryvalues2d
  public :: mhd_coeffVectorBdr2d_sim
  public :: mhd_hadaptCallbackScalar2d
  public :: mhd_hadaptCallbackBlock2d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGal2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the standard Galerkin
    ! discretisation in 2D.
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

  ! Number of edges
  integer, intent(in) :: nedges
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
    real(DP), dimension(NVAR2D) :: Fxi,Fxj,Fyi,Fyj
#else
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
#endif
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx

    
    do idx = 1, nedges

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
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
      FLUX_MHD_2T_XDIR_2D(Fxi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_2D(Fxj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      FLUX_MHD_2T_YDIR_2D(Fyi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_YDIR_2D(Fyj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
      
      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*Fyj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*Fyi )
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_2D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_2D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
           
      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij)
#endif
    end do

  end subroutine mhd_calcFluxGal2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGalNoBdr2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the TVD discretisation
    ! in 2D. The symmetric boundary contributions are neglected and
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

  ! Number of edges
  integer, intent(in) :: nedges
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
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
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
      FLUXDIFF_MHD_2T_XDIR_2D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_2D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
      
      ! Assemble symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale *&
          (0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))*Fx_ij+&
           0.5*(DmatrixCoeffsAtEdge(2,1,idx)-DmatrixCoeffsAtEdge(2,2,idx))*Fy_ij)
      DfluxesAtEdge(:,2,idx) = DfluxesAtEdge(:,1,idx)
    end do

  end subroutine mhd_calcFluxGalNoBdr2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxScDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 2D using scalar artificial viscosities proportional to the
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

  ! Number of edges
  integer, intent(in) :: nedges
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
    real(DP), dimension(NVAR2D) :: Fxi,Fxj,Fyi,Fyj
#else
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx

    
    do idx = 1, nedges

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
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
      FLUX_MHD_2T_XDIR_2D(Fxi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_2D(Fxj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      FLUX_MHD_2T_YDIR_2D(Fyi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_YDIR_2D(Fyj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_2D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_2D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
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
                                         DmatrixCoeffsAtEdge(2,2,idx)*Fyj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*Fyi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij + Diff)
#endif
    end do

  end subroutine mhd_calcFluxScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxScDissDiSp2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)
    

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 2D using scalar artificial viscosities proportional to the
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

  ! Number of edges
  integer, intent(in) :: nedges
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
    real(DP), dimension(NVAR2D) :: Fxi,Fxj,Fyi,Fyj
#else
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: anorm
    integer :: idx
    
    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
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
      FLUX_MHD_2T_XDIR_2D(Fxi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_2D(Fxj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      FLUX_MHD_2T_YDIR_2D(Fyi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_YDIR_2D(Fyj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_2D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_2D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
#endif
      
      !-------------------------------------------------------------------------
      ! Evaluate the scalar artificial dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      Diff = 0

      !!! TODO !!!

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*Fyj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*Fyi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij + Diff)
#endif
    end do

  end subroutine mhd_calcFluxScDissDiSp2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRoeDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 2D using tensorial artificial viscosities of Roe-type.
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

  ! Number of edges
  integer, intent(in) :: nedges
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
    real(DP), dimension(NVAR2D) :: Fxi,Fxj,Fyi,Fyj
#else
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: anorm
    integer :: idx


    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
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
      FLUX_MHD_2T_XDIR_2D(Fxi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_2D(Fxj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      FLUX_MHD_2T_YDIR_2D(Fyi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_YDIR_2D(Fyj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_2D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_2D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

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
                                           DmatrixCoeffsAtEdge(2,2,idx)*Fyj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*Fyi + Diff)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij + Diff)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij + Diff)
#endif
      else
        
#ifdef MHD_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*Fyj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*Fyi)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij)
#endif
      end if
    end do

  end subroutine mhd_calcFluxRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRoeDissDiSp2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)
    

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 2D using tensorial artificial viscosities of Roe-type, whereby
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

  ! Number of edges
  integer, intent(in) :: nedges
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
    real(DP), dimension(NVAR2D) :: Fxi,Fxj,Fyi,Fyj
#else
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
#endif
    real(DP), dimension(NVAR2D) :: DiffX,DiffY
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: anorm
    integer :: idx


    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
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
      FLUX_MHD_2T_XDIR_2D(Fxi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_2D(Fxj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      FLUX_MHD_2T_YDIR_2D(Fyi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_YDIR_2D(Fyj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_2D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_2D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      if (anorm .gt. SYS_EPSREAL) then

        DiffX = 0; DiffY = 0

        !!! TODO !!!

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef MHD_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*Fyj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*Fyi+&
                                           DiffX+DiffY)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                            DiffX+DiffY)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                            DiffX+DiffY)
#endif
      else

#ifdef MHD_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*Fyj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*Fyi)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij)
#endif
      end if
    end do

  end subroutine mhd_calcFluxRoeDissDiSp2d_sim
 
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRusDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 2D using scalar artificial viscosities of Rusanov-type.
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

  ! Number of edges
  integer, intent(in) :: nedges
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
    real(DP), dimension(NVAR2D) :: Fxi,Fxj,Fyi,Fyj
#else
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: ca1i,ca1j,ca2i,ca2j,cf1i,cf1j,cf2i,cf2j,d_ij
    real(DP) :: aPow2i,aPow2j,astPow2i,astPow2j
    integer :: idx


    do idx = 1, nedges

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
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
      FLUX_MHD_2T_XDIR_2D(Fxi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_2D(Fxj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      FLUX_MHD_2T_YDIR_2D(Fyi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_YDIR_2D(Fyj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_2D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_2D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
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
      aPow2i = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      aPow2j = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
#else
#error "Speed of sound must be implemented!"
#endif

      ! Compute auxiliary quantities
      astPow2i = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)/\
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + aPow2i
      astPow2j = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)/\
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + aPow2j

      ! Compute the speed of the fast waves in x-direction
      cf1i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca1i**2)))
      cf1j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca1j**2)))

      ! Compute the speed of the fast waves in y-direction
      cf2i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca2i**2)))
      cf2j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca2j**2)))
      
      ! Compute scalar dissipation for the Rusanov flux
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
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*Fyj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*Fyi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij + Diff)
#endif
    end do

  end subroutine mhd_calcFluxRusDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRusDissDiSp2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 2D using scalar artificial viscosities of Rusanov-type, whereby
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

  ! Number of edges
  integer, intent(in) :: nedges
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
    real(DP), dimension(NVAR2D) :: Fxi,Fxj,Fyi,Fyj
#else
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: ca1i,ca1j,ca2i,ca2j,cf1i,cf1j,cf2i,cf2j,d_ij
    real(DP) :: aPow2i,aPow2j,astPow2i,astPow2j
    integer :: idx

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
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
      FLUX_MHD_2T_XDIR_2D(Fxi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_2D(Fxj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      FLUX_MHD_2T_YDIR_2D(Fyi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_YDIR_2D(Fyj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_2D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_MHD_2T_YDIR_2D(Fy_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
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
      aPow2i = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      aPow2j = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
#else
#error "Speed of sound must be implemented!"
#endif

      ! Compute auxiliary quantities
      astPow2i = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)/\
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + aPow2i
      astPow2j = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)/\
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + aPow2j

      ! Compute the speed of the fast waves in x-direction
      cf1i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca1i**2)))
      cf1j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca1j**2)))

      ! Compute the speed of the fast waves in y-direction
      cf2i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca2i**2)))
      cf2j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca2j**2)))

      ! Compute scalar dissipation with dimensional splitting
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
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*Fyj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*Fyi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij + Diff)
#endif
    end do

  end subroutine mhd_calcFluxRusDissDiSp2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiagMatD2d_sim(DdataAtNode, DmatrixCoeffsAtNode,&
      IverticesAtNode, dscale, nnodes, DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 2D.
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

  ! Number of nodes
  integer, intent(in) :: nnodes
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


  end subroutine mhd_calcMatDiagMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiag2d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 2D.
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

  ! Number of nodes
  integer, intent(in) :: nnodes
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


  end subroutine mhd_calcMatDiag2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGalMatD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 2D.
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

  ! Number of edges
  integer, intent(in) :: nedges
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
  

  end subroutine mhd_calcMatGalMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGal2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 2D.
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

  ! Number of edges
  integer, intent(in) :: nedges
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

      
  end subroutine mhd_calcMatGal2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDissMatD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 2D and applies scalar artificial viscosities proportional to
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

  ! Number of edges
  integer, intent(in) :: nedges
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


  end subroutine mhd_calcMatScDissMatD2d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDiss2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 2D and applies
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

  ! Number of edges
  integer, intent(in) :: nedges
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

   
  end subroutine mhd_calcMatScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDissMatD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 2D and applies
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

  ! Number of edges
  integer, intent(in) :: nedges
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

  
  end subroutine mhd_calcMatRoeDissMatD2d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDiss2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 2D and applies
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

  ! Number of edges
  integer, intent(in) :: nedges
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


  end subroutine mhd_calcMatRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDissMatD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 2D and applies the scalar artificial viscosities of Rusanov-type.
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

  ! Number of edges
  integer, intent(in) :: nedges
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


  end subroutine mhd_calcMatRusDissMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDiss2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 2D applies
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

  ! Number of edges
  integer, intent(in) :: nedges
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

  
  end subroutine mhd_calcMatRusDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcCharacteristics2d_sim(Dweight, DdataAtEdge,&
      nedges, DcharVariablesAtEdge, DeigenvaluesAtEdge,&
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

    ! Number of edges
    integer, intent(in) :: nedges
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
    

  end subroutine mhd_calcCharacteristics2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTScDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for FCT
    ! algorithms in 2D using scalar dissipation proportional to the
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

  ! Number of edges
  integer, intent(in) :: nedges
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
  

  end subroutine mhd_calcFluxFCTScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTRoeDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes forFCT
    ! algorithms in 2D using tensorial dissipation of Roe-type.
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

  ! Number of edges
  integer, intent(in) :: nedges
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


  end subroutine mhd_calcFluxFCTRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTRusDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for FCT
    ! algorithms in 2D using scalar dissipation of Rusanov-type.
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

  ! Number of edges
  integer, intent(in) :: nedges
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
    real(DP) :: ca1i,ca1j,ca2i,ca2j,cf1i,cf1j,cf2i,cf2j,d_ij
    real(DP) :: aPow2i,aPow2j,astPow2i,astPow2j
    integer :: idx
    

    do idx = 1, nedges

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Compute the speed of the Alfven waves in x-direction
      ca1i = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx))
      ca1j = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))

      ! Compute the speed of the Alfven waves in y-direction
      ca2i = abs(Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx))
      ca2j = abs(Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx))

! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      aPow2i = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      aPow2j = GAMMA*PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
#else
#error "Speed of sound must be implemented!"
#endif

      ! Compute auxiliary quantities
      astPow2i = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)/\
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx) + aPow2i
      astPow2j = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)/\
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx) + aPow2j

      ! Compute the speed of the fast waves in x-direction
      cf1i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca1i**2)))
      cf1j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca1j**2)))

      ! Compute the speed of the fast waves in y-direction
      cf2i = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*ca2i**2)))
      cf2j = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*ca2j**2)))
      
      ! Compute scalar dissipation for the Rusanov flux
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj+&
                      DmatrixCoeffsAtEdge(2,1,idx)*vj)+&
                 sqrt((DmatrixCoeffsAtEdge(1,1,idx)**2)*cf1j+&
                      (DmatrixCoeffsAtEdge(2,1,idx)**2)*cf2j),&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui+&
                      DmatrixCoeffsAtEdge(2,2,idx)*vi)+&
                 sqrt((DmatrixCoeffsAtEdge(1,2,idx)**2)*cf1i+&
                      (DmatrixCoeffsAtEdge(2,2,idx)**2)*cf2i) )

      ! Compute conservative fluxes
      DfluxesAtEdge(:,idx) = dscale*d_ij*(DdataAtEdge(:,1,idx)-DdataAtEdge(:,2,idx))
    end do

  end subroutine mhd_calcFluxFCTRusDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDensity2d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative variables to fluxes for the density in 2D
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

    ! Number of edges
    integer, intent(in) :: nedges
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
    
    do idx = 1, nedges
      
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
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density in 2D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

    ! Number of edges
    integer, intent(in) :: nedges
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
    
    do idx = 1, nedges
      
      ! Transformed density difference
      DtransformedDataAtEdge(1,idx) =&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
    end do

  end subroutine mhd_trafoDiffDensity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxEnergy2d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the energy in 2D
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

    ! Number of edges
    integer, intent(in) :: nedges
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
    
    do idx = 1, nedges
      
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
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the energy in 2D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

    ! Number of edges
    integer, intent(in) :: nedges
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
    
    do idx = 1, nedges
      
      ! Transformed total density difference
      DtransformedDataAtEdge(1,idx) =&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
    end do

  end subroutine mhd_trafoDiffEnergy2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxPressure2d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the pressure in 2D
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

    ! Number of edges
    integer, intent(in) :: nedges
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

    do idx = 1, nedges
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      
      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(1,1,idx) = (GAMMA-1.0)*&
          (0.5*(ui*ui+vi*vi+wi*wi)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vi*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          wi*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
          X_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
          Y_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
          Z_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))
      DtransformedFluxesAtEdge(1,2,idx) =-(GAMMA-1.0)*&
          (0.5*(uj*uj+vj*vj+wj*wj)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vj*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          wj*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
          X_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
          Y_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
          Z_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do

  end subroutine mhd_trafoFluxPressure2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffPressure2d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the pressure in 2D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

    ! Number of edges
    integer, intent(in) :: nedges
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

    do idx = 1, nedges
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(1,idx) =&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
    end do

  end subroutine mhd_trafoDiffPressure2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxVelocity2d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the velocity in 2D
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

    ! Number of edges
    integer, intent(in) :: nedges
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

    do idx = 1, nedges
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)


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

      ! Transformed velocity fluxes in z-direction
      DtransformedFluxesAtEdge(3,1,idx) =&
          (Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
           wi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      DtransformedFluxesAtEdge(3,2,idx) =&
         -(Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
           wj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
    end do
    
  end subroutine mhd_trafoFluxVelocity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffVelocity2d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the velocity in 2D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

    ! Number of edges
    integer, intent(in) :: nedges
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
    
    do idx = 1, nedges

      ! Transformed velocity difference in x-direction
      DtransformedDataAtEdge(1,idx) =&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)

      ! Transformed velocity difference in y-direction
      DtransformedDataAtEdge(2,idx) =&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)

      ! Transformed velocity difference in z-direction
      DtransformedDataAtEdge(3,idx) =&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
    end do
    
  end subroutine mhd_trafoDiffVelocity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxMomentum2d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the momentum in 2D
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

    ! Number of edges
    integer, intent(in) :: nedges
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
    
    do idx = 1, nedges
      
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

      ! Transformed momentum fluxes in z-direction
      DtransformedFluxesAtEdge(3,1,idx) =&
          Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
      DtransformedFluxesAtEdge(3,2,idx) =&
         -Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
    end do
    
  end subroutine mhd_trafoFluxMomentum2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffMomentum2d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the momentum in 2D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

    ! Number of edges
    integer, intent(in) :: nedges
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
    
    do idx = 1, nedges
      
      ! Transformed momentum difference in x-direction
      DtransformedDataAtEdge(1,idx) =&
          X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)

      ! Transformed momentum difference in y-direction
      DtransformedDataAtEdge(2,idx) =&
          Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)

      ! Transformed momentum difference in z-direction
      DtransformedDataAtEdge(3,idx) =&
          Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
    end do
    
  end subroutine mhd_trafoDiffMomentum2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDenEng2d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density and energy in 2D
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

    ! Number of edges
    integer, intent(in) :: nedges
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
    
    do idx = 1, nedges
      
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
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density and energy in 2D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

    ! Number of edges
    integer, intent(in) :: nedges
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
    
    do idx = 1, nedges

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
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density and energy in 2D
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

    ! Number of edges
    integer, intent(in) :: nedges
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

    do idx = 1, nedges
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
           
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
      
      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
     DtransformedFluxesAtEdge(2,1,idx) = (GAMMA-1.0)*&
          (0.5*(ui*ui+vi*vi+wi*wi)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vi*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          wi*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
          X_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
          Y_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
          Z_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))
      DtransformedFluxesAtEdge(2,2,idx) =-(GAMMA-1.0)*&
          (0.5*(uj*uj+vj*vj+wj*wj)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vj*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          wj*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
          X_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
          Y_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
          Z_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do

  end subroutine mhd_trafoFluxDenPre2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDenPre2d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density and energy in 2D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

    ! Number of edges
    integer, intent(in) :: nedges
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

    do idx = 1, nedges
      
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
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density, pressure 
    ! and velocity in 2D
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

    ! Number of edges
    integer, intent(in) :: nedges
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

    do idx = 1, nedges
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)
      
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

      ! Transformed velocity fluxes in z-direction
      DtransformedFluxesAtEdge(4,1,idx) =&
          (Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
           wi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      DtransformedFluxesAtEdge(4,2,idx) =&
         -(Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
           wj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)

      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(5,1,idx) = (GAMMA-1.0)*&
          (0.5*(ui*ui+vi*vi+wi*wi)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vi*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          wi*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
          X_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
          Y_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)*&
          Z_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))
      DtransformedFluxesAtEdge(5,2,idx) =-(GAMMA-1.0)*&
          (0.5*(uj*uj+vj*vj+wj*wj)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          vj*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          wj*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
          X_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
          Y_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)*&
          Z_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)-&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do

  end subroutine mhd_trafoFluxDenPreVel2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine mhd_trafoDiffDenPreVel2d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density,
    ! pressure and velocity in 2D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Number of edges
    integer, intent(in) :: nedges
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

    do idx = 1, nedges
      
      ! Transformed density difference
      DtransformedDataAtEdge(1,idx) =&
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

      ! Transformed velocity difference in z-direction
      DtransformedDataAtEdge(4,idx) =&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(5,idx) =&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
    end do

  end subroutine mhd_trafoDiffDenPreVel2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxMagfield2d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the magnetic field in 2D
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

    ! Number of edges
    integer, intent(in) :: nedges
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

    do idx = 1, nedges

      ! Transformed magnetic field fluxes in x-direction
      DtransformedFluxesAtEdge(1,1,idx) =&
          X_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -X_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)

      ! Transformed magnetic field fluxes in y-direction
      DtransformedFluxesAtEdge(2,1,idx) =&
          Y_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -Y_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)

      ! Transformed magnetic field fluxes in z-direction
      DtransformedFluxesAtEdge(3,1,idx) =&
          Z_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
      DtransformedFluxesAtEdge(3,2,idx) =&
         -Z_MAGNETICFIELD_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR2D,idx)
    end do

  end subroutine mhd_trafoFluxMagfield2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffMagfield2d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative convervative to differences for the magnetic 
    ! field in 2D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

    ! Number of edges
    integer, intent(in) :: nedges
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

    do idx = 1, nedges
      
      ! Transformed magnetic field difference in x-direction
      DtransformedDataAtEdge(1,idx) =&
          X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          X_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)

      ! Transformed magnetic field difference in y-direction
      DtransformedDataAtEdge(2,idx) =&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)

      ! Transformed magnetic field difference in z-direction
      DtransformedDataAtEdge(3,idx) =&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,2,idx)-&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR(DdataAtEdge,NVAR2D,1,idx)
    end do

  end subroutine mhd_trafoDiffMagfield2d_sim

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

    ! An array accepting the DOF`s on all elements in the test space.
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