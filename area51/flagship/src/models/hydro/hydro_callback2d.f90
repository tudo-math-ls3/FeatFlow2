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
!#     -> Computes fluxes for standard Galerkin scheme
!#
!# 2.) hydro_calcFluxGalNoBdr2d_sim
!#     -> Computes fluxes for standard Galerkin scheme without
!#        assembling the symmetric boundary contribution
!#
!# 3.) hydro_calcFluxScDiss2d_sim
!#     -> Computes fluxes for low-order discretisation adopting
!#        scalar artificial viscosities
!#
!# 4.) hydro_calcFluxScDissDiSp2d_sim
!#     -> Computes fluxes for low-order discretisation adopting
!#        scalar artificial viscosities based on dimensional
!#        splitting approach
!#
!# 5.) hydro_calcFluxRoeDiss2d_sim
!#     -> Computes fluxes for low-order discretisation adopting
!#       tensorial artificial viscosities of Roe-type
!#
!# 6.) hydro_calcFluxRoeDissDiSp2d_sim
!#     -> Computes fluxes for low-order discretisation adopting
!#        tensorial artificial viscosities of Roe-type based on
!#        dimensional splitting approach
!#
!# 7.) hydro_calcFluxRusDiss2d_sim
!#     -> Computes fluxes for low-order discretisation adopting
!#        scalar artificial diffusion of Rusanov-type
!#
!# 8.) hydro_calcFluxRusDissDiSp2d_sim
!#     -> Computes fluxes for low-order discretisation adopting
!#        scalar artificial diffusion of Rusanov-type based on
!#        dimensional splitting approach
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
!#         adopting tensorial artificial viscosities of Roe-type
!#
!# 16.) hydro_calcMatRoeDiss2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities of Roe-type
!#
!# 17.) hydro_calcMatRusDissMatD2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities of Rusanov-type
!#
!# 18.) hydro_calcMatRusDiss2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities of Rusanov-type
!#
!# 19.) hydro_calcCharacteristics2d_sim
!#      -> Computes characteristic variables
!#
!# 20.) hydro_calcFluxFCTScDiss2d_sim
!#      -> Computes fluxes for FCT algorithm adopting scalar
!#         artificial viscosities
!#
!# 21.) hydro_calcFluxFCTRoeDiss2d_sim
!#      -> Computes fluxes for FCT algorithm adopting tensorial
!#         artificial viscosities of Roe-type
!#
!# 22.) hydro_calcFluxFCTRusDiss2d_sim
!#      -> Computes fluxes for FCT algorithm adopting scalar
!#         artificial viscosities of Rusanov-type
!#
!# 23.) hydro_trafoFluxDensity2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 24.) hydro_trafoDiffDensity2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 25. hydro_trafoNodalDensity2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density values
!#
!# 26.) hydro_trafoFluxEnergy2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 27.) hydro_trafoDiffEnergy2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 28. hydro_trafoNodalEnergy2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal energy values
!#
!# 29.) hydro_trafoFluxPressure2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 30.) hydro_trafoDiffPressure2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the pressure
!#
!# 31. hydro_trafoNodalPressure2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal pressure values
!#
!# 32.) hydro_trafoFluxVelocity2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 33.) hydro_trafoDiffVelocity2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 34.) hydro_trafoNodalVelocity2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal velocity values
!#
!# 35.) hydro_trafoFluxMomentum2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 36.) hydro_trafoDiffMomentum2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 37.) hydro_trafoNodalMomentum2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal momentum values
!#
!# 38.) hydro_trafoFluxDenEng2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 39.) hydro_trafoDiffDenEng2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 40.) hydro_trafoNodalDenEng2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density and energy values
!#
!# 41.) hydro_trafoFluxDenPre2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 42.) hydro_trafoDiffDenPre2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 43.) hydro_trafoNodalDenPre2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density and pressure values
!#
!# 44.) hydro_trafoFluxDenPreVel2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 45.) hydro_trafoDiffDenPreVel2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure 
!#         and the velocity
!#
!# 46.) hydro_trafoNodalDenPreVel2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density, pressure and velocity values
!#
!# 47.) hydro_calcBoundaryvalues2d
!#      -> Computes the boundary values for a given node
!#
!# 48.) hydro_hadaptCallbackScalar2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in interleave format
!#
!# 49.) hydro_hadaptCallbackBlock2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in block format
!#
!# 50.) hydro_coeffVectorBdr2d_sim
!#      -> Calculates the coefficients for the linear form in 2D
!#
!# </purpose>
!##############################################################################

module hydro_callback2d

#define HYDRO_NDIM 2
#include "hydro.h"

  use boundary
  use boundaryaux
  use boundarycondaux
  use collection
  use cubature
  use derivatives
  use domainintegration
  use element
  use feevaluation
  use flagship_callback
  use fparser
  use fsystem
  use genoutput
  use graph
  use groupfemsystem
  use hadaptaux
  use hydro_basic
  use linearsystemblock
  use linearsystemscalar
  use mprimitives
  use problem
  use scalarpde
  use solveraux
  use spatialdiscretisation
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
  public :: hydro_calcFluxFCTScDiss2d_sim
  public :: hydro_calcFluxFCTRoeDiss2d_sim
  public :: hydro_calcFluxFCTRusDiss2d_sim
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
  public :: hydro_trafoNodalDensity2d_sim
  public :: hydro_trafoNodalEnergy2d_sim
  public :: hydro_trafoNodalPressure2d_sim
  public :: hydro_trafoNodalVelocity2d_sim
  public :: hydro_trafoNodalMomentum2d_sim
  public :: hydro_trafoNodalDenEng2d_sim
  public :: hydro_trafoNodalDenPre2d_sim
  public :: hydro_trafoNodalDenPreVel2d_sim
  public :: hydro_calcBoundaryvalues2d
  public :: hydro_coeffVectorBdr2d_sim
  public :: hydro_hadaptCallbackScalar2d
  public :: hydro_hadaptCallbackBlock2d

  public :: hydro_calcDivVecScDiss2d_cuda

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGal2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
    real(DP), dimension(NVAR2D) :: Fxi,Fxj,Fyi,Fyj
#else
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
#endif
    real(DP) :: pi,pj,ui,uj,vi,vj
    integer :: idx

    
    do idx = 1, nedges

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute pressures
      pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      Fxi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)

      Fxj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)

      ! Compute fluxes for x-direction
      Fyi(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)

      Fyj(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)

      ! Assemble skew-symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
           IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fyj-&
           IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
           IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fyi )
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
      ! Compute flux difference for x-direction
      Fx_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
                        
      ! Compute flux difference for y-direction
      Fy_ij(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      
      ! Assemble fluxes
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
           IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
           IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij)
#endif
    end do

  end subroutine hydro_calcFluxGal2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGalNoBdr2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
    real(DP) :: pi,pj,ui,uj,vi,vj
    integer :: idx

    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute pressures
      pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute flux difference for x-direction
      Fx_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
                        
      ! Compute flux difference for y-direction
      Fy_ij(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      
      ! Assemble symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (0.5_DP*(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)-&
                   IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0))*Fx_ij+&
           0.5_DP*(IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)-&
                   IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0))*Fy_ij)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
    end do

  end subroutine hydro_calcFluxGalNoBdr2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxScDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
    real(DP), dimension(NVAR2D) :: Fxi,Fxj,Fyi,Fyj
#else
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj
    real(DP) :: H_ij,anorm,aux,c_ij,d_ij,q_ij,u_ij,v_ij,vel_ij
    integer :: idx

    
    do idx = 1, nedges

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute pressures
      pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)

      ! Compute fluxes for x-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                      INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                      INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                      INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                      INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                      INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                      INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                      INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                      INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient
      a = RCONST(0.5)*(IDX3(DmatrixCoeffsAtEdge,:,1,idx,0,0,0)-&
                       IDX3(DmatrixCoeffsAtEdge,:,2,idx,0,0,0))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      ! Compute densities
      ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute enthalpies
      hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
      hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(ri,rj)
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      H_ij = ROE_MEAN_VALUE(hi,hj,aux)
      
      ! Compute auxiliary variables
      vel_ij = u_ij*a(1) + v_ij*a(2)
      q_ij   = RCONST(0.5)*(u_ij*u_ij+v_ij*v_ij)

      ! Compute the speed of sound
      c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))
      
      ! Compute scalar dissipation
      d_ij = abs(vel_ij) + anorm*c_ij

      ! Multiply the solution difference by the scalar dissipation
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
                   IDX3(DdataAtEdge,:,1,idx,0,0,0))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
           IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fyj-&
           IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
           IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fyi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
           IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
           IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxScDiss2d_sim

  ! ***************************************************************************

!<subroutine>

  subroutine hydro_calcDivVecScDiss2d_cuda(rafcstab, rx, ry, dscale,&
      bclear, rcollection)

    use afcstabilisation
    use collection
    use fsystem
    use linearsystemblock

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 2D using scalar artificial viscosities proportional to the
    ! spectral radius (largest eigenvalue) of the Roe-matrix.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! divergence vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT
    
    ! local variables
    integer, dimension(:), pointer :: p_IverticesAtEdgeIdx
    integer :: IEDGEmax,IEDGEset,igroup
    integer(I64) :: p_DmatrixCoeffsAtEdge
    integer(I64) :: p_IverticesAtEdge
    integer(I64) :: p_Dx, p_Dy
    
    
    ! Check if edge structure is available on device and copy it otherwise
    if (storage_getMemoryAddress(rafcstab%h_IverticesAtEdge) .eq. 0_I64)&
        call afcstab_copyH2D_IverticesAtEdge(rafcstab)
    p_IverticesAtEdge = storage_getMemoryAddress(rafcstab%h_IverticesAtEdge)
    
    ! Check if matrix coefficients are available on device and copy it otherwise
    if (storage_getMemoryAddress(rafcstab%h_DmatrixCoeffsAtEdge) .eq. 0_I64)&
        call afcstab_copyH2D_DmatCoeffAtEdge(rafcstab)
    p_DmatrixCoeffsAtEdge = storage_getMemoryAddress(rafcstab%h_DmatrixCoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    p_Dx = storage_getMemoryAddress(rx%h_Ddata)
    if (p_Dx .eq. 0_I64) then
      call lsysbl_copyH2D_Vector(rx, .false.)
      p_Dx = storage_getMemoryAddress(rx%h_Ddata)
    end if
   
    ! Make sure that the destination vector ry exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Vector(ry, .true.)
    p_Dy = storage_getMemoryAddress(ry%h_Ddata)
   
    ! Set pointer
    call afcstab_getbase_IvertAtEdgeIdx(rafcstab, p_IverticesAtEdgeIdx)
    
    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronize memory access
    do igroup = 1, size(p_IverticesAtEdgeIdx)-1
      
      ! Do nothing for empty groups
      if (p_IverticesAtEdgeIdx(igroup+1)-p_IverticesAtEdgeIdx(igroup) .le. 0) cycle

      ! Get position of first edge in group
      IEDGEset = p_IverticesAtEdgeIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IverticesAtEdgeIdx(igroup+1)-1
      
      ! Use callback function to compute internodal fluxes
      call hydro_calcFluxScDiss2d_cuda(p_DmatrixCoeffsAtEdge, p_IverticesAtEdge,&
          p_Dx, p_Dy, dscale, rafcstab%NEQ, rafcstab%NVAR, rafcstab%NEDGE,&
          rafcstab%nmatCoeff, IEDGEmax-IEDGEset+1, IEDGEset)
    end do

    ! Transfer destination vector back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Vector(ry, bclear)

#else

    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivVecScDiss2d_cuda')
    call sys_halt()
    
#endif

  end subroutine hydro_calcDivVecScDiss2d_cuda

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxScDissDiSp2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
    real(DP), dimension(NVAR2D) :: Fxi,Fxj,Fyi,Fyj
#else
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj
    real(DP) :: H_ij,aux,c_ij,d_ij,q_ij,u_ij,v_ij
    integer :: idx
    
    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute pressures
      pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      Fxi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)

      Fxj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)

      ! Compute fluxes for x-direction
      Fyi(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)

      Fyj(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
#else
      ! Compute flux difference for x-direction
      Fx_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
                        
      ! Compute flux difference for y-direction
      Fy_ij(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
#endif
      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient
      a = 0.5_DP*(IDX3(DmatrixCoeffsAtEdge,:,1,idx,0,0,0)-&
                  IDX3(DmatrixCoeffsAtEdge,:,2,idx,0,0,0))

      ! Compute densities
      ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute enthalpies
      hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
      hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(ri,rj)
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      H_ij = ROE_MEAN_VALUE(hi,hj,aux)
      
      ! Compute auxiliary variable
      q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)

! Compute the speed of sound
      c_ij = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(H_ij-q_ij), SYS_EPSREAL_DP))
      
      ! Compute scalar dissipation with dimensional splitting
      d_ij = ( abs(a(1)*u_ij) + abs(a(1))*c_ij +&
               abs(a(2)*v_ij) + abs(a(2))*c_ij )

      ! Multiply the solution difference by the scalar dissipation
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
                   IDX3(DdataAtEdge,:,1,idx,0,0,0))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
           IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fyj-&
           IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
           IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fyi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
           IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
           IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxScDissDiSp2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRoeDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
    real(DP), dimension(NVAR2D) :: Fxi,Fxj,Fyi,Fyj
#else
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj
    real(DP) :: H_ij,c2_ij,c_ij,q_ij,u_ij,v_ij,vel_ij
    real(DP) :: anorm,aux,aux1,aux2
    real(DP) :: l1,l2,l3,l4,w1,w2,w3,w4
    integer :: idx


    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute pressures
      pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      Fxi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)

      Fxj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)

      ! Compute fluxes for x-direction
      Fyi(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)

      Fyj(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
#else
      ! Compute flux difference for x-direction
      Fx_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
                        
      ! Compute flux difference for y-direction
      Fy_ij(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = 0.5_DP*(IDX3(DmatrixCoeffsAtEdge,:,1,idx,0,0,0)-&
                  IDX3(DmatrixCoeffsAtEdge,:,2,idx,0,0,0))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Normalize the skew-symmetric coefficient
        a = a/anorm
        
        ! Compute densities
        ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        
        ! Compute enthalpies
        hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
        hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2)
        q_ij   = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)

        ! Compute the speed of sound
        c2_ij = max(((HYDRO_GAMMA)-1.0_DP)*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij  = sqrt(c2_ij)
        
        ! Compute eigenvalues
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)
        
        ! Compute solution difference U_j-U_i
        Diff = IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
               IDX3(DdataAtEdge,:,1,idx,0,0,0)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-1.0_DP)*(q_ij*Diff(1)&
                                      -u_ij*Diff(2)&
                                      -v_ij*Diff(3)&
                                           +Diff(4))/2.0_DP/c2_ij
        aux2 = (vel_ij*Diff(1)&
                 -a(1)*Diff(2)&
                 -a(2)*Diff(3))/2.0_DP/c_ij
        
        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((1.0_DP-((HYDRO_GAMMA)-1.0_DP)*q_ij/c2_ij)*Diff(1)&
                                +((HYDRO_GAMMA)-1.0_DP)*(u_ij*Diff(2)&
                                                        +v_ij*Diff(3)&
                                                             -Diff(4))/c2_ij)
        w3 = l3 * (aux1 - aux2)
        w4 = l4 * ((a(1)*v_ij-a(2)*u_ij)*Diff(1)&
                                   +a(2)*Diff(2)&
                                   -a(1)*Diff(3))
        
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        Diff(1) = anorm * ( w1 + w2 + w3 )
        Diff(2) = anorm * ( (u_ij-c_ij*a(1))*w1 + u_ij*w2 +&
                            (u_ij+c_ij*a(1))*w3 + a(2)*w4 )
        Diff(3) = anorm * ( (v_ij-c_ij*a(2))*w1 + v_ij*w2 +&
                            (v_ij+c_ij*a(2))*w3 - a(1)*w4 )
        Diff(4) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 +&
                            (H_ij+c_ij*vel_ij)*w3 + (u_ij*a(2)-v_ij*a(1))*w4 )

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
            (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
             IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fyj-&
             IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
             IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fyi + Diff)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
            (IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
             IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij + Diff)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
            (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
             IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij + Diff)
#endif
      else
        
#ifdef HYDRO_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
            (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
             IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fyj-&
             IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
             IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fyi)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
            (IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
             IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
            (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
             IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij)
#endif
      end if
    end do

  end subroutine hydro_calcFluxRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRoeDissDiSp2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
    real(DP), dimension(NVAR2D) :: Fxi,Fxj,Fyi,Fyj
#else
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
#endif
    real(DP), dimension(NVAR2D) :: DiffX,DiffY
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj
    real(DP) :: H_ij,c2_ij,c_ij,q_ij,u_ij,v_ij
    real(DP) :: anorm,aux,aux1,aux2
    real(DP) :: l1,l2,l3,l4,w1,w2,w3,w4
    integer :: idx


    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute pressures
      pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      Fxi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)

      Fxj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)

      ! Compute fluxes for x-direction
      Fyi(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)

      Fyj(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
#else
      ! Compute flux difference for x-direction
      Fx_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
                        
      ! Compute flux difference for y-direction
      Fy_ij(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
#endif
      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = 0.5_DP*(IDX3(DmatrixCoeffsAtEdge,:,1,idx,0,0,0)-&
                  IDX3(DmatrixCoeffsAtEdge,:,2,idx,0,0,0))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      if (anorm .gt. SYS_EPSREAL_DP) then

        ! Compute the absolute value
        a = abs(a)
        
        ! Compute densities
        ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        
        ! Compute enthalpies
        hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
        hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variable
        q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)

        ! Compute the speed of sound
        c2_ij = max(((HYDRO_GAMMA)-1.0_DP)*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij  = sqrt(c2_ij)

        !-----------------------------------------------------------------------
        ! Dimensional splitting: x-direction
        !-----------------------------------------------------------------------
        
        ! Compute eigenvalues
        l1 = abs(u_ij-c_ij)
        l2 = abs(u_ij)
        l3 = abs(u_ij+c_ij)
        l4 = abs(u_ij)
        
        ! Compute solution difference U_j-U_i
        DiffX = IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
                IDX3(DdataAtEdge,:,1,idx,0,0,0)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-1.0_DP)*(q_ij*DiffX(1)&
                                      -u_ij*DiffX(2)&
                                      -v_ij*DiffX(3)&
                                           +DiffX(4))/2.0_DP/c2_ij
        aux2 = (u_ij*DiffX(1)&
                    -DiffX(2))/2.0_DP/c_ij
        
        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((1.0_DP-((HYDRO_GAMMA)-1.0_DP)*q_ij/c2_ij)*DiffX(1)&
                                +((HYDRO_GAMMA)-1.0_DP)*(u_ij*DiffX(2)&
                                                        +v_ij*DiffX(3)&
                                                             -DiffX(4))/c2_ij)
        w3 = l3 * (aux1 - aux2)
        w4 = l4 * (v_ij*DiffX(1)&
                       -DiffX(3))
        
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        DiffX(1) = a(1) * ( w1 + w2 + w3 )
        DiffX(2) = a(1) * ( (u_ij-c_ij)*w1 + u_ij*w2 + (u_ij+c_ij)*w3 )
        DiffX(3) = a(1) * (        v_ij*w1 + v_ij*w2 +        v_ij*w3 - w4 )
        DiffX(4) = a(1) * ( (H_ij-c_ij*u_ij)*w1 + q_ij*w2 +&
                            (H_ij+c_ij*u_ij)*w3 - v_ij*w4 )
        
        !-----------------------------------------------------------------------
        ! Dimensional splitting: y-direction
        !-----------------------------------------------------------------------

        ! Compute eigenvalues
        l1 = abs(v_ij-c_ij)
        l2 = abs(v_ij)
        l3 = abs(v_ij+c_ij)
        l4 = abs(v_ij)
        
        ! Compute solution difference U_j-U_i
        DiffY = IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
                IDX3(DdataAtEdge,:,1,idx,0,0,0)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-1.0_DP)*(q_ij*DiffY(1)&
                                      -u_ij*DiffY(2)&
                                      -v_ij*DiffY(3)&
                                           +DiffY(4))/2.0_DP/c2_ij
        aux2 = (v_ij*DiffY(1)&
                    -DiffY(3))/2.0_DP/c_ij

        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((1.0_DP-((HYDRO_GAMMA)-1.0_DP)*q_ij/c2_ij)*DiffY(1)&
                                +((HYDRO_GAMMA)-1.0_DP)*(u_ij*DiffY(2)&
                                                        +v_ij*DiffY(3)&
                                                             -DiffY(4))/c2_ij)
        w3 = l3 * (aux1 - aux2)
        w4 = l4 * (-u_ij*DiffY(1)&
                        +DiffY(2))
        
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        DiffY(1) = a(2) * ( w1 + w2 + w3 )
        DiffY(2) = a(2) * (        u_ij*w1 + u_ij*w2 +        u_ij*w3 + w4 )
        DiffY(3) = a(2) * ( (v_ij-c_ij)*w1 + v_ij*w2 + (v_ij+c_ij)*w3 )
        DiffY(4) = a(2) * ( (H_ij-c_ij*v_ij)*w1 + q_ij*w2 +&
                            (H_ij+c_ij*v_ij)*w3 + u_ij*w4 )

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
            (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
             IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fyj-&
             IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
             IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fyi+&
             DiffX + DiffY)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
            (IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
             IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij+&
             DiffX + DiffY)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
            (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
             IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij+&
             DiffX + DiffY)
#endif
      else
        
#ifdef HYDRO_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
            (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
             IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fyj-&
             IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
             IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fyi)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
            (IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
             IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
            (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
             IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij)
#endif
      end if
    end do

  end subroutine hydro_calcFluxRoeDissDiSp2d_sim
 
  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRusDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
    real(DP), dimension(NVAR2D) :: Fxi,Fxj,Fyi,Fyj
#else
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: Ei,Ej,ci,cj,pi,pj,ui,uj,vi,vj
    real(DP) :: d_ij
    integer :: idx


    do idx = 1, nedges

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute pressures
      pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      Fxi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)

      Fxj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)

      ! Compute fluxes for x-direction
      Fyi(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)

      Fyj(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
#else
      ! Compute flux difference for x-direction
      Fx_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
                        
      ! Compute flux difference for y-direction
      Fy_ij(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !---------------------------------------------------------------------------
      
      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL_DP))
      
#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(0.5_DP*(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0))*uj+&
                      0.5_DP*(IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0))*vj)+&
                 0.5_DP*sqrt((IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0))**2+&
                             (IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0))**2)*cj,&
                  abs(0.5_DP*(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0))*ui+&
                      0.5_DP*(IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0))*vi)+&
                 0.5_DP*sqrt((IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0))**2+&
                             (IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0))**2)*ci )
#else
      ! Compute scalar dissipation
      d_ij = max( abs(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*uj+&
                      IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*vj)+&
                 sqrt(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)**2+&
                      IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)**2)*cj,&
                  abs(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*ui+&
                      IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*vi)+&
                 sqrt(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)**2+&
                      IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)**2)*ci )
#endif

      ! Multiply the solution difference by the scalar dissipation
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
                   IDX3(DdataAtEdge,:,1,idx,0,0,0))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------
      
#ifdef HYDRO_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
           IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fyj-&
           IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
           IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fyi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
           IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
           IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxRusDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRusDissDiSp2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
    real(DP), dimension(NVAR2D) :: Fxi,Fxj,Fyi,Fyj
#else
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: Ei,Ej,ci,cj,pi,pj,ui,uj,vi,vj
    real(DP) :: d_ij
    integer :: idx
    

    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute pressures
      pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      Fxi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fxi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)

      Fxj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fxj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)

      ! Compute fluxes for x-direction
      Fyi(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)
      Fyi(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)

      Fyj(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fyj(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
#else
      ! Compute flux difference for x-direction
      Fx_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fx_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                 INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
                        
      ! Compute flux difference for y-direction
      Fy_ij(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
      Fy_ij(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,vi,pi)-&
                 INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,vj,pj)
#endif
      
      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !-------------------------------------------------------------------------

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL_DP))

#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation with dimensional splitting based on
      ! the skew-symmetric part which does not include the symmetric
      ! boundary contribution
      d_ij = max( abs(0.5_DP*(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0))*uj)+&
                  abs(0.5_DP*(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)))*cj,&
                  abs(0.5_DP*(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0))*ui)+&
                  abs(0.5_DP*(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)))*ci )&
           + max( abs(0.5_DP*(IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0))*vj)+&
                  abs(0.5_DP*(IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)))*cj,&
                  abs(0.5_DP*(IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0))*vi)+&
                  abs(0.5_DP*(IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)))*ci )
#else      
      ! Compute scalar dissipation with dimensional splitting
      d_ij = max( abs(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*uj)+&
                  abs(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0))*cj,&
                  abs(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*ui)+&
                  abs(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0))*ci )&
           + max( abs(IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*vj)+&
                  abs(IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0))*cj,&
                  abs(IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*vi)+&
                  abs(IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0))*ci )
#endif

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
                   IDX3(DdataAtEdge,:,1,idx,0,0,0))
      
      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
           IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fyj-&
           IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
           IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fyi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
           IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
           IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxRusDissDiSp2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiagMatD2d_sim(DdataAtNode, DmatrixCoeffsAtNode,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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


    do inode = 1, nnodes
      
      ! Compute auxiliary variables
      ui = XVELOCITY2(DdataAtNode,IDX2,inode,0,0)
      vi = YVELOCITY2(DdataAtNode,IDX2,inode,0,0)
      
#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ii = diag(A_i)*C_{ii}$
      IDX3(DcoefficientsAtNode,1,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,_)
      IDX3(DcoefficientsAtNode,2,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,_)
      IDX3(DcoefficientsAtNode,3,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,_)
      IDX3(DcoefficientsAtNode,4,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,_)
#else
      ! Compute Galerkin coefficient $K_ii = -diag(A_i)*C_{ii}$
      IDX3(DcoefficientsAtNode,1,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,_)
      IDX3(DcoefficientsAtNode,2,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,_)
      IDX3(DcoefficientsAtNode,3,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,_)
      IDX3(DcoefficientsAtNode,4,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,_)
#endif
    end do

  end subroutine hydro_calcMatDiagMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiag2d_sim(DdataAtNode,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all nodes under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP) :: Ei,ui,vi
    integer :: inode


    do inode = 1, nnodes
      
      ! Compute auxiliary variables
      ui = XVELOCITY2(DdataAtNode,IDX2,inode,0,0)
      vi = YVELOCITY2(DdataAtNode,IDX2,inode,0,0)
      Ei = SPECIFICTOTALENERGY2(DdataAtNode,IDX2,inode,0,0)

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      IDX3(DcoefficientsAtNode,1,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,2,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,3,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,4,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX41(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,5,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,6,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,7,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,8,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX42(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,9,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,10,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,11,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,12,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX43(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,13,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX14(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,14,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX24(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,15,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX34(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,16,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                    IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
#else
      ! Compute Galerkin coefficient $K_ii = -A_i*C_{ii}$
      IDX3(DcoefficientsAtNode,1,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,2,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,3,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,4,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX41(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,5,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,6,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,7,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,8,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX42(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,9,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,10,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,11,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,12,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX43(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,13,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX14(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,14,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX24(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,15,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX34(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtNode,16,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX2(DmatrixCoeffsAtNode,1,inode,0,0),
                                     IDX2(DmatrixCoeffsAtNode,2,inode,0,0),ui,vi,Ei)
#endif
    end do

  end subroutine hydro_calcMatDiag2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGalMatD2d_sim(DdataAtEdge,&
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

  ! Number of edges,
  integer, intent(in) :: nedges
!</input>

!<inputoutput>
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Nullify dissipation tensor
      IDX3(DcoefficientsAtEdge,:,1,idx,0,0,0) = 0.0_DP

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,_)
#endif
    end do

  end subroutine hydro_calcMatGalMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGal2d_sim(DdataAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP) :: Ei,Ej,ui,uj,vi,vj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Nullify dissipation tensor
      IDX3(DcoefficientsAtEdge,:,1,idx,0,0,0) = 0.0_DP

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,10,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,11,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,12,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,13,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,14,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,15,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,16,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,5,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,6,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,7,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,8,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,9,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,10,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,11,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,12,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,13,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,14,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,15,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,16,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,10,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,11,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,12,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,13,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,14,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,15,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,16,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,5,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,6,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,7,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,8,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,9,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,10,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,11,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,12,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,13,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,14,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,15,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,16,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
#endif
    end do
      
  end subroutine hydro_calcMatGal2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDissMatD2d_sim(DdataAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj
    real(DP) :: H_ij,anorm,aux,c_ij,q_ij,u_ij,v_ij,vel_ij
    integer :: idx

    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vii,_)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,_)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5_DP*(IDX3(DmatrixCoeffsAtEdge,:,1,idx,0,0,0)-&
                  IDX3(DmatrixCoeffsAtEdge,:,2,idx,0,0,0))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Compute densities
        ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        
        ! Compute pressures
        pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
        pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

        ! Compute enthalpies
        hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
        hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
               
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2)
        q_ij   = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)

        ! Compute the speed of sound
        c_ij = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(H_ij-q_ij), SYS_EPSREAL_DP))

        ! Compute scalar dissipation
        IDX3(DcoefficientsAtEdge,:,1,idx,0,0,0) = dscale * (abs(vel_ij) + anorm*c_ij)

      else
        
        ! Nullify dissipation tensor
        IDX3(DcoefficientsAtEdge,:,1,idx,0,0,0) = 0.0_DP

      end if
    end do

  end subroutine hydro_calcMatScDissMatD2d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDiss2d_sim(DdataAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: Ei,Ej,hi,hj,ri,rj,pi,pj,ui,uj,vi,vj
    real(DP) :: H_ij,anorm,aux,c_ij,q_ij,u_ij,v_ij,vel_ij
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,10,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,11,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,12,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,13,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,14,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,15,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,16,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,5,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,6,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,7,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,8,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,9,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,10,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,11,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,12,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,13,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,14,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,15,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,16,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,10,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,11,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,12,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,13,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,14,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,15,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,16,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,5,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,6,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,7,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,8,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,9,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,10,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,11,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,12,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,13,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,14,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,15,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,16,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5_DP*(IDX3(DmatrixCoeffsAtEdge,:,1,idx,0,0,0)-&
                  IDX3(DmatrixCoeffsAtEdge,:,2,idx,0,0,0))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Compute densities
        ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        
        ! Compute pressures
        pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
        pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

        ! Compute enthalpies
        hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
        hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2)
        q_ij   = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)

        ! Compute the speed of sound
        c_ij = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(H_ij-q_ij), SYS_EPSREAL_DP))

        ! Compute scalar dissipation
        aux = dscale * (abs(vel_ij) + anorm*c_ij)
        
        IDX3(DcoefficientsAtEdge, :,1,idx,0,0,0) = 0.0_DP
        IDX3(DcoefficientsAtEdge, 1,1,idx,0,0,0) = aux
        IDX3(DcoefficientsAtEdge, 6,1,idx,0,0,0) = aux
        IDX3(DcoefficientsAtEdge,11,1,idx,0,0,0) = aux
        IDX3(DcoefficientsAtEdge,16,1,idx,0,0,0) = aux

      else

        ! Nullify dissipation tensor
        IDX3(DcoefficientsAtEdge,:,1,idx,0,0,0) = 0.0_DP

      end if
    end do

  end subroutine hydro_calcMatScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDissMatD2d_sim(DdataAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj
    real(DP) :: H_ij,anorm,aux,cPow2_ij,c_ij,q_ij,u_ij,v_ij,vel_ij
    real(DP) :: l1,l2,l3,l4
    integer :: idx


    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vii,_)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,_)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !---------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5_DP*(IDX3(DmatrixCoeffsAtEdge,:,1,idx,0,0,0)-&
                  IDX3(DmatrixCoeffsAtEdge,:,2,idx,0,0,0))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL_DP) then

        ! Normalize the skew-symmetric coefficient
        a = a/anorm

        ! Compute densities
        ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        
        ! Compute pressures
        pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
        pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

        ! Compute enthalpies
        hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
        hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1)+v_ij*a(2)
        q_ij   = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
        
        ! Compute speed of sound
        cPow2_ij = max(((HYDRO_GAMMA)-1.0_DP)*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij     = sqrt(cPow2_ij)
        
        ! Diagonal matrix of eigenvalues
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)

        ! Matrix of right eigenvectors
        R_ij(1,1) =  l1
        R_ij(2,1) =  l1*(u_ij-c_ij*a(1))
        R_ij(3,1) =  l1*(v_ij-c_ij*a(2))
        R_ij(4,1) =  l1*(H_ij-c_ij*vel_ij)
        
        R_ij(1,2) =  l2
        R_ij(2,2) =  l2*u_ij
        R_ij(3,2) =  l2*v_ij
        R_ij(4,2) =  l2*q_ij
        
        R_ij(1,3) =  l3
        R_ij(2,3) =  l3*(u_ij+c_ij*a(1))
        R_ij(3,3) =  l3*(v_ij+c_ij*a(2))
        R_ij(4,3) =  l3*(H_ij+c_ij*vel_ij)
        
        R_ij(1,4) =  0.0_DP
        R_ij(2,4) =  l4*a(2)
        R_ij(3,4) = -l4*a(1)
        R_ij(4,4) =  l4*(u_ij*a(2)-v_ij*a(1))
        
        ! Matrix of left eigenvectors
        L_ij(1,1) =  0.5_DP*(((HYDRO_GAMMA)-1.0_DP)*q_ij+c_ij*vel_ij)/cPow2_ij
        L_ij(2,1) =  (cPow2_ij-((HYDRO_GAMMA)-1.0_DP)*q_ij)/cPow2_ij
        L_ij(3,1) =  0.5_DP*(((HYDRO_GAMMA)-1.0_DP)*q_ij-c_ij*vel_ij)/cPow2_ij
        L_ij(4,1) =  v_ij*a(1)-u_ij*a(2)
        
        L_ij(1,2) =  0.5_DP*(-((HYDRO_GAMMA)-1.0_DP)*u_ij-c_ij*a(1))/cPow2_ij
        L_ij(2,2) =  ((HYDRO_GAMMA)-1.0_DP)*u_ij/cPow2_ij
        L_ij(3,2) =  0.5_DP*(-((HYDRO_GAMMA)-1.0_DP)*u_ij+c_ij*a(1))/cPow2_ij
        L_ij(4,2) =  a(2)

        L_ij(1,3) =  0.5_DP*(-((HYDRO_GAMMA)-1.0_DP)*v_ij-c_ij*a(2))/cPow2_ij
        L_ij(2,3) =  ((HYDRO_GAMMA)-1.0_DP)*v_ij/cPow2_ij
        L_ij(3,3) =  0.5_DP*(-((HYDRO_GAMMA)-1.0_DP)*v_ij+c_ij*a(2))/cPow2_ij
        L_ij(4,3) = -a(1)
        
        L_ij(1,4) =  0.5_DP*((HYDRO_GAMMA)-1.0_DP)/cPow2_ij
        L_ij(2,4) = -((HYDRO_GAMMA)-1.0_DP)/cPow2_ij
        L_ij(3,4) =  0.5_DP*((HYDRO_GAMMA)-1.0_DP)/cPow2_ij
        L_ij(4,4) =  0.0_DP
        
        ! Include scaling parameter
        anorm = dscale*anorm
        
        ! Compute tensorial dissipation D_ij = diag(R_ij*|Lbd_ij|*L_ij)*I
        IDX3(DcoefficientsAtEdge,:,1,idx,0,0,0) = 0.0_DP
        IDX3(DcoefficientsAtEdge,1,1,idx,0,0,0) = anorm*( R_ij(1,1)*L_ij(1,1)+&
                                                           R_ij(1,2)*L_ij(2,1)+&
                                                           R_ij(1,3)*L_ij(3,1)+&
                                                           R_ij(1,4)*L_ij(4,1)  )
        IDX3(DcoefficientsAtEdge,2,1,idx,0,0,0) = anorm*( R_ij(2,1)*L_ij(1,2)+&
                                                           R_ij(2,2)*L_ij(2,2)+&
                                                           R_ij(2,3)*L_ij(3,2)+&
                                                           R_ij(2,4)*L_ij(4,2)  )
        IDX3(DcoefficientsAtEdge,3,1,idx,0,0,0) = anorm*( R_ij(3,1)*L_ij(1,3)+&
                                                           R_ij(3,2)*L_ij(2,3)+&
                                                           R_ij(3,3)*L_ij(3,3)+&
                                                           R_ij(3,4)*L_ij(4,3)  )
        IDX3(DcoefficientsAtEdge,4,1,idx,0,0,0) = anorm*( R_ij(4,1)*L_ij(1,4)+&
                                                           R_ij(4,2)*L_ij(2,4)+&
                                                           R_ij(4,3)*L_ij(3,4)+&
                                                           R_ij(4,4)*L_ij(4,4)  )
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0_DP
        
      end if
    end do

  end subroutine hydro_calcMatRoeDissMatD2d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDiss2d_sim(DdataAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
    real(DP) :: Ei,Ej,hi,hj,pi,pj,ri,rj,ui,uj,vi,vj
    real(DP) :: H_ij,anorm,aux,cPow2_ij,c_ij,q_ij,u_ij,v_ij,vel_ij
    real(DP) :: l1,l2,l3,l4
    integer :: idx,i,j,k


    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,10,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,11,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,12,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,13,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,14,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,15,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,16,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,5,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,6,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,7,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,8,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,9,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,10,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,11,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,12,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,13,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,14,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,15,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,16,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,10,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,11,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,12,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,13,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,14,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,15,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,16,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,5,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,6,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,7,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,8,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,9,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,10,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,11,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,12,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,13,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,14,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,15,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,16,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5_DP*(IDX3(DmatrixCoeffsAtEdge,:,1,idx,0,0,0)-&
                  IDX3(DmatrixCoeffsAtEdge,:,2,idx,0,0,0))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Normalize the skew-symmetric coefficient
        a = a/anorm

        ! Compute densities
        ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        
        ! Compute pressures
        pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
        pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

        ! Compute enthalpies
        hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
        hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1)+v_ij*a(2)
        q_ij   = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)

        ! Compute speed of sound
        cPow2_ij = max(((HYDRO_GAMMA)-1.0_DP)*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij     = sqrt(cPow2_ij)
        
        ! Diagonal matrix of eigenvalues
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)
        
        ! Matrix of right eigenvectors
        R_ij(1,1) =  l1
        R_ij(2,1) =  l1*(u_ij-c_ij*a(1))
        R_ij(3,1) =  l1*(v_ij-c_ij*a(2))
        R_ij(4,1) =  l1*(H_ij-c_ij*vel_ij)
        
        R_ij(1,2) =  l2
        R_ij(2,2) =  l2*u_ij
        R_ij(3,2) =  l2*v_ij
        R_ij(4,2) =  l2*q_ij
        
        R_ij(1,3) =  l3
        R_ij(2,3) =  l3*(u_ij+c_ij*a(1))
        R_ij(3,3) =  l3*(v_ij+c_ij*a(2))
        R_ij(4,3) =  l3*(H_ij+c_ij*vel_ij)
        
        R_ij(1,4) =  0.0_DP
        R_ij(2,4) =  l4*a(2)
        R_ij(3,4) = -l4*a(1)
        R_ij(4,4) =  l4*(u_ij*a(2)-v_ij*a(1))
        
        ! Matrix of left eigenvectors
        L_ij(1,1) =  0.5_DP*(((HYDRO_GAMMA)-1.0_DP)*q_ij+c_ij*vel_ij)/cPow2_ij
        L_ij(2,1) =  (cPow2_ij-((HYDRO_GAMMA)-1.0_DP)*q_ij)/cPow2_ij
        L_ij(3,1) =  0.5_DP*(((HYDRO_GAMMA)-1.0_DP)*q_ij-c_ij*vel_ij)/cPow2_ij
        L_ij(4,1) =  v_ij*a(1)-u_ij*a(2)
        
        L_ij(1,2) =  0.5_DP*(-((HYDRO_GAMMA)-1.0_DP)*u_ij-c_ij*a(1))/cPow2_ij
        L_ij(2,2) =  ((HYDRO_GAMMA)-1.0_DP)*u_ij/cPow2_ij
        L_ij(3,2) =  0.5_DP*(-((HYDRO_GAMMA)-1.0_DP)*u_ij+c_ij*a(1))/cPow2_ij
        L_ij(4,2) =  a(2)
        
        L_ij(1,3) =  0.5_DP*(-((HYDRO_GAMMA)-1.0_DP)*v_ij-c_ij*a(2))/cPow2_ij
        L_ij(2,3) =  ((HYDRO_GAMMA)-1.0_DP)*v_ij/cPow2_ij
        L_ij(3,3) =  0.5_DP*(-((HYDRO_GAMMA)-1.0_DP)*v_ij+c_ij*a(2))/cPow2_ij
        L_ij(4,3) = -a(1)
        
        L_ij(1,4) =  ((HYDRO_GAMMA)-1.0_DP)/2.0_DP/cPow2_ij
        L_ij(2,4) = -((HYDRO_GAMMA)-1.0_DP)/cPow2_ij
        L_ij(3,4) =  ((HYDRO_GAMMA)-1.0_DP)/2.0_DP/cPow2_ij
        L_ij(4,4) =  0.0_DP
        
        ! Include scaling parameter
        anorm = dscale*anorm

        ! Compute tensorial dissipation D_ij = R_ij*|Lbd_ij|*L_ij
        do i = 1, NVAR2D
          do j = 1, NVAR2D
            aux = 0.0_DP
            do k = 1, NVAR2D
              aux = aux + R_ij(i,k)*L_ij(k,j)
            end do
            IDX3(DcoefficientsAtEdge,NVAR2D*(j-1)+i,1,idx,0,0,0) = anorm*aux
          end do
        end do
        
      else
        
        ! Nullify dissipation tensor
        IDX3(DcoefficientsAtEdge,:,1,idx,0,0,0) = 0.0_DP
        
      end if
    end do

  end subroutine hydro_calcMatRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDissMatD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 2D and applies the scalar artificial viscosities of
    ! Rusanov-type.
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP) :: Ei,Ej,ci,cj,ui,uj,vi,vj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vii,_)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,_)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,_)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate scalar artificial dissipation of Rusanov-type
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL_DP))
      
      ! Compute dissipation tensor
      IDX3(DcoefficientsAtEdge,:,1,idx,0,0,0) = dscale *&
          max( abs(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*uj+&
                   IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*vj) +&
                   sqrt(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)**2+&
                        IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)**2)*cj,&
               abs(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*ui+&
                   IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*vi) +&
                   sqrt(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)**2+&
                        IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)**2)*ci )
    end do

  end subroutine hydro_calcMatRusDissMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDiss2d_sim(DdataAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP) :: Ei,Ej,aux,ci,cj,ui,uj,vi,vj
    integer :: idx


    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,10,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,11,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,12,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,13,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,14,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,15,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,16,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,5,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,6,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,7,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,8,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,9,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,10,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,11,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,12,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,13,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,14,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,15,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,16,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                    IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),ui,vi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DcoefficientsAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,10,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,11,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,12,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)

      IDX3(DcoefficientsAtEdge,13,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,14,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,15,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      IDX3(DcoefficientsAtEdge,16,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DcoefficientsAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX41(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,5,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,6,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,7,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,8,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX42(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,9,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,10,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,11,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,12,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX43(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)

      IDX3(DcoefficientsAtEdge,13,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX14(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,14,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX24(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,15,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX34(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
      IDX3(DcoefficientsAtEdge,16,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX44(-dscale,IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),
                                     IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),ui,vi,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate scalar artificial dissipation of Rusanov-type
      !---------------------------------------------------------------------------

      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL_DP))
      
      ! Compute dissipation tensor
      aux = dscale *&
          max( abs(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*uj+&
                   IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*vj) +&
                   sqrt(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)**2+&
                        IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)**2)*cj,&
               abs(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*ui+&
                   IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*vi) +&
                   sqrt(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)**2+&
                        IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)**2)*ci )

      IDX3(DcoefficientsAtEdge, :,1,idx,0,0,0) = 0.0_DP
      IDX3(DcoefficientsAtEdge, 1,1,idx,0,0,0) = aux
      IDX3(DcoefficientsAtEdge, 6,1,idx,0,0,0) = aux
      IDX3(DcoefficientsAtEdge,11,1,idx,0,0,0) = aux
      IDX3(DcoefficientsAtEdge,16,1,idx,0,0,0) = aux
    end do

  end subroutine hydro_calcMatRusDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcCharacteristics2d_sim(Dweight, DdataAtEdge,&
      nedges, DcharVariablesAtEdge, DeigenvaluesAtEdge,&
      DrightEigenvectorsAtEdge, DleftEigenvectorsAtEdge, rcollection)

!<description>
    ! This subroutine computes the characteristic variables in 2D.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! OPTIONAL: Characteristic variables for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each edge
    real(DP), dimension(:,:), intent(out), optional :: DcharVariablesAtEdge
    
    ! OPTIONAL: Eigenvalues for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each edge
    real(DP), dimension(:,:), intent(out), optional :: DeigenvaluesAtEdge
    
    ! OPTIONAL: Matrices of left eigenvectors for all edges under consideration
    !   DIMENSION(nvar*nvar,nedges)
    ! with nvar the number of variables at each edge
    real(DP), dimension(:,:), intent(out), optional :: DleftEigenvectorsAtEdge
    
    ! OPTIONAL: Matrices of right eigenvectors for all edges under consideration
    !   DIMENSION(nvar*nvar,nedges)
    ! with nvar the number of variables at each edge
    real(DP), dimension(:,:), intent(out), optional :: DrightEigenvectorsAtEdge
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(NVAR2D) :: a,Diff
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj
    real(DP) :: H_ij,anorm,aux,aux1,aux2,cPow2_ij,c_ij,q_ij,u_ij,v_ij
    integer :: idx


    ! Compute norm of weighting coefficient
    anorm = sqrt(Dweight(1)*Dweight(1)+Dweight(2)*Dweight(2))

    ! Check if weighting coefficient is zero
    if (anorm .le. SYS_EPSREAL_DP) then
      if (present(DcharVariablesAtEdge))     DcharVariablesAtEdge     = 0.0_DP
      if (present(DeigenvaluesAtEdge))       DeigenvaluesAtEdge       = 0.0_DP
      if (present(DrightEigenvectorsAtEdge)) DrightEigenvectorsAtEdge = 0.0_DP
      if (present(DleftEigenvectorsAtEdge))  DleftEigenvectorsAtEdge  = 0.0_DP

      ! That`s it
      return
    end if

    ! Compute normalised weighting coefficient
    a = Dweight/anorm

    ! Do we have to compute characteristic variables
    if (present(DcharVariablesAtEdge)) then
      do idx = 1, nedges
        
        ! Compute velocities
        ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        
        uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

        ! Compute densities
        ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        
        ! Compute pressures
        pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
        pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

        ! Compute enthalpies
        hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
        hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        q_ij  = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)

        ! Compute speed of sound
        cPow2_ij = max(((HYDRO_GAMMA)-1.0_DP)*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij = sqrt(cPow2_ij)
        aux  = u_ij*a(1)+v_ij*a(2)

        ! Compute solution difference U_j-U_i
        Diff = IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
               IDX3(DdataAtEdge,:,1,idx,0,0,0)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-1.0_DP)*(q_ij*Diff(1)&
                                      -u_ij*Diff(2)&
                                      -v_ij*Diff(3)&
                                           +Diff(4))/2.0_DP/cPow2_ij
        aux2 = (aux*Diff(1)&
              -a(1)*Diff(2)&
              -a(2)*Diff(3))/2.0_DP/c_ij

        ! Compute characteristic variables
        IDX2(DcharVariablesAtEdge,1,idx,0,0) = anorm * (aux1 + aux2)
        IDX2(DcharVariablesAtEdge,2,idx,0,0) = anorm *&
            ((1.0_DP-((HYDRO_GAMMA)-1.0_DP)*q_ij/cPow2_ij)*Diff(1)+&
                              ((HYDRO_GAMMA)-1.0_DP)*(u_ij*Diff(2)+&
                                                      v_ij*Diff(3)-&
                                                          Diff(4))/cPow2_ij)
        IDX2(DcharVariablesAtEdge,3,idx,0,0) = anorm * (aux1 - aux2)
        IDX2(DcharVariablesAtEdge,4,idx,0,0) = anorm * ((a(1)*v_ij-a(2)*u_ij)*Diff(1)+&
                                                                          a(2)*Diff(2)-&
                                                                          a(1)*Diff(3))
      end do
    end if

    ! Do we have to compute eigenvalues
    if (present(DeigenvaluesAtEdge)) then
      do idx = 1, nedges

        ! Compute velocities
        ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        
        uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

        ! Compute densities
        ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        
        ! Compute pressures
        pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
        pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

        ! Compute enthalpies
        hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
        hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
        c_ij = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(H_ij-q_ij), SYS_EPSREAL_DP))
        aux  = a(1)*u_ij+a(2)*v_ij

        ! Compute eigenvalues
        IDX2(DeigenvaluesAtEdge,1,idx,0,0) = aux-c_ij
        IDX2(DeigenvaluesAtEdge,2,idx,0,0) = aux
        IDX2(DeigenvaluesAtEdge,3,idx,0,0) = aux+c_ij
        IDX2(DeigenvaluesAtEdge,4,idx,0,0) = aux
      end do
    end if

    ! Do we have to compute right eigenvectors
    if (present(DrightEigenvectorsAtEdge)) then
      do idx = 1, nedges

        ! Compute velocities
        ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        
        uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

        ! Compute densities
        ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        
        ! Compute pressures
        pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
        pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

        ! Compute enthalpies
        hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
        hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
        c_ij = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(H_ij-q_ij), SYS_EPSREAL_DP))
        aux  = a(1)*u_ij+a(2)*v_ij

        ! Compute right eigenvectors
        IDX2(DrightEigenvectorsAtEdge, 1,idx,0,0) =  1.0_DP
        IDX2(DrightEigenvectorsAtEdge, 2,idx,0,0) =  u_ij-c_ij*a(1)
        IDX2(DrightEigenvectorsAtEdge, 3,idx,0,0) =  v_ij-c_ij*a(2)
        IDX2(DrightEigenvectorsAtEdge, 4,idx,0,0) =  H_ij-c_ij*aux

        IDX2(DrightEigenvectorsAtEdge, 5,idx,0,0) =  1.0_DP
        IDX2(DrightEigenvectorsAtEdge, 6,idx,0,0) =  u_ij
        IDX2(DrightEigenvectorsAtEdge, 7,idx,0,0) =  v_ij
        IDX2(DrightEigenvectorsAtEdge, 8,idx,0,0) =  q_ij

        IDX2(DrightEigenvectorsAtEdge, 9,idx,0,0) =  1.0_DP
        IDX2(DrightEigenvectorsAtEdge,10,idx,0,0) =  u_ij+c_ij*a(1)
        IDX2(DrightEigenvectorsAtEdge,11,idx,0,0) =  v_ij+c_ij*a(2)
        IDX2(DrightEigenvectorsAtEdge,12,idx,0,0) =  H_ij+c_ij*aux

        IDX2(DrightEigenvectorsAtEdge,13,idx,0,0) =  0.0_DP
        IDX2(DrightEigenvectorsAtEdge,14,idx,0,0) =  a(2)
        IDX2(DrightEigenvectorsAtEdge,15,idx,0,0) = -a(1)
        IDX2(DrightEigenvectorsAtEdge,16,idx,0,0) =  u_ij*a(2)-v_ij*a(1)
      end do
    end if

    ! Do we have to compute left eigenvectors
    if (present(DleftEigenvectorsAtEdge)) then
      do idx = 1, nedges

        ! Compute velocities
        ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        
        uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

        ! Compute densities
        ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        
        ! Compute pressures
        pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
        pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

        ! Compute enthalpies
        hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
        hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        q_ij     = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
        cPow2_ij = max(((HYDRO_GAMMA)-1.0_DP)*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij     = sqrt(cPow2_ij)
        aux      = a(1)*u_ij+a(2)*v_ij

        ! Compute left eigenvectors
        IDX2(DleftEigenvectorsAtEdge, 1,idx,0,0) =  0.5_DP*(((HYDRO_GAMMA)-1.0_DP)*q_ij+c_ij*aux)/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge, 2,idx,0,0) = (cPow2_ij-((HYDRO_GAMMA)-1.0_DP)*q_ij)/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge, 3,idx,0,0) =  0.5_DP*(((HYDRO_GAMMA)-1.0_DP)*q_ij-c_ij*aux)/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge, 4,idx,0,0) =  v_ij*a(1)-u_ij*a(2)

        IDX2(DleftEigenvectorsAtEdge, 5,idx,0,0) =  0.5_DP*(-((HYDRO_GAMMA)-1.0_DP)*u_ij-c_ij*a(1))/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge, 6,idx,0,0) =  ((HYDRO_GAMMA)-1.0_DP)*u_ij/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge, 7,idx,0,0) =  0.5_DP*(-((HYDRO_GAMMA)-1.0_DP)*u_ij+c_ij*a(1))/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge, 8,idx,0,0) =  a(2)

        IDX2(DleftEigenvectorsAtEdge, 9,idx,0,0) =  0.5_DP*(-((HYDRO_GAMMA)-1.0_DP)*v_ij-c_ij*a(2))/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge,10,idx,0,0) =  ((HYDRO_GAMMA)-1.0_DP)*v_ij/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge,11,idx,0,0) =  0.5_DP*(-((HYDRO_GAMMA)-1.0_DP)*v_ij+c_ij*a(2))/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge,12,idx,0,0) = -a(1)

        IDX2(DleftEigenvectorsAtEdge,13,idx,0,0) =  ((HYDRO_GAMMA)-1.0_DP)/2.0_DP/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge,14,idx,0,0) = -((HYDRO_GAMMA)-1.0_DP)/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge,15,idx,0,0) =  ((HYDRO_GAMMA)-1.0_DP)/2.0_DP/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge,16,idx,0,0) =  0.0_DP
      end do
    end if

  end subroutine hydro_calcCharacteristics2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTScDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj
    real(DP) :: H_ij,anorm,aux,c_ij,d_ij,q_ij,u_ij,v_ij,vel_ij
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute pressures
      pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute skew-symmetric coefficient
      a = 0.5_DP*(IDX3(DmatrixCoeffsAtEdge,:,1,idx,0,0,0)-&
                  IDX3(DmatrixCoeffsAtEdge,:,2,idx,0,0,0))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      ! Compute densities
      ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute enthalpies
      hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
      hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(ri,rj)
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      H_ij = ROE_MEAN_VALUE(hi,hj,aux)
      
      ! Compute auxiliary variables
      vel_ij = u_ij*a(1) + v_ij*a(2)
      q_ij   = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)

      ! Compute the speed of sound
      c_ij = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(H_ij-q_ij), SYS_EPSREAL_DP))
      
      ! Compute scalar dissipation
      d_ij = abs(vel_ij) + anorm*c_ij

      ! Compute conservative fluxes
      IDX2(DfluxesAtEdge,:,idx,0,0) = dscale*d_ij*&
          (IDX3(DdataAtEdge,:,1,idx,0,0,0)-IDX3(DdataAtEdge,:,2,idx,0,0,0))

    end do

  end subroutine hydro_calcFluxFCTScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRoeDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)


!<description>
    ! This subroutine computes the raw antidiffusive fluxes for FCT
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj
    real(DP) :: H_ij,c2_ij,c_ij,q_ij,u_ij,v_ij,vel_ij
    real(DP) :: anorm,aux,aux1,aux2
    real(DP) :: l1,l2,l3,l4,w1,w2,w3,w4
    integer :: idx


    do idx = 1, nedges

      ! Compute the skew-symmetric coefficient and its norm
      a = 0.5_DP*(IDX3(DmatrixCoeffsAtEdge,:,1,idx,0,0,0)-&
                  IDX3(DmatrixCoeffsAtEdge,:,2,idx,0,0,0))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      if (anorm .gt. SYS_EPSREAL_DP) then

        ! Normalize the skew-symmetric coefficient
        a = a/anorm

        ! Compute velocities
        ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        
        uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

        ! Compute pressures
        pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
        pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
        ! Compute densities
        ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        
        ! Compute enthalpies
        hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
        hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2)
        q_ij   = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)

        ! Compute the speed of sound
        c2_ij = max(((HYDRO_GAMMA)-1.0_DP)*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij  = sqrt(c2_ij)
        
        ! Compute eigenvalues
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)
        
        ! Compute solution difference U_i-U_j
        Diff = IDX3(DdataAtEdge,:,1,idx,0,0,0)-&
               IDX3(DdataAtEdge,:,2,idx,0,0,0)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-1.0_DP)*(q_ij*Diff(1)&
                                     -u_ij*Diff(2)&
                                     -v_ij*Diff(3)&
                                          +Diff(4))/2.0_DP/c2_ij

        aux2 = (vel_ij*Diff(1)&
                 -a(1)*Diff(2)&
                 -a(2)*Diff(3))/2.0_DP/c_ij
        
        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((1.0_DP-((HYDRO_GAMMA)-1.0_DP)*q_ij/c2_ij)*Diff(1)&
                                +((HYDRO_GAMMA)-1.0_DP)*(u_ij*Diff(2)&
                                                        +v_ij*Diff(3)&
                                                             -Diff(4))/c2_ij)
        w3 = l3 * (aux1 - aux2)
        w4 = l4 * ((a(1)*v_ij-a(2)*u_ij)*Diff(1)&
                                   +a(2)*Diff(2)&
                                   -a(1)*Diff(3))
        
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        Diff(1) = anorm * ( w1 + w2 + w3 )
        Diff(2) = anorm * ( (u_ij-c_ij*a(1))*w1 + u_ij*w2 +&
                            (u_ij+c_ij*a(1))*w3 + a(2)*w4 )
        Diff(3) = anorm * ( (v_ij-c_ij*a(2))*w1 + v_ij*w2 +&
                            (v_ij+c_ij*a(2))*w3 - a(1)*w4 )
        Diff(4) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 +&
                            (H_ij+c_ij*vel_ij)*w3 + (u_ij*a(2)-v_ij*a(1))*w4 )

        ! Compute conservative flux
        IDX2(DfluxesAtEdge,:,idx,0,0) = dscale*Diff
      else
        ! Cancel conservative flux
        IDX2(DfluxesAtEdge,:,idx,0,0) = 0
      end if
    end do

  end subroutine hydro_calcFluxFCTRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRusDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
    real(DP) :: Ei,Ej,ci,cj,ui,uj,vi,vj
    real(DP) :: d_ij
    integer :: idx
    
    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-1.0_DP)*(HYDRO_GAMMA)*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL_DP))
      
#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(0.5_DP*(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0))*uj+&
                      0.5_DP*(IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0))*vj)+&
                 0.5_DP*sqrt((IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0))**2+&
                             (IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0))**2)*cj,&
                  abs(0.5_DP*(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0))*ui+&
                      0.5_DP*(IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0))*vi)+&
                 0.5_DP*sqrt((IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0))**2+&
                             (IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)-&
                              IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0))**2)*ci )
#else
      ! Compute scalar dissipation
      d_ij = max( abs(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*uj+&
                      IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*vj)+&
                 sqrt(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)**2+&
                      IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)**2)*cj,&
                  abs(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*ui+&
                      IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*vi)+&
                 sqrt(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)**2+&
                      IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)**2)*ci )
#endif
      
      ! Compute conservative flux
      IDX2(DfluxesAtEdge,:,idx,0,0) = dscale*d_ij*&
          (IDX3(DdataAtEdge,:,1,idx,0,0,0)-IDX3(DdataAtEdge,:,2,idx,0,0,0))
    end do

  end subroutine hydro_calcFluxFCTRusDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDensity2d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density in 2D.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)
    end do

  end subroutine hydro_trafoFluxDensity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDensity2d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density in 2D.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
          DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
    end do

  end subroutine hydro_trafoDiffDensity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalDensity2d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the density in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DdataAtNode

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtNode
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, nnodes
      
      ! Transformed density values
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          DENSITY2(DdataAtNode,IDX2,idx,0,0)
    end do

  end subroutine hydro_trafoNodalDensity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxEnergy2d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the energy in 2D.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0)
    end do

  end subroutine hydro_trafoFluxEnergy2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffEnergy2d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the energy in 2D.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
          TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
    end do

  end subroutine hydro_trafoDiffEnergy2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalEnergy2d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the energy in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DdataAtNode

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtNode
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, nnodes
      
      ! Transformed energy values
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          TOTALENERGY2(DdataAtNode,IDX2,idx,0,0)
    end do

  end subroutine hydro_trafoNodalEnergy2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxPressure2d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the pressure in 2D.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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

    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          ((HYDRO_GAMMA)-1.0_DP)*(0.5_DP*(ui*ui+vi*vi)*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
                                  ui*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                                  vi*YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)+&
                                   TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -((HYDRO_GAMMA)-1.0_DP)*(0.5_DP*(uj*uj+vj*vj)*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
                                  uj*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                                  vj*YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)+&
                                   TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
    end do

  end subroutine hydro_trafoFluxPressure2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffPressure2d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the pressure in 2D.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
          PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
    end do

  end subroutine hydro_trafoDiffPressure2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalPressure2d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the pressure in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DdataAtNode

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtNode
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, nnodes
      
      ! Transformed pressure values
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          PRESSURE2(DdataAtNode,IDX2,idx,0,0)
    end do

  end subroutine hydro_trafoNodalPressure2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxVelocity2d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the x- and y-velocity.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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

    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Transformed velocity fluxes in x-direction
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          (XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
           ui*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0))/&
          DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -(XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
           uj*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0))/&
          DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Transformed velocity fluxes in y-direction
      IDX3(DtransformedFluxesAtEdge,2,1,idx,0,0,0) =&
          (YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
           vi*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0))/&
          DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,0,0,0) =&
         -(YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
           vj*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0))/&
          DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
    end do
    
  end subroutine hydro_trafoFluxVelocity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffVelocity2d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the x- and y-velocity.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
          XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      ! Transformed velocity difference in y-direction
      IDX2(DtransformedDataAtEdge,2,idx,0,0) =&
          YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
    end do
    
  end subroutine hydro_trafoDiffVelocity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalVelocity2d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the velocity in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DdataAtNode

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtNode
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, nnodes
      
      ! Transformed x-velocity values
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          XVELOCITY2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed y-velocity values
      IDX2(DtransformedDataAtNode,2,idx,0,0) =&
          YVELOCITY2(DdataAtNode,IDX2,idx,0,0)
    end do

  end subroutine hydro_trafoNodalVelocity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxMomentum2d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the x- and y-momentum.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)

      ! Transformed momentum fluxes in y-direction
      IDX3(DtransformedFluxesAtEdge,2,1,idx,0,0,0) =&
          YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,0,0,0) =&
         -YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)
    end do
    
  end subroutine hydro_trafoFluxMomentum2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffMomentum2d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the x- and y-momentum.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
          XMOMENTUM3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          XMOMENTUM3(DdataAtEdge,IDX3,1,idx,0,0,0)

      ! Transformed momentum difference in y-direction
      IDX2(DtransformedDataAtEdge,2,idx,0,0) =&
          YMOMENTUM3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          YMOMENTUM3(DdataAtEdge,IDX3,1,idx,0,0,0)
    end do
    
  end subroutine hydro_trafoDiffMomentum2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalMomentum2d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the momentum in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DdataAtNode

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtNode
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, nnodes
      
      ! Transformed x-momentum values
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          XMOMENTUM2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed y-momentum values
      IDX2(DtransformedDataAtNode,2,idx,0,0) =&
          YMOMENTUM2(DdataAtNode,IDX2,idx,0,0)
    end do

  end subroutine hydro_trafoNodalMomentum2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenEng2d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density and energy in 2D.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)

      ! Transformed total energy fluxes
      IDX3(DtransformedFluxesAtEdge,2,1,idx,0,0,0) =&
          TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,0,0,0) =&
         -TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0)
    end do

  end subroutine hydro_trafoFluxDenEng2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenEng2d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density and energy in 2D.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
          DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      ! Transformed total energy difference
      IDX2(DtransformedDataAtEdge,2,idx,0,0) =&
          TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
    end do

  end subroutine hydro_trafoDiffDenEng2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalDenEng2d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the density and
    ! energy in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DdataAtNode

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtNode
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, nnodes
      
      ! Transformed density values
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          DENSITY2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed energy values
      IDX2(DtransformedDataAtNode,2,idx,0,0) =&
          TOTALENERGY2(DdataAtNode,IDX2,idx,0,0)
    end do

  end subroutine hydro_trafoNodalDenEng2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenPre2d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density and energy in 2D.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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

    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Transformed density fluxes
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)
      
      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,2,1,idx,0,0,0) =&
          ((HYDRO_GAMMA)-1.0_DP)*(0.5_DP*(ui*ui+vi*vi)*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
                                  ui*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                                  vi*YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)+&
                                  TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
      IDX3(DtransformedFluxesAtEdge,2,2,idx,0,0,0) =&
         -((HYDRO_GAMMA)-1.0_DP)*(0.5_DP*(uj*uj+vj*vj)*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
                                  uj*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                                  vj*YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)+&
                                  TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
    end do

  end subroutine hydro_trafoFluxDenPre2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenPre2d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density and energy in 2D.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
          DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      
      ! Transformed pressure difference
      IDX2(DtransformedDataAtEdge,2,idx,0,0) =&
          PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
    end do
    
  end subroutine hydro_trafoDiffDenPre2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalDenPre2d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the density and
    ! pressure in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DdataAtNode

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtNode
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, nnodes
      
      ! Transformed density values
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          DENSITY2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed pressure values
      IDX2(DtransformedDataAtNode,2,idx,0,0) =&
          PRESSURE2(DdataAtNode,IDX2,idx,0,0)
    end do

  end subroutine hydro_trafoNodalDenPre2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenPreVel2d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation
    ! of the given flux into primitive variables in 2D.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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

    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Transformed density fluxes
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)

      ! Transformed velocity fluxes in x-direction
      IDX3(DtransformedFluxesAtEdge,2,1,idx,0,0,0) =&
          (XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
           ui*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0))/&
          DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,0,0,0) =&
         -(XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
           uj*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0))/&
          DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Transformed velocity fluxes in y-direction
      IDX3(DtransformedFluxesAtEdge,3,1,idx,0,0,0) =&
          (YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
           vi*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0))/&
          DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      IDX3(DtransformedFluxesAtEdge,3,2,idx,0,0,0) =&
         -(YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
           vj*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0))/&
          DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,4,1,idx,0,0,0) =&
          ((HYDRO_GAMMA)-1.0_DP)*(0.5_DP*(ui*ui+vi*vi)*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
                                  ui*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                                  vi*YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)+&
                                  TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
      IDX3(DtransformedFluxesAtEdge,4,2,idx,0,0,0) =&
         -((HYDRO_GAMMA)-1.0_DP)*(0.5_DP*(uj*uj+vj*vj)*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
                                  uj*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                                  vj*YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)+&
                                  TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
    end do

  end subroutine hydro_trafoFluxDenPreVel2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenPreVel2d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density, pressure
    ! and velocity in 2D.
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
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
      IDX2(DtransformedDataAtEdge,2,idx,0,0) =&
          DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      ! Transformed velocity difference in x-direction
      IDX2(DtransformedDataAtEdge,2,idx,0,0) =&
          XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      
      ! Transformed velocity difference in y-direction
      IDX2(DtransformedDataAtEdge,3,idx,0,0) =&
          YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      ! Transformed pressure difference
      IDX2(DtransformedDataAtEdge,4,idx,0,0) =&
          PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
    end do

  end subroutine hydro_trafoDiffDenPreVel2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalDenPreVel2d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the density,
    ! pressure and velocity in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DdataAtNode

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed solution values for all nodes under consideration
    !   DIMENSION(nvar,nnodes)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtNode
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, nnodes
      
      ! Transformed density values
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          DENSITY2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed x-velocity values
      IDX2(DtransformedDataAtNode,2,idx,0,0) =&
          XVELOCITY2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed y-velocity values
      IDX2(DtransformedDataAtNode,3,idx,0,0) =&
          YVELOCITY2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed pressure values
      IDX2(DtransformedDataAtNode,4,idx,0,0) =&
          PRESSURE2(DdataAtNode,IDX2,idx,0,0)
    end do

  end subroutine hydro_trafoNodalDenPreVel2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcBoundaryvalues2d(DbdrNormal, DpointNormal,&
      DbdrValue, ibdrCondType, Du, Du0, istatus)

!<description>
    ! This subroutine computes the boundary values for a given node in 2D.
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
      p   = ((HYDRO_GAMMA)-1.0_DP)*rho*(E-0.5_DP*(v1*v1+v2*v2))
      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))

      ! Precompute auxiliary data
      dnxy = DbdrNormal(1)*DbdrNormal(2)
      dnx2 = DbdrNormal(1)*DbdrNormal(1)
      dny2 = DbdrNormal(2)*DbdrNormal(2)

      if (ibdrCondType .eq. BDRC_FREESLIP) then
        ! Compute reflected velocities at the boundary
        v1_b = (dny2-dnx2)*v1 - 2.0_DP*dnxy*v2
        v2_b = (dnx2-dny2)*v2 - 2.0_DP*dnxy*v1
      else
        ! Compute initial velocity from previous time step
        v1_0 = Du0(2)/Du0(1)
        v2_0 = Du0(3)/Du0(1)

        ! Compute semi-reflected velocities at the boundary
        v1_b = (dny2-dnx2)*v1-2.0_DP*dnxy*v2 + 2*DbdrValue(1)*(v1_0*dnx2+v2_0*dnxy)
        v2_b = (dnx2-dny2)*v2-2.0_DP*dnxy*v1 + 2*DbdrValue(1)*(v2_0*dny2+v1_0*dnxy)
      end if

      ! Compute normal velocities at the boundary and the ghost state
      ! w.r.t. the numerical/approximate  outward unit normal vector
      vn   = DpointNormal(1)*v1   + DpointNormal(2)*v2
      vn_b = DpointNormal(1)*v1_b + DpointNormal(2)*v2_b

      ! Compute the tangential velocity depending on the sign of N*v
      if (vn .gt. 0.0_DP) then
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
      if (2.0_DP*(2.0_DP/((HYDRO_GAMMA)-1.0_DP))*c .le. vn_b-vn) then
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
      ppv  = p+0.5_DP*(vn-vn_b)*cup
      ppv  = max(0.0_DP, ppv)

      if (ppv .eq. p) then

        ! Select guessed pressure from PVRS Riemann solver
        pstar = ppv
      else
        if (ppv .lt. p) then

          ! Guess pressure from the Two-Rarefaction Riemann solver
          vm    = 0.5_DP*(vn+vn_b)
          ptl   = 1.0_DP + ((HYDRO_GAMMA)-1.0_DP)/2.0_DP*(vn-vm)/c
          ptr   = 1.0_DP + ((HYDRO_GAMMA)-1.0_DP)/2.0_DP*(vm-vn_b)/c
          pstar = 0.5_DP*(p*ptl + p*ptr)**(2.0_DP*(HYDRO_GAMMA)/((HYDRO_GAMMA)-1.0_DP))
        else

          ! Guess pressure from the Two-Shock Riemann solver
          ! with PVRS as estimated pressure value
          ge    = sqrt((2.0_DP/((HYDRO_GAMMA)+1.0_DP)/rho)/(((HYDRO_GAMMA)-1.0_DP)/((HYDRO_GAMMA)+1.0_DP)*p+ppv))
          pstar = p - 0.5_DP*(vn_b-vn)/ge
        end if
      end if

      ! Initialize solution difference and pressure
      vdiff = (vn_b-vn)/2.0_DP
      pold  = pstar

      newton: do ite = 1, 100

        ! Compute pressure function f(pold) and its derivative f1(pold)
        if (pold .le. p) then

          ! Rarefaction wave
          prat = pold/p

          f  = (2.0_DP/((HYDRO_GAMMA)-1.0_DP))*c*(prat**(((HYDRO_GAMMA)-1.0_DP)/(2.0_DP*(HYDRO_GAMMA))) - 1.0_DP)
          fd = (1.0_DP/(rho*c))*prat**(-((HYDRO_GAMMA)+1.0_DP)/(2.0_DP*(HYDRO_GAMMA)))
        else

          ! Shock wave
          auxA = 2.0_DP/((HYDRO_GAMMA)+1.0_DP)/rho
          auxB = ((HYDRO_GAMMA)-1.0_DP)/((HYDRO_GAMMA)+1.0_DP)*p
          qrt  = sqrt(auxA/(auxB + pold))

          f  = (pold-p)*qrt
          fd = (1.0_DP - 0.5_DP*(pold - p)/(auxB + pold))*qrt
        end if

        pstar = pold - (f+vdiff)/fd
        if (pstar .lt. 0.0_DP) then
          pold = 1.0E-6
          cycle newton
        end if

        aux = 2.0_DP*abs((pstar-pold)/(pstar+pold))
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
      vn = 0.5_DP*(vn+vn_b)


      !-------------------------------------------------------------------------
      ! Calculate the density in the star region
      !-------------------------------------------------------------------------

      if (pstar .le. p) then

        ! Rarefaction wave
        rho = rho*(pstar/p)**(1.0_DP/(HYDRO_GAMMA))
      else

        ! Shock wave
        rho = rho*(pstar/p+((HYDRO_GAMMA)-1.0_DP)/((HYDRO_GAMMA)+1.0_DP))/&
                  (((HYDRO_GAMMA)-1.0_DP)/((HYDRO_GAMMA)+1.0_DP)*(pstar/p)+1.0_DP)
      end if

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DpointNormal(2)*vt+DpointNormal(1)*vn)
      Du(3) = rho*(-DpointNormal(1)*vt+DpointNormal(2)*vn)
      Du(4) = pstar/((HYDRO_GAMMA)-1.0_DP)+0.5_DP*rho*(vn*vn+vt*vt)


    case(BDRC_VISCOUSWALL)
      !-------------------------------------------------------------------------

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      v2  = Du(3)/rho
      E   = Du(4)/rho
      p   = ((HYDRO_GAMMA)-1.0_DP)*rho*(E-0.5_DP*(v1*v1+v2*v2))

      ! Update the solution vector and let vn:=0 and vt:=0
      Du(2) = 0.0_DP
      Du(3) = 0.0_DP
      Du(4) = p/((HYDRO_GAMMA)-1.0_DP)


    case(BDRC_SUPERINLET)
      !-------------------------------------------------------------------------

      ! The free stream primitive variables are Deval=[rho,v1,v2,p]
      rho = DbdrValue(1)
      p   = DbdrValue(4)
      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))
      vn  = DbdrNormal(1)*DbdrValue(2)+DbdrNormal(2)*DbdrValue(3)
      vt  = DbdrNormal(2)*DbdrValue(2)-DbdrNormal(1)*DbdrValue(3)

      ! Compute Riemann invariants based on the free stream values
      W(1) = vn-2*c/((HYDRO_GAMMA)-1.0_DP)
      W(2) = vn+2*c/((HYDRO_GAMMA)-1.0_DP)
      W(3) = p/(rho**(HYDRO_GAMMA))
      W(4) = vt

      ! Transform back into conservative variables
      vn   = 0.5_DP*(W(1)+W(2))
      c    = 0.25*((HYDRO_GAMMA)-1.0_DP)*(W(2)-W(1))
      rho  = (c*c/(HYDRO_GAMMA)/W(3))**(1.0_DP/((HYDRO_GAMMA)-1.0_DP))
      p    = rho*c*c/(HYDRO_GAMMA)
      vt   = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/((HYDRO_GAMMA)-1.0_DP)+0.5_DP*rho*(vn*vn+vt*vt)


    case(BDRC_FREESTREAM)
      !-------------------------------------------------------------------------

      ! The free stream primitive variables are Deval=[rho,v1,v2,p]
      rho = DbdrValue(1)
      p   = DbdrValue(4)
      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))
      vn  = DbdrNormal(1)*DbdrValue(2)+DbdrNormal(2)*DbdrValue(3)
      vt  = DbdrNormal(2)*DbdrValue(2)-DbdrNormal(1)*DbdrValue(3)

      ! Compute Riemann invariants based on the free stream values
      Winf(1) = vn-2.0_DP*c/((HYDRO_GAMMA)-1.0_DP)
      Winf(2) = vn+2.0_DP*c/((HYDRO_GAMMA)-1.0_DP)
      Winf(3) = p/(rho**(HYDRO_GAMMA))
      Winf(4) = vt

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      v2  = Du(3)/rho
      E   = Du(4)/rho
      p   = ((HYDRO_GAMMA)-1.0_DP)*rho*(E-0.5_DP*(v1*v1+v2*v2))

      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))
      vn  = DbdrNormal(1)*v1+DbdrNormal(2)*v2
      vt  = DbdrNormal(2)*v1-DbdrNormal(1)*v2

      ! Compute Riemann invariants based on the solution values
      Wu(1) = vn-2.0_DP*c/((HYDRO_GAMMA)-1.0_DP)
      Wu(2) = vn+2.0_DP*c/((HYDRO_GAMMA)-1.0_DP)
      Wu(3) = p/(rho**(HYDRO_GAMMA))
      Wu(4) = vt

      ! Adopt free stream/computed values depending on the sign of the eigenvalue
      W(1) = merge(Winf(1), Wu(1), vn <  c)
      W(2) = merge(Winf(2), Wu(2), vn < -c)
      W(3) = merge(Winf(3), Wu(3), vn <  SYS_EPSREAL_DP)
      W(4) = merge(Winf(4), Wu(4), vn <  SYS_EPSREAL_DP)

      ! Transform back into conservative variables
      vn   = 0.5_DP*(W(1)+W(2))
      c    = 0.25*((HYDRO_GAMMA)-1.0_DP)*(W(2)-W(1))
      rho  = (c*c/(HYDRO_GAMMA)/W(3))**(1.0_DP/((HYDRO_GAMMA)-1.0_DP))
      p    = rho*c*c/(HYDRO_GAMMA)
      vt   = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/((HYDRO_GAMMA)-1.0_DP)+0.5_DP*rho*(vn*vn+vt*vt)


    case(BDRC_SUBINLET)
      !-------------------------------------------------------------------------

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      v2  = Du(3)/rho
      E   = Du(4)/rho
      p   = ((HYDRO_GAMMA)-1.0_DP)*rho*(E-0.5_DP*(v1*v1+v2*v2))

      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))
      vn  = DbdrNormal(1)*v1+DbdrNormal(2)*v2
      vt  = DbdrNormal(2)*v1-DbdrNormal(1)*v2

      ! The specified density and pressure is Deval=[rho,p]
      rho = DbdrValue(1)
      p   = DbdrValue(2)

      ! Compute Riemann invariants
      W(1) = vn-2.0_DP*c/((HYDRO_GAMMA)-1.0_DP)
      W(2) = vn+2.0_DP*c/((HYDRO_GAMMA)-1.0_DP)
      W(3) = p/(rho**(HYDRO_GAMMA))
      W(4) = vt

      ! Transform back into conservative variables
      vn   = 0.5_DP*(W(1)+W(2))
      c    = 0.25*((HYDRO_GAMMA)-1.0_DP)*(W(2)-W(1))
      rho  = (c*c/(HYDRO_GAMMA)/W(3))**(1.0_DP/((HYDRO_GAMMA)-1.0_DP))
      p    = rho*c*c/(HYDRO_GAMMA)
      vt   = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/((HYDRO_GAMMA)-1.0_DP)+0.5_DP*rho*(vn*vn+vt*vt)


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
      p   = ((HYDRO_GAMMA)-1.0_DP)*rho*(E-0.5_DP*(v1*v1+v2*v2))

      vn  = DbdrNormal(1)*v1+DbdrNormal(2)*v2
      vt  = DbdrNormal(2)*v1-DbdrNormal(1)*v2
      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))

      ! Compute Riemann invariants based on the solution values and prescribed exit pressure
      W(2) = 2*c/((HYDRO_GAMMA)-1.0_DP)-vn
      W(3) = p/(rho**(HYDRO_GAMMA))
      W(4) = vt
      W(1) = 4/((HYDRO_GAMMA)-1.0_DP)*sqrt(max((HYDRO_GAMMA)*ps/rho*(p/ps)**(1.0_DP/(HYDRO_GAMMA)), SYS_EPSREAL_DP))-W(2)

      ! Transform back into conservative variables
      vn  = 0.5_DP*(W(1)-W(2))
      c   = 0.25*((HYDRO_GAMMA)-1.0_DP)*(W(1)+W(2))
      rho = (c*c/(HYDRO_GAMMA)/W(3))**(1.0_DP/((HYDRO_GAMMA)-1.0_DP))
      p   = rho*c*c/(HYDRO_GAMMA)
      vt  = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/((HYDRO_GAMMA)-1.0_DP)+0.5_DP*rho*(vn*vn+vt*vt)


    case default
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
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   IquickAccess(3):     maximum number of expressions
    !   IquickAccess(4):     cubature rule
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
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
    type(t_boundaryRegion), pointer :: p_rboundaryRegionMirror
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(:), pointer :: Domega,DlocalData
    real(DP), dimension(:,:), pointer :: Daux,Dnx,Dny,DpointParMirror
    real(DP), dimension(:,:), pointer :: DcubPtsRef,Dbas,Dflux,Ddiff
    real(DP), dimension(:,:,:), pointer :: DstateI,DstateM,Dcoords
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dminParam,dmaxParam,dminParamMirror,dmaxParamMirror
    real(DP) :: dtime,dscale,cI,cM,dvnI,dvnM,dvtI,dvtM
    real(DP) :: hI,hM,l1,l2,l3,l4,pI,pM,rI,rM,uI,uM,vI,vM,w1,w2,w3,w4
    real(DP) :: aux,aux1,aux2,u_IM,v_IM,H_IM,vel_IM,q_IM,c_IM,c2_IM
    integer :: ibdrtype,isegment,nmaxExpr,ccubType
    integer :: iel,icubp,ipoint,npoints,ivar,nvar,iexpr,ivt,nve,neq
    

#ifndef HYDRO_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DHYDRO_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'hydro_coeffVectorBdr2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first quick access vector
    ! points to the solution vector
    p_rsolution => rcollection%p_rvectorQuickAccess1

#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
    ! Set pointer
    call lsysbl_getbase_double(p_rsolution, p_Ddata)
#endif

    ! Check if the solution is given in block or interleaved format
    if (p_rsolution%nblocks .eq. 1) then
      nvar = p_rsolution%RvectorBlock(1)%NVAR
      neq  = p_rsolution%NEQ/nvar
    else
      nvar = p_rsolution%nblocks
      neq  = p_rsolution%NEQ/nvar
    end if
    
    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first three quick access integer values hold:
    ! - the type of boundary condition
    ! - the segment number
    ! - the maximum number of expressions
    ! - the cubature rule
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)
    nmaxExpr = rcollection%IquickAccess(3)
    ccubType = rcollection%IquickAccess(4)
    
#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
    ! Evaluate one-dimensional basis functions on the boundary edge
    if (npointsPerElement .ne. cub_igetNumPts(ccubType)) then
      call output_line('Type of cubature rule at boundary mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_coeffVectorBdr2d_sim')
      call sys_halt()
    else
      ! How many DOFs are located at the boundary? This should be made
      ! more flexible by checking the type of element. For the time
      ! being, only linear and bilinear finite elements are supported
      npoints = 2
      
      ! How many vertices per element do we have?
      nve = elem_igetNVE(rdomainIntSubset%celement)
      
      ! Allocate temporal memory for one-dimensional
      ! cubature along the boundary edge
      allocate(Dbas(npoints,npointsPerElement))
      allocate(Domega(npointsPerElement))
      allocate(DcubPtsRef(1,npointsPerElement))

      ! Get the coordinates of the cubature points and the
      ! corresponding weights for the given cubature rule
      call cub_getCubature(ccubType, DcubPtsRef, Domega)
      
      ! Evaluate the one-dimensional basis functions 
      ! in the cubature points on the boundary
      do icubp = 1, npointsPerElement
        Dbas(1,icubp) = 0.5_DP*(1.0_DP-DcubPtsRef(1,icubp))
        Dbas(2,icubp) = 0.5_DP*(1.0_DP+DcubPtsRef(1,icubp))
      end do

      ! Deallocate temporal memory which is no longer required
      deallocate(DcubPtsRef,Domega)
    end if
#else
    ! Boundary values are evaluated directly at the cubature points
    npoints = npointsPerElement
#endif
    
    ! Allocate temporal memory for normal vectors, the coordinates and
    ! the solution state vectors in the DOFs located at the boundary
    allocate(Dnx(npoints,nelements), Dny(npoints,nelements))
    allocate(IDX3_ALLOCATE(DstateI,nvar,npoints,nelements))
    allocate(IDX3_ALLOCATE(DstateM,nvar,npoints,nelements))
#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
    allocate(Dcoords(NDIM2D,npoints,nelements))
#endif
    
    ! Get coordinates and internal state vector ...
    if (p_rsolution%nblocks .eq. 1) then

      ! ... for solutions stored in interleaved format

#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
      do iel = 1, nelements
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        ivt    = IdofsTest(ipoint,iel)

        ! Store internal state vector
        IDX3(DstateI,1,1,iel,0,0,0) = p_Ddata(nvar*(ivt-1)+1)
        IDX3(DstateI,2,1,iel,0,0,0) = p_Ddata(nvar*(ivt-1)+2)
        IDX3(DstateI,3,1,iel,0,0,0) = p_Ddata(nvar*(ivt-1)+3)
        IDX3(DstateI,4,1,iel,0,0,0) = p_Ddata(nvar*(ivt-1)+4)

        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)

        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store internal state vector
        IDX3(DstateI,1,2,iel,0,0,0) = p_Ddata(nvar*(ivt-1)+1)
        IDX3(DstateI,2,2,iel,0,0,0) = p_Ddata(nvar*(ivt-1)+2)
        IDX3(DstateI,3,2,iel,0,0,0) = p_Ddata(nvar*(ivt-1)+3)
        IDX3(DstateI,4,2,iel,0,0,0) = p_Ddata(nvar*(ivt-1)+4)

        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
      end do
#else
      ! Allocate temporal memory
      allocate(Daux(npointsPerElement*nvar, nelements))
      
      ! Evaluate the solution in the cubature points on the boundary
      call fevl_evaluate_sim(DER_FUNC2D, Daux, p_rsolution%RvectorBlock(1),&
          Dpoints, rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Distribute solution values to the internal state vector
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          IDX3(DstateI,1,ipoint,iel,0,0,0) = Daux((ipoint-1)*NVAR2D+1,iel)
          IDX3(DstateI,2,ipoint,iel,0,0,0) = Daux((ipoint-1)*NVAR2D+2,iel)
          IDX3(DstateI,3,ipoint,iel,0,0,0) = Daux((ipoint-1)*NVAR2D+3,iel)
          IDX3(DstateI,4,ipoint,iel,0,0,0) = Daux((ipoint-1)*NVAR2D+4,iel)
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Daux)
#endif

    else

      ! ... for solutions stored in block format

#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
      do iel = 1, nelements
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store internal state vector
        IDX3(DstateI,1,1,iel,0,0,0) = p_Ddata(      ivt)
        IDX3(DstateI,2,1,iel,0,0,0) = p_Ddata(neq  +ivt)
        IDX3(DstateI,3,1,iel,0,0,0) = p_Ddata(neq*2+ivt)
        IDX3(DstateI,4,1,iel,0,0,0) = p_Ddata(neq*3+ivt)

        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)

        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store internal state vector
        IDX3(DstateI,1,2,iel,0,0,0) = p_Ddata(      ivt)
        IDX3(DstateI,2,2,iel,0,0,0) = p_Ddata(neq  +ivt)
        IDX3(DstateI,3,2,iel,0,0,0) = p_Ddata(neq*2+ivt)
        IDX3(DstateI,4,2,iel,0,0,0) = p_Ddata(neq*3+ivt)

        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
      end do
#else
      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements))
      
      ! Evaluate the solution in the cubature points on the boundary
      do ivar = 1, nvar
        call fevl_evaluate_sim(DER_FUNC2D, Daux,&
            p_rsolution%RvectorBlock(ivar), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
        ! Distribute solution values to the internal state vector
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            IDX3(DstateI,ivar,ipoint,iel,0,0,0) = Daux(ipoint,iel)
          end do
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Daux)
#endif
    end if
    
#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
    ! Calculate the normal vectors in DOFs on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else
    ! Calculate the normal vectors in cubature on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))
      
    case (BDRC_FREESTREAM)
      !-----------------------------------------------------------------------
      ! Free-stream boundary conditions:
      !
      ! Compute the Riemann invariants based on the computed (internal)
      ! state vector and the given freestream state vector and select
      ! the Riemann invariant for each characteristic field based on the
      ! sign of the corresponding eigenvalue.
      
      ! Initialize values for function parser
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Set values for function parser
#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
          Dvalue(1:NDIM2D) = Dcoords(1:NDIM2D,ipoint,iel)
#else
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
#endif
          
          ! Compute free stream values from function parser given in
          ! term of the primitive variables [rho,v1,v2,p]
          do iexpr = 1, 4
            call fparser_evalFunction(p_rfparser,&
                nmaxExpr*(isegment-1)+iexpr,&
                Dvalue, IDX3(DstateM,iexpr,ipoint,iel,0,0,0))
          end do
          
          ! Compute auxiliary quantities based on free stream state vector
          rM = IDX3(DstateM,1,ipoint,iel,0,0,0)
          pM = IDX3(DstateM,4,ipoint,iel,0,0,0)
          cM = sqrt((HYDRO_GAMMA)*pM/rM)
          dvnM =  Dnx(ipoint,iel)*IDX3(DstateM,2,ipoint,iel,0,0,0)+&
                  Dny(ipoint,iel)*IDX3(DstateM,3,ipoint,iel,0,0,0)
          dvtM = -Dny(ipoint,iel)*IDX3(DstateM,2,ipoint,iel,0,0,0)+&
                  Dnx(ipoint,iel)*IDX3(DstateM,3,ipoint,iel,0,0,0)
          
          ! Compute auxiliary quantities based on internal state vector
          pI = PRESSURE3(DstateI,IDX3,ipoint,iel,0,0,0)
          cI = sqrt(max((HYDRO_GAMMA)*pI/&
               DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0), SYS_EPSREAL_DP))
          
          ! Compute the normal and tangential velocities based
          ! on internal state vector
          dvnI = ( Dnx(ipoint,iel)*&
                   XMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)+&
                   Dny(ipoint,iel)*&
                   YMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0))/&
                   DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0)
          dvtI = (-Dny(ipoint,iel)*&
                   XMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)+&
                   Dnx(ipoint,iel)*&
                   YMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0))/&
                   DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0)

          ! Select free stream or computed Riemann invariant depending
          ! on the sign of the corresponding eigenvalue
          if (dvnI .lt. cI) then
            w1 = dvnM-2.0_DP*cM/((HYDRO_GAMMA)-1.0_DP)
          else
            w1 = dvnI-2.0_DP*cI/((HYDRO_GAMMA)-1.0_DP)
          end if
          
          if (dvnI .lt. 0.0_DP) then
            w2 = pM/(rM**(HYDRO_GAMMA))
            w3 = dvtM
          else
            w2 = pI/DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0)**(HYDRO_GAMMA)
            w3 = dvtI
          end if
          
          if (dvnI .lt. -cI) then
            w4 = dvnM+2.0_DP*cM/((HYDRO_GAMMA)-1.0_DP)
          else
            w4 = dvnI+2.0_DP*cI/((HYDRO_GAMMA)-1.0_DP)
          end if
          
          ! Convert Riemann invariants into conservative state variables
          cM = 0.25*((HYDRO_GAMMA)-1.0_DP)*(w4-w1)
          rM = (cM*cM/(HYDRO_GAMMA)/w2)**(1.0_DP/((HYDRO_GAMMA)-1.0_DP))
          pM = rM*cM*cM/(HYDRO_GAMMA)
          dvnM = 0.5_DP*(w1+w4)
          dvtM = w3
          
          ! Calculate the state vector based on Riemann invariants
          IDX3(DstateM,1,ipoint,iel,0,0,0) = rM
          IDX3(DstateM,2,ipoint,iel,0,0,0) = rM*(Dnx(ipoint,iel)*dvnM-Dny(ipoint,iel)*dvtM)
          IDX3(DstateM,3,ipoint,iel,0,0,0) = rM*(Dny(ipoint,iel)*dvnM+Dnx(ipoint,iel)*dvtM)
          IDX3(DstateM,4,ipoint,iel,0,0,0) = pM/((HYDRO_GAMMA)-1.0_DP)+0.5_DP*rM*(dvnM**2+dvtM**2)
        end do
      end do

      
    case (BDRC_FREESLIP)
      !-----------------------------------------------------------------------
      ! Free-slip boundary condition:
      !
      ! Compute the mirrored state vector based on the values of the
      ! computed state vector and use an approximate Riemann solver
      
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Compute the normal and tangential velocities based
          ! on the internal state vector
          dvnI = ( Dnx(ipoint,iel)*&
                   XMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)+&
                   Dny(ipoint,iel)*&
                   YMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0))/&
                   DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0)
          dvtI = (-Dny(ipoint,iel)*&
                   XMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)+&
                   Dnx(ipoint,iel)*&
                   YMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0))/&
                   DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0)

          ! Compute the mirrored state vector
          IDX3(DstateM,1,ipoint,iel,0,0,0) = DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0)
          IDX3(DstateM,2,ipoint,iel,0,0,0) = DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0)*&
                                              (-dvnI*Dnx(ipoint,iel)-dvtI*Dny(ipoint,iel))
          IDX3(DstateM,3,ipoint,iel,0,0,0) = DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0)*&
                                              (-dvnI*Dny(ipoint,iel)+dvtI*Dnx(ipoint,iel))
          IDX3(DstateM,4,ipoint,iel,0,0,0) = TOTALENERGY3(DstateI,IDX3,ipoint,iel,0,0,0)
        end do
      end do


    case (BDRC_SUPERINLET)
      !-----------------------------------------------------------------------
      ! Supersonic inlet boundary conditions:
      !
      ! Prescribe the state vector in conservative variables
      
      ! Initialize values for function parser
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Set values for function parser
#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
          Dvalue(1:NDIM2D) = Dcoords(1:NDIM2D,ipoint,iel)
#else
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
#endif
          
          ! Compute boundary values from function parser given in
          ! term of the primitive variables [rho,v1,v2,p]
          do iexpr = 1, 4
            call fparser_evalFunction(p_rfparser,&
                nmaxExpr*(isegment-1)+iexpr,&
                Dvalue, IDX3(DstateM,iexpr,ipoint,iel,0,0,0))
          end do
          
          ! Compute convervative variables
          IDX3(DstateM,4,ipoint,iel,0,0,0) = IDX3(DstateM,4,ipoint,iel,0,0,0)/((HYDRO_GAMMA)-1.0_DP)&
              + IDX3(DstateM,1,ipoint,iel,0,0,0)*0.5_DP*(IDX3(DstateM,2,ipoint,iel,0,0,0)**2+&
                                                          IDX3(DstateM,3,ipoint,iel,0,0,0)**2)
          IDX3(DstateM,2,ipoint,iel,0,0,0) = IDX3(DstateM,1,ipoint,iel,0,0,0)*IDX3(DstateM,2,ipoint,iel,0,0,0)
          IDX3(DstateM,3,ipoint,iel,0,0,0) = IDX3(DstateM,1,ipoint,iel,0,0,0)*IDX3(DstateM,3,ipoint,iel,0,0,0)
        end do
      end do

        
    case (BDRC_SUPEROUTLET)
      !-----------------------------------------------------------------------
      ! Supersonic outlet boundary conditions:
      !
      ! Evaluate the boundary fluxes based on the computed state
      ! vector; since no Riemann problem is solved at the boundary we can
      ! treat this case in a special way and leave this routine immediately.
      
      ! Allocate temporal memory
      allocate(Dflux(nvar,npoints))
#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
      allocate(DlocalData(nvar))
#endif
      
      do iel = 1, nelements
        
        ! Loop over the DOFs and evaluate the Galerkin fluxes at DOFs
        do ipoint = 1, npoints
          
          ! Compute velocities and pressure
          uI = XVELOCITY3(DstateI,IDX3,ipoint,iel,0,0,0)
          vI = YVELOCITY3(DstateI,IDX3,ipoint,iel,0,0,0)
          pI = PRESSURE3(DstateI,IDX3,ipoint,iel,0,0,0)

          ! Calculate normal flux: ${\bf n}\cdot{\bf F}(U)$
          Dflux(1,ipoint) = Dnx(ipoint,iel)*&
                            XMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)&
                          + Dny(ipoint,iel)*&
                            YMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)
          Dflux(2,ipoint) = Dnx(ipoint,iel)*&
                            (XMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)*uI+pI)&
                          + Dny(ipoint,iel)*&
                            (XMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)*vI)
          Dflux(3,ipoint) = Dnx(ipoint,iel)*&
                            (YMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)*uI)&
                          + Dny(ipoint,iel)*&
                            (YMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)*vI+pI)
          Dflux(4,ipoint) = Dnx(ipoint,iel)*&
                            (TOTALENERGY3(DstateI,IDX3,ipoint,iel,0,0,0)+pI)*uI&
                          + Dny(ipoint,iel)*&
                            (TOTALENERGY3(DstateI,IDX3,ipoint,iel,0,0,0)+pI)*vI
        end do
        
#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Galerkin
        ! fluxes from the DOFs to the cubature points, where they are
        ! needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          DlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Galerkin fluxes
          do ipoint = 1, npoints
            DlocalData = DlocalData + Dbas(ipoint,icubp)*Dflux(:,ipoint)
          end do
          
          ! Store flux in the cubature points
          Dcoefficients(:,1,icubp,iel) = dscale*DlocalData
        end do
#else
        ! Loop over the cubature points and store the fluxes
        do ipoint = 1, npointsPerElement
          Dcoefficients(:,1,ipoint,iel) = dscale*Dflux(:,ipoint)
        end do
#endif
      end do
      
      ! Deallocate temporal memory
      deallocate(Dnx,Dny,DstateI,DstateM,Dflux)
#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
      deallocate(Dcoords,Dbas,DlocalData)
#endif

      ! That`s it
      return

      
    case (BDRC_SUBINLET)
      !-----------------------------------------------------------------------
      ! Subsonic pressure-density inlet boundary conditions:
      !
      ! Prescribe the density, pressure and tangential velocity at the inlet
      
      ! Initialize values for function parser
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
          
          ! Compute boundary values from function parser given in
          ! terms of the density, pressure and tangential velocity
          do iexpr = 1, 3
            call fparser_evalFunction(p_rfparser,&
                nmaxExpr*(isegment-1)+iexpr,&
                Dvalue, IDX3(DstateM,iexpr,ipoint,iel,0,0,0))
          end do
          
          ! Compute auxiliary quantities based on prescribed boundary values
          rM   = IDX3(DstateM,1,ipoint,iel,0,0,0)
          pM   = IDX3(DstateM,2,ipoint,iel,0,0,0)
          dvtM = IDX3(DstateM,3,ipoint,iel,0,0,0)
          cM   = sqrt((HYDRO_GAMMA)*pM/rM)
          
          ! Compute auxiliary quantities based on internal state vector
          pI = PRESSURE3(DstateI,IDX3,ipoint,iel,0,0,0)
          cI = sqrt(max((HYDRO_GAMMA)*pI/&
               DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0), SYS_EPSREAL_DP))
          
          ! Compute the normal velocity based on the internal state vector
          dvnI = ( Dnx(ipoint,iel)*&
                   XMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)+&
                   Dny(ipoint,iel)*&
                   YMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0))/&
                   DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0)
          
          ! Compute fourth Riemann invariant based on the internal state vector
          w4 = dvnI+2.0_DP*cI/((HYDRO_GAMMA)-1.0_DP)
          
          ! Compute the first Riemann invariant based on the fourth Riemann
          ! invariant and the prescribed boundary values
          w1 = w4-4.0*cM/((HYDRO_GAMMA)-1.0_DP)
          
          ! Compute the normal velocity based on the first and fourth Riemann
          ! invarient
          dvnM = 0.5_DP*(w1+w4)
          
          ! Setup the state vector based on Rimann invariants
          IDX3(DstateM,1,ipoint,iel,0,0,0) = rM
          IDX3(DstateM,2,ipoint,iel,0,0,0) = rM*(Dnx(ipoint,iel)*dvnM-Dny(ipoint,iel)*dvtM)
          IDX3(DstateM,3,ipoint,iel,0,0,0) = rM*(Dny(ipoint,iel)*dvnM+Dnx(ipoint,iel)*dvtM)
          IDX3(DstateM,4,ipoint,iel,0,0,0) = pM/((HYDRO_GAMMA)-1.0_DP)+0.5_DP*rM*(dvnM**2+dvtM**2)
        end do
      end do


    case (BDRC_SUBOUTLET)
      !-----------------------------------------------------------------------
      ! Subsonic pressure outlet boundary condition:
      !
      ! Prescribe the pressure at the outlet
      
      ! Initialize values for function parser
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Set values for function parser
#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
          Dvalue(1:NDIM2D) = Dcoords(1:NDIM2D,ipoint,iel)
#else
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
#endif
          
          ! Compute pressure value from function parser
          call fparser_evalFunction(p_rfparser,&
              nmaxExpr*(isegment-1)+1, Dvalue, pM)
          
          ! Compute auxiliary quantities based on internal state vector
          pI = PRESSURE3(DstateI,IDX3,ipoint,iel,0,0,0)
          cI = sqrt(max((HYDRO_GAMMA)*pI/&
               DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0), SYS_EPSREAL_DP))
          
          ! Compute the normal and tangential velocities based
          ! on internal state vector
          dvnI = ( Dnx(ipoint,iel)*&
                   XMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)+&
                   Dny(ipoint,iel)*&
                   YMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0))/&
                   DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0)
          dvtI = (-Dny(ipoint,iel)*&
                   XMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)+&
                   Dnx(ipoint,iel)*&
                   YMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0))/&
                   DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0)
          
          ! Compute three Riemann invariants based on internal state vector
          w2 = pI/DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0)**(HYDRO_GAMMA)
          w3 = dvtI
          w4 = dvnI+2.0_DP*cI/((HYDRO_GAMMA)-1.0_DP)
          
          ! Compute first Riemann invariant based on fourth Riemann invariant,
          ! the computed density and pressure and the prescribed exit pressure
          w1 = w4-4.0/((HYDRO_GAMMA)-1.0_DP)*sqrt((HYDRO_GAMMA)*pM/&
               DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0)*(pI/pM)**(1.0_DP/(HYDRO_GAMMA)))
          
          ! Convert Riemann invariants into conservative state variables
          cM = 0.25*((HYDRO_GAMMA)-1.0_DP)*(w4-w1)
          rM = (cM*cM/(HYDRO_GAMMA)/w2)**(1.0_DP/((HYDRO_GAMMA)-1.0_DP))
          pM = rM*cM*cM/(HYDRO_GAMMA)
          dvnM = 0.5_DP*(w1+w4)
          dvtM = w3
          
          ! Setup the state vector based on Riemann invariants
          IDX3(DstateM,1,ipoint,iel,0,0,0) = rM
          IDX3(DstateM,2,ipoint,iel,0,0,0) = rM*(Dnx(ipoint,iel)*dvnM-Dny(ipoint,iel)*dvtM)
          IDX3(DstateM,3,ipoint,iel,0,0,0) = rM*(Dny(ipoint,iel)*dvnM+Dnx(ipoint,iel)*dvtM)
          IDX3(DstateM,4,ipoint,iel,0,0,0) = pM/((HYDRO_GAMMA)-1.0_DP)+0.5_DP*rM*(dvnM**2+dvtM**2)
        end do
      end do


!!$      case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
!!$        !-----------------------------------------------------------------------
!!$        ! Periodic boundary conditions:
!!$        !
!!$        ! Compute the Riemann invariants based on the computed
!!$        ! (internal) state vector and on the state vector evaluated at
!!$        ! the mirror boundary and select the Riemann invariant for
!!$        ! each characteristic field based on the sign of the
!!$        ! corresponding eigenvalue.
!!$
!!$        ! Get mirrored boundary region from collection structure
!!$        p_rboundaryRegionMirror => collct_getvalue_bdreg(rcollection,&
!!$            'rboundaryRegionMirror')
!!$        
!!$        ! Get minimum/maximum parameter values from collection structure
!!$        dminParam = rcollection%DquickAccess(3)
!!$        dmaxParam = rcollection%DquickAccess(4)
!!$        dminParamMirror = rcollection%DquickAccess(5)
!!$        dmaxParamMirror = rcollection%DquickAccess(6)
!!$
!!$        ! Allocate temporal memory
!!$        allocate(DpointParMirror(npointsPerElement,nelements))
!!$        allocate(Daux3(npointsPerElement*nvar, nelements))
!!$
!!$        ! Rescale parameter values DpointPar on the boundary segment
!!$        ! where to compute the boundary conditions into parameter
!!$        ! values on the mirror boundary region
!!$        if (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_PERIODIC) then
!!$          call mprim_linearRescale(DpointPar, dminParam, dmaxParam,&
!!$              dmaxParamMirror, dminParamMirror, DpointParMirror)
!!$        else
!!$          call mprim_linearRescale(DpointPar, dminParam, dmaxParam,&
!!$              dmaxParamMirror, dminParamMirror, DpointParMirror)
!!$        end if
!!$
!!$        ! Evaluate the solution in the cubature points on the mirrored boundary
!!$        call doEvaluateAtBdrScalar(DER_FUNC2D, npointsPerElement*nelements*nvar,&
!!$            Daux3, p_rsolution%RvectorBlock(1), npointsPerElement*nelements,&
!!$            DpointParMirror, ibct, BDR_PAR_LENGTH, p_rboundaryRegionMirror)
!!$
!!$        do iel = 1, nelements
!!$          do ipoint = 1, npointsPerElement
!!$        
!!$            ! Compute auxiliary quantities based on the internal state
!!$            ! vector evaluated on the boundary
!!$            pI = ((HYDRO_GAMMA)-1.0_DP)*(Daux1((ipoint-1)*NVAR2D+4,iel)-&
!!$                              0.5_DP*(Daux1((ipoint-1)*NVAR2D+2,iel)**2+&
!!$                                   Daux1((ipoint-1)*NVAR2D+3,iel)**2)/&
!!$                                   Daux1((ipoint-1)*NVAR2D+1,iel))
!!$            cI = sqrt(max((HYDRO_GAMMA)*pI/Daux1((ipoint-1)*NVAR2D+1,iel), SYS_EPSREAL_DP))
!!$
!!$            ! Compute the normal and tangential velocities based on
!!$            ! the internal state vector evaluated on the boundary
!!$            dvnI = ( Dnx(ipoint,iel)*Daux1((ipoint-1)*NVAR2D+2,iel)+&
!!$                     Dny(ipoint,iel)*Daux1((ipoint-1)*NVAR2D+3,iel) )/&
!!$                     Daux1((ipoint-1)*NVAR2D+1,iel)
!!$            dvtI = (-Dny(ipoint,iel)*Daux1((ipoint-1)*NVAR2D+2,iel)+&
!!$                     Dnx(ipoint,iel)*Daux1((ipoint-1)*NVAR2D+3,iel) )/&
!!$                     Daux1((ipoint-1)*NVAR2D+1,iel)
!!$
!!$            ! Compute auxiliary quantities based on state vector
!!$            ! evaluated on the mirrored boundary
!!$            pM = ((HYDRO_GAMMA)-1.0_DP)*(Daux3((ipoint-1)*NVAR2D+4,iel)-&
!!$                              0.5_DP*(Daux3((ipoint-1)*NVAR2D+2,iel)**2+&
!!$                                   Daux3((ipoint-1)*NVAR2D+3,iel)**2)/&
!!$                                   Daux3((ipoint-1)*NVAR2D+1,iel))
!!$            cM = sqrt(max((HYDRO_GAMMA)*pM/Daux3((ipoint-1)*NVAR2D+1,iel), SYS_EPSREAL_DP))
!!$
!!$            ! Compute the normal and tangential velocities based on
!!$            ! state vector evaluated on the mirrored boundary
!!$            dvnM = ( Dnx(ipoint,iel)*Daux3((ipoint-1)*NVAR2D+2,iel)+&
!!$                     Dny(ipoint,iel)*Daux3((ipoint-1)*NVAR2D+3,iel) )/&
!!$                     Daux3((ipoint-1)*NVAR2D+1,iel)
!!$            dvtM = (-Dny(ipoint,iel)*Daux3((ipoint-1)*NVAR2D+2,iel)+&
!!$                     Dnx(ipoint,iel)*Daux3((ipoint-1)*NVAR2D+3,iel) )/&
!!$                     Daux3((ipoint-1)*NVAR2D+1,iel)
!!$
!!$            ! Select internal or mirrored Riemann invariant depending
!!$            ! on the sign of the corresponding eigenvalue
!!$            if (dvnI .lt. cI) then
!!$              DstateM(1) = dvnM-2.0_DP*cM/((HYDRO_GAMMA)-1.0_DP)
!!$            else
!!$              DstateM(1) = dvnI-2.0_DP*cI/((HYDRO_GAMMA)-1.0_DP)
!!$            end if
!!$
!!$            if (dvnI .lt. SYS_EPSREAL_DP) then
!!$              DstateM(2) = pM/(Daux3((ipoint-1)*NVAR2D+1,iel)**(HYDRO_GAMMA))
!!$              DstateM(3) = dvtM
!!$            else
!!$              DstateM(2) = pI/(Daux1((ipoint-1)*NVAR2D+1,iel)**(HYDRO_GAMMA))
!!$              DstateM(3) = dvtI
!!$            end if
!!$
!!$            if (dvnI .lt. -cI) then
!!$              DstateM(4) = dvnM+2.0_DP*cM/((HYDRO_GAMMA)-1.0_DP)
!!$            else
!!$              DstateM(4) = dvnI+2.0_DP*cI/((HYDRO_GAMMA)-1.0_DP)
!!$            end if
!!$            
!!$            ! Convert Riemann invariants into conservative state variables
!!$            cM = 0.25*((HYDRO_GAMMA)-1.0_DP)*(DstateM(4)-DstateM(1))
!!$            rM = (cM*cM/(HYDRO_GAMMA)/DstateM(2))**(1.0_DP/((HYDRO_GAMMA)-1.0_DP))
!!$            pM = rM*cM*cM/(HYDRO_GAMMA)
!!$            dvnM = 0.5_DP*(DstateM(1)+DstateM(4))
!!$            dvtM = DstateM(3)
!!$
!!$            ! Setup the state vector based on Riemann invariants
!!$            DstateM(1) = rM
!!$            DstateM(2) = rM*(Dnx(ipoint,iel)*dvnM-Dny(ipoint,iel)*dvtM)
!!$            DstateM(3) = rM*(Dny(ipoint,iel)*dvnM+Dnx(ipoint,iel)*dvtM)
!!$            DstateM(4) = pM/((HYDRO_GAMMA)-1.0_DP) + 0.5_DP*rM*(dvnM**2+dvtM**2)
!!$            
!!$            ! Setup the computed internal state vector
!!$            DstateI(1) = Daux1((ipoint-1)*NVAR2D+1,iel)
!!$            DstateI(2) = Daux1((ipoint-1)*NVAR2D+2,iel)
!!$            DstateI(3) = Daux1((ipoint-1)*NVAR2D+3,iel)
!!$            DstateI(4) = Daux1((ipoint-1)*NVAR2D+4,iel)
!!$            
!!$            ! Invoke Riemann solver
!!$            call doRiemannSolver(DstateI, DstateM,&
!!$                Dnx(ipoint,iel), Dny(ipoint,iel), Dflux, Ddiff)
!!$            
!!$            ! Store flux in the cubature points
!!$            Dcoefficients(:,1,ipoint,iel) = dscale*0.5_DP*(Dflux-Ddiff)
!!$          end do
!!$        end do
!!$
!!$        ! Deallocate temporal memory
!!$        deallocate(DpointParMirror, Daux3)

    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_coeffVectorBdr2d_sim')
      call sys_halt()
      
    end select


    ! Allocate temporal memory
    allocate(Dflux(nvar,npoints), Ddiff(nvar,npoints))
#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
    allocate(DlocalData(nvar))
#endif
    
    do iel = 1, nelements
      
      ! Loop over the DOFs and evaluate the Galerkin fluxes at DOFs
      do ipoint = 1, npoints
        
        !-----------------------------------------------------------------------
        ! Solve the boundary Riemann problem by Roe`s approximate Riemann solver
        !-----------------------------------------------------------------------

        ! Compute velocities and pressure from internal state
        uI = XVELOCITY3(DstateI,IDX3,ipoint,iel,0,0,0)
        vI = YVELOCITY3(DstateI,IDX3,ipoint,iel,0,0,0)
        pI = PRESSURE3(DstateI,IDX3,ipoint,iel,0,0,0)
        rI = DENSITY3(DstateI,IDX3,ipoint,iel,0,0,0)
        hI = (TOTALENERGY3(DstateI,IDX3,ipoint,iel,0,0,0)+pI)/rI
        
        ! Compute velocities and pressure from mirrored state
        uM = XVELOCITY3(DstateM,IDX3,ipoint,iel,0,0,0)
        vM = YVELOCITY3(DstateM,IDX3,ipoint,iel,0,0,0)
        pM = PRESSURE3(DstateM,IDX3,ipoint,iel,0,0,0)
        rM = DENSITY3(DstateM,IDX3,ipoint,iel,0,0,0)
        hM = (TOTALENERGY3(DstateM,IDX3,ipoint,iel,0,0,0)+pM)/rM

        ! Calculate normal flux: $\frac12{\bf n}\cdot[{\bf F}(U_I)+{\bf F}(U_M)]$
        Dflux(1,ipoint) = Dnx(ipoint,iel)*&
                          (XMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)+&
                           XMOMENTUM3(DstateM,IDX3,ipoint,iel,0,0,0))&
                        + Dny(ipoint,iel)*&
                          (YMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)+&
                           YMOMENTUM3(DstateM,IDX3,ipoint,iel,0,0,0))
        Dflux(2,ipoint) = Dnx(ipoint,iel)*&
                          (XMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)*uI+pI+&
                           XMOMENTUM3(DstateM,IDX3,ipoint,iel,0,0,0)*uM+pM)&
                        + Dny(ipoint,iel)*&
                          (XMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)*vI+&
                           XMOMENTUM3(DstateM,IDX3,ipoint,iel,0,0,0)*vM)
        Dflux(3,ipoint) = Dnx(ipoint,iel)*&
                          (YMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)*uI+&
                           YMOMENTUM3(DstateM,IDX3,ipoint,iel,0,0,0)*uM)&
                        + Dny(ipoint,iel)*&
                          (YMOMENTUM3(DstateI,IDX3,ipoint,iel,0,0,0)*vI+pI+&
                           YMOMENTUM3(DstateM,IDX3,ipoint,iel,0,0,0)*vM+pM)
        Dflux(4,ipoint) = Dnx(ipoint,iel)*&
                          ((TOTALENERGY3(DstateI,IDX3,ipoint,iel,0,0,0)+pI)*uI+&
                           (TOTALENERGY3(DstateM,IDX3,ipoint,iel,0,0,0)+pM)*uM)&
                        + Dny(ipoint,iel)*&
                          ((TOTALENERGY3(DstateI,IDX3,ipoint,iel,0,0,0)+pI)*vI+&
                           (TOTALENERGY3(DstateM,IDX3,ipoint,iel,0,0,0)+pM)*vM)

        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(rI,rM)
        u_IM = ROE_MEAN_VALUE(uI,uM,aux)
        v_IM = ROE_MEAN_VALUE(vI,vM,aux)
        H_IM = ROE_MEAN_VALUE(hI,hM,aux)
      
        ! Compute auxiliary variables
        vel_IM = Dnx(ipoint,iel)*u_IM + Dny(ipoint,iel)*v_IM
        q_IM   = 0.5_DP*(u_IM*u_IM+v_IM*v_IM)

        ! Compute the speed of sound
        c2_IM = max(((HYDRO_GAMMA)-1.0_DP)*(H_IM-q_IM), SYS_EPSREAL_DP)
        c_IM  = sqrt(c2_IM)
        
        ! Compute eigenvalues
        l1 = abs(vel_IM-c_IM)
        l2 = abs(vel_IM)
        l3 = abs(vel_IM+c_IM)
        l4 = abs(vel_IM)
      
        ! Compute solution difference U_M-U_I
        Ddiff(:,ipoint) = IDX3(DstateM,:,ipoint,iel,0,0,0)-&
                          IDX3(DstateI,:,ipoint,iel,0,0,0)
      
        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-1.0_DP)*(q_IM*Ddiff(1,ipoint)&
                                      -u_IM*Ddiff(2,ipoint)&
                                      -v_IM*Ddiff(3,ipoint)&
                                           +Ddiff(4,ipoint))/2.0_DP/c2_IM
        aux2 =        (vel_IM*Ddiff(1,ipoint)&
             -Dnx(ipoint,iel)*Ddiff(2,ipoint)&
             -Dny(ipoint,iel)*Ddiff(3,ipoint))/2.0_DP/c_IM
      
        ! Compute characteristic variables multiplied by the
        ! corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((1.0_DP-((HYDRO_GAMMA)-1.0_DP)*q_IM/c2_IM)*Ddiff(1,ipoint)&
                                +((HYDRO_GAMMA)-1.0_DP)*(u_IM*Ddiff(2,ipoint)&
                                                        +v_IM*Ddiff(3,ipoint)&
                                                             -Ddiff(4,ipoint))/c2_IM)
        w3 = l3 * (aux1 - aux2)
        w4 = l4 * ((Dnx(ipoint,iel)*v_IM-Dny(ipoint,iel)*u_IM)*Ddiff(1,ipoint)&
                                              +Dny(ipoint,iel)*Ddiff(2,ipoint)&
                                              -Dnx(ipoint,iel)*Ddiff(3,ipoint))
      
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        Ddiff(1,ipoint) = w1 + w2 + w3
        Ddiff(2,ipoint) = (u_IM-c_IM*Dnx(ipoint,iel))*w1 + u_IM*w2 +&
                          (u_IM+c_IM*Dnx(ipoint,iel))*w3 + Dny(ipoint,iel)*w4
        Ddiff(3,ipoint) = (v_IM-c_IM*Dny(ipoint,iel))*w1 + v_IM*w2 +&
                          (v_IM+c_IM*Dny(ipoint,iel))*w3 - Dnx(ipoint,iel)*w4
        Ddiff(4,ipoint) = (H_IM-c_IM*vel_IM)*w1 + q_IM*w2 + (H_IM+c_IM*vel_IM)*w3 +&
                          (u_IM*Dny(ipoint,iel)-v_IM*Dnx(ipoint,iel))*w4
      end do
        
#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
      ! Loop over the cubature points and interpolate the Galerkin
      ! fluxes from the DOFs to the cubature points, where they are
      ! needed by the linear form assembly routine
      do icubp = 1, npointsPerElement
        
        DlocalData = 0.0_DP
        
        ! Loop over the DOFs and interpolate the Galerkin fluxes
        do ipoint = 1, npoints
          DlocalData = DlocalData + Dbas(ipoint,icubp)*0.5_DP*(Dflux(:,ipoint)-Ddiff(:,ipoint))
        end do
        
        ! Store flux in the cubature points
        Dcoefficients(:,1,icubp,iel) = dscale*DlocalData
      end do
#else
      ! Loop over the cubature points and store the fluxes
      do ipoint = 1, npointsPerElement
        Dcoefficients(:,1,ipoint,iel) = dscale*0.5_DP*(Dflux(:,ipoint)-Ddiff(:,ipoint))
      end do
#endif
    end do

    ! Deallocate temporal memory
    deallocate(Dnx,Dny,DstateI,DstateM,Dflux,Ddiff)
#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
    deallocate(Dcoords,Dbas,DlocalData)
#endif

  contains
    
    !***************************************************************************
    ! Evaluate th solution vector at some boundary points given in
    ! terms of their parameter values. This ugly trick is necessary
    ! since we have to pass the 2d-array Dvalues and DpointsPar to a
    ! subroutine which accepts only 1d-arrays.
    ! ***************************************************************************

    subroutine doEvaluateAtBdrScalar(iderType, n, Dvalues, rvectorScalar,&
        m, DpointsPar, ibdc, cparType, rboundaryRegion)
      
      integer, intent(in) :: iderType,ibdc,cparType,n,m
      real(DP), dimension(m), intent(in) :: DpointsPar
      type(t_vectorScalar), intent(in) :: rvectorScalar
      type(t_boundaryRegion), intent(in) :: rboundaryRegion
      
      real(DP), dimension(n), intent(out) :: Dvalues
      
      call fevl_evaluateBdr2D(iderType, Dvalues, rvectorScalar,&
          DpointsPar, ibdc, cparType, rboundaryRegion)
      
    end subroutine doEvaluateAtBdrScalar

    !***************************************************************************
    ! Evaluate th solution vector at some boundary points given in
    ! terms of their parameter values. This ugly trick is necessary
    ! since we have to pass the 2d-array Dvalues and DpointsPar to a
    ! subroutine which accepts only 1d-arrays.
    ! ***************************************************************************

    subroutine doEvaluateAtBdrBlock(iderType, n1, n2, Dvalues, rvectorBlock,&
        m, DpointsPar, ibdc, cparType, rboundaryRegion)
      
      integer, intent(in) :: iderType,ibdc,cparType,n1,n2,m
      real(DP), dimension(m), intent(in) :: DpointsPar
      type(t_vectorBlock), intent(in) :: rvectorBlock
      type(t_boundaryRegion), intent(in) :: rboundaryRegion
      
      real(DP), dimension(n1,1,n2), intent(out) :: Dvalues
      
      call fevl_evaluateBdr2D((/iderType/), Dvalues, rvectorBlock,&
          DpointsPar, ibdc, cparType, rboundaryRegion)
      
    end subroutine doEvaluateAtBdrBlock
    
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
    ! A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   IquickAccess(1):     NEQ or ivt
    !   IquickAccess(2:5):   ivt1,...,ivt5
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
            0.5_DP*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR2D+ivar)+&
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
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR2D+ivar) = 0.0_DP
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case default
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
    ! A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   IquickAccess(1):     NEQ or ivt
    !   IquickAccess(2:5):   ivt1,...,ivt5
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
            0.5_DP*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
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
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = 0.0_DP
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)

    end select

  end subroutine hydro_hadaptCallbackBlock2d

end module hydro_callback2d
