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
!# 1.) hydro_calcFluxGalerkin3d_sim
!#     -> Computes fluxes for standard Galerkin scheme
!#
!# 2.) hydro_calcFluxGalNoBdr3d_sim
!#     -> Computes fluxes for standard Galerkin scheme without
!#        assembling the symmetric boundary contribution
!#
!# 3.) hydro_calcFluxScDiss3d_sim
!#     -> Computes fluxes for low-order discretisation adopting
!#        scalar artificial viscosities
!#
!# 4.) hydro_calcFluxScDissDiSp3d_sim
!#     -> Computes fluxes for low-order discretisation adopting
!#        scalar artificial viscosities based on dimensional
!#        splitting approach
!#
!# 5.) hydro_calcFluxRoeDiss3d_sim
!#     -> Computes fluxes for low-order discretisation adopting
!#        tensorial artificial viscosities of Roe-type
!#
!# 6.) hydro_calcFluxRoeDissDiSp3d_sim
!#     -> Computes fluxes for low-order discretisation adopting
!#        tensorial artificial viscosities of Roe-type based on
!#        dimensional splitting approach
!#
!# 7.) hydro_calcFluxRusDiss3d_sim
!#     -> Computes fluxes for low-order discretisation adopting
!#        scalar artificial diffusion of Rusanov-type
!#
!# 8.) hydro_calcFluxRusDissDiSp3d_sim
!#     -> Computes fluxes for low-order discretisation adopting
!#        scalar artificial diffusion of Rusanov-type based on
!#        dimensional splitting approach
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
!# 12.) hydro_calcMatGalerkin3d_sim
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
!#         adopting tensorial artificial viscosities of Roe-type
!#
!# 16.) hydro_calcMatRoeDiss3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities of Roe-type
!#
!# 17.) hydro_calcMatRusDissMatD3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities of Rusanov-type
!#
!# 18.) hydro_calcMatRusDiss3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities of Rusanov-type
!#
!# 19.) hydro_calcCharacteristics3d_sim
!#      -> Computes characteristic variables
!#
!# 20.) hydro_calcFluxFCTScDiss3d_sim
!#      -> Computes fluxes for FCT algorithm adopting scalar
!#         artificial viscosities
!#
!# 21.) hydro_calcFluxFCTRoeDiss3d_sim
!#      -> Computes fluxes for FCT algorithm adopting tensorial
!#         artificial viscosities of Roe-type
!#
!# 22.) hydro_calcFluxFCTRusDiss3d_sim
!#      -> Computes fluxes for FCT algorithm adopting scalar
!#         artificial viscosities of Rusanov-type
!#
!# 23.) hydro_trafoFluxDensity3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 24.) hydro_trafoDiffDensity3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 25. hydro_trafoNodalDensity3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density values
!#
!# 26.) hydro_trafoFluxEnergy3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 27.) hydro_trafoDiffEnergy3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 28. hydro_trafoNodalEnergy3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal energy values
!#
!# 29.) hydro_trafoFluxPressure3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 30.) hydro_trafoDiffPressure3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the pressure
!#
!# 31. hydro_trafoNodalPressure3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal pressure values
!#
!# 32.) hydro_trafoFluxVelocity3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 33.) hydro_trafoDiffVelocity3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 34.) hydro_trafoNodalVelocity3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal velocity values
!#
!# 35.) hydro_trafoFluxMomentum3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 36.) hydro_trafoDiffMomentum3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 37.) hydro_trafoNodalMomentum3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal momentum values
!#
!# 38.) hydro_trafoFluxDenEng3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 39.) hydro_trafoDiffDenEng3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 40.) hydro_trafoNodalDenEng3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density and energy values
!#
!# 41.) hydro_trafoFluxDenPre3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 42.) hydro_trafoDiffDenPre3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 43.) hydro_trafoNodalDenPre3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density and pressure values
!#
!# 44.) hydro_trafoFluxDenPreVel3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 45.) hydro_trafoDiffDenPreVel3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure
!#         and the velocity
!#
!# 46.) hydro_trafoNodalDenPreVel3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density, pressure and velocity values
!#
!# 47.) hydro_calcBoundaryvalues3d
!#      -> Computes the boundary values for a given node
!#
!# </purpose>
!##############################################################################

module hydro_callback3d

#include "../../flagship.h"
#include "hydro.h"

!$use omp_lib
  use basicgeometry
  use collection
  use fsystem
  use genoutput
  use graph
  use groupfemsystem
  use linearsystemblock
  use linearsystemscalar
  use problem
  use solveraux
  use storage

  ! Modules from hydrodynamic model
  use hydro_basic

  implicit none

  private

  public :: hydro_calcFluxGalerkin3d_sim
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
  public :: hydro_calcMatGalerkin3d_sim
  public :: hydro_calcMatScDissMatD3d_sim
  public :: hydro_calcMatScDiss3d_sim
  public :: hydro_calcMatRoeDissMatD3d_sim
  public :: hydro_calcMatRoeDiss3d_sim
  public :: hydro_calcMatRusDissMatD3d_sim
  public :: hydro_calcMatRusDiss3d_sim
  public :: hydro_calcCharacteristics3d_sim
  public :: hydro_calcFluxFCTScDiss3d_sim
  public :: hydro_calcFluxFCTRoeDiss3d_sim
  public :: hydro_calcFluxFCTRusDiss3d_sim
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
  public :: hydro_trafoNodalDensity3d_sim
  public :: hydro_trafoNodalEnergy3d_sim
  public :: hydro_trafoNodalPressure3d_sim
  public :: hydro_trafoNodalVelocity3d_sim
  public :: hydro_trafoNodalMomentum3d_sim
  public :: hydro_trafoNodalDenEng3d_sim
  public :: hydro_trafoNodalDenPre3d_sim
  public :: hydro_trafoNodalDenPreVel3d_sim
  public :: hydro_calcBoundaryvalues3d

  ! CUDA wrapper routines
  public :: hydro_calcDivVecGalerkin3d_cuda
  public :: hydro_calcDivVecScDiss3d_cuda
  public :: hydro_calcDivVecScDissDiSp3d_cuda
  public :: hydro_calcDivVecRoeDiss3d_cuda
  public :: hydro_calcDivVecRoeDissDiSp3d_cuda
  public :: hydro_calcDivVecRusDiss3d_cuda
  public :: hydro_calcDivVecRusDissDiSp3d_cuda

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGalerkin3d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

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
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IedgeList

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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP) :: pi,pj,ui,uj,vi,vj,wi,wj
    integer :: idx

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)

      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)

      ! Assemble skew-symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj+&
           IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fzj-&
           IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi-&
           IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fzi )
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)

      ! Assemble fluxes
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fz_ij)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fz_ij)
#endif
    end do

  end subroutine hydro_calcFluxGalerkin3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGalNoBdr3d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

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
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IedgeList

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
  real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
  real(DP) :: pi,pj,ui,uj,vi,vj,wi,wj
  integer :: idx
  
    
    do idx = 1, nedges

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)

      ! Assemble symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
          (RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                        IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*Fx_ij+&
           RCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                        IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*Fy_ij+&
           RCONST(0.5)*(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)-&
                        IDX3(DcoeffsAtEdge,3,2,idx,_,_,_))*Fz_ij)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
    end do

  end subroutine hydro_calcFluxGalNoBdr3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxScDiss3d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)
    
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
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IedgeList

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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,anorm,aux,c_ij,d_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    integer :: idx

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)

      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      ! Compute densities
      ri = DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      rj = DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute enthalpies
      hi = (TOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
      hj = (TOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(ri,rj)
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      w_ij = ROE_MEAN_VALUE(wi,wj,aux)
      H_ij = ROE_MEAN_VALUE(hi,hj,aux)

      ! Compute auxiliary variables
      vel_ij = u_ij*a(1) + v_ij*a(2) + w_ij*a(3)
      q_ij   = RCONST(0.5)*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

      ! Compute the speed of sound
      c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))

      ! Compute scalar dissipation
      d_ij = abs(vel_ij) + anorm*c_ij

      ! Multiply the solution difference by the scalar dissipation
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,_,_,_)-&
                   IDX3(DdataAtEdge,:,1,idx,_,_,_))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj+&
           IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fzj-&
           IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi-&
           IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fzi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fz_ij+ Diff)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fz_ij+ Diff)
#endif
    end do

  end subroutine hydro_calcFluxScDiss3d_sim

  !***************************************************************************

!<subroutine>

  subroutine hydro_calcDivVecGalerkin3d_cuda(rgroupFEMSet, rx, ry, dscale,&
      bclear, fcb_calcFluxSys_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine computes the fluxes for the standard Galerkin
    ! discretisation in 3D.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: callback function to compute local fluxes
    include '../../../../../kernel/PDEOperators/intf_calcFluxSys_sim.inc'
    optional :: fcb_calcFluxSys_sim
!</input>

!<inputoutput>
    ! Destination vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT
    
    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: IEDGEmax,IEDGEset,igroup
    type(C_PTR) :: cptr_DcoeffsAtEdge
    type(C_PTR) :: cptr_IedgeList
    type(C_PTR) :: cptr_Dx, cptr_Dy
    integer(I64) :: istream

    ! Create CUDA stream
    call coproc_createStream(istream)
    
    ! Get pointers to memory blocks on device
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IedgeList, cptr_IedgeList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtEdge, cptr_DcoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    if (.not.storage_isAssociated(cptr_Dx)) then
      call lsysbl_copyH2D_Vector(rx, .false., .false., istream)
      call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    end if
   
    ! Make sure that the destination vector ry exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Vector(ry, .true., .false., istream)
    call storage_getMemPtrOnDevice(ry%h_Ddata, cptr_Dy)
   
    ! Set pointer
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
    
    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronise memory access
    do igroup = 1, size(p_IedgeListIdx)-1
      
      ! Do nothing for empty groups
      if (p_IedgeListIdx(igroup+1)-p_IedgeListIdx(igroup) .le. 0) cycle

      ! Get position of first edge in group
      IEDGEset = p_IedgeListIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IedgeListIdx(igroup+1)-1
      
      ! Use callback function to compute internodal fluxes
      call hydro_calcFluxGalerkin3d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
          cptr_Dx, cptr_Dy, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
          rgroupFEMSet%NEDGE, rgroupFEMSet%ncoeffsAtEdge, IEDGEmax-IEDGEset+1,&
          IEDGEset, istream)
    end do

    ! Transfer destination vector back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Vector(ry, bclear, .false., istream)

    ! Ensure data consistency
    call coproc_synchronizeStream(istream)
    call coproc_destroyStream(istream)

#else

    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivVecGalerkin3d_cuda')
    call sys_halt()
    
#endif

  end subroutine hydro_calcDivVecGalerkin3d_cuda

  !***************************************************************************

!<subroutine>

  subroutine hydro_calcDivVecScDiss3d_cuda(rgroupFEMSet, rx, ry, dscale,&
      bclear, fcb_calcFluxSys_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 3D using scalar artificial viscosities proportional to the
    ! spectral radius (largest eigenvalue) of the Roe-matrix.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: callback function to compute local fluxes
    include '../../../../../kernel/PDEOperators/intf_calcFluxSys_sim.inc'
    optional :: fcb_calcFluxSys_sim
!</input>

!<inputoutput>
    ! Destination vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT
    
    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: IEDGEmax,IEDGEset,igroup
    type(C_PTR) :: cptr_DcoeffsAtEdge
    type(C_PTR) :: cptr_IedgeList
    type(C_PTR) :: cptr_Dx, cptr_Dy
    integer(I64) :: istream

    ! Create CUDA stream
    call coproc_createStream(istream)
    
    ! Get pointers to memory blocks on device
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IedgeList, cptr_IedgeList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtEdge, cptr_DcoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    if (.not.storage_isAssociated(cptr_Dx)) then
      call lsysbl_copyH2D_Vector(rx, .false., .false., istream)
      call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    end if
   
    ! Make sure that the destination vector ry exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Vector(ry, .true., .false., istream)
    call storage_getMemPtrOnDevice(ry%h_Ddata, cptr_Dy)
   
    ! Set pointer
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
    
    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronise memory access
    do igroup = 1, size(p_IedgeListIdx)-1
      
      ! Do nothing for empty groups
      if (p_IedgeListIdx(igroup+1)-p_IedgeListIdx(igroup) .le. 0) cycle

      ! Get position of first edge in group
      IEDGEset = p_IedgeListIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IedgeListIdx(igroup+1)-1
      
      ! Use callback function to compute internodal fluxes
      call hydro_calcFluxScDiss3d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
          cptr_Dx, cptr_Dy, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
          rgroupFEMSet%NEDGE, rgroupFEMSet%ncoeffsAtEdge, IEDGEmax-IEDGEset+1,&
          IEDGEset, istream)
    end do

    ! Transfer destination vector back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Vector(ry, bclear, .false., istream)

    ! Ensure data consistency
    call coproc_synchronizeStream(istream)
    call coproc_destroyStream(istream)

#else

    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivVecScDiss3d_cuda')
    call sys_halt()
    
#endif

  end subroutine hydro_calcDivVecScDiss3d_cuda

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxScDissDiSp3d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)
    

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
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IedgeList

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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,aux,c_ij,d_ij,q_ij,u_ij,v_ij,w_ij
    integer :: idx

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)

      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      
      ! Compute densities
      ri = DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      rj = DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute enthalpies
      hi = (TOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
      hj = (TOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(ri,rj)
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      w_ij = ROE_MEAN_VALUE(wi,wj,aux)
      H_ij = ROE_MEAN_VALUE(hi,hj,aux)
      
      ! Compute auxiliary variable
      q_ij = RCONST(0.5)*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

      ! Compute the speed of sound
      c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))
      
      ! Compute scalar dissipation with dimensional splitting
      d_ij = ( abs(a(1)*u_ij) + abs(a(1))*c_ij +&
               abs(a(2)*v_ij) + abs(a(2))*c_ij +&
               abs(a(3)*w_ij) + abs(a(3))*c_ij )

      ! Multiply the solution difference by the scalar dissipation
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,_,_,_)-&
                   IDX3(DdataAtEdge,:,1,idx,_,_,_))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj+&
           IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fzj-&
           IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi-&
           IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fzi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fz_ij+ Diff)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fz_ij+ Diff)
#endif
    end do

  end subroutine hydro_calcFluxScDissDiSp3d_sim

  !***************************************************************************

!<subroutine>

  subroutine hydro_calcDivVecScDissDiSp3d_cuda(rgroupFEMSet, rx, ry, dscale,&
      bclear, fcb_calcFluxSys_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 3D using scalar artificial viscosities proportional to the
    ! spectral radius (largest eigenvalue) of the Roe-matrix, whereby
    ! dimensional splitting is employed.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: callback function to compute local fluxes
    include '../../../../../kernel/PDEOperators/intf_calcFluxSys_sim.inc'
    optional :: fcb_calcFluxSys_sim
!</input>

!<inputoutput>
    ! Destination vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT
    
    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: IEDGEmax,IEDGEset,igroup
    type(C_PTR) :: cptr_DcoeffsAtEdge
    type(C_PTR) :: cptr_IedgeList
    type(C_PTR) :: cptr_Dx, cptr_Dy
    integer(I64) :: istream

    ! Create CUDA stream
    call coproc_createStream(istream)
    
    ! Get pointers to memory blocks on device
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IedgeList, cptr_IedgeList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtEdge, cptr_DcoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    if (.not.storage_isAssociated(cptr_Dx)) then
      call lsysbl_copyH2D_Vector(rx, .false., .false., istream)
      call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    end if
   
    ! Make sure that the destination vector ry exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Vector(ry, .true., .false., istream)
    call storage_getMemPtrOnDevice(ry%h_Ddata, cptr_Dy)
   
    ! Set pointer
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
    
    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronise memory access
    do igroup = 1, size(p_IedgeListIdx)-1
      
      ! Do nothing for empty groups
      if (p_IedgeListIdx(igroup+1)-p_IedgeListIdx(igroup) .le. 0) cycle

      ! Get position of first edge in group
      IEDGEset = p_IedgeListIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IedgeListIdx(igroup+1)-1
      
      ! Use callback function to compute internodal fluxes
      call hydro_calcFluxScDissDiSp3d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
          cptr_Dx, cptr_Dy, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
          rgroupFEMSet%NEDGE, rgroupFEMSet%ncoeffsAtEdge, IEDGEmax-IEDGEset+1,&
          IEDGEset, istream)
    end do

    ! Transfer destination vector back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Vector(ry, bclear, .false., istream)

    ! Ensure data consistency
    call coproc_synchronizeStream(istream)
    call coproc_destroyStream(istream)

#else

    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivVecScDissDiSp3d_cuda')
    call sys_halt()
    
#endif

  end subroutine hydro_calcDivVecScDissDiSp3d_cuda

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRoeDiss3d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

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
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IedgeList

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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,c2_ij,c_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    real(DP) :: anorm,aux,aux1,aux2
    real(DP) :: l1,l2,l3,l4,l5,w1,w2,w3,w4,w5
    integer :: idx
#if defined(HYDRO_USE_ENTROPYFIX) && (HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX)
    real(DP) :: dtol,ci,cj,veli,velj
#endif

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)

      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Normalise the skew-symmetric coefficient
        a = a/anorm
        
        ! Compute densities
        ri = DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute enthalpies
        hi = (TOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)

        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2) + w_ij*a(3)
        q_ij   = RCONST(0.5)*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

        ! Compute the speed of sound
        c2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij  = sqrt(c2_ij)

        ! Compute eigenvalues
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)
        l5 = abs(vel_ij)

#if defined(HYDRO_USE_ENTROPYFIX)

#if HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX

        ! Entropy-fix by Harten and Hyman
        ci = SOUNDSPEED3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        cj = SOUNDSPEED3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
        veli = ui*a(1) + vi*a(2) + wi*a(3)
        velj = uj*a(1) + vj*a(2) + wj*a(3)

        dtol = max(RCONST(0.0), (vel_ij-c_ij) - (veli-ci), (velj-cj) - (vel_ij-c_ij) )
        
        if (l1 .lt. dtol)&
            l1 = RCONST(0.5)*((l1*l1)/dtol + dtol)
        
        dtol = max(RCONST(0.0), (vel_ij+c_ij) - (veli+ci), (velj+cj) - (vel_ij+c_ij) )

        if (l3 .lt. dtol)&
            l3 = RCONST(0.5)*((l3*l3)/dtol + dtol)

#elif HYDRO_USE_ENTROPYFIX == HARTEN_ENTROPYFIX

#ifndef HYDRO_HARTEN_ENTROPYFIX
#error "Value HYDRO_HARTEN_ENTROPYFIX is required!"
#else
        ! Entropy-fix by Harten
        if (l1 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l1 = RCONST(0.5)*((l1*l1)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))

        if (l3 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l3 = RCONST(0.5)*((l3*l3)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))
#endif
#else
#error "Invalid type of entropy fix!"
#endif
#endif

        ! Compute solution difference U_j-U_i
        Diff = IDX3(DdataAtEdge,:,2,idx,_,_,_)-&
               IDX3(DdataAtEdge,:,1,idx,_,_,_)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-RCONST(1.0))/RCONST(2.0)/c2_ij*(q_ij*Diff(1)&
                                                             -u_ij*Diff(2)&
                                                             -v_ij*Diff(3)&
                                                             -w_ij*Diff(4)&
                                                                  +Diff(5))
        aux2 = (vel_ij*Diff(1)&
                 -a(1)*Diff(2)&
                 -a(2)*Diff(3)&
                 -a(3)*Diff(4))/RCONST(2.0)/c_ij

        ! Get the dimension with largest coefficient
        select case(maxloc(a,1))
        case(1)
          ! Compute characteristic variables multiplied by the corresponding eigenvalue
          w1 = l1 * (aux1 + aux2)
          w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*Diff(1)&
                                      +((HYDRO_GAMMA)-RCONST(1.0))*( u_ij*Diff(2)&
                                                                    +v_ij*Diff(3)&
                                                                    +w_ij*Diff(4)&
                                                                         -Diff(5))/c2_ij)
          w3 = l3 * (aux1 - aux2)
          w4 = l4 * ( (v_ij-vel_ij*a(2))/a(1)*Diff(1)&
                                        +a(2)*Diff(2)&
                +(a(2)*a(2)-RCONST(1.0))/a(1)*Diff(3)&
                              +a(2)*a(3)/a(1)*Diff(4))
          w5 = l5 * ( (vel_ij*a(3)-w_ij)/a(1)*Diff(1)&
                                        -a(3)*Diff(2)&
                              -a(2)*a(3)/a(1)*Diff(3)&
                +(RCONST(1.0)-a(3)*a(3))/a(1)*Diff(4))

          ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
          Diff(1) = anorm * ( w1 + w2 + w3 )
          Diff(2) = anorm * ( (u_ij-c_ij*a(1))*w1 + u_ij*w2 +&
                              (u_ij+c_ij*a(1))*w3 + a(2)*w4 - a(3)*w5 )
          Diff(3) = anorm * ( (v_ij-c_ij*a(2))*w1 + v_ij*w2 +&
                              (v_ij+c_ij*a(2))*w3 - a(1)*w4 )
          Diff(4) = anorm * ( (w_ij-c_ij*a(3))*w1 + w_ij*w2 +&
                              (w_ij+c_ij*a(3))*w3 + a(1)*w5 )
          Diff(5) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3&
                              + (u_ij*a(2)-v_ij*a(1))*w4 + (w_ij*a(1)-u_ij*a(3))*w5 )

        case(2)
          ! Compute characteristic variables multiplied by the corresponding eigenvalue
          w1 = l1 * (aux1 + aux2)
          w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*Diff(1)&
                                       +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*Diff(2)&
                                                                    +v_ij*Diff(3)&
                                                                    +w_ij*Diff(4)&
                                                                         -Diff(5))/c2_ij)
          w3 = l3 * (aux1 - aux2)
          w4 = l4 * ( (vel_ij*a(1)-u_ij)/a(2)*Diff(1)&
                +(RCONST(1.0)-a(1)*a(1))/a(2)*Diff(2)&
                                        -a(1)*Diff(3)&
                              -a(1)*a(3)/a(2)*Diff(4))
          w5 = l5 * ( (w_ij-vel_ij*a(3))/a(2)*Diff(1)&
                              +a(1)*a(3)/a(2)*Diff(2)&
                                        +a(3)*Diff(3)&
                +(a(3)*a(3)-RCONST(1.0))/a(2)*Diff(4))

          ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
          Diff(1) = anorm * ( w1 + w2 + w3 )
          Diff(2) = anorm * ( (u_ij-c_ij*a(1))*w1 + u_ij*w2 +&
                              (u_ij+c_ij*a(1))*w3 + a(2)*w4 )
          Diff(3) = anorm * ( (v_ij-c_ij*a(2))*w1 + v_ij*w2 +&
                              (v_ij+c_ij*a(2))*w3 - a(1)*w4 + a(3)*w5 )
          Diff(4) = anorm * ( (w_ij-c_ij*a(3))*w1 + w_ij*w2 +&
                              (w_ij+c_ij*a(3))*w3 - a(2)*w5 )
          Diff(5) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3&
                              + (u_ij*a(2)-v_ij*a(1))*w4 + (v_ij*a(3)-w_ij*a(2))*w5 )

        case(3)
          ! Compute characteristic variables multiplied by the corresponding eigenvalue
          w1 = l1 * (aux1 + aux2)
          w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*Diff(1)&
                                       +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*Diff(2)&
                                                                    +v_ij*Diff(3)&
                                                                    +w_ij*Diff(4)&
                                                                         -Diff(5))/c2_ij)
          w3 = l3 * (aux1 - aux2)
          w4 = l4 * ( (u_ij-vel_ij*a(1))/a(3)*Diff(1)&
                +(a(1)*a(1)-RCONST(1.0))/a(3)*Diff(2)&
                              +a(1)*a(2)/a(3)*Diff(3)&
                                        +a(1)*Diff(4) )
          w5 = l5 * ( (vel_ij*a(2)-v_ij)/a(3)*Diff(1)&
                              -a(1)*a(2)/a(3)*Diff(2)&
                +(RCONST(1.0)-a(2)*a(2))/a(3)*Diff(3)&
                                        -a(2)*Diff(4) )

          ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
          Diff(1) = anorm * ( w1 + w2 + w3 )
          Diff(2) = anorm * ( (u_ij-c_ij*a(1))*w1 + u_ij*w2 +&
                              (u_ij+c_ij*a(1))*w3 - a(3)*w4 )
          Diff(3) = anorm * ( (v_ij-c_ij*a(2))*w1 + v_ij*w2 +&
                              (v_ij+c_ij*a(2))*w3 + a(3)*w5 )
          Diff(4) = anorm * ( (w_ij-c_ij*a(3))*w1 + w_ij*w2 +&
                              (w_ij+c_ij*a(3))*w3 + a(1)*w4 - a(2)*w5)
          Diff(5) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3&
                              + (w_ij*a(1)-u_ij*a(3))*w4 + (v_ij*a(3)-w_ij*a(2))*w5 )
        end select

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
             IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj+&
             IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fzj-&
             IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
             IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi-&
             IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fzi + Diff)
        IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
        IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
            (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fz_ij + Diff)
        IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fz_ij + Diff)
#endif
      else
#ifdef HYDRO_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
             IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj+&
             IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fzj-&
             IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
             IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi-&
             IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fzi)
        IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
        IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
            (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fz_ij)
        IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fz_ij)
#endif
      end if
    end do

  end subroutine hydro_calcFluxRoeDiss3d_sim

  !***************************************************************************

!<subroutine>

  subroutine hydro_calcDivVecRoeDiss3d_cuda(rgroupFEMSet, rx, ry, dscale, bclear,&
      fcb_calcFluxSys_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 3D using scalar artificial viscosities of Roe-type.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
    
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: callback function to compute local fluxes
    include '../../../../../kernel/PDEOperators/intf_calcFluxSys_sim.inc'
    optional :: fcb_calcFluxSys_sim
!</input>

!<inputoutput>
    ! Destination vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT
    
    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: IEDGEmax,IEDGEset,igroup
    type(C_PTR) :: cptr_DcoeffsAtEdge
    type(C_PTR) :: cptr_IedgeList
    type(C_PTR) :: cptr_Dx, cptr_Dy
    integer(I64) :: istream

    ! Create CUDA stream
    call coproc_createStream(istream)
    
    ! Get pointers to memory blocks on device
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IedgeList, cptr_IedgeList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtEdge, cptr_DcoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    if (.not.storage_isAssociated(cptr_Dx)) then
      call lsysbl_copyH2D_Vector(rx, .false., .false., istream)
      call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    end if
   
    ! Make sure that the destination vector ry exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Vector(ry, .true., .false., istream)
    call storage_getMemPtrOnDevice(ry%h_Ddata, cptr_Dy)
   
    ! Set pointer
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
    
    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronise memory access
    do igroup = 1, size(p_IedgeListIdx)-1
      
      ! Do nothing for empty groups
      if (p_IedgeListIdx(igroup+1)-p_IedgeListIdx(igroup) .le. 0) cycle

      ! Get position of first edge in group
      IEDGEset = p_IedgeListIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IedgeListIdx(igroup+1)-1
      
      ! Use callback function to compute internodal fluxes
      call hydro_calcFluxRoeDiss3d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
          cptr_Dx, cptr_Dy, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
          rgroupFEMSet%NEDGE, rgroupFEMSet%ncoeffsAtEdge, IEDGEmax-IEDGEset+1,&
          IEDGEset, istream)
    end do

    ! Transfer destination vector back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Vector(ry, bclear, .false., istream)

    ! Ensure data consistency
    call coproc_synchronizeStream(istream)
    call coproc_destroyStream(istream)

#else

    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivVecRoeDiss3d_cuda')
    call sys_halt()
    
#endif

  end subroutine hydro_calcDivVecRoeDiss3d_cuda

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRoeDissDiSp3d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)
    

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
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IedgeList

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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: DiffX,DiffY,DiffZ
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,c2_ij,c_ij,q_ij,u_ij,v_ij,w_ij
    real(DP) :: anorm,aux,aux1,aux2
    real(DP) :: l1,l2,l3,l4,l5,w1,w2,w3,w4,w5
    integer :: idx
#if defined(HYDRO_USE_ENTROPYFIX) && (HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX)
    real(DP) :: dtol,ci,cj
#endif

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)

      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
#endif
      
      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Compute the absolute value
        a = abs(a)
        
        ! Compute densities
        ri = DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute enthalpies
        hi = (TOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)

        ! Compute auxiliary variable
        q_ij = RCONST(0.5)*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

        ! Compute the speed of sound
        c2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij  = sqrt(c2_ij)

        !-----------------------------------------------------------------------
        ! Dimensional splitting: x-direction
        !-----------------------------------------------------------------------

        ! Compute eigenvalues
        l1 = abs(u_ij-c_ij)
        l2 = abs(u_ij)
        l3 = abs(u_ij+c_ij)
        l4 = abs(u_ij)
        l5 = abs(u_ij)
        
#if defined(HYDRO_USE_ENTROPYFIX)

#if HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX

        ! Entropy-fix by Harten and Hyman
        ci = SOUNDSPEED3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        cj = SOUNDSPEED3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

        dtol = max(RCONST(0.0), (u_ij-c_ij) - (ui-ci), (uj-cj) - (u_ij-c_ij) )
        
        if (l1 .lt. dtol)&
            l1 = RCONST(0.5)*((l1*l1)/dtol + dtol)
        
        dtol = max(RCONST(0.0), (u_ij+c_ij) - (ui+ci), (uj+cj) - (u_ij+c_ij) )

        if (l3 .lt. dtol)&
            l3 = RCONST(0.5)*((l3*l3)/dtol + dtol)

#elif HYDRO_USE_ENTROPYFIX == HARTEN_ENTROPYFIX

#ifndef HYDRO_HARTEN_ENTROPYFIX
#error "Value HYDRO_HARTEN_ENTROPYFIX is required!"
#else
        ! Entropy-fix by Harten
        if (l1 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l1 = RCONST(0.5)*((l1*l1)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))

        if (l3 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l3 = RCONST(0.5)*((l3*l3)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))
#endif
#else
#error "Invalid type of entropy fix!"
#endif
#endif

        ! Compute solution difference U_j-U_i
        DiffX = IDX3(DdataAtEdge,:,2,idx,_,_,_)-&
                IDX3(DdataAtEdge,:,1,idx,_,_,_)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-RCONST(1.0))/RCONST(2.0)/c2_ij*(q_ij*DiffX(1)&
                                                             -u_ij*DiffX(2)&
                                                             -v_ij*DiffX(3)&
                                                             -w_ij*DiffX(4)&
                                                                  +DiffX(5))
        aux2 = (u_ij*DiffX(1)-DiffX(2))/RCONST(2.0)/c_ij

        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*DiffX(1)&
                                    +((HYDRO_GAMMA)-RCONST(1.0))*( u_ij*DiffX(2)&
                                                                  +v_ij*DiffX(3)&
                                                                  +w_ij*DiffX(4)&
                                                                       -DiffX(5))/c2_ij)
        w3 = l3 * (aux1 - aux2)
        w4 = l4 * ( v_ij*DiffX(1)-DiffX(3))
        w5 = l5 * (-w_ij*DiffX(1)+DiffX(4))

        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        DiffX(1) = a(1) * ( w1 + w2 + w3 )
        DiffX(2) = a(1) * ( (u_ij-c_ij)*w1 + u_ij*w2 + (u_ij+c_ij)*w3 )
        DiffX(3) = a(1) * ( v_ij*(w1 + w2 + w3) - w4 )
        DiffX(4) = a(1) * ( w_ij*(w1 + w2 + w3) + w5 )
        DiffX(5) = a(1) * ( (H_ij-c_ij*u_ij)*w1 + q_ij*w2 + (H_ij+c_ij*u_ij)*w3&
                             -v_ij*w4 + w_ij*w5 )

        !-----------------------------------------------------------------------
        ! Dimensional splitting: y-direction
        !-----------------------------------------------------------------------

        ! Compute eigenvalues
        l1 = abs(v_ij-c_ij)
        l2 = abs(v_ij)
        l3 = abs(v_ij+c_ij)
        l4 = abs(v_ij)
        l5 = abs(v_ij)

#if defined(HYDRO_USE_ENTROPYFIX)

#if HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX

        ! Entropy-fix by Harten and Hyman
        dtol = max(RCONST(0.0), (v_ij-c_ij) - (vi-ci), (vj-cj) - (v_ij-c_ij) )
        
        if (l1 .lt. dtol)&
            l1 = RCONST(0.5)*((l1*l1)/dtol + dtol)
        
        dtol = max(RCONST(0.0), (v_ij+c_ij) - (vi+ci), (vj+cj) - (v_ij+c_ij) )

        if (l3 .lt. dtol)&
            l3 = RCONST(0.5)*((l3*l3)/dtol + dtol)

#elif HYDRO_USE_ENTROPYFIX == HARTEN_ENTROPYFIX

#ifndef HYDRO_HARTEN_ENTROPYFIX
#error "Value HYDRO_HARTEN_ENTROPYFIX is required!"
#else
        ! Entropy-fix by Harten
        if (l1 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l1 = RCONST(0.5)*((l1*l1)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))

        if (l3 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l3 = RCONST(0.5)*((l3*l3)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))
#endif
#else
#error "Invalid type of entropy fix!"
#endif
#endif

        ! Compute solution difference U_j-U_i
        DiffY = IDX3(DdataAtEdge,:,2,idx,_,_,_)-&
                IDX3(DdataAtEdge,:,1,idx,_,_,_)

        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-RCONST(1.0))/RCONST(2.0)/c2_ij*(q_ij*DiffY(1)&
                                                             -u_ij*DiffY(2)&
                                                             -v_ij*DiffY(3)&
                                                             -w_ij*DiffY(4)&
                                                                  +DiffY(5))
        aux2 = (v_ij*DiffY(1)-DiffY(3))/RCONST(2.0)/c_ij

        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*DiffY(1)&
                                     +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*DiffY(2)&
                                                                  +v_ij*DiffY(3)&
                                                                  +w_ij*DiffY(4)&
                                                                       -DiffY(5))/c2_ij)
        w3 = l3 * (aux1 - aux2)
        w4 = l4 * (-u_ij*DiffY(1)+DiffY(2))
        w5 = l5 * ( w_ij*DiffY(1)-DiffY(4))

        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        DiffY(1) = a(2) * ( w1 + w2 + w3 )
        DiffY(2) = a(2) * ( u_ij*(w1 + w2 + w3) + w4 )
        DiffY(3) = a(2) * ( (v_ij-c_ij)*w1 + v_ij*w2 + (v_ij+c_ij)*w3 )
        DiffY(4) = a(2) * ( w_ij*(w1 + w2 + w3) - w5 )
        DiffY(5) = a(2) * ( (H_ij-c_ij*v_ij)*w1 + q_ij*w2 + (H_ij+c_ij*v_ij)*w3&
                            + u_ij*w4 -w_ij*w5 )

        !-----------------------------------------------------------------------
        ! Dimensional splitting: z-direction
        !-----------------------------------------------------------------------
        
        ! Compute eigenvalues
        l1 = abs(w_ij-c_ij)
        l2 = abs(w_ij)
        l3 = abs(w_ij+c_ij)
        l4 = abs(w_ij)
        l5 = abs(w_ij)

#if defined(HYDRO_USE_ENTROPYFIX)

#if HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX

        ! Entropy-fix by Harten and Hyman
        dtol = max(RCONST(0.0), (w_ij-c_ij) - (wi-ci), (wj-cj) - (w_ij-c_ij) )
        
        if (l1 .lt. dtol)&
            l1 = RCONST(0.5)*((l1*l1)/dtol + dtol)
        
        dtol = max(RCONST(0.0), (w_ij+c_ij) - (wi+ci), (wj+cj) - (w_ij+c_ij) )

        if (l3 .lt. dtol)&
            l3 = RCONST(0.5)*((l3*l3)/dtol + dtol)

#elif HYDRO_USE_ENTROPYFIX == HARTEN_ENTROPYFIX

#ifndef HYDRO_HARTEN_ENTROPYFIX
#error "Value HYDRO_HARTEN_ENTROPYFIX is required!"
#else
        ! Entropy-fix by Harten
        if (l1 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l1 = RCONST(0.5)*((l1*l1)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))

        if (l3 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l3 = RCONST(0.5)*((l3*l3)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))
#endif
#else
#error "Invalid type of entropy fix!"
#endif
#endif

        ! Compute solution difference U_j-U_i
        DiffZ = IDX3(DdataAtEdge,:,2,idx,_,_,_)-&
                IDX3(DdataAtEdge,:,1,idx,_,_,_)

        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-RCONST(1.0))/RCONST(2.0)/c2_ij*(q_ij*DiffZ(1)&
                                                             -u_ij*DiffZ(2)&
                                                             -v_ij*DiffZ(3)&
                                                             -w_ij*DiffZ(4)&
                                                                  +DiffZ(5))
        aux2 = (w_ij*DiffZ(1)-DiffZ(3))/RCONST(2.0)/c_ij

        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*DiffZ(1)&
                                     +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*DiffZ(2)&
                                                                  +v_ij*DiffZ(3)&
                                                                  +w_ij*DiffZ(4)&
                                                                       -DiffZ(5))/c2_ij)
        w3 = l3 * (aux1 - aux2)
        w4 = l4 * ( u_ij*DiffZ(1)-DiffZ(2))
        w5 = l5 * (-v_ij*DiffZ(1)+DiffZ(3))
        
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        DiffZ(1) = a(3) * ( w1 + w2 + w3 )
        DiffZ(2) = a(3) * ( u_ij*(w1 + w2 + w3) - w4 )
        DiffZ(3) = a(3) * ( v_ij*(w1 + w2 + w3) + w5 )
        DiffZ(4) = a(3) * ( (w_ij-c_ij)*w1 + w_ij*w2 + (w_ij+c_ij)*w3 )
        DiffZ(5) = a(3) * ( (H_ij-c_ij*w_ij)*w1 + q_ij*w2 + (H_ij+c_ij*w_ij)*w3&
                            -u_ij*w4 + v_ij*w5 )

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
             IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj+&
             IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fzj-&
             IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
             IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi-&
             IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fzi+&
             DiffX+DiffY+DiffZ)
        IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
        IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
            (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fz_ij+&
             DiffX+DiffY+DiffZ)
        IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fz_ij+&
             DiffX+DiffY+DiffZ)
#endif
      else
 #ifdef HYDRO_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
             IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj+&
             IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fzj-&
             IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
             IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi-&
             IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fzi)
        IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
        IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
            (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fz_ij)
        IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fz_ij)
#endif
      end if
    end do

  end subroutine hydro_calcFluxRoeDissDiSp3d_sim

  !***************************************************************************

!<subroutine>

  subroutine hydro_calcDivVecRoeDissDiSp3d_cuda(rgroupFEMSet, rx, ry, dscale,&
      bclear, fcb_calcFluxSys_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 3D using scalar artificial viscosities of Roe-type, whereby
    ! dimensional splitting is employed.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: callback function to compute local fluxes
    include '../../../../../kernel/PDEOperators/intf_calcFluxSys_sim.inc'
    optional :: fcb_calcFluxSys_sim
!</input>

!<inputoutput>
    ! Destination vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT
    
    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: IEDGEmax,IEDGEset,igroup
    type(C_PTR) :: cptr_DcoeffsAtEdge
    type(C_PTR) :: cptr_IedgeList
    type(C_PTR) :: cptr_Dx, cptr_Dy
    integer(I64) :: istream

    ! Create CUDA stream
    call coproc_createStream(istream)
    
    ! Get pointers to memory blocks on device
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IedgeList, cptr_IedgeList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtEdge, cptr_DcoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    if (.not.storage_isAssociated(cptr_Dx)) then
      call lsysbl_copyH2D_Vector(rx, .false., .false., istream)
      call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    end if
   
    ! Make sure that the destination vector ry exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Vector(ry, .true., .false., istream)
    call storage_getMemPtrOnDevice(ry%h_Ddata, cptr_Dy)
   
    ! Set pointer
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
    
    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronise memory access
    do igroup = 1, size(p_IedgeListIdx)-1
      
      ! Do nothing for empty groups
      if (p_IedgeListIdx(igroup+1)-p_IedgeListIdx(igroup) .le. 0) cycle

      ! Get position of first edge in group
      IEDGEset = p_IedgeListIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IedgeListIdx(igroup+1)-1
      
      ! Use callback function to compute internodal fluxes
      call hydro_calcFluxRoeDissDiSp3d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
          cptr_Dx, cptr_Dy, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
          rgroupFEMSet%NEDGE, rgroupFEMSet%ncoeffsAtEdge, IEDGEmax-IEDGEset+1,&
          IEDGEset, istream)
    end do

    ! Transfer destination vector back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Vector(ry, bclear, .false., istream)

    ! Ensure data consistency
    call coproc_synchronizeStream(istream)
    call coproc_destroyStream(istream)

#else

    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivVecRoeDissDiSp3d_cuda')
    call sys_halt()
    
#endif

  end subroutine hydro_calcDivVecRoeDissDiSp3d_cuda

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRusDiss3d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)
    
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
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IedgeList

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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP) :: Ei,Ej,ci,cj,pi,pj,ui,uj,vi,vj,wi,wj
    real(DP) :: d_ij
    integer :: idx

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)

      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !---------------------------------------------------------------------------
      
      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*&
          (Ei-RCONST(0.5)*(ui*ui+vi*vi+wi*wi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*&
          (Ej-RCONST(0.5)*(uj*uj+vj*vj+wj*wj)), SYS_EPSREAL_DP))

#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      ! (COMPARE WITH 2D CASE!!!)
      d_ij = max( abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*uj+&
                      RCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*vj+&
                      RCONST(0.5)*(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,3,2,idx,_,_,_))*wj)+&
              RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)+&
                               POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2)+&
                               POW(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),2))*cj,&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_))*ui+&
                      RCONST(0.5)*(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_))*vi+&
                      RCONST(0.5)*(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,3,1,idx,_,_,_))*wi)+&
              RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)+&
                               POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2)+&
                               POW(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),2))*ci )
#else
      ! Compute scalar dissipation (COMPARE WITH 2D CASE!!!)
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*uj+&
                      IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*vj+&
                      IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*wj)+&
             sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)+&
                  POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2)+&
                  POW(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),2))*cj,&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*ui+&
                      IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*vi+&
                      IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*wi)+&
             sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)+&
                  POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2)+&
                  POW(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),2))*ci )
#endif

      ! Multiply the solution difference by the scalar dissipation
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,_,_,_)-&
                   IDX3(DdataAtEdge,:,1,idx,_,_,_))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------
      
#ifdef HYDRO_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj+&
           IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fzj-&
           IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi-&
           IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fzi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fz_ij + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fz_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxRusDiss3d_sim

  !***************************************************************************

!<subroutine>

  subroutine hydro_calcDivVecRusDiss3d_cuda(rgroupFEMSet, rx, ry, dscale, bclear,&
      fcb_calcFluxSys_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 3D using scalar artificial viscosities of Rusanov-type.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: callback function to compute local fluxes
    include '../../../../../kernel/PDEOperators/intf_calcFluxSys_sim.inc'
    optional :: fcb_calcFluxSys_sim
!</input>

!<inputoutput>
    ! Destination vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT
    
    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: IEDGEmax,IEDGEset,igroup
    type(C_PTR) :: cptr_DcoeffsAtEdge
    type(C_PTR) :: cptr_IedgeList
    type(C_PTR) :: cptr_Dx, cptr_Dy
    integer(I64) :: istream

    ! Create CUDA stream
    call coproc_createStream(istream)
    
    ! Get pointers to memory blocks on device
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IedgeList, cptr_IedgeList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtEdge, cptr_DcoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    if (.not.storage_isAssociated(cptr_Dx)) then
      call lsysbl_copyH2D_Vector(rx, .false., .false., istream)
      call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    end if
   
    ! Make sure that the destination vector ry exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Vector(ry, .true., .false., istream)
    call storage_getMemPtrOnDevice(ry%h_Ddata, cptr_Dy)
   
    ! Set pointer
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
    
    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronise memory access
    do igroup = 1, size(p_IedgeListIdx)-1
      
      ! Do nothing for empty groups
      if (p_IedgeListIdx(igroup+1)-p_IedgeListIdx(igroup) .le. 0) cycle

      ! Get position of first edge in group
      IEDGEset = p_IedgeListIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IedgeListIdx(igroup+1)-1
      
      ! Use callback function to compute internodal fluxes
      call hydro_calcFluxRusDiss3d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
          cptr_Dx, cptr_Dy, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
          rgroupFEMSet%NEDGE, rgroupFEMSet%ncoeffsAtEdge, IEDGEmax-IEDGEset+1,&
          IEDGEset, istream)
    end do

    ! Transfer destination vector back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Vector(ry, bclear, .false., istream)

    ! Ensure data consistency
    call coproc_synchronizeStream(istream)
    call coproc_destroyStream(istream)

#else

    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivVecRusDiss3d_cuda')
    call sys_halt()
    
#endif

  end subroutine hydro_calcDivVecRusDiss3d_cuda

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRusDissDiSp3d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

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
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IedgeList

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
    real(DP), dimension(NVAR3D) :: Fxi,Fyi,Fxj,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP) :: Ei,Ej,ci,cj,pi,pj,ui,uj,vi,vj,wi,wj
    real(DP) :: d_ij
    integer :: idx

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)

      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,wi,pi)-&
                      INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,wj,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !-------------------------------------------------------------------------

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*&
          (Ei-RCONST(0.5)*(ui*ui+vi*vi+wi*wi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*&
          (Ej-RCONST(0.5)*(uj*uj+vj*vj+wj*wj)), SYS_EPSREAL_DP))
      
#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation with dimensional splitting based on
      ! the skew-symmetric part which does not include the symmetric
      ! boundary contribution
      d_ij = max( abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*uj)+&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)))*cj,&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_))*ui)+&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)))*ci )&
           + max( abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*vj)+&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)))*cj,&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_))*vi)+&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)))*ci )&
           + max( abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,3,2,idx,_,_,_))*wj)+&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)))*cj,&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,3,1,idx,_,_,_))*wi)+&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)))*ci )
#else
      ! Compute scalar dissipation with dimensional splitting
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*uj)+&
                  abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_))*cj,&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*ui)+&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*ci )&
           + max( abs(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*vj)+&
                  abs(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_))*cj,&
                  abs(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*vi)+&
                  abs(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*ci )&
           + max( abs(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*wj)+&
                  abs(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_))*cj,&
                  abs(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*wi)+&
                  abs(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_))*ci )
#endif

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,_,_,_)-&
                   IDX3(DdataAtEdge,:,1,idx,_,_,_))
      
      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj+&
           IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fzj-&
           IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi-&
           IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fzi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*Fz_ij + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fz_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxRusDissDiSp3d_sim

  !***************************************************************************

!<subroutine>

  subroutine hydro_calcDivVecRusDissDiSp3d_cuda(rgroupFEMSet, rx, ry, dscale,&
      bclear, fcb_calcFluxSys_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 3D using scalar artificial viscosities of Rusanov-type, whereby
    ! dimensional splitting is employed.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
    
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: callback function to compute local fluxes
    include '../../../../../kernel/PDEOperators/intf_calcFluxSys_sim.inc'
    optional :: fcb_calcFluxSys_sim
!</input>

!<inputoutput>
    ! Destination vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT
    
    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: IEDGEmax,IEDGEset,igroup
    type(C_PTR) :: cptr_DcoeffsAtEdge
    type(C_PTR) :: cptr_IedgeList
    type(C_PTR) :: cptr_Dx, cptr_Dy
    integer(I64) :: istream

    ! Create CUDA stream
    call coproc_createStream(istream)
    
    ! Get pointers to memory blocks on device
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IedgeList, cptr_IedgeList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtEdge, cptr_DcoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    if (.not.storage_isAssociated(cptr_Dx)) then
      call lsysbl_copyH2D_Vector(rx, .false., .false., istream)
      call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    end if
   
    ! Make sure that the destination vector ry exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Vector(ry, .true., .false., istream)
    call storage_getMemPtrOnDevice(ry%h_Ddata, cptr_Dy)
   
    ! Set pointer
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
    
    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronise memory access
    do igroup = 1, size(p_IedgeListIdx)-1
      
      ! Do nothing for empty groups
      if (p_IedgeListIdx(igroup+1)-p_IedgeListIdx(igroup) .le. 0) cycle

      ! Get position of first edge in group
      IEDGEset = p_IedgeListIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IedgeListIdx(igroup+1)-1
      
      ! Use callback function to compute internodal fluxes
      call hydro_calcFluxRusDissDiSp3d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
          cptr_Dx, cptr_Dy, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
          rgroupFEMSet%NEDGE, rgroupFEMSet%ncoeffsAtEdge, IEDGEmax-IEDGEset+1,&
          IEDGEset, istream)
    end do

    ! Transfer destination vector back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Vector(ry, bclear, .false., istream)

    ! Ensure data consistency
    call coproc_synchronizeStream(istream)
    call coproc_destroyStream(istream)

#else

    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivVecRusDissDiSp3d_cuda')
    call sys_halt()
    
#endif

  end subroutine hydro_calcDivVecRusDissDiSp3d_cuda

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiagMatD3d_sim(DdataAtNode, DcoeffsAtNode,&
      InodeList, dscale, nnodes, DmatrixAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 3D.
!</description>

!<input>
  ! Nodal solution values for all nodes under consideration
  real(DP), dimension(:,:), intent(in) :: DdataAtNode

  ! Entries of the coefficient matrices for all nodes under consideration
  real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode

  ! Numbers of vertices and matrix entries for all nodes under consideration
  integer, dimension(:,:), intent(in) :: InodeList

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
  real(DP), dimension(:,:,:), intent(out) :: DmatrixAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP) :: ui,vi,wi
    integer :: inode


    do inode = 1, nnodes
      
      ! Compute auxiliary variables
      ui = XVELOCITY2_3D(DdataAtNode,IDX2,inode,_,_)
      vi = YVELOCITY2_3D(DdataAtNode,IDX2,inode,_,_)
      wi = ZVELOCITY2_3D(DdataAtNode,IDX2,inode,_,_)
      
#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ii = diag(A_i)*C_{ii}$
      IDX3(DmatrixAtNode,1,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtNode,2,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtNode,3,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtNode,4,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtNode,5,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,_)
#else
      ! Compute Galerkin coefficient $K_ii = -diag(A_i)*C_{ii}$
      IDX3(DmatrixAtNode,1,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtNode,2,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtNode,3,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtNode,4,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtNode,5,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,_)
#endif
    end do

  end subroutine hydro_calcMatDiagMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiag3d_sim(DdataAtNode,&
      DcoeffsAtNode, InodeList, dscale, nnodes,&
      DmatrixAtNode, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 3D.
!</description>

!<input>
  ! Nodal solution values for all nodes under consideration
  real(DP), dimension(:,:), intent(in) :: DdataAtNode

  ! Entries of the coefficient matrices for all nodes under consideration
  real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode

  ! Numbers of vertices and matrix entries for all nodes under consideration
  integer, dimension(:,:), intent(in) :: InodeList

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
  real(DP), dimension(:,:,:), intent(out) :: DmatrixAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP) :: Ei,ui,vi,wi
    integer :: inode


    do inode = 1, nnodes
      
      ! Compute auxiliary variables
      ui = XVELOCITY2_3D(DdataAtNode,IDX2,inode,_,_)
      vi = YVELOCITY2_3D(DdataAtNode,IDX2,inode,_,_)
      wi = ZVELOCITY2_3D(DdataAtNode,IDX2,inode,_,_)
      Ei = SPECIFICTOTALENERGY2_3D(DdataAtNode,IDX2,inode,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      IDX3(DmatrixAtNode,1,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,2,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,3,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,4,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,5,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,6,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,7,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,8,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,9,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,10,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,11,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,12,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,13,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,14,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,15,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,16,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,17,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,18,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,19,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,20,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,21,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,22,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,23,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,24,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,25,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                    IDX2(DcoeffsAtNode,2,inode,_,_),
                                    IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)

#else
      ! Compute Galerkin coefficient $K_ii = -A_i*C_{ii}$
      IDX3(DmatrixAtNode,1,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,2,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,3,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,4,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,5,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,6,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,7,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,8,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,9,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,10,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,11,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,12,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,13,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,14,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,15,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,16,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,17,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,18,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,19,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,20,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,21,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,22,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,23,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,24,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtNode,25,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),
                                     IDX2(DcoeffsAtNode,2,inode,_,_),
                                     IDX2(DcoeffsAtNode,3,inode,_,_),ui,vi,wi,Ei)
#endif
    end do

  end subroutine hydro_calcMatDiag3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGalMatD3d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 3D.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IedgeList

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
  real(DP), dimension(:,:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DmatrixAtEdge,1,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,2,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,3,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,4,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,5,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DmatrixAtEdge,1,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,2,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,3,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,4,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,5,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
#endif
    end do

  end subroutine hydro_calcMatGalMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGalerkin3d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IedgeList

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
  real(DP), dimension(:,:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP) :: Ei,Ej,ui,uj,vi,vj,wi,wj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DmatrixAtEdge,1,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,2,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,3,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,4,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,5,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,6,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,7,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,8,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,9,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,10,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,11,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,12,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,13,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,14,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,15,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,16,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,17,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,18,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,19,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,20,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,21,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,22,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,23,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,24,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,25,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,17,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,18,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,19,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,20,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,21,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,22,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,23,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,24,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,25,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DmatrixAtEdge,1,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,2,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,3,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,4,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,5,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,6,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,7,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,8,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,9,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,10,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,11,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,12,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,13,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,14,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,15,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,16,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,17,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,18,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,19,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,20,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,21,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,22,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,23,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,24,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,25,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,17,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,18,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,19,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,20,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,21,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,22,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,23,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,24,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,25,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
#endif
    end do

  end subroutine hydro_calcMatGalerkin3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDissMatD3d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 3D and applies scalar artificial viscosities proportional to
    ! the spectral radius (largest eigenvalue) of the Roe-matrix.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IedgeList

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
  real(DP), dimension(:,:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,anorm,aux,c_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    integer :: idx

    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      
      if (anorm .gt. SYS_EPSREAL_DP) then

        ! Compute densities
        ri = DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute pressures
        pi = PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        pj = PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute enthalpies
        hi = (TOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2) + w_ij*a(3)
        q_ij   = RCONST(0.5)*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

        ! Compute the speed of sound
        c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))

        ! Compute scalar dissipation
        IDX3(DmatrixAtEdge,:,1,idx,_,_,_) = dscale * (abs(vel_ij) + anorm*c_ij)

      else
        
        ! Nullify dissipation tensor
        IDX3(DmatrixAtEdge,:,1,idx,_,_,_) = RCONST(0.0)

      end if
    end do

  end subroutine hydro_calcMatScDissMatD3d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDiss3d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D and applies
    ! scalar artificial viscosities proportional to the spectral
    ! radius (largest eigenvalue) of the Roe-matrix.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IedgeList

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
  real(DP), dimension(:,:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: Ei,Ej,hi,hj,pi,pj,ri,rj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,anorm,aux,c_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,17,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,18,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,19,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,20,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,21,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,22,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,23,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,24,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,25,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,6,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,7,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,9,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,10,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,11,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,12,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,13,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,14,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,15,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,16,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,17,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,18,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,19,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,20,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,21,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,22,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,23,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,24,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,25,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,17,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,18,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,19,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,20,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,21,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,22,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,23,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,24,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,25,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,6,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,7,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,9,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,10,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,11,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,12,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,13,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,14,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,15,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,16,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,17,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,18,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,19,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,20,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,21,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,22,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,23,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,24,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,25,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      
      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Compute densities
        ri = DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute pressures
        pi = PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        pj = PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute enthalpies
        hi = (TOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2) + w_ij*a(3)
        q_ij   = RCONST(0.5)*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

        ! Compute the speed of sound
        c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))

        ! Compute scalar dissipation
        aux = dscale * (abs(vel_ij) + anorm*c_ij)

        IDX3(DmatrixAtEdge,:,1,idx,_,_,_) = RCONST(0.0)
        IDX3(DmatrixAtEdge, 1,1,idx,_,_,_) = aux
        IDX3(DmatrixAtEdge, 6,1,idx,_,_,_) = aux
        IDX3(DmatrixAtEdge,11,1,idx,_,_,_) = aux
        IDX3(DmatrixAtEdge,16,1,idx,_,_,_) = aux

      else

        ! Nullify dissipation tensor
        IDX3(DmatrixAtEdge,:,1,idx,_,_,_) = RCONST(0.0)
        
      end if
    end do

  end subroutine hydro_calcMatScDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDissMatD3d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D and applies
    ! tensorial artificial viscosities of Roe-type.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IedgeList

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
  real(DP), dimension(:,:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NVAR3D,NVAR3D) :: R_ij,L_ij
    real(DP), dimension(NVAR3D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,anorm,aux,cPow2_ij,c_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    real(DP) :: l1,l2,l3,l4,l5
    integer :: idx
#if defined(HYDRO_USE_ENTROPYFIX) && (HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX)
    real(DP) :: dtol,ci,cj,veli,velj
#endif


    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
#endif
        
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      if (anorm .gt. SYS_EPSREAL_DP) then

        ! Normalise the skew-symmetric coefficient
        a = a/anorm

        ! Compute densities
        ri = DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute pressures
        pi = PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        pj = PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute enthalpies
        hi = (TOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)

        ! Compute auxiliary values
        vel_ij = u_ij*a(1)+v_ij*a(2)+w_ij*a(3)
        q_ij   = RCONST(0.5)*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)
        
        ! Compute speed of sound
        cPow2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij     = sqrt(cPow2_ij)

        ! Diagonal matrix of eigenvectors
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)
        l5 = abs(vel_ij)

#if defined(HYDRO_USE_ENTROPYFIX)

#if HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX

        ! Entropy-fix by Harten and Hyman
        ci = SOUNDSPEED3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        cj = SOUNDSPEED3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
        veli = ui*a(1) + vi*a(2) + wi*a(3)
        velj = uj*a(1) + vj*a(2) + wj*a(3)

        dtol = max(RCONST(0.0), (vel_ij-c_ij) - (veli-ci), (velj-cj) - (vel_ij-c_ij) )
        
        if (l1 .lt. dtol)&
            l1 = RCONST(0.5)*((l1*l1)/dtol + dtol)
        
        dtol = max(RCONST(0.0), (vel_ij+c_ij) - (veli+ci), (velj+cj) - (vel_ij+c_ij) )

        if (l3 .lt. dtol)&
            l3 = RCONST(0.5)*((l3*l3)/dtol + dtol)

#elif HYDRO_USE_ENTROPYFIX == HARTEN_ENTROPYFIX

#ifndef HYDRO_HARTEN_ENTROPYFIX
#error "Value HYDRO_HARTEN_ENTROPYFIX is required!"
#else
        ! Entropy-fix by Harten
        if (l1 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l1 = RCONST(0.5)*((l1*l1)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))

        if (l3 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l3 = RCONST(0.5)*((l3*l3)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))
#endif
#else
#error "Invalid type of entropy fix!"
#endif
#endif

        ! Get the dimension with largest coefficient
        select case(maxloc(a,1))
        case(1)
          ! Matrix of right eigenvectors
          R_ij(1,1) =  l1
          R_ij(2,1) =  l1*(u_ij-c_ij*a(1))
          R_ij(3,1) =  l1*(v_ij-c_ij*a(2))
          R_ij(4,1) =  l1*(w_ij-c_ij*a(3))
          R_ij(5,1) =  l1*(H_ij-c_ij*vel_ij)
          
          R_ij(1,2) =  l2
          R_ij(2,2) =  l2*u_ij
          R_ij(3,2) =  l2*v_ij
          R_ij(4,2) =  l2*w_ij
          R_ij(5,2) =  l2*q_ij
        
          R_ij(1,3) =  l3
          R_ij(2,3) =  l3*(u_ij+c_ij*a(1))
          R_ij(3,3) =  l3*(v_ij+c_ij*a(2))
          R_ij(4,3) =  l3*(v_ij+c_ij*a(3))
          R_ij(5,3) =  l3*(H_ij+c_ij*vel_ij)
        
          R_ij(1,4) =  RCONST(0.0)
          R_ij(2,4) =  l4*a(2)
          R_ij(3,4) = -l4*a(1)
          R_ij(4,4) =  RCONST(0.0)
          R_ij(5,4) =  l4*(u_ij*a(2)-v_ij*a(1))

          R_ij(1,5) =  RCONST(0.0)
          R_ij(2,5) = -l5*a(3)
          R_ij(3,5) =  RCONST(0.0)
          R_ij(4,5) =  l5*a(1)
          R_ij(5,5) =  l5*(w_ij*a(1)-u_ij*a(3))

          ! Matrix of left eigenvectors
          L_ij(1,1) =  RCONST(0.5)*(((HYDRO_GAMMA)-RCONST(1.0))*q_ij+c_ij*vel_ij)/cPow2_ij
          L_ij(2,1) =  (cPow2_ij-((HYDRO_GAMMA)-RCONST(1.0))*q_ij)/cPow2_ij
          L_ij(3,1) =  RCONST(0.5)*(((HYDRO_GAMMA)-RCONST(1.0))*q_ij-c_ij*vel_ij)/cPow2_ij
          L_ij(4,1) =  (v_ij-vel_ij*a(2))/a(1)
          L_ij(5,1) =  (vel_ij*a(3)-w_ij)/a(1)
          
          L_ij(1,2) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*u_ij-c_ij*a(1))/cPow2_ij
          L_ij(2,2) =  ((HYDRO_GAMMA)-RCONST(1.0))*u_ij/cPow2_ij
          L_ij(3,2) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*u_ij+c_ij*a(1))/cPow2_ij
          L_ij(4,2) =  a(2)
          L_ij(5,2) = -a(3)
          
          L_ij(1,3) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*v_ij-c_ij*a(2))/cPow2_ij
          L_ij(2,3) =  ((HYDRO_GAMMA)-RCONST(1.0))*v_ij/cPow2_ij
          L_ij(3,3) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*v_ij+c_ij*a(2))/cPow2_ij
          L_ij(4,3) =  (a(2)*a(2)-RCONST(1.0))/a(1)
          L_ij(5,3) = -a(2)*a(3)/a(1)
          
          L_ij(1,4) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*w_ij-c_ij*a(3))/cPow2_ij
          L_ij(2,4) =  ((HYDRO_GAMMA)-RCONST(1.0))*w_ij/cPow2_ij
          L_ij(3,4) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*w_ij+c_ij*a(3))/cPow2_ij
          L_ij(4,4) =  a(2)*a(3)/a(1)
          L_ij(5,4) = -(1-a(3)*a(3))/a(1)

          L_ij(1,5) =  RCONST(0.5)*((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(2,5) = -((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(3,5) =  RCONST(0.5)*((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(4,5) =  RCONST(0.0)
          L_ij(5,5) =  RCONST(0.0)

        case(2)
          ! Matrix of right eigenvectors
          R_ij(1,1) =  l1
          R_ij(2,1) =  l1*(u_ij-c_ij*a(1))
          R_ij(3,1) =  l1*(v_ij-c_ij*a(2))
          R_ij(4,1) =  l1*(w_ij-c_ij*a(3))
          R_ij(5,1) =  l1*(H_ij-c_ij*vel_ij)
          
          R_ij(1,2) =  l2
          R_ij(2,2) =  l2*u_ij
          R_ij(3,2) =  l2*v_ij
          R_ij(4,2) =  l2*w_ij
          R_ij(5,2) =  l2*q_ij
        
          R_ij(1,3) =  l3
          R_ij(2,3) =  l3*(u_ij+c_ij*a(1))
          R_ij(3,3) =  l3*(v_ij+c_ij*a(2))
          R_ij(4,3) =  l3*(v_ij+c_ij*a(3))
          R_ij(5,3) =  l3*(H_ij+c_ij*vel_ij)
        
          R_ij(1,4) =  RCONST(0.0)
          R_ij(2,4) =  l4*a(2)
          R_ij(3,4) = -l4*a(1)
          R_ij(4,4) =  RCONST(0.0)
          R_ij(5,4) =  l4*(u_ij*a(2)-v_ij*a(1))

          R_ij(1,5) =  RCONST(0.0)
          R_ij(2,5) =  RCONST(0.0)
          R_ij(3,5) =  l5*a(3)
          R_ij(4,5) = -l5*a(2)
          R_ij(5,5) =  l5*(v_ij*a(3)-w_ij*a(2))

          ! Matrix of left eigenvectors
          L_ij(1,1) =  RCONST(0.5)*(((HYDRO_GAMMA)-RCONST(1.0))*q_ij+c_ij*vel_ij)/cPow2_ij
          L_ij(2,1) =  (cPow2_ij-((HYDRO_GAMMA)-RCONST(1.0))*q_ij)/cPow2_ij
          L_ij(3,1) =  RCONST(0.5)*(((HYDRO_GAMMA)-RCONST(1.0))*q_ij-c_ij*vel_ij)/cPow2_ij
          L_ij(4,1) =  (vel_ij*a(1)-u_ij)/a(2)
          L_ij(5,1) =  (w_ij-vel_ij*a(3))/a(2)
          
          L_ij(1,2) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*u_ij-c_ij*a(1))/cPow2_ij
          L_ij(2,2) =  ((HYDRO_GAMMA)-RCONST(1.0))*u_ij/cPow2_ij
          L_ij(3,2) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*u_ij+c_ij*a(1))/cPow2_ij
          L_ij(4,2) =  (RCONST(1.0)-a(1)*a(1))/a(2)
          L_ij(5,2) =  a(1)*a(3)/a(2)
          
          L_ij(1,3) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*v_ij-c_ij*a(2))/cPow2_ij
          L_ij(2,3) =  ((HYDRO_GAMMA)-RCONST(1.0))*v_ij/cPow2_ij
          L_ij(3,3) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*v_ij+c_ij*a(2))/cPow2_ij
          L_ij(4,3) = -a(1)
          L_ij(5,3) =  a(3)
          
          L_ij(1,4) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*w_ij-c_ij*a(3))/cPow2_ij
          L_ij(2,4) =  ((HYDRO_GAMMA)-RCONST(1.0))*w_ij/cPow2_ij
          L_ij(3,4) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*w_ij+c_ij*a(3))/cPow2_ij
          L_ij(4,4) = -a(1)*a(3)/a(2)
          L_ij(5,4) =  (a(3)*a(3)-RCONST(1.0))/a(2)

          L_ij(1,5) =  RCONST(0.5)*((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(2,5) = -((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(3,5) =  RCONST(0.5)*((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(4,5) =  RCONST(0.0)
          L_ij(5,5) =  RCONST(0.0)
          
        case(3)
          ! Matrix of right eigenvectors
          R_ij(1,1) =  l1
          R_ij(2,1) =  l1*(u_ij-c_ij*a(1))
          R_ij(3,1) =  l1*(v_ij-c_ij*a(2))
          R_ij(4,1) =  l1*(w_ij-c_ij*a(3))
          R_ij(5,1) =  l1*(H_ij-c_ij*vel_ij)
          
          R_ij(1,2) =  l2
          R_ij(2,2) =  l2*u_ij
          R_ij(3,2) =  l2*v_ij
          R_ij(4,2) =  l2*w_ij
          R_ij(5,2) =  l2*q_ij
        
          R_ij(1,3) =  l3
          R_ij(2,3) =  l3*(u_ij+c_ij*a(1))
          R_ij(3,3) =  l3*(v_ij+c_ij*a(2))
          R_ij(4,3) =  l3*(v_ij+c_ij*a(3))
          R_ij(5,3) =  l3*(H_ij+c_ij*vel_ij)
        
          R_ij(1,4) =  RCONST(0.0)
          R_ij(2,4) = -l4*a(3)
          R_ij(3,4) =  RCONST(0.0)
          R_ij(4,4) =  l4*a(1)
          R_ij(5,4) =  l4*(w_ij*a(1)-u_ij*a(3))

          R_ij(1,5) =  RCONST(0.0)
          R_ij(2,5) =  RCONST(0.0)
          R_ij(3,5) =  l5*a(3)
          R_ij(4,5) = -l5*a(2)
          R_ij(5,5) =  l5*(v_ij*a(3)-w_ij*a(2))

          ! Matrix of left eigenvectors
          L_ij(1,1) =  RCONST(0.5)*(((HYDRO_GAMMA)-RCONST(1.0))*q_ij+c_ij*vel_ij)/cPow2_ij
          L_ij(2,1) =  (cPow2_ij-((HYDRO_GAMMA)-RCONST(1.0))*q_ij)/cPow2_ij
          L_ij(3,1) =  RCONST(0.5)*(((HYDRO_GAMMA)-RCONST(1.0))*q_ij-c_ij*vel_ij)/cPow2_ij
          L_ij(4,1) =  (u_ij-vel_ij*a(1))/a(3)
          L_ij(5,1) =  (vel_ij*a(2)-v_ij)/a(3)
          
          L_ij(1,2) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*u_ij-c_ij*a(1))/cPow2_ij
          L_ij(2,2) =  ((HYDRO_GAMMA)-RCONST(1.0))*u_ij/cPow2_ij
          L_ij(3,2) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*u_ij+c_ij*a(1))/cPow2_ij
          L_ij(4,2) =  (a(1)*a(1)-RCONST(1.0))/a(3)
          L_ij(5,2) = -a(1)*a(2)/a(3)
          
          L_ij(1,3) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*v_ij-c_ij*a(2))/cPow2_ij
          L_ij(2,3) =  ((HYDRO_GAMMA)-RCONST(1.0))*v_ij/cPow2_ij
          L_ij(3,3) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*v_ij+c_ij*a(2))/cPow2_ij
          L_ij(4,3) =  a(1)*a(2)/a(3)
          L_ij(5,3) =  (RCONST(1.0)-a(2)*a(2))/a(3)
          
          L_ij(1,4) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*w_ij-c_ij*a(3))/cPow2_ij
          L_ij(2,4) =  ((HYDRO_GAMMA)-RCONST(1.0))*w_ij/cPow2_ij
          L_ij(3,4) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*w_ij+c_ij*a(3))/cPow2_ij
          L_ij(4,4) =  a(1)
          L_ij(5,4) = -a(2)

          L_ij(1,5) =  RCONST(0.5)*((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(2,5) = -((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(3,5) =  RCONST(0.5)*((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(4,5) =  RCONST(0.0)
          L_ij(5,5) =  RCONST(0.0)

        end select

        ! Include scaling parameter
        anorm = dscale*anorm
        
        ! Compute tensorial dissipation D_ij = diag(R_ij*|Lbd_ij|*L_ij)*I
        IDX3(DmatrixAtEdge,:,1,idx,_,_,_) = RCONST(0.0)
        IDX3(DmatrixAtEdge,1,1,idx,_,_,_) = anorm*( R_ij(1,1)*L_ij(1,1)+&
                                                          R_ij(1,2)*L_ij(2,1)+&
                                                          R_ij(1,3)*L_ij(3,1)+&
                                                          R_ij(1,4)*L_ij(4,1)+&
                                                          R_ij(1,5)*L_ij(5,1)  )
        IDX3(DmatrixAtEdge,2,1,idx,_,_,_) = anorm*( R_ij(2,1)*L_ij(1,2)+&
                                                          R_ij(2,2)*L_ij(2,2)+&
                                                          R_ij(2,3)*L_ij(3,2)+&
                                                          R_ij(2,4)*L_ij(4,2)+&
                                                          R_ij(2,5)*L_ij(5,2)  )
        IDX3(DmatrixAtEdge,3,1,idx,_,_,_) = anorm*( R_ij(3,1)*L_ij(1,3)+&
                                                          R_ij(3,2)*L_ij(2,3)+&
                                                          R_ij(3,3)*L_ij(3,3)+&
                                                          R_ij(3,4)*L_ij(4,3)+&
                                                          R_ij(3,5)*L_ij(5,3)  )
        IDX3(DmatrixAtEdge,4,1,idx,_,_,_) = anorm*( R_ij(4,1)*L_ij(1,4)+&
                                                          R_ij(4,2)*L_ij(2,4)+&
                                                          R_ij(4,3)*L_ij(3,4)+&
                                                          R_ij(4,4)*L_ij(4,4)+&
                                                          R_ij(4,5)*L_ij(5,4)  )
        IDX3(DmatrixAtEdge,5,1,idx,_,_,_) = anorm*( R_ij(5,1)*L_ij(1,5)+&
                                                          R_ij(5,2)*L_ij(2,5)+&
                                                          R_ij(5,3)*L_ij(3,5)+&
                                                          R_ij(5,4)*L_ij(4,5)+&
                                                          R_ij(5,5)*L_ij(5,5)  )
      else

        ! Nullify dissipation tensor
        IDX3(DmatrixAtEdge,:,1,idx,_,_,_) = RCONST(0.0)
        
      end if
    end do

  end subroutine hydro_calcMatRoeDissMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDiss3d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D and applies
    ! tensorial artificial viscosities of Roe-type.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IedgeList

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
  real(DP), dimension(:,:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NVAR3D,NVAR3D) :: R_ij,L_ij
    real(DP), dimension(NVAR3D) :: a
    real(DP) :: Ei,Ej,hi,hj,pi,pj,ri,rj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,anorm,aux,cPow2_ij,c_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    real(DP) :: l1,l2,l3,l4,l5
    integer :: idx,i,j,k
#if defined(HYDRO_USE_ENTROPYFIX) && (HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX)
    real(DP) :: dtol,ci,cj,veli,velj
#endif


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,17,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,18,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,19,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,20,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,21,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,22,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,23,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,24,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,25,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,6,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,7,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,9,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,10,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,11,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,12,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,13,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,14,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,15,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,16,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,17,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,18,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,19,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,20,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,21,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,22,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,23,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,24,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,25,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,17,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,18,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,19,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,20,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,21,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,22,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,23,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,24,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,25,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,6,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,7,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,9,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,10,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,11,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,12,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,13,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,14,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,15,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,16,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,17,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,18,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,19,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,20,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,21,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,22,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,23,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,24,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,25,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      
      if (anorm .gt. SYS_EPSREAL_DP) then

        ! Compute densities
        ri = DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute pressures
        pi = PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        pj = PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute enthalpies
        hi = (TOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)

        ! Compute auxiliary values
        vel_ij = u_ij*a(1)+v_ij*a(2)+w_ij*a(3)
        q_ij   = RCONST(0.5)*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)
        
        ! Compute speed of sound
        cPow2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij     = sqrt(cPow2_ij)

        ! Diagonal matrix of eigenvectors
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)
        l5 = abs(vel_ij)

#if defined(HYDRO_USE_ENTROPYFIX)

#if HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX

        ! Entropy-fix by Harten and Hyman
        ci = SOUNDSPEED3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        cj = SOUNDSPEED3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
        veli = ui*a(1) + vi*a(2) + wi*a(3)
        velj = uj*a(1) + vj*a(2) + wj*a(3)

        dtol = max(RCONST(0.0), (vel_ij-c_ij) - (veli-ci), (velj-cj) - (vel_ij-c_ij) )
        
        if (l1 .lt. dtol)&
            l1 = RCONST(0.5)*((l1*l1)/dtol + dtol)
        
        dtol = max(RCONST(0.0), (vel_ij+c_ij) - (veli+ci), (velj+cj) - (vel_ij+c_ij) )

        if (l3 .lt. dtol)&
            l3 = RCONST(0.5)*((l3*l3)/dtol + dtol)

#elif HYDRO_USE_ENTROPYFIX == HARTEN_ENTROPYFIX

#ifndef HYDRO_HARTEN_ENTROPYFIX
#error "Value HYDRO_HARTEN_ENTROPYFIX is required!"
#else
        ! Entropy-fix by Harten
        if (l1 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l1 = RCONST(0.5)*((l1*l1)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))

        if (l3 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l3 = RCONST(0.5)*((l3*l3)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))
#endif
#else
#error "Invalid type of entropy fix!"
#endif
#endif
        
        ! Get the dimension with largest coefficient
        select case(maxloc(a,1))
        case(1)
          ! Matrix of right eigenvectors
          R_ij(1,1) =  l1
          R_ij(2,1) =  l1*(u_ij-c_ij*a(1))
          R_ij(3,1) =  l1*(v_ij-c_ij*a(2))
          R_ij(4,1) =  l1*(w_ij-c_ij*a(3))
          R_ij(5,1) =  l1*(H_ij-c_ij*vel_ij)
          
          R_ij(1,2) =  l2
          R_ij(2,2) =  l2*u_ij
          R_ij(3,2) =  l2*v_ij
          R_ij(4,2) =  l2*w_ij
          R_ij(5,2) =  l2*q_ij
        
          R_ij(1,3) =  l3
          R_ij(2,3) =  l3*(u_ij+c_ij*a(1))
          R_ij(3,3) =  l3*(v_ij+c_ij*a(2))
          R_ij(4,3) =  l3*(v_ij+c_ij*a(3))
          R_ij(5,3) =  l3*(H_ij+c_ij*vel_ij)
        
          R_ij(1,4) =  RCONST(0.0)
          R_ij(2,4) =  l4*a(2)
          R_ij(3,4) = -l4*a(1)
          R_ij(4,4) =  RCONST(0.0)
          R_ij(5,4) =  l4*(u_ij*a(2)-v_ij*a(1))

          R_ij(1,5) =  RCONST(0.0)
          R_ij(2,5) = -l5*a(3)
          R_ij(3,5) =  RCONST(0.0)
          R_ij(4,5) =  l5*a(1)
          R_ij(5,5) =  l5*(w_ij*a(1)-u_ij*a(3))

          ! Matrix of left eigenvectors
          L_ij(1,1) =  RCONST(0.5)*(((HYDRO_GAMMA)-RCONST(1.0))*q_ij+c_ij*vel_ij)/cPow2_ij
          L_ij(2,1) =  (cPow2_ij-((HYDRO_GAMMA)-RCONST(1.0))*q_ij)/cPow2_ij
          L_ij(3,1) =  RCONST(0.5)*(((HYDRO_GAMMA)-RCONST(1.0))*q_ij-c_ij*vel_ij)/cPow2_ij
          L_ij(4,1) =  (v_ij-vel_ij*a(2))/a(1)
          L_ij(5,1) =  (vel_ij*a(3)-w_ij)/a(1)
          
          L_ij(1,2) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*u_ij-c_ij*a(1))/cPow2_ij
          L_ij(2,2) =  ((HYDRO_GAMMA)-RCONST(1.0))*u_ij/cPow2_ij
          L_ij(3,2) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*u_ij+c_ij*a(1))/cPow2_ij
          L_ij(4,2) =  a(2)
          L_ij(5,2) = -a(3)
          
          L_ij(1,3) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*v_ij-c_ij*a(2))/cPow2_ij
          L_ij(2,3) =  ((HYDRO_GAMMA)-RCONST(1.0))*v_ij/cPow2_ij
          L_ij(3,3) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*v_ij+c_ij*a(2))/cPow2_ij
          L_ij(4,3) =  (a(2)*a(2)-RCONST(1.0))/a(1)
          L_ij(5,3) = -a(2)*a(3)/a(1)
          
          L_ij(1,4) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*w_ij-c_ij*a(3))/cPow2_ij
          L_ij(2,4) =  ((HYDRO_GAMMA)-RCONST(1.0))*w_ij/cPow2_ij
          L_ij(3,4) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*w_ij+c_ij*a(3))/cPow2_ij
          L_ij(4,4) =  a(2)*a(3)/a(1)
          L_ij(5,4) = -(1-a(3)*a(3))/a(1)

          L_ij(1,5) =  RCONST(0.5)*((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(2,5) = -((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(3,5) =  RCONST(0.5)*((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(4,5) =  RCONST(0.0)
          L_ij(5,5) =  RCONST(0.0)

        case(2)
          ! Matrix of right eigenvectors
          R_ij(1,1) =  l1
          R_ij(2,1) =  l1*(u_ij-c_ij*a(1))
          R_ij(3,1) =  l1*(v_ij-c_ij*a(2))
          R_ij(4,1) =  l1*(w_ij-c_ij*a(3))
          R_ij(5,1) =  l1*(H_ij-c_ij*vel_ij)
          
          R_ij(1,2) =  l2
          R_ij(2,2) =  l2*u_ij
          R_ij(3,2) =  l2*v_ij
          R_ij(4,2) =  l2*w_ij
          R_ij(5,2) =  l2*q_ij
        
          R_ij(1,3) =  l3
          R_ij(2,3) =  l3*(u_ij+c_ij*a(1))
          R_ij(3,3) =  l3*(v_ij+c_ij*a(2))
          R_ij(4,3) =  l3*(v_ij+c_ij*a(3))
          R_ij(5,3) =  l3*(H_ij+c_ij*vel_ij)
        
          R_ij(1,4) =  RCONST(0.0)
          R_ij(2,4) =  l4*a(2)
          R_ij(3,4) = -l4*a(1)
          R_ij(4,4) =  RCONST(0.0)
          R_ij(5,4) =  l4*(u_ij*a(2)-v_ij*a(1))

          R_ij(1,5) =  RCONST(0.0)
          R_ij(2,5) =  RCONST(0.0)
          R_ij(3,5) =  l5*a(3)
          R_ij(4,5) = -l5*a(2)
          R_ij(5,5) =  l5*(v_ij*a(3)-w_ij*a(2))

          ! Matrix of left eigenvectors
          L_ij(1,1) =  RCONST(0.5)*(((HYDRO_GAMMA)-RCONST(1.0))*q_ij+c_ij*vel_ij)/cPow2_ij
          L_ij(2,1) =  (cPow2_ij-((HYDRO_GAMMA)-RCONST(1.0))*q_ij)/cPow2_ij
          L_ij(3,1) =  RCONST(0.5)*(((HYDRO_GAMMA)-RCONST(1.0))*q_ij-c_ij*vel_ij)/cPow2_ij
          L_ij(4,1) =  (vel_ij*a(1)-u_ij)/a(2)
          L_ij(5,1) =  (w_ij-vel_ij*a(3))/a(2)
          
          L_ij(1,2) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*u_ij-c_ij*a(1))/cPow2_ij
          L_ij(2,2) =  ((HYDRO_GAMMA)-RCONST(1.0))*u_ij/cPow2_ij
          L_ij(3,2) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*u_ij+c_ij*a(1))/cPow2_ij
          L_ij(4,2) =  (RCONST(1.0)-a(1)*a(1))/a(2)
          L_ij(5,2) =  a(1)*a(3)/a(2)
          
          L_ij(1,3) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*v_ij-c_ij*a(2))/cPow2_ij
          L_ij(2,3) =  ((HYDRO_GAMMA)-RCONST(1.0))*v_ij/cPow2_ij
          L_ij(3,3) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*v_ij+c_ij*a(2))/cPow2_ij
          L_ij(4,3) = -a(1)
          L_ij(5,3) =  a(3)
          
          L_ij(1,4) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*w_ij-c_ij*a(3))/cPow2_ij
          L_ij(2,4) =  ((HYDRO_GAMMA)-RCONST(1.0))*w_ij/cPow2_ij
          L_ij(3,4) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*w_ij+c_ij*a(3))/cPow2_ij
          L_ij(4,4) = -a(1)*a(3)/a(2)
          L_ij(5,4) =  (a(3)*a(3)-RCONST(1.0))/a(2)

          L_ij(1,5) =  RCONST(0.5)*((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(2,5) = -((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(3,5) =  RCONST(0.5)*((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(4,5) =  RCONST(0.0)
          L_ij(5,5) =  RCONST(0.0)
          
        case(3)
          ! Matrix of right eigenvectors
          R_ij(1,1) =  l1
          R_ij(2,1) =  l1*(u_ij-c_ij*a(1))
          R_ij(3,1) =  l1*(v_ij-c_ij*a(2))
          R_ij(4,1) =  l1*(w_ij-c_ij*a(3))
          R_ij(5,1) =  l1*(H_ij-c_ij*vel_ij)
          
          R_ij(1,2) =  l2
          R_ij(2,2) =  l2*u_ij
          R_ij(3,2) =  l2*v_ij
          R_ij(4,2) =  l2*w_ij
          R_ij(5,2) =  l2*q_ij
        
          R_ij(1,3) =  l3
          R_ij(2,3) =  l3*(u_ij+c_ij*a(1))
          R_ij(3,3) =  l3*(v_ij+c_ij*a(2))
          R_ij(4,3) =  l3*(v_ij+c_ij*a(3))
          R_ij(5,3) =  l3*(H_ij+c_ij*vel_ij)
        
          R_ij(1,4) =  RCONST(0.0)
          R_ij(2,4) = -l4*a(3)
          R_ij(3,4) =  RCONST(0.0)
          R_ij(4,4) =  l4*a(1)
          R_ij(5,4) =  l4*(w_ij*a(1)-u_ij*a(3))

          R_ij(1,5) =  RCONST(0.0)
          R_ij(2,5) =  RCONST(0.0)
          R_ij(3,5) =  l5*a(3)
          R_ij(4,5) = -l5*a(2)
          R_ij(5,5) =  l5*(v_ij*a(3)-w_ij*a(2))

          ! Matrix of left eigenvectors
          L_ij(1,1) =  RCONST(0.5)*(((HYDRO_GAMMA)-RCONST(1.0))*q_ij+c_ij*vel_ij)/cPow2_ij
          L_ij(2,1) =  (cPow2_ij-((HYDRO_GAMMA)-RCONST(1.0))*q_ij)/cPow2_ij
          L_ij(3,1) =  RCONST(0.5)*(((HYDRO_GAMMA)-RCONST(1.0))*q_ij-c_ij*vel_ij)/cPow2_ij
          L_ij(4,1) =  (u_ij-vel_ij*a(1))/a(3)
          L_ij(5,1) =  (vel_ij*a(2)-v_ij)/a(3)
          
          L_ij(1,2) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*u_ij-c_ij*a(1))/cPow2_ij
          L_ij(2,2) =  ((HYDRO_GAMMA)-RCONST(1.0))*u_ij/cPow2_ij
          L_ij(3,2) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*u_ij+c_ij*a(1))/cPow2_ij
          L_ij(4,2) =  (a(1)*a(1)-RCONST(1.0))/a(3)
          L_ij(5,2) = -a(1)*a(2)/a(3)
          
          L_ij(1,3) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*v_ij-c_ij*a(2))/cPow2_ij
          L_ij(2,3) =  ((HYDRO_GAMMA)-RCONST(1.0))*v_ij/cPow2_ij
          L_ij(3,3) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*v_ij+c_ij*a(2))/cPow2_ij
          L_ij(4,3) =  a(1)*a(2)/a(3)
          L_ij(5,3) =  (RCONST(1.0)-a(2)*a(2))/a(3)
          
          L_ij(1,4) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*w_ij-c_ij*a(3))/cPow2_ij
          L_ij(2,4) =  ((HYDRO_GAMMA)-RCONST(1.0))*w_ij/cPow2_ij
          L_ij(3,4) =  RCONST(0.5)*(-((HYDRO_GAMMA)-RCONST(1.0))*w_ij+c_ij*a(3))/cPow2_ij
          L_ij(4,4) =  a(1)
          L_ij(5,4) = -a(2)

          L_ij(1,5) =  RCONST(0.5)*((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(2,5) = -((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(3,5) =  RCONST(0.5)*((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
          L_ij(4,5) =  RCONST(0.0)
          L_ij(5,5) =  RCONST(0.0)

        end select

        ! Include scaling parameter
        anorm = dscale*anorm

        ! Compute tensorial dissipation D_ij = R_ij*|Lbd_ij|*L_ij
        do i = 1, NVAR3D
          do j = 1, NVAR3D
            aux = RCONST(0.0)
            do k = 1, NVAR3D
              aux = aux + R_ij(i,k)*L_ij(k,j)
            end do
            IDX3(DmatrixAtEdge,NVAR3D*(j-1)+i,1,idx,_,_,_) = anorm*aux
          end do
        end do
        
      else
        
        ! Nullify dissipation tensor
        IDX3(DmatrixAtEdge,:,1,idx,_,_,_) = RCONST(0.0)
        
      end if
    end do

  end subroutine hydro_calcMatRoeDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDissMatD3d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 3D and applies the scalar artificial viscosities of
    ! Rusanov-type.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IedgeList

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
  real(DP), dimension(:,:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP) :: Ei,Ej,ci,cj,ui,uj,vi,vj,wi,wj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,_)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate scalar artificial dissipation of Rusanov-type
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*&
          (Ei-RCONST(0.5)*(ui*ui+vi*vi+wi*wi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*&
          (Ej-RCONST(0.5)*(uj*uj+vj*vj+wj*wj)), SYS_EPSREAL_DP))
      
      ! Compute dissipation tensor
      IDX3(DmatrixAtEdge,:,1,idx,_,_,_) = dscale *&
          max( abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*uj+&
                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*vj+&
                   IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*wj) +&
          sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)+&
               POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2)+&
               POW(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),2))*cj,&
               abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*ui+&
                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*vi+&
                   IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*wi) +&
          sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)+&
               POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2)+&
               POW(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),2))*ci )
    end do

  end subroutine hydro_calcMatRusDissMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDiss3d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D applies
    ! scalar artificial viscosities of Rusanov-type.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IedgeList

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
  real(DP), dimension(:,:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP) :: Ei,Ej,aux,ci,cj,ui,uj,vi,vj,wi,wj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,17,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,18,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,19,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,20,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,21,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,22,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,23,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,24,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,25,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,6,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,7,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,9,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,10,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,11,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,12,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,13,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,14,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,15,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,16,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,17,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,18,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,19,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,20,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,21,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,22,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,23,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,24,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,25,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),ui,vi,wi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,17,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                    IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,18,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,19,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,20,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)

      IDX3(DmatrixAtEdge,21,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,22,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,23,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,24,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      IDX3(DmatrixAtEdge,25,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX51_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,6,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,7,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,9,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,10,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX52_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,11,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,12,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,13,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,14,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,15,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX53_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,16,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,17,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,18,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,19,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,20,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX54_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)

      IDX3(DmatrixAtEdge,21,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX15_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,22,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX25_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,23,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX35_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,24,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX45_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
      IDX3(DmatrixAtEdge,25,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX55_3D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),
                                     IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),ui,vi,wi,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate scalar artificial dissipation of Rusanov-type
      !---------------------------------------------------------------------------

      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*&
          (Ei-RCONST(0.5)*(ui*ui+vi*vi+wi*wi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*&
          (Ej-RCONST(0.5)*(uj*uj+vj*vj+wj*wj)), SYS_EPSREAL_DP))
      
      ! Compute dissipation tensor
      aux = dscale *&
          max( abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*uj+&
                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*vj+&
                   IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*wj) +&
          sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)+&
               POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2)+&
               POW(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),2))*cj,&
               abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*ui+&
                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*vi+&
                   IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*wi) +&
          sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)+&
               POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2)+&
               POW(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),2))*ci )

      IDX3(DmatrixAtEdge, :,1,idx,_,_,_) = RCONST(0.0)
      IDX3(DmatrixAtEdge, 1,1,idx,_,_,_) = aux
      IDX3(DmatrixAtEdge, 6,1,idx,_,_,_) = aux
      IDX3(DmatrixAtEdge,11,1,idx,_,_,_) = aux
      IDX3(DmatrixAtEdge,16,1,idx,_,_,_) = aux
    end do
  
  end subroutine hydro_calcMatRusDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcCharacteristics3d_sim(Dweight, DdataAtEdge,&
      nedges, DcharVariablesAtEdge, DeigenvaluesAtEdge,&
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

  pure subroutine hydro_calcFluxFCTScDiss3d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)


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
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IedgeList

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
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,anorm,aux,c_ij,d_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    integer :: idx


    do idx = 1, nedges

      ! Compute skew-symmetric coefficient
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                  IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute densities
      ri = DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      rj = DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute enthalpies
      hi = (TOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
      hj = (TOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(ri,rj)
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      w_ij = ROE_MEAN_VALUE(wi,wj,aux)
      H_ij = ROE_MEAN_VALUE(hi,hj,aux)

      ! Compute auxiliary variables
      vel_ij = u_ij*a(1) + v_ij*a(2) + w_ij*a(3)
      q_ij   = RCONST(0.5)*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

      ! Compute the speed of sound
      c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))

      ! Compute scalar dissipation
      d_ij = abs(vel_ij) + anorm*c_ij

      ! Compute antidiffusive flux
       IDX2(DfluxesAtEdge,:,idx,_,_) = dscale*d_ij*&
           (IDX3(DdataAtEdge,:,1,idx,_,_,_)-IDX3(DdataAtEdge,:,2,idx,_,_,_))
     end do

   end subroutine hydro_calcFluxFCTScDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRoeDiss3d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for FCT
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
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IedgeList

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
    real(DP), dimension(NVAR3D) :: Diff
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,c2_ij,c_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    real(DP) :: anorm,aux,aux1,aux2
    real(DP) :: l1,l2,l3,l4,l5,w1,w2,w3,w4,w5
    integer :: idx


    do idx = 1, nedges
      
      ! Compute the skew-symmetric coefficient and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Normalise the skew-symmetric coefficient
        a = a/anorm

        ! Compute velocities
        ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        
        uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
        vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
        wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute pressures
        pi = PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        pj = PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute densities
        ri = DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute enthalpies
        hi = (TOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2) + w_ij*a(3)
        q_ij   = RCONST(0.5)*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

        ! Compute the speed of sound
        c2_ij  = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij   = sqrt(c2_ij)

        ! Compute eigenvalues
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)
        l5 = abs(vel_ij)

        ! Compute solution difference U_i-U_j
        Diff = IDX3(DdataAtEdge,:,1,idx,_,_,_)-&
               IDX3(DdataAtEdge,:,2,idx,_,_,_)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-RCONST(1.0))/RCONST(2.0)/c2_ij*(q_ij*Diff(1)&
                                                             -u_ij*Diff(2)&
                                                             -v_ij*Diff(3)&
                                                             -w_ij*Diff(4)&
                                                                  +Diff(5))
        aux2 = (vel_ij*Diff(1)&
                 -a(1)*Diff(2)&
                 -a(2)*Diff(3)&
                 -a(3)*Diff(4))/RCONST(2.0)/c_ij

        ! Get the dimension with largest coefficient
        select case(maxloc(a,1))
        case(1)
          ! Compute characteristic variables multiplied by the corresponding eigenvalue
          w1 = l1 * (aux1 + aux2)
          w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*Diff(1)&
                                      +((HYDRO_GAMMA)-RCONST(1.0))*( u_ij*Diff(2)&
                                                                    +v_ij*Diff(3)&
                                                                    +w_ij*Diff(4)&
                                                                         -Diff(5))/c2_ij)
          w3 = l3 * (aux1 - aux2)
          w4 = l4 * ( (v_ij-vel_ij*a(2))/a(1)*Diff(1)&
                                        +a(2)*Diff(2)&
                +(a(2)*a(2)-RCONST(1.0))/a(1)*Diff(3)&
                              +a(2)*a(3)/a(1)*Diff(4))
          w5 = l5 * ( (vel_ij*a(3)-w_ij)/a(1)*Diff(1)&
                                        -a(3)*Diff(2)&
                              -a(2)*a(3)/a(1)*Diff(3)&
                +(RCONST(1.0)-a(3)*a(3))/a(1)*Diff(4))

          ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
          Diff(1) = anorm * ( w1 + w2 + w3 )
          Diff(2) = anorm * ( (u_ij-c_ij*a(1))*w1 + u_ij*w2 +&
                              (u_ij+c_ij*a(1))*w3 + a(2)*w4 - a(3)*w5 )
          Diff(3) = anorm * ( (v_ij-c_ij*a(2))*w1 + v_ij*w2 +&
                              (v_ij+c_ij*a(2))*w3 - a(1)*w4 )
          Diff(4) = anorm * ( (w_ij-c_ij*a(3))*w1 + w_ij*w2 +&
                              (w_ij+c_ij*a(3))*w3 + a(1)*w5 )
          Diff(5) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3&
                              + (u_ij*a(2)-v_ij*a(1))*w4 + (w_ij*a(1)-u_ij*a(3))*w5 )

        case(2)
          ! Compute characteristic variables multiplied by the corresponding eigenvalue
          w1 = l1 * (aux1 + aux2)
          w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*Diff(1)&
                                       +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*Diff(2)&
                                                                    +v_ij*Diff(3)&
                                                                    +w_ij*Diff(4)&
                                                                         -Diff(5))/c2_ij)
          w3 = l3 * (aux1 - aux2)
          w4 = l4 * ( (vel_ij*a(1)-u_ij)/a(2)*Diff(1)&
                +(RCONST(1.0)-a(1)*a(1))/a(2)*Diff(2)&
                                        -a(1)*Diff(3)&
                              -a(1)*a(3)/a(2)*Diff(4))
          w5 = l5 * ( (w_ij-vel_ij*a(3))/a(2)*Diff(1)&
                              +a(1)*a(3)/a(2)*Diff(2)&
                                        +a(3)*Diff(3)&
                +(a(3)*a(3)-RCONST(1.0))/a(2)*Diff(4))

          ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
          Diff(1) = anorm * ( w1 + w2 + w3 )
          Diff(2) = anorm * ( (u_ij-c_ij*a(1))*w1 + u_ij*w2 +&
                              (u_ij+c_ij*a(1))*w3 + a(2)*w4 )
          Diff(3) = anorm * ( (v_ij-c_ij*a(2))*w1 + v_ij*w2 +&
                              (v_ij+c_ij*a(2))*w3 - a(1)*w4 + a(3)*w5 )
          Diff(4) = anorm * ( (w_ij-c_ij*a(3))*w1 + w_ij*w2 +&
                              (w_ij+c_ij*a(3))*w3 - a(2)*w5 )
          Diff(5) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3&
                              + (u_ij*a(2)-v_ij*a(1))*w4 + (v_ij*a(3)-w_ij*a(2))*w5 )

        case(3)
          ! Compute characteristic variables multiplied by the corresponding eigenvalue
          w1 = l1 * (aux1 + aux2)
          w2 = l2 * ((RCONST(1.0)-((HYDRO_GAMMA)-RCONST(1.0))*q_ij/c2_ij)*Diff(1)&
                                       +((HYDRO_GAMMA)-RCONST(1.0))*(u_ij*Diff(2)&
                                                                    +v_ij*Diff(3)&
                                                                    +w_ij*Diff(4)&
                                                                         -Diff(5))/c2_ij)
          w3 = l3 * (aux1 - aux2)
          w4 = l4 * ( (u_ij-vel_ij*a(1))/a(3)*Diff(1)&
                +(a(1)*a(1)-RCONST(1.0))/a(3)*Diff(2)&
                              +a(1)*a(2)/a(3)*Diff(3)&
                                        +a(1)*Diff(4) )
          w5 = l5 * ( (vel_ij*a(2)-v_ij)/a(3)*Diff(1)&
                              -a(1)*a(2)/a(3)*Diff(2)&
                +(RCONST(1.0)-a(2)*a(2))/a(3)*Diff(3)&
                                        -a(2)*Diff(4) )

          ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
          Diff(1) = anorm * ( w1 + w2 + w3 )
          Diff(2) = anorm * ( (u_ij-c_ij*a(1))*w1 + u_ij*w2 +&
                              (u_ij+c_ij*a(1))*w3 - a(3)*w4 )
          Diff(3) = anorm * ( (v_ij-c_ij*a(2))*w1 + v_ij*w2 +&
                              (v_ij+c_ij*a(2))*w3 + a(3)*w5 )
          Diff(4) = anorm * ( (w_ij-c_ij*a(3))*w1 + w_ij*w2 +&
                              (w_ij+c_ij*a(3))*w3 + a(1)*w4 - a(2)*w5)
          Diff(5) = anorm * ( (H_ij-c_ij*vel_ij)*w1 + q_ij*w2 + (H_ij+c_ij*vel_ij)*w3&
                              + (w_ij*a(1)-u_ij*a(3))*w4 + (v_ij*a(3)-w_ij*a(2))*w5 )
        end select

        ! Compute antidiffusive flux
        IDX2(DfluxesAtEdge,:,idx,_,_) = dscale*Diff
      else
        ! Cancel antidiffusive flux
        IDX2(DfluxesAtEdge,:,idx,_,_) = 0
      end if
    end do

  end subroutine hydro_calcFluxFCTRoeDiss3d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRusDiss3d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

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
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IedgeList

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
    real(DP) :: Ei,Ej,ci,cj,ui,uj,vi,vj,wi,wj
    real(DP) :: d_ij
    integer :: idx

    
    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*&
          (Ei-RCONST(0.5)*(ui*ui+vi*vi+wi*wi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(HYDRO_GAMMA)*&
          (Ej-RCONST(0.5)*(uj*uj+vj*vj+wj*wj)), SYS_EPSREAL_DP))
      
#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*uj+&
                      RCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*vj+&
                      RCONST(0.5)*(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,3,2,idx,_,_,_))*wj)+&
              RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)+&
                               POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2)+&
                               POW(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),2))*cj,&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_))*ui+&
                      RCONST(0.5)*(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_))*vi+&
                      RCONST(0.5)*(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,3,1,idx,_,_,_))*wi)+&
              RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)+&
                               POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2)+&
                               POW(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),2))*ci )
#else
      ! Compute scalar dissipation
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*uj+&
                      IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*vj+&
                      IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*wj)+&
             sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)+&
                  POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2)+&
                  POW(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),2))*cj,&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*ui+&
                      IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*vi+&
                      IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*wi)+&
             sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)+&
                  POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2)+&
                  POW(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),2))*ci )
#endif

      ! Compute antidiffusive flux
      IDX2(DfluxesAtEdge,:,idx,_,_) = dscale*d_ij*&
          (IDX3(DdataAtEdge,:,1,idx,_,_,_)-IDX3(DdataAtEdge,:,2,idx,_,_,_))
    end do

  end subroutine hydro_calcFluxFCTRusDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDensity3d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density in 3D.
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
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_)
    end do

  end subroutine hydro_trafoFluxDensity3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDensity3d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density in 3D.
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
      IDX2(DtransformedDataAtEdge,1,idx,_,_) =&
          DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do

  end subroutine hydro_trafoDiffDensity3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalDensity3d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the density in 3D.
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          DENSITY2_3D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine hydro_trafoNodalDensity3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxEnergy3d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the energy in 3D.
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
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          TOTALENERGY2_3D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -TOTALENERGY2_3D(DfluxesAtEdge,IDX2,idx,_,_)
    end do

  end subroutine hydro_trafoFluxEnergy3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffEnergy3d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the energy in 3D.
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
      IDX2(DtransformedDataAtEdge,1,idx,_,_) =&
          TOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          TOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do

  end subroutine hydro_trafoDiffEnergy3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalEnergy3d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the energy in 3D.
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          TOTALENERGY2_3D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine hydro_trafoNodalEnergy3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxPressure3d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the pressure in 3D.
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
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          ((HYDRO_GAMMA)-RCONST(1.0))*&
          (RCONST(0.5)*(ui*ui+vi*vi+wi*wi)*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                      ui*XMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                      vi*YMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                      wi*ZMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)+&
                                       TOTALENERGY2_3D(DfluxesAtEdge,IDX2,idx,_,_))
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -((HYDRO_GAMMA)-RCONST(1.0))*&
         (RCONST(0.5)*(uj*uj+vj*vj+wj*wj)*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                     uj*XMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                     vj*YMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                     wj*ZMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)+&
                                      TOTALENERGY2_3D(DfluxesAtEdge,IDX2,idx,_,_))
    end do

  end subroutine hydro_trafoFluxPressure3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffPressure3d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the pressure in 3D.
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
      IDX2(DtransformedDataAtEdge,1,idx,_,_) =&
          PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do

  end subroutine hydro_trafoDiffPressure3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalPressure3d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the pressure in 3D.
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          PRESSURE2_3D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine hydro_trafoNodalPressure3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxVelocity3d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the x-, y- and z-velocity
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
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      
      ! Transformed velocity fluxes in x-direction
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          (XMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
          ui*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_))/&
          DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -(XMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
          uj*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_))/&
          DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed velocity fluxes in y-direction
      IDX3(DtransformedFluxesAtEdge,2,1,idx,_,_,_) =&
          (YMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
          vi*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_))/&
          DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,_,_,_) =&
         -(YMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
          vj*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_))/&
          DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed velocity fluxes in z-direction
      IDX3(DtransformedFluxesAtEdge,3,1,idx,_,_,_) =&
          (ZMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
          wi*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_))/&
          DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      IDX3(DtransformedFluxesAtEdge,3,2,idx,_,_,_) =&
         -(ZMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
          wj*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_))/&
          DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
    end do
    
  end subroutine hydro_trafoFluxVelocity3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffVelocity3d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the x-, y- and z-velocity
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
      IDX2(DtransformedDataAtEdge,1,idx,_,_) =&
          XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed velocity difference in y-direction
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed velocity difference in z-direction
      IDX2(DtransformedDataAtEdge,3,idx,_,_) =&
          ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do
    
  end subroutine hydro_trafoDiffVelocity3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalVelocity3d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the velocity in 3D.
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          XVELOCITY2_3D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed y-velocity values
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          YVELOCITY2_3D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed z-velocity values
      IDX2(DtransformedDataAtNode,3,idx,_,_) =&
          ZVELOCITY2_3D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine hydro_trafoNodalVelocity3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxMomentum3d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the x-, y- and z-momentum
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
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          XMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -XMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)

      ! Transformed momentum fluxes in y-direction
      IDX3(DtransformedFluxesAtEdge,2,1,idx,_,_,_) =&
          YMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,_,_,_) =&
         -YMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)

      ! Transformed momentum fluxes in z-direction
      IDX3(DtransformedFluxesAtEdge,3,1,idx,_,_,_) =&
          ZMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,3,2,idx,_,_,_) =&
         -ZMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)
    end do
    
  end subroutine hydro_trafoFluxMomentum3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffMomentum3d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the x-, y- and z-momentum
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
      IDX2(DtransformedDataAtEdge,1,idx,_,_) =&
          XMOMENTUM3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          XMOMENTUM3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed momentum difference in y-direction
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          YMOMENTUM3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          YMOMENTUM3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed momentum difference in z-direction
      IDX2(DtransformedDataAtEdge,3,idx,_,_) =&
          ZMOMENTUM3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          ZMOMENTUM3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do
    
  end subroutine hydro_trafoDiffMomentum3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalMomentum3d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the momentum in 3D.
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          XMOMENTUM2_3D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed y-momentum values
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          YMOMENTUM2_3D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed z-momentum values
      IDX2(DtransformedDataAtNode,3,idx,_,_) =&
          ZMOMENTUM2_3D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine hydro_trafoNodalMomentum3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenEng3d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density and energy in 3D.
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
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_)

      ! Transformed total energy fluxes
      IDX3(DtransformedFluxesAtEdge,2,1,idx,_,_,_) =&
          TOTALENERGY2_3D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,_,_,_) =&
         -TOTALENERGY2_3D(DfluxesAtEdge,IDX2,idx,_,_)
    end do

  end subroutine hydro_trafoFluxDenEng3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenEng3d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density and energy in 3D.
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
      IDX2(DtransformedDataAtEdge,1,idx,_,_) =&
          DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed total energy difference
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          TOTALENERGY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          TOTALENERGY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do

  end subroutine hydro_trafoDiffDenEng3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalDenEng3d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the density and
    ! energy in 3D.
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          DENSITY2_3D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed energy values
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          TOTALENERGY2_3D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine hydro_trafoNodalDenEng3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenPre3d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density and energy in 3D.
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
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed density fluxes
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_)
      
      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,2,1,idx,_,_,_) =&
          ((HYDRO_GAMMA)-RCONST(1.0))*&
          (RCONST(0.5)*(ui*ui+vi*vi+wi*wi)*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                      ui*XMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                      vi*YMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                      wi*ZMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)+&
                                       TOTALENERGY2_3D(DfluxesAtEdge,IDX2,idx,_,_))
      IDX3(DtransformedFluxesAtEdge,2,2,idx,_,_,_) =&
         -((HYDRO_GAMMA)-RCONST(1.0))*&
         (RCONST(0.5)*(uj*uj+vj*vj+wj*wj)*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                     uj*XMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                     vj*YMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                     wj*ZMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)+&
                                      TOTALENERGY2_3D(DfluxesAtEdge,IDX2,idx,_,_))
    end do

  end subroutine hydro_trafoFluxDenPre3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenPre3d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density and energy in 3D.
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
      IDX2(DtransformedDataAtEdge,1,idx,_,_) =&
          DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      
      ! Transformed pressure difference
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do
    
  end subroutine hydro_trafoDiffDenPre3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalDenPre3d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the density and
    ! pressure in 3D.
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          DENSITY2_3D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed pressure values
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          PRESSURE2_3D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine hydro_trafoNodalDenPre3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenPreVel3d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation
    ! of the given flux into primitive variables in 3D.
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
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Transformed density fluxes
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_)

      ! Transformed velocity fluxes in x-direction
      IDX3(DtransformedFluxesAtEdge,2,1,idx,_,_,_) =&
          (XMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
          ui*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_))/&
          DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,_,_,_) =&
         -(XMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
          uj*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_))/&
          DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed velocity fluxes in y-direction
      IDX3(DtransformedFluxesAtEdge,3,1,idx,_,_,_) =&
          (YMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
          vi*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_))/&
          DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      IDX3(DtransformedFluxesAtEdge,3,2,idx,_,_,_) =&
         -(YMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
          vj*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_))/&
          DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed velocity fluxes in z-direction
      IDX3(DtransformedFluxesAtEdge,4,1,idx,_,_,_) =&
          (ZMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
          wi*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_))/&
          DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      IDX3(DtransformedFluxesAtEdge,4,2,idx,_,_,_) =&
         -(ZMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
          wj*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_))/&
          DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,5,1,idx,_,_,_) =&
          ((HYDRO_GAMMA)-RCONST(1.0))*&
          (RCONST(0.5)*(ui*ui+vi*vi+wi*wi)*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                      ui*XMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                      vi*YMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                      wi*ZMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)+&
                                       TOTALENERGY2_3D(DfluxesAtEdge,IDX2,idx,_,_))
      IDX3(DtransformedFluxesAtEdge,5,2,idx,_,_,_) =&
         -((HYDRO_GAMMA)-RCONST(1.0))*&
         (RCONST(0.5)*(uj*uj+vj*vj+wj*wj)*DENSITY2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                     uj*XMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                     vj*YMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)-&
                                     wj*ZMOMENTUM2_3D(DfluxesAtEdge,IDX2,idx,_,_)+&
                                      TOTALENERGY2_3D(DfluxesAtEdge,IDX2,idx,_,_))
    end do

  end subroutine hydro_trafoFluxDenPreVel3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenPreVel3d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density, pressure and velocity in 3D.
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
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed velocity difference in x-direction
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      
      ! Transformed velocity difference in y-direction
      IDX2(DtransformedDataAtEdge,3,idx,_,_) =&
          YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed velocity difference in z-direction
      IDX2(DtransformedDataAtEdge,4,idx,_,_) =&
          ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed pressure difference
      IDX2(DtransformedDataAtEdge,5,idx,_,_) =&
          PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do

  end subroutine hydro_trafoDiffDenPreVel3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalDenPreVel3d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the density,
    ! pressure and velocity in 3D.
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          DENSITY2_3D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed x-velocity values
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          XVELOCITY2_3D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed y-velocity values
      IDX2(DtransformedDataAtNode,3,idx,_,_) =&
          YVELOCITY2_3D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed z-velocity values
      IDX2(DtransformedDataAtNode,4,idx,_,_) =&
          ZVELOCITY2_3D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed pressure values
      IDX2(DtransformedDataAtNode,5,idx,_,_) =&
          PRESSURE2_3D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine hydro_trafoNodalDenPreVel3d_sim

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcBoundaryvalues3d(DbdrNormal, DpointNormal,&
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

  end subroutine hydro_calcBoundaryvalues3d

end module hydro_callback3d
