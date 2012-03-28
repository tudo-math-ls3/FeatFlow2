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
!# 1.) mhd_calcFluxGalerkin3d_sim
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
!# 12.) mhd_calcMatGalerkin3d_sim
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
!# 25.) mhd_trafoNodalDensity3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density values
!#
!# 26.) mhd_trafoFluxEnergy3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 27.) mhd_trafoDiffEnergy3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 28.) mhd_trafoNodalEnergy3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal energy values
!#
!# 29.) mhd_trafoFluxPressure3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 30.) mhd_trafoDiffPressure3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the presure
!#
!# 31.) mhd_trafoNodalPressure3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal pressure values
!#
!# 32.) mhd_trafoFluxVelocity3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 33.) mhd_trafoDiffVelocity3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 34.) mhd_trafoNodalVelocity3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal velocity values
!#
!# 35.) mhd_trafoFluxMomentum3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 36.) mhd_trafoDiffMomentum3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 37.) mhd_trafoNodalMomentum3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal momentum values
!#
!# 38.) mhd_trafoFluxDenEng3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 39.) mhd_trafoDiffDenEng3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 40.) mhd_trafoNodalDenEng3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density and energy values
!#
!# 41.) mhd_trafoFluxDenPre3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 42.) mhd_trafoDiffDenPre3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 43.) mhd_trafoNodalDenPre3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density and pressure values
!#
!# 44.) mhd_trafoFluxDenPreVel3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 45.) mhd_trafoDiffDenPreVel3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure
!#         and the velocity
!#
!# 46.) mhd_trafoNodalDenPreVel3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density, pressure and velocity values
!#
!# 47.) mhd_trafoDiffMagfield3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the magnetic field
!#
!# 48.) mhd_trafoFluxMagfield3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the magnetic field
!#
!# 49.) mhd_trafoNodalMagfield3d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal values for the magnetic field
!#
!# 50.) mhd_calcBoundaryvalues3d
!#      -> Computes the boundary values for a given node
!#
!# </purpose>
!##############################################################################

module mhd_callback3d

#define MHD_NDIM 3
#include "mhd.h"

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

  ! Modules from MHD model
  use mhd_basic

  ! The following 3D-routine coincide with their 2D-versions
  use mhd_callback2d, only : &
      mhd_trafoFluxDensity3d_sim => mhd_trafoFluxDensity2d_sim,&
      mhd_trafoDiffDensity3d_sim => mhd_trafoDiffDensity2d_sim,&
      mhd_trafoFluxEnergy3d_sim => mhd_trafoFluxEnergy2d_sim,&
      mhd_trafoDiffEnergy3d_sim => mhd_trafoDiffEnergy2d_sim,&
      mhd_trafoFluxPressure3d_sim => mhd_trafoFluxPressure2d_sim,&
      mhd_trafoDiffPressure3d_sim => mhd_trafoDiffPressure2d_sim,&
      mhd_trafoFluxVelocity3d_sim => mhd_trafoFluxVelocity2d_sim,&
      mhd_trafoDiffVelocity3d_sim => mhd_trafoDiffVelocity2d_sim,&
      mhd_trafoFluxMomentum3d_sim => mhd_trafoFluxMomentum2d_sim,&
      mhd_trafoDiffMomentum3d_sim => mhd_trafoDiffMomentum2d_sim,&
      mhd_trafoFluxDenEng3d_sim => mhd_trafoFluxDenEng2d_sim,&
      mhd_trafoDiffDenEng3d_sim => mhd_trafoDiffDenEng2d_sim,&
      mhd_trafoFluxDenPre3d_sim => mhd_trafoFluxDenPre2d_sim,&
      mhd_trafoDiffDenPre3d_sim => mhd_trafoDiffDenPre2d_sim,&
      mhd_trafoFluxDenPreVel3d_sim => mhd_trafoFluxDenPreVel2d_sim,&
      mhd_trafoDiffDenPreVel3d_sim => mhd_trafoDiffDenPreVel2d_sim,&
      mhd_trafoDiffMagfield3d_sim => mhd_trafoDiffMagfield2d_sim,&
      mhd_trafoFluxMagfield3d_sim => mhd_trafoFluxMagfield2d_sim,&
      mhd_trafoNodalDensity3d_sim => mhd_trafoNodalDensity2d_sim,&
      mhd_trafoNodalEnergy3d_sim => mhd_trafoNodalEnergy2d_sim,&
      mhd_trafoNodalPressure3d_sim => mhd_trafoNodalPressure2d_sim,&
      mhd_trafoNodalVelocity3d_sim => mhd_trafoNodalVelocity2d_sim,&
      mhd_trafoNodalMomentum3d_sim => mhd_trafoNodalMomentum2d_sim,&
      mhd_trafoNodalDenEng3d_sim => mhd_trafoNodalDenEng2d_sim,&
      mhd_trafoNodalDenPre3d_sim => mhd_trafoNodalDenPre2d_sim,&
      mhd_trafoNodalDenPreVel3d_sim => mhd_trafoNodalDenPreVel2d_sim,&
      mhd_trafoNodalMagfield3d_sim => mhd_trafoNodalMagfield2d_sim

  implicit none

  private
  public :: mhd_calcFluxGalerkin3d_sim
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
  public :: mhd_calcMatGalerkin3d_sim
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
  public :: mhd_trafoNodalDensity3d_sim
  public :: mhd_trafoNodalEnergy3d_sim
  public :: mhd_trafoNodalPressure3d_sim
  public :: mhd_trafoNodalVelocity3d_sim
  public :: mhd_trafoNodalMomentum3d_sim
  public :: mhd_trafoNodalDenEng3d_sim
  public :: mhd_trafoNodalDenPre3d_sim
  public :: mhd_trafoNodalDenPreVel3d_sim
  public :: mhd_trafoNodalMagfield3d_sim
  public :: mhd_calcBoundaryvalues3d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGalerkin3d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      wi = ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      wj = ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = TOTALPRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi)
      qj = MAG_DOT_VEL3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      Fxi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fxj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      Fyi(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fyj(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      Fzi(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fzj(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      
      ! Assemble skew-symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
           IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fyj+&
           IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*Fzj-&
           IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
           IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fyi-&
           IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fzi )
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
      ! Compute flux difference for x-direction
      Fx_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      Fy_ij(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      Fz_ij(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Assemble fluxes
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)*Fz_ij)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*Fz_ij)
#endif
    end do

  end subroutine mhd_calcFluxGalerkin3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGalNoBdr3d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
  real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
  integer :: idx
  
    
    do idx = 1, nedges

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      wi = ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      wj = ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = TOTALPRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi)
      qj = MAG_DOT_VEL3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj)

      ! Compute flux difference for x-direction
      Fx_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      Fy_ij(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      Fz_ij(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Assemble fluxes
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
          (RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                        IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*Fx_ij+&
           RCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)-&
                        IDX3(DcoeffsAtEdge,2,2,idx,0,0,0))*Fy_ij+&
           RCONST(0.5)*(IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)-&
                        IDX3(DcoeffsAtEdge,3,2,idx,0,0,0))*Fz_ij)
      DfluxesAtEdge(:,2,idx) = DfluxesAtEdge(:,1,idx)
    end do

  end subroutine mhd_calcFluxGalNoBdr3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxScDiss3d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      wi = ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      wj = ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = TOTALPRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi)
      qj = MAG_DOT_VEL3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj)


#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      Fxi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fxj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      Fyi(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fyj(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      Fzi(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fzj(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      Fx_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      Fy_ij(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      Fz_ij(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar artificial dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------

      Diff = RCONST(0.0)

      !!! TODO !!!

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
           IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fyj+&
           IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*Fzj-&
           IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
           IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fyi-&
           IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)*Fzi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)*Fz_ij+ Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*Fz_ij+ Diff)
#endif
    end do

  end subroutine mhd_calcFluxScDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxScDissDiSp3d_sim(DdataAtEdge, DcoeffsAtEdge,&
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

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      wi = ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      wj = ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = TOTALPRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi)
      qj = MAG_DOT_VEL3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      Fxi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fxj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      Fyi(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fyj(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      Fzi(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fzj(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      Fx_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      Fy_ij(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      Fz_ij(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar artificial dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,0,0,0)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,0,0,0))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      Diff = RCONST(0.0)

      !!! TODO !!!

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
           IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fyj+&
           IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*Fzj-&
           IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
           IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fyi-&
           IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)*Fzi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)*Fz_ij+ Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij+&
           IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*Fz_ij+ Diff)
#endif
    end do

  end subroutine mhd_calcFluxScDissDiSp3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRoeDiss3d_sim(DdataAtEdge, DcoeffsAtEdge,&
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

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      wi = ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      wj = ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = TOTALPRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi)
      qj = MAG_DOT_VEL3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      Fxi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fxj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      Fyi(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fyj(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      Fzi(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fzj(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      Fx_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      Fy_ij(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      Fz_ij(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,0,0,0)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,0,0,0))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Normalise the skew-symmetric coefficient
        a = a/anorm
        
        Diff = RCONST(0.0)

        !! TODO !!

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef MHD_USE_IBP
       IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
             IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fyj+&
             IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*Fzj-&
             IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
             IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fyi-&
             IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)*Fzi + Diff)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
            (IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)*Fz_ij + Diff)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*Fz_ij + Diff)
#endif
      else
        
#ifdef MHD_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
             IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fyj+&
             IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*Fzj-&
             IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
             IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fyi-&
             IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)*Fzi)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
            (IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)*Fz_ij)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*Fz_ij)
#endif
      end if
    end do

  end subroutine mhd_calcFluxRoeDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRoeDissDiSp3d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
    

    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      wi = ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      wj = ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = TOTALPRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi)
      qj = MAG_DOT_VEL3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      Fxi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fxj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      Fyi(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fyj(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      Fzi(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fzj(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      Fx_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      Fy_ij(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      Fz_ij(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#endif
      
      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,0,0,0)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,0,0,0))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Normalise the skew-symmetric coefficient
        a = a/anorm
        
        DiffX = RCONST(0.0); DiffY = RCONST(0.0); DiffZ = RCONST(0.0)

        !!! TODO !!!

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef MHD_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
             IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fyj+&
             IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*Fzj-&
             IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
             IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fyi-&
             IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)*Fzi+&
             DiffX+DiffY+DiffZ)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
            (IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)*Fz_ij+&
             DiffX+DiffY+DiffZ)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*Fz_ij+&
             DiffX+DiffY+DiffZ)
#endif
      else
        
#ifdef MHD_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
             IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fyj+&
             IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*Fzj-&
             IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
             IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fyi-&
             IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)*Fzi)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
            (IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)*Fz_ij)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij+&
             IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*Fz_ij)
#endif
      end if
    end do

  end subroutine mhd_calcFluxRoeDissDiSp3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRusDiss3d_sim(DdataAtEdge, DcoeffsAtEdge,&
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

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      wi = ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      wj = ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = TOTALPRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi)
      qj = MAG_DOT_VEL3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      Fxi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fxj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      Fyi(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fyj(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      Fzi(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fzj(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      Fx_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      Fy_ij(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      Fz_ij(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
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
      ca1i = abs(XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0))
      ca1j = abs(XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0))

      ! Compute the speed of the Alfven waves in y-direction
      ca2i = abs(YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0))
      ca2j = abs(YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0))

      ! Compute the speed of the Alfven waves in z-direction
      ca3i = abs(ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0))
      ca3j = abs(ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0))

      ! Compute the speed of sound
      aPow2i = RCONST(MAGNETOHYDRODYN_GAMMA)*&
               PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)/&
               DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      aPow2j = RCONST(MAGNETOHYDRODYN_GAMMA)*&
               PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)/&
               DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute auxiliary quantities
      astPow2i = MAGFIELDMAGNITUDE3(DdataAtEdge,IDX3,1,idx,0,0,0)/&
                 DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0) + aPow2i
      astPow2j = MAGFIELDMAGNITUDE3(DdataAtEdge,IDX3,2,idx,0,0,0)/&
                 DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0) + aPow2j

      ! Compute the speed of the fast waves in x-direction
      cf1i = sqrt(RCONST(0.5)*(astPow2i+&
                  sqrt(POW(astPow2i,2)-RCONST(4.0)*aPow2i*POW(ca1i,2))))
      cf1j = sqrt(RCONST(0.5)*(astPow2j+&
                  sqrt(POW(astPow2j,2)-RCONST(4.0)*aPow2j*POW(ca1j,2))))

      ! Compute the speed of the fast waves in y-direction
      cf2i = sqrt(RCONST(0.5)*(astPow2i+&
                  sqrt(POW(astPow2i,2)-RCONST(4.0)*aPow2i*POW(ca2i,2))))
      cf2j = sqrt(RCONST(0.5)*(astPow2j+&
                  sqrt(POW(astPow2j,2)-RCONST(4.0)*aPow2j*POW(ca2j,2))))

      ! Compute the speed of the fast waves in z-direction
      cf3i = sqrt(RCONST(0.5)*(astPow2i+&
                  sqrt(POW(astPow2i,2)-RCONST(4.0)*aPow2i*POW(ca3i,2))))
      cf3j = sqrt(RCONST(0.5)*(astPow2j+&
                  sqrt(POW(astPow2j,2)-RCONST(4.0)*aPow2j*POW(ca3j,2))))

#ifdef MHD_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*uj+&
                      RCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,2,2,idx,0,0,0))*vj+&
                      RCONST(0.5)*(IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,3,2,idx,0,0,0))*wj)+&
                 RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-
                                      IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),2)*cf1j+&
                                  POW(IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)-
                                      IDX3(DcoeffsAtEdge,2,2,idx,0,0,0),2)*cf2j+&
                                  POW(IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)-
                                      IDX3(DcoeffsAtEdge,3,2,idx,0,0,0),2)*cf3j),&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,0,0,0))*ui+&
                      RCONST(0.5)*(IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,2,1,idx,0,0,0))*vi+&
                      RCONST(0.5)*(IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,3,1,idx,0,0,0))*wi)+&
                 RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)-
                                      IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),2)*cf1i+&
                                  POW(IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)-
                                      IDX3(DcoeffsAtEdge,2,1,idx,0,0,0),2)*cf2i+&
                                  POW(IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)-
                                      IDX3(DcoeffsAtEdge,3,1,idx,0,0,0),2)*cf3i) )
#else
       ! Compute scalar dissipation
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*uj+&
                      IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*vj+&
                      IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)*wj)+&
                 sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),2)*cf1j+&
                      POW(IDX3(DcoeffsAtEdge,2,1,idx,0,0,0),2)*cf2j+&
                      POW(IDX3(DcoeffsAtEdge,3,1,idx,0,0,0),2)*cf3j),&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*ui+&
                      IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*vi+&
                      IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*wi)+&
                 sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),2)*cf1i+&
                      POW(IDX3(DcoeffsAtEdge,2,2,idx,0,0,0),2)*cf2i+&
                      POW(IDX3(DcoeffsAtEdge,3,2,idx,0,0,0),2)*cf3i) )
#endif

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
                   IDX3(DdataAtEdge,:,1,idx,0,0,0))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------
      
#ifdef MHD_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
           IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fyj-&
           IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
           IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fyi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij + Diff)
#endif
    end do

  end subroutine mhd_calcFluxRusDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRusDissDiSp3d_sim(DdataAtEdge, DcoeffsAtEdge,&
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

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      wi = ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      wj = ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = TOTALPRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi)
      qj = MAG_DOT_VEL3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      Fxi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fxi(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fxj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fxj(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      Fyi(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fyi(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fyj(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fyj(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      Fzi(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fzi(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fzj(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fzj(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      Fx_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fx_ij(8) = INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      Fy_ij(1) = INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(2) = INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(3) = INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(4) = INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(5) = INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(6) = INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(7) = INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fy_ij(8) = INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_YDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      Fz_ij(1) = INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(2) = INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(3) = INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(4) = INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(5) = INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(6) = INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(7) = INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fz_ij(8) = INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_ZDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
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
      ca1i = abs(XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0))
      ca1j = abs(XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0))

      ! Compute the speed of the Alfven waves in y-direction
      ca2i = abs(YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0))
      ca2j = abs(YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0))

      ! Compute the speed of the Alfven waves in z-direction
      ca3i = abs(ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0))
      ca3j = abs(ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0))

      ! Compute the speed of sound
      aPow2i = RCONST(MAGNETOHYDRODYN_GAMMA)*&
               PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)/&
               DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      aPow2j = RCONST(MAGNETOHYDRODYN_GAMMA)*&
               PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)/&
               DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute auxiliary quantities
      astPow2i = MAGFIELDMAGNITUDE3(DdataAtEdge,IDX3,1,idx,0,0,0)/&
                 DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0) + aPow2i
      astPow2j = MAGFIELDMAGNITUDE3(DdataAtEdge,IDX3,2,idx,0,0,0)/&
                 DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0) + aPow2j

      ! Compute the speed of the fast waves in x-direction
      cf1i = sqrt(RCONST(0.5)*(astPow2i+&
                  sqrt(POW(astPow2i,2)-RCONST(4.0)*aPow2i*POW(ca1i,2))))
      cf1j = sqrt(RCONST(0.5)*(astPow2j+&
                  sqrt(POW(astPow2j,2)-RCONST(4.0)*aPow2j*POW(ca1j,2))))

      ! Compute the speed of the fast waves in y-direction
      cf2i = sqrt(RCONST(0.5)*(astPow2i+&
                  sqrt(POW(astPow2i,2)-RCONST(4.0)*aPow2i*POW(ca2i,2))))
      cf2j = sqrt(RCONST(0.5)*(astPow2j+&
                  sqrt(POW(astPow2j,2)-RCONST(4.0)*aPow2j*POW(ca2j,2))))

      ! Compute the speed of the fast waves in z-direction
      cf3i = sqrt(RCONST(0.5)*(astPow2i+&
                  sqrt(POW(astPow2i,2)-RCONST(4.0)*aPow2i*POW(ca3i,2))))
      cf3j = sqrt(RCONST(0.5)*(astPow2j+&
                  sqrt(POW(astPow2j,2)-RCONST(4.0)*aPow2j*POW(ca3j,2))))
      
#ifdef MHD_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*uj)+&
              RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-
                                   IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),2))*cf1j+&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,2,2,idx,0,0,0))*vj)+&
              RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)-
                                   IDX3(DcoeffsAtEdge,2,2,idx,0,0,0),2))*cf2j+&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,3,2,idx,0,0,0))*wj)+&
              RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)-
                                   IDX3(DcoeffsAtEdge,3,2,idx,0,0,0),2))*cf3j,&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,0,0,0))*ui)+&
              RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)-
                                   IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),2))*cf1i+&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,2,1,idx,0,0,0))*vi)+&
              RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)-
                                   IDX3(DcoeffsAtEdge,2,1,idx,0,0,0),2))*cf2i+&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,3,1,idx,0,0,0))*wi)+&
              RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)-
                                   IDX3(DcoeffsAtEdge,3,1,idx,0,0,0),2))*cf3i)
#else
      ! Compute scalar dissipation with dimensional splitting
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*uj)+&
                  abs(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0))*cf1j,&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*ui)+&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*cf1i )&
           + max( abs(IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*vj)+&
                  abs(IDX3(DcoeffsAtEdge,2,1,idx,0,0,0))*cf2j,&
                  abs(IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*vi)+&
                  abs(IDX3(DcoeffsAtEdge,2,2,idx,0,0,0))*cf2i )&
           + max( abs(IDX3(DcoeffsAtEdge,3,1,idx,0,0,0)*wj)+&
                  abs(IDX3(DcoeffsAtEdge,3,1,idx,0,0,0))*cf3j,&
                  abs(IDX3(DcoeffsAtEdge,3,2,idx,0,0,0)*wi)+&
                  abs(IDX3(DcoeffsAtEdge,3,2,idx,0,0,0))*cf3i )
#endif

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))
      
      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
           IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fyj-&
           IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
           IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fyi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,1,idx,0,0,0)*Fy_ij + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,2,idx,0,0,0)*Fy_ij + Diff)
#endif
    end do

  end subroutine mhd_calcFluxRusDissDiSp3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiagMatD3d_sim(DdataAtNode, DcoeffsAtNode,&
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


  end subroutine mhd_calcMatDiagMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiag3d_sim(DdataAtNode,&
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


  end subroutine mhd_calcMatDiag3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGalMatD3d_sim(DdataAtEdge,&
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


  end subroutine mhd_calcMatGalMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGalerkin3d_sim(DdataAtEdge,&
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

  end subroutine mhd_calcMatGalerkin3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDissMatD3d_sim(DdataAtEdge,&
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


  end subroutine mhd_calcMatScDissMatD3d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDiss3d_sim(DdataAtEdge,&
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


  end subroutine mhd_calcMatScDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDissMatD3d_sim(DdataAtEdge,&
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


  end subroutine mhd_calcMatRoeDissMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDiss3d_sim(DdataAtEdge,&
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


  end subroutine mhd_calcMatRoeDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDissMatD3d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 3D and applies the scalar artificial viscosities of Rusanov-type.
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


  end subroutine mhd_calcMatRusDissMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDiss3d_sim(DdataAtEdge,&
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


  end subroutine mhd_calcMatRusDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcCharacteristics3d_sim(Dweight, DdataAtEdge,&
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


  end subroutine mhd_calcCharacteristics3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTScDiss3d_sim(DdataAtEdge, DcoeffsAtEdge,&
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


  end subroutine mhd_calcFluxFCTScDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTRoeDiss3d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

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

  end subroutine mhd_calcFluxFCTRoeDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTRusDiss3d_sim(DdataAtEdge, DcoeffsAtEdge,&
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

end module mhd_callback3d
