!##############################################################################
!# ****************************************************************************
!# <name> mhd_callback3d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible ideal MHD equations in 3D.
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

#include "flagship.h"
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
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Assemble skew-symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj+&
           IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*Fzj-&
           IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi-&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fzi )
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

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
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)

      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Assemble fluxes
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
          (DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                        IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*Fx_ij+&
           DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                        IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*Fy_ij+&
           DCONST(0.5)*(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)-&
                        IDX3(DcoeffsAtEdge,3,2,idx,_,_,_))*Fz_ij)
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
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)


#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar artificial dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------

      Diff = DCONST(0.0)

      !!! TODO !!!

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
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
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar artificial dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = DCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      Diff = DCONST(0.0)

      !!! TODO !!!

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
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
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = DCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Normalise the skew-symmetric coefficient
        a = a/anorm
        
        Diff = DCONST(0.0)

        !! TODO !!

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef MHD_USE_IBP
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
        
#ifdef MHD_USE_IBP
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
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#endif
      
      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = DCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Normalise the skew-symmetric coefficient
        a = a/anorm
        
        DiffX = DCONST(0.0); DiffY = DCONST(0.0); DiffZ = DCONST(0.0)

        !!! TODO !!!

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef MHD_USE_IBP
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
        
#ifdef MHD_USE_IBP
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
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
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
      ca1i = abs(XMAGFIELD3_3D(DdataAtEdge,IDX3,1,idx,_,_,_))
      ca1j = abs(XMAGFIELD3_3D(DdataAtEdge,IDX3,2,idx,_,_,_))

      ! Compute the speed of the Alfven waves in y-direction
      ca2i = abs(YMAGFIELD3_3D(DdataAtEdge,IDX3,1,idx,_,_,_))
      ca2j = abs(YMAGFIELD3_3D(DdataAtEdge,IDX3,2,idx,_,_,_))

      ! Compute the speed of the Alfven waves in z-direction
      ca3i = abs(ZMAGFIELD3_3D(DdataAtEdge,IDX3,1,idx,_,_,_))
      ca3j = abs(ZMAGFIELD3_3D(DdataAtEdge,IDX3,2,idx,_,_,_))

      ! Compute the speed of sound
      aPow2i = DCONST(MAGNETOHYDRODYN_GAMMA)*&
               PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)/&
               DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      aPow2j = DCONST(MAGNETOHYDRODYN_GAMMA)*&
               PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)/&
               DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      astPow2i = MAGFIELDMAGNITUDE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)/&
                 DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_) + aPow2i
      astPow2j = MAGFIELDMAGNITUDE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)/&
                 DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_) + aPow2j

      ! Compute the speed of the fast waves in x-direction
      cf1i = sqrt(DCONST(0.5)*(astPow2i+&
                  sqrt(POW(astPow2i,2)-DCONST(4.0)*aPow2i*POW(ca1i,2))))
      cf1j = sqrt(DCONST(0.5)*(astPow2j+&
                  sqrt(POW(astPow2j,2)-DCONST(4.0)*aPow2j*POW(ca1j,2))))

      ! Compute the speed of the fast waves in y-direction
      cf2i = sqrt(DCONST(0.5)*(astPow2i+&
                  sqrt(POW(astPow2i,2)-DCONST(4.0)*aPow2i*POW(ca2i,2))))
      cf2j = sqrt(DCONST(0.5)*(astPow2j+&
                  sqrt(POW(astPow2j,2)-DCONST(4.0)*aPow2j*POW(ca2j,2))))

      ! Compute the speed of the fast waves in z-direction
      cf3i = sqrt(DCONST(0.5)*(astPow2i+&
                  sqrt(POW(astPow2i,2)-DCONST(4.0)*aPow2i*POW(ca3i,2))))
      cf3j = sqrt(DCONST(0.5)*(astPow2j+&
                  sqrt(POW(astPow2j,2)-DCONST(4.0)*aPow2j*POW(ca3j,2))))

#ifdef MHD_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*uj+&
                      DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*vj+&
                      DCONST(0.5)*(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,3,2,idx,_,_,_))*wj)+&
                 DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-
                                      IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)*cf1j+&
                                  POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-
                                      IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2)*cf2j+&
                                  POW(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)-
                                      IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),2)*cf3j),&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_))*ui+&
                      DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_))*vi+&
                      DCONST(0.5)*(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,3,1,idx,_,_,_))*wi)+&
                 DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-
                                      IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)*cf1i+&
                                  POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-
                                      IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2)*cf2i+&
                                  POW(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)-
                                      IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),2)*cf3i) )
#else
       ! Compute scalar dissipation
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*uj+&
                      IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*vj+&
                      IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*wj)+&
                 sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)*cf1j+&
                      POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2)*cf2j+&
                      POW(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),2)*cf3j),&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*ui+&
                      IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*vi+&
                      IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*wi)+&
                 sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)*cf1i+&
                      POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2)*cf2i+&
                      POW(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),2)*cf3i) )
#endif

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,_,_,_)-&
                   IDX3(DdataAtEdge,:,1,idx,_,_,_))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------
      
#ifdef MHD_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj-&
           IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij + Diff)
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
      ui = XVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for z-direction
      IDX1(Fzi,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fzi,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fzj,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fzj,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute flux difference for z-direction
      IDX1(Fz_ij,1) = INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX1_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,2) = INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX2_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,3) = INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX3_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,4) = INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX4_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,5) = INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX5_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,6) = INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX6_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,7) = INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX7_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fz_ij,8) = INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                 INVISCIDFLUX8_ZDIR3_3D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
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
      ca1i = abs(XMAGFIELD3_3D(DdataAtEdge,IDX3,1,idx,_,_,_))
      ca1j = abs(XMAGFIELD3_3D(DdataAtEdge,IDX3,2,idx,_,_,_))

      ! Compute the speed of the Alfven waves in y-direction
      ca2i = abs(YMAGFIELD3_3D(DdataAtEdge,IDX3,1,idx,_,_,_))
      ca2j = abs(YMAGFIELD3_3D(DdataAtEdge,IDX3,2,idx,_,_,_))

      ! Compute the speed of the Alfven waves in z-direction
      ca3i = abs(ZMAGFIELD3_3D(DdataAtEdge,IDX3,1,idx,_,_,_))
      ca3j = abs(ZMAGFIELD3_3D(DdataAtEdge,IDX3,2,idx,_,_,_))

      ! Compute the speed of sound
      aPow2i = DCONST(MAGNETOHYDRODYN_GAMMA)*&
               PRESSURE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)/&
               DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)
      aPow2j = DCONST(MAGNETOHYDRODYN_GAMMA)*&
               PRESSURE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)/&
               DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      astPow2i = MAGFIELDMAGNITUDE3_3D(DdataAtEdge,IDX3,1,idx,_,_,_)/&
                 DENSITY3_3D(DdataAtEdge,IDX3,1,idx,_,_,_) + aPow2i
      astPow2j = MAGFIELDMAGNITUDE3_3D(DdataAtEdge,IDX3,2,idx,_,_,_)/&
                 DENSITY3_3D(DdataAtEdge,IDX3,2,idx,_,_,_) + aPow2j

      ! Compute the speed of the fast waves in x-direction
      cf1i = sqrt(DCONST(0.5)*(astPow2i+&
                  sqrt(POW(astPow2i,2)-DCONST(4.0)*aPow2i*POW(ca1i,2))))
      cf1j = sqrt(DCONST(0.5)*(astPow2j+&
                  sqrt(POW(astPow2j,2)-DCONST(4.0)*aPow2j*POW(ca1j,2))))

      ! Compute the speed of the fast waves in y-direction
      cf2i = sqrt(DCONST(0.5)*(astPow2i+&
                  sqrt(POW(astPow2i,2)-DCONST(4.0)*aPow2i*POW(ca2i,2))))
      cf2j = sqrt(DCONST(0.5)*(astPow2j+&
                  sqrt(POW(astPow2j,2)-DCONST(4.0)*aPow2j*POW(ca2j,2))))

      ! Compute the speed of the fast waves in z-direction
      cf3i = sqrt(DCONST(0.5)*(astPow2i+&
                  sqrt(POW(astPow2i,2)-DCONST(4.0)*aPow2i*POW(ca3i,2))))
      cf3j = sqrt(DCONST(0.5)*(astPow2j+&
                  sqrt(POW(astPow2j,2)-DCONST(4.0)*aPow2j*POW(ca3j,2))))
      
#ifdef MHD_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*uj)+&
              DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2))*cf1j+&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*vj)+&
              DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2))*cf2j+&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,3,2,idx,_,_,_))*wj)+&
              DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,3,2,idx,_,_,_),2))*cf3j,&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_))*ui)+&
              DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2))*cf1i+&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_))*vi)+&
              DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2))*cf2i+&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,3,1,idx,_,_,_))*wi)+&
              DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)-
                                   IDX3(DcoeffsAtEdge,3,1,idx,_,_,_),2))*cf3i)
#else
      ! Compute scalar dissipation with dimensional splitting
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*uj)+&
                  abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_))*cf1j,&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*ui)+&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*cf1i )&
           + max( abs(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*vj)+&
                  abs(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_))*cf2j,&
                  abs(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*vi)+&
                  abs(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*cf2i )&
           + max( abs(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_)*wj)+&
                  abs(IDX3(DcoeffsAtEdge,3,1,idx,_,_,_))*cf3j,&
                  abs(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_)*wi)+&
                  abs(IDX3(DcoeffsAtEdge,3,2,idx,_,_,_))*cf3i )
#endif

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))
      
      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj-&
           IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij + Diff)
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
