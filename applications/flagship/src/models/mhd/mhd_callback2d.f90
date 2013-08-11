!##############################################################################
!# ****************************************************************************
!# <name> mhd_callback2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible ideal MHD equations in 2D.
!#
!# The following callback functions are available:
!#
!# 1.) mhd_calcFluxGalerkin2d_sim
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
!# 12.) mhd_calcMatGalerkin2d_sim
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
!# 25.) mhd_trafoNodalDensity2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density values
!#
!# 26.) mhd_trafoFluxEnergy2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 27.) mhd_trafoDiffEnergy2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 28.) mhd_trafoNodalEnergy2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal energy values
!#
!# 29.) mhd_trafoFluxPressure2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 30.) mhd_trafoDiffPressure2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the presure
!#
!# 31.) mhd_trafoNodalPressure2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal pressure values
!#
!# 32.) mhd_trafoFluxVelocity2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 33.) mhd_trafoDiffVelocity2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 34.) mhd_trafoNodalVelocity2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal velocity values
!#
!# 35.) mhd_trafoFluxMomentum2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 36.) mhd_trafoDiffMomentum2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 37.) mhd_trafoNodalMomentum2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal momentum values
!#
!# 38.) mhd_trafoFluxDenEng2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 39.) mhd_trafoDiffDenEng2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 40.) mhd_trafoNodalDenEng2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density and energy values
!#
!# 41.) mhd_trafoFluxDenPre2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 42.) mhd_trafoDiffDenPre2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 43.) mhd_trafoNodalDenPre2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density and pressure values
!#
!# 44.) mhd_trafoFluxDenPreVel2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 45.) mhd_trafoDiffDenPreVel2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure
!#         and the velocity
!#
!# 46.) mhd_trafoNodalDenPreVel2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density, pressure and velocity values
!#
!# 47.) mhd_trafoDiffMagfield2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the magnetic field
!#
!# 48.) mhd_trafoFluxMagfield2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the magnetic field
!#
!# 49.) mhd_trafoNodalMagfield2d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal values for the magnetic field
!#
!# 50.) mhd_calcBoundaryvalues2d
!#      -> Computes the boundary values for a given node
!#
!# 51.) mhd_calcBilfBdrCond2d
!#      -> Calculates the bilinear form arising from the weak
!#         imposition of boundary conditions in 2D
!#
!# 52.) mhd_calcLinfBdrCond2d
!#      -> Calculates the linear form arising from the weak
!#         imposition of boundary conditions in 2D
!#
!# 53.) mhd_coeffVectorBdr2d_sim
!#      -> Calculates the coefficients for the linear form in 2D
!#
!# 54.) mhd_coeffMatrixBdr2d_sim
!#      -> Calculates the coefficients for the bilinear form in 2D
!# </purpose>
!##############################################################################

module mhd_callback2d

#include "flagship.h"
#include "mhd.h"
#include "kernel/System/fmath.h"

!$use omp_lib
  use basicgeometry
  use boundary
  use boundarycondaux
  use collection
  use derivatives
  use flagship_callback
  use fsystem
  use genoutput
  use graph
  use groupfemsystem
  use linearformevaluation
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use problem
  use scalarpde
  use solveraux
  use storage

  ! Modules from MHD model
  use mhd_basic

  implicit none

  private
  public :: mhd_calcFluxGalerkin2d_sim
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
  public :: mhd_calcMatGalerkin2d_sim
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
  public :: mhd_trafoNodalDensity2d_sim
  public :: mhd_trafoNodalEnergy2d_sim
  public :: mhd_trafoNodalPressure2d_sim
  public :: mhd_trafoNodalVelocity2d_sim
  public :: mhd_trafoNodalMomentum2d_sim
  public :: mhd_trafoNodalDenEng2d_sim
  public :: mhd_trafoNodalDenPre2d_sim
  public :: mhd_trafoNodalDenPreVel2d_sim
  public :: mhd_trafoNodalMagfield2d_sim
  public :: mhd_calcBoundaryvalues2d
  public :: mhd_calcBilfBdrCond2d
  public :: mhd_calcLinfBdrCond2d
  public :: mhd_coeffVectorBdr2d_sim

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGalerkin2d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)
      
#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Assemble skew-symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj-&
           IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi )
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
           
      ! Assemble fluxes
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij)
#endif
    end do

  end subroutine mhd_calcFluxGalerkin2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGalNoBdr2d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

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
    real(DP), dimension(NVAR2D) :: Fx_ij,Fy_ij
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)

      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Assemble symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
          (DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                        IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*Fx_ij+&
           DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                        IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*Fy_ij)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
    end do

  end subroutine mhd_calcFluxGalNoBdr2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxScDiss2d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)
    
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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
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

  end subroutine mhd_calcFluxScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxScDissDiSp2d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)
    

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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#endif
      
      !-------------------------------------------------------------------------
      ! Evaluate the scalar artificial dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = DCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      Diff = DCONST(0.0)

      !!! TODO !!!

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

  end subroutine mhd_calcFluxScDissDiSp2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRoeDiss2d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = DCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

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
      else
        
#ifdef MHD_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
             IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj-&
             IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
             IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi)
        IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
        IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
            (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij)
        IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij)
#endif
      end if
    end do

  end subroutine mhd_calcFluxRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRoeDissDiSp2d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)
    

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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = DCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      if (anorm .gt. SYS_EPSREAL_DP) then

        DiffX = DCONST(0.0); DiffY = DCONST(0.0)

        !!! TODO !!!

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef MHD_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
             IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj-&
             IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
             IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi+&
             DiffX + DiffY)
        IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
        IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
            (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij+&
             DiffX + DiffY)
        IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij+&
             DiffX + DiffY)
#endif
      else
        
#ifdef MHD_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
             IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj-&
             IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
             IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi)
        IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
        IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
            (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij)
        IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
             IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij)
#endif
      end if
    end do

  end subroutine mhd_calcFluxRoeDissDiSp2d_sim
 
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRusDiss2d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)
    
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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
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
      ca1i = abs(XMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_))
      ca1j = abs(XMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_))

      ! Compute the speed of the Alfven waves in y-direction
      ca2i = abs(YMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_))
      ca2j = abs(YMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_))

      ! Compute the speed of sound
      aPow2i = DCONST(MAGNETOHYDRODYN_GAMMA)*&
               PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)/&
               DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      aPow2j = DCONST(MAGNETOHYDRODYN_GAMMA)*&
               PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)/&
               DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      astPow2i = MAGFIELDMAGNITUDE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)/&
                 DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_) + aPow2i
      astPow2j = MAGFIELDMAGNITUDE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)/&
                 DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_) + aPow2j

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

#ifdef MHD_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*uj+&
                      DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*vj)+&
                 DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-MYNEWLINE \
                                      IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)*cf1j+&
                                  POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-MYNEWLINE \
                                      IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2)*cf2j),&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_))*ui+&
                      DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_))*vi)+&
                 DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-MYNEWLINE \
                                      IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)*cf1i+&
                                  POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-MYNEWLINE \
                                      IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2)*cf2i) )
#else
       ! Compute scalar dissipation
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*uj+&
                      IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*vj)+&
                 sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)*cf1j+&
                      POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2)*cf2j),&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*ui+&
                      IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*vi)+&
                 sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)*cf1i+&
                      POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2)*cf2i) )
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

  end subroutine mhd_calcFluxRusDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRusDissDiSp2d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wi = ZVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wj = ZVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Compute total pressures
      pi = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = TOTALPRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      qi = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi)
      qj = MAG_DOT_VEL3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj)

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fxi,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fxj,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)

      ! Compute fluxes for y-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)
      IDX1(Fyi,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fyj,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,5) = INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,6) = INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,7) = INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fx_ij,8) = INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,5) = INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX5_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,6) = INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX6_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,7) = INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX7_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
      IDX1(Fy_ij,8) = INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,vi,wi,pi,qi)-&
                      INVISCIDFLUX8_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,vj,wj,pj,qj)
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
      ca1i = abs(XMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_))
      ca1j = abs(XMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_))

      ! Compute the speed of the Alfven waves in y-direction
      ca2i = abs(YMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_))
      ca2j = abs(YMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_))

      ! Compute the speed of sound
      aPow2i = DCONST(MAGNETOHYDRODYN_GAMMA)*&
               PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)/&
               DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      aPow2j = DCONST(MAGNETOHYDRODYN_GAMMA)*&
               PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)/&
               DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute auxiliary quantities
      astPow2i = MAGFIELDMAGNITUDE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)/&
                 DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_) + aPow2i
      astPow2j = MAGFIELDMAGNITUDE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)/&
                 DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_) + aPow2j

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

#ifdef MHD_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*uj)+&
              DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-MYNEWLINE \
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2))*cf1j+&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*vj)+&
              DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-MYNEWLINE \
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2))*cf2j,&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_))*ui)+&
              DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-MYNEWLINE \
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2))*cf1i+&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_))*vi)+&
              DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-MYNEWLINE \
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2))*cf2i)
#else
      ! Compute scalar dissipation with dimensional splitting
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*uj)+&
                  abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_))*cf1j,&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*ui)+&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*cf1i )&
           + max( abs(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*vj)+&
                  abs(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_))*cf2j,&
                  abs(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*vi)+&
                  abs(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*cf2i )
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

  end subroutine mhd_calcFluxRusDissDiSp2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiagMatD2d_sim(DdataAtNode, DcoeffsAtNode,&
      InodeList, dscale, nnodes, DmatrixAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 2D.
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


  end subroutine mhd_calcMatDiagMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiag2d_sim(DdataAtNode,&
      DcoeffsAtNode, InodeList, dscale, nnodes,&
      DmatrixAtNode, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 2D.
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


  end subroutine mhd_calcMatDiag2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGalMatD2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 2D.
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
  

  end subroutine mhd_calcMatGalMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGalerkin2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 2D.
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

      
  end subroutine mhd_calcMatGalerkin2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDissMatD2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 2D and applies scalar artificial viscosities proportional to
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


  end subroutine mhd_calcMatScDissMatD2d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDiss2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 2D and applies
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

   
  end subroutine mhd_calcMatScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDissMatD2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 2D and applies
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

  
  end subroutine mhd_calcMatRoeDissMatD2d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDiss2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 2D and applies
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


  end subroutine mhd_calcMatRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDissMatD2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 2D and applies the scalar artificial viscosities of Rusanov-type.
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


  end subroutine mhd_calcMatRusDissMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDiss2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 2D applies
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
    

  end subroutine mhd_calcCharacteristics2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTScDiss2d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

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
  

  end subroutine mhd_calcFluxFCTScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTRoeDiss2d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

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


  end subroutine mhd_calcFluxFCTRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTRusDiss2d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

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
          DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)
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
          DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do

  end subroutine mhd_trafoDiffDensity2d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalDensity2d_sim(DdataAtNode,&
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          DENSITY2_2D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine mhd_trafoNodalDensity2d_sim

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
          TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_)
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
          TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do

  end subroutine mhd_trafoDiffEnergy2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalEnergy2d_sim(DdataAtNode,&
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          TOTALENERGY2_2D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine mhd_trafoNodalEnergy2d_sim

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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wi = ZVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wj = ZVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          ((MAGNETOHYDRODYN_GAMMA)-DCONST(1.0))*(DCONST(0.5)*&
          (ui*ui+vi*vi+wi*wi)*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         ui*XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         vi*YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         wi*ZMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            XMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)*&
                            XMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            YMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)*&
                            YMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            ZMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)*&
                            ZMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                          TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_))
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
          -((MAGNETOHYDRODYN_GAMMA)-DCONST(1.0))*(DCONST(0.5)*&
          (uj*uj+vj*vj+wj*wj)*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         uj*XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         vj*YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         wj*ZMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            XMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)*&
                            XMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            YMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)*&
                            YMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            ZMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)*&
                            ZMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                          TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_))
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
          PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do

  end subroutine mhd_trafoDiffPressure2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalPressure2d_sim(DdataAtNode,&
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          PRESSURE2_2D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine mhd_trafoNodalPressure2d_sim
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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wi = ZVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wj = ZVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed velocity fluxes in x-direction
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          (XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
          ui*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_))/&
             DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -(XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
          uj*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_))/&
             DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed velocity fluxes in y-direction
      IDX3(DtransformedFluxesAtEdge,2,1,idx,_,_,_) =&
          (YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
          vi*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_))/&
             DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,_,_,_) =&
         -(YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
          vj*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_))/&
             DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed velocity fluxes in z-direction
      IDX3(DtransformedFluxesAtEdge,3,1,idx,_,_,_) =&
          (ZMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
          wi*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_))/&
             DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      IDX3(DtransformedFluxesAtEdge,3,2,idx,_,_,_) =&
         -(ZMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
          wj*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_))/&
             DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
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
          XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed velocity difference in y-direction
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed velocity difference in z-direction
      IDX2(DtransformedDataAtEdge,3,idx,_,_) =&
          ZVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          ZVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do
    
  end subroutine mhd_trafoDiffVelocity2d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalVelocity2d_sim(DdataAtNode,&
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          XVELOCITY2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed y-velocity values
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          YVELOCITY2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed z-velocity values
      IDX2(DtransformedDataAtNode,3,idx,_,_) =&
          ZVELOCITY2_2D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine mhd_trafoNodalVelocity2d_sim

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
          XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)

      ! Transformed momentum fluxes in y-direction
      IDX3(DtransformedFluxesAtEdge,2,1,idx,_,_,_) =&
          YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,_,_,_) =&
         -YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)

      ! Transformed momentum fluxes in z-direction
      IDX3(DtransformedFluxesAtEdge,3,1,idx,_,_,_) =&
          ZMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,3,2,idx,_,_,_) =&
         -ZMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)
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
          XMOMENTUM3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          XMOMENTUM3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed momentum difference in y-direction
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          YMOMENTUM3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          YMOMENTUM3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed momentum difference in z-direction
      IDX2(DtransformedDataAtEdge,3,idx,_,_) =&
          ZMOMENTUM3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          ZMOMENTUM3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do
    
  end subroutine mhd_trafoDiffMomentum2d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalMomentum2d_sim(DdataAtNode,&
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          XMOMENTUM2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed y-momentum values
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          YMOMENTUM2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed z-momentum values
      IDX2(DtransformedDataAtNode,3,idx,_,_) =&
          ZMOMENTUM2_2D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine mhd_trafoNodalMomentum2d_sim

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
          DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)

      ! Transformed total energy fluxes
      IDX3(DtransformedFluxesAtEdge,2,1,idx,_,_,_) =&
          TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,_,_,_) =&
         -TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_)
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
          DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed total energy difference
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do

  end subroutine mhd_trafoDiffDenEng2d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalDenEng2d_sim(DdataAtNode,&
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          DENSITY2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed energy values
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          TOTALENERGY2_2D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine mhd_trafoNodalDenEng2d_sim

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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wi = ZVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wj = ZVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed density fluxes
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)

      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,2,1,idx,_,_,_) =&
          ((MAGNETOHYDRODYN_GAMMA)-DCONST(1.0))*(DCONST(0.5)*&
          (ui*ui+vi*vi+wi*wi)*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         ui*XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         vi*YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         wi*ZMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            XMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)*&
                            XMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            YMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)*&
                            YMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            ZMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)*&
                            ZMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                          TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_))
      IDX3(DtransformedFluxesAtEdge,2,2,idx,_,_,_) =&
          -((MAGNETOHYDRODYN_GAMMA)-DCONST(1.0))*(DCONST(0.5)*&
          (uj*uj+vj*vj+wj*wj)*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         uj*XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         vj*YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         wj*ZMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            XMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)*&
                            XMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            YMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)*&
                            YMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            ZMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)*&
                            ZMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                          TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_))
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
      IDX2(DtransformedDataAtEdge,1,idx,_,_) = &
          DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      
      ! Transformed pressure difference
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do
    
  end subroutine mhd_trafoDiffDenPre2d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalDenPre2d_sim(DdataAtNode,&
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          DENSITY2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed pressure values
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          PRESSURE2_2D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine mhd_trafoNodalDenPre2d_sim

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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      wi = ZVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      wj = ZVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Transformed density fluxes
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)

      ! Transformed velocity fluxes in x-direction
      IDX3(DtransformedFluxesAtEdge,2,1,idx,_,_,_) =&
          (XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
          ui*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_))/&
             DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,_,_,_) =&
         -(XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
          uj*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_))/&
             DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed velocity fluxes in y-direction
      IDX3(DtransformedFluxesAtEdge,3,1,idx,_,_,_) =&
          (YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
          vi*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_))/&
             DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      IDX3(DtransformedFluxesAtEdge,3,2,idx,_,_,_) =&
         -(YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
          vj*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_))/&
             DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed velocity fluxes in z-direction
      IDX3(DtransformedFluxesAtEdge,4,1,idx,_,_,_) =&
          (ZMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
          wi*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_))/&
             DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      IDX3(DtransformedFluxesAtEdge,4,2,idx,_,_,_) =&
         -(ZMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
          wj*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_))/&
             DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,5,1,idx,_,_,_) =&
          ((MAGNETOHYDRODYN_GAMMA)-DCONST(1.0))*(DCONST(0.5)*&
          (ui*ui+vi*vi+wi*wi)*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         ui*XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         vi*YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         wi*ZMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            XMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)*&
                            XMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            YMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)*&
                            YMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            ZMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)*&
                            ZMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                          TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_))
      IDX3(DtransformedFluxesAtEdge,5,2,idx,_,_,_) =&
          -((MAGNETOHYDRODYN_GAMMA)-DCONST(1.0))*(DCONST(0.5)*&
          (uj*uj+vj*vj+wj*wj)*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         uj*XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         vj*YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                         wj*ZMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            XMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)*&
                            XMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            YMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)*&
                            YMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                            ZMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)*&
                            ZMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
                          TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_))
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
          DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed velocity difference in x-direction
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed velocity difference in y-direction
      IDX2(DtransformedDataAtEdge,3,idx,_,_) =&
          YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed velocity difference in z-direction
      IDX2(DtransformedDataAtEdge,4,idx,_,_) =&
          ZVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          ZVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed pressure difference
      IDX2(DtransformedDataAtEdge,5,idx,_,_) =&
          PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do

  end subroutine mhd_trafoDiffDenPreVel2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalDenPreVel2d_sim(DdataAtNode,&
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          DENSITY2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed x-velocity values
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          XVELOCITY2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed y-velocity values
      IDX2(DtransformedDataAtNode,3,idx,_,_) =&
          YVELOCITY2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed z-velocity values
      IDX2(DtransformedDataAtNode,4,idx,_,_) =&
          ZVELOCITY2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed pressure values
      IDX2(DtransformedDataAtNode,5,idx,_,_) =&
          PRESSURE2_2D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine mhd_trafoNodalDenPreVel2d_sim

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

      ! Transformed magnetic field fluxes in x-direction
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          XMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -XMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)

      ! Transformed magnetic field fluxes in y-direction
      IDX3(DtransformedFluxesAtEdge,2,1,idx,_,_,_) =&
          YMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,_,_,_) =&
         -YMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)

      ! Transformed magnetic field fluxes in z-direction
      IDX3(DtransformedFluxesAtEdge,3,1,idx,_,_,_) =&
          ZMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,3,2,idx,_,_,_) =&
         -ZMAGFIELD2_2D(DfluxesAtEdge,IDX2,idx,_,_)
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

      ! Transformed magnetic field difference in x-direction
      IDX2(DtransformedDataAtEdge,1,idx,_,_) =&
          XMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          XMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed magnetic field difference in y-direction
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          YMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          YMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed magnetic field difference in z-direction
      IDX2(DtransformedDataAtEdge,3,idx,_,_) =&
          ZMAGFIELD3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          ZMAGFIELD3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
    end do

  end subroutine mhd_trafoDiffMagfield2d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalMagfield2d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the magnetic field in 2D.
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
      
      ! Transformed x-component of the magnetic field
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          XMAGFIELD2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed y-component of the magnetic field
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          YMAGFIELD2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed z-component of the magnetic field
      IDX2(DtransformedDataAtNode,3,idx,_,_) =&
          ZMAGFIELD2_2D(DdataAtNode,IDX2,idx,_,_)
    end do

  end subroutine mhd_trafoNodalMagfield2d_sim
  
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

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcBilfBdrCond2D(rproblemLevel, rboundaryCondition,&
      rsolution, dtime, dscale, ssectionName, fcoeff_buildMatrixScBdr2D_sim,&
      rmatrix, rcollection, cconstrType)

!<description>
    ! This subroutine computes the bilinear form arising from the weak
    ! imposition of boundary conditions in 2D.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! boundary condition
    type(t_boundaryCondition), intent(in) :: rboundaryCondition
    
    ! solution vector
    type(t_vectorBlock), intent(in), target :: rsolution

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! callback routine for nonconstant coefficient matrices.
    include '../../../../../kernel/DOFMaintenance/intf_coefficientMatrixScBdr2D.inc'

    ! OPTIONAL: One of the BILF_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! BILF_MATC_ELEMENTBASED is used.
    integer, intent(in), optional :: cconstrType
!</intput>

!<inputoutput>
    ! matrix
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! At the moment, nothing is done in this subroutine and it should
    ! not be called. It may be necessary to assemble some bilinear
    ! forms at the boundary in future.

  end subroutine mhd_calcBilfBdrCond2D

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcLinfBdrCond2D(rproblemLevel, rboundaryCondition,&
      rsolution, dtime, dscale, ssectionName, fcoeff_buildVectorBlBdr2D_sim,&
      rvector, rcollection)

!<description>
    ! This subroutine computes the linear form arising from the weak
    ! imposition of boundary conditions in 2D.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! boundary condition
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! solution vector
    type(t_vectorBlock), intent(in), target :: rsolution

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! callback routine for nonconstant coefficient vectors.
    include '../../../../../kernel/DOFMaintenance/intf_coefficientVectorBlBdr2D.inc'
!</intput>

!<inputoutput>
    ! residual/right-hand side vector
    type(t_vectorBlock), intent(inout) :: rvector

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_collection) :: rcollectionTmp
    type(t_boundaryRegion) :: rboundaryRegion,rboundaryRegionMirror,rregion
    type(t_linearForm) :: rform
    integer, dimension(:), pointer :: p_IbdrCondCpIdx, p_IbdrCondType
    integer, dimension(:), pointer :: p_IbdrCompPeriodic, p_IbdrCondPeriodic
    integer :: ibct, isegment
    integer(I32) :: ccubTypeBdr

    ! Evaluate linear form for boundary integral and return if
    ! there are no weak boundary conditions available
    if (.not.rboundaryCondition%bWeakBdrCond) return
    
    ! Check if we are in 2D
    if (rproblemLevel%rtriangulation%ndim .ne. NDIM2D) then
      call output_line('Spatial dimension must be 2D!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcLinfBdrCond2D')
      call sys_halt()
    end if

    ! Get pointer to parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
          'rparlist', ssectionName=ssectionName)
    
    ! Get parameters from parameter list
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'ccubTypeBdr', ccubTypeBdr)

    ! Initialise temporal collection structure
    call collct_init(rcollectionTmp)
    
    ! Prepare quick access arrays of temporal collection structure
    rcollectionTmp%SquickAccess(1) = ''
    rcollectionTmp%SquickAccess(2) = 'rfparser'
    rcollectionTmp%DquickAccess(1) = dtime
    rcollectionTmp%DquickAccess(2) = dscale
    rcollectionTmp%IquickAccess(4) = int(ccubTypeBdr)

    ! Attach user-defined collection structure to temporal collection
    ! structure (may be required by the callback function)
    rcollectionTmp%p_rnextCollection => rcollection

    ! Attach solution vector to first quick access vector of the
    ! temporal collection structure
    rcollectionTmp%p_rvectorQuickAccess1 => rsolution
    
    ! Attach function parser from boundary conditions to collection
    ! structure and specify its name in quick access string array
    call collct_setvalue_pars(rcollectionTmp, 'rfparser',&
        rboundaryCondition%rfparser, .true.)
    
    
    ! Set pointers
    call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx,&
        p_IbdrCondCpIdx)
    call storage_getbase_int(rboundaryCondition%h_IbdrCondType,&
        p_IbdrCondType)

    ! Set additional pointers for periodic boundary conditions
    if (rboundaryCondition%bPeriodic) then
      call storage_getbase_int(rboundaryCondition%h_IbdrCompPeriodic,&
          p_IbdrCompPeriodic)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondPeriodic,&
          p_IbdrCondPeriodic)
    end if
    
    ! Loop over all boundary components
    do ibct = 1, rboundaryCondition%iboundarycount
      
      ! Loop over all boundary segments
      do isegment = p_IbdrCondCpIdx(ibct), p_IbdrCondCpIdx(ibct+1)-1
        
        ! Check if this segment has weak boundary conditions
        if (iand(p_IbdrCondType(isegment), BDRC_WEAK) .ne. BDRC_WEAK) cycle
        
        ! Prepare further quick access arrays of temporal collection
        ! structure with boundary component, type and maximum expressions
        rcollectionTmp%IquickAccess(1) = p_IbdrCondType(isegment)
        rcollectionTmp%IquickAccess(2) = isegment
        rcollectionTmp%IquickAccess(3) = rboundaryCondition%nmaxExpressions
        
        ! Initialise the linear form
        rform%itermCount = 1
        rform%Idescriptors(1) = DER_FUNC
        
        ! Create boundary segment in 01-parametrisation
        call bdrc_createRegion(rboundaryCondition, ibct,&
            isegment-p_IbdrCondCpIdx(ibct)+1, rboundaryRegion)
        
        ! Check if special treatment of mirror boundary condition is required
        if ((iand(p_IbdrCondType(isegment), BDRC_TYPEMASK) .eq. BDRC_PERIODIC) .or.&
            (iand(p_IbdrCondType(isegment), BDRC_TYPEMASK) .eq. BDRC_ANTIPERIODIC)) then
          
          ! Create boundary region for mirror boundary in 01-parametrisation
          call bdrc_createRegion(rboundaryCondition, p_IbdrCompPeriodic(isegment),&
              p_IbdrCondPeriodic(isegment)-p_IbdrCondCpIdx(p_IbdrCompPeriodic(isegment))+1,&
              rboundaryRegionMirror)
          
          ! Attach boundary regin to temporal collection structure
          call collct_setvalue_bdreg(rcollectionTmp, 'rboundaryRegionMirror',&
              rboundaryRegionMirror, .true.)
          ! In the callback-function, the minimum/maximum parameter
          ! values of the boundary region and its mirrored
          ! counterpartqq are required in length parametrisation to
          ! determine the parameter values of the mirrored cubature
          ! points. Therefore, we make a copy of both boundary
          ! regions, convert them to length parametrisation and attach
          ! the minimum/maximum parameter values to the quick access
          ! arrays of the temporal collection structure.
          rregion = rboundaryRegion
          call boundary_convertRegion(&
              rvector%RvectorBlock(1)%p_rspatialDiscr%p_rboundary,&
              rregion, BDR_PAR_LENGTH)

          ! Prepare quick access array of temporal collection structure
          rcollectionTmp%DquickAccess(3) = rregion%dminParam
          rcollectionTmp%DquickAccess(4) = rregion%dmaxParam

          rregion = rboundaryRegionMirror
          call boundary_convertRegion(&
              rvector%RvectorBlock(1)%p_rspatialDiscr%p_rboundary,&
              rregion, BDR_PAR_LENGTH)
          
          ! Prepare quick access array of temporal collection structure
          rcollectionTmp%DquickAccess(5) = rregion%dminParam
          rcollectionTmp%DquickAccess(6) = rregion%dmaxParam
        end if

        ! Assemble the linear form
        if (rvector%nblocks .eq. 1) then
          call linf_buildVecIntlScalarBdr2d(rform, ccubTypeBdr, .false.,&
              rvector%RvectorBlock(1), fcoeff_buildVectorBlBdr2D_sim,&
              rboundaryRegion, rcollectionTmp)
        else
          call linf_buildVectorBlockBdr2d(rform, ccubTypeBdr, .false.,&
              rvector, fcoeff_buildVectorBlBdr2D_sim,&
              rboundaryRegion, rcollectionTmp)
        end if
        
      end do ! isegment
    end do ! ibct
        
    ! Release temporal collection structure
    call collct_done(rcollectionTmp)

  end subroutine mhd_calcLinfBdrCond2D

  !***************************************************************************

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
    ! OPTIONAL: A collection structure to provide additional
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

end module mhd_callback2d
