
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
!# 51.) mhd_hadaptCallbackScalar2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in interleave format
!#
!# 52.) mhd_hadaptCallbackBlock2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in block format
!#
!# 53.) mhd_coeffVectorBdr2d_sim
!#      -> Calculates the coefficients for the linear form in 2D
!#
!# 54.) mhd_coeffMatrixBdr2d_sim
!#      -> Calculates the coefficients for the bilinear form in 2D
!# </purpose>
!##############################################################################

module mhd_callback2d

#define MHD_NDIM 2
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
      
      ! Assemble skew-symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*Fxj+&
           IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*Fyj-&
           IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*Fxi-&
           IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*Fyi )
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
           
      ! Assemble fluxes
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale * IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*F_ij
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale * IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*F_ij
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
      
      ! Assemble symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (RCONST(0.5)*(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)-&
                        IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0))*Fx_ij+&
           RCONST(0.5)*(IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)-&
                        IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0))*Fy_ij)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
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
#endif
      
      !-------------------------------------------------------------------------
      ! Evaluate the scalar artificial dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = RCONST(0.5)*(IDX3(DmatrixCoeffsAtEdge,:,1,idx,0,0,0)-&
                       IDX3(DmatrixCoeffsAtEdge,:,2,idx,0,0,0))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      Diff = RCONST(0.0)

      !!! TODO !!!

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
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
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = RCONST(0.5)*(IDX3(DmatrixCoeffsAtEdge,:,1,idx,0,0,0)-&
                       IDX3(DmatrixCoeffsAtEdge,:,2,idx,0,0,0))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Normalize the skew-symmetric coefficient
        a = a/anorm
        
        Diff = RCONST(0.0)

        !! TODO !!

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef MHD_USE_IBP
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
        
#ifdef MHD_USE_IBP
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
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = RCONST(0.5)*(IDX3(DmatrixCoeffsAtEdge,:,1,idx,0,0,0)-&
                       IDX3(DmatrixCoeffsAtEdge,:,2,idx,0,0,0))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      if (anorm .gt. SYS_EPSREAL_DP) then

        DiffX = RCONST(0.0); DiffY = RCONST(0.0)

        !!! TODO !!!

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef MHD_USE_IBP
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
        
#ifdef MHD_USE_IBP
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

#ifdef MHD_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(RCONST(0.5)*(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)-&
                                   IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0))*uj+&
                      RCONST(0.5)*(IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)-&
                                   IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0))*vj)+&
                 RCONST(0.5)*sqrt(POW(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)-
                                      IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),2)*cf1j+&
                                  POW(IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)-
                                      IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),2)*cf2j),&
                  abs(RCONST(0.5)*(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)-&
                                   IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0))*ui+&
                      RCONST(0.5)*(IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)-&
                                   IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0))*vi)+&
                 RCONST(0.5)*sqrt(POW(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)-
                                      IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),2)*cf1i+&
                                  POW(IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)-
                                      IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),2)*cf2i) )
#else
       ! Compute scalar dissipation
      d_ij = max( abs(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*uj+&
                      IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*vj)+&
                 sqrt(POW(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),2)*cf1j+&
                      POW(IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),2)*cf2j),&
                  abs(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*ui+&
                      IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*vi)+&
                 sqrt(POW(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),2)*cf1i+&
                      POW(IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),2)*cf2i) )
#endif

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
                   IDX3(DdataAtEdge,:,1,idx,0,0,0))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------
      
#ifdef MHD_USE_IBP
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

#ifdef MHD_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(RCONST(0.5)*(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)-&
                                   IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0))*uj)+&
              RCONST(0.5)*sqrt(POW(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)-
                                   IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0),2))*cf1j+&
                  abs(RCONST(0.5)*(IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)-&
                                   IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0))*vj)+&
              RCONST(0.5)*sqrt(POW(IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)-
                                   IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0),2))*cf2j,&
                  abs(RCONST(0.5)*(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)-&
                                   IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0))*ui)+&
              RCONST(0.5)*sqrt(POW(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)-
                                   IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0),2))*cf1i+&
                  abs(RCONST(0.5)*(IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)-&
                                   IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0))*vi)+&
              RCONST(0.5)*sqrt(POW(IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)-
                                   IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0),2))*cf2i)
#else   
      ! Compute scalar dissipation with dimensional splitting
      d_ij = max( abs(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0)*uj)+&
                  abs(IDX3(DmatrixCoeffsAtEdge,1,1,idx,0,0,0))*cf1j,&
                  abs(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0)*ui)+&
                  abs(IDX3(DmatrixCoeffsAtEdge,1,2,idx,0,0,0))*cf1i )&
           + max( abs(IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0)*vj)+&
                  abs(IDX3(DmatrixCoeffsAtEdge,2,1,idx,0,0,0))*cf2j,&
                  abs(IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0)*vi)+&
                  abs(IDX3(DmatrixCoeffsAtEdge,2,2,idx,0,0,0))*cf2i )
#endif

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
                   IDX3(DdataAtEdge,:,1,idx,0,0,0))
      
      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
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
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
          DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
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
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          DENSITY2(DdataAtNode,IDX2,idx,0,0)
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
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0)
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
          TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
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
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          TOTALENERGY2(DdataAtNode,IDX2,idx,0,0)
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
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      wi = ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      wj = ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          ((MAGNETOHYDRODYN_GAMMA)-RCONST(1.0))*(RCONST(0.5)*&
          (ui*ui+vi*vi+wi*wi)*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         ui*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         vi*YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         wi*ZMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)*&
                            XMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)*&
                            YMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)*&
                            ZMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                          TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
          -((MAGNETOHYDRODYN_GAMMA)-RCONST(1.0))*(RCONST(0.5)*&
          (uj*uj+vj*vj+wj*wj)*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         uj*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         vj*YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         wj*ZMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0)*&
                            XMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0)*&
                            YMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0)*&
                            ZMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                          TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
          PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
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
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          PRESSURE2(DdataAtNode,IDX2,idx,0,0)
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
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      wi = ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      wj = ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

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

      ! Transformed velocity fluxes in z-direction
      IDX3(DtransformedFluxesAtEdge,3,1,idx,0,0,0) =&
          (ZMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
          wi*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0))/&
             DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      IDX3(DtransformedFluxesAtEdge,3,2,idx,0,0,0) =&
         -(ZMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
          wj*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0))/&
             DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
          XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      ! Transformed velocity difference in y-direction
      IDX2(DtransformedDataAtEdge,2,idx,0,0) =&
          YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      ! Transformed velocity difference in z-direction
      IDX2(DtransformedDataAtEdge,3,idx,0,0) =&
          ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
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
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          XVELOCITY2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed y-velocity values
      IDX2(DtransformedDataAtNode,2,idx,0,0) =&
          YVELOCITY2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed z-velocity values
      IDX2(DtransformedDataAtNode,3,idx,0,0) =&
          ZVELOCITY2(DdataAtNode,IDX2,idx,0,0)
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
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)

      ! Transformed momentum fluxes in y-direction
      IDX3(DtransformedFluxesAtEdge,2,1,idx,0,0,0) =&
          YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,0,0,0) =&
         -YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)

      ! Transformed momentum fluxes in z-direction
      IDX3(DtransformedFluxesAtEdge,3,1,idx,0,0,0) =&
          ZMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,3,2,idx,0,0,0) =&
         -ZMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
          XMOMENTUM3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          XMOMENTUM3(DdataAtEdge,IDX3,1,idx,0,0,0)

      ! Transformed momentum difference in y-direction
      IDX2(DtransformedDataAtEdge,2,idx,0,0) =&
          YMOMENTUM3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          YMOMENTUM3(DdataAtEdge,IDX3,1,idx,0,0,0)

      ! Transformed momentum difference in z-direction
      IDX2(DtransformedDataAtEdge,3,idx,0,0) =&
          ZMOMENTUM3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          ZMOMENTUM3(DdataAtEdge,IDX3,1,idx,0,0,0)
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
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          XMOMENTUM2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed y-momentum values
      IDX2(DtransformedDataAtNode,2,idx,0,0) =&
          YMOMENTUM2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed z-momentum values
      IDX2(DtransformedDataAtNode,3,idx,0,0) =&
          ZMOMENTUM2(DdataAtNode,IDX2,idx,0,0)
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
          DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      ! Transformed total energy difference
      IDX2(DtransformedDataAtEdge,2,idx,0,0) =&
          TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
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
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          DENSITY2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed energy values
      IDX2(DtransformedDataAtNode,2,idx,0,0) =&
          TOTALENERGY2(DdataAtNode,IDX2,idx,0,0)
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
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      wi = ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      wj = ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Transformed density fluxes
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)

      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,2,1,idx,0,0,0) =&
          ((MAGNETOHYDRODYN_GAMMA)-RCONST(1.0))*(RCONST(0.5)*&
          (ui*ui+vi*vi+wi*wi)*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         ui*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         vi*YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         wi*ZMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)*&
                            XMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)*&
                            YMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)*&
                            ZMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                          TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
      IDX3(DtransformedFluxesAtEdge,2,2,idx,0,0,0) =&
          -((MAGNETOHYDRODYN_GAMMA)-RCONST(1.0))*(RCONST(0.5)*&
          (uj*uj+vj*vj+wj*wj)*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         uj*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         vj*YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         wj*ZMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0)*&
                            XMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0)*&
                            YMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0)*&
                            ZMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                          TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) = &
          DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      
      ! Transformed pressure difference
      IDX2(DtransformedDataAtEdge,2,idx,0,0) =&
          PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
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
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          DENSITY2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed pressure values
      IDX2(DtransformedDataAtNode,2,idx,0,0) =&
          PRESSURE2(DdataAtNode,IDX2,idx,0,0)
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
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      wi = ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      wj = ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
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

      ! Transformed velocity fluxes in z-direction
      IDX3(DtransformedFluxesAtEdge,4,1,idx,0,0,0) =&
          (ZMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
          wi*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0))/&
             DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      IDX3(DtransformedFluxesAtEdge,4,2,idx,0,0,0) =&
         -(ZMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
          wj*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0))/&
             DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,5,1,idx,0,0,0) =&
          ((MAGNETOHYDRODYN_GAMMA)-RCONST(1.0))*(RCONST(0.5)*&
          (ui*ui+vi*vi+wi*wi)*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         ui*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         vi*YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         wi*ZMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)*&
                            XMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)*&
                            YMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)*&
                            ZMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                          TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
      IDX3(DtransformedFluxesAtEdge,5,2,idx,0,0,0) =&
          -((MAGNETOHYDRODYN_GAMMA)-RCONST(1.0))*(RCONST(0.5)*&
          (uj*uj+vj*vj+wj*wj)*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         uj*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         vj*YMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                         wj*ZMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0)*&
                            XMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0)*&
                            YMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                            ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0)*&
                            ZMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)-&
                          TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
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

      ! Transformed velocity difference in z-direction
      IDX2(DtransformedDataAtEdge,4,idx,0,0) =&
          ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)

      ! Transformed pressure difference
      IDX2(DtransformedDataAtEdge,5,idx,0,0) =&
          PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
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
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          DENSITY2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed x-velocity values
      IDX2(DtransformedDataAtNode,2,idx,0,0) =&
          XVELOCITY2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed y-velocity values
      IDX2(DtransformedDataAtNode,3,idx,0,0) =&
          YVELOCITY2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed z-velocity values
      IDX2(DtransformedDataAtNode,4,idx,0,0) =&
          ZVELOCITY2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed pressure values
      IDX2(DtransformedDataAtNode,5,idx,0,0) =&
          PRESSURE2(DdataAtNode,IDX2,idx,0,0)
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
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          XMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -XMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)

      ! Transformed magnetic field fluxes in y-direction
      IDX3(DtransformedFluxesAtEdge,2,1,idx,0,0,0) =&
          YMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,0,0,0) =&
         -YMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)

      ! Transformed magnetic field fluxes in z-direction
      IDX3(DtransformedFluxesAtEdge,3,1,idx,0,0,0) =&
          ZMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,3,2,idx,0,0,0) =&
         -ZMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)
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
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
          XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)

      ! Transformed magnetic field difference in y-direction
      IDX2(DtransformedDataAtEdge,2,idx,0,0) =&
          YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)

      ! Transformed magnetic field difference in z-direction
      IDX2(DtransformedDataAtEdge,3,idx,0,0) =&
          ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)
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
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          XMAGFIELD2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed y-component of the magnetic field
      IDX2(DtransformedDataAtNode,2,idx,0,0) =&
          YMAGFIELD2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed z-component of the magnetic field
      IDX2(DtransformedDataAtNode,3,idx,0,0) =&
          ZMAGFIELD2(DdataAtNode,IDX2,idx,0,0)
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


    case default
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


    case default
      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)

    end select

  end subroutine mhd_hadaptCallbackBlock2d

end module mhd_callback2d
