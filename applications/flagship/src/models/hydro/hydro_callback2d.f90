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
!# 1.) hydro_calcFluxGalerkin2d_sim
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
!# 12.) hydro_calcMatGalerkin2d_sim
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
!# 48.) hydro_calcBilfBdrCond2d
!#      -> Calculates the bilinear form arising from the weak
!#         imposition of boundary conditions in 2D
!#
!# 49.) hydro_calcLinfBdrCond2d
!#      -> Calculates the linear form arising from the weak
!#         imposition of boundary conditions in 2D
!#
!# 50.) hydro_coeffVectorBdr2d_sim
!#      -> Calculates the coefficients for the linear form in 2D
!#
!# </purpose>
!##############################################################################

module hydro_callback2d

#include "flagship.h"
#include "hydro.h"
#include "kernel/System/fmath.h"
#include "kernel/feat2constants.h"

!$use omp_lib
  use basicgeometry
  use boundary
  use boundaryaux
  use boundarycondaux
  use collection
  use cubature
  use derivatives
  use domainintegration
  use element
  use feevaluation
  use fparser
  use fsystem
  use genoutput
  use graph
  use groupfemsystem
  use linearformevaluation
  use linearsystemblock
  use linearsystemscalar
  use mprimitives
  use paramlist
  use problem
  use scalarpde
  use solveraux
  use spatialdiscretisation
  use storage

  ! Modules from hydrodynamic model
  use hydro_basic

  implicit none

  private

  public :: hydro_calcFluxGalerkin2d_sim
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
  public :: hydro_calcMatGalerkin2d_sim
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
  public :: hydro_calcBilfBdrCond2d
  public :: hydro_calcLinfBdrCond2d
  public :: hydro_coeffVectorBdr2d_sim
  
  ! CUDA wrapper routines
  public :: hydro_calcDivVecGalerkin2d_cuda
  public :: hydro_calcDivVecScDiss2d_cuda
  public :: hydro_calcDivVecScDissDiSp2d_cuda
  public :: hydro_calcDivVecRoeDiss2d_cuda
  public :: hydro_calcDivVecRoeDissDiSp2d_cuda
  public :: hydro_calcDivVecRusDiss2d_cuda
  public :: hydro_calcDivVecRusDissDiSp2d_cuda

  public :: hydro_calcDivMatGalMatD2d_cuda
  public :: hydro_calcDivMatGalerkin2d_cuda
  public :: hydro_calcDivMatScDissMatD2d_cuda
  public :: hydro_calcDivMatScDiss2d_cuda
  public :: hydro_calcDivMatRoeDissMatD2d_cuda
  public :: hydro_calcDivMatRoeDiss2d_cuda
  public :: hydro_calcDivMatRusDissMatD2d_cuda
  public :: hydro_calcDivMatRusDiss2d_cuda

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGalerkin2d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)

      ! Compute fluxes for x-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)

      ! Assemble skew-symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fxj+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fyj-&
           IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fxi-&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fyi )
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      
      ! Assemble fluxes
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*Fy_ij)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*Fx_ij+&
           IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*Fy_ij)
#endif
    end do

  end subroutine hydro_calcFluxGalerkin2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGalNoBdr2d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
    real(DP) :: pi,pj,ui,uj,vi,vj
    integer :: idx

    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------
      
      ! Compute velocities
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      
      ! Assemble symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,_,_,_) = dscale *&
          (DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                        IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*Fx_ij+&
           DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                        IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*Fy_ij)
      IDX3(DfluxesAtEdge,:,2,idx,_,_,_) = IDX3(DfluxesAtEdge,:,1,idx,_,_,_)
    end do

  end subroutine hydro_calcFluxGalNoBdr2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxScDiss2d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)

      ! Compute fluxes for x-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient
      a = DCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      ! Compute densities
      ri = DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      rj = DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute enthalpies
      hi = (TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
      hj = (TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(ri,rj)
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      H_ij = ROE_MEAN_VALUE(hi,hj,aux)
      
      ! Compute auxiliary variables
      vel_ij = u_ij*a(1) + v_ij*a(2)
      q_ij   = DCONST(0.5)*(u_ij*u_ij+v_ij*v_ij)

      ! Compute the speed of sound
      c_ij = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))
      
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

  end subroutine hydro_calcFluxScDiss2d_sim

  !***************************************************************************

!<subroutine>

  subroutine hydro_calcDivVecGalerkin2d_cuda(rgroupFEMSet, rx, ry, dscale,&
      bclear, fcb_calcFluxSys_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine computes the fluxes for the standard Galerkin
    ! discretisation in 2D.
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
      call hydro_calcFluxGalerkin2d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
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
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivVecGalerkin2d_cuda')
    call sys_halt()
    
#endif

  end subroutine hydro_calcDivVecGalerkin2d_cuda

  !***************************************************************************

!<subroutine>

  subroutine hydro_calcDivVecScDiss2d_cuda(rgroupFEMSet, rx, ry, dscale,&
      bclear, fcb_calcFluxSys_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 2D using scalar artificial viscosities proportional to the
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
      call hydro_calcFluxScDiss2d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
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
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivVecScDiss2d_cuda')
    call sys_halt()
    
#endif

  end subroutine hydro_calcDivVecScDiss2d_cuda 

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxScDissDiSp2d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)

      ! Compute fluxes for x-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient
      a = DCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))

      ! Compute densities
      ri = DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      rj = DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute enthalpies
      hi = (TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
      hj = (TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(ri,rj)
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      H_ij = ROE_MEAN_VALUE(hi,hj,aux)
      
      ! Compute auxiliary variable
      q_ij = DCONST(0.5)*(u_ij*u_ij+v_ij*v_ij)

      ! Compute the speed of sound
      c_ij = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))
      
      ! Compute scalar dissipation with dimensional splitting
      d_ij = ( abs(a(1)*u_ij) + abs(a(1))*c_ij +&
               abs(a(2)*v_ij) + abs(a(2))*c_ij )

      ! Multiply the solution difference by the scalar dissipation
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,_,_,_)-&
                   IDX3(DdataAtEdge,:,1,idx,_,_,_))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
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

  end subroutine hydro_calcFluxScDissDiSp2d_sim

  !***************************************************************************

!<subroutine>

  subroutine hydro_calcDivVecScDissDiSp2d_cuda(rgroupFEMSet, rx, ry, dscale,&
      bclear, fcb_calcFluxSys_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 2D using scalar artificial viscosities proportional to the
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
      call hydro_calcFluxScDissDiSp2d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
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
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivVecScDissDiSp2d_cuda')
    call sys_halt()
    
#endif

  end subroutine hydro_calcDivVecScDissDiSp2d_cuda

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRoeDiss2d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
#if defined(HYDRO_USE_ENTROPYFIX) && (HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX)
    real(DP) :: dtol,ci,cj,veli,velj
#endif

    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)

      ! Compute fluxes for x-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
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
        
        ! Compute densities
        ri = DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute enthalpies
        hi = (TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2)
        q_ij   = DCONST(0.5)*(u_ij*u_ij+v_ij*v_ij)

        ! Compute the speed of sound
        c2_ij = max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij  = sqrt(c2_ij)
        
        ! Compute eigenvalues
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)

#if defined(HYDRO_USE_ENTROPYFIX)

#if HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX

        ! Entropy-fix by Harten and Hyman
        ci = SOUNDSPEED3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        cj = SOUNDSPEED3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        veli = ui*a(1) + vi*a(2)
        velj = uj*a(1) + vj*a(2)

        dtol = max(DCONST(0.0), (vel_ij-c_ij) - (veli-ci), (velj-cj) - (vel_ij-c_ij) )
        
        if (l1 .lt. dtol)&
            l1 = DCONST(0.5)*((l1*l1)/dtol + dtol)
        
        dtol = max(DCONST(0.0), (vel_ij+c_ij) - (veli+ci), (velj+cj) - (vel_ij+c_ij) )

        if (l3 .lt. dtol)&
            l3 = DCONST(0.5)*((l3*l3)/dtol + dtol)

#elif HYDRO_USE_ENTROPYFIX == HARTEN_ENTROPYFIX

#ifndef HYDRO_HARTEN_ENTROPYFIX
#error "Value HYDRO_HARTEN_ENTROPYFIX is required!"
#else
        ! Entropy-fix by Harten
        if (l1 .lt. DCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l1 = DCONST(0.5)*((l1*l1)/DCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + DCONST(HYDRO_HARTEN_ENTROPYFIX))

        if (l3 .lt. DCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l3 = DCONST(0.5)*((l3*l3)/DCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + DCONST(HYDRO_HARTEN_ENTROPYFIX))
#endif
#else
#error "Invalid type of entropy fix!"
#endif
#endif
        
        ! Compute solution difference U_j-U_i
        Diff = IDX3(DdataAtEdge,:,2,idx,_,_,_)-&
               IDX3(DdataAtEdge,:,1,idx,_,_,_)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-DCONST(1.0))*(q_ij*Diff(1)&
                                           -u_ij*Diff(2)&
                                           -v_ij*Diff(3)&
                                                +Diff(4))/DCONST(2.0)/c2_ij
        aux2 = (vel_ij*Diff(1)&
                 -a(1)*Diff(2)&
                 -a(2)*Diff(3))/DCONST(2.0)/c_ij
        
        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((DCONST(1.0)-((HYDRO_GAMMA)-DCONST(1.0))*q_ij/c2_ij)*Diff(1)&
                                     +((HYDRO_GAMMA)-DCONST(1.0))*(u_ij*Diff(2)&
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
        
#ifdef HYDRO_USE_IBP
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

  end subroutine hydro_calcFluxRoeDiss2d_sim

  !***************************************************************************

!<subroutine>

  subroutine hydro_calcDivVecRoeDiss2d_cuda(rgroupFEMSet, rx, ry, dscale, bclear,&
      fcb_calcFluxSys_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 2D using scalar artificial viscosities of Roe-type.
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
      call hydro_calcFluxRoeDiss2d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
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
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivVecRoeDiss2d_cuda')
    call sys_halt()
    
#endif

  end subroutine hydro_calcDivVecRoeDiss2d_cuda

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRoeDissDiSp2d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
#if defined(HYDRO_USE_ENTROPYFIX) && (HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX)
    real(DP) :: dtol,ci,cj
#endif


    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)

      ! Compute fluxes for x-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
#endif
      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = DCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      if (anorm .gt. SYS_EPSREAL_DP) then

        ! Compute the absolute value
        a = abs(a)
        
        ! Compute densities
        ri = DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute enthalpies
        hi = (TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variable
        q_ij = DCONST(0.5)*(u_ij*u_ij+v_ij*v_ij)

        ! Compute the speed of sound
        c2_ij = max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij  = sqrt(c2_ij)

        !-----------------------------------------------------------------------
        ! Dimensional splitting: x-direction
        !-----------------------------------------------------------------------
        
        ! Compute eigenvalues
        l1 = abs(u_ij-c_ij)
        l2 = abs(u_ij)
        l3 = abs(u_ij+c_ij)
        l4 = abs(u_ij)

#if defined(HYDRO_USE_ENTROPYFIX)

#if HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX

        ! Entropy-fix by Harten and Hyman
        ci = SOUNDSPEED3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        cj = SOUNDSPEED3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

        dtol = max(DCONST(0.0), (u_ij-c_ij) - (ui-ci), (uj-cj) - (u_ij-c_ij) )
        
        if (l1 .lt. dtol)&
            l1 = DCONST(0.5)*((l1*l1)/dtol + dtol)
        
        dtol = max(DCONST(0.0), (u_ij+c_ij) - (ui+ci), (uj+cj) - (u_ij+c_ij) )

        if (l3 .lt. dtol)&
            l3 = DCONST(0.5)*((l3*l3)/dtol + dtol)

#elif HYDRO_USE_ENTROPYFIX == HARTEN_ENTROPYFIX

#ifndef HYDRO_HARTEN_ENTROPYFIX
#error "Value HYDRO_HARTEN_ENTROPYFIX is required!"
#else
        ! Entropy-fix by Harten
        if (l1 .lt. DCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l1 = DCONST(0.5)*((l1*l1)/DCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + DCONST(HYDRO_HARTEN_ENTROPYFIX))

        if (l3 .lt. DCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l3 = DCONST(0.5)*((l3*l3)/DCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + DCONST(HYDRO_HARTEN_ENTROPYFIX))
#endif
#else
#error "Invalid type of entropy fix!"
#endif
#endif

        ! Compute solution difference U_j-U_i
        DiffX = IDX3(DdataAtEdge,:,2,idx,_,_,_)-&
                IDX3(DdataAtEdge,:,1,idx,_,_,_)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-DCONST(1.0))*(q_ij*DiffX(1)&
                                      -u_ij*DiffX(2)&
                                      -v_ij*DiffX(3)&
                                           +DiffX(4))/DCONST(2.0)/c2_ij
        aux2 = (u_ij*DiffX(1)&
                    -DiffX(2))/DCONST(2.0)/c_ij
        
        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((DCONST(1.0)-((HYDRO_GAMMA)-DCONST(1.0))*q_ij/c2_ij)*DiffX(1)&
                                     +((HYDRO_GAMMA)-DCONST(1.0))*(u_ij*DiffX(2)&
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
        
#if defined(HYDRO_USE_ENTROPYFIX)

#if HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX

        ! Entropy-fix by Harten and Hyman
        dtol = max(DCONST(0.0), (v_ij-c_ij) - (vi-ci), (vj-cj) - (v_ij-c_ij) )
        
        if (l1 .lt. dtol)&
            l1 = DCONST(0.5)*((l1*l1)/dtol + dtol)
        
        dtol = max(DCONST(0.0), (v_ij+c_ij) - (vi+ci), (vj+cj) - (v_ij+c_ij) )

        if (l3 .lt. dtol)&
            l3 = DCONST(0.5)*((l3*l3)/dtol + dtol)

#elif HYDRO_USE_ENTROPYFIX == HARTEN_ENTROPYFIX

#ifndef HYDRO_HARTEN_ENTROPYFIX
#error "Value HYDRO_HARTEN_ENTROPYFIX is required!"
#else
        ! Entropy-fix by Harten
        if (l1 .lt. DCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l1 = DCONST(0.5)*((l1*l1)/DCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + DCONST(HYDRO_HARTEN_ENTROPYFIX))

        if (l3 .lt. DCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l3 = DCONST(0.5)*((l3*l3)/DCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + DCONST(HYDRO_HARTEN_ENTROPYFIX))
#endif
#else
#error "Invalid type of entropy fix!"
#endif
#endif

        ! Compute solution difference U_j-U_i
        DiffY = IDX3(DdataAtEdge,:,2,idx,_,_,_)-&
                IDX3(DdataAtEdge,:,1,idx,_,_,_)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-DCONST(1.0))*(q_ij*DiffY(1)&
                                      -u_ij*DiffY(2)&
                                      -v_ij*DiffY(3)&
                                           +DiffY(4))/DCONST(2.0)/c2_ij
        aux2 = (v_ij*DiffY(1)&
                    -DiffY(3))/DCONST(2.0)/c_ij

        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((DCONST(1.0)-((HYDRO_GAMMA)-DCONST(1.0))*q_ij/c2_ij)*DiffY(1)&
                                     +((HYDRO_GAMMA)-DCONST(1.0))*(u_ij*DiffY(2)&
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
        
#ifdef HYDRO_USE_IBP
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

  end subroutine hydro_calcFluxRoeDissDiSp2d_sim
 
  !***************************************************************************

!<subroutine>

  subroutine hydro_calcDivVecRoeDissDiSp2d_cuda(rgroupFEMSet, rx, ry, dscale,&
      bclear, fcb_calcFluxSys_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 2D using scalar artificial viscosities of Roe-type, whereby
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
      call lsysbl_copyH2D_Vector(rx, .false., .false.,  istream)
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
      call hydro_calcFluxRoeDissDiSp2d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
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
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivVecRoeDissDiSp2d_cuda')
    call sys_halt()
    
#endif

  end subroutine hydro_calcDivVecRoeDissDiSp2d_cuda

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRusDiss2d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)

      ! Compute fluxes for x-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !---------------------------------------------------------------------------
      
      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(HYDRO_GAMMA)*&
          (Ei-DCONST(0.5)*(ui*ui+vi*vi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(HYDRO_GAMMA)*&
          (Ej-DCONST(0.5)*(uj*uj+vj*vj)), SYS_EPSREAL_DP))
      
#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*uj+&
                      DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*vj)+&
                 DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-MYNEWLINE \
                                      IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)+&
                                  POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-MYNEWLINE \
                                      IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2))*cj,&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_))*ui+&
                      DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_))*vi)+&
                 DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-MYNEWLINE \
                                      IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)+&
                                  POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-MYNEWLINE \
                                      IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2))*ci )
#else
      ! Compute scalar dissipation
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*uj+&
                      IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*vj)+&
                 sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)+&
                      POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2))*cj,&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*ui+&
                      IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*vi)+&
                 sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)+&
                      POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2))*ci )
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

  end subroutine hydro_calcFluxRusDiss2d_sim

  !***************************************************************************

!<subroutine>

  subroutine hydro_calcDivVecRusDiss2d_cuda(rgroupFEMSet, rx, ry, dscale, bclear,&
      fcb_calcFluxSys_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 2D using scalar artificial viscosities of Rusanov-type.
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
      call hydro_calcFluxRusDiss2d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
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
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivVecRusDiss2d_cuda')
    call sys_halt()
    
#endif

  end subroutine hydro_calcDivVecRusDiss2d_cuda

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRusDissDiSp2d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      IDX1(Fxi,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)
      IDX1(Fxi,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)

      IDX1(Fxj,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fxj,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)

      ! Compute fluxes for x-direction
      IDX1(Fyi,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)
      IDX1(Fyi,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)

      IDX1(Fyj,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fyj,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
#else
      ! Compute flux difference for x-direction
      IDX1(Fx_ij,1) = INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX1_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,2) = INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX2_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,3) = INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX3_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
      IDX1(Fx_ij,4) = INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,ui,pi)-&
                      INVISCIDFLUX4_XDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,uj,pj)
                        
      ! Compute flux difference for y-direction
      IDX1(Fy_ij,1) = INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX1_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,2) = INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX2_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,3) = INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX3_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
      IDX1(Fy_ij,4) = INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,1,idx,_,_,_,vi,pi)-&
                      INVISCIDFLUX4_YDIR3_2D(DdataAtEdge,IDX3,2,idx,_,_,_,vj,pj)
#endif
      
      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !-------------------------------------------------------------------------

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(HYDRO_GAMMA)*&
          (Ei-DCONST(0.5)*(ui*ui+vi*vi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(HYDRO_GAMMA)*&
          (Ej-DCONST(0.5)*(uj*uj+vj*vj)), SYS_EPSREAL_DP))

#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation with dimensional splitting based on
      ! the skew-symmetric part which does not include the symmetric
      ! boundary contribution
      d_ij = max( abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*uj)+&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)))*cj,&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_))*ui)+&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)))*ci )&
           + max( abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*vj)+&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)))*cj,&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_))*vi)+&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)))*ci )
#else
      ! Compute scalar dissipation with dimensional splitting
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*uj)+&
                  abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_))*cj,&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*ui)+&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*ci )&
           + max( abs(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*vj)+&
                  abs(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_))*cj,&
                  abs(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*vi)+&
                  abs(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*ci )
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

  end subroutine hydro_calcFluxRusDissDiSp2d_sim

  !***************************************************************************

!<subroutine>

  subroutine hydro_calcDivVecRusDissDiSp2d_cuda(rgroupFEMSet, rx, ry, dscale,&
      bclear, fcb_calcFluxSys_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 2D using scalar artificial viscosities of Rusanov-type, whereby
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
      call hydro_calcFluxRusDissDiSp2d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
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
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivVecRusDissDiSp2d_cuda')
    call sys_halt()
    
#endif

  end subroutine hydro_calcDivVecRusDissDiSp2d_cuda

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiagMatD2d_sim(DdataAtNode, DcoeffsAtNode,&
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

    ! local variable
    real(DP) :: ui,vi
    integer :: inode


    do inode = 1, nnodes
      
      ! Compute auxiliary variables
      ui = XVELOCITY2_2D(DdataAtNode,IDX2,inode,_,_)
      vi = YVELOCITY2_2D(DdataAtNode,IDX2,inode,_,_)
      
#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ii = diag(A_i)*C_{ii}$
      IDX3(DmatrixAtNode,1,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,_)
      IDX3(DmatrixAtNode,2,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,_)
      IDX3(DmatrixAtNode,3,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,_)
      IDX3(DmatrixAtNode,4,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,_)
#else
      ! Compute Galerkin coefficient $K_ii = -diag(A_i)*C_{ii}$
      IDX3(DmatrixAtNode,1,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,_)
      IDX3(DmatrixAtNode,2,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,_)
      IDX3(DmatrixAtNode,3,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,_)
      IDX3(DmatrixAtNode,4,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,_)
#endif
    end do

  end subroutine hydro_calcMatDiagMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiag2d_sim(DdataAtNode,&
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

    ! local variable
    real(DP) :: Ei,ui,vi
    integer :: inode


    do inode = 1, nnodes
      
      ! Compute auxiliary variables
      ui = XVELOCITY2_2D(DdataAtNode,IDX2,inode,_,_)
      vi = YVELOCITY2_2D(DdataAtNode,IDX2,inode,_,_)
      Ei = SPECIFICTOTALENERGY2_2D(DdataAtNode,IDX2,inode,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      IDX3(DmatrixAtNode,1,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,2,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,3,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,4,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,5,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,6,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,7,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,8,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,9,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,10,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,11,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,12,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,13,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,14,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,15,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,16,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                               IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
#else
      ! Compute Galerkin coefficient $K_ii = -A_i*C_{ii}$
      IDX3(DmatrixAtNode,1,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,2,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,3,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,4,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,5,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,6,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,7,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,8,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,9,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,10,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,11,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,12,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,13,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,14,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,15,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
      IDX3(DmatrixAtNode,16,1,inode,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX2(DcoeffsAtNode,1,inode,_,_),MYNEWLINE \
                                                IDX2(DcoeffsAtNode,2,inode,_,_),ui,vi,Ei)
#endif
    end do

  end subroutine hydro_calcMatDiag2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGalMatD2d_sim(DdataAtEdge,&
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
  real(DP), dimension(:,:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>
  
    ! local variable
    real(DP) :: ui,uj,vi,vj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DmatrixAtEdge,1,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,2,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,3,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,4,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DmatrixAtEdge,1,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,2,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,3,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,4,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,_)
#endif
    end do

  end subroutine hydro_calcMatGalMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGalerkin2d_sim(DdataAtEdge,&
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

    ! local variable
    real(DP) :: Ei,Ej,ui,uj,vi,vj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DmatrixAtEdge,1,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,2,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,3,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,4,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,5,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,6,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,7,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,8,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,9,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,10,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,11,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,12,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,13,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,14,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,15,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,16,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DmatrixAtEdge,1,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,2,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,3,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,4,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,5,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,6,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,7,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,8,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,9,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,10,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,11,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,12,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,13,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,14,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,15,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,16,1,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
#endif
    end do
      
  end subroutine hydro_calcMatGalerkin2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDissMatD2d_sim(DdataAtEdge,&
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

    ! local variable
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj
    real(DP) :: H_ij,anorm,aux,c_ij,q_ij,u_ij,v_ij,vel_ij
    integer :: idx

    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vii,_)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,_)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = DCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Compute densities
        ri = DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute pressures
        pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute enthalpies
        hi = (TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
               
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2)
        q_ij   = DCONST(0.5)*(u_ij*u_ij+v_ij*v_ij)

        ! Compute the speed of sound
        c_ij = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))

        ! Compute scalar dissipation
        IDX3(DmatrixAtEdge,:,1,idx,_,_,_) = dscale * (abs(vel_ij) + anorm*c_ij)

      else
        
        ! Nullify dissipation tensor
        IDX3(DmatrixAtEdge,:,1,idx,_,_,_) = DCONST(0.0)

      end if
    end do

  end subroutine hydro_calcMatScDissMatD2d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDiss2d_sim(DdataAtEdge,&
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

    ! local variable
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: Ei,Ej,hi,hj,ri,rj,pi,pj,ui,uj,vi,vj
    real(DP) :: H_ij,anorm,aux,c_ij,q_ij,u_ij,v_ij,vel_ij
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,6,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,7,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,9,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,10,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,11,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,12,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,13,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,14,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,15,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,16,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,6,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,7,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,9,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,10,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,11,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,12,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,13,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,14,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,15,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,16,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = DCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Compute densities
        ri = DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute pressures
        pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute enthalpies
        hi = (TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2)
        q_ij   = DCONST(0.5)*(u_ij*u_ij+v_ij*v_ij)

        ! Compute the speed of sound
        c_ij = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))

        ! Compute scalar dissipation
        aux = dscale * (abs(vel_ij) + anorm*c_ij)
        
        IDX3(DmatrixAtEdge, :,1,idx,_,_,_) = DCONST(0.0)
        IDX3(DmatrixAtEdge, 1,1,idx,_,_,_) = aux
        IDX3(DmatrixAtEdge, 6,1,idx,_,_,_) = aux
        IDX3(DmatrixAtEdge,11,1,idx,_,_,_) = aux
        IDX3(DmatrixAtEdge,16,1,idx,_,_,_) = aux

      else

        ! Nullify dissipation tensor
        IDX3(DmatrixAtEdge,:,1,idx,_,_,_) = DCONST(0.0)

      end if
    end do

  end subroutine hydro_calcMatScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDissMatD2d_sim(DdataAtEdge,&
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

    ! local variable
    real(DP), dimension(NVAR2D,NVAR2D) :: R_ij,L_ij
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj
    real(DP) :: H_ij,anorm,aux,cPow2_ij,c_ij,q_ij,u_ij,v_ij,vel_ij
    real(DP) :: l1,l2,l3,l4
    integer :: idx
#if defined(HYDRO_USE_ENTROPYFIX) && (HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX)
    real(DP) :: dtol,ci,cj,veli,velj
#endif


    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vii,_)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,_)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !---------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = DCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL_DP) then

        ! Normalise the skew-symmetric coefficient
        a = a/anorm

        ! Compute densities
        ri = DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute pressures
        pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute enthalpies
        hi = (TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1)+v_ij*a(2)
        q_ij   = DCONST(0.5)*(u_ij*u_ij+v_ij*v_ij)
        
        ! Compute speed of sound
        cPow2_ij = max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij     = sqrt(cPow2_ij)
        
        ! Diagonal matrix of eigenvalues
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)

#if defined(HYDRO_USE_ENTROPYFIX)

#if HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX

        ! Entropy-fix by Harten and Hyman
        ci = SOUNDSPEED3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        cj = SOUNDSPEED3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        veli = ui*a(1) + vi*a(2)
        velj = uj*a(1) + vj*a(2)

        dtol = max(DCONST(0.0), (vel_ij-c_ij) - (veli-ci), (velj-cj) - (vel_ij-c_ij) )
        
        if (l1 .lt. dtol)&
            l1 = DCONST(0.5)*((l1*l1)/dtol + dtol)
        
        dtol = max(DCONST(0.0), (vel_ij+c_ij) - (veli+ci), (velj+cj) - (vel_ij+c_ij) )

        if (l3 .lt. dtol)&
            l3 = DCONST(0.5)*((l3*l3)/dtol + dtol)

#elif HYDRO_USE_ENTROPYFIX == HARTEN_ENTROPYFIX

#ifndef HYDRO_HARTEN_ENTROPYFIX
#error "Value HYDRO_HARTEN_ENTROPYFIX is required!"
#else
        ! Entropy-fix by Harten
        if (l1 .lt. DCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l1 = DCONST(0.5)*((l1*l1)/DCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + DCONST(HYDRO_HARTEN_ENTROPYFIX))

        if (l3 .lt. DCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l3 = DCONST(0.5)*((l3*l3)/DCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + DCONST(HYDRO_HARTEN_ENTROPYFIX))
#endif
#else
#error "Invalid type of entropy fix!"
#endif
#endif

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
        
        R_ij(1,4) =  DCONST(0.0)
        R_ij(2,4) =  l4*a(2)
        R_ij(3,4) = -l4*a(1)
        R_ij(4,4) =  l4*(u_ij*a(2)-v_ij*a(1))
        
        ! Matrix of left eigenvectors
        L_ij(1,1) =  DCONST(0.5)*(((HYDRO_GAMMA)-DCONST(1.0))*q_ij+c_ij*vel_ij)/cPow2_ij
        L_ij(2,1) =  (cPow2_ij-((HYDRO_GAMMA)-DCONST(1.0))*q_ij)/cPow2_ij
        L_ij(3,1) =  DCONST(0.5)*(((HYDRO_GAMMA)-DCONST(1.0))*q_ij-c_ij*vel_ij)/cPow2_ij
        L_ij(4,1) =  v_ij*a(1)-u_ij*a(2)
        
        L_ij(1,2) =  DCONST(0.5)*(-((HYDRO_GAMMA)-DCONST(1.0))*u_ij-c_ij*a(1))/cPow2_ij
        L_ij(2,2) =  ((HYDRO_GAMMA)-DCONST(1.0))*u_ij/cPow2_ij
        L_ij(3,2) =  DCONST(0.5)*(-((HYDRO_GAMMA)-DCONST(1.0))*u_ij+c_ij*a(1))/cPow2_ij
        L_ij(4,2) =  a(2)

        L_ij(1,3) =  DCONST(0.5)*(-((HYDRO_GAMMA)-DCONST(1.0))*v_ij-c_ij*a(2))/cPow2_ij
        L_ij(2,3) =  ((HYDRO_GAMMA)-DCONST(1.0))*v_ij/cPow2_ij
        L_ij(3,3) =  DCONST(0.5)*(-((HYDRO_GAMMA)-DCONST(1.0))*v_ij+c_ij*a(2))/cPow2_ij
        L_ij(4,3) = -a(1)
        
        L_ij(1,4) =  DCONST(0.5)*((HYDRO_GAMMA)-DCONST(1.0))/cPow2_ij
        L_ij(2,4) = -((HYDRO_GAMMA)-DCONST(1.0))/cPow2_ij
        L_ij(3,4) =  DCONST(0.5)*((HYDRO_GAMMA)-DCONST(1.0))/cPow2_ij
        L_ij(4,4) =  DCONST(0.0)
        
        ! Include scaling parameter
        anorm = dscale*anorm
        
        ! Compute tensorial dissipation D_ij = diag(R_ij*|Lbd_ij|*L_ij)*I
        IDX3(DmatrixAtEdge,:,1,idx,_,_,_) = DCONST(0.0)
        IDX3(DmatrixAtEdge,1,1,idx,_,_,_) = anorm*( R_ij(1,1)*L_ij(1,1)+&
                                                          R_ij(1,2)*L_ij(2,1)+&
                                                          R_ij(1,3)*L_ij(3,1)+&
                                                          R_ij(1,4)*L_ij(4,1)  )
        IDX3(DmatrixAtEdge,2,1,idx,_,_,_) = anorm*( R_ij(2,1)*L_ij(1,2)+&
                                                          R_ij(2,2)*L_ij(2,2)+&
                                                          R_ij(2,3)*L_ij(3,2)+&
                                                          R_ij(2,4)*L_ij(4,2)  )
        IDX3(DmatrixAtEdge,3,1,idx,_,_,_) = anorm*( R_ij(3,1)*L_ij(1,3)+&
                                                          R_ij(3,2)*L_ij(2,3)+&
                                                          R_ij(3,3)*L_ij(3,3)+&
                                                          R_ij(3,4)*L_ij(4,3)  )
        IDX3(DmatrixAtEdge,4,1,idx,_,_,_) = anorm*( R_ij(4,1)*L_ij(1,4)+&
                                                          R_ij(4,2)*L_ij(2,4)+&
                                                          R_ij(4,3)*L_ij(3,4)+&
                                                          R_ij(4,4)*L_ij(4,4)  )
      else
        
        ! Nullify dissipation tensor
        DmatrixAtEdge(:,1,idx) = DCONST(0.0)
        
      end if
    end do

  end subroutine hydro_calcMatRoeDissMatD2d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDiss2d_sim(DdataAtEdge,&
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

    ! local variable
    real(DP), dimension(NVAR2D,NVAR2D) :: R_ij,L_ij
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: Ei,Ej,hi,hj,pi,pj,ri,rj,ui,uj,vi,vj
    real(DP) :: H_ij,anorm,aux,cPow2_ij,c_ij,q_ij,u_ij,v_ij,vel_ij
    real(DP) :: l1,l2,l3,l4
    integer :: idx,i,j,k
#if defined(HYDRO_USE_ENTROPYFIX) && (HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX)
    real(DP) :: dtol,ci,cj,veli,velj
#endif


    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,6,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,7,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,9,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,10,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,11,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,12,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,13,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,14,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,15,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,16,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,6,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,7,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,9,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,10,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,11,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,12,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,13,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,14,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,15,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,16,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = DCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Normalise the skew-symmetric coefficient
        a = a/anorm

        ! Compute densities
        ri = DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute pressures
        pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute enthalpies
        hi = (TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1)+v_ij*a(2)
        q_ij   = DCONST(0.5)*(u_ij*u_ij+v_ij*v_ij)

        ! Compute speed of sound
        cPow2_ij = max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij     = sqrt(cPow2_ij)
        
        ! Diagonal matrix of eigenvalues
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)
        
#if defined(HYDRO_USE_ENTROPYFIX)

#if HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX

        ! Entropy-fix by Harten and Hyman
        ci = SOUNDSPEED3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        cj = SOUNDSPEED3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        veli = ui*a(1) + vi*a(2)
        velj = uj*a(1) + vj*a(2)

        dtol = max(DCONST(0.0), (vel_ij-c_ij) - (veli-ci), (velj-cj) - (vel_ij-c_ij) )
        
        if (l1 .lt. dtol)&
            l1 = DCONST(0.5)*((l1*l1)/dtol + dtol)
        
        dtol = max(DCONST(0.0), (vel_ij+c_ij) - (veli+ci), (velj+cj) - (vel_ij+c_ij) )

        if (l3 .lt. dtol)&
            l3 = DCONST(0.5)*((l3*l3)/dtol + dtol)

#elif HYDRO_USE_ENTROPYFIX == HARTEN_ENTROPYFIX

#ifndef HYDRO_HARTEN_ENTROPYFIX
#error "Value HYDRO_HARTEN_ENTROPYFIX is required!"
#else
        ! Entropy-fix by Harten
        if (l1 .lt. DCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l1 = DCONST(0.5)*((l1*l1)/DCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + DCONST(HYDRO_HARTEN_ENTROPYFIX))

        if (l3 .lt. DCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l3 = DCONST(0.5)*((l3*l3)/DCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + DCONST(HYDRO_HARTEN_ENTROPYFIX))
#endif
#else
#error "Invalid type of entropy fix!"
#endif
#endif

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
        
        R_ij(1,4) =  DCONST(0.0)
        R_ij(2,4) =  l4*a(2)
        R_ij(3,4) = -l4*a(1)
        R_ij(4,4) =  l4*(u_ij*a(2)-v_ij*a(1))
        
        ! Matrix of left eigenvectors
        L_ij(1,1) =  DCONST(0.5)*(((HYDRO_GAMMA)-DCONST(1.0))*q_ij+c_ij*vel_ij)/cPow2_ij
        L_ij(2,1) =  (cPow2_ij-((HYDRO_GAMMA)-DCONST(1.0))*q_ij)/cPow2_ij
        L_ij(3,1) =  DCONST(0.5)*(((HYDRO_GAMMA)-DCONST(1.0))*q_ij-c_ij*vel_ij)/cPow2_ij
        L_ij(4,1) =  v_ij*a(1)-u_ij*a(2)
        
        L_ij(1,2) =  DCONST(0.5)*(-((HYDRO_GAMMA)-DCONST(1.0))*u_ij-c_ij*a(1))/cPow2_ij
        L_ij(2,2) =  ((HYDRO_GAMMA)-DCONST(1.0))*u_ij/cPow2_ij
        L_ij(3,2) =  DCONST(0.5)*(-((HYDRO_GAMMA)-DCONST(1.0))*u_ij+c_ij*a(1))/cPow2_ij
        L_ij(4,2) =  a(2)
        
        L_ij(1,3) =  DCONST(0.5)*(-((HYDRO_GAMMA)-DCONST(1.0))*v_ij-c_ij*a(2))/cPow2_ij
        L_ij(2,3) =  ((HYDRO_GAMMA)-DCONST(1.0))*v_ij/cPow2_ij
        L_ij(3,3) =  DCONST(0.5)*(-((HYDRO_GAMMA)-DCONST(1.0))*v_ij+c_ij*a(2))/cPow2_ij
        L_ij(4,3) = -a(1)
        
        L_ij(1,4) =  ((HYDRO_GAMMA)-DCONST(1.0))/DCONST(2.0)/cPow2_ij
        L_ij(2,4) = -((HYDRO_GAMMA)-DCONST(1.0))/cPow2_ij
        L_ij(3,4) =  ((HYDRO_GAMMA)-DCONST(1.0))/DCONST(2.0)/cPow2_ij
        L_ij(4,4) =  DCONST(0.0)
        
        ! Include scaling parameter
        anorm = dscale*anorm

        ! Compute tensorial dissipation D_ij = R_ij*|Lbd_ij|*L_ij
        do i = 1, NVAR2D
          do j = 1, NVAR2D
            aux = DCONST(0.0)
            do k = 1, NVAR2D
              aux = aux + R_ij(i,k)*L_ij(k,j)
            end do
            IDX3(DmatrixAtEdge,NVAR2D*(j-1)+i,1,idx,_,_,_) = anorm*aux
          end do
        end do
        
      else
        
        ! Nullify dissipation tensor
        IDX3(DmatrixAtEdge,:,1,idx,_,_,_) = DCONST(0.0)
        
      end if
    end do

  end subroutine hydro_calcMatRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDissMatD2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 2D and applies the scalar artificial viscosities of
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
    real(DP) :: Ei,Ej,ci,cj,ui,uj,vi,vj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vii,_)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,_)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,_)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate scalar artificial dissipation of Rusanov-type
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(HYDRO_GAMMA)*&
          (Ei-DCONST(0.5)*(ui*ui+vi*vi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(HYDRO_GAMMA)*&
          (Ej-DCONST(0.5)*(uj*uj+vj*vj)), SYS_EPSREAL_DP))
      
      ! Compute dissipation tensor
      IDX3(DmatrixAtEdge,:,1,idx,_,_,_) = dscale *&
          max( abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*uj+&
                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*vj) +&
                   sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)+&
                        POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2))*cj,&
               abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*ui+&
                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*vi) +&
                   sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)+&
                        POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2))*ci )
    end do

  end subroutine hydro_calcMatRusDissMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDiss2d_sim(DdataAtEdge,&
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

    ! local variable
    real(DP) :: Ei,Ej,aux,ci,cj,ui,uj,vi,vj
    integer :: idx


    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,6,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,7,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,9,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,10,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,11,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,12,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,13,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,14,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,15,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,16,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                               IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),ui,vi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,4,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,5,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,6,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,7,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,9,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,10,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,11,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,12,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)

      IDX3(DmatrixAtEdge,13,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,14,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,15,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      IDX3(DmatrixAtEdge,16,2,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),uj,vj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX11_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX21_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX31_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,4,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX41_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,5,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX12_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,6,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX22_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,7,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX32_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX42_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,9,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX13_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,10,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX23_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,11,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX33_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,12,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX43_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)

      IDX3(DmatrixAtEdge,13,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX14_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,14,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX24_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,15,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX34_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
      IDX3(DmatrixAtEdge,16,3,idx,_,_,_) =&
          INVISCIDFLUXJACOBIMATRIX44_2D(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),MYNEWLINE \
                                                IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),ui,vi,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate scalar artificial dissipation of Rusanov-type
      !---------------------------------------------------------------------------

      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(HYDRO_GAMMA)*&
          (Ei-DCONST(0.5)*(ui*ui+vi*vi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(HYDRO_GAMMA)*&
          (Ej-DCONST(0.5)*(uj*uj+vj*vj)), SYS_EPSREAL_DP))
      
      ! Compute dissipation tensor
      aux = dscale *&
          max( abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*uj+&
                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*vj) +&
                   sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)+&
                        POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2))*cj,&
               abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*ui+&
                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*vi) +&
                   sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)+&
                        POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2))*ci )

      IDX3(DmatrixAtEdge, :,1,idx,_,_,_) = DCONST(0.0)
      IDX3(DmatrixAtEdge, 1,1,idx,_,_,_) = aux
      IDX3(DmatrixAtEdge, 6,1,idx,_,_,_) = aux
      IDX3(DmatrixAtEdge,11,1,idx,_,_,_) = aux
      IDX3(DmatrixAtEdge,16,1,idx,_,_,_) = aux
    end do

  end subroutine hydro_calcMatRusDiss2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcDivMatGalMatD2d_cuda(rgroupFEMSet, rx, rmatrix,&
      dscale, bclear, fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
      cconstrType, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine assembles the block-diagonal part of the
    ! standard Galerkin operator in 2D.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the matrix entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

    ! OPTIONAL: callback function to compute local matrices
    include '../../../../../kernel/PDEOperators/intf_calcMatrixDiagSys_sim.inc'
    include '../../../../../kernel/PDEOperators/intf_calcMatrixSys_sim.inc'
    optional :: fcb_calcMatrixDiagSys_sim
    optional :: fcb_calcMatrixSys_sim
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixBlock), intent(inout) :: rmatrix
    
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT

    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: IEDGEmax,IEDGEset,igroup,ccType
    type(C_PTR), dimension(:,:), pointer :: cptr_Da
    type(C_PTR) :: cptr_DcoeffsAtEdge, cptr_DcoeffsAtDiag
    type(C_PTR) :: cptr_IedgeList, cptr_IdiagList
    type(C_PTR) :: cptr_Dx
    integer(I64) :: istream

    ccType = GFEM_MATC_CONSISTENT
    if (present(cconstrType)) ccType = cconstrType

    ! Create CUDA stream
    call coproc_createStream(istream)

    ! Get pointers to memory blocks on device
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IdiagList, cptr_IdiagList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IedgeList, cptr_IedgeList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtDiag, cptr_DcoeffsAtDiag)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtEdge, cptr_DcoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    if (.not.storage_isAssociated(cptr_Dx)) then
      call lsysbl_copyH2D_Vector(rx, .false., .false., istream)
      call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    end if
   
    ! Make sure that the destination matrix rmatrix exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Matrix(rmatrix, .true., .false., istream, cptr_Da)
    
    ! Set pointer
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)

    ! Use callback function to compute diagonal matrix entries
    call hydro_calcMatDiagMatD2d_cuda(cptr_DcoeffsAtDiag, cptr_IdiagList,&
        cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
        rgroupFEMSet%NA, rgroupFEMSet%ncoeffsAtDiag, istream)

    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronise memory access
    do igroup = 1, size(p_IedgeListIdx)-1
      
      ! Do nothing for empty groups
      if (p_IedgeListIdx(igroup+1)-p_IedgeListIdx(igroup) .le. 0) cycle
      
      ! Get position of first edge in group
      IEDGEset = p_IedgeListIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IedgeListIdx(igroup+1)-1
      
      ! Use callback function to compute off-diagonal matrix entries
      call hydro_calcMatGalMatD2d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
          cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
          rgroupFEMSet%NA, rgroupFEMSet%NEDGE, rgroupFEMSet%ncoeffsAtEdge,&
          IEDGEmax-IEDGEset+1, IEDGEset, ccType, istream)
    end do
    
    deallocate(cptr_Da)

    ! Transfer destination matrix back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Matrix(rmatrix, bclear, .false., istream)
       
    ! Ensure data consistency
    call coproc_synchronizeStream(istream)
    call coproc_destroyStream(istream)

#else
    
    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivMatGalMatD2d_cuda')
    call sys_halt()
    
#endif
    
  end subroutine hydro_calcDivMatGalMatD2d_cuda

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcDivMatGalerkin2d_cuda(rgroupFEMSet, rx, rmatrix,&
      dscale, bclear, fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
      cconstrType, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine assembles the standard Galerkin operator in 2D
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the matrix entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

    ! OPTIONAL: callback function to compute local matrices
    include '../../../../../kernel/PDEOperators/intf_calcMatrixDiagSys_sim.inc'
    include '../../../../../kernel/PDEOperators/intf_calcMatrixSys_sim.inc'
    optional :: fcb_calcMatrixDiagSys_sim
    optional :: fcb_calcMatrixSys_sim
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT

    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: IEDGEmax,IEDGEset,igroup,ccType
    type(C_PTR), dimension(:,:), pointer :: cptr_Da
    type(C_PTR) :: cptr_DcoeffsAtEdge, cptr_DcoeffsAtDiag
    type(C_PTR) :: cptr_IedgeList, cptr_IdiagList
    type(C_PTR) :: cptr_Dx
    integer(I64) :: istream

    ccType = GFEM_MATC_CONSISTENT
    if (present(cconstrType)) ccType = cconstrType
    
    ! Create CUDA stream
    call coproc_createStream(istream)

    ! Get pointers to memory blocks on device
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IdiagList, cptr_IdiagList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IedgeList, cptr_IedgeList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtDiag, cptr_DcoeffsAtDiag)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtEdge, cptr_DcoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    if (.not.storage_isAssociated(cptr_Dx)) then
      call lsysbl_copyH2D_Vector(rx, .false., .false., istream)
      call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    end if
   
    ! Make sure that the destination matrix rmatrix exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Matrix(rmatrix, .true., .false., istream, cptr_Da)

    ! Set pointer
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
    
    ! Use callback function to compute diagonal matrix entries
    call hydro_calcMatDiag2d_cuda(cptr_DcoeffsAtDiag, cptr_IdiagList,&
        cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
        rgroupFEMSet%NA, rgroupFEMSet%ncoeffsAtDiag, istream)

    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronise memory access
    do igroup = 1, size(p_IedgeListIdx)-1
      
      ! Do nothing for empty groups
      if (p_IedgeListIdx(igroup+1)-p_IedgeListIdx(igroup) .le. 0) cycle
      
      ! Get position of first edge in group
      IEDGEset = p_IedgeListIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IedgeListIdx(igroup+1)-1
      
      ! Use callback function to compute off-diagonal matrix entries
      call hydro_calcMatGalerkin2d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
          cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
          rgroupFEMSet%NA, rgroupFEMSet%NEDGE, rgroupFEMSet%ncoeffsAtEdge,&
          IEDGEmax-IEDGEset+1, IEDGEset, ccType, istream)
    end do

    deallocate(cptr_Da)
    
    ! Transfer destination matrix back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Matrix(rmatrix, bclear, .false., istream)

    ! Ensure data consistency
    call coproc_synchronizeStream(istream)
    call coproc_destroyStream(istream)
#else
    
    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivMatGalerkin2d_cuda')
    call sys_halt()
    
#endif
    
  end subroutine hydro_calcDivMatGalerkin2d_cuda

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcDivMatScDissMatD2d_cuda(rgroupFEMSet, rx, rmatrix,&
      dscale, bclear, fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
      cconstrType, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine assembles the block-diagonal part of the
    ! discrete low-order operator in 2D using scalar artificial
    ! viscosities proportional to the spectral radius (largest
    ! eigenvalue) of the Roe-matrix.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the matrix entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

    ! OPTIONAL: callback function to compute local matrices
    include '../../../../../kernel/PDEOperators/intf_calcMatrixDiagSys_sim.inc'
    include '../../../../../kernel/PDEOperators/intf_calcMatrixSys_sim.inc'
    optional :: fcb_calcMatrixDiagSys_sim
    optional :: fcb_calcMatrixSys_sim
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT

    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: IEDGEmax,IEDGEset,igroup,ccType
    type(C_PTR), dimension(:,:), pointer :: cptr_Da
    type(C_PTR) :: cptr_DcoeffsAtEdge, cptr_DcoeffsAtDiag
    type(C_PTR) :: cptr_IedgeList, cptr_IdiagList
    type(C_PTR) :: cptr_Dx
    integer(I64) :: istream

    ccType = GFEM_MATC_CONSISTENT
    if (present(cconstrType)) ccType = cconstrType

    ! Create CUDA stream
    call coproc_createStream(istream)

    ! Get pointers to memory blocks on device
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IdiagList, cptr_IdiagList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IedgeList, cptr_IedgeList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtDiag, cptr_DcoeffsAtDiag)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtEdge, cptr_DcoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    if (.not.storage_isAssociated(cptr_Dx)) then
      call lsysbl_copyH2D_Vector(rx, .false., .false., istream)
      call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    end if
   
    ! Make sure that the destination matrix rmatrix exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Matrix(rmatrix, .true., .false., istream, cptr_Da)
    
    ! Set pointer
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)

    ! Use callback function to compute diagonal matrix entries
    call hydro_calcMatDiagMatD2d_cuda(cptr_DcoeffsAtDiag, cptr_IdiagList,&
        cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
        rgroupFEMSet%NA, rgroupFEMSet%ncoeffsAtDiag, istream)

    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronise memory access
    do igroup = 1, size(p_IedgeListIdx)-1
      
      ! Do nothing for empty groups
      if (p_IedgeListIdx(igroup+1)-p_IedgeListIdx(igroup) .le. 0) cycle
      
      ! Get position of first edge in group
      IEDGEset = p_IedgeListIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IedgeListIdx(igroup+1)-1
      
      ! Use callback function to compute off-diagonal matrix entries
      call hydro_calcMatScDissMatD2d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
          cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
          rgroupFEMSet%NA, rgroupFEMSet%NEDGE, rgroupFEMSet%ncoeffsAtEdge,&
          IEDGEmax-IEDGEset+1, IEDGEset, ccType, istream)
    end do
    
    deallocate(cptr_Da)

    ! Transfer destination matrix back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Matrix(rmatrix, bclear, .false., istream)
       
    ! Ensure data consistency
    call coproc_synchronizeStream(istream)
    call coproc_destroyStream(istream)

#else
    
    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivMatScDissMatD2d_cuda')
    call sys_halt()
    
#endif
    
  end subroutine hydro_calcDivMatScDissMatD2d_cuda

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcDivMatScDiss2d_cuda(rgroupFEMSet, rx, rmatrix,&
      dscale, bclear, fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
      cconstrType, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine assembles the discrete low-order operator in 2D
    ! using scalar artificial viscosities proportional to the spectral
    ! radius (largest eigenvalue) of the Roe-matrix.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the matrix entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

    ! OPTIONAL: callback function to compute local matrices
    include '../../../../../kernel/PDEOperators/intf_calcMatrixDiagSys_sim.inc'
    include '../../../../../kernel/PDEOperators/intf_calcMatrixSys_sim.inc'
    optional :: fcb_calcMatrixDiagSys_sim
    optional :: fcb_calcMatrixSys_sim
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT

    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: IEDGEmax,IEDGEset,igroup,ccType
    type(C_PTR), dimension(:,:), pointer :: cptr_Da
    type(C_PTR) :: cptr_DcoeffsAtEdge, cptr_DcoeffsAtDiag
    type(C_PTR) :: cptr_IedgeList, cptr_IdiagList
    type(C_PTR) :: cptr_Dx
    integer(I64) :: istream

    ccType = GFEM_MATC_CONSISTENT
    if (present(cconstrType)) ccType = cconstrType

    ! Create CUDA stream
    call coproc_createStream(istream)

    ! Get pointers to memory blocks on device
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IdiagList, cptr_IdiagList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IedgeList, cptr_IedgeList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtDiag, cptr_DcoeffsAtDiag)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtEdge, cptr_DcoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    if (.not.storage_isAssociated(cptr_Dx)) then
      call lsysbl_copyH2D_Vector(rx, .false., .false., istream)
      call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    end if
   
    ! Make sure that the destination matrix rmatrix exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Matrix(rmatrix, .true., .false., istream, cptr_Da)
    
    ! Set pointer
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)

    ! Use callback function to compute diagonal matrix entries
    call hydro_calcMatDiag2d_cuda(cptr_DcoeffsAtDiag, cptr_IdiagList,&
        cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
        rgroupFEMSet%NA, rgroupFEMSet%ncoeffsAtDiag, istream)

    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronise memory access
    do igroup = 1, size(p_IedgeListIdx)-1
      
      ! Do nothing for empty groups
      if (p_IedgeListIdx(igroup+1)-p_IedgeListIdx(igroup) .le. 0) cycle
      
      ! Get position of first edge in group
      IEDGEset = p_IedgeListIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IedgeListIdx(igroup+1)-1
      
      ! Use callback function to compute off-diagonal matrix entries
      call hydro_calcMatScDiss2d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
          cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
          rgroupFEMSet%NA, rgroupFEMSet%NEDGE, rgroupFEMSet%ncoeffsAtEdge,&
          IEDGEmax-IEDGEset+1, IEDGEset, ccType, istream)
    end do
    
    deallocate(cptr_Da)

    ! Transfer destination matrix back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Matrix(rmatrix, bclear, .false., istream)
       
    ! Ensure data consistency
    call coproc_synchronizeStream(istream)
    call coproc_destroyStream(istream)

#else
    
    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivMatScDiss2d_cuda')
    call sys_halt()
    
#endif
    
  end subroutine hydro_calcDivMatScDiss2d_cuda

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcDivMatRoeDissMatD2d_cuda(rgroupFEMSet, rx, rmatrix,&
      dscale, bclear, fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
      cconstrType, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine assembles the block-diagonal part of the
    ! discrete low-order operator in 2D using tensorial artificial
    ! viscosities of Roe-type.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the matrix entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

    ! OPTIONAL: callback function to compute local matrices
    include '../../../../../kernel/PDEOperators/intf_calcMatrixDiagSys_sim.inc'
    include '../../../../../kernel/PDEOperators/intf_calcMatrixSys_sim.inc'
    optional :: fcb_calcMatrixDiagSys_sim
    optional :: fcb_calcMatrixSys_sim
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT
    
    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: IEDGEmax,IEDGEset,igroup,ccType
    type(C_PTR), dimension(:,:), pointer :: cptr_Da
    type(C_PTR) :: cptr_DcoeffsAtEdge, cptr_DcoeffsAtDiag
    type(C_PTR) :: cptr_IedgeList, cptr_IdiagList
    type(C_PTR) :: cptr_Dx
    integer(I64) :: istream

    ccType = GFEM_MATC_CONSISTENT
    if (present(cconstrType)) ccType = cconstrType

    ! Create CUDA stream
    call coproc_createStream(istream)

    ! Get pointers to memory blocks on device
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IdiagList, cptr_IdiagList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IedgeList, cptr_IedgeList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtDiag, cptr_DcoeffsAtDiag)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtEdge, cptr_DcoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    if (.not.storage_isAssociated(cptr_Dx)) then
      call lsysbl_copyH2D_Vector(rx, .false., .false., istream)
      call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    end if
   
    ! Make sure that the destination matrix rmatrix exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Matrix(rmatrix, .true., .false., istream, cptr_Da)
    
    ! Set pointer
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)

    ! Use callback function to compute diagonal matrix entries
    call hydro_calcMatDiagMatD2d_cuda(cptr_DcoeffsAtDiag, cptr_IdiagList,&
        cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
        rgroupFEMSet%NA, rgroupFEMSet%ncoeffsAtDiag, istream)

    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronise memory access
    do igroup = 1, size(p_IedgeListIdx)-1
      
      ! Do nothing for empty groups
      if (p_IedgeListIdx(igroup+1)-p_IedgeListIdx(igroup) .le. 0) cycle
      
      ! Get position of first edge in group
      IEDGEset = p_IedgeListIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IedgeListIdx(igroup+1)-1
      
      ! Use callback function to compute off-diagonal matrix entries
      call hydro_calcMatRoeDissMatD2d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
          cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
          rgroupFEMSet%NA, rgroupFEMSet%NEDGE, rgroupFEMSet%ncoeffsAtEdge,&
          IEDGEmax-IEDGEset+1, IEDGEset, ccType, istream)
    end do
    
    deallocate(cptr_Da)

    ! Transfer destination matrix back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Matrix(rmatrix, bclear, .false., istream)
       
    ! Ensure data consistency
    call coproc_synchronizeStream(istream)
    call coproc_destroyStream(istream)

#else
    
    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivMatRoeDissMatD2d_cuda')
    call sys_halt()
    
#endif
    
  end subroutine hydro_calcDivMatRoeDissMatD2d_cuda

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcDivMatRoeDiss2d_cuda(rgroupFEMSet, rx, rmatrix,&
      dscale, bclear, fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
      cconstrType, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine assembles the discrete low-order operator in 2D
    ! using tensorial artificial viscosities of Roe-type.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the matrix entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

    ! OPTIONAL: callback function to compute local matrices
    include '../../../../../kernel/PDEOperators/intf_calcMatrixDiagSys_sim.inc'
    include '../../../../../kernel/PDEOperators/intf_calcMatrixSys_sim.inc'
    optional :: fcb_calcMatrixDiagSys_sim
    optional :: fcb_calcMatrixSys_sim
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT
    
    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: IEDGEmax,IEDGEset,igroup,ccType
    type(C_PTR), dimension(:,:), pointer :: cptr_Da
    type(C_PTR) :: cptr_DcoeffsAtEdge, cptr_DcoeffsAtDiag
    type(C_PTR) :: cptr_IedgeList, cptr_IdiagList
    type(C_PTR) :: cptr_Dx
    integer(I64) :: istream

    ccType = GFEM_MATC_CONSISTENT
    if (present(cconstrType)) ccType = cconstrType

    ! Create CUDA stream
    call coproc_createStream(istream)

    ! Get pointers to memory blocks on device
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IdiagList, cptr_IdiagList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IedgeList, cptr_IedgeList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtDiag, cptr_DcoeffsAtDiag)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtEdge, cptr_DcoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    if (.not.storage_isAssociated(cptr_Dx)) then
      call lsysbl_copyH2D_Vector(rx, .false., .false., istream)
      call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    end if
   
    ! Make sure that the destination matrix rmatrix exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Matrix(rmatrix, .true., .false., istream, cptr_Da)
    
    ! Set pointer
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)

    ! Use callback function to compute diagonal matrix entries
    call hydro_calcMatDiag2d_cuda(cptr_DcoeffsAtDiag, cptr_IdiagList,&
        cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
        rgroupFEMSet%NA, rgroupFEMSet%ncoeffsAtDiag, istream)

    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronise memory access
    do igroup = 1, size(p_IedgeListIdx)-1
      
      ! Do nothing for empty groups
      if (p_IedgeListIdx(igroup+1)-p_IedgeListIdx(igroup) .le. 0) cycle
      
      ! Get position of first edge in group
      IEDGEset = p_IedgeListIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IedgeListIdx(igroup+1)-1
      
      ! Use callback function to compute off-diagonal matrix entries
      call hydro_calcMatRoeDiss2d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
          cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
          rgroupFEMSet%NA, rgroupFEMSet%NEDGE, rgroupFEMSet%ncoeffsAtEdge,&
          IEDGEmax-IEDGEset+1, IEDGEset, ccType, istream)
    end do
    
    deallocate(cptr_Da)

    ! Transfer destination matrix back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Matrix(rmatrix, bclear, .false., istream)
       
    ! Ensure data consistency
    call coproc_synchronizeStream(istream)
    call coproc_destroyStream(istream)

#else
    
    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivMatRoeDiss2d_cuda')
    call sys_halt()
    
#endif
    
  end subroutine hydro_calcDivMatRoeDiss2d_cuda

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcDivMatRusDissMatD2d_cuda(rgroupFEMSet, rx, rmatrix,&
      dscale, bclear, fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
      cconstrType, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine assembles the block-diagonal part of the
    ! discrete low-order operator in 2D using scalar artificial
    ! viscosities of Rusanov-type.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the matrix entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

    ! OPTIONAL: callback function to compute local matrices
    include '../../../../../kernel/PDEOperators/intf_calcMatrixDiagSys_sim.inc'
    include '../../../../../kernel/PDEOperators/intf_calcMatrixSys_sim.inc'
    optional :: fcb_calcMatrixDiagSys_sim
    optional :: fcb_calcMatrixSys_sim
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT

    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: IEDGEmax,IEDGEset,igroup,ccType
    type(C_PTR), dimension(:,:), pointer :: cptr_Da
    type(C_PTR) :: cptr_DcoeffsAtEdge, cptr_DcoeffsAtDiag
    type(C_PTR) :: cptr_IedgeList, cptr_IdiagList
    type(C_PTR) :: cptr_Dx
    integer(I64) :: istream

    ccType = GFEM_MATC_CONSISTENT
    if (present(cconstrType)) ccType = cconstrType

    ! Create CUDA stream
    call coproc_createStream(istream)

    ! Get pointers to memory blocks on device
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IdiagList, cptr_IdiagList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IedgeList, cptr_IedgeList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtDiag, cptr_DcoeffsAtDiag)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtEdge, cptr_DcoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    if (.not.storage_isAssociated(cptr_Dx)) then
      call lsysbl_copyH2D_Vector(rx, .false., .false., istream)
      call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    end if
   
    ! Make sure that the destination matrix rmatrix exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Matrix(rmatrix, .true., .false., istream, cptr_Da)
    
    ! Set pointer
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)

    ! Use callback function to compute diagonal matrix entries
    call hydro_calcMatDiagMatD2d_cuda(cptr_DcoeffsAtDiag, cptr_IdiagList,&
        cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
        rgroupFEMSet%NA, rgroupFEMSet%ncoeffsAtDiag, istream)

    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronise memory access
    do igroup = 1, size(p_IedgeListIdx)-1
      
      ! Do nothing for empty groups
      if (p_IedgeListIdx(igroup+1)-p_IedgeListIdx(igroup) .le. 0) cycle
      
      ! Get position of first edge in group
      IEDGEset = p_IedgeListIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IedgeListIdx(igroup+1)-1
      
      ! Use callback function to compute off-diagonal matrix entries
      call hydro_calcMatRusDissMatD2d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
          cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
          rgroupFEMSet%NA, rgroupFEMSet%NEDGE, rgroupFEMSet%ncoeffsAtEdge,&
          IEDGEmax-IEDGEset+1, IEDGEset, ccType, istream)
    end do
    
    deallocate(cptr_Da)

    ! Transfer destination matrix back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Matrix(rmatrix, bclear, .false., istream)
       
    ! Ensure data consistency
    call coproc_synchronizeStream(istream)
    call coproc_destroyStream(istream)

#else
    
    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivMatRusDissMatD2d_cuda')
    call sys_halt()
    
#endif
    
  end subroutine hydro_calcDivMatRusDissMatD2d_cuda

!*****************************************************************************

!<subroutine>

  subroutine hydro_calcDivMatRusDiss2d_cuda(rgroupFEMSet, rx, rmatrix,&
      dscale, bclear, fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
      cconstrType, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine assembles the discrete low-order operator in 2D
    ! using scalar artificial viscosities of Rusanov-type.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the matrix entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

    ! OPTIONAL: callback function to compute local matrices
    include '../../../../../kernel/PDEOperators/intf_calcMatrixDiagSys_sim.inc'
    include '../../../../../kernel/PDEOperators/intf_calcMatrixSys_sim.inc'
    optional :: fcb_calcMatrixDiagSys_sim
    optional :: fcb_calcMatrixSys_sim
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

#ifdef ENABLE_COPROCESSOR_SUPPORT

    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: IEDGEmax,IEDGEset,igroup,ccType
    type(C_PTR), dimension(:,:), pointer :: cptr_Da
    type(C_PTR) :: cptr_DcoeffsAtEdge, cptr_DcoeffsAtDiag
    type(C_PTR) :: cptr_IedgeList, cptr_IdiagList
    type(C_PTR) :: cptr_Dx
    integer(I64) :: istream

    ccType = GFEM_MATC_CONSISTENT
    if (present(cconstrType)) ccType = cconstrType

    ! Create CUDA stream
    call coproc_createStream(istream)

    ! Get pointers to memory blocks on device
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IdiagList, cptr_IdiagList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_IedgeList, cptr_IedgeList)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtDiag, cptr_DcoeffsAtDiag)
    call storage_getMemPtrOnDevice(rgroupFEMSet%h_CoeffsAtEdge, cptr_DcoeffsAtEdge)

    ! In the very first call to this routine, the source vector may be
    ! uninitialised on the device. In this case, we have to do it here.
    call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    if (.not.storage_isAssociated(cptr_Dx)) then
      call lsysbl_copyH2D_Vector(rx, .false., .false., istream)
      call storage_getMemPtrOnDevice(rx%h_Ddata, cptr_Dx)
    end if
   
    ! Make sure that the destination matrix rmatrix exists on the
    ! coprocessor device and is initialised by zeros
    call lsysbl_copyH2D_Matrix(rmatrix, .true., .false., istream, cptr_Da)
    
    ! Set pointer
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)

    ! Use callback function to compute diagonal matrix entries
    call hydro_calcMatDiag2d_cuda(cptr_DcoeffsAtDiag, cptr_IdiagList,&
        cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
        rgroupFEMSet%NA, rgroupFEMSet%ncoeffsAtDiag, istream)

    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronise memory access
    do igroup = 1, size(p_IedgeListIdx)-1
      
      ! Do nothing for empty groups
      if (p_IedgeListIdx(igroup+1)-p_IedgeListIdx(igroup) .le. 0) cycle
      
      ! Get position of first edge in group
      IEDGEset = p_IedgeListIdx(igroup)
      
      ! Get position of last edge in group
      IEDGEmax = p_IedgeListIdx(igroup+1)-1
      
      ! Use callback function to compute off-diagonal matrix entries
      call hydro_calcMatRusDiss2d_cuda(cptr_DcoeffsAtEdge, cptr_IedgeList,&
          cptr_Dx, cptr_Da, dscale, rx%nblocks, rgroupFEMSet%NEQ,&
          rgroupFEMSet%NA, rgroupFEMSet%NEDGE, rgroupFEMSet%ncoeffsAtEdge,&
          IEDGEmax-IEDGEset+1, IEDGEset, ccType, istream)
    end do
    
    deallocate(cptr_Da)

    ! Transfer destination matrix back to host memory. If bclear is
    ! .TRUE. then the content of the host memory can be overwritten;
    ! otherwise we need to copy-add the content from device memory
    call lsysbl_copyD2H_Matrix(rmatrix, bclear, .false., istream)
       
    ! Ensure data consistency
    call coproc_synchronizeStream(istream)
    call coproc_destroyStream(istream)

#else
    
    call output_line('Coprocessor support is disabled!',&
        OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivMatRusDiss2d_cuda')
    call sys_halt()
    
#endif
    
  end subroutine hydro_calcDivMatRusDiss2d_cuda

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
      if (present(DcharVariablesAtEdge))     DcharVariablesAtEdge     = DCONST(0.0)
      if (present(DeigenvaluesAtEdge))       DeigenvaluesAtEdge       = DCONST(0.0)
      if (present(DrightEigenvectorsAtEdge)) DrightEigenvectorsAtEdge = DCONST(0.0)
      if (present(DleftEigenvectorsAtEdge))  DleftEigenvectorsAtEdge  = DCONST(0.0)

      ! That`s it
      return
    end if

    ! Compute normalised weighting coefficient
    a = Dweight/anorm

    ! Do we have to compute characteristic variables
    if (present(DcharVariablesAtEdge)) then
      do idx = 1, nedges
        
        ! Compute velocities
        ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        
        uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute densities
        ri = DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute pressures
        pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute enthalpies
        hi = (TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        q_ij  = DCONST(0.5)*(u_ij*u_ij+v_ij*v_ij)

        ! Compute speed of sound
        cPow2_ij = max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij = sqrt(cPow2_ij)
        aux  = u_ij*a(1)+v_ij*a(2)

        ! Compute solution difference U_j-U_i
        Diff = IDX3(DdataAtEdge,:,2,idx,_,_,_)-&
               IDX3(DdataAtEdge,:,1,idx,_,_,_)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-DCONST(1.0))*(q_ij*Diff(1)&
                                           -u_ij*Diff(2)&
                                           -v_ij*Diff(3)&
                                                +Diff(4))/DCONST(2.0)/cPow2_ij
        aux2 = (aux*Diff(1)&
              -a(1)*Diff(2)&
              -a(2)*Diff(3))/DCONST(2.0)/c_ij

        ! Compute characteristic variables
        IDX2(DcharVariablesAtEdge,1,idx,_,_) = anorm * (aux1 + aux2)
        IDX2(DcharVariablesAtEdge,2,idx,_,_) = anorm *&
            ((DCONST(1.0)-((HYDRO_GAMMA)-DCONST(1.0))*q_ij/cPow2_ij)*Diff(1)+&
                                   ((HYDRO_GAMMA)-DCONST(1.0))*(u_ij*Diff(2)+&
                                                                v_ij*Diff(3)-&
                                                                     Diff(4))/cPow2_ij)
        IDX2(DcharVariablesAtEdge,3,idx,_,_) = anorm * (aux1 - aux2)
        IDX2(DcharVariablesAtEdge,4,idx,_,_) = anorm * ((a(1)*v_ij-a(2)*u_ij)*Diff(1)+&
                                                                         a(2)*Diff(2)-&
                                                                         a(1)*Diff(3))
      end do
    end if

    ! Do we have to compute eigenvalues
    if (present(DeigenvaluesAtEdge)) then
      do idx = 1, nedges

        ! Compute velocities
        ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        
        uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute densities
        ri = DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute pressures
        pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute enthalpies
        hi = (TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        q_ij = DCONST(0.5)*(u_ij*u_ij+v_ij*v_ij)
        c_ij = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))
        aux  = a(1)*u_ij+a(2)*v_ij

        ! Compute eigenvalues
        IDX2(DeigenvaluesAtEdge,1,idx,_,_) = aux-c_ij
        IDX2(DeigenvaluesAtEdge,2,idx,_,_) = aux
        IDX2(DeigenvaluesAtEdge,3,idx,_,_) = aux+c_ij
        IDX2(DeigenvaluesAtEdge,4,idx,_,_) = aux
      end do
    end if

    ! Do we have to compute right eigenvectors
    if (present(DrightEigenvectorsAtEdge)) then
      do idx = 1, nedges

        ! Compute velocities
        ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        
        uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute densities
        ri = DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute pressures
        pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute enthalpies
        hi = (TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        q_ij = DCONST(0.5)*(u_ij*u_ij+v_ij*v_ij)
        c_ij = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))
        aux  = a(1)*u_ij+a(2)*v_ij

        ! Compute right eigenvectors
        IDX2(DrightEigenvectorsAtEdge, 1,idx,_,_) =  DCONST(1.0)
        IDX2(DrightEigenvectorsAtEdge, 2,idx,_,_) =  u_ij-c_ij*a(1)
        IDX2(DrightEigenvectorsAtEdge, 3,idx,_,_) =  v_ij-c_ij*a(2)
        IDX2(DrightEigenvectorsAtEdge, 4,idx,_,_) =  H_ij-c_ij*aux

        IDX2(DrightEigenvectorsAtEdge, 5,idx,_,_) =  DCONST(1.0)
        IDX2(DrightEigenvectorsAtEdge, 6,idx,_,_) =  u_ij
        IDX2(DrightEigenvectorsAtEdge, 7,idx,_,_) =  v_ij
        IDX2(DrightEigenvectorsAtEdge, 8,idx,_,_) =  q_ij

        IDX2(DrightEigenvectorsAtEdge, 9,idx,_,_) =  DCONST(1.0)
        IDX2(DrightEigenvectorsAtEdge,10,idx,_,_) =  u_ij+c_ij*a(1)
        IDX2(DrightEigenvectorsAtEdge,11,idx,_,_) =  v_ij+c_ij*a(2)
        IDX2(DrightEigenvectorsAtEdge,12,idx,_,_) =  H_ij+c_ij*aux

        IDX2(DrightEigenvectorsAtEdge,13,idx,_,_) =  DCONST(0.0)
        IDX2(DrightEigenvectorsAtEdge,14,idx,_,_) =  a(2)
        IDX2(DrightEigenvectorsAtEdge,15,idx,_,_) = -a(1)
        IDX2(DrightEigenvectorsAtEdge,16,idx,_,_) =  u_ij*a(2)-v_ij*a(1)
      end do
    end if

    ! Do we have to compute left eigenvectors
    if (present(DleftEigenvectorsAtEdge)) then
      do idx = 1, nedges

        ! Compute velocities
        ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        
        uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute densities
        ri = DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute pressures
        pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute enthalpies
        hi = (TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        q_ij     = DCONST(0.5)*(u_ij*u_ij+v_ij*v_ij)
        cPow2_ij = max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij     = sqrt(cPow2_ij)
        aux      = a(1)*u_ij+a(2)*v_ij

        ! Compute left eigenvectors
        IDX2(DleftEigenvectorsAtEdge, 1,idx,_,_) =  DCONST(0.5)*(((HYDRO_GAMMA)-DCONST(1.0))*q_ij+c_ij*aux)/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge, 2,idx,_,_) = (cPow2_ij-((HYDRO_GAMMA)-DCONST(1.0))*q_ij)/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge, 3,idx,_,_) =  DCONST(0.5)*(((HYDRO_GAMMA)-DCONST(1.0))*q_ij-c_ij*aux)/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge, 4,idx,_,_) =  v_ij*a(1)-u_ij*a(2)

        IDX2(DleftEigenvectorsAtEdge, 5,idx,_,_) =  DCONST(0.5)*(-((HYDRO_GAMMA)-DCONST(1.0))*u_ij-c_ij*a(1))/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge, 6,idx,_,_) =  ((HYDRO_GAMMA)-DCONST(1.0))*u_ij/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge, 7,idx,_,_) =  DCONST(0.5)*(-((HYDRO_GAMMA)-DCONST(1.0))*u_ij+c_ij*a(1))/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge, 8,idx,_,_) =  a(2)

        IDX2(DleftEigenvectorsAtEdge, 9,idx,_,_) =  DCONST(0.5)*(-((HYDRO_GAMMA)-DCONST(1.0))*v_ij-c_ij*a(2))/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge,10,idx,_,_) =  ((HYDRO_GAMMA)-DCONST(1.0))*v_ij/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge,11,idx,_,_) =  DCONST(0.5)*(-((HYDRO_GAMMA)-DCONST(1.0))*v_ij+c_ij*a(2))/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge,12,idx,_,_) = -a(1)

        IDX2(DleftEigenvectorsAtEdge,13,idx,_,_) =  ((HYDRO_GAMMA)-DCONST(1.0))/DCONST(2.0)/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge,14,idx,_,_) = -((HYDRO_GAMMA)-DCONST(1.0))/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge,15,idx,_,_) =  ((HYDRO_GAMMA)-DCONST(1.0))/DCONST(2.0)/cPow2_ij
        IDX2(DleftEigenvectorsAtEdge,16,idx,_,_) =  DCONST(0.0)
      end do
    end if

  end subroutine hydro_calcCharacteristics2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTScDiss2d_sim(DdataAtEdge, DcoeffsAtEdge,&
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

    ! local variables
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj
    real(DP) :: H_ij,anorm,aux,c_ij,d_ij,q_ij,u_ij,v_ij,vel_ij
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute pressures
      pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute skew-symmetric coefficient
      a = DCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      ! Compute densities
      ri = DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      rj = DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute enthalpies
      hi = (TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
      hj = (TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(ri,rj)
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      H_ij = ROE_MEAN_VALUE(hi,hj,aux)
      
      ! Compute auxiliary variables
      vel_ij = u_ij*a(1) + v_ij*a(2)
      q_ij   = DCONST(0.5)*(u_ij*u_ij+v_ij*v_ij)

      ! Compute the speed of sound
      c_ij = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))
      
      ! Compute scalar dissipation
      d_ij = abs(vel_ij) + anorm*c_ij

      ! Compute conservative fluxes
      IDX2(DfluxesAtEdge,:,idx,_,_) = dscale*d_ij*&
          (IDX3(DdataAtEdge,:,1,idx,_,_,_)-IDX3(DdataAtEdge,:,2,idx,_,_,_))

    end do

  end subroutine hydro_calcFluxFCTScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRoeDiss2d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)


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
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj
    real(DP) :: H_ij,c2_ij,c_ij,q_ij,u_ij,v_ij,vel_ij
    real(DP) :: anorm,aux,aux1,aux2
    real(DP) :: l1,l2,l3,l4,w1,w2,w3,w4
    integer :: idx


    do idx = 1, nedges

      ! Compute the skew-symmetric coefficient and its norm
      a = DCONST(0.5)*(IDX3(DcoeffsAtEdge,:,1,idx,_,_,_)-&
                       IDX3(DcoeffsAtEdge,:,2,idx,_,_,_))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      if (anorm .gt. SYS_EPSREAL_DP) then

        ! Normalise the skew-symmetric coefficient
        a = a/anorm

        ! Compute velocities
        ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        
        uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

        ! Compute pressures
        pi = PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        pj = PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
        ! Compute densities
        ri = DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
        rj = DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
        
        ! Compute enthalpies
        hi = (TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)+pi)/ri
        hj = (TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2)
        q_ij   = DCONST(0.5)*(u_ij*u_ij+v_ij*v_ij)

        ! Compute the speed of sound
        c2_ij = max(((HYDRO_GAMMA)-DCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij  = sqrt(c2_ij)
        
        ! Compute eigenvalues
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)
        
        ! Compute solution difference U_i-U_j
        Diff = IDX3(DdataAtEdge,:,1,idx,_,_,_)-&
               IDX3(DdataAtEdge,:,2,idx,_,_,_)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = ((HYDRO_GAMMA)-DCONST(1.0))*(q_ij*Diff(1)&
                                           -u_ij*Diff(2)&
                                           -v_ij*Diff(3)&
                                                +Diff(4))/DCONST(2.0)/c2_ij

        aux2 = (vel_ij*Diff(1)&
                 -a(1)*Diff(2)&
                 -a(2)*Diff(3))/DCONST(2.0)/c_ij
        
        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((DCONST(1.0)-((HYDRO_GAMMA)-DCONST(1.0))*q_ij/c2_ij)*Diff(1)&
                                     +((HYDRO_GAMMA)-DCONST(1.0))*(u_ij*Diff(2)&
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
        IDX2(DfluxesAtEdge,:,idx,_,_) = dscale*Diff
      else
        ! Cancel conservative flux
        IDX2(DfluxesAtEdge,:,idx,_,_) = 0
      end if
    end do

  end subroutine hydro_calcFluxFCTRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRusDiss2d_sim(DdataAtEdge, DcoeffsAtEdge,&
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

    ! local variables
    real(DP) :: Ei,Ej,ci,cj,ui,uj,vi,vj
    real(DP) :: d_ij
    integer :: idx
    
    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute specific energies
      Ei = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      Ej = SPECIFICTOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(HYDRO_GAMMA)*&
          (Ei-DCONST(0.5)*(ui*ui+vi*vi)), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(HYDRO_GAMMA)*&
          (Ej-DCONST(0.5)*(uj*uj+vj*vj)), SYS_EPSREAL_DP))
      
#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,_,_,_))*uj+&
                      DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,2,idx,_,_,_))*vj)+&
                 DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)-MYNEWLINE \
                                      IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)+&
                                  POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)-MYNEWLINE \
                                      IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2))*cj,&
                  abs(DCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,_,_,_))*ui+&
                      DCONST(0.5)*(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-&
                                   IDX3(DcoeffsAtEdge,2,1,idx,_,_,_))*vi)+&
                 DCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)-MYNEWLINE \
                                      IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)+&
                                  POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)-MYNEWLINE \
                                      IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2))*ci )
#else
      ! Compute scalar dissipation
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_)*uj+&
                      IDX3(DcoeffsAtEdge,2,1,idx,_,_,_)*vj)+&
                 sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,_,_,_),2)+&
                      POW(IDX3(DcoeffsAtEdge,2,1,idx,_,_,_),2))*cj,&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_)*ui+&
                      IDX3(DcoeffsAtEdge,2,2,idx,_,_,_)*vi)+&
                 sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,_,_,_),2)+&
                      POW(IDX3(DcoeffsAtEdge,2,2,idx,_,_,_),2))*ci )
#endif
      
      ! Compute conservative flux
      IDX2(DfluxesAtEdge,:,idx,_,_) = dscale*d_ij*&
          (IDX3(DdataAtEdge,:,1,idx,_,_,_)-IDX3(DdataAtEdge,:,2,idx,_,_,_))
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
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)
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
      IDX2(DtransformedDataAtEdge,1,idx,_,_) =&
          DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          DENSITY2_2D(DdataAtNode,IDX2,idx,_,_)
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
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_)
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
      IDX2(DtransformedDataAtEdge,1,idx,_,_) =&
          TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          TOTALENERGY2_2D(DdataAtNode,IDX2,idx,_,_)
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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          ((HYDRO_GAMMA)-DCONST(1.0))*(&
            DCONST(0.5)*(ui*ui+vi*vi)*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
            ui*XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
            vi*YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)+&
            TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_))
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -((HYDRO_GAMMA)-DCONST(1.0))*(&
            DCONST(0.5)*(uj*uj+vj*vj)*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
            uj*XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
            vj*YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)+&
            TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_))
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
      IDX2(DtransformedDataAtEdge,1,idx,_,_) =&
          PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          PRESSURE2_2D(DdataAtNode,IDX2,idx,_,_)
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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

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
      IDX2(DtransformedDataAtEdge,1,idx,_,_) =&
          XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed velocity difference in y-direction
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          XVELOCITY2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed y-velocity values
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          YVELOCITY2_2D(DdataAtNode,IDX2,idx,_,_)
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
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)

      ! Transformed momentum fluxes in y-direction
      IDX3(DtransformedFluxesAtEdge,2,1,idx,_,_,_) =&
          YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,_,_,_) =&
         -YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)
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
      IDX2(DtransformedDataAtEdge,1,idx,_,_) =&
          XMOMENTUM3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          XMOMENTUM3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed momentum difference in y-direction
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          YMOMENTUM3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          YMOMENTUM3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          XMOMENTUM2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed y-momentum values
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          YMOMENTUM2_2D(DdataAtNode,IDX2,idx,_,_)
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
      IDX2(DtransformedDataAtEdge,1,idx,_,_) =&
          DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      ! Transformed total energy difference
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          TOTALENERGY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          TOTALENERGY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          DENSITY2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed energy values
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          TOTALENERGY2_2D(DdataAtNode,IDX2,idx,_,_)
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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)

      ! Transformed density fluxes
      IDX3(DtransformedFluxesAtEdge,1,1,idx,_,_,_) =&
          DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,_,_,_) =&
         -DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)
      
      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,2,1,idx,_,_,_) =&
          ((HYDRO_GAMMA)-DCONST(1.0))*(&
            DCONST(0.5)*(ui*ui+vi*vi)*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
            ui*XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
            vi*YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)+&
            TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_))
      IDX3(DtransformedFluxesAtEdge,2,2,idx,_,_,_) =&
         -((HYDRO_GAMMA)-DCONST(1.0))*(&
            DCONST(0.5)*(uj*uj+vj*vj)*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
            uj*XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
            vj*YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)+&
            TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_))
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
      IDX2(DtransformedDataAtEdge,1,idx,_,_) =&
          DENSITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          DENSITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      
      ! Transformed pressure difference
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
          PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          DENSITY2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed pressure values
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          PRESSURE2_2D(DdataAtNode,IDX2,idx,_,_)
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
      ui = XVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
      vi = YVELOCITY3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)

      uj = XVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      vj = YVELOCITY3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)
      
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

      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,4,1,idx,_,_,_) =&
          ((HYDRO_GAMMA)-DCONST(1.0))*(&
            DCONST(0.5)*(ui*ui+vi*vi)*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
            ui*XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
            vi*YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)+&
            TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_))
      IDX3(DtransformedFluxesAtEdge,4,2,idx,_,_,_) =&
         -((HYDRO_GAMMA)-DCONST(1.0))*(&
            DCONST(0.5)*(uj*uj+vj*vj)*DENSITY2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
            uj*XMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)-&
            vj*YMOMENTUM2_2D(DfluxesAtEdge,IDX2,idx,_,_)+&
            TOTALENERGY2_2D(DfluxesAtEdge,IDX2,idx,_,_))
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
      IDX2(DtransformedDataAtEdge,2,idx,_,_) =&
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

      ! Transformed pressure difference
      IDX2(DtransformedDataAtEdge,4,idx,_,_) =&
          PRESSURE3_2D(DdataAtEdge,IDX3,2,idx,_,_,_)-&
          PRESSURE3_2D(DdataAtEdge,IDX3,1,idx,_,_,_)
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
      IDX2(DtransformedDataAtNode,1,idx,_,_) =&
          DENSITY2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed x-velocity values
      IDX2(DtransformedDataAtNode,2,idx,_,_) =&
          XVELOCITY2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed y-velocity values
      IDX2(DtransformedDataAtNode,3,idx,_,_) =&
          YVELOCITY2_2D(DdataAtNode,IDX2,idx,_,_)

      ! Transformed pressure values
      IDX2(DtransformedDataAtNode,4,idx,_,_) =&
          PRESSURE2_2D(DdataAtNode,IDX2,idx,_,_)
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
      p   = ((HYDRO_GAMMA)-DCONST(1.0))*rho*(E-DCONST(0.5)*(v1*v1+v2*v2))
      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))

      ! Precompute auxiliary data
      dnxy = DbdrNormal(1)*DbdrNormal(2)
      dnx2 = DbdrNormal(1)*DbdrNormal(1)
      dny2 = DbdrNormal(2)*DbdrNormal(2)

      if (ibdrCondType .eq. BDRC_FREESLIP) then
        ! Compute reflected velocities at the boundary
        v1_b = (dny2-dnx2)*v1 - DCONST(2.0)*dnxy*v2
        v2_b = (dnx2-dny2)*v2 - DCONST(2.0)*dnxy*v1
      else
        ! Compute initial velocity from previous time step
        v1_0 = Du0(2)/Du0(1)
        v2_0 = Du0(3)/Du0(1)

        ! Compute semi-reflected velocities at the boundary
        v1_b = (dny2-dnx2)*v1-DCONST(2.0)*dnxy*v2 + 2*DbdrValue(1)*(v1_0*dnx2+v2_0*dnxy)
        v2_b = (dnx2-dny2)*v2-DCONST(2.0)*dnxy*v1 + 2*DbdrValue(1)*(v2_0*dny2+v1_0*dnxy)
      end if

      ! Compute normal velocities at the boundary and the ghost state
      ! w.r.t. the numerical/approximate  outward unit normal vector
      vn   = DpointNormal(1)*v1   + DpointNormal(2)*v2
      vn_b = DpointNormal(1)*v1_b + DpointNormal(2)*v2_b

      ! Compute the tangential velocity depending on the sign of N*v
      if (vn .gt. DCONST(0.0)) then
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
      if (DCONST(2.0)*(DCONST(2.0)/((HYDRO_GAMMA)-DCONST(1.0)))*c .le. vn_b-vn) then
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
      ppv  = p+DCONST(0.5)*(vn-vn_b)*cup
      ppv  = max(DCONST(0.0), ppv)

      if (ppv .eq. p) then

        ! Select guessed pressure from PVRS Riemann solver
        pstar = ppv
      else
        if (ppv .lt. p) then

          ! Guess pressure from the Two-Rarefaction Riemann solver
          vm    = DCONST(0.5)*(vn+vn_b)
          ptl   = DCONST(1.0) + ((HYDRO_GAMMA)-DCONST(1.0))/DCONST(2.0)*(vn-vm)/c
          ptr   = DCONST(1.0) + ((HYDRO_GAMMA)-DCONST(1.0))/DCONST(2.0)*(vm-vn_b)/c
          pstar = DCONST(0.5)*(p*ptl + p*ptr)**(DCONST(2.0)*(HYDRO_GAMMA)/((HYDRO_GAMMA)-DCONST(1.0)))
        else

          ! Guess pressure from the Two-Shock Riemann solver
          ! with PVRS as estimated pressure value
          ge    = sqrt((DCONST(2.0)/((HYDRO_GAMMA)+DCONST(1.0))/rho)/(((HYDRO_GAMMA)-&
                       DCONST(1.0))/((HYDRO_GAMMA)+DCONST(1.0))*p+ppv))
          pstar = p - DCONST(0.5)*(vn_b-vn)/ge
        end if
      end if

      ! Initialise solution difference and pressure
      vdiff = (vn_b-vn)/DCONST(2.0)
      pold  = pstar

      newton: do ite = 1, 100

        ! Compute pressure function f(pold) and its derivative f1(pold)
        if (pold .le. p) then

          ! Rarefaction wave
          prat = pold/p

          f  = (DCONST(2.0)/((HYDRO_GAMMA)-DCONST(1.0)))*c*(prat**(((HYDRO_GAMMA)-&
                DCONST(1.0))/(DCONST(2.0)*(HYDRO_GAMMA))) - DCONST(1.0))
          fd = (DCONST(1.0)/(rho*c))*prat**(-((HYDRO_GAMMA)+DCONST(1.0))/(DCONST(2.0)*(HYDRO_GAMMA)))
        else

          ! Shock wave
          auxA = DCONST(2.0)/((HYDRO_GAMMA)+DCONST(1.0))/rho
          auxB = ((HYDRO_GAMMA)-DCONST(1.0))/((HYDRO_GAMMA)+DCONST(1.0))*p
          qrt  = sqrt(auxA/(auxB + pold))

          f  = (pold-p)*qrt
          fd = (DCONST(1.0) - DCONST(0.5)*(pold - p)/(auxB + pold))*qrt
        end if

        pstar = pold - (f+vdiff)/fd
        if (pstar .lt. DCONST(0.0)) then
          pold = 1.0E-6
          cycle newton
        end if

        aux = DCONST(2.0)*abs((pstar-pold)/(pstar+pold))
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
      vn = DCONST(0.5)*(vn+vn_b)


      !-------------------------------------------------------------------------
      ! Calculate the density in the star region
      !-------------------------------------------------------------------------

      if (pstar .le. p) then

        ! Rarefaction wave
        rho = rho*(pstar/p)**(DCONST(1.0)/(HYDRO_GAMMA))
      else

        ! Shock wave
        rho = rho*(pstar/p+((HYDRO_GAMMA)-DCONST(1.0))/((HYDRO_GAMMA)+DCONST(1.0)))/&
                  (((HYDRO_GAMMA)-DCONST(1.0))/((HYDRO_GAMMA)+DCONST(1.0))*(pstar/p)+DCONST(1.0))
      end if

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DpointNormal(2)*vt+DpointNormal(1)*vn)
      Du(3) = rho*(-DpointNormal(1)*vt+DpointNormal(2)*vn)
      Du(4) = pstar/((HYDRO_GAMMA)-DCONST(1.0))+DCONST(0.5)*rho*(vn*vn+vt*vt)


    case(BDRC_VISCOUSWALL)
      !-------------------------------------------------------------------------

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      v2  = Du(3)/rho
      E   = Du(4)/rho
      p   = ((HYDRO_GAMMA)-DCONST(1.0))*rho*(E-DCONST(0.5)*(v1*v1+v2*v2))

      ! Update the solution vector and let vn:=0 and vt:=0
      Du(2) = DCONST(0.0)
      Du(3) = DCONST(0.0)
      Du(4) = p/((HYDRO_GAMMA)-DCONST(1.0))


    case(BDRC_SUPERINLET)
      !-------------------------------------------------------------------------

      ! The free stream primitive variables are Deval=[rho,v1,v2,p]
      rho = DbdrValue(1)
      p   = DbdrValue(4)
      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))
      vn  = DbdrNormal(1)*DbdrValue(2)+DbdrNormal(2)*DbdrValue(3)
      vt  = DbdrNormal(2)*DbdrValue(2)-DbdrNormal(1)*DbdrValue(3)

      ! Compute Riemann invariants based on the free stream values
      W(1) = vn-2*c/((HYDRO_GAMMA)-DCONST(1.0))
      W(2) = vn+2*c/((HYDRO_GAMMA)-DCONST(1.0))
      W(3) = p/(rho**(HYDRO_GAMMA))
      W(4) = vt

      ! Transform back into conservative variables
      vn   = DCONST(0.5)*(W(1)+W(2))
      c    = DCONST(0.25)*((HYDRO_GAMMA)-DCONST(1.0))*(W(2)-W(1))
      rho  = (c*c/(HYDRO_GAMMA)/W(3))**(DCONST(1.0)/((HYDRO_GAMMA)-DCONST(1.0)))
      p    = rho*c*c/(HYDRO_GAMMA)
      vt   = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/((HYDRO_GAMMA)-DCONST(1.0))+DCONST(0.5)*rho*(vn*vn+vt*vt)


    case(BDRC_FREESTREAM)
      !-------------------------------------------------------------------------

      ! The free stream primitive variables are Deval=[rho,v1,v2,p]
      rho = DbdrValue(1)
      p   = DbdrValue(4)
      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))
      vn  = DbdrNormal(1)*DbdrValue(2)+DbdrNormal(2)*DbdrValue(3)
      vt  = DbdrNormal(2)*DbdrValue(2)-DbdrNormal(1)*DbdrValue(3)

      ! Compute Riemann invariants based on the free stream values
      Winf(1) = vn-DCONST(2.0)*c/((HYDRO_GAMMA)-DCONST(1.0))
      Winf(2) = vn+DCONST(2.0)*c/((HYDRO_GAMMA)-DCONST(1.0))
      Winf(3) = p/(rho**(HYDRO_GAMMA))
      Winf(4) = vt

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      v2  = Du(3)/rho
      E   = Du(4)/rho
      p   = ((HYDRO_GAMMA)-DCONST(1.0))*rho*(E-DCONST(0.5)*(v1*v1+v2*v2))

      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))
      vn  = DbdrNormal(1)*v1+DbdrNormal(2)*v2
      vt  = DbdrNormal(2)*v1-DbdrNormal(1)*v2

      ! Compute Riemann invariants based on the solution values
      Wu(1) = vn-DCONST(2.0)*c/((HYDRO_GAMMA)-DCONST(1.0))
      Wu(2) = vn+DCONST(2.0)*c/((HYDRO_GAMMA)-DCONST(1.0))
      Wu(3) = p/(rho**(HYDRO_GAMMA))
      Wu(4) = vt

      ! Adopt free stream/computed values depending on the sign of the eigenvalue
      W(1) = merge(Winf(1), Wu(1), vn <  c)
      W(2) = merge(Winf(2), Wu(2), vn < -c)
      W(3) = merge(Winf(3), Wu(3), vn <  SYS_EPSREAL_DP)
      W(4) = merge(Winf(4), Wu(4), vn <  SYS_EPSREAL_DP)

      ! Transform back into conservative variables
      vn   = DCONST(0.5)*(W(1)+W(2))
      c    = DCONST(0.25)*((HYDRO_GAMMA)-DCONST(1.0))*(W(2)-W(1))
      rho  = (c*c/(HYDRO_GAMMA)/W(3))**(DCONST(1.0)/((HYDRO_GAMMA)-DCONST(1.0)))
      p    = rho*c*c/(HYDRO_GAMMA)
      vt   = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/((HYDRO_GAMMA)-DCONST(1.0))+DCONST(0.5)*rho*(vn*vn+vt*vt)


    case(BDRC_SUBINLET)
      !-------------------------------------------------------------------------

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      v2  = Du(3)/rho
      E   = Du(4)/rho
      p   = ((HYDRO_GAMMA)-DCONST(1.0))*rho*(E-DCONST(0.5)*(v1*v1+v2*v2))

      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))
      vn  = DbdrNormal(1)*v1+DbdrNormal(2)*v2
      vt  = DbdrNormal(2)*v1-DbdrNormal(1)*v2

      ! The specified density and pressure is Deval=[rho,p]
      rho = DbdrValue(1)
      p   = DbdrValue(2)

      ! Compute Riemann invariants
      W(1) = vn-DCONST(2.0)*c/((HYDRO_GAMMA)-DCONST(1.0))
      W(2) = vn+DCONST(2.0)*c/((HYDRO_GAMMA)-DCONST(1.0))
      W(3) = p/(rho**(HYDRO_GAMMA))
      W(4) = vt

      ! Transform back into conservative variables
      vn   = DCONST(0.5)*(W(1)+W(2))
      c    = DCONST(0.25)*((HYDRO_GAMMA)-DCONST(1.0))*(W(2)-W(1))
      rho  = (c*c/(HYDRO_GAMMA)/W(3))**(DCONST(1.0)/((HYDRO_GAMMA)-DCONST(1.0)))
      p    = rho*c*c/(HYDRO_GAMMA)
      vt   = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/((HYDRO_GAMMA)-DCONST(1.0))+DCONST(0.5)*rho*(vn*vn+vt*vt)


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
      p   = ((HYDRO_GAMMA)-DCONST(1.0))*rho*(E-DCONST(0.5)*(v1*v1+v2*v2))

      vn  = DbdrNormal(1)*v1+DbdrNormal(2)*v2
      vt  = DbdrNormal(2)*v1-DbdrNormal(1)*v2
      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))

      ! Compute Riemann invariants based on the solution values and prescribed exit pressure
      W(2) = 2*c/((HYDRO_GAMMA)-DCONST(1.0))-vn
      W(3) = p/(rho**(HYDRO_GAMMA))
      W(4) = vt
      W(1) = 4/((HYDRO_GAMMA)-DCONST(1.0))*sqrt(&
          max((HYDRO_GAMMA)*ps/rho*(p/ps)**(DCONST(1.0)/(HYDRO_GAMMA)), SYS_EPSREAL_DP))-W(2)

      ! Transform back into conservative variables
      vn  = DCONST(0.5)*(W(1)-W(2))
      c   = DCONST(0.25)*((HYDRO_GAMMA)-DCONST(1.0))*(W(1)+W(2))
      rho = (c*c/(HYDRO_GAMMA)/W(3))**(DCONST(1.0)/((HYDRO_GAMMA)-DCONST(1.0)))
      p   = rho*c*c/(HYDRO_GAMMA)
      vt  = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/((HYDRO_GAMMA)-DCONST(1.0))+DCONST(0.5)*rho*(vn*vn+vt*vt)

    case (BDRC_OPEN)
      !-------------------------------------------------------------------------
      ! The open boundary conditions prescribes the computed solution
      ! values at the boundary as freestream values
      
      ! That`s it
      
    case default
      call output_line('Unsupported type of boundary condition!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcBoundaryvalues2d')
      call sys_halt()
    end select

  end subroutine hydro_calcBoundaryvalues2d

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcBilfBdrCond2D(rproblemLevel, rboundaryCondition,&
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

  end subroutine hydro_calcBilfBdrCond2D

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcLinfBdrCond2D(rproblemLevel, rboundaryCondition,&
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
    integer, dimension(:), pointer :: p_IbdrCondCpIdx,p_IbdrCondType
    integer, dimension(:), pointer :: p_IbdrCompPeriodic,p_IbdrCondPeriodic
    integer :: ibct,isegment,idissipationtype
    integer(I32) :: ccubTypeBdr

    ! Evaluate linear form for boundary integral and return if
    ! there are no weak boundary conditions available
    if (.not.rboundaryCondition%bWeakBdrCond) return

    ! Check if we are in 2D
    if (rproblemLevel%rtriangulation%ndim .ne. NDIM2D) then
      call output_line('Spatial dimension must be 2D!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcLinfBdrCond2D')
      call sys_halt()
    end if

    ! Get pointer to parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)
    
    ! Get parameters from parameter list
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'ccubTypeBdr', ccubTypeBdr)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'idissipationtype', idissipationtype)

    ! Initialise temporal collection structure
    call collct_init(rcollectionTmp)

    ! Prepare quick access arrays of temporal collection structure
    rcollectionTmp%SquickAccess(1) = ''
    rcollectionTmp%SquickAccess(2) = 'rfparser'
    rcollectionTmp%DquickAccess(1) = dtime
    rcollectionTmp%DquickAccess(2) = dscale
    rcollectionTmp%IquickAccess(1) = idissipationtype
    rcollectionTmp%IquickAccess(2) = int(ccubTypeBdr)

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
        rcollectionTmp%IquickAccess(3) = p_IbdrCondType(isegment)
        rcollectionTmp%IquickAccess(4) = isegment
        rcollectionTmp%IquickAccess(5) = rboundaryCondition%nmaxExpressions
        
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
    
  end subroutine hydro_calcLinfBdrCond2D

  !***************************************************************************

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
    !   IquickAccess(1):     disipation type
    !   IquickAccess(2):     cubature rule
    !   IquickAccess(3):     boundary type
    !   IquickAccess(4):     segment number
    !   IquickAccess(5):     maximum number of expressions
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
    real(DP) :: dmaxParam,dmaxParamMirror,dminParam,dminParamMirror
    real(DP) :: dscale,dtime,dvnI,dvnM,dvtI,dvtM
    real(DP) :: eI,eM,cI,cM,hI,hM,pI,pM,rI,rM,uI,uM,vI,vM
    real(DP) :: aux,aux1,aux2,l1,l2,l3,l4,w1,w2,w3,w4
    real(DP) :: H_IM,c2_IM,c_IM,d_IM,q_IM,u_IM,v_IM,vel_IM
    integer :: ibdrtype,idissipationtype,isegment,nmaxExpr
    integer :: icubp,iel,iexpr,ipoint,ivar,ivt,neq,npoints,nvar,nve
    integer(I32) :: ccubType
    

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
    ! - the type of the dissipation
    ! - the type of boundary condition
    ! - the segment number
    ! - the maximum number of expressions
    ! - the cubature rule
    idissipationtype = rcollection%IquickAccess(1)
    ccubType         = int(rcollection%IquickAccess(2),I32)
    ibdrtype         = rcollection%IquickAccess(3)
    isegment         = rcollection%IquickAccess(4)
    nmaxExpr         = rcollection%IquickAccess(5)
    
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
        Dbas(1,icubp) = DCONST(0.5)*(DCONST(1.0)-DcubPtsRef(1,icubp))
        Dbas(2,icubp) = DCONST(0.5)*(DCONST(1.0)+DcubPtsRef(1,icubp))
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
        IDX3(DstateI,1,1,iel,_,_,_) = p_Ddata(nvar*(ivt-1)+1)
        IDX3(DstateI,2,1,iel,_,_,_) = p_Ddata(nvar*(ivt-1)+2)
        IDX3(DstateI,3,1,iel,_,_,_) = p_Ddata(nvar*(ivt-1)+3)
        IDX3(DstateI,4,1,iel,_,_,_) = p_Ddata(nvar*(ivt-1)+4)

        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)

        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store internal state vector
        IDX3(DstateI,1,2,iel,_,_,_) = p_Ddata(nvar*(ivt-1)+1)
        IDX3(DstateI,2,2,iel,_,_,_) = p_Ddata(nvar*(ivt-1)+2)
        IDX3(DstateI,3,2,iel,_,_,_) = p_Ddata(nvar*(ivt-1)+3)
        IDX3(DstateI,4,2,iel,_,_,_) = p_Ddata(nvar*(ivt-1)+4)

        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
      end do
#else
      ! Allocate temporal memory
      allocate(Daux(npointsPerElement*nvar, nelements))
      
      ! Evaluate the solution in the cubature points on the boundary
      call fevl_evaluate_sim(DER_FUNC, Daux, p_rsolution%RvectorBlock(1),&
          Dpoints, rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Distribute solution values to the internal state vector
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          IDX3(DstateI,1,ipoint,iel,_,_,_) = Daux((ipoint-1)*NVAR2D+1,iel)
          IDX3(DstateI,2,ipoint,iel,_,_,_) = Daux((ipoint-1)*NVAR2D+2,iel)
          IDX3(DstateI,3,ipoint,iel,_,_,_) = Daux((ipoint-1)*NVAR2D+3,iel)
          IDX3(DstateI,4,ipoint,iel,_,_,_) = Daux((ipoint-1)*NVAR2D+4,iel)
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
        IDX3(DstateI,1,1,iel,_,_,_) = p_Ddata(      ivt)
        IDX3(DstateI,2,1,iel,_,_,_) = p_Ddata(neq  +ivt)
        IDX3(DstateI,3,1,iel,_,_,_) = p_Ddata(neq*2+ivt)
        IDX3(DstateI,4,1,iel,_,_,_) = p_Ddata(neq*3+ivt)

        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)

        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store internal state vector
        IDX3(DstateI,1,2,iel,_,_,_) = p_Ddata(      ivt)
        IDX3(DstateI,2,2,iel,_,_,_) = p_Ddata(neq  +ivt)
        IDX3(DstateI,3,2,iel,_,_,_) = p_Ddata(neq*2+ivt)
        IDX3(DstateI,4,2,iel,_,_,_) = p_Ddata(neq*3+ivt)

        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
      end do
#else
      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements))
      
      ! Evaluate the solution in the cubature points on the boundary
      do ivar = 1, nvar
        call fevl_evaluate_sim(DER_FUNC, Daux,&
            p_rsolution%RvectorBlock(ivar), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
        ! Distribute solution values to the internal state vector
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            IDX3(DstateI,ivar,ipoint,iel,_,_,_) = Daux(ipoint,iel)
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
      
      ! Initialise values for function parser
      Dvalue = DCONST(0.0)
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
                Dvalue, IDX3(DstateM,iexpr,ipoint,iel,_,_,_))
          end do
          
          ! Compute auxiliary quantities based on free stream state vector
          rM = IDX3(DstateM,1,ipoint,iel,_,_,_)
          pM = IDX3(DstateM,4,ipoint,iel,_,_,_)
          cM = sqrt((HYDRO_GAMMA)*pM/rM)
          dvnM =  Dnx(ipoint,iel)*IDX3(DstateM,2,ipoint,iel,_,_,_)+&
                  Dny(ipoint,iel)*IDX3(DstateM,3,ipoint,iel,_,_,_)
          dvtM = -Dny(ipoint,iel)*IDX3(DstateM,2,ipoint,iel,_,_,_)+&
                  Dnx(ipoint,iel)*IDX3(DstateM,3,ipoint,iel,_,_,_)
          
          ! Compute auxiliary quantities based on internal state vector
          pI = PRESSURE3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          cI = sqrt(max((HYDRO_GAMMA)*pI/&
               DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_), SYS_EPSREAL_DP))
          
          ! Compute the normal and tangential velocities based
          ! on internal state vector
          dvnI = ( Dnx(ipoint,iel)*&
                   XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                   Dny(ipoint,iel)*&
                   YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_))/&
                   DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          dvtI = (-Dny(ipoint,iel)*&
                   XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                   Dnx(ipoint,iel)*&
                   YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_))/&
                   DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)

          ! Select free stream or computed Riemann invariant depending
          ! on the sign of the corresponding eigenvalue
          if (dvnI .lt. cI) then
            w1 = dvnM-DCONST(2.0)*cM/((HYDRO_GAMMA)-DCONST(1.0))
          else
            w1 = dvnI-DCONST(2.0)*cI/((HYDRO_GAMMA)-DCONST(1.0))
          end if
          
          if (dvnI .lt. DCONST(0.0)) then
            w2 = pM/(rM**(HYDRO_GAMMA))
            w3 = dvtM
          else
            w2 = pI/DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)**(HYDRO_GAMMA)
            w3 = dvtI
          end if
          
          if (dvnI .lt. -cI) then
            w4 = dvnM+DCONST(2.0)*cM/((HYDRO_GAMMA)-DCONST(1.0))
          else
            w4 = dvnI+DCONST(2.0)*cI/((HYDRO_GAMMA)-DCONST(1.0))
          end if
          
          ! Convert Riemann invariants into conservative state variables
          cM = DCONST(0.25)*((HYDRO_GAMMA)-DCONST(1.0))*(w4-w1)
          rM = (cM*cM/(HYDRO_GAMMA)/w2)**(DCONST(1.0)/((HYDRO_GAMMA)-DCONST(1.0)))
          pM = rM*cM*cM/(HYDRO_GAMMA)
          dvnM = DCONST(0.5)*(w1+w4)
          dvtM = w3
          
          ! Calculate the state vector based on Riemann invariants
          IDX3(DstateM,1,ipoint,iel,_,_,_) = rM
          IDX3(DstateM,2,ipoint,iel,_,_,_) = rM*(Dnx(ipoint,iel)*dvnM-Dny(ipoint,iel)*dvtM)
          IDX3(DstateM,3,ipoint,iel,_,_,_) = rM*(Dny(ipoint,iel)*dvnM+Dnx(ipoint,iel)*dvtM)
          IDX3(DstateM,4,ipoint,iel,_,_,_) = pM/((HYDRO_GAMMA)-DCONST(1.0))+DCONST(0.5)*rM*(dvnM*dvnM+dvtM*dvtM)
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
                   XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                   Dny(ipoint,iel)*&
                   YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_))/&
                   DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          dvtI = (-Dny(ipoint,iel)*&
                   XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                   Dnx(ipoint,iel)*&
                   YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_))/&
                   DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)

          ! Compute the mirrored state vector
          IDX3(DstateM,1,ipoint,iel,_,_,_) = DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          IDX3(DstateM,2,ipoint,iel,_,_,_) = DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*&
                                              (-dvnI*Dnx(ipoint,iel)-dvtI*Dny(ipoint,iel))
          IDX3(DstateM,3,ipoint,iel,_,_,_) = DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*&
                                              (-dvnI*Dny(ipoint,iel)+dvtI*Dnx(ipoint,iel))
          IDX3(DstateM,4,ipoint,iel,_,_,_) = TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
        end do
      end do
      
    case (BDRC_SUPERINLET)
      !-----------------------------------------------------------------------
      ! Supersonic inlet boundary conditions:
      !
      ! Prescribe the state vector in conservative variables
      
      ! Initialise values for function parser
      Dvalue = DCONST(0.0)
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
                Dvalue, IDX3(DstateM,iexpr,ipoint,iel,_,_,_))
          end do
          
          ! Compute convervative variables
          IDX3(DstateM,4,ipoint,iel,_,_,_) = IDX3(DstateM,4,ipoint,iel,_,_,_)/((HYDRO_GAMMA)-DCONST(1.0))&
              + IDX3(DstateM,1,ipoint,iel,_,_,_)*DCONST(0.5)*(POW(IDX3(DstateM,2,ipoint,iel,_,_,_),2)+&
                                                              POW(IDX3(DstateM,3,ipoint,iel,_,_,_),2))
          IDX3(DstateM,2,ipoint,iel,_,_,_) = IDX3(DstateM,1,ipoint,iel,_,_,_)*IDX3(DstateM,2,ipoint,iel,_,_,_)
          IDX3(DstateM,3,ipoint,iel,_,_,_) = IDX3(DstateM,1,ipoint,iel,_,_,_)*IDX3(DstateM,3,ipoint,iel,_,_,_)
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
          uI = XVELOCITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          vI = YVELOCITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          pI = PRESSURE3_2D(DstateI,IDX3,ipoint,iel,_,_,_)

          ! Calculate normal flux: ${\bf n}\cdot{\bf F}(U)$
          Dflux(1,ipoint) = Dnx(ipoint,iel)*&
                            XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)&
                          + Dny(ipoint,iel)*&
                            YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          Dflux(2,ipoint) = Dnx(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*uI+pI)&
                          + Dny(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*vI)
          Dflux(3,ipoint) = Dnx(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*uI)&
                          + Dny(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*vI+pI)
          Dflux(4,ipoint) = Dnx(ipoint,iel)*&
                            (TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+pI)*uI&
                          + Dny(ipoint,iel)*&
                            (TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+pI)*vI
        end do
        
#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Galerkin
        ! fluxes from the DOFs to the cubature points, where they are
        ! needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          DlocalData = DCONST(0.0)
          
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
      
      ! Initialise values for function parser
      Dvalue = DCONST(0.0)
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
                Dvalue, IDX3(DstateM,iexpr,ipoint,iel,_,_,_))
          end do
          
          ! Compute auxiliary quantities based on prescribed boundary values
          rM   = IDX3(DstateM,1,ipoint,iel,_,_,_)
          pM   = IDX3(DstateM,2,ipoint,iel,_,_,_)
          dvtM = IDX3(DstateM,3,ipoint,iel,_,_,_)
          cM   = sqrt((HYDRO_GAMMA)*pM/rM)
          
          ! Compute auxiliary quantities based on internal state vector
          pI = PRESSURE3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          cI = sqrt(max((HYDRO_GAMMA)*pI/&
               DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_), SYS_EPSREAL_DP))
          
          ! Compute the normal velocity based on the internal state vector
          dvnI = ( Dnx(ipoint,iel)*&
                   XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                   Dny(ipoint,iel)*&
                   YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_))/&
                   DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          
          ! Compute fourth Riemann invariant based on the internal state vector
          w4 = dvnI+DCONST(2.0)*cI/((HYDRO_GAMMA)-DCONST(1.0))
          
          ! Compute the first Riemann invariant based on the fourth Riemann
          ! invariant and the prescribed boundary values
          w1 = w4-DCONST(4.0)*cM/((HYDRO_GAMMA)-DCONST(1.0))
          
          ! Compute the normal velocity based on the first and fourth Riemann
          ! invarient
          dvnM = DCONST(0.5)*(w1+w4)
          
          ! Setup the state vector based on Rimann invariants
          IDX3(DstateM,1,ipoint,iel,_,_,_) = rM
          IDX3(DstateM,2,ipoint,iel,_,_,_) = rM*(Dnx(ipoint,iel)*dvnM-Dny(ipoint,iel)*dvtM)
          IDX3(DstateM,3,ipoint,iel,_,_,_) = rM*(Dny(ipoint,iel)*dvnM+Dnx(ipoint,iel)*dvtM)
          IDX3(DstateM,4,ipoint,iel,_,_,_) = pM/((HYDRO_GAMMA)-DCONST(1.0))+DCONST(0.5)*rM*(dvnM*dvnM+dvtM*dvtM)
        end do
      end do

    case (BDRC_SUBOUTLET)
      !-----------------------------------------------------------------------
      ! Subsonic pressure outlet boundary condition:
      !
      ! Prescribe the pressure at the outlet
      
      ! Initialise values for function parser
      Dvalue = DCONST(0.0)
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
          pI = PRESSURE3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          cI = sqrt(max((HYDRO_GAMMA)*pI/&
               DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_), SYS_EPSREAL_DP))
          
          ! Compute the normal and tangential velocities based
          ! on internal state vector
          dvnI = ( Dnx(ipoint,iel)*&
                   XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                   Dny(ipoint,iel)*&
                   YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_))/&
                   DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          dvtI = (-Dny(ipoint,iel)*&
                   XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                   Dnx(ipoint,iel)*&
                   YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_))/&
                   DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          
          ! Compute three Riemann invariants based on internal state vector
          w2 = pI/POW(DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_),HYDRO_GAMMA)
          w3 = dvtI
          w4 = dvnI+DCONST(2.0)*cI/((HYDRO_GAMMA)-DCONST(1.0))
          
          ! Compute first Riemann invariant based on fourth Riemann invariant,
          ! the computed density and pressure and the prescribed exit pressure
          w1 = w4-DCONST(4.0)/((HYDRO_GAMMA)-DCONST(1.0))*sqrt((HYDRO_GAMMA)*pM/&
               DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*(pI/pM)**(DCONST(1.0)/(HYDRO_GAMMA)))
          
          ! Convert Riemann invariants into conservative state variables
          cM = DCONST(0.25)*((HYDRO_GAMMA)-DCONST(1.0))*(w4-w1)
          rM = (cM*cM/(HYDRO_GAMMA)/w2)**(DCONST(1.0)/((HYDRO_GAMMA)-DCONST(1.0)))
          pM = rM*cM*cM/(HYDRO_GAMMA)
          dvnM = DCONST(0.5)*(w1+w4)
          dvtM = w3
          
          ! Setup the state vector based on Riemann invariants
          IDX3(DstateM,1,ipoint,iel,_,_,_) = rM
          IDX3(DstateM,2,ipoint,iel,_,_,_) = rM*(Dnx(ipoint,iel)*dvnM-Dny(ipoint,iel)*dvtM)
          IDX3(DstateM,3,ipoint,iel,_,_,_) = rM*(Dny(ipoint,iel)*dvnM+Dnx(ipoint,iel)*dvtM)
          IDX3(DstateM,4,ipoint,iel,_,_,_) = pM/((HYDRO_GAMMA)-DCONST(1.0))+DCONST(0.5)*rM*(dvnM*dvnM+dvtM*dvtM)
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
!!$        call doEvaluateAtBdrScalar(DER_FUNC, npointsPerElement*nelements*nvar,&
!!$            Daux3, p_rsolution%RvectorBlock(1), npointsPerElement*nelements,&
!!$            DpointParMirror, ibct, BDR_PAR_LENGTH, p_rboundaryRegionMirror)
!!$
!!$        do iel = 1, nelements
!!$          do ipoint = 1, npointsPerElement
!!$
!!$            ! Compute auxiliary quantities based on the internal state
!!$            ! vector evaluated on the boundary
!!$            pI = ((HYDRO_GAMMA)-DCONST(1.0))*(Daux1((ipoint-1)*NVAR2D+4,iel)-&
!!$                              DCONST(0.5)*(Daux1((ipoint-1)*NVAR2D+2,iel)**2+&
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
!!$            pM = ((HYDRO_GAMMA)-DCONST(1.0))*(Daux3((ipoint-1)*NVAR2D+4,iel)-&
!!$                              DCONST(0.5)*(Daux3((ipoint-1)*NVAR2D+2,iel)**2+&
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
!!$              DstateM(1) = dvnM-DCONST(2.0)*cM/((HYDRO_GAMMA)-DCONST(1.0))
!!$            else
!!$              DstateM(1) = dvnI-DCONST(2.0)*cI/((HYDRO_GAMMA)-DCONST(1.0))
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
!!$              DstateM(4) = dvnM+DCONST(2.0)*cM/((HYDRO_GAMMA)-DCONST(1.0))
!!$            else
!!$              DstateM(4) = dvnI+DCONST(2.0)*cI/((HYDRO_GAMMA)-DCONST(1.0))
!!$            end if
!!$
!!$            ! Convert Riemann invariants into conservative state variables
!!$            cM = DCONST(0.25)*((HYDRO_GAMMA)-DCONST(1.0))*(DstateM(4)-DstateM(1))
!!$            rM = (cM*cM/(HYDRO_GAMMA)/DstateM(2))**(DCONST(1.0)/((HYDRO_GAMMA)-DCONST(1.0)))
!!$            pM = rM*cM*cM/(HYDRO_GAMMA)
!!$            dvnM = DCONST(0.5)*(DstateM(1)+DstateM(4))
!!$            dvtM = DstateM(3)
!!$
!!$            ! Setup the state vector based on Riemann invariants
!!$            DstateM(1) = rM
!!$            DstateM(2) = rM*(Dnx(ipoint,iel)*dvnM-Dny(ipoint,iel)*dvtM)
!!$            DstateM(3) = rM*(Dny(ipoint,iel)*dvnM+Dnx(ipoint,iel)*dvtM)
!!$            DstateM(4) = pM/((HYDRO_GAMMA)-DCONST(1.0)) + DCONST(0.5)*rM*(dvnM**2+dvtM**2)
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
!!$            Dcoefficients(:,1,ipoint,iel) = dscale*DCONST(0.5)*(Dflux-Ddiff)
!!$          end do
!!$        end do
!!$
!!$        ! Deallocate temporal memory
!!$        deallocate(DpointParMirror, Daux3)

    case (BDRC_OPEN)
      !-------------------------------------------------------------------------
      ! The open boundary conditions prescribes the computed solution
      ! values at the boundary as freestream values

      ! Set the computed external state vector to the computed
      ! internal state vector. That`s it!
      DstateM = DstateI

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

    ! What type of dissipation are we?
    select case(idissipationtype)

    case (DISSIPATION_SCALAR)

      !-------------------------------------------------------------------------
      ! Solve the boundary Riemann problem at the boundary using
      ! scalar dissipation proportional to the spectral radius of the
      ! Roe matrix
      
      ! Loop over the DOFs and evaluate the Galerkin fluxes at DOFs
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Compute velocities and pressure from internal state
          uI = XVELOCITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          vI = YVELOCITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          pI = PRESSURE3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          rI = DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          hI = (TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+pI)/rI
          
          ! Compute velocities and pressure from mirrored state
          uM = XVELOCITY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          vM = YVELOCITY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          pM = PRESSURE3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          rM = DENSITY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          hM = (TOTALENERGY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)+pM)/rM
          
          ! Calculate normal flux: $\frac12{\bf n}\cdot[{\bf F}(U_I)+{\bf F}(U_M)]$
          Dflux(1,ipoint) = Dnx(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                             XMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_))&
                          + Dny(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                             YMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_))
          Dflux(2,ipoint) = Dnx(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*uI+pI+&
                             XMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*uM+pM)&
                          + Dny(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*vI+&
                             XMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*vM)
          Dflux(3,ipoint) = Dnx(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*uI+&
                             YMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*uM)&
                          + Dny(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*vI+pI+&
                             YMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*vM+pM)
          Dflux(4,ipoint) = Dnx(ipoint,iel)*&
                            ((TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+pI)*uI+&
                             (TOTALENERGY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)+pM)*uM)&
                          + Dny(ipoint,iel)*&
                            ((TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+pI)*vI+&
                             (TOTALENERGY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)+pM)*vM)

          ! Compute Roe mean values
          aux  = ROE_MEAN_RATIO(rI,rM)
          u_IM = ROE_MEAN_VALUE(uI,uM,aux)
          v_IM = ROE_MEAN_VALUE(vI,vM,aux)
          H_IM = ROE_MEAN_VALUE(hI,hM,aux)

          ! Compute auxiliary variables
          vel_IM = Dnx(ipoint,iel)*u_IM + Dny(ipoint,iel)*v_IM
          q_IM   = DCONST(0.5)*(u_IM*u_IM+v_IM*v_IM)
          
          ! Compute the speed of sound
          c_IM = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(H_IM-q_IM), SYS_EPSREAL_DP))

          ! Compute scalar dissipation
          d_IM = abs(vel_IM) + c_IM

          ! Multiply solution difference U_M-U_I by the scalar dissipation
          Ddiff(:,ipoint) = d_IM*(IDX3(DstateM,:,ipoint,iel,_,_,_)-&
                                  IDX3(DstateI,:,ipoint,iel,_,_,_))
        end do

#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Galerkin
        ! fluxes from the DOFs to the cubature points, where they are
        ! needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          DlocalData = DCONST(0.0)
          ! Loop over the DOFs and interpolate the Galerkin fluxes
          do ipoint = 1, npoints
            DlocalData = DlocalData&
                + Dbas(ipoint,icubp)*DCONST(0.5)*(Dflux(:,ipoint)-Ddiff(:,ipoint))
          end do
          
          ! Store flux in the cubature points
          Dcoefficients(:,1,icubp,iel) = dscale*DlocalData
        end do
#else
        ! Loop over the cubature points and store the fluxes
        do ipoint = 1, npointsPerElement
          Dcoefficients(:,1,ipoint,iel) =&
              dscale*DCONST(0.5)*(Dflux(:,ipoint)-Ddiff(:,ipoint))
        end do
#endif
      end do

    case (DISSIPATION_SCALAR_DSPLIT)

      !-------------------------------------------------------------------------
      ! Solve the boundary Riemann problem at the boundary using
      ! scalar dissipation proportional to the spectral radius of the
      ! Roe matrix, whereby dimensional splitting is employed
      
      ! Loop over the DOFs and evaluate the Galerkin fluxes at DOFs
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Compute velocities and pressure from internal state
          uI = XVELOCITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          vI = YVELOCITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          pI = PRESSURE3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          rI = DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          hI = (TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+pI)/rI
          
          ! Compute velocities and pressure from mirrored state
          uM = XVELOCITY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          vM = YVELOCITY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          pM = PRESSURE3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          rM = DENSITY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          hM = (TOTALENERGY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)+pM)/rM
          
          ! Calculate normal flux: $\frac12{\bf n}\cdot[{\bf F}(U_I)+{\bf F}(U_M)]$
          Dflux(1,ipoint) = Dnx(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                             XMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_))&
                          + Dny(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                             YMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_))
          Dflux(2,ipoint) = Dnx(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*uI+pI+&
                             XMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*uM+pM)&
                          + Dny(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*vI+&
                             XMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*vM)
          Dflux(3,ipoint) = Dnx(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*uI+&
                             YMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*uM)&
                          + Dny(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*vI+pI+&
                             YMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*vM+pM)
          Dflux(4,ipoint) = Dnx(ipoint,iel)*&
                            ((TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+pI)*uI+&
                             (TOTALENERGY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)+pM)*uM)&
                          + Dny(ipoint,iel)*&
                            ((TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+pI)*vI+&
                             (TOTALENERGY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)+pM)*vM)

          ! Compute Roe mean values
          aux  = ROE_MEAN_RATIO(rI,rM)
          u_IM = ROE_MEAN_VALUE(uI,uM,aux)
          v_IM = ROE_MEAN_VALUE(vI,vM,aux)
          H_IM = ROE_MEAN_VALUE(hI,hM,aux)

          ! Compute auxiliary variables
          q_IM   = DCONST(0.5)*(u_IM*u_IM+v_IM*v_IM)
          
          ! Compute the speed of sound
          c_IM = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(H_IM-q_IM), SYS_EPSREAL_DP))

          ! Compute scalar dissipation with dimensional splitting
          d_IM = ( abs(Dnx(ipoint,iel)*u_IM) + abs(Dnx(ipoint,iel))*c_IM +&
                   abs(Dny(ipoint,iel)*v_IM) + abs(Dny(ipoint,iel))*c_IM )

          ! Multiply solution difference U_M-U_I by the scalar dissipation
          Ddiff(:,ipoint) = d_IM*(IDX3(DstateM,:,ipoint,iel,_,_,_)-&
                                  IDX3(DstateI,:,ipoint,iel,_,_,_))
        end do

#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Galerkin
        ! fluxes from the DOFs to the cubature points, where they are
        ! needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          DlocalData = DCONST(0.0)
          ! Loop over the DOFs and interpolate the Galerkin fluxes
          do ipoint = 1, npoints
            DlocalData = DlocalData&
                + Dbas(ipoint,icubp)*DCONST(0.5)*(Dflux(:,ipoint)-Ddiff(:,ipoint))
          end do
          
          ! Store flux in the cubature points
          Dcoefficients(:,1,icubp,iel) = dscale*DlocalData
        end do
#else
        ! Loop over the cubature points and store the fluxes
        do ipoint = 1, npointsPerElement
          Dcoefficients(:,1,ipoint,iel) =&
              dscale*DCONST(0.5)*(Dflux(:,ipoint)-Ddiff(:,ipoint))
        end do
#endif
      end do

    case (DISSIPATION_ROE)
      !-------------------------------------------------------------------------
      ! Solve the boundary Riemann problem at the boundary using
      ! tensorial dissipation of Roe-type
      
      ! Loop over the DOFs and evaluate the Galerkin fluxes at DOFs
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Compute velocities and pressure from internal state
          uI = XVELOCITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          vI = YVELOCITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          pI = PRESSURE3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          rI = DENSITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          hI = (TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+pI)/rI
          
          ! Compute velocities and pressure from mirrored state
          uM = XVELOCITY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          vM = YVELOCITY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          pM = PRESSURE3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          rM = DENSITY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          hM = (TOTALENERGY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)+pM)/rM
          
          ! Calculate normal flux: $\frac12{\bf n}\cdot[{\bf F}(U_I)+{\bf F}(U_M)]$
          Dflux(1,ipoint) = Dnx(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                             XMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_))&
                          + Dny(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                             YMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_))
          Dflux(2,ipoint) = Dnx(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*uI+pI+&
                             XMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*uM+pM)&
                          + Dny(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*vI+&
                             XMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*vM)
          Dflux(3,ipoint) = Dnx(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*uI+&
                             YMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*uM)&
                          + Dny(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*vI+pI+&
                             YMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*vM+pM)
          Dflux(4,ipoint) = Dnx(ipoint,iel)*&
                            ((TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+pI)*uI+&
                             (TOTALENERGY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)+pM)*uM)&
                          + Dny(ipoint,iel)*&
                            ((TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+pI)*vI+&
                             (TOTALENERGY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)+pM)*vM)

          ! Compute Roe mean values
          aux  = ROE_MEAN_RATIO(rI,rM)
          u_IM = ROE_MEAN_VALUE(uI,uM,aux)
          v_IM = ROE_MEAN_VALUE(vI,vM,aux)
          H_IM = ROE_MEAN_VALUE(hI,hM,aux)
          
          ! Compute auxiliary variables
          vel_IM = Dnx(ipoint,iel)*u_IM + Dny(ipoint,iel)*v_IM
          q_IM   = DCONST(0.5)*(u_IM*u_IM+v_IM*v_IM)
          
          ! Compute the speed of sound
          c2_IM = max(((HYDRO_GAMMA)-DCONST(1.0))*(H_IM-q_IM), SYS_EPSREAL_DP)
          c_IM  = sqrt(c2_IM)
          
          ! Compute eigenvalues
          l1 = abs(vel_IM-c_IM)
          l2 = abs(vel_IM)
          l3 = abs(vel_IM+c_IM)
          l4 = abs(vel_IM)
          
          ! Compute solution difference U_M-U_I
          Ddiff(:,ipoint) = IDX3(DstateM,:,ipoint,iel,_,_,_)-&
                            IDX3(DstateI,:,ipoint,iel,_,_,_)
          
          ! Compute auxiliary quantities for characteristic variables
          aux1 = ((HYDRO_GAMMA)-DCONST(1.0))*(q_IM*Ddiff(1,ipoint)&
                                             -u_IM*Ddiff(2,ipoint)&
                                             -v_IM*Ddiff(3,ipoint)&
                                                  +Ddiff(4,ipoint))/DCONST(2.0)/c2_IM
          aux2 =        (vel_IM*Ddiff(1,ipoint)&
               -Dnx(ipoint,iel)*Ddiff(2,ipoint)&
               -Dny(ipoint,iel)*Ddiff(3,ipoint))/DCONST(2.0)/c_IM
          
          ! Compute characteristic variables multiplied by the
          ! corresponding eigenvalue
          w1 = l1 * (aux1 + aux2)
          w2 = l2 * ((DCONST(1.0)-((HYDRO_GAMMA)-DCONST(1.0))*q_IM/c2_IM)*Ddiff(1,ipoint)&
                                       +((HYDRO_GAMMA)-DCONST(1.0))*(u_IM*Ddiff(2,ipoint)&
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
          DlocalData = DCONST(0.0)
          ! Loop over the DOFs and interpolate the Galerkin fluxes
          do ipoint = 1, npoints
            DlocalData = DlocalData&
                + Dbas(ipoint,icubp)*DCONST(0.5)*(Dflux(:,ipoint)-Ddiff(:,ipoint))
          end do
          
          ! Store flux in the cubature points
          Dcoefficients(:,1,icubp,iel) = dscale*DlocalData
        end do
#else
        ! Loop over the cubature points and store the fluxes
        do ipoint = 1, npointsPerElement
          Dcoefficients(:,1,ipoint,iel) =&
              dscale*DCONST(0.5)*(Dflux(:,ipoint)-Ddiff(:,ipoint))
        end do
#endif
      end do

    case (DISSIPATION_RUSANOV)

      !-------------------------------------------------------------------------
      ! Solve the boundary Riemann problem at the boundary using
      ! scalar dissipation of Rusanov-type
      
      ! Loop over the DOFs and evaluate the Galerkin fluxes at DOFs
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Compute velocities and pressure from internal state
          uI = XVELOCITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          vI = YVELOCITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          pI = PRESSURE3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          
          ! Compute velocities and pressure from mirrored state
          uM = XVELOCITY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          vM = YVELOCITY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          pM = PRESSURE3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          
          ! Calculate normal flux: $\frac12{\bf n}\cdot[{\bf F}(U_I)+{\bf F}(U_M)]$
          Dflux(1,ipoint) = Dnx(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                             XMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_))&
                          + Dny(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                             YMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_))
          Dflux(2,ipoint) = Dnx(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*uI+pI+&
                             XMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*uM+pM)&
                          + Dny(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*vI+&
                             XMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*vM)
          Dflux(3,ipoint) = Dnx(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*uI+&
                             YMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*uM)&
                          + Dny(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*vI+pI+&
                             YMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*vM+pM)
          Dflux(4,ipoint) = Dnx(ipoint,iel)*&
                            ((TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+pI)*uI+&
                             (TOTALENERGY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)+pM)*uM)&
                          + Dny(ipoint,iel)*&
                            ((TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+pI)*vI+&
                             (TOTALENERGY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)+pM)*vM)

          ! Compute specific energies
          eI = SPECIFICTOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          eM = SPECIFICTOTALENERGY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          
          ! Compute the speed of sound
          cI = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(HYDRO_GAMMA)*&
              (eI-DCONST(0.5)*(uI*uI+vI*vI)), SYS_EPSREAL_DP))
          cM = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(HYDRO_GAMMA)*&
              (eM-DCONST(0.5)*(uM*uM+vM*vM)), SYS_EPSREAL_DP))
          
          ! Compute scalar dissipation
          d_IM = max( abs(Dnx(ipoint,iel)*uI) + abs(Dny(ipoint,iel)*vI) + cI,&
                      abs(Dnx(ipoint,iel)*uM) + abs(Dny(ipoint,iel)*vM) + cM )

          ! Multiply solution difference U_M-U_I by the scalar dissipation
          Ddiff(:,ipoint) = d_IM*(IDX3(DstateM,:,ipoint,iel,_,_,_)-&
                                  IDX3(DstateI,:,ipoint,iel,_,_,_))
        end do

#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Galerkin
        ! fluxes from the DOFs to the cubature points, where they are
        ! needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          DlocalData = DCONST(0.0)
          ! Loop over the DOFs and interpolate the Galerkin fluxes
          do ipoint = 1, npoints
            DlocalData = DlocalData&
                + Dbas(ipoint,icubp)*DCONST(0.5)*(Dflux(:,ipoint)-Ddiff(:,ipoint))
          end do
          
          ! Store flux in the cubature points
          Dcoefficients(:,1,icubp,iel) = dscale*DlocalData
        end do
#else
        ! Loop over the cubature points and store the fluxes
        do ipoint = 1, npointsPerElement
          Dcoefficients(:,1,ipoint,iel) =&
              dscale*DCONST(0.5)*(Dflux(:,ipoint)-Ddiff(:,ipoint))
        end do
#endif
      end do

    case (DISSIPATION_RUSANOV_DSPLIT)

      !-------------------------------------------------------------------------
      ! Solve the boundary Riemann problem at the boundary using
      ! scalar dissipation of Rusanov-type, whereby dimensional
      ! splitting is employed
      
      ! Loop over the DOFs and evaluate the Galerkin fluxes at DOFs
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Compute velocities and pressure from internal state
          uI = XVELOCITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          vI = YVELOCITY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          pI = PRESSURE3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          
          ! Compute velocities and pressure from mirrored state
          uM = XVELOCITY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          vM = YVELOCITY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          pM = PRESSURE3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          
          ! Calculate normal flux: $\frac12{\bf n}\cdot[{\bf F}(U_I)+{\bf F}(U_M)]$
          Dflux(1,ipoint) = Dnx(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                             XMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_))&
                          + Dny(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+&
                             YMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_))
          Dflux(2,ipoint) = Dnx(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*uI+pI+&
                             XMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*uM+pM)&
                          + Dny(ipoint,iel)*&
                            (XMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*vI+&
                             XMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*vM)
          Dflux(3,ipoint) = Dnx(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*uI+&
                             YMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*uM)&
                          + Dny(ipoint,iel)*&
                            (YMOMENTUM3_2D(DstateI,IDX3,ipoint,iel,_,_,_)*vI+pI+&
                             YMOMENTUM3_2D(DstateM,IDX3,ipoint,iel,_,_,_)*vM+pM)
          Dflux(4,ipoint) = Dnx(ipoint,iel)*&
                            ((TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+pI)*uI+&
                             (TOTALENERGY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)+pM)*uM)&
                          + Dny(ipoint,iel)*&
                            ((TOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)+pI)*vI+&
                             (TOTALENERGY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)+pM)*vM)

          ! Compute specific energies
          eI = SPECIFICTOTALENERGY3_2D(DstateI,IDX3,ipoint,iel,_,_,_)
          eM = SPECIFICTOTALENERGY3_2D(DstateM,IDX3,ipoint,iel,_,_,_)
          
          ! Compute the speed of sound
          cI = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(HYDRO_GAMMA)*&
              (eI-DCONST(0.5)*(uI*uI+vI*vI)), SYS_EPSREAL_DP))
          cM = sqrt(max(((HYDRO_GAMMA)-DCONST(1.0))*(HYDRO_GAMMA)*&
              (eM-DCONST(0.5)*(uM*uM+vM*vM)), SYS_EPSREAL_DP))
          
          ! Compute scalar dissipation
          d_IM = max( abs(Dnx(ipoint,iel)*uI) + abs(Dnx(ipoint,iel))*cI,&
                      abs(Dnx(ipoint,iel)*uM) + abs(Dnx(ipoint,iel))*cM )&
               + max( abs(Dny(ipoint,iel)*vI) + abs(Dny(ipoint,iel))*cI,&
                      abs(Dny(ipoint,iel)*vM) + abs(Dny(ipoint,iel))*cM )

          ! Multiply solution difference U_M-U_I by the scalar dissipation
          Ddiff(:,ipoint) = d_IM*(IDX3(DstateM,:,ipoint,iel,_,_,_)-&
                                  IDX3(DstateI,:,ipoint,iel,_,_,_))
        end do

#ifdef HYDRO_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Galerkin
        ! fluxes from the DOFs to the cubature points, where they are
        ! needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          DlocalData = DCONST(0.0)
          ! Loop over the DOFs and interpolate the Galerkin fluxes
          do ipoint = 1, npoints
            DlocalData = DlocalData&
                + Dbas(ipoint,icubp)*DCONST(0.5)*(Dflux(:,ipoint)-Ddiff(:,ipoint))
          end do
          
          ! Store flux in the cubature points
          Dcoefficients(:,1,icubp,iel) = dscale*DlocalData
        end do
#else
        ! Loop over the cubature points and store the fluxes
        do ipoint = 1, npointsPerElement
          Dcoefficients(:,1,ipoint,iel) =&
              dscale*DCONST(0.5)*(Dflux(:,ipoint)-Ddiff(:,ipoint))
        end do
#endif
      end do

    case default
      call output_line('Invalid type of dissipation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_coeffVectorBdr2d_sim')
      call sys_halt()
    end select

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
    !***************************************************************************

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
    !***************************************************************************

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

end module hydro_callback2d
