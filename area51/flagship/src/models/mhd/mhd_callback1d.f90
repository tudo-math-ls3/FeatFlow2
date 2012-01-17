!##############################################################################
!# ****************************************************************************
!# <name> mhd_callback1d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible MHDequations in 1D.
!#
!# The following callback functions are available:
!#
!# 1.) mhd_calcFluxGal1d_sim
!#     -> Computes fluxes for standard Galerkin scheme
!#
!# 2.) mhd_calcFluxGalNoBdr1d_sim
!#     -> Computes fluxes for standard Galerkin scheme without
!#        assembling the symmetric boundary contribution
!#
!# 3.) mhd_calcFluxScDiss1d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting scalar artificial viscosities
!#
!# 4.) mhd_calcFluxRoeDiss1d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities
!#
!# 5.) mhd_calcFluxRusDiss1d_sim
!#     -> Computes fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion
!#
!# 6.) mhd_calcMatDiagMatD1d_sim
!#     -> Computes local matrix for diagonal entry
!#
!# 7.) mhd_calcMatDiag1d_sim
!#     -> Computes local matrix for diagonal entry
!#
!# 8.) mhd_calcMatGalMatD1d_sim
!#     -> Computes local matrices for standard Galerkin scheme
!#
!# 9.) mhd_calcMatGal1d_sim
!#     -> Computes local matrices for standard Galerkin scheme
!#
!# 10.) mhd_calcMatScDissMatD1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 11.) mhd_calcMatScDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 12.) mhd_calcMatRoeDissMatD1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities of Roe-type
!#
!# 13.) mhd_calcMatRoeDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities of Roe-type
!#
!# 14.) mhd_calcMatRusDissMatD1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities of Rusanov-type
!#
!# 15.) mhd_calcMatRusDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities of Rusanov-type
!#
!# 16.) mhd_calcCharacteristics1d_sim
!#      -> Computes characteristic variables
!#
!# 17.) mhd_calcFluxFCTScDiss1d_sim
!#      -> Computes fluxes for FCT algorithm adopting scalar
!#         artificial viscosities
!#
!# 18.) mhd_calcFluxFCTRoeDiss1d_sim
!#      -> Computes fluxes for FCT algorithm adopting tensorial
!#         artificial viscosities of Roe-type
!#
!# 19.) mhd_calcFluxFCTRusDiss1d_sim
!#      -> Computes fluxes for FCT algorithm adopting scalar
!#         artificial viscosities of Rusanov-type
!#
!# 20.) mhd_trafoFluxDensity1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 21.) mhd_trafoDiffDensity1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 22.) mhd_trafoNodalDensity1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density values
!#
!# 23.) mhd_trafoFluxEnergy1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 24.) mhd_trafoDiffEnergy1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 25.) mhd_trafoNodalEnergy1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal energy values
!#
!# 26.) mhd_trafoFluxPressure1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 27.) mhd_trafoDiffPressure1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the pressure
!#
!# 28.) mhd_trafoNodalPressure1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal pressure values
!#
!# 29.) mhd_trafoFluxVelocity1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 30.) mhd_trafoDiffVelocity1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 31.) mhd_trafoNodalVelocity1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal velocity values
!#
!# 32.) mhd_trafoFluxMomentum1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 33.) mhd_trafoDiffMomentum1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 34.) mhd_trafoNodalMomentum1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal momentum values
!#
!# 35.) mhd_trafoFluxDenEng1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 36.) mhd_trafoDiffDenEng1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 37.) mhd_trafoNodalDenEng1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density and energy values
!#
!# 38.) mhd_trafoFluxDenPre1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 39.) mhd_trafoDiffDenPre1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 40.) mhd_trafoNodalDenPre1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density and pressure values
!#
!# 41.) mhd_trafoFluxDenPreVel1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 42.) mhd_trafoDiffDenPreVel1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure
!#         and the velocity
!#
!# 43.) mhd_trafoNodalDenPreVel1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density, pressure and velocity values
!#
!# 44.) mhd_trafoDiffMagfield1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the magnetic field
!#
!# 45.) mhd_trafoFluxMagfield1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the magnetic field
!#
!# 46.) mhd_trafoNodalMagfield1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal values for the magnetic field
!#
!# 47.) mhd_calcBoundaryvalues1d
!#      -> Computes the boundary values for a given node
!#
!# 48.) mhd_calcBilfBdrCond1d
!#      -> Calculates the bilinear form arising from the weak
!#         imposition of boundary conditions in 1D
!#
!# 49.) mhd_calcLinfBdrCond1d
!#      -> Calculates the linear form arising from the weak
!#         imposition of boundary conditions in 1D
!#
!# 50.) mhd_coeffVectorBdr1d_sim
!#      -> Calculates the coefficients for the linear form in 1D
!#
!# </purpose>
!##############################################################################

module mhd_callback1d

#define MHD_NDIM 1
#include "mhd.h"

!$use omp_lib
  use basicgeometry
  use boundarycondaux
  use collection
  use derivatives
  use domainintegration
  use feevaluation
  use fparser
  use fsystem
  use genoutput
  use graph
  use groupfemsystem
  use linearformevaluation
  use linearsystemblock
  use linearsystemscalar
  use problem
  use scalarpde
  use solveraux
  use spatialdiscretisation
  use storage

  ! Modules from MHD model
  use mhd_basic

  implicit none

  private
  public :: mhd_calcFluxGal1d_sim
  public :: mhd_calcFluxGalNoBdr1d_sim
  public :: mhd_calcFluxScDiss1d_sim
  public :: mhd_calcFluxRoeDiss1d_sim
  public :: mhd_calcFluxRusDiss1d_sim
  public :: mhd_calcMatDiagMatD1d_sim
  public :: mhd_calcMatDiag1d_sim
  public :: mhd_calcMatGalMatD1d_sim
  public :: mhd_calcMatGal1d_sim
  public :: mhd_calcMatScDissMatD1d_sim
  public :: mhd_calcMatScDiss1d_sim
  public :: mhd_calcMatRoeDissMatD1d_sim
  public :: mhd_calcMatRoeDiss1d_sim
  public :: mhd_calcMatRusDissMatD1d_sim
  public :: mhd_calcMatRusDiss1d_sim
  public :: mhd_calcCharacteristics1d_sim
  public :: mhd_calcFluxFCTScDiss1d_sim
  public :: mhd_calcFluxFCTRoeDiss1d_sim
  public :: mhd_calcFluxFCTRusDiss1d_sim
  public :: mhd_trafoFluxDensity1d_sim
  public :: mhd_trafoFluxEnergy1d_sim
  public :: mhd_trafoFluxPressure1d_sim
  public :: mhd_trafoFluxVelocity1d_sim
  public :: mhd_trafoFluxMomentum1d_sim
  public :: mhd_trafoFluxDenEng1d_sim
  public :: mhd_trafoFluxDenPre1d_sim
  public :: mhd_trafoFluxDenPreVel1d_sim
  public :: mhd_trafoFluxMagfield1d_sim
  public :: mhd_trafoDiffDensity1d_sim
  public :: mhd_trafoDiffEnergy1d_sim
  public :: mhd_trafoDiffPressure1d_sim
  public :: mhd_trafoDiffVelocity1d_sim
  public :: mhd_trafoDiffMomentum1d_sim
  public :: mhd_trafoDiffDenEng1d_sim
  public :: mhd_trafoDiffDenPre1d_sim
  public :: mhd_trafoDiffDenPreVel1d_sim
  public :: mhd_trafoDiffMagfield1d_sim
  public :: mhd_trafoNodalDensity1d_sim
  public :: mhd_trafoNodalEnergy1d_sim
  public :: mhd_trafoNodalPressure1d_sim
  public :: mhd_trafoNodalVelocity1d_sim
  public :: mhd_trafoNodalMomentum1d_sim
  public :: mhd_trafoNodalDenEng1d_sim
  public :: mhd_trafoNodalDenPre1d_sim
  public :: mhd_trafoNodalDenPreVel1d_sim
  public :: mhd_trafoNodalMagfield1d_sim
  public :: mhd_calcBoundaryvalues1d
  public :: mhd_calcBilfBdrCond1d
  public :: mhd_calcLinfBdrCond1d
  public :: mhd_coeffVectorBdr1d_sim

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGal1d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the standard Galerkin
    ! discretisation in 1D.
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
    real(DP), dimension(NVAR1D) :: Fi,Fj
#else
    real(DP), dimension(NVAR1D) :: F_ij
#endif
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx


    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
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
      Fi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)

      ! Assemble skew-symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fj-&
           IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fi )
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
      ! Compute flux difference for x-direction
      F_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
                 
      ! Assemble fluxes
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale * IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*F_ij
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale * IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*F_ij
#endif
    end do

  end subroutine mhd_calcFluxGal1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGalNoBdr1d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the TVD discretisation
    ! in 1D. The symmetric boundary contributions are neglected and
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
    real(DP), dimension(NVAR1D) :: F_ij
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx


    do idx = 1, nedges

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
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
      F_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      
      ! Assemble symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                       IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*F_ij
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
    end do

  end subroutine mhd_calcFluxGalNoBdr1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxScDiss1d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 1D using scalar artificial viscosities proportional to the
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
    real(DP), dimension(NVAR1D) :: Fi,Fj
#else
    real(DP), dimension(NVAR1D) :: F_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: hi,hj,pi,pj,qi,qj,ri,rj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,X_ij,aux,cf_ij,d_ij,q_ij,rho_ij,u_ij
    real(DP) :: aPow2_ij,astPow2_ij,bPow2_ij,bxPow2_ij
    integer :: idx


    do idx = 1, nedges
          
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
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
      Fi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      F_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar artificial dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !
      ! There are seven eigenvalues
      !
      !   u-cf, u-ca, u-cs, u, u+cs, u+ca, u+cf,
      !
      ! where u is the x-velocity component and ca, cs and cf are the
      ! velocities of the Alfveen waves, the slow and fast waves.
      !
      ! The largest in magnitude eigenvalue is |u|+cf
      !-------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                       IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))

      ! Compute densities
      ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute enthalpies
      hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
      hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(ri,rj)
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      H_ij = ROE_MEAN_VALUE(hi,hj,aux)

      ! Compute the Roe-averaged density with left and right states interchanged!
      rho_ij = ROE_MEAN_VALUE(rj,ri,aux)
   
      ! Compute the square of the Roe-averaged speed of the Alfven waves.
      ! Note that left and right states are interchanged!
      bxPow2_ij = (POW(ROE_MEAN_VALUE(\
                       XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                       XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2))/rho_ij
    
      ! Compute the density-averaged magnetic field
      X_ij = (POW(XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)-XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),2)+&
              POW(YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)-YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),2)+&
              POW(ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)-ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),2))/&
              RCONST(2.0)*POW(sqrt(ri)+sqrt(rj),2)

      ! Compute the square of the Roe-averaged magnetic field.
      ! Note that left and right states are interchanged!
      bPow2_ij = (POW(ROE_MEAN_VALUE(\
                      XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                      XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2)+\
                  POW(ROE_MEAN_VALUE(\
                      YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                      YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2)+\
                  POW(ROE_MEAN_VALUE(\
                      ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                      ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2))/rho_ij

      ! Compute the magnitude of the Roe-averaged velocity
      q_ij = (POW(ROE_MEAN_VALUE(ui,uj,aux),2)+&
              POW(ROE_MEAN_VALUE(vi,vj,aux),2)+&
              POW(ROE_MEAN_VALUE(wi,wj,aux),2))*RCONST(0.5)

      ! Compute the Roe-averaged speed of sound
      aPow2_ij = (RCONST(2.0)-(MAGNETOHYDRODYN_GAMMA))*X_ij +&
          (MAGNETOHYDRODYN_GAMMA-RCONST(1.0))*(H_ij-q_ij-bPow2_ij)
      
      ! Compute auxiliary variables
      astPow2_ij = aPow2_ij+bPow2_ij
      aux        = sqrt(POW(astPow2_ij,2)-RCONST(4.0)*aPow2_ij*bxPow2_ij)
            
      ! Compute the Roe-averagred speed of the fast waves
      cf_ij = sqrt(RCONST(0.5)*(astPow2_ij+aux))

      ! Scalar dissipation
      d_ij = abs(u_ij*a(1)) + abs(a(1))*cf_ij
      
      ! Multiply the solution difference by the scalar dissipation
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
                   IDX3(DdataAtEdge,:,1,idx,0,0,0))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------
      
#ifdef MHD_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fj-&
           IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*F_ij + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*F_ij + Diff)
#endif
    end do

  end subroutine mhd_calcFluxScDiss1d_sim

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcFluxRoeDiss1d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 1D using tensorial artificial viscosities of Roe-type.
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
    real(DP), dimension(NVAR1D) :: Fi,Fj
#else
    real(DP), dimension(NVAR1D) :: F_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP), dimension(7,7) :: Reig
    real(DP) :: hi,hj,pi,pj,qi,qj,ri,rj,ui,uj,vi,vj,wi,wj
    real(DP) :: caPow2_ij,ca_ij,cfPow2_ij,cf_ij,csPow2_ij,cs_ij
    real(DP) :: S,aux,auxf,auxs,auxsqr,auxy,auxz,anorm
    real(DP) :: H_ij,X_ij,q_ij,rho_ij,u_ij,v_ij,w_ij
    real(DP) :: aPow2_ij,astPow2_ij,bPow2_ij,bxPow2_ij
    real(DP) :: l1,l2,l3,l4,l5,l6,l7,w1,w2,w3,w4,w5,w6,w7
    integer, dimension(7) :: Ipiv
    integer :: idx,info

    
    do idx = 1, nedges
          
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
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
      Fi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      F_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !
      ! There are seven eigenvalues
      !
      !   u-cf, u-ca, u-cs, u, u+cs, u+ca, u+cf,
      !
      ! where u is the x-velocity component and ca, cs and cf are the
      ! velocities of the Alfveen waves, the slow and fast waves.
      !-------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                       IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))
      anorm = abs(a(1)) ! = sqrt(a(1)*a(1))

      if (anorm .gt. SYS_EPSREAL_DP) then

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
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)

        ! Compute the Roe-averaged density with left and right states interchanged!
        rho_ij = ROE_MEAN_VALUE(rj,ri,aux)
        
        ! Compute the square of the Roe-averaged speed of the Alfven waves.
        ! Note that left and right states are interchanged!
        bxPow2_ij = (POW(ROE_MEAN_VALUE(\
                       XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                       XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2))/rho_ij
        ca_ij = sqrt(bxPow2_ij)
    
        ! Compute the density-averaged magnetic field
        X_ij = (POW(XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)-XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),2)+&
                POW(YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)-YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),2)+&
                POW(ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)-ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),2))/&
                RCONST(2.0)*POW(sqrt(ri)+sqrt(rj),2)

        ! Compute the square of the Roe-averaged magnetic field.
        ! Note that left and right states are interchanged!
        bPow2_ij = (POW(ROE_MEAN_VALUE(\
                      XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                      XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2)+\
                  POW(ROE_MEAN_VALUE(\
                      YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                      YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2)+\
                  POW(ROE_MEAN_VALUE(\
                      ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                      ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2))/rho_ij

        ! Compute the magnitude of the Roe-averaged velocity
        q_ij = (POW(ROE_MEAN_VALUE(ui,uj,aux),2)+&
              POW(ROE_MEAN_VALUE(vi,vj,aux),2)+&
              POW(ROE_MEAN_VALUE(wi,wj,aux),2))*RCONST(0.5)

        ! Compute the Roe-averaged speed of sound
        aPow2_ij = (RCONST(2.0)-(MAGNETOHYDRODYN_GAMMA))*X_ij +&
            (MAGNETOHYDRODYN_GAMMA-RCONST(1.0))*(H_ij-q_ij-bPow2_ij)

        ! Compute auxiliary quantities
        astPow2_ij = aPow2_ij+bPow2_ij
        auxsqr     = sqrt(POW(astPow2_ij,2)-RCONST(4.0)*aPow2_ij*bxPow2_ij)

        ! Compute the Roe-averagred speed of the slow and fast waves
        cfPow2_ij = RCONST(0.5)*(astPow2_ij+auxsqr); cf_ij=sqrt(cfPow2_ij)
        csPow2_ij = RCONST(0.5)*(astPow2_ij-auxsqr); cs_ij=sqrt(csPow2_ij)

        ! Compute eigenvalues
        l1 = abs(u_ij-cf_ij)
        l2 = abs(u_ij-ca_ij)
        l3 = abs(u_ij-cs_ij)
        l4 = abs(u_ij)
        l5 = abs(u_ij+cs_ij)
        l6 = abs(u_ij+ca_ij)
        l7 = abs(u_ij+cf_ij)
        
        ! Compute solution difference U_j-U_i
        Diff = IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
               IDX3(DdataAtEdge,:,1,idx,0,0,0)

        ! The variable names are adopted from J.M. Stone, T.A. Gardiner,
        ! P. Teuben, J.F. Hawlay and J.B. Simon ATHENA: A New Code for
        ! Astrophysical MHD. The Astrophysocal Journal Supplement Series
        ! (ISSN 0067-0049), vol. 178, September 2008, p. 137-177.

        ! Compute the "alpha_f,s" values
        auxf = sqrt((aPow2_ij-csPow2_ij)/(cfPow2_ij-csPow2_ij))
        auxs = sqrt((cfPow2_ij-aPow2_ij)/(cfPow2_ij-csPow2_ij))

        ! Compute the "beta_"y,z" values (with left and right states interchanged!)
        auxy = ROE_MEAN_VALUE(\
                YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux)
        auxz = ROE_MEAN_VALUE(\
                ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux)
        aux  = sqrt(POW(auxy,2)+POW(auxz,2))
        auxy = auxy/aux
        auxz = auxz/aux
        
        ! Compute the sign if the magnetic field
        S = sign(RCONST(1.0), RCONST(MHD_XMAGFIELD_CONST))

        ! Compute auxiliary square root
        auxsqr = sqrt(rho_ij*aPow2_ij)

        ! Compute matrix of right eigenvectors
        Reig(1,1) = auxf/aPow2_ij
        Reig(1,2) = RCONST(0.0)
        Reig(1,3) = auxs/aPow2_ij
        Reig(1,4) = RCONST(1.0)/aPow2_ij
        Reig(1,5) = Reig(1,3)
        Reig(1,6) = RCONST(0.0)
        Reig(1,7) = Reig(1,1)

        Reig(2,1) = auxf*(u_ij-cf_ij)/aPow2_ij
        Reig(2,2) = RCONST(0.0)
        Reig(2,3) = auxs*(u_ij-cs_ij)/aPow2_ij
        Reig(2,4) =              u_ij/aPow2_ij
        Reig(2,5) = auxs*(u_ij+cs_ij)/aPow2_ij
        Reig(2,6) = RCONST(0.0)
        Reig(2,7) = auxf*(u_ij+cf_ij)/aPow2_ij

        Reig(3,1) = (auxf*v_ij+auxs*cs_ij*auxy*S)/aPow2_ij
        Reig(3,2) = -rho_ij*auxz
        Reig(3,3) = (auxs*v_ij-auxf*cf_ij*auxy*S)/aPow2_ij
        Reig(3,4) =                          v_ij/aPow2_ij
        Reig(3,5) = (auxs*v_ij+auxf*cf_ij*auxy*S)/aPow2_ij
        Reig(3,6) =  rho_ij*auxz
        Reig(3,7) = (auxf*v_ij-auxs*cs_ij*auxy*S)/aPow2_ij

        Reig(4,1) = (auxf*w_ij+auxs*cs_ij*auxz*S)/aPow2_ij
        Reig(4,2) =  rho_ij*auxy
        Reig(4,3) = (auxs*w_ij-auxf*cf_ij*auxz*S)/aPow2_ij
        Reig(4,4) =                          w_ij/aPow2_ij
        Reig(4,5) = (auxs*w_ij+auxf*cf_ij*auxz*S)/aPow2_ij
        Reig(4,6) = -rho_ij*auxy
        Reig(4,7) = (auxf*w_ij-auxs*cs_ij*auxz*S)/aPow2_ij

        Reig(5,1) =  auxs*auxy/auxsqr
        Reig(5,2) = -S*sqrt(rho_ij)*auxz
        Reig(5,3) = -auxf*auxy/auxsqr
        Reig(5,4) =  RCONST(0.0)
        Reig(5,5) =  Reig(5,3)
        Reig(5,6) =  Reig(5,2)
        Reig(5,7) =  Reig(5,1)

        Reig(6,1) =  auxs*auxz/auxsqr
        Reig(6,2) =  S*sqrt(rho_ij)*auxy
        Reig(6,3) = -auxf*auxz/auxsqr
        Reig(6,4) =  RCONST(0.0)
        Reig(6,5) =  Reig(6,3)
        Reig(6,6) =  Reig(6,2)
        Reig(6,7) =  Reig(6,1)

        Reig(7,1) = (auxf*(H_ij-bPow2_ij-u_ij*cf_ij)+&
                           auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-auxs*aux/auxsqr
        Reig(7,2) = -rho_ij*(v_ij*auxz-w_ij*auxy)
        Reig(7,3) = (auxs*(H_ij-bPow2_ij-u_ij*cs_ij)-&
                           auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-auxf*aux/auxsqr
        Reig(7,4) =  (q_ij+(RCONST(MAGNETOHYDRODYN_GAMMA)-RCONST(2.0))/&
                           (RCONST(MAGNETOHYDRODYN_GAMMA)-RCONST(1.0))*X_ij)/aPow2_ij
        Reig(7,5) = (auxs*(H_ij-bPow2_ij+u_ij*cs_ij)+&
                           auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-auxf*aux/auxsqr
        Reig(7,6) = -Reig(7,2)
        Reig(7,7) = (auxf*(H_ij-bPow2_ij+u_ij*cf_ij)-&
                           auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-auxs*aux/auxsqr

        ! Compute characteristic variables by "solving" R_ij * dW = dU
        call dgesv(7, 1, Reig, 7, Ipiv, Diff, 7, info)

        ! Multiply characteristic variables by the corresponding eigenvalue
        w1 = l1 * Diff(1)
        w2 = l2 * Diff(2)
        w3 = l3 * Diff(3)
        w4 = l4 * Diff(4)
        w5 = l5 * Diff(5)
        w6 = l6 * Diff(6)
        w7 = l7 * Diff(7)

        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        Diff(1) = anorm * ( auxf*(w1+w7) + auxs*(w3+w5) + w4 )/aPow2_ij
        Diff(2) = anorm * ( auxf*(u_ij-cf_ij)*w1 + auxs*(u_ij-cs_ij)*w3 + u_ij*w4 +&
                            auxs*(u_ij+cs_ij)*w5 + auxf*(u_ij+cf_ij)*w7 )/aPow2_ij
        Diff(3) = anorm * ( ((auxf*v_ij+auxs*cs_ij*auxy*S)*w1 +&
                             (auxs*v_ij-auxf*cf_ij*auxy*S)*w3 +&
                                                      v_ij*w4 +&
                             (auxs*v_ij+auxf*cf_ij*auxy*S)*w5 +&
                             (auxf*v_ij-auxs*cs_ij*auxy*S)*w7)/aPow2_ij +&
                            rho_ij*auxz*(-w2+w6) )
        Diff(4) = anorm * ( ((auxf*w_ij+auxs*cs_ij*auxz*S)*w1 +&
                             (auxs*w_ij-auxf*cf_ij*auxz*S)*w3 +&
                                                      w_ij*w4 +&
                             (auxs*w_ij+auxf*cf_ij*auxz*S)*w5 +&
                             (auxf*w_ij-auxs*cs_ij*auxz*S)*w7)/aPow2_ij +&
                            rho_ij*auxy*(w2-w6) )
        Diff(5) = anorm * ( (auxs*auxy*(w1+w7) -&
                             auxf*auxy*(w3+w5))/sqrt(rho_ij*aPow2_ij) -&
                            S*sqrt(rho_ij)*auxz*(w2+w6) )
        Diff(6) = anorm * ( (auxs*auxz*(w1+w7) -&
                             auxf*auxz*(w3+w5))/sqrt(rho_ij*aPow2_ij) +&
                            S*sqrt(rho_ij)*auxy*(w2+w6) )
        Diff(7) = anorm * ( ((auxf*(H_ij-bPow2_ij-u_ij*cf_ij)+&
                              auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))*w1 +&
                             (auxf*(H_ij-bPow2_ij+u_ij*cf_ij)-&
                              auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))*w7)/aPow2_ij -&
                             auxs*aux*(w1+w7)/sqrt(rho_ij*aPow2_ij) +&
                             rho_ij*(v_ij*auxz-w_ij*auxy)*(-w2+w6) +&
                            ((auxs*(H_ij-bPow2_ij-u_ij*cs_ij)-&
                              auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))*w3 +&
                             (auxs*(H_ij-bPow2_ij+u_ij*cs_ij)+&
                              auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))*w5)/aPow2_ij -&
                             auxf*aux*(w3+w5)/sqrt(rho_ij*aPow2_ij) +&
                             (q_ij+(RCONST(MAGNETOHYDRODYN_GAMMA)-RCONST(2.0))/&
                                   (RCONST(MAGNETOHYDRODYN_GAMMA)-RCONST(1.0))*X_ij)*w4/aPow2_ij )

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fj-&
             IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fi + Diff)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
            (IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*F_ij + Diff)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*F_ij + Diff)
#endif
      else

#ifdef MHD_USE_IBP
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fj-&
             IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fi )
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
        IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
            (IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*F_ij)
        IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
            (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*F_ij)
#endif
      end if
    end do

  end subroutine mhd_calcFluxRoeDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxRusDiss1d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)


!<description>
    ! This subroutine computes the fluxes for the low-order scheme in
    ! 1D using scalar artificial viscosities of Rusanov-type.
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
    real(DP), dimension(NVAR1D) :: Fi,Fj
#else
    real(DP), dimension(NVAR1D) :: F_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: cai,caj,cfi,cfj,d_ij
    real(DP) :: aPow2i,aPow2j,astPow2i,astPow2j
    integer :: idx

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
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
      Fi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)
      Fi(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)

      Fj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      Fj(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      F_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(4) = INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX4_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(5) = INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX5_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(6) = INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX6_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
      F_ij(7) = INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,vi,wi,pi,qi)-&
                INVISCIDFLUX7_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,vj,wj,pj,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !
      ! There are seven eigenvalues
      !
      !   u-cf, u-ca, u-cs, u, u+cs, u+ca, u+cf,
      !
      ! where u is the x-velocity component and ca, cs and cf are the
      ! velocities of the Alfveen waves, the slow and fast waves. Since
      !
      !   cf >= ca >= cs >= 0
      !
      ! it suffices to consider only the two eigenvalues
      !
      !   u-cf and u+cf
      !
      ! to construct the Rusanov fluxes
      ! -------------------------------------------------------------------------
      
      ! Compute the speed of the Alfven waves
      cai = abs(XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0))
      caj = abs(XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0))

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

      ! Compute the speed of the fast waves
      cfi = sqrt(RCONST(0.5)*(astPow2i+&
                 sqrt(POW(astPow2i,2)-RCONST(4.0)*aPow2i*POW(cai,2))))
      cfj = sqrt(RCONST(0.5)*(astPow2j+&
                 sqrt(POW(astPow2j,2)-RCONST(4.0)*aPow2j*POW(caj,2))))
            
#ifdef MHD_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*uj)+&
                 RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-\
                                      IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),2))*cfj,&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,0,0,0))*ui)+&
                 RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)-\
                                      IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),2))*cfi )
#else
      ! Compute scalar dissipation
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*uj)+&
                  abs(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0))*cfj,&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*ui)+&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*cfi )
#endif
      
      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
                   IDX3(DdataAtEdge,:,1,idx,0,0,0))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fj-&
           IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fi + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale *&
          (IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*F_ij + Diff)
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*F_ij + Diff)
#endif
    end do

  end subroutine mhd_calcFluxRusDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiagMatD1d_sim(DdataAtNode, DcoeffsAtNode,&
      InodeList, dscale, nnodes, DmatrixAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 1D.
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
    integer :: inode


    do inode = 1, nnodes
      
      ! Set coefficient to zero
      DmatrixAtNode(:,:,inode) = 0.0
    end do

  end subroutine mhd_calcMatDiagMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiag1d_sim(DdataAtNode,&
      DcoeffsAtNode, InodeList, dscale, nnodes,&
      DmatrixAtNode, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 1D.
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
    integer :: inode


    do inode = 1, nnodes
      
      ! Set coefficient to zero
      DmatrixAtNode(:,:,inode) = 0.0
    end do
  
  end subroutine mhd_calcMatDiag1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGalMatD1d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 1D.
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
    integer :: idx


    do idx = 1, nedges
      
      ! Set coefficient to zero
      DmatrixAtEdge(:,:,idx) = 0.0
    end do

  end subroutine mhd_calcMatGalMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGal1d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 1D.
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
    integer :: idx


    do idx = 1, nedges
      
      ! Set coefficient to zero
      DmatrixAtEdge(:,:,idx) = 0.0
    end do
   
  end subroutine mhd_calcMatGal1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDissMatD1d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 1D and applies scalar artificial viscosities proportional to
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
    integer :: idx


    do idx = 1, nedges
      
      ! Set coefficient to zero
      DmatrixAtEdge(:,:,idx) = 0.0
    end do

  end subroutine mhd_calcMatScDissMatD1d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDiss1d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 1D and applies
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
    integer :: idx


    do idx = 1, nedges
      
      ! Set coefficient to zero
      DmatrixAtEdge(:,:,idx) = 0.0
    end do
   
  end subroutine mhd_calcMatScDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDissMatD1d_sim(DdataAtEdge,&
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
    integer :: idx


    do idx = 1, nedges
      
      ! Set coefficient to zero
      DmatrixAtEdge(:,:,idx) = 0.0
    end do

  end subroutine mhd_calcMatRoeDissMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDiss1d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 1D and applies
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
    integer :: idx


    do idx = 1, nedges
      
      ! Set coefficient to zero
      DmatrixAtEdge(:,:,idx) = 0.0
    end do
  
  end subroutine mhd_calcMatRoeDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDissMatD1d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 1D and applies the scalar artificial viscosities of Rusanov-type.
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
    integer :: idx


    do idx = 1, nedges
      
      ! Set coefficient to zero
      DmatrixAtEdge(:,:,idx) = 0.0
    end do
    
  end subroutine mhd_calcMatRusDissMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDiss1d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 1D applies
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
    integer :: idx


    do idx = 1, nedges
      
      ! Set coefficient to zero
      DmatrixAtEdge(:,:,idx) = 0.0
    end do
  
  end subroutine mhd_calcMatRusDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcCharacteristics1d_sim(Dweight, DdataAtEdge,&
      nedges, DcharVariablesAtEdge, DeigenvaluesAtEdge,&
      DrightEigenvectorsAtEdge, DleftEigenvectorsAtEdge, rcollection)

!<description>
    ! This subroutine computes the characteristic variables in 1D.
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

    
  end subroutine mhd_calcCharacteristics1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTScDiss1d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the raw antidiffusive fluxes for FCT
    ! algorithms in 1D using scalar dissipation proportional to the
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
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,X_ij,aux,cf_ij,d_ij,q_ij,rho_ij,u_ij
    real(DP) :: aPow2_ij,astPow2_ij,bPow2_ij,bxPow2_ij
    integer :: idx


    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      wi = ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      wj = ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute skew-symmetric coefficient and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                       IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))

      ! Compute densities
      ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute total pressures
      pi = TOTALPRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = TOTALPRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute enthalpies
      hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
      hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(ri,rj)
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      H_ij = ROE_MEAN_VALUE(hi,hj,aux)
    
      ! Compute the Roe-averaged density with left and right states interchanged!
      rho_ij = ROE_MEAN_VALUE(rj,ri,aux)

      ! Compute the square of the Roe-averaged speed of the Alfven waves.
      ! Note that left and right states are interchanged!
      bxPow2_ij = (POW(ROE_MEAN_VALUE(\
                       XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                       XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2))/rho_ij
    
      ! Compute the density-averaged magnetic field
      X_ij = (POW(XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)-XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),2)+&
              POW(YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)-YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),2)+&
              POW(ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)-ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),2))/&
              RCONST(2.0)*POW(sqrt(ri)+sqrt(rj),2)

      ! Compute the square of the Roe-averaged magnetic field.
      ! Note that left and right states are interchanged!
      bPow2_ij = (POW(ROE_MEAN_VALUE(\
                      XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                      XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2)+\
                  POW(ROE_MEAN_VALUE(\
                      YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                      YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2)+\
                  POW(ROE_MEAN_VALUE(\
                      ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                      ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2))/rho_ij

      ! Compute the magnitude of the Roe-averaged velocity
      q_ij = (POW(ROE_MEAN_VALUE(ui,uj,aux),2)+&
              POW(ROE_MEAN_VALUE(vi,vj,aux),2)+&
              POW(ROE_MEAN_VALUE(wi,wj,aux),2))*RCONST(0.5)


      ! Compute the Roe-averaged speed of sound
      aPow2_ij = (RCONST(2.0)-(MAGNETOHYDRODYN_GAMMA))*X_ij +&
          (MAGNETOHYDRODYN_GAMMA-RCONST(1.0))*(H_ij-q_ij-bPow2_ij)
      
      ! Compute auxiliary variables
      astPow2_ij = aPow2_ij+bPow2_ij
      aux        = sqrt(POW(astPow2_ij,2)-RCONST(4.0)*aPow2_ij*bxPow2_ij)
            
      ! Compute the Roe-averagred speed of the fast waves
      cf_ij = sqrt(RCONST(0.5)*(astPow2_ij+aux))

      ! Scalar dissipation
      d_ij = abs(u_ij*a(1)) + abs(a(1))*cf_ij
      
      ! Compute conservative fluxes
      IDX2(DfluxesAtEdge,:,idx,0,0) = dscale*&
          d_ij*(IDX3(DdataAtEdge,:,1,idx,0,0,0)-&
                IDX3(DdataAtEdge,:,2,idx,0,0,0))
    end do

  end subroutine mhd_calcFluxFCTScDiss1d_sim

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcFluxFCTRoeDiss1d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes forFCT
    ! algorithms in 1D using tensorial dissipation of Roe-type.
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
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP), dimension(7,7) :: Reig
    real(DP) :: hi,hj,pi,pj,qi,qj,ri,rj,ui,uj,vi,vj,wi,wj
    real(DP) :: caPow2_ij,ca_ij,cfPow2_ij,cf_ij,csPow2_ij,cs_ij
    real(DP) :: S,aux,auxf,auxs,auxsqr,auxy,auxz
    real(DP) :: H_ij,X_ij,q_ij,rho_ij,u_ij,v_ij,w_ij
    real(DP) :: aPow2_ij,astPow2_ij,bPow2_ij,bxPow2_ij
    real(DP) :: l1,l2,l3,l4,l5,l6,l7,w1,w2,w3,w4,w5,w6,w7
    real(DP) :: anorm
    integer, dimension(7) :: Ipiv
    integer :: idx,info

    
    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      wi = ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      wj = ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !
      ! There are seven eigenvalues
      !
      !   u-cf, u-ca, u-cs, u, u+cs, u+ca, u+cf,
      !
      ! where u is the x-velocity component and ca, cs and cf are the
      ! velocities of the Alfveen waves, the slow and fast waves.
      !-------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                       IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))
      anorm = abs(a(1)) ! = sqrt(a(1)*a(1))

      if (anorm .gt. SYS_EPSREAL_DP) then

        ! Compute densities
        ri = DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        rj = DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        
        ! Compute total pressures
        pi = TOTALPRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
        pj = TOTALPRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)
        
        ! Compute enthalpies
        hi = (TOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)+pi)/ri
        hj = (TOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)+pj)/rj
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(ri,rj)
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)

        ! Compute the Roe-averaged density with left and right states interchanged!
        rho_ij = ROE_MEAN_VALUE(rj,ri,aux)
        
        ! Compute the square of the Roe-averaged speed of the Alfven waves.
        ! Note that left and right states are interchanged!
        bxPow2_ij = (POW(ROE_MEAN_VALUE(\
                       XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                       XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2))/rho_ij
        ca_ij = sqrt(bxPow2_ij)
    
        ! Compute the density-averaged magnetic field
        X_ij = (POW(XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)-XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),2)+&
                POW(YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)-YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),2)+&
                POW(ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)-ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),2))/&
                RCONST(2.0)*POW(sqrt(ri)+sqrt(rj),2)

        ! Compute the square of the Roe-averaged magnetic field.
        bPow2_ij = (POW(ROE_MEAN_VALUE(\
                      XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                      XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2)+\
                  POW(ROE_MEAN_VALUE(\
                      YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                      YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2)+\
                  POW(ROE_MEAN_VALUE(\
                      ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                      ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux),2))/rho_ij

        ! Compute the magnitude of the Roe-averaged velocity
        q_ij = (POW(ROE_MEAN_VALUE(ui,uj,aux),2)+&
              POW(ROE_MEAN_VALUE(vi,vj,aux),2)+&
              POW(ROE_MEAN_VALUE(wi,wj,aux),2))*RCONST(0.5)

        ! Compute the Roe-averaged speed of sound
        aPow2_ij = (RCONST(2.0)-(MAGNETOHYDRODYN_GAMMA))*X_ij +&
            (MAGNETOHYDRODYN_GAMMA-RCONST(1.0))*(H_ij-q_ij-bPow2_ij)

        ! Compute auxiliary quantities
        astPow2_ij = aPow2_ij+bPow2_ij
        auxsqr     = sqrt(POW(astPow2_ij,2)-RCONST(4.0)*aPow2_ij*bxPow2_ij)

        ! Compute the Roe-averagred speed of the slow and fast waves
        cfPow2_ij = RCONST(0.5)*(astPow2_ij+auxsqr); cf_ij=sqrt(cfPow2_ij)
        csPow2_ij = RCONST(0.5)*(astPow2_ij-auxsqr); cs_ij=sqrt(csPow2_ij)

        ! Compute eigenvalues
        l1 = abs(u_ij-cf_ij)
        l2 = abs(u_ij-ca_ij)
        l3 = abs(u_ij-cs_ij)
        l4 = abs(u_ij)
        l5 = abs(u_ij+cs_ij)
        l6 = abs(u_ij+ca_ij)
        l7 = abs(u_ij+cf_ij)
        
        ! Compute solution difference U_i-U_j
        Diff = IDX3(DdataAtEdge,:,1,idx,0,0,0)-&
               IDX3(DdataAtEdge,:,2,idx,0,0,0)

        ! The variable names are adopted from J.M. Stone, T.A. Gardiner,
        ! P. Teuben, J.F. Hawlay and J.B. Simon ATHENA: A New Code for
        ! Astrophysical MHD. The Astrophysocal Journal Supplement Series
        ! (ISSN 0067-0049), vol. 178, September 2008, p. 137-177.

        ! Compute the "alpha_f,s" values
        auxf = sqrt((aPow2_ij-csPow2_ij)/(cfPow2_ij-csPow2_ij))
        auxs = sqrt((cfPow2_ij-aPow2_ij)/(cfPow2_ij-csPow2_ij))

        ! Compute the "beta_"y,z" values (with left and right states interchanged!)
        auxy = ROE_MEAN_VALUE(\
                YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux)
        auxz = ROE_MEAN_VALUE(\
                ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0),\
                ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0),aux)
        aux  = sqrt(POW(auxy,2)+POW(auxz,2))
        auxy = auxy/aux
        auxz = auxz/aux
        
        ! Compute the sign if the magnetic field
        S = sign(RCONST(1.0), RCONST(MHD_XMAGFIELD_CONST))

        ! Compute auxiliary square root
        auxsqr = sqrt(rho_ij*aPow2_ij)

        ! Compute matrix of right eigenvectors
        Reig(1,1) = auxf/aPow2_ij
        Reig(1,2) = RCONST(0.0)
        Reig(1,3) = auxs/aPow2_ij
        Reig(1,4) = RCONST(1.0)/aPow2_ij
        Reig(1,5) = Reig(1,3)
        Reig(1,6) = RCONST(0.0)
        Reig(1,7) = Reig(1,1)

        Reig(2,1) = auxf*(u_ij-cf_ij)/aPow2_ij
        Reig(2,2) = RCONST(0.0)
        Reig(2,3) = auxs*(u_ij-cs_ij)/aPow2_ij
        Reig(2,4) =              u_ij/aPow2_ij
        Reig(2,5) = auxs*(u_ij+cs_ij)/aPow2_ij
        Reig(2,6) = RCONST(0.0)
        Reig(2,7) = auxf*(u_ij+cf_ij)/aPow2_ij

        Reig(3,1) = (auxf*v_ij+auxs*cs_ij*auxy*S)/aPow2_ij
        Reig(3,2) = -rho_ij*auxz
        Reig(3,3) = (auxs*v_ij-auxf*cf_ij*auxy*S)/aPow2_ij
        Reig(3,4) =                          v_ij/aPow2_ij
        Reig(3,5) = (auxs*v_ij+auxf*cf_ij*auxy*S)/aPow2_ij
        Reig(3,6) =  rho_ij*auxz
        Reig(3,7) = (auxf*v_ij-auxs*cs_ij*auxy*S)/aPow2_ij

        Reig(4,1) = (auxf*w_ij+auxs*cs_ij*auxz*S)/aPow2_ij
        Reig(4,2) =  rho_ij*auxy
        Reig(4,3) = (auxs*w_ij-auxf*cf_ij*auxz*S)/aPow2_ij
        Reig(4,4) =                          w_ij/aPow2_ij
        Reig(4,5) = (auxs*w_ij+auxf*cf_ij*auxz*S)/aPow2_ij
        Reig(4,6) = -rho_ij*auxy
        Reig(4,7) = (auxf*w_ij-auxs*cs_ij*auxz*S)/aPow2_ij

        Reig(5,1) =  auxs*auxy/auxsqr
        Reig(5,2) = -S*sqrt(rho_ij)*auxz
        Reig(5,3) = -auxf*auxy/auxsqr
        Reig(5,4) =  RCONST(0.0)
        Reig(5,5) =  Reig(5,3)
        Reig(5,6) =  Reig(5,2)
        Reig(5,7) =  Reig(5,1)

        Reig(6,1) =  auxs*auxz/auxsqr
        Reig(6,2) =  S*sqrt(rho_ij)*auxy
        Reig(6,3) = -auxf*auxz/auxsqr
        Reig(6,4) =  RCONST(0.0)
        Reig(6,5) =  Reig(6,3)
        Reig(6,6) =  Reig(6,2)
        Reig(6,7) =  Reig(6,1)

        Reig(7,1) = (auxf*(H_ij-bPow2_ij-u_ij*cf_ij)+&
                           auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-auxs*aux/auxsqr
        Reig(7,2) = -rho_ij*(v_ij*auxz-w_ij*auxy)
        Reig(7,3) = (auxs*(H_ij-bPow2_ij-u_ij*cs_ij)-&
                           auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-auxf*aux/auxsqr
        Reig(7,4) =  (q_ij+(RCONST(MAGNETOHYDRODYN_GAMMA)-RCONST(2.0))/&
                           (RCONST(MAGNETOHYDRODYN_GAMMA)-RCONST(1.0))*X_ij)/aPow2_ij
        Reig(7,5) = (auxs*(H_ij-bPow2_ij+u_ij*cs_ij)+&
                           auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-auxf*aux/auxsqr
        Reig(7,6) = -Reig(7,2)
        Reig(7,7) = (auxf*(H_ij-bPow2_ij+u_ij*cf_ij)-&
                           auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-auxs*aux/auxsqr

        ! Compute characteristic variables by "solving" R_ij * dW = dU
        call dgesv(7, 1, Reig, 7, Ipiv, Diff, 7, info)

        ! Multiply characteristic variables by the corresponding eigenvalue
        w1 = l1 * Diff(1)
        w2 = l2 * Diff(2)
        w3 = l3 * Diff(3)
        w4 = l4 * Diff(4)
        w5 = l5 * Diff(5)
        w6 = l6 * Diff(6)
        w7 = l7 * Diff(7)

        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        Diff(1) = anorm * ( auxf*(w1+w7) + auxs*(w3+w5) + w4 )/aPow2_ij
        Diff(2) = anorm * ( auxf*(u_ij-cf_ij)*w1 + auxs*(u_ij-cs_ij)*w3 + u_ij*w4 +&
                            auxs*(u_ij+cs_ij)*w5 + auxf*(u_ij+cf_ij)*w7 )/aPow2_ij
        Diff(3) = anorm * ( ((auxf*v_ij+auxs*cs_ij*auxy*S)*w1 +&
                             (auxs*v_ij-auxf*cf_ij*auxy*S)*w3 +&
                                                      v_ij*w4 +&
                             (auxs*v_ij+auxf*cf_ij*auxy*S)*w5 +&
                             (auxf*v_ij-auxs*cs_ij*auxy*S)*w7)/aPow2_ij +&
                            rho_ij*auxz*(-w2+w6) )
        Diff(4) = anorm * ( ((auxf*w_ij+auxs*cs_ij*auxz*S)*w1 +&
                             (auxs*w_ij-auxf*cf_ij*auxz*S)*w3 +&
                                                      w_ij*w4 +&
                             (auxs*w_ij+auxf*cf_ij*auxz*S)*w5 +&
                             (auxf*w_ij-auxs*cs_ij*auxz*S)*w7)/aPow2_ij +&
                            rho_ij*auxy*(w2-w6) )
        Diff(5) = anorm * ( (auxs*auxy*(w1+w7) -&
                             auxf*auxy*(w3+w5))/sqrt(rho_ij*aPow2_ij) -&
                            S*sqrt(rho_ij)*auxz*(w2+w6) )
        Diff(6) = anorm * ( (auxs*auxz*(w1+w7) -&
                             auxf*auxz*(w3+w5))/sqrt(rho_ij*aPow2_ij) +&
                            S*sqrt(rho_ij)*auxy*(w2+w6) )
        Diff(7) = anorm * ( ((auxf*(H_ij-bPow2_ij-u_ij*cf_ij)+&
                              auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))*w1 +&
                             (auxf*(H_ij-bPow2_ij+u_ij*cf_ij)-&
                              auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))*w7)/aPow2_ij -&
                             auxs*aux*(w1+w7)/sqrt(rho_ij*aPow2_ij) +&
                             rho_ij*(v_ij*auxz-w_ij*auxy)*(-w2+w6) +&
                            ((auxs*(H_ij-bPow2_ij-u_ij*cs_ij)-&
                              auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))*w3 +&
                             (auxs*(H_ij-bPow2_ij+u_ij*cs_ij)+&
                              auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))*w5)/aPow2_ij -&
                             auxf*aux*(w3+w5)/sqrt(rho_ij*aPow2_ij) +&
                             (q_ij+(RCONST(MAGNETOHYDRODYN_GAMMA)-RCONST(2.0))/&
                                   (RCONST(MAGNETOHYDRODYN_GAMMA)-RCONST(1.0))*X_ij)*w4/aPow2_ij )

        ! Compute antidiffusive flux
        IDX2(DfluxesAtEdge,:,idx,0,0) = dscale*Diff
      else
        ! Clear antidiffusive flux
        IDX2(DfluxesAtEdge,:,idx,0,0) = 0
      end if
    end do

  end subroutine mhd_calcFluxFCTRoeDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTRusDiss1d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for FCT
    ! algorithms in 1D using scalar dissipation of Rusanov-type.
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
    real(DP) :: ui,uj,vi,vj,wi,wj
    real(DP) :: cai,caj,cfi,cfj,d_ij
    real(DP) :: aPow2i,aPow2j,astPow2i,astPow2j
    integer :: idx
    

    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      vi = YVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      wi = ZVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      vj = YVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      wj = ZVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !
      ! There are seven eigenvalues
      !
      !   u-cf, u-ca, u-cs, u, u+cs, u+ca, u+cf,
      !
      ! where u is the x-velocity component and ca, cs and cf are the
      ! velocities of the Alfveen waves, the slow and fast waves. Since
      !
      !   cf >= ca >= cs >= 0
      !
      ! it suffices to consider only the two eigenvalues
      !
      !   u-cf and u+cf
      !
      ! to construct the Rusanov fluxes
      ! -------------------------------------------------------------------------

      ! Compute the speed of the Alfven waves
      cai = abs(XMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0))
      caj = abs(XMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0))

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

      ! Compute the speed of the fast waves
      cfi = sqrt(RCONST(0.5)*(astPow2i+&
                 sqrt(POW(astPow2i,2)-RCONST(4.0)*aPow2i*POW(cai,2))))
      cfj = sqrt(RCONST(0.5)*(astPow2j+&
                 sqrt(POW(astPow2j,2)-RCONST(4.0)*aPow2j*POW(caj,2))))
            
#ifdef MHD_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*uj)+&
                 RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-\
                                      IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),2))*cfj,&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,0,0,0))*ui)+&
                 RCONST(0.5)*sqrt(POW(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)-\
                                      IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),2))*cfi )
#else
      ! Compute scalar dissipation
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*uj)+&
                  abs(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0))*cfj,&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*ui)+&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*cfi )
#endif

      ! Compute conservative fluxes
      IDX2(DfluxesAtEdge,:,idx,0,0) = dscale*&
          d_ij*(IDX3(DdataAtEdge,:,1,idx,0,0,0)-&
                IDX3(DdataAtEdge,:,2,idx,0,0,0))
    end do

  end subroutine mhd_calcFluxFCTRusDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDensity1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density in 1D.
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

  end subroutine mhd_trafoFluxDensity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDensity1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density in 1D.
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

  end subroutine mhd_trafoDiffDensity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalDensity1d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the density in 1D.
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

  end subroutine mhd_trafoNodalDensity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxEnergy1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the energy in 1D.
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

  end subroutine mhd_trafoFluxEnergy1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffEnergy1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the energy in 1D.
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

  end subroutine mhd_trafoDiffEnergy1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalEnergy1d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the energy in 1D.
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

  end subroutine mhd_trafoNodalEnergy1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxPressure1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the pressure in 1D.
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
    
  end subroutine mhd_trafoFluxPressure1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffPressure1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the pressure in 1D.
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

  end subroutine mhd_trafoDiffPressure1d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalPressure1d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the pressure in 1D.
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

  end subroutine mhd_trafoNodalPressure1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxVelocity1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative fluxes to fluxes for the velocity in 1D.
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
    
  end subroutine mhd_trafoFluxVelocity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffVelocity1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the velocity in 1D.
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

  end subroutine mhd_trafoDiffVelocity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalVelocity1d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the velocity in 1D.
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

  end subroutine mhd_trafoNodalVelocity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxMomentum1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the momentum in 1D.
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
    
  end subroutine mhd_trafoFluxMomentum1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffMomentum1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the momentum in 1D.
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
    
  end subroutine mhd_trafoDiffMomentum1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalMomentum1d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the momentum in 1D.
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

  end subroutine mhd_trafoNodalMomentum1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDenEng1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density and energy in 1D.
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

  end subroutine mhd_trafoFluxDenEng1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDenEng1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density and energy in 1D.
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

  end subroutine mhd_trafoDiffDenEng1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalDenEng1d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the density and
    ! energy in 1D.
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

  end subroutine mhd_trafoNodalDenEng1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDenPre1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density and energy in 1D.
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

  end subroutine mhd_trafoFluxDenPre1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDenPre1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density and energy in 1D.
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

  end subroutine mhd_trafoDiffDenPre1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalDenPre1d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the density and
    ! pressure in 1D.
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

  end subroutine mhd_trafoNodalDenPre1d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDenPreVel1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density, pressure
    ! and velocity in 1D.
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
      
  end subroutine mhd_trafoFluxDenPreVel1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDenPreVel1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density,
    ! pressure and velocity in 1D.
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

  end subroutine mhd_trafoDiffDenPreVel1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalDenPreVel1d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the density,
    ! pressure and velocity in 1D.
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

  end subroutine mhd_trafoNodalDenPreVel1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxMagfield1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the magnetic field in 1D.
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

      ! Transformed magnetic field fluxes in y-direction
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          YMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -YMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)

      ! Transformed magnetic field fluxes in z-direction
      IDX3(DtransformedFluxesAtEdge,2,1,idx,0,0,0) =&
          ZMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,2,2,idx,0,0,0) =&
         -ZMAGFIELD2(DfluxesAtEdge,IDX2,idx,0,0)
    end do

  end subroutine mhd_trafoFluxMagfield1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffMagfield1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative convervative to differences for the magnetic
    ! field in 1D.
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
      
      ! Transformed magnetic field difference in y-direction
      IDX2(DtransformedDataAtEdge,1,idx,0,0) =&
          YMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          YMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)

      ! Transformed magnetic field difference in z-direction
      IDX2(DtransformedDataAtEdge,2,idx,0,0) =&
          ZMAGFIELD3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          ZMAGFIELD3(DdataAtEdge,IDX3,1,idx,0,0,0)
    end do

  end subroutine mhd_trafoDiffMagfield1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoNodalMagfield1d_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to nodal values for the magnetic field in 1D.
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
      
      ! Transformed y-component of the magnetic field
      IDX2(DtransformedDataAtNode,1,idx,0,0) =&
          YMAGFIELD2(DdataAtNode,IDX2,idx,0,0)

      ! Transformed z-component of the magnetic field
      IDX2(DtransformedDataAtNode,2,idx,0,0) =&
          ZMAGFIELD2(DdataAtNode,IDX2,idx,0,0)
    end do

  end subroutine mhd_trafoNodalMagfield1d_sim

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcBoundaryvalues1d(DbdrNormal, DpointNormal,&
      DbdrValue, ibdrCondType, Du, Du0, istatus)

!<description>
    ! This subroutine computes the boundary values for a given node in 1D.
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

  end subroutine mhd_calcBoundaryvalues1d

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcBilfBdrCond1D(rproblemLevel, rboundaryCondition,&
      rsolution, dtime, dscale, ssectionName, fcoeff_buildMatrixScBdr1D_sim,&
      rmatrix, rcollection, cconstrType)

!<description>
    ! This subroutine computes the bilinear form arising from the weak
    ! imposition of boundary conditions in 1D.
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
    include '../../../../../kernel/DOFMaintenance/intf_coefficientMatrixScBdr1D.inc'

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

  end subroutine mhd_calcBilfBdrCond1D

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcLinfBdrCond1D(rproblemLevel, rboundaryCondition,&
      rsolution, dtime, dscale, ssectionName, fcoeff_buildVectorBlBdr1D_sim,&
      rvector, rcollection)

!<description>
    ! This subroutine computes the linear form arising from the weak
    ! imposition of boundary conditions in 1D.
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
    include '../../../../../kernel/DOFMaintenance/intf_coefficientVectorBlBdr1D.inc'
!</intput>

!<inputoutput>
    ! residual/right-hand side vector
    type(t_vectorBlock), intent(inout) :: rvector

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>
!</subroutine>

     ! local variables
    type(t_collection) :: rcollectionTmp
    type(t_linearForm) :: rform
    integer, dimension(:), pointer :: p_IbdrCondType
    integer :: ibct

    ! Evaluate linear form for boundary integral and return if
    ! there are no weak boundary conditions available
    if (.not.rboundaryCondition%bWeakBdrCond) return

    ! Check if we are in 1D
    if (rproblemLevel%rtriangulation%ndim .ne. NDIM1D) then
      call output_line('Spatial dimension must be 1D!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcLinfBdrCond1D')
      call sys_halt()
    end if

    ! Initialise temporal collection structure
    call collct_init(rcollectionTmp)

    ! Prepare quick access arrays of temporal collection structure
    rcollectionTmp%SquickAccess(1) = ''
    rcollectionTmp%SquickAccess(2) = 'rfparser'
    rcollectionTmp%DquickAccess(1) = dtime
    rcollectionTmp%DquickAccess(2) = dscale

    ! Attach user-defined collection structure to temporal collection
    ! structure (may be required by the callback function)
    rcollectionTmp%p_rnextCollection => rcollection
    
    ! Attach solution vector to temporal collection structure
    rcollectionTmp%p_rvectorQuickAccess1 => rsolution
    
    ! Attach function parser from boundary conditions to collection
    ! structure and specify its name in quick access string array
    call collct_setvalue_pars(rcollectionTmp, 'rfparser',&
        rboundaryCondition%rfparser, .true.)
    
    
    ! Set pointers
    call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)
    
    ! Loop over all boundary components
    do ibct = 1, rboundaryCondition%iboundarycount
      
      ! Check if this component has weak boundary conditions
      if (iand(p_IbdrCondType(ibct), BDRC_WEAK) .ne. BDRC_WEAK) cycle

      ! Prepare further quick access arrays of temporal collection
      ! structure with boundary component, type and maximum expressions
      rcollectionTmp%IquickAccess(1) = p_IbdrCondType(ibct)
      rcollectionTmp%IquickAccess(2) = ibct
      rcollectionTmp%IquickAccess(3) = rboundaryCondition%nmaxExpressions
      
      ! Initialise the linear form
      rform%itermCount = 1
      rform%Idescriptors(1) = DER_FUNC
      
      ! Assemble the linear form
      if (rvector%nblocks .eq. 1) then
        call linf_buildVecIntlScalarBdr1d(rform, .false.,&
            rvector%RvectorBlock(1), fcoeff_buildVectorBlBdr1D_sim,&
            ibct, rcollectionTmp)
      else
        call linf_buildVectorBlockBdr1d(rform, .false.,&
            rvector, fcoeff_buildVectorBlBdr1D_sim,&
            ibct, rcollectionTmp)
      end if
      
    end do ! ibct
    
    ! Release temporal collection structure
    call collct_done(rcollectionTmp)

  end subroutine mhd_calcLinfBdrCond1D

  ! ***************************************************************************

!<subroutine>

  subroutine mhd_coeffVectorBdr1d_sim(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
    !
    ! The routine accepts a set of points (cubature points) in real
    ! coordinates.  According to the terms in the linear form, the
    ! routine has to compute simultaneously for all these points and
    ! all the terms in the linear form the corresponding coefficients
    ! in front of the terms.
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

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
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
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(nbocks,itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:,:), intent(out) :: Dcoefficients
!</output>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution
    real(DP), dimension(:,:), pointer :: Daux1
    real(DP), dimension(:,:,:), pointer :: Daux2
    real(DP), dimension(NVAR1D) :: DstateI,DstateM,Dflux,Diff
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dnx,dtime,dscale,pI,cI,rM,pM,cM,dvnI,dvnM,w1,w3
    integer :: ibdrtype,isegment,iel,ipoint,ndim,ivar,nvar,iexpr,nmaxExpr

#ifndef MHD_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DMHD_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'mhd_coeffVectorBdr1d_sim')
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

    ! Check if the solution is given in block or interleaved format
    if (p_rsolution%nblocks .eq. 1) then
      nvar = p_rsolution%RvectorBlock(1)%NVAR
    else
      nvar = p_rsolution%nblocks
    end if
    
    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first three quick access integer values hold:
    ! - the type of boundary condition
    ! - the segment number
    ! - the maximum number of expressions
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)
    nmaxExpr = rcollection%IquickAccess(3)

    if (p_rsolution%nblocks .eq. 1) then

      !-------------------------------------------------------------------------
      ! Solution is stored in interleaved format
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(Daux1(npointsPerElement*nvar, nelements))
      
      ! Evaluate the solution in the cubature points on the boundary
      call fevl_evaluate_sim(DER_FUNC, Daux1, p_rsolution%RvectorBlock(1),&
          Dpoints, rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! What type of boundary conditions are we?
      select case(iand(ibdrtype, BDRC_TYPEMASK))

      case (BDRC_SUPEROUTLET)
        !-----------------------------------------------------------------------
        ! Supersonic outlet boundary conditions:
        !
        ! Evaluate the boundary fluxes based on the computed state vector
        
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Get the normal vector in the point from the boundary
            dnx = merge(1.0_DP, -1.0_DP, mod(ibct,2) .eq. 0)
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux1((ipoint-1)*NVAR1D+1,iel)
            DstateI(2) = Daux1((ipoint-1)*NVAR1D+2,iel)
            DstateI(3) = Daux1((ipoint-1)*NVAR1D+3,iel)
            DstateI(4) = Daux1((ipoint-1)*NVAR1D+4,iel)
            DstateI(5) = Daux1((ipoint-1)*NVAR1D+5,iel)
            DstateI(6) = Daux1((ipoint-1)*NVAR1D+6,iel)
            DstateI(7) = Daux1((ipoint-1)*NVAR1D+7,iel)
            
            ! Assemble Galerkin fluxes at the boundary
            call doGalerkinFlux(DstateI, dnx, Dflux)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*Dflux
          end do
        end do

      case default
        call output_line('Invalid type of boundary conditions!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mhd_coeffVectorBdr1d_sim')
        call sys_halt()
        
      end select
      
      ! Deallocate temporal memory
      deallocate(Daux1)

    else

      !-------------------------------------------------------------------------
      ! Solution is stored in block format
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(Daux2(npointsPerElement*nvar, nelements, nvar))
      
      ! Evaluate the solution in the cubature points on the boundary
      do ivar = 1, nvar
        call fevl_evaluate_sim(DER_FUNC, Daux2(:,:,ivar),&
            p_rsolution%RvectorBlock(ivar), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      end do
      
      ! What type of boundary conditions are we?
      select case(iand(ibdrtype, BDRC_TYPEMASK))
        
      case (BDRC_SUPEROUTLET)
        !-----------------------------------------------------------------------
        ! Supersonic outlet boundary conditions:
        !
        ! Evaluate the boundary fluxes based on the computed state vector

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement

            ! Get the normal vector in the point from the boundary
            dnx = merge(1.0_DP, -1.0_DP, mod(ibct,2) .eq. 0)
        
            ! Setup the computed internal state vector
            DstateI(1) = Daux2(ipoint,iel,1)
            DstateI(2) = Daux2(ipoint,iel,2)
            DstateI(3) = Daux2(ipoint,iel,3)
            DstateI(4) = Daux2(ipoint,iel,4)
            DstateI(5) = Daux2(ipoint,iel,5)
            DstateI(6) = Daux2(ipoint,iel,6)
            DstateI(7) = Daux2(ipoint,iel,7)

            ! Assemble Galerkin fluxes at the boundary
            call doGalerkinFlux(DstateI, dnx, Dflux)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*Dflux
          end do
        end do

      case default
        call output_line('Invalid type of boundary conditions!',&
            OU_CLASS_ERROR,OU_MODE_STD,'c')
        call sys_halt()
        
      end select

      ! Deallocate temporal memory
      deallocate(Daux2)

    end if

  contains

    ! Here come the working routines
    
    !***************************************************************************
    ! Approximate Riemann solver along the outward unit normal
    !***************************************************************************
    
    subroutine doRiemannSolver(DstateI, DstateM, dnx, Dflux, Diff)
      
      ! input parameters
      real(DP), dimension(NVAR1D), intent(in) :: DstateI, DstateM
      real(DP), intent(in) :: dnx
      
      ! output parameters
      real(DP), dimension(NVAR1D), intent(out) :: Dflux, Diff

      print *, "NOT IMPLEMENTED"
      stop

    end subroutine doRiemannSolver

    !***************************************************************************
    ! Compute the Galerkin flux (used for supersonic outflow)
    !***************************************************************************

    subroutine doGalerkinFlux(Dstate, dnx, Dflux)

      ! input parameters
      real(DP), dimension(NVAR1D), intent(in) :: Dstate
      real(DP), intent(in) :: dnx

      ! output parameters
      real(DP), dimension(NVAR1D), intent(out) :: Dflux

      ! local variables
      real(DP) :: u,v,w,p,q
      
      
      ! Compute auxiliary quantities
      u = XVELOCITY1(Dstate,IDX1)
      v = YVELOCITY1(Dstate,IDX1)
      w = ZVELOCITY1(Dstate,IDX1)
      p = TOTALPRESSURE1(Dstate,IDX1)
      q = MAG_DOT_VEL1(Dstate,IDX1,u,v,w)
      
      ! Calculate ${\bf n}\cdot{\bf F}(U)$
      Dflux(1) = dnx*(XMOMENTUM1(Dstate,IDX1))
      Dflux(2) = dnx*(XMOMENTUM1(Dstate,IDX1)*u+p)
      Dflux(3) = dnx*(YMOMENTUM1(Dstate,IDX1)*u-&
                      XMAGFIELD1(Dstate,IDX1)*&
                      YMAGFIELD1(Dstate,IDX1))
      Dflux(4) = dnx*(ZMOMENTUM1(Dstate,IDX1)*u-&
                      XMAGFIELD1(Dstate,IDX1)*&
                      ZMAGFIELD1(Dstate,IDX1))
      Dflux(5) = dnx*(YMAGFIELD1(Dstate,IDX1)*u-&
                      XMAGFIELD1(Dstate,IDX1)*v)
      Dflux(6) = dnx*(ZMAGFIELD1(Dstate,IDX1)*u-&
                      XMAGFIELD1(Dstate,IDX1)*w)
      Dflux(7) = dnx*((TOTALENERGY1(Dstate,IDX1)+p)*u-&
                      XMAGFIELD1(Dstate,IDX1)*q)

    end subroutine doGalerkinFlux

  end subroutine mhd_coeffVectorBdr1d_sim

end module mhd_callback1d
