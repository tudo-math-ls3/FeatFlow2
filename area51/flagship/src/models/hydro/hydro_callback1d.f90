!##############################################################################
!# ****************************************************************************
!# <name> hydro_callback1d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible Euler/Navier-Stokes equations in 1D.
!#
!# The following callback functions are available:
!#
!# 1.) hydro_calcFluxGal1d_sim
!#     -> Computes fluxes for standard Galerkin scheme
!#
!# 2.) hydro_calcFluxGalNoBdr1d_sim
!#     -> Computes fluxes for standard Galerkin scheme without
!#        assembling the symmetric boundary contribution
!#
!# 3.) hydro_calcFluxScDiss1d_sim
!#     -> Computes fluxes for low-order discretisation adoption
!#        scalar artificial viscosities
!#
!# 4.) hydro_calcFluxRoeDiss1d_sim
!#     -> Computes fluxes for low-order discretisation adopting
!#        tensorial artificial viscosities of Roe-type
!#
!# 5.) hydro_calcFluxRusDiss1d_sim
!#     -> Computes fluxes for low-order discretisation adoption
!#        scalar artificial diffusion of Rusanov-type
!#
!# 6.) hydro_calcMatDiagMatD1d_sim
!#     -> Computes local matrix for diagonal entry
!#
!# 7.) hydro_calcMatDiag1d_sim
!#     -> Computes local matrix for diagonal entry
!#
!# 8.) hydro_calcMatGalMatD1d_sim
!#     -> Computes local matrices for standard Galerkin scheme
!#
!# 9.) hydro_calcMatGal1d_sim
!#     -> Computes local matrices for standard Galerkin scheme
!#
!# 10.) hydro_calcMatScDissMatD1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 11.) hydro_calcMatScDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 12.) hydro_calcMatRoeDissMatD1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities of Roe-type
!#
!# 13.) hydro_calcMatRoeDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities of Roe-type
!#
!# 14.) hydro_calcMatRusDissMatD1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities of Rusanov-type
!#
!# 15.) hydro_calcMatRusDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities of Rusanov-type
!#
!# 16.) hydro_calcCharacteristics1d_sim
!#      -> Computes characteristic variables
!#
!# 17.) hydro_calcFluxFCTScDiss1d_sim
!#      -> Computes fluxes for FCT algorithm adopting scalar
!#         artificial viscosities
!#
!# 18.) hydro_calcFluxFCTRoeDiss1d_sim
!#      -> Computes fluxes for FCT algorithm adopting tensorial
!#         artificial viscosities of Roe-type
!#
!# 19.) hydro_calcFluxFCTRusDiss1d_sim
!#      -> Computes fluxes for FCT algorithm adopting scalar
!#         artificial viscosities of Rusanov-type
!#
!# 20.) hydro_trafoFluxDensity1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 21.) hydro_trafoDiffDensity1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 22.) hydro_trafoNodalDensity1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density values
!#
!# 23.) hydro_trafoFluxEnergy1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 24.) hydro_trafoDiffEnergy1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 25.) hydro_trafoNodalEnergy1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal energy values
!#
!# 26.) hydro_trafoFluxPressure1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 27.) hydro_trafoDiffPressure1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the pressure
!#
!# 28.) hydro_trafoNodalPressure1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal pressure values
!#
!# 29.) hydro_trafoFluxVelocity1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 30.) hydro_trafoDiffVelocity1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 31.) hydro_trafoNodalVelocity1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal velocity values
!#
!# 32.) hydro_trafoFluxMomentum1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 33.) hydro_trafoDiffMomentum1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 34.) hydro_trafoNodalMomentum1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal momentum values
!#
!# 35.) hydro_trafoFluxDenEng1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 36.) hydro_trafoDiffDenEng1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 37.) hydro_trafoNodalDenEng1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density and energy values
!#
!# 38.) hydro_trafoFluxDenPre1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 39.) hydro_trafoDiffDenPre1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 40.) hydro_trafoNodalDenPre1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density and pressure values
!#
!# 41.) hydro_trafoFluxDenPreVel1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 42.) hydro_trafoDiffDenPreVel1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure
!#         and the velocity
!#
!# 43.) hydro_trafoNodalDenPreVel1d_sim
!#      -> Computes the transformation from conservative solution
!#         values to nodal density, pressure and velocity values
!#
!# 44.) hydro_calcBoundaryvalues1d
!#      -> Computes the boundary values for a given node
!#
!# 45.) hydro_hadaptCallbackScalar1d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 1D, whereby the vector is stored in interleave format
!#
!# 46.) hydro_hadaptCallbackBlock1d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 1D, whereby the vector is stored in block format
!#
!# 47.) hydro_coeffVectorBdr1d_sim
!#      -> Calculates the coefficients for the linear form in 1D
!#
!# </purpose>
!##############################################################################

module hydro_callback1d

#define HYDRO_NDIM 1
#include "hydro.h"

  use boundarycondaux
  use collection
  use derivatives
  use domainintegration
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
  use problem
  use scalarpde
  use solveraux
  use spatialdiscretisation
  use storage

  implicit none

  private
  public :: hydro_calcFluxGal1d_sim
  public :: hydro_calcFluxGalNoBdr1d_sim
  public :: hydro_calcFluxScDiss1d_sim
  public :: hydro_calcFluxRoeDiss1d_sim
  public :: hydro_calcFluxRusDiss1d_sim
  public :: hydro_calcMatDiagMatD1d_sim
  public :: hydro_calcMatDiag1d_sim
  public :: hydro_calcMatGalMatD1d_sim
  public :: hydro_calcMatGal1d_sim
  public :: hydro_calcMatScDissMatD1d_sim
  public :: hydro_calcMatScDiss1d_sim
  public :: hydro_calcMatRoeDissMatD1d_sim
  public :: hydro_calcMatRoeDiss1d_sim
  public :: hydro_calcMatRusDissMatD1d_sim
  public :: hydro_calcMatRusDiss1d_sim
  public :: hydro_calcCharacteristics1d_sim
  public :: hydro_calcFluxFCTScDiss1d_sim
  public :: hydro_calcFluxFCTRoeDiss1d_sim
  public :: hydro_calcFluxFCTRusDiss1d_sim
  public :: hydro_trafoFluxDensity1d_sim
  public :: hydro_trafoFluxEnergy1d_sim
  public :: hydro_trafoFluxPressure1d_sim
  public :: hydro_trafoFluxVelocity1d_sim
  public :: hydro_trafoFluxMomentum1d_sim
  public :: hydro_trafoFluxDenEng1d_sim
  public :: hydro_trafoFluxDenPre1d_sim
  public :: hydro_trafoFluxDenPreVel1d_sim
  public :: hydro_trafoDiffDensity1d_sim
  public :: hydro_trafoDiffEnergy1d_sim
  public :: hydro_trafoDiffPressure1d_sim
  public :: hydro_trafoDiffVelocity1d_sim
  public :: hydro_trafoDiffMomentum1d_sim
  public :: hydro_trafoDiffDenEng1d_sim
  public :: hydro_trafoDiffDenPre1d_sim
  public :: hydro_trafoDiffDenPreVel1d_sim
  public :: hydro_trafoNodalDensity1d_sim
  public :: hydro_trafoNodalEnergy1d_sim
  public :: hydro_trafoNodalPressure1d_sim
  public :: hydro_trafoNodalVelocity1d_sim
  public :: hydro_trafoNodalMomentum1d_sim
  public :: hydro_trafoNodalDenEng1d_sim
  public :: hydro_trafoNodalDenPre1d_sim
  public :: hydro_trafoNodalDenPreVel1d_sim
  public :: hydro_calcBoundaryvalues1d
  public :: hydro_coeffVectorBdr1d_sim
  public :: hydro_hadaptCallbackScalar1d
  public :: hydro_hadaptCallbackBlock1d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGal1d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
#ifdef HYDRO_USE_IBP
    real(DP), dimension(NVAR1D) :: Fi,Fj
#else
    real(DP), dimension(NVAR1D) :: F_ij
#endif
    real(DP) :: pi,pj,ui,uj
    integer :: idx


    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute pressures
      pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      Fi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)

      Fj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      
      ! Assemble skew-symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          (IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*Fj-&
           IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*Fi )
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
#else
      ! Compute flux difference for x-direction
      F_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      F_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      F_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
                 
      ! Assemble fluxes
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) =  dscale * IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*F_ij
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = -dscale * IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*F_ij
#endif
    end do

  end subroutine hydro_calcFluxGal1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGalNoBdr1d_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the fluxes for the TVD discretisation
    ! in 1D. The symmetric boundary contributions are neglected and
    ! incorporated into the antidiffusive flux. Hence, this is simply
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
    real(DP) :: pi,pj,ui,uj
    integer :: idx


    do idx = 1, nedges

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !-------------------------------------------------------------------------
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute pressures
      pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute flux difference for x-direction
      F_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      F_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      F_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      
      ! Assemble symmetric fluxes
      IDX3(DfluxesAtEdge,:,1,idx,0,0,0) = dscale *&
          RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                       IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*F_ij
      IDX3(DfluxesAtEdge,:,2,idx,0,0,0) = IDX3(DfluxesAtEdge,:,1,idx,0,0,0)
    end do

  end subroutine hydro_calcFluxGalNoBdr1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxScDiss1d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
#ifdef HYDRO_USE_IBP
    real(DP), dimension(NVAR1D) :: Fi,Fj
#else
    real(DP), dimension(NVAR1D) :: F_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj
    real(DP) :: H_ij,aux,c_ij,d_ij,q_ij,u_ij
    integer :: idx


    do idx = 1, nedges
          
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute pressures
      pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      Fi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)

      Fj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
#else
      ! Compute flux difference for x-direction
      F_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      F_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      F_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
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
      
      ! Compute auxiliary variable
      q_ij = RCONST(0.5)*(u_ij*u_ij)

      ! Compute the speed of sound
      c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))
      
      ! Compute scalar dissipation
      d_ij = abs(a(1)*u_ij) + abs(a(1))*c_ij
      
      ! Multiply the solution difference by the scalar dissipation
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
                   IDX3(DdataAtEdge,:,1,idx,0,0,0))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------
      
#ifdef HYDRO_USE_IBP
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

  end subroutine hydro_calcFluxScDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRoeDiss1d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
#ifdef HYDRO_USE_IBP
    real(DP), dimension(NVAR1D) :: Fi,Fj
#else
    real(DP), dimension(NVAR1D) :: F_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj
    real(DP) :: H_ij,anorm,aux,b1,b2,cPow2_ij,c_ij,q_ij,u_ij
    real(DP) :: l1,l2,l3,w1,w2,w3
    integer :: idx
#if defined(HYDRO_USE_ENTROPYFIX) && (HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX)
    real(DP) :: dtol,ci,cj
#endif

    
    do idx = 1, nedges
          
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute pressures
      pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      Fi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)

      Fj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
#else
      ! Compute flux difference for x-direction
      F_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      F_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      F_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
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
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
                
        ! Compute auxiliary variable
        q_ij = RCONST(0.5)*(u_ij*u_ij)

        ! Compute the speed of sound
        cPow2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij     = sqrt(cPow2_ij)
        
        ! Compute eigenvalues
        l1 = abs(u_ij-c_ij)
        l2 = abs(u_ij)
        l3 = abs(u_ij+c_ij)

#if defined(HYDRO_USE_ENTROPYFIX)

#if HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX

        ! Entropy-fix by Harten and Hyman
        ci = SOUNDSPEED3_1D(DdataAtEdge,IDX3,1,idx,0,0,0)
        cj = SOUNDSPEED3_1D(DdataAtEdge,IDX3,2,idx,0,0,0)

        dtol = max(RCONST(0.0), (u_ij-c_ij) - (ui-ci), (uj-cj) - (u_ij-c_ij) )
        
        if (l1 .lt. dtol)&
            l1 = RCONST(0.5)*((l1**2)/dtol + dtol)
        
        dtol = max(RCONST(0.0), (u_ij+c_ij) - (ui+ci), (uj+cj) - (u_ij+c_ij) )

        if (l3 .lt. dtol)&
            l3 = RCONST(0.5)*((l3**2)/dtol + dtol)

#elif HYDRO_USE_ENTROPYFIX == HARTEN_ENTROPYFIX

#ifndef HYDRO_HARTEN_ENTROPYFIX
#error "Value HYDRO_HARTEN_ENTROPYFIX is required!"
#else
        ! Entropy-fix by Harten
        if (l1 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l1 = RCONST(0.5)*((l1**2)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))

        if (l3 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l3 = RCONST(0.5)*((l3**2)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))
#endif
#else
#error "Invalid type of entropy fix!"
#endif
#endif
        
        ! Compute solution difference U_j-U_i
        Diff = IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
               IDX3(DdataAtEdge,:,1,idx,0,0,0)
        
        ! Compute auxiliary quantities for characteristic variables
        b2 = ((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij; b1 = b2*q_ij
        
        ! Compute characteristic variables multiplied by the
        ! corresponding eigenvalue
        w1 = l1 * RCONST(0.5) * (       (b1+u_ij/c_ij)*Diff(1) -&
                            (b2*u_ij+RCONST(1.0)/c_ij)*Diff(2) +&
                                                    b2*Diff(3) )
        w2 = l2 * (                   (RCONST(1.0)-b1)*Diff(1) +&
                                               b2*u_ij*Diff(2) -&
                                                    b2*Diff(3) )
        w3 = l3 * RCONST(0.5) * (       (b1-u_ij/c_ij)*Diff(1) -&
                            (b2*u_ij-RCONST(1.0)/c_ij)*Diff(2) +&
                                                    b2*Diff(3) )
        
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        Diff(1) = anorm * ( w1 + w2 + w3 )
        Diff(2) = anorm * ( (u_ij-c_ij)*w1 + u_ij*w2 + (u_ij+c_ij)*w3 )
        Diff(3) = anorm * ( (H_ij-u_ij*c_ij)*w1 + q_ij*w2 + (H_ij+u_ij*c_ij)*w3 )
        
        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
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

#ifdef HYDRO_USE_IBP
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

  end subroutine hydro_calcFluxRoeDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRusDiss1d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
#ifdef HYDRO_USE_IBP
    real(DP), dimension(NVAR1D) :: Fi,Fj
#else
    real(DP), dimension(NVAR1D) :: F_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP) :: Ei,Ej,ci,cj,pi,pj,ui,uj
    real(DP) :: d_ij
    integer :: idx

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute pressures
      pi = PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
      pj = PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      Fi(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fi(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)
      Fi(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)

      Fj(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fj(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      Fj(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
#else
      ! Compute flux difference for x-direction
      F_ij(1) = INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                INVISCIDFLUX1_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      F_ij(2) = INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                INVISCIDFLUX2_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
      F_ij(3) = INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,1,idx,0,0,0,ui,pi)-&
                INVISCIDFLUX3_XDIR3(DdataAtEdge,IDX3,2,idx,0,0,0,uj,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !-------------------------------------------------------------------------
      
      ! Compute specific total energies
      Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*&
          (HYDRO_GAMMA)*(Ei-RCONST(0.5)*ui*ui), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*&
          (HYDRO_GAMMA)*(Ej-RCONST(0.5)*uj*uj), SYS_EPSREAL_DP))

#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*uj)+&
                 RCONST(0.5)*sqrt((IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))**2)*cj,&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,0,0,0))*ui)+&
                 RCONST(0.5)*sqrt((IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,0,0,0))**2)*ci )
#else
      ! Compute scalar dissipation
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*uj)+&
                  abs(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0))*cj,&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*ui)+&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*ci )
#endif

      ! Multiply the solution difference by the scalar dissipation
      Diff = d_ij*(IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
                   IDX3(DdataAtEdge,:,1,idx,0,0,0))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
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

  end subroutine hydro_calcFluxRusDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiagMatD1d_sim(DdataAtNode, DcoeffsAtNode,&
      InodeList, dscale, nnodes, DmatrixAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 1D
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
    real(DP) :: ui
    integer :: inode


    do inode = 1, nnodes
      
      ! Compute velocity
      ui = XVELOCITY2(DdataAtNode,IDX2,inode,0,0)
      
#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ii = diag(A_i)*C_{ii}$
      IDX3(DmatrixAtNode,1,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,_)
      IDX3(DmatrixAtNode,2,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,_)
      IDX3(DmatrixAtNode,3,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,_)
#else
      ! Compute Galerkin coefficient $K_ii = -diag(A_i)*C_{ii}$
      IDX3(DmatrixAtNode,1,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,_)
      IDX3(DmatrixAtNode,2,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,_)
      IDX3(DmatrixAtNode,3,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,_)
#endif
    end do

  end subroutine hydro_calcMatDiagMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiag1d_sim(DdataAtNode,&
      DcoeffsAtNode, InodeList, dscale, nnodes,&
      DmatrixAtNode, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 1D
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
    real(DP) :: Ei,ui
    integer :: inode


    do inode = 1, nnodes
      
      ! Compute auxiliary variables
      ui = XVELOCITY2(DdataAtNode,IDX2,inode,0,0)
      Ei = SPECIFICTOTALENERGY2(DdataAtNode,IDX2,inode,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      IDX3(DmatrixAtNode,1,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,2,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,3,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,4,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,5,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,6,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,7,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,8,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,9,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
#else
      ! Compute Galerkin coefficient $K_ii = -A_i*C_{ii}$
      IDX3(DmatrixAtNode,1,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,2,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,3,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,4,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,5,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,6,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,7,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,8,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
      IDX3(DmatrixAtNode,9,1,inode,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX2(DcoeffsAtNode,1,inode,0,0),ui,Ei)
#endif
    end do
    
  end subroutine hydro_calcMatDiag1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGalMatD1d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 1D
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
    real(DP) :: ui,uj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DmatrixAtEdge,1,1,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,2,1,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,3,1,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DmatrixAtEdge,1,1,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,2,1,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,3,1,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,_)
#endif
    end do

  end subroutine hydro_calcMatGalMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGal1d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 1D
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
    real(DP) :: Ei,Ej,ui,uj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute auxiliary variables
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      
      IDX3(DmatrixAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)

      IDX3(DmatrixAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,5,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,6,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,7,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,9,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DmatrixAtEdge,1,1,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,2,1,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,3,1,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)

      IDX3(DmatrixAtEdge,4,1,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,5,1,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,6,1,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)

      IDX3(DmatrixAtEdge,7,1,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,8,1,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,9,1,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
#endif
    end do

  end subroutine hydro_calcMatGal1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDissMatD1d_sim(DdataAtEdge,&
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
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj
    real(DP) :: H_ij,anorm,aux,c_ij,q_ij,u_ij,vel_ij
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
     
#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,_)
#endif
      
      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ji}-C_{ij})$ and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                       IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))
      anorm = abs(a(1))
      
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
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary values
        vel_ij = u_ij*a(1)
        q_ij   = RCONST(0.5)*(u_ij*u_ij)
        
        ! Compute the speed of sound
        c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))

        ! Compute scalar dissipation
        IDX3(DmatrixAtEdge,:,1,idx,0,0,0) = dscale * (abs(vel_ij) + anorm*c_ij)

      else
        
        ! Nullify dissipation tensor
        IDX3(DmatrixAtEdge,:,1,idx,0,0,0) = RCONST(0.0)
        
      end if
    end do

  end subroutine hydro_calcMatScDissMatD1d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDiss1d_sim(DdataAtEdge,&
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
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: Ei,Ej,hi,hj,pi,pj,ri,rj,ui,uj
    real(DP) :: H_ij,anorm,aux,c_ij,q_ij,u_ij,vel_ij
    integer :: idx


    do idx = 1, nedges
    
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute specific total energies
      Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      
      IDX3(DmatrixAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)

      IDX3(DmatrixAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,5,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,6,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,7,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,9,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)

      IDX3(DmatrixAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)

      IDX3(DmatrixAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,5,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,6,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,7,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,9,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                       IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))
      anorm = abs(a(1)) ! = sqrt(a(1)*a(1))
      
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
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary values
        vel_ij = u_ij*a(1)
        q_ij   = RCONST(0.5)*(u_ij*u_ij)
        
        ! Compute the speed of sound
        c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))
        
        ! Compute scalar dissipation
        aux = dscale * (abs(vel_ij) + anorm*c_ij)

        IDX3(DmatrixAtEdge,:,1,idx,0,0,0) = RCONST(0.0)
        IDX3(DmatrixAtEdge,1,1,idx,0,0,0) = aux
        IDX3(DmatrixAtEdge,5,1,idx,0,0,0) = aux
        IDX3(DmatrixAtEdge,9,1,idx,0,0,0) = aux

      else

        ! Nullify dissipation tensor
        IDX3(DmatrixAtEdge,:,1,idx,0,0,0) = RCONST(0.0)

      end if
    end do
    
  end subroutine hydro_calcMatScDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDissMatD1d_sim(DdataAtEdge,&
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
    real(DP), dimension(NVAR1D,NVAR1D) :: R_ij,L_ij
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj
    real(DP) :: H_ij,anorm,aux,b1,b2,cPow2_ij,c_ij,q_ij,u_ij
    real(DP) :: l1,l2,l3
    integer :: idx
#if defined(HYDRO_USE_ENTROPYFIX) && (HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX)
    real(DP) :: dtol,ci,cj
#endif

    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,_)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                       IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))
      anorm = abs(a(1))
      
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
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
                
        ! Compute auxiliary variable
        q_ij = RCONST(0.5)*(u_ij*u_ij)

        ! Compute speed of sound
        cPow2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij     = sqrt(cPow2_ij)
        
        b2 = ((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
        b1 = b2*q_ij

        ! Diagonal matrix of eigenvalues
        l1 = abs(u_ij-c_ij)
        l2 = abs(u_ij)
        l3 = abs(u_ij+c_ij)

#if defined(HYDRO_USE_ENTROPYFIX)

#if HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX

        ! Entropy-fix by Harten and Hyman
        ci = SOUNDSPEED3_1D(DdataAtEdge,IDX3,1,idx,0,0,0)
        cj = SOUNDSPEED3_1D(DdataAtEdge,IDX3,2,idx,0,0,0)

        dtol = max(RCONST(0.0), (u_ij-c_ij) - (ui-ci), (uj-cj) - (u_ij-c_ij) )
        
        if (l1 .lt. dtol)&
            l1 = RCONST(0.5)*((l1**2)/dtol + dtol)
        
        dtol = max(RCONST(0.0), (u_ij+c_ij) - (ui+ci), (uj+cj) - (u_ij+c_ij) )

        if (l3 .lt. dtol)&
            l3 = RCONST(0.5)*((l3**2)/dtol + dtol)

#elif HYDRO_USE_ENTROPYFIX == HARTEN_ENTROPYFIX

#ifndef HYDRO_HARTEN_ENTROPYFIX
#error "Value HYDRO_HARTEN_ENTROPYFIX is required!"
#else
        ! Entropy-fix by Harten
        if (l1 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l1 = RCONST(0.5)*((l1**2)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))

        if (l3 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l3 = RCONST(0.5)*((l3**2)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))
#endif
#else
#error "Invalid type of entropy fix!"
#endif
#endif

        ! Matrix of right eigenvectors
        R_ij(1,1) =  l1
        R_ij(2,1) =  l1*(u_ij-c_ij)
        R_ij(3,1) =  l1*(H_ij-c_ij*u_ij)
        
        R_ij(1,2) =  l2
        R_ij(2,2) =  l2*u_ij
        R_ij(3,2) =  l2*q_ij
        
        R_ij(1,3) =  l3
        R_ij(2,3) =  l3*(u_ij+c_ij)
        R_ij(3,3) =  l3*(H_ij+c_ij*u_ij)
        
        ! Matrix of left eigenvectors
        L_ij(1,1) = RCONST(0.5)*(b1+u_ij/c_ij)
        L_ij(2,1) = RCONST(1.0)-b1
        L_ij(3,1) = RCONST(0.5)*(b1-u_ij/c_ij)
        
        L_ij(1,2) = -RCONST(0.5)*(b2*u_ij+1/c_ij)
        L_ij(2,2) =  b2*u_ij
        L_ij(3,2) = -RCONST(0.5)*(b2*u_ij-1/c_ij)
        
        L_ij(1,3) = RCONST(0.5)*b2
        L_ij(2,3) = -b2
        L_ij(3,3) = RCONST(0.5)*b2

        ! Include scaling parameter
        anorm = dscale*anorm

        ! Compute tensorial dissipation D_ij = diag(R_ij*|Lbd_ij|*L_ij)*I
        IDX3(DmatrixAtEdge,:,1,idx,0,0,0) = RCONST(0.0)
        IDX3(DmatrixAtEdge,1,1,idx,0,0,0) = anorm*( R_ij(1,1)*L_ij(1,1)+&
                                                          R_ij(1,2)*L_ij(2,1)+&
                                                          R_ij(1,3)*L_ij(3,1)  )
        IDX3(DmatrixAtEdge,2,1,idx,0,0,0) = anorm*( R_ij(2,1)*L_ij(1,2)+&
                                                          R_ij(2,2)*L_ij(2,2)+&
                                                          R_ij(2,3)*L_ij(3,2)  )
        IDX3(DmatrixAtEdge,3,1,idx,0,0,0) = anorm*( R_ij(3,1)*L_ij(1,3)+&
                                                          R_ij(3,2)*L_ij(2,3)+&
                                                          R_ij(3,3)*L_ij(3,3)  )
      else
        
        ! Nullify dissipation tensor
        IDX3(DmatrixAtEdge,:,1,idx,0,0,0) = RCONST(0.0)

      end if
    end do

  end subroutine hydro_calcMatRoeDissMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDiss1d_sim(DdataAtEdge,&
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
    real(DP), dimension(NVAR1D,NVAR1D) :: R_ij,L_ij
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: Ei,Ej,hi,hj,pi,pj,ri,rj,ui,uj
    real(DP) :: H_ij,anorm,aux,b1,b2,cPow2_ij,c_ij,q_ij,u_ij
    real(DP) :: l1,l2,l3
    integer :: idx,i,j,k
#if defined(HYDRO_USE_ENTROPYFIX) && (HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX)
    real(DP) :: dtol,ci,cj
#endif

    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute specific total energies
      Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      
      IDX3(DmatrixAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)

      IDX3(DmatrixAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,5,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,6,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,7,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,9,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)

      IDX3(DmatrixAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)

      IDX3(DmatrixAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,5,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,6,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,7,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,9,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = RCONST(0.5)*(DcoeffsAtEdge(1,1,idx)-DcoeffsAtEdge(1,2,idx))
      anorm = abs(a(1))

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
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
                
        ! Compute auxiliary values
        q_ij = RCONST(0.5)*u_ij*u_ij
        cPow2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij     = sqrt(cPow2_ij)

        b2 = ((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
        b1 = b2*q_ij
        
        ! Diagonal matrix of eigenvalues
        l1 = abs(u_ij-c_ij)
        l2 = abs(u_ij)
        l3 = abs(u_ij+c_ij)
        
#if defined(HYDRO_USE_ENTROPYFIX)

#if HYDRO_USE_ENTROPYFIX == HARTEN_HYMAN_ENTROPYFIX

        ! Entropy-fix by Harten and Hyman
        ci = SOUNDSPEED3_1D(DdataAtEdge,IDX3,1,idx,0,0,0)
        cj = SOUNDSPEED3_1D(DdataAtEdge,IDX3,2,idx,0,0,0)

        dtol = max(RCONST(0.0), (u_ij-c_ij) - (ui-ci), (uj-cj) - (u_ij-c_ij) )
        
        if (l1 .lt. dtol)&
            l1 = RCONST(0.5)*((l1**2)/dtol + dtol)
        
        dtol = max(RCONST(0.0), (u_ij+c_ij) - (ui+ci), (uj+cj) - (u_ij+c_ij) )

        if (l3 .lt. dtol)&
            l3 = RCONST(0.5)*((l3**2)/dtol + dtol)

#elif HYDRO_USE_ENTROPYFIX == HARTEN_ENTROPYFIX

#ifndef HYDRO_HARTEN_ENTROPYFIX
#error "Value HYDRO_HARTEN_ENTROPYFIX is required!"
#else
        ! Entropy-fix by Harten
        if (l1 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l1 = RCONST(0.5)*((l1**2)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))

        if (l3 .lt. RCONST(HYDRO_HARTEN_ENTROPYFIX))&
            l3 = RCONST(0.5)*((l3**2)/RCONST(HYDRO_HARTEN_ENTROPYFIX)&
               + RCONST(HYDRO_HARTEN_ENTROPYFIX))
#endif
#else
#error "Invalid type of entropy fix!"
#endif
#endif

        ! Matrix of right eigenvectors
        R_ij(1,1) =  l1
        R_ij(2,1) =  l1*(u_ij-c_ij)
        R_ij(3,1) =  l1*(H_ij-c_ij*u_ij)
        
        R_ij(1,2) =  l2
        R_ij(2,2) =  l2*u_ij
        R_ij(3,2) =  l2*q_ij
        
        R_ij(1,3) =  l3
        R_ij(2,3) =  l3*(u_ij+c_ij)
        R_ij(3,3) =  l3*(H_ij+c_ij*u_ij)
        
        ! Matrix of left eigenvectors
        L_ij(1,1) = RCONST(0.5)*(b1+u_ij/c_ij)
        L_ij(2,1) = RCONST(1.0)-b1
        L_ij(3,1) = RCONST(0.5)*(b1-u_ij/c_ij)
        
        L_ij(1,2) = -RCONST(0.5)*(b2*u_ij+1/c_ij)
        L_ij(2,2) =  b2*u_ij
        L_ij(3,2) = -RCONST(0.5)*(b2*u_ij-1/c_ij)
        
        L_ij(1,3) = RCONST(0.5)*b2
        L_ij(2,3) = -b2
        L_ij(3,3) = RCONST(0.5)*b2
        
        ! Include scaling parameter
        anorm = dscale*anorm

        ! Compute tensorial dissipation D_ij = R_ij*|Lbd_ij|*L_ij
        do i = 1, NVAR1D
          do j = 1, NVAR1D
            aux = RCONST(0.0)
            do k = 1, NVAR1D
              aux = aux + R_ij(i,k)*L_ij(k,j)
            end do
            IDX3(DmatrixAtEdge,NVAR1D*(j-1)+i,1,idx,0,0,0) = anorm*aux
          end do
        end do
        
      else
        
        ! Nullify dissipation tensor
        IDX3(DmatrixAtEdge,:,1,idx,0,0,0) = RCONST(0.0)
        
      end if
    end do

  end subroutine hydro_calcMatRoeDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDissMatD1d_sim(DdataAtEdge,&
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
    real(DP) :: Ei,Ej,ci,cj,ui,uj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute specific total energies
      Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,_)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,_)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,_)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,_)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,_)
      IDX3(DmatrixAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,_)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate scalar artificial dissipation of Rusanov-type
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*&
          (HYDRO_GAMMA)*(Ei-RCONST(0.5)*ui*ui), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*&
          (HYDRO_GAMMA)*(Ej-RCONST(0.5)*uj*uj), SYS_EPSREAL_DP))
      
      ! Compute dissipation tensor
      IDX3(DmatrixAtEdge,:,1,idx,0,0,0) = dscale *&
          max( abs(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*uj)+&
               abs(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0))*cj,&
               abs(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*ui)+&
               abs(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*ci )
    end do
    
  end subroutine hydro_calcMatRusDissMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDiss1d_sim(DdataAtEdge,&
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
    real(DP) :: Ei,Ej,aux,ci,cj,ui,uj
    integer :: idx


    do idx = 1, nedges

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Compute specific total energies
      Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      
      IDX3(DmatrixAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)

      IDX3(DmatrixAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      IDX3(DmatrixAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,5,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,6,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,7,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,9,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),ui,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      IDX3(DmatrixAtEdge,1,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,2,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,3,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)

      IDX3(DmatrixAtEdge,4,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,5,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,6,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)

      IDX3(DmatrixAtEdge,7,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,8,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      IDX3(DmatrixAtEdge,9,2,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,1,idx,0,0,0),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      IDX3(DmatrixAtEdge,1,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX11(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,2,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX21(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,3,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX31(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,4,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX12(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,5,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX22(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,6,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX32(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)

      IDX3(DmatrixAtEdge,7,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX13(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,8,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX23(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
      IDX3(DmatrixAtEdge,9,3,idx,0,0,0) =&
          FLUXJACOBIMATRIX33(-dscale,IDX3(DcoeffsAtEdge,1,2,idx,0,0,0),ui,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate scalar artificial dissipation of Rusanov-type
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*&
          (HYDRO_GAMMA)*(Ei-RCONST(0.5)*ui*ui), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*&
          (HYDRO_GAMMA)*(Ej-RCONST(0.5)*uj*uj), SYS_EPSREAL_DP))

      ! Compute dissipation tensor
      aux = dscale * max( abs(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*uj)+&
                          abs(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0))*cj,&
                          abs(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*ui)+&
                          abs(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*ci )

      IDX3(DmatrixAtEdge,:,1,idx,0,0,0) = RCONST(0.0)
      IDX3(DmatrixAtEdge,1,1,idx,0,0,0) = aux
      IDX3(DmatrixAtEdge,5,1,idx,0,0,0) = aux
      IDX3(DmatrixAtEdge,9,1,idx,0,0,0) = aux
    end do
    
  end subroutine hydro_calcMatRusDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcCharacteristics1d_sim(Dweight, DdataAtEdge,&
      nedges, DcharVariablesAtEdge, DeigenvaluesAtEdge,&
      DrightEigenvectorsAtEdge, DleftEigenvectorsAtEdge, rcollection)

!<description>
    ! This subroutine computes the characteristic variables in 1D
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

    ! local variables
    real(DP), dimension(NVAR1D) :: Diff
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj
    real(DP) :: H_ij,anorm,aux,b1_ij,b2_ij,cPow2_ij,c_ij,q_ij,u_ij
    integer :: idx


    ! Compute norm of weighting coefficient
    anorm = abs(Dweight(1))

    ! Check if weighting coefficient is zero
    if (anorm .le. SYS_EPSREAL_DP) then
      if (present(DcharVariablesAtEdge))     DcharVariablesAtEdge     = RCONST(0.0)
      if (present(DeigenvaluesAtEdge))       DeigenvaluesAtEdge       = RCONST(0.0)
      if (present(DrightEigenvectorsAtEdge)) DrightEigenvectorsAtEdge = RCONST(0.0)
      if (present(DleftEigenvectorsAtEdge))  DleftEigenvectorsAtEdge  = RCONST(0.0)

      ! That's it
      return
    end if
    
    
    ! Do we have to compute characteristic variables
    if (present(DcharVariablesAtEdge)) then
      do idx = 1, nedges
        
        ! Compute velocities
        ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
        
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
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)

        ! Compute auxiliary variable
        q_ij = RCONST(0.5)*(u_ij*u_ij)

        ! Compute speed of sound
        cPow2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij = sqrt(cPow2_ij)

        b2_ij = ((HYDRO_GAMMA)-1.0)/cPow2_ij
        b1_ij  = b2_ij*q_ij
        
        ! Compute solution difference U_j-U_i
        Diff = IDX3(DdataAtEdge,:,2,idx,0,0,0)-&
               IDX3(DdataAtEdge,:,1,idx,0,0,0)
          
        ! Compute characteristic variables
        IDX2(DcharVariablesAtEdge,1,idx,0,0) = anorm * RCONST(0.5) *&
            (            (b1_ij+u_ij/c_ij)*Diff(1)-&
             (b2_ij*u_ij+RCONST(1.0)/c_ij)*Diff(2)+&
                                     b2_ij*Diff(3))
        IDX2(DcharVariablesAtEdge,2,idx,0,0) = anorm *&
            (          (RCONST(1.0)-b1_ij)*Diff(1)+&
                                b2_ij*u_ij*Diff(2)-&
                                     b2_ij*Diff(3) )
        IDX2(DcharVariablesAtEdge,3,idx,0,0) = anorm * RCONST(0.5) *&
            (            (b1_ij-u_ij/c_ij)*Diff(1)-&
             (b2_ij*u_ij-RCONST(1.0)/c_ij)*Diff(2)+&
                                     b2_ij*Diff(3) )
      end do
    end if


    ! Do we have to compute eigenvalues
    if (present(DeigenvaluesAtEdge)) then
      do idx = 1, nedges
        
        ! Compute velocities
        ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

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
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variable
        c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*&
            (H_ij-RCONST(0.5)*(u_ij*u_ij)), SYS_EPSREAL_DP))

        ! Compute eigenvalues
        IDX2(DeigenvaluesAtEdge,1,idx,0,0) = u_ij-c_ij
        IDX2(DeigenvaluesAtEdge,2,idx,0,0) = u_ij
        IDX2(DeigenvaluesAtEdge,3,idx,0,0) = u_ij+c_ij
      end do
    end if


    ! Do we have to compute right eigenvectors
    if (present(DrightEigenvectorsAtEdge)) then
      do idx = 1, nedges
        
        ! Compute velocities
        ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

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
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        q_ij = RCONST(0.5)*(u_ij*u_ij)
        c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))

        ! Compute right eigenvectors
        IDX2(DrightEigenvectorsAtEdge,1,idx,0,0) =  RCONST(1.0)
        IDX2(DrightEigenvectorsAtEdge,2,idx,0,0) =  u_ij-c_ij
        IDX2(DrightEigenvectorsAtEdge,3,idx,0,0) =  H_ij-u_ij*c_ij

        IDX2(DrightEigenvectorsAtEdge,4,idx,0,0) =  RCONST(1.0)
        IDX2(DrightEigenvectorsAtEdge,5,idx,0,0) =  u_ij
        IDX2(DrightEigenvectorsAtEdge,6,idx,0,0) =  q_ij

        IDX2(DrightEigenvectorsAtEdge,7,idx,0,0) =  RCONST(1.0)
        IDX2(DrightEigenvectorsAtEdge,8,idx,0,0) =  u_ij+c_ij
        IDX2(DrightEigenvectorsAtEdge,9,idx,0,0) =  H_ij+u_ij*c_ij
      end do
    end if


    ! Do we have to compute left eigenvectors
    if (present(DleftEigenvectorsAtEdge)) then
      do idx = 1, nedges
        
        ! Compute velocities
        ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

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
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
        
        ! Compute auxiliary variables
        q_ij = RCONST(0.5)*(u_ij*u_ij)
        cPow2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij = sqrt(cPow2_ij)

        b2_ij = ((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij
        b1_ij = b2_ij*q_ij

        ! Compute left eigenvectors
        IDX2(DleftEigenvectorsAtEdge,1,idx,0,0) = RCONST(0.5) * (b1_ij+u_ij/c_ij)
        IDX2(DleftEigenvectorsAtEdge,2,idx,0,0) =               (RCONST(1.0)-b1_ij)
        IDX2(DleftEigenvectorsAtEdge,3,idx,0,0) = RCONST(0.5) * (b1_ij-u_ij/c_ij)

        IDX2(DleftEigenvectorsAtEdge,4,idx,0,0) =-RCONST(0.5) * (b2_ij*u_ij+RCONST(1.0)/c_ij)
        IDX2(DleftEigenvectorsAtEdge,5,idx,0,0) =               (b2_ij*u_ij)
        IDX2(DleftEigenvectorsAtEdge,6,idx,0,0) =-RCONST(0.5) * (b2_ij*u_ij-RCONST(1.0)/c_ij)

        IDX2(DleftEigenvectorsAtEdge,7,idx,0,0) = RCONST(0.5)*b2_ij
        IDX2(DleftEigenvectorsAtEdge,8,idx,0,0) =            -b2_ij
        IDX2(DleftEigenvectorsAtEdge,9,idx,0,0) = RCONST(0.5)*b2_ij
      end do
    end if
    
  end subroutine hydro_calcCharacteristics1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTScDiss1d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj
    real(DP) :: H_ij,aux,c_ij,d_ij,q_ij,u_ij,vel_ij
    integer :: idx


    do idx = 1, nedges

      ! Compute skew-symmetric coefficient and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                       IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))

      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
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
      H_ij = ROE_MEAN_VALUE(hi,hj,aux)

      ! Compute auxiliary variables
      vel_ij = u_ij*a(1)
      q_ij   = RCONST(0.5)*(u_ij*u_ij)

      ! Compute the speed of sound
      c_ij = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP))
      
      ! Compute scalar dissipation
      d_ij = abs(vel_ij) + abs(a(1))*c_ij

      ! Compute conservative fluxes
      IDX2(DfluxesAtEdge,:,idx,0,0) = dscale*&
          d_ij*(IDX3(DdataAtEdge,:,1,idx,0,0,0)-&
                IDX3(DdataAtEdge,:,2,idx,0,0,0))
    end do

  end subroutine hydro_calcFluxFCTScDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRoeDiss1d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
    real(DP) :: hi,hj,pi,pj,ri,rj,ui,uj
    real(DP) :: H_ij,anorm,aux,b1,b2,cPow2_ij,c_ij,q_ij,u_ij
    real(DP) :: l1,l2,l3,w1,w2,w3
    integer :: idx

    do idx = 1, nedges

      ! Compute skew-symmetric coefficient and its norm
      a = RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                       IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))
      anorm = abs(a(1)) ! = sqrt(a(1)*a(1))

      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Compute velocities
        ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
        uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

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
        H_ij = ROE_MEAN_VALUE(hi,hj,aux)
         
        ! Compute auxiliary variables
        q_ij = RCONST(0.5)*(u_ij*u_ij)

        ! Compute the speed of sound
        cPow2_ij = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_ij-q_ij), SYS_EPSREAL_DP)
        c_ij     = sqrt(cPow2_ij)
        
        ! Compute eigenvalues
        l1 = abs(u_ij-c_ij)
        l2 = abs(u_ij)
        l3 = abs(u_ij+c_ij)
        
        ! Compute solution difference U_i-U_j
        Diff = IDX3(DdataAtEdge,:,1,idx,0,0,0)-&
               IDX3(DdataAtEdge,:,2,idx,0,0,0)
        
        ! Compute auxiliary quantities for characteristic variables
        b2 = ((HYDRO_GAMMA)-RCONST(1.0))/cPow2_ij; b1 = b2*q_ij
        
        ! Compute characteristic variables multiplied by the
        ! corresponding eigenvalue
        w1 = l1 * RCONST(0.5) * (       (b1+u_ij/c_ij)*Diff(1) -&
                            (b2*u_ij+RCONST(1.0)/c_ij)*Diff(2) +&
                                                    b2*Diff(3) )
        w2 = l2 * (                   (RCONST(1.0)-b1)*Diff(1) +&
                                               b2*u_ij*Diff(2) -&
                                                    b2*Diff(3) )
        w3 = l3 * RCONST(0.5) * (       (b1-u_ij/c_ij)*Diff(1) -&
                            (b2*u_ij-RCONST(1.0)/c_ij)*Diff(2) +&
                                                    b2*Diff(3) )
        
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        Diff(1) = anorm * ( w1 + w2 + w3 )
        Diff(2) = anorm * ( (u_ij-c_ij)*w1 + u_ij*w2 + (u_ij+c_ij)*w3 )
        Diff(3) = anorm * ( (H_ij-u_ij*c_ij)*w1 + q_ij*w2 + (H_ij+u_ij*c_ij)*w3 )

        ! Compute conservative fluxes
        IDX2(DfluxesAtEdge,:,idx,0,0) = dscale*Diff
      else
        ! Ccalcel conservative fluxes
        IDX2(DfluxesAtEdge,:,idx,0,0) = 0
      end if
    end do
    
  end subroutine hydro_calcFluxFCTRoeDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRusDiss1d_sim(DdataAtEdge, DcoeffsAtEdge,&
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
    real(DP) :: Ei,Ej,ci,cj,ui,uj
    real(DP) :: d_ij
    integer :: idx

    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute specific total energies
      Ei = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      Ej = SPECIFICTOTALENERGY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Compute the speed of sound
      ci = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*&
          (HYDRO_GAMMA)*(Ei-RCONST(0.5)*ui*ui), SYS_EPSREAL_DP))
      cj = sqrt(max(((HYDRO_GAMMA)-RCONST(1.0))*&
          (HYDRO_GAMMA)*(Ej-RCONST(0.5)*uj*uj), SYS_EPSREAL_DP))
     
#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*uj)+&
                 RCONST(0.5)*sqrt((IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))**2)*cj,&
                  abs(RCONST(0.5)*(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,0,0,0))*ui)+&
                 RCONST(0.5)*sqrt((IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)-&
                                   IDX3(DcoeffsAtEdge,1,1,idx,0,0,0))**2)*ci )
#else
      ! Compute scalar dissipation
      d_ij = max( abs(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0)*uj)+&
                  abs(IDX3(DcoeffsAtEdge,1,1,idx,0,0,0))*cj,&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0)*ui)+&
                  abs(IDX3(DcoeffsAtEdge,1,2,idx,0,0,0))*ci )
#endif
      
      ! Compute conservative fluxes
      IDX2(DfluxesAtEdge,:,idx,0,0) = dscale*&
          d_ij*(IDX3(DdataAtEdge,:,1,idx,0,0,0)-&
                IDX3(DdataAtEdge,:,2,idx,0,0,0))
    end do

  end subroutine hydro_calcFluxFCTRusDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDensity1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density in 1D
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

  end subroutine hydro_trafoFluxDensity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDensity1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density in 1D
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

  end subroutine hydro_trafoDiffDensity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalDensity1d_sim(DdataAtNode,&
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

  end subroutine hydro_trafoNodalDensity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxEnergy1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the energy in 1D
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

  end subroutine hydro_trafoFluxEnergy1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffEnergy1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the energy in 1D
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

  end subroutine hydro_trafoDiffEnergy1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalEnergy1d_sim(DdataAtNode,&
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

  end subroutine hydro_trafoNodalEnergy1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxPressure1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the pressure in 1D
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
    real(DP) :: ui,uj
    integer :: idx
    
    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          ((HYDRO_GAMMA)-RCONST(1.0))*(RCONST(0.5)*&
          ui*ui*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
           ui*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)+&
            TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
          -((HYDRO_GAMMA)-RCONST(1.0))*(RCONST(0.5)*&
          uj*uj*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
           uj*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)+&
            TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
    end do
    
  end subroutine hydro_trafoFluxPressure1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffPressure1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the pressure in 1D
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

  end subroutine hydro_trafoDiffPressure1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalPressure1d_sim(DdataAtNode,&
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

  end subroutine hydro_trafoNodalPressure1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxVelocity1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the x-velocity
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
    real(DP) :: ui,uj

    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)

      ! Transformed velocity fluxes in x-direction
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          (XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
          ui*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0))/&
             DENSITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -(XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)-&
          uj*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0))/&
             DENSITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
    end do
    
  end subroutine hydro_trafoFluxVelocity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffVelocity1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the x-velocity
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
    end do

  end subroutine hydro_trafoDiffVelocity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalVelocity1d_sim(DdataAtNode,&
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
    end do

  end subroutine hydro_trafoNodalVelocity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxMomentum1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the x-momentum
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
    end do
    
  end subroutine hydro_trafoFluxMomentum1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffMomentum1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the x-momentum
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
    end do
    
  end subroutine hydro_trafoDiffMomentum1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalMomentum1d_sim(DdataAtNode,&
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
    end do

  end subroutine hydro_trafoNodalMomentum1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenEng1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density and energy in 1D
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

  end subroutine hydro_trafoFluxDenEng1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenEng1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density and energy in 1D
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

  end subroutine hydro_trafoDiffDenEng1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalDenEng1d_sim(DdataAtNode,&
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

  end subroutine hydro_trafoNodalDenEng1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenPre1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density and energy in 1D
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
    real(DP) :: ui,uj
    integer :: idx

    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
           
      ! Transformed density fluxes
      IDX3(DtransformedFluxesAtEdge,1,1,idx,0,0,0) =&
          DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)
      IDX3(DtransformedFluxesAtEdge,1,2,idx,0,0,0) =&
         -DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)

      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,2,1,idx,0,0,0) =&
          ((HYDRO_GAMMA)-RCONST(1.0))*(RCONST(0.5)*&
          ui*ui*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
           ui*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)+&
            TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
      IDX3(DtransformedFluxesAtEdge,2,2,idx,0,0,0) =&
          -((HYDRO_GAMMA)-RCONST(1.0))*(RCONST(0.5)*&
          uj*uj*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
           uj*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)+&
            TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
    end do

  end subroutine hydro_trafoFluxDenPre1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenPre1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density and energy in 1D
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

  end subroutine hydro_trafoDiffDenPre1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalDenPre1d_sim(DdataAtNode,&
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

  end subroutine hydro_trafoNodalDenPre1d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenPreVel1d_sim(DdataAtEdge,&
      DfluxesAtEdge, nedges, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to fluxes for the density, pressure and velocity in 1D
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
    real(DP) :: ui,uj
    integer :: idx

    do idx = 1, nedges
      
      ! Compute velocities
      ui = XVELOCITY3(DdataAtEdge,IDX3,1,idx,0,0,0)
      uj = XVELOCITY3(DdataAtEdge,IDX3,2,idx,0,0,0)
      
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

      ! Transformed pressure fluxes
      IDX3(DtransformedFluxesAtEdge,3,1,idx,0,0,0) =&
          ((HYDRO_GAMMA)-RCONST(1.0))*(RCONST(0.5)*&
          ui*ui*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
           ui*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)+&
            TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
      IDX3(DtransformedFluxesAtEdge,3,2,idx,0,0,0) =&
          -((HYDRO_GAMMA)-RCONST(1.0))*(RCONST(0.5)*&
          uj*uj*DENSITY2(DfluxesAtEdge,IDX2,idx,0,0)-&
           uj*XMOMENTUM2(DfluxesAtEdge,IDX2,idx,0,0)+&
            TOTALENERGY2(DfluxesAtEdge,IDX2,idx,0,0))
    end do
      
  end subroutine hydro_trafoFluxDenPreVel1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenPreVel1d_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density, pressure and velocity in 1D
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
      
      ! Transformed pressure difference
      IDX2(DtransformedDataAtEdge,3,idx,0,0) =&
          PRESSURE3(DdataAtEdge,IDX3,2,idx,0,0,0)-&
          PRESSURE3(DdataAtEdge,IDX3,1,idx,0,0,0)
    end do

  end subroutine hydro_trafoDiffDenPreVel1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoNodalDenPreVel1d_sim(DdataAtNode,&
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

      ! Transformed pressure values
      IDX2(DtransformedDataAtNode,3,idx,0,0) =&
          PRESSURE2(DdataAtNode,IDX2,idx,0,0)
    end do

  end subroutine hydro_trafoNodalDenPreVel1d_sim

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcBoundaryvalues1d(DbdrNormal, DpointNormal,&
      DbdrValue, ibdrCondType, Du, Du0, istatus)

!<description>
    ! This subroutine computes the boundary values for a given node in 1D
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
    real(DP), dimension(NVAR1D) :: W,Wu,Winf    ! Riemann invariants, eigenvalues, etc.
    real(DP) :: rho,v1,p,E,c                    ! primitive variables
    real(DP) :: vn_b,vn,pstar,ps                ! velocities and boundary values
    real(DP) :: cup,f,fd,ge,qrt                 ! auxiliary variables ...
    real(DP) :: pold,ppv,prat,ptl,ptr,vdiff,vm  ! ... for the Riemann solver
    real(DP) :: auxA,auxB,aux
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
      ! See the 2D-version of this routine for details.

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      E   = Du(3)/rho
      p   = ((HYDRO_GAMMA)-RCONST(1.0))*rho*(E-RCONST(0.5)*v1*v1)
      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))

      ! Compute normal velocities at the boundary and the ghost state
      ! w.r.t. the numerical/approximate  outward unit normal vector
      vn   =  DbdrNormal(1)*v1
      vn_b = -DbdrNormal(1)*v1

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
      if (RCONST(4.0)/((HYDRO_GAMMA)-RCONST(1.0))*c .le. vn_b-vn) then
        if (present(istatus)) then
          istatus = -ibdrCondType
          return
        else
          call output_line('Riemann solver failed due to vacuum',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcBoundaryvalues1d')
          call sys_halt()
        end if
      end if

      ! Provide a guess value for pressure in the "star region"
      ! by using the PVRS Riemann solver as suggested by Toro

      cup  = rho*c
      ppv  = p+RCONST(0.5)*(vn-vn_b)*cup
      ppv  = max(RCONST(0.0), ppv)

      if (ppv .eq. p) then

        ! Select guessed pressure from PVRS Riemann solver
        pstar = ppv
      else
        if (ppv .lt. p) then

          ! Guess pressure from the Two-Rarefaction Riemann solver
          vm    = RCONST(0.5)*(vn+vn_b)
          ptl   = RCONST(1.0)+((HYDRO_GAMMA)-RCONST(1.0))/RCONST(2.0)*(vn-vm)/c
          ptr   = RCONST(1.0)+((HYDRO_GAMMA)-RCONST(1.0))/RCONST(2.0)*(vm-vn_b)/c
          pstar = RCONST(0.5)*(p*ptl + p*ptr)**(RCONST(2.0)*(HYDRO_GAMMA)/((HYDRO_GAMMA)-RCONST(1.0)))
        else

          ! Guess pressure from the Two-Shock Riemann solver
          ! with PVRS as estimated pressure value
          ge    = sqrt((RCONST(2.0)/((HYDRO_GAMMA)+RCONST(1.0))/rho)/&
              (((HYDRO_GAMMA)-RCONST(1.0))/((HYDRO_GAMMA)+RCONST(1.0))*p+ppv))
          pstar = p - RCONST(0.5)*(vn_b-vn)/ge
        end if
      end if

      ! Initialize solution difference and pressure
      vdiff = (vn_b-vn)/RCONST(2.0)
      pold  = pstar

      newton: do ite = 1, 100

        ! Compute pressure function f(pold) and its derivative f1(pold)
        if (pold .le. p) then

          ! Rarefaction wave
          prat = pold/p

          f  = RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0))*c*&
              (prat**(((HYDRO_GAMMA)-RCONST(1.0))/(RCONST(2.0)*(HYDRO_GAMMA)))-RCONST(1.0))
          fd = (RCONST(1.0)/(rho*c))*prat**(-((HYDRO_GAMMA)+RCONST(1.0))/(RCONST(2.0)*(HYDRO_GAMMA)))
        else

          ! Shock wave
          auxA = RCONST(2.0)/((HYDRO_GAMMA)+RCONST(1.0))/rho
          auxB = ((HYDRO_GAMMA)-RCONST(1.0))/((HYDRO_GAMMA)+RCONST(1.0))*p
          qrt  = sqrt(auxA/(auxB + pold))

          f  = (pold-p)*qrt
          fd = (RCONST(1.0) - RCONST(0.5)*(pold - p)/(auxB + pold))*qrt
        end if

        pstar = pold - (f+vdiff)/fd
        if (pstar .lt. RCONST(0.0)) then
          pold = 1.0E-6
          cycle newton
        end if

        aux = RCONST(2.0)*abs((pstar-pold)/(pstar+pold))
        if (aux .le. RCONST(1.0E-6))  exit newton

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
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcBoundaryvalues1d')
          call sys_halt()
        end if
      end if

      !-------------------------------------------------------------------------
      ! Calculate the velocity in the star region
      !-------------------------------------------------------------------------

      ! Note that the contribution fR-fL vanishes due to constant states
      vn = RCONST(0.5)*(vn+vn_b)


      !-------------------------------------------------------------------------
      ! Calculate the density in the star region
      !-------------------------------------------------------------------------

      if (pstar .le. p) then

        ! Rarefaction wave
        rho = rho*(pstar/p)**(RCONST(1.0)/(HYDRO_GAMMA))
      else

        ! Shock wave
        rho = rho*(pstar/p+((HYDRO_GAMMA)-RCONST(1.0))/((HYDRO_GAMMA)+RCONST(1.0)))/&
                  (((HYDRO_GAMMA)-RCONST(1.0))/((HYDRO_GAMMA)+RCONST(1.0))*(pstar/p)+RCONST(1.0))
      end if

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*DpointNormal(1)*vn
      Du(3) = pstar/((HYDRO_GAMMA)-RCONST(1.0))+RCONST(0.5)*rho*(vn*vn)

    case(BDRC_VISCOUSWALL)
      !-------------------------------------------------------------------------

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      E   = Du(3)/rho
      p   = ((HYDRO_GAMMA)-RCONST(1.0))*rho*(E-RCONST(0.5)*v1*v1)

      ! Update the solution vector and let vn:=0
      Du(2) = RCONST(0.0)
      Du(3) = p/((HYDRO_GAMMA)-RCONST(1.0))

    case(BDRC_SUPERINLET)
      !-------------------------------------------------------------------------

      ! The free stream primitive variables are Deval=[rho,v1,p]
      rho = DbdrValue(1)
      p   = DbdrValue(3)
      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))
      vn  = DbdrNormal(1)*DbdrValue(2)

      ! Compute Riemann invariants based on the free stream values
      W(1) = vn-2*c/((HYDRO_GAMMA)-RCONST(1.0))
      W(2) = vn+2*c/((HYDRO_GAMMA)-RCONST(1.0))
      W(3) = p/(rho**(HYDRO_GAMMA))

      ! Transform back into conservative variables
      vn   = RCONST(0.5)*(W(1)+W(2))
      c    = RCONST(0.25)*((HYDRO_GAMMA)-RCONST(1.0))*(W(2)-W(1))
      rho  = (c*c/(HYDRO_GAMMA)/W(3))**(RCONST(1.0)/((HYDRO_GAMMA)-RCONST(1.0)))
      p    = rho*c*c/(HYDRO_GAMMA)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*DbdrNormal(1)*vn
      Du(3) = p/((HYDRO_GAMMA)-RCONST(1.0))+RCONST(0.5)*rho*(vn*vn)


    case(BDRC_FREESTREAM)
      !-------------------------------------------------------------------------

      ! The free stream primitive variables are Deval=[rho,v1,p]
      rho = DbdrValue(1)
      p   = DbdrValue(3)
      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))
      vn  = DbdrNormal(1)*DbdrValue(2)

      ! Compute Riemann invariants based on the free stream values
      Winf(1) = vn-RCONST(2.0)*c/((HYDRO_GAMMA)-RCONST(1.0))
      Winf(2) = vn+RCONST(2.0)*c/((HYDRO_GAMMA)-RCONST(1.0))
      Winf(3) = p/(rho**(HYDRO_GAMMA))

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      E   = Du(3)/rho
      p   = ((HYDRO_GAMMA)-RCONST(1.0))*rho*(E-RCONST(0.5)*v1*v1)

      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))
      vn  = DbdrNormal(1)*v1

      ! Compute Riemann invariants based on the solution values
      Wu(1) = vn-RCONST(2.0)*c/((HYDRO_GAMMA)-RCONST(1.0))
      Wu(2) = vn+RCONST(2.0)*c/((HYDRO_GAMMA)-RCONST(1.0))
      Wu(3) = p/(rho**(HYDRO_GAMMA))

      ! Adopt free stream/computed values depending on the sign of the eigenvalue
      W(1) = merge(Winf(1), Wu(1), vn <  c)
      W(2) = merge(Winf(2), Wu(2), vn < -c)
      W(3) = merge(Winf(3), Wu(3), vn <  SYS_EPSREAL_DP)

      ! Transform back into conservative variables
      vn   = RCONST(0.5)*(W(1)+W(2))
      c    = RCONST(0.25)*((HYDRO_GAMMA)-RCONST(1.0))*(W(2)-W(1))
      rho  = (c*c/(HYDRO_GAMMA)/W(3))**(RCONST(1.0)/((HYDRO_GAMMA)-RCONST(1.0)))
      p    = rho*c*c/(HYDRO_GAMMA)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*DbdrNormal(1)*vn
      Du(3) = p/((HYDRO_GAMMA)-RCONST(1.0))+RCONST(0.5)*rho*(vn*vn)


    case(BDRC_SUBINLET)
      !-------------------------------------------------------------------------

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      E   = Du(3)/rho
      p   = ((HYDRO_GAMMA)-RCONST(1.0))*rho*(E-RCONST(0.5)*v1*v1)

      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))
      vn  = DbdrNormal(1)*v1

      ! The specified density and pressure is Deval=[rho,p]
      rho = DbdrValue(1)
      p   = DbdrValue(2)

      ! Compute Riemann invariants
      W(1) = vn-RCONST(2.0)*c/((HYDRO_GAMMA)-RCONST(1.0))
      W(2) = vn+RCONST(2.0)*c/((HYDRO_GAMMA)-RCONST(1.0))
      W(3) = p/(rho**(HYDRO_GAMMA))

      ! Transform back into conservative variables
      vn   = RCONST(0.5)*(W(1)+W(2))
      c    = RCONST(0.25)*((HYDRO_GAMMA)-RCONST(1.0))*(W(2)-W(1))
      rho  = (c*c/(HYDRO_GAMMA)/W(3))**(RCONST(1.0)/((HYDRO_GAMMA)-RCONST(1.0)))
      p    = rho*c*c/(HYDRO_GAMMA)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*DbdrNormal(1)*vn
      Du(3) = p/((HYDRO_GAMMA)-RCONST(1.0))+RCONST(0.5)*rho*(vn*vn)


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
      E   = Du(3)/rho
      p   = ((HYDRO_GAMMA)-RCONST(1.0))*rho*(E-RCONST(0.5)*v1*v1)

      vn  = DbdrNormal(1)*v1
      c   = sqrt(max((HYDRO_GAMMA)*p/rho, SYS_EPSREAL_DP))

      ! Compute Riemann invariants based on the solution values and prescribed exit pressure
      W(2) = 2*c/((HYDRO_GAMMA)-RCONST(1.0))-vn
      W(3) = p/(rho**(HYDRO_GAMMA))
      W(1) = 4/((HYDRO_GAMMA)-RCONST(1.0))*sqrt(max((HYDRO_GAMMA)*ps/&
          rho*(p/ps)**(RCONST(1.0)/(HYDRO_GAMMA)), SYS_EPSREAL_DP))-W(2)

      ! Transform back into conservative variables
      vn  = RCONST(0.5)*(W(1)-W(2))
      c   = RCONST(0.25)*((HYDRO_GAMMA)-RCONST(1.0))*(W(1)+W(2))
      rho = (c*c/(HYDRO_GAMMA)/W(3))**(RCONST(1.0)/((HYDRO_GAMMA)-RCONST(1.0)))
      p   = rho*c*c/(HYDRO_GAMMA)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*DbdrNormal(1)*vn
      Du(3) = p/((HYDRO_GAMMA)-RCONST(1.0))+RCONST(0.5)*rho*(vn*vn)


    case default
      call output_line('Unsupported type of boundary condition!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcBoundaryvalues1d')
      call sys_halt()
    end select

  end subroutine hydro_calcBoundaryvalues1d
  
  ! ***************************************************************************

!<subroutine>

  subroutine hydro_coeffVectorBdr1d_sim(rdiscretisation, rform,&
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

#ifndef HYDRO_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DHYDRO_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'hydro_coeffVectorBdr1d_sim')
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
      call fevl_evaluate_sim(DER_FUNC1D, Daux1, p_rsolution%RvectorBlock(1),&
          Dpoints, rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! What type of boundary conditions are we?
      select case(iand(ibdrtype, BDRC_TYPEMASK))
        
      case (BDRC_FREESTREAM)
        !-----------------------------------------------------------------------
        ! Free-stream boundary conditions:
        !
        ! Compute the Riemann invariants based on the computed (internal)
        ! state vector and the given freestream state vector and select
        ! the Riemman invariant for each characteristic fields based on
        ! the sign of the corresponding eigenvalue.
        
        ! Initialize values for function parser
        Dvalue = RCONST(0.0)
        Dvalue(NDIM3D+1) = dtime

        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)
        
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Get the normal vector in the point from the boundary
            dnx = merge(RCONST(1.0), -RCONST(1.0), mod(ibct,2) .eq. 0)
            
            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)
            
            ! Compute free stream values from function parser given in
            ! term of the primitive variables [rho,v1,p]
            do iexpr = 1, 3
              call fparser_evalFunction(p_rfparser,&
                  nmaxExpr*(isegment-1)+iexpr, Dvalue, DstateM(iexpr))
            end do
            
            ! Compute auxiliary quantities based on free stream state vector
            rM = DstateM(1)
            pM = DstateM(3)
            cM = sqrt(max((HYDRO_GAMMA)*pM/rM, SYS_EPSREAL_DP))
            dvnM = dnx*DstateM(2)
            
            ! Compute auxiliary quantities based on internal state vector
            pI = ((HYDRO_GAMMA)-RCONST(1.0))*(Daux1((ipoint-1)*NVAR1D+3,iel)-RCONST(0.5)*&
                (Daux1((ipoint-1)*NVAR1D+2,iel)**2))/Daux1((ipoint-1)*NVAR1D+1,iel)
            cI = sqrt(max((HYDRO_GAMMA)*pI/Daux1((ipoint-1)*NVAR1D+1,iel), SYS_EPSREAL_DP))

            ! Compute the normal velocity based on internal state vector
            dvnI = dnx*Daux1((ipoint-1)*NVAR1D+2,iel)/Daux1((ipoint-1)*NVAR1D+1,iel)

            ! Select free stream or computed Riemann invariant depending
            ! on the sign of the corresponding eigenvalue
            if (dvnI .lt. cI) then
              DstateM(1) = dvnM-(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*cM
            else
              DstateM(1) = dvnI-(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*cI
            end if

            if (dvnI .lt. SYS_EPSREAL_DP) then
              DstateM(2) = pM/(rM**(HYDRO_GAMMA))
            else
              DstateM(2) = pI/(Daux1((ipoint-1)*NVAR1D+1,iel)**(HYDRO_GAMMA))
            end if

            if (dvnI .lt. -cI) then
              DstateM(3) = dvnM+(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*cM
            else
              DstateM(3) = dvnI+(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*cI
            end if
            
            ! Convert Riemann invariants into conservative state variables
            cM = RCONST(0.25)*((HYDRO_GAMMA)-RCONST(1.0))*(DstateM(3)-DstateM(1))
            rM = (cM*cM/((HYDRO_GAMMA)*DstateM(2)))**(RCONST(1.0)/((HYDRO_GAMMA)-RCONST(1.0)))
            pM = rM*cM*cM/(HYDRO_GAMMA)
            dvnM = RCONST(0.5)*(DstateM(1)+DstateM(3))
            
            ! Setup the state vector based on Riemann invariants
            DstateM(1) = rM
            DstateM(2) = rM*dnx*dvnM
            DstateM(3) = pM/((HYDRO_GAMMA)-RCONST(1.0))+RCONST(0.5)*dvnM*dvnM
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux1((ipoint-1)*NVAR1D+1,iel)
            DstateI(2) = Daux1((ipoint-1)*NVAR1D+2,iel)
            DstateI(3) = Daux1((ipoint-1)*NVAR1D+3,iel)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*RCONST(0.5)*(Dflux-Diff)
          end do
        end do

        
      case (BDRC_FREESLIP)
        !-----------------------------------------------------------------------
        ! Free-slip boundary condition:
        !
        ! Compute the mirrored state vector based on the values of the
        ! computed state vector and use an approximate Riemann solver
        
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Get the normal vector in the point from the boundary
            dnx = merge(RCONST(1.0), -RCONST(1.0), mod(ibct,2) .eq. 0)
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux1((ipoint-1)*NVAR1D+1,iel)
            DstateI(2) = Daux1((ipoint-1)*NVAR1D+2,iel)
            DstateI(3) = Daux1((ipoint-1)*NVAR1D+3,iel)

            ! Compute the normal velocity based on the internal state vector
            dvnI = dnx*Daux1((ipoint-1)*NVAR1D+2,iel)/Daux1((ipoint-1)*NVAR1D+1,iel)

            ! Compute the mirrored state vector
            DstateM(1) = DstateI(1)
            DstateM(2) = DstateM(1)*(-dvnI*dnx)
            DstateM(3) = DstateI(3)

            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*RCONST(0.5)*(Dflux-Diff)
          end do
        end do

        
      case (BDRC_SUPERINLET)
        !-----------------------------------------------------------------------
        ! Supersonic inlet boundary conditions:
        !
        ! Prescribe the state vector in conservative variables
        
        ! Initialize values for function parser
        Dvalue = RCONST(0.0)
        Dvalue(NDIM3D+1) = dtime
        
        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)
        
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Get the normal vector in the point from the boundary
            dnx = merge(RCONST(1.0), -RCONST(1.0), mod(ibct,2) .eq. 0)
            
            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute boundary values from function parser given in
            ! term of the primitive variables [rho,v,p]
            do iexpr = 1, 3
              call fparser_evalFunction(p_rfparser,&
                  nmaxExpr*(isegment-1)+iexpr, Dvalue, DstateM(iexpr))
            end do

            ! Compute convervative variables
            DstateM(3) = DstateM(3)*(RCONST(1.0)/((HYDRO_GAMMA)-RCONST(1.0)))&
                       + DstateM(1)*RCONST(0.5)*(DstateM(2)**2)
            DstateM(2) = DstateM(1)*DstateM(2)

            ! Setup the computed internal state vector
            DstateI(1) = Daux1((ipoint-1)*NVAR1D+1,iel)
            DstateI(2) = Daux1((ipoint-1)*NVAR1D+2,iel)
            DstateI(3) = Daux1((ipoint-1)*NVAR1D+3,iel)

            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*RCONST(0.5)*(Dflux-Diff)
          end do
        end do


      case (BDRC_SUPEROUTLET)
        !-----------------------------------------------------------------------
        ! Supersonic outlet boundary conditions:
        !
        ! Evaluate the boundary fluxes based on the computed state vector

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement

            ! Get the normal vector in the point from the boundary
            dnx = merge(RCONST(1.0), -RCONST(1.0), mod(ibct,2) .eq. 0)
        
            ! Setup the computed internal state vector
            DstateI(1) = Daux1((ipoint-1)*NVAR1D+1,iel)
            DstateI(2) = Daux1((ipoint-1)*NVAR1D+2,iel)
            DstateI(3) = Daux1((ipoint-1)*NVAR1D+3,iel)

            ! Assemble Galerkin fluxes at the boundary
            call doGalerkinFlux(DstateI, dnx, Dflux)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*Dflux
          end do
        end do


      case (BDRC_SUBINLET)
        !-----------------------------------------------------------------------
        ! Subsonic pressure-density inlet boundary conditions:
        !
        ! Prescribe the density, pressure and tangential velocity at the inlet
        
        ! Initialize values for function parser
        Dvalue = RCONST(0.0)
        Dvalue(NDIM3D+1) = dtime

        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement

            ! Get the normal vector in the point from the boundary
            dnx = merge(RCONST(1.0), -RCONST(1.0), mod(ibct,2) .eq. 0)
            
            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute boundary values from function parser given in
            ! terms of the density and pressure
            do iexpr = 1, 2
              call fparser_evalFunction(p_rfparser,&
                  nmaxExpr*(isegment-1)+iexpr, Dvalue, DstateM(iexpr))
            end do

            ! Compute auxiliary quantities based on prescribed boundary values
            rM = DstateM(1)
            pM = DstateM(2)
            cM = sqrt(max((HYDRO_GAMMA)*pM/rM, SYS_EPSREAL_DP))

            ! Compute the normal velocity based on the internal state vector
            dvnI = dnx*Daux1((ipoint-1)*NVAR1D+2,iel)/Daux1((ipoint-1)*NVAR1D+1,iel)

            ! Compute the speed of sound based on the internal state vector
            cI = sqrt(max((HYDRO_GAMMA)*pI/Daux1((ipoint-1)*NVAR1D+1,iel), SYS_EPSREAL_DP))

            ! Compute fourth Riemann invariant based on the internal state vector
            w3 = dvnI+(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*cI

            ! Compute the first Riemann invariant based on the third Riemann
            ! invariant and the prescribed boundary values
            w1 = w3-2*(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*cM

            ! Setup the state vector based on Rimann invariants
            DstateM(1) = rM
            DstateM(2) = rM*dnx*RCONST(0.5)*(w1+w3)
            DstateM(3) = (RCONST(1.0)/&
                ((HYDRO_GAMMA)-RCONST(1.0)))*pM+RCONST(0.5)*rM*((dnx*RCONST(0.5)*(w1+w3))**2)

            ! Setup the computed internal state vector
            DstateI(1) = Daux1((ipoint-1)*NVAR1D+1,iel)
            DstateI(2) = Daux1((ipoint-1)*NVAR1D+2,iel)
            DstateI(3) = Daux1((ipoint-1)*NVAR1D+3,iel)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*RCONST(0.5)*(Dflux-Diff)
          end do
        end do

      case (BDRC_SUBOUTLET)
        !-----------------------------------------------------------------------
        ! Subsonic pressure outlet boundary condition:
        !
        ! Prescribe the pressure at the outlet

        ! Initialize values for function parser
        Dvalue = RCONST(0.0)
        Dvalue(NDIM3D+1) = dtime
        
        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement

            ! Get the normal vector in the point from the boundary
            dnx = merge(RCONST(1.0), -RCONST(1.0), mod(ibct,2) .eq. 0)
            
            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute pressure value from function parser
            call fparser_evalFunction(p_rfparser,&
                nmaxExpr*(isegment-1)+1, Dvalue, pM)

            ! Compute auxiliary quantities based on internal state vector
            pI = ((HYDRO_GAMMA)-RCONST(1.0))*(Daux1((ipoint-1)*NVAR1D+3,iel)-RCONST(0.5)*&
                (Daux1((ipoint-1)*NVAR1D+2,iel)**2)/Daux1((ipoint-1)*NVAR1D+1,iel))
            cI = sqrt(max((HYDRO_GAMMA)*pI/Daux1((ipoint-1)*NVAR1D+1,iel), SYS_EPSREAL_DP))

            ! Compute the normal velocity based on internal state vector
            dvnI = dnx*Daux1((ipoint-1)*NVAR1D+2,iel)/Daux1((ipoint-1)*NVAR1D+1,iel)
            
            ! Compute three Riemann invariants based on internal state vector
            DstateM(2) = pI/(Daux1((ipoint-1)*NVAR1D+1,iel)**(HYDRO_GAMMA))
            DstateM(3) = dvnI+(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*cI
            
            ! Compute first Riemann invariant based on third Riemann invariant,
            ! the computed density and pressure and the prescribed exit pressure
            DstateM(1) = DstateM(3)-2*(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*sqrt(max(SYS_EPSREAL_DP,&
                (HYDRO_GAMMA)*pM/Daux1((ipoint-1)*NVAR1D+1,iel)*(pI/pM)**(RCONST(1.0)/(HYDRO_GAMMA))))

            ! Convert Riemann invariants into conservative state variables
            cM = RCONST(0.25)*((HYDRO_GAMMA)-RCONST(1.0))*(DstateM(3)-DstateM(1))
            rM = (cM*cM/((HYDRO_GAMMA)*DstateM(2)))**(RCONST(1.0)/((HYDRO_GAMMA)-RCONST(1.0)))
            pM = rM*cM*cM/(HYDRO_GAMMA)
            dvnM = RCONST(0.5)*(DstateM(1)+DstateM(3))

            ! Setup the state vector based on Riemann invariants
            DstateM(1) = rM
            DstateM(2) = rM*dnx*dvnM
            DstateM(3) = pM/((HYDRO_GAMMA)-RCONST(1.0))+RCONST(0.5)*dvnM*dvnM
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux1((ipoint-1)*NVAR1D+1,iel)
            DstateI(2) = Daux1((ipoint-1)*NVAR1D+2,iel)
            DstateI(3) = Daux1((ipoint-1)*NVAR1D+3,iel)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*RCONST(0.5)*(Dflux-Diff)
          end do
        end do
                
      case default
        call output_line('Invalid type of boundary conditions!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_coeffVectorBdr1d_sim')
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
        call fevl_evaluate_sim(DER_FUNC1D, Daux2(:,:,ivar),&
            p_rsolution%RvectorBlock(ivar), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      end do
      
      ! What type of boundary conditions are we?
      select case(iand(ibdrtype, BDRC_TYPEMASK))
        
      case (BDRC_FREESTREAM)
        !-----------------------------------------------------------------------
        ! Free-stream boundary conditions:
        !
        ! Compute the Riemann invariants based on the computed (internal)
        ! state vector and the given freestream state vector and select
        ! the Riemman invariant for each characteristic fields based on
        ! the sign of the corresponding eigenvalue.
        
        ! Initialize values for function parser
        Dvalue = RCONST(0.0)
        Dvalue(NDIM3D+1) = dtime

        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)
        
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Get the normal vector in the point from the boundary
            dnx = merge(RCONST(1.0), -RCONST(1.0), mod(ibct,2) .eq. 0)
            
            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)
            
            ! Compute free stream values from function parser given in
            ! term of the primitive variables [rho,v1,p]
            do iexpr = 1, 3
              call fparser_evalFunction(p_rfparser,&
                  nmaxExpr*(isegment-1)+iexpr, Dvalue, DstateM(iexpr))
            end do
            
            ! Compute auxiliary quantities based on free stream state vector
            rM = DstateM(1)
            pM = DstateM(3)
            cM = sqrt(max((HYDRO_GAMMA)*pM/rM, SYS_EPSREAL_DP))
            dvnM = dnx*DstateM(2)
            
            ! Compute auxiliary quantities based on internal state vector
            pI = ((HYDRO_GAMMA)-RCONST(1.0))*(Daux2(ipoint,iel,3)-RCONST(0.5)*&
                (Daux2(ipoint,iel,2)**2))/Daux2(ipoint,iel,1)
            cI = sqrt(max((HYDRO_GAMMA)*pI/Daux2(ipoint,iel,1), SYS_EPSREAL_DP))

            ! Compute the normal velocity based on internal state vector
            dvnI = dnx*Daux2(ipoint,iel,2)/Daux2(ipoint,iel,1)

            ! Select free stream or computed Riemann invariant depending
            ! on the sign of the corresponding eigenvalue
            if (dvnI .lt. cI) then
              DstateM(1) = dvnM-(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*cM
            else
              DstateM(1) = dvnI-(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*cI
            end if

            if (dvnI .lt. SYS_EPSREAL_DP) then
              DstateM(2) = pM/(rM**(HYDRO_GAMMA))
            else
              DstateM(2) = pI/(Daux2(ipoint,iel,1)**(HYDRO_GAMMA))
            end if

            if (dvnI .lt. -cI) then
              DstateM(3) = dvnM+(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*cM
            else
              DstateM(3) = dvnI+(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*cI
            end if
            
            ! Convert Riemann invariants into conservative state variables
            cM = RCONST(0.25)*((HYDRO_GAMMA)-RCONST(1.0))*(DstateM(3)-DstateM(1))
            rM = (cM*cM/((HYDRO_GAMMA)*DstateM(2)))**(RCONST(1.0)/((HYDRO_GAMMA)-RCONST(1.0)))
            pM = rM*cM*cM/(HYDRO_GAMMA)
            dvnM = RCONST(0.5)*(DstateM(1)+DstateM(3))
            
            ! Setup the state vector based on Riemann invariants
            DstateM(1) = rM
            DstateM(2) = rM*dnx*dvnM
            DstateM(3) = pM/((HYDRO_GAMMA)-RCONST(1.0))+RCONST(0.5)*dvnM*dvnM
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux2(ipoint,iel,1)
            DstateI(2) = Daux2(ipoint,iel,2)
            DstateI(3) = Daux2(ipoint,iel,3)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*RCONST(0.5)*(Dflux-Diff)
          end do
        end do

        
      case (BDRC_FREESLIP)
        !-----------------------------------------------------------------------
        ! Free-slip boundary condition:
        !
        ! Compute the mirrored state vector based on the values of the
        ! computed state vector and use an approximate Riemann solver
        
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Get the normal vector in the point from the boundary
            dnx = merge(RCONST(1.0), -RCONST(1.0), mod(ibct,2) .eq. 0)
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux2(ipoint,iel,1)
            DstateI(2) = Daux2(ipoint,iel,2)
            DstateI(3) = Daux2(ipoint,iel,3)

            ! Compute the normal velocity based on the internal state vector
            dvnI = dnx*Daux2(ipoint,iel,2)/Daux2(ipoint,iel,1)

            ! Compute the mirrored state vector
            DstateM(1) = DstateI(1)
            DstateM(2) = DstateM(1)*(-dvnI*dnx)
            DstateM(3) = DstateI(3)

            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*RCONST(0.5)*(Dflux-Diff)
          end do
        end do

        
      case (BDRC_SUPERINLET)
        !-----------------------------------------------------------------------
        ! Supersonic inlet boundary conditions:
        !
        ! Prescribe the state vector in conservative variables
        
        ! Initialize values for function parser
        Dvalue = RCONST(0.0)
        Dvalue(NDIM3D+1) = dtime
        
        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)
        
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Get the normal vector in the point from the boundary
            dnx = merge(RCONST(1.0), -RCONST(1.0), mod(ibct,2) .eq. 0)
            
            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute boundary values from function parser given in
            ! term of the primitive variables [rho,v,p]
            do iexpr = 1, 3
              call fparser_evalFunction(p_rfparser,&
                  nmaxExpr*(isegment-1)+iexpr, Dvalue, DstateM(iexpr))
            end do

            ! Compute convervative variables
            DstateM(3) = DstateM(3)*(RCONST(1.0)/((HYDRO_GAMMA)-RCONST(1.0)))&
                       + DstateM(1)*RCONST(0.5)*(DstateM(2)**2)
            DstateM(2) = DstateM(1)*DstateM(2)

            ! Setup the computed internal state vector
            DstateI(1) = Daux2(ipoint,iel,1)
            DstateI(2) = Daux2(ipoint,iel,2)
            DstateI(3) = Daux2(ipoint,iel,33)

            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*RCONST(0.5)*(Dflux-Diff)
          end do
        end do


      case (BDRC_SUPEROUTLET)
        !-----------------------------------------------------------------------
        ! Supersonic outlet boundary conditions:
        !
        ! Evaluate the boundary fluxes based on the computed state vector

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement

            ! Get the normal vector in the point from the boundary
            dnx = merge(RCONST(1.0), -RCONST(1.0), mod(ibct,2) .eq. 0)
        
            ! Setup the computed internal state vector
            DstateI(1) = Daux2(ipoint,iel,1)
            DstateI(2) = Daux2(ipoint,iel,2)
            DstateI(3) = Daux2(ipoint,iel,3)

            ! Assemble Galerkin fluxes at the boundary
            call doGalerkinFlux(DstateI, dnx, Dflux)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*Dflux
          end do
        end do


      case (BDRC_SUBINLET)
        !-----------------------------------------------------------------------
        ! Subsonic pressure-density inlet boundary conditions:
        !
        ! Prescribe the density, pressure and tangential velocity at the inlet
        
        ! Initialize values for function parser
        Dvalue = RCONST(0.0)
        Dvalue(NDIM3D+1) = dtime

        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement

            ! Get the normal vector in the point from the boundary
            dnx = merge(RCONST(1.0), -RCONST(1.0), mod(ibct,2) .eq. 0)
            
            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute boundary values from function parser given in
            ! terms of the density and pressure
            do iexpr = 1, 2
              call fparser_evalFunction(p_rfparser,&
                  nmaxExpr*(isegment-1)+iexpr, Dvalue, DstateM(iexpr))
            end do

            ! Compute auxiliary quantities based on prescribed boundary values
            rM = DstateM(1)
            pM = DstateM(2)
            cM = sqrt(max((HYDRO_GAMMA)*pM/rM, SYS_EPSREAL_DP))

            ! Compute the normal velocity based on the internal state vector
            dvnI = dnx*Daux2(ipoint,iel,2)/Daux2(ipoint,iel,1)

            ! Compute the speed of sound based on the internal state vector
            cI = sqrt(max((HYDRO_GAMMA)*pI/Daux2(ipoint,iel,1), SYS_EPSREAL_DP))

            ! Compute fourth Riemann invariant based on the internal state vector
            w3 = dvnI+(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*cI

            ! Compute the first Riemann invariant based on the third Riemann
            ! invariant and the prescribed boundary values
            w1 = w3-2*(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*cM

            ! Setup the state vector based on Rimann invariants
            DstateM(1) = rM
            DstateM(2) = rM*dnx*RCONST(0.5)*(w1+w3)
            DstateM(3) = (RCONST(1.0)/((HYDRO_GAMMA)-RCONST(1.0)))*pM+RCONST(0.5)*rM*((dnx*RCONST(0.5)*(w1+w3))**2)

            ! Setup the computed internal state vector
            DstateI(1) = Daux2(ipoint,iel,1)
            DstateI(2) = Daux2(ipoint,iel,2)
            DstateI(3) = Daux2(ipoint,iel,3)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*RCONST(0.5)*(Dflux-Diff)
          end do
        end do

      case (BDRC_SUBOUTLET)
        !-----------------------------------------------------------------------
        ! Subsonic pressure outlet boundary condition:
        !
        ! Prescribe the pressure at the outlet

        ! Initialize values for function parser
        Dvalue = RCONST(0.0)
        Dvalue(NDIM3D+1) = dtime
        
        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement

            ! Get the normal vector in the point from the boundary
            dnx = merge(RCONST(1.0), -RCONST(1.0), mod(ibct,2) .eq. 0)
            
            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute pressure value from function parser
            call fparser_evalFunction(p_rfparser,&
                nmaxExpr*(isegment-1)+1, Dvalue, pM)

            ! Compute auxiliary quantities based on internal state vector
            pI = ((HYDRO_GAMMA)-RCONST(1.0))*(Daux2(ipoint,iel,3)-RCONST(0.5)*&
                (Daux2(ipoint,iel,2)**2)/Daux2(ipoint,iel,1))
            cI = sqrt(max((HYDRO_GAMMA)*pI/Daux2(ipoint,iel,1), SYS_EPSREAL_DP))

            ! Compute the normal velocity based on internal state vector
            dvnI = dnx*Daux2(ipoint,iel,2)/Daux2(ipoint,iel,11)
            
            ! Compute three Riemann invariants based on internal state vector
            DstateM(2) = pI/(Daux2(ipoint,iel,1)**(HYDRO_GAMMA))
            DstateM(3) = dvnI+(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*cI
            
            ! Compute first Riemann invariant based on third Riemann invariant,
            ! the computed density and pressure and the prescribed exit pressure
            DstateM(1) = DstateM(3)-2*(RCONST(2.0)/((HYDRO_GAMMA)-RCONST(1.0)))*sqrt(max(SYS_EPSREAL_DP,&
                (HYDRO_GAMMA)*pM/Daux2(ipoint,iel,1)*(pI/pM)**(RCONST(1.0)/(HYDRO_GAMMA))))

            ! Convert Riemann invariants into conservative state variables
            cM = RCONST(0.25)*((HYDRO_GAMMA)-RCONST(1.0))*(DstateM(3)-DstateM(1))
            rM = (cM*cM/((HYDRO_GAMMA)*DstateM(2)))**(RCONST(1.0)/((HYDRO_GAMMA)-RCONST(1.0)))
            pM = rM*cM*cM/(HYDRO_GAMMA)
            dvnM = RCONST(0.5)*(DstateM(1)+DstateM(3))

            ! Setup the state vector based on Riemann invariants
            DstateM(1) = rM
            DstateM(2) = rM*dnx*dvnM
            DstateM(3) = pM/((HYDRO_GAMMA)-RCONST(1.0))+RCONST(0.5)*dvnM*dvnM
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux2(ipoint,iel,1)
            DstateI(2) = Daux2(ipoint,iel,2)
            DstateI(3) = Daux2(ipoint,iel,3)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*RCONST(0.5)*(Dflux-Diff)
          end do
        end do
                
      case default
        call output_line('Invalid type of boundary conditions!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_coeffVectorBdr1d_sim')
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

      ! local variables
      real(DP) :: hI,hM,uI,pI,uM,pM,rI,rM
      real(DP) :: cPow2,c_IM,H_IM,q_IM,u_IM
      real(DP) :: l1,l2,l3,w1,w2,w3,aux,b1,b2,dveln
      
      ! Compute auxiliary quantities
      uI = XVELOCITY1(DstateI,IDX1)
      pI = PRESSURE1(DstateI,IDX1)
      rI = DENSITY1(DstateI,IDX1)
      hI = (TOTALENERGY1(DstateI,IDX1)+pI)/rI
      
      ! Compute auxiliary quantities
      uM = XVELOCITY1(DstateM,IDX1)
      pM = PRESSURE1(DstateM,IDX1)
      rM = DENSITY1(DstateM,IDX1)
      hM = (TOTALENERGY1(DstateM,IDX1)+pM)/rM

      ! Calculate $\frac12{\bf n}\cdot[{\bf F}(U_I)+{\bf F}(U_M)]$
      Dflux(1) = dnx*(DstateI(2) + DstateM(2))
      Dflux(2) = dnx*(DstateI(2)*uI+pI + DstateM(2)*uM+pM)
      Dflux(3) = dnx*((DstateI(3)+pI)*uI + (DstateM(3)+pM)*uM)

      
      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(rI,rM)
      u_IM = ROE_MEAN_VALUE(uI,uM,aux)
      H_IM = ROE_MEAN_VALUE(hI,hM,aux)
      
      ! Compute auxiliary variable
      q_IM  = RCONST(0.5)*u_IM**2

      ! Compute the speed of sound
      cPow2 = max(((HYDRO_GAMMA)-RCONST(1.0))*(H_IM-q_IM), SYS_EPSREAL_DP)
      c_IM  = sqrt(cPow2)

      ! Compute normal velocity
      dveln = dnx*uI

      ! Compute eigenvalues
      l1 = abs(dveln-c_IM)
      l2 = abs(dveln)
      l3 = abs(dveln+c_IM)
      
      ! Compute solution difference U_M-U_I
      Diff = DstateM-DstateI
      
      ! Compute auxiliary quantities for characteristic variables
      b2 = ((HYDRO_GAMMA)-RCONST(1.0))/cPow2; b1 = b2*q_IM
      
      ! Compute characteristic variables multiplied by the
      ! corresponding eigenvalue
      w1 = l1 * RCONST(0.5) * (   (b1+dveln/c_IM)*Diff(1) -&
                          (b2*u_IM+dnx/c_IM)*Diff(2) +&
                                          b2*Diff(3) )
      w2 = l2 * (                (RCONST(1.0)-b1)*Diff(1) +&
                                     b2*u_IM*Diff(2) -&
                                          b2*Diff(3) )
      w3 = l3 * RCONST(0.5) * (   (b1-dveln/c_IM)*Diff(1) -&
                          (b2*u_IM-dnx/c_IM)*Diff(2) +&
                                          b2*Diff(3) )
            
      ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
      Diff(1) = w1 + w2 + w3
      Diff(2) = (u_IM-c_IM*dnx)*w1 + u_IM*w2 +&
                (u_IM+c_IM*dnx)*w3
      Diff(3) = (H_IM-c_IM*dveln)*w1 + q_IM*w2 +&
                (H_IM+c_IM*dveln)*w3

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
      real(DP) :: u,p
      
      
      ! Compute auxiliary quantities
      u = XVELOCITY1(Dstate,IDX1)
      p = PRESSURE1(Dstate,IDX1)
      
      ! Calculate ${\bf n}\cdot{\bf F}(U)$
      Dflux(1) = dnx*Dstate(2)
      Dflux(2) = dnx*(Dstate(2)*u+p)
      Dflux(3) = dnx*(Dstate(3)+p)*u

    end subroutine doGalerkinFlux

  end subroutine hydro_coeffVectorBdr1d_sim

  !*****************************************************************************

!<subroutine>

  subroutine hydro_hadaptCallbackScalar1d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 1D. The solution vector is assumed
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
    !   IquickAccess(2:3):   ivt1,ivt2
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
            OU_CLASS_WARNING,OU_MODE_STD,'hydro_hadaptCallbackScalar1d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR1D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR1D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR1D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR1D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      do ivar = 1, NVAR1D
        p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR1D+ivar) = &
            RCONST(0.5)*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR1D+ivar)+&
                    p_Dsolution((rcollection%IquickAccess(3)-1)*NVAR1D+ivar))
      end do

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        do ivar = 1, NVAR1D
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR1D+ivar) = &
              p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR1D+ivar)
        end do
      else
        do ivar = 1, NVAR1D
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR1D+ivar) = RCONST(0.0)
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)

    end select

  end subroutine hydro_hadaptCallbackScalar1d

  !*****************************************************************************

!<subroutine>

  subroutine hydro_hadaptCallbackBlock1d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 1D. The solution vector is assumed
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
    !   IquickAccess(2:3):   ivt1,ivt2
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
      if (rsolution%nblocks .ne. NVAR1D) then
        call output_line('Vector is not in block format!',&
            OU_CLASS_WARNING,OU_MODE_STD,'hydro_hadaptCallbackBlock1d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR1D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR1D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR1D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR1D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      neq = rsolution%NEQ/NVAR1D
      do ivar = 1, NVAR1D
        p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = &
            RCONST(0.5)*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
                    p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(3)) )
      end do

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        neq = rsolution%NEQ/NVAR1D
        do ivar = 1, NVAR1D
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = &
              p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))
        end do
      else
        neq = rsolution%NEQ/NVAR1D
        do ivar = 1, NVAR1D
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = RCONST(0.0)
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)

    end select

  end subroutine hydro_hadaptCallbackBlock1d

end module hydro_callback1d
