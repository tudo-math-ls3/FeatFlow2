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
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#
!# 2.) hydro_calcFluxGalNoBdr1d_sim
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#        without assembling the symmetric boundary contribution
!#
!# 3.) hydro_calcFluxScDiss1d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting scalar artificial viscosities
!#
!# 4.) hydro_calcFluxRoeDiss1d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities
!#
!# 5.) hydro_calcFluxRusDiss1d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion
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
!#         adopting tensorial artificial viscosities
!#
!# 13.) hydro_calcMatRoeDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 14.) hydro_calcMatRusDissMatD1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov artificial viscosities
!#
!# 15.) hydro_calcMatRusDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov flux artificial viscosities
!#
!# 16.) hydro_calcCharacteristics1d_sim
!#      -> Computes characteristic variables
!#
!# 17.) hydro_calcFluxFCTScalarDiss1d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting scalar artificial viscosities
!#
!# 18.) hydro_calcFluxFCTRoeDiss1d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting tensorial artificial viscosities
!#
!# 19.) hydro_calcFluxFCTRusanov1d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting the Rusanov artificial viscosities
!#
!# 20.) hydro_trafoFluxDensity1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 21.) hydro_trafoDiffDensity1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 22.) hydro_trafoFluxEnergy1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 23.) hydro_trafoDiffEnergy1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 24.) hydro_trafoFluxPressure1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 25.) hydro_trafoDiffPressure1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the pressure
!#
!# 26.) hydro_trafoFluxVelocity1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 27.) hydro_trafoDiffVelocity1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 28.) hydro_trafoFluxMomentum1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 29.) hydro_trafoDiffMomentum1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 30.) hydro_trafoFluxDenEng1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 31.) hydro_trafoDiffDenEng1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 32.) hydro_trafoFluxDenPre1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 33.) hydro_trafoDiffDenPre1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 34.) hydro_trafoFluxDenPreVel1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 35.) hydro_trafoDiffDenPreVel1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure 
!#         and the velocity
!#
!# 36.) hydro_calcBoundaryvalues1d
!#      -> Computes the boundary values for a given node
!#
!# 37.) hydro_hadaptCallbackScalar1d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 1D, whereby the vector is stored in interleave format
!#
!# 38.) hydro_hadaptCallbackBlock1d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 1D, whereby the vector is stored in block format
!#
!# </purpose>
!##############################################################################

module hydro_callback1d

#include "hydro.h"

  use boundarycondaux
  use collection
  use flagship_callback
  use fsystem
  use genoutput
  use graph
  use groupfemsystem
  use hadaptaux
  use hydro_basic
  use linearsystemblock
  use linearsystemscalar
  use problem
  use solveraux
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
  public :: hydro_calcFluxFCTScalarDiss1d
  public :: hydro_calcFluxFCTRoeDiss1d
  public :: hydro_calcFluxFCTRusanov1d
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
  public :: hydro_calcBoundaryvalues1d
  public :: hydro_hadaptCallbackScalar1d
  public :: hydro_hadaptCallbackBlock1d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGal1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the standard
    ! Galerkin discretisation in 1D.
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
#ifdef HYDRO_USE_IBP
    real(DP), dimension(NVAR1D) :: Fi,Fj
#else
    real(DP), dimension(NVAR1D) :: F_ij
#endif
    real(DP) :: pi,pj,ui,uj
    integer :: idx


    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !
      !     / rho*u         \
      ! F = | rho*u*u + p   |
      !     \ rho*E*u + p*u /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_HYDRO_2T_XDIR_1D(Fi,DdataAtEdge,1,idx,ui,pi)
      FLUX_HYDRO_2T_XDIR_1D(Fj,DdataAtEdge,2,idx,uj,pj)
      
      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fi )
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_HYDRO_2T_XDIR_1D(F_ij,DdataAtEdge,1,2,idx,ui,uj,pi,pj)
                 
      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * DmatrixCoeffsAtEdge(1,1,idx)*F_ij
      DfluxesAtEdge(:,2,idx) = -dscale * DmatrixCoeffsAtEdge(1,2,idx)*F_ij      
#endif
    end do

  end subroutine hydro_calcFluxGal1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGalNoBdr1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the TVD
    ! discretisation in 1D. The symmetric boundary contributions
    ! are neglected and incorporated in the antidiffusive flux.
    ! Hence, this is simply the standard Galerkin flux for the
    ! skew-symmetric internal contributions.
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
    real(DP), dimension(NVAR1D) :: F_ij
    real(DP) :: pi,pj,ui,uj
    integer :: idx


    do idx = 1, size(DfluxesAtEdge,3)

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !
      !     / rho*u         \
      ! F = | rho*u*u + p   |
      !     \ rho*E*u + p*u /
      !-------------------------------------------------------------------------
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute flux difference for x-direction
      FLUXDIFF_HYDRO_2T_XDIR_1D(F_ij,DdataAtEdge,1,2,idx,ui,uj,pi,pj)

      ! Assemble symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale *&
          0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))*F_ij
      DfluxesAtEdge(:,2,idx) = DfluxesAtEdge(:,1,idx)
    end do

  end subroutine hydro_calcFluxGalNoBdr1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxScDiss1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 1D using scalar dissipation.
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
#ifdef HYDRO_USE_IBP
    real(DP), dimension(NVAR1D) :: Fi,Fj
#else
    real(DP), dimension(NVAR1D) :: F_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: pi,pj,ui,uj
    real(DP) :: aux,anorm,vel_ij,c_ij,d_ij,H_ij,q_ij,u_ij
    integer :: idx

    do idx = 1, size(DfluxesAtEdge,3)
          
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !
      !     / rho*u         \
      ! F = | rho*u*u + p   |
      !     \ rho*E*u + p*u /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_HYDRO_2T_XDIR_1D(Fi,DdataAtEdge,1,idx,ui,pi)
      FLUX_HYDRO_2T_XDIR_1D(Fj,DdataAtEdge,2,idx,uj,pj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_HYDRO_2T_XDIR_1D(F_ij,DdataAtEdge,1,2,idx,ui,uj,pi,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral radius
      !-------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))
      anorm = abs(a(1)) ! = sqrt(a(1)*a(1))

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      H_ij = ROE_MEAN_VALUE(\
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+pi)/\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+pj)/\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)
      
      ! Compute auxiliary variables
      vel_ij = u_ij*a(1)
      q_ij   = 0.5*(u_ij*u_ij)

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Scalar dissipation
      d_ij = abs(vel_ij) + anorm*c_ij
      
      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------
      
#ifdef HYDRO_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*F_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*F_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxScDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRoeDiss1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 1D using tensorial dissipation.
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
#ifdef HYDRO_USE_IBP
    real(DP), dimension(NVAR1D) :: Fi,Fj
#else
    real(DP), dimension(NVAR1D) :: F_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: pi,pj,ui,uj,b1,b2,aux,H_ij,q_ij,u_ij
    real(DP) :: anorm,l1,l2,l3,w1,w2,w3,c_ij,cPow2_ij
    integer :: idx
    
    
    do idx = 1, size(DfluxesAtEdge,3)
          
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !
      !     / rho*u         \
      ! F = | rho*u*u + p   |
      !     \ rho*E*u + p*u /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_HYDRO_2T_XDIR_1D(Fi,DdataAtEdge,1,idx,ui,pi)
      FLUX_HYDRO_2T_XDIR_1D(Fj,DdataAtEdge,2,idx,uj,pj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_HYDRO_2T_XDIR_1D(F_ij,DdataAtEdge,1,2,idx,ui,uj,pi,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor by Roe
      !-------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))
      anorm = abs(a(1)) ! = sqrt(a(1)*a(1))

      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      H_ij = ROE_MEAN_VALUE(\
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+pi)/\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+pj)/\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)
        
        ! Compute auxiliary variables
        q_ij = 0.5*(u_ij*u_ij)

        ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
        cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(cPow2_ij)
        
        ! Compute eigenvalues
        l1 = abs(u_ij-c_ij)
        l2 = abs(u_ij)
        l3 = abs(u_ij+c_ij)
        
        ! Compute solution difference U_j-U_i
        Diff = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
        
        ! Compute auxiliary quantities for characteristic variables
        b2 = (GAMMA-1.0)/cPow2_ij; b1 = b2*q_ij
        
        ! Compute characteristic variables multiplied by the
        ! corresponding eigenvalue
        w1 = l1 * 0.5 * (    (b1+u_ij/c_ij)*Diff(1) -&
                         (b2*u_ij+1.0/c_ij)*Diff(2) +&
                                         b2*Diff(3) )
        w2 = l2 * (                  (1-b1)*Diff(1) +&
                                    b2*u_ij*Diff(2) -&
                                         b2*Diff(3) )
        w3 = l3 * 0.5 * (    (b1-u_ij/c_ij)*Diff(1) -&
                         (b2*u_ij-1.0/c_ij)*Diff(2) +&
                                         b2*Diff(3) )
        
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        Diff(1) = anorm * ( w1 + w2 + w3 )
        Diff(2) = anorm * ( (u_ij-c_ij)*w1 + u_ij*w2 + (u_ij+c_ij)*w3 )
        Diff(3) = anorm * ( (H_ij-u_ij*c_ij)*w1 + q_ij*w2 + (H_ij+u_ij*c_ij)*w3 )
        
        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fi + Diff)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else        
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*F_ij + Diff)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*F_ij + Diff)
#endif
      else

#ifdef HYDRO_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fi )
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*F_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*F_ij)
#endif
      end if
    end do

  end subroutine hydro_calcFluxRoeDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRusDiss1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)


!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 1D using the Rusanov dissipation.
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
#ifdef HYDRO_USE_IBP
    real(DP), dimension(NVAR1D) :: Fi,Fj
#else
    real(DP), dimension(NVAR1D) :: F_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP) :: pi,pj,ui,uj,d_ij,ci,cj,Ei,Ej
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !
      !     / rho*u         \
      ! F = | rho*u*u + p   |
      !     \ rho*E*u + p*u /
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_HYDRO_2T_XDIR_1D(Fi,DdataAtEdge,1,idx,ui,pi)
      FLUX_HYDRO_2T_XDIR_1D(Fj,DdataAtEdge,2,idx,uj,pj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_HYDRO_2T_XDIR_1D(F_ij,DdataAtEdge,1,2,idx,ui,uj,pi,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !-------------------------------------------------------------------------
      
      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*ui*ui), SYS_EPSREAL))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*uj*uj), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Scalar dissipation for the Rusanov flux
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj)+&
                  abs(DmatrixCoeffsAtEdge(1,1,idx))*cj,&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui)+&
                  abs(DmatrixCoeffsAtEdge(1,2,idx))*ci )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*F_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*F_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxRusDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiagMatD1d_sim(DdataAtNode, DmatrixCoeffsAtNode,&
      IverticesAtNode, dscale, DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 1D
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

    ! local variable
    real(DP) :: ui
    integer :: inode

    do inode = 1, size(DcoefficientsAtNode,3)
      
      ! Compute velocity
      ui = X_VELOCITY_1T_FROM_CONSVAR(DdataAtNode,NVAR1D,inode)
      
#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ii = diag(A_i)*C_{ii}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtNode,1,inode,\
        dscale,DmatrixCoeffsAtNode(1,inode),ui)
#else
      ! Compute Galerkin coefficient $K_ii = -diag(A_i)*C_{ii}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtNode,1,inode,\
        -dscale,DmatrixCoeffsAtNode(1,inode),ui)
#endif
    end do

  end subroutine hydro_calcMatDiagMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiag1d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 1D
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

    ! local variable
    real(DP) :: ui,Ei
    integer :: inode

    do inode = 1, size(DcoefficientsAtNode,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_1T_FROM_CONSVAR(DdataAtNode,NVAR1D,inode)
      Ei = TOTAL_ENERGY_1T_FROM_CONSVAR(DdataAtNode,NVAR1D,inode)

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtNode,1,inode,\
        dscale,DmatrixCoeffsAtNode(1,inode),ui,Ei)
#else
      ! Compute Galerkin coefficient $K_ii = -A_i*C_{ii}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtNode,1,inode,\
        -dscale,DmatrixCoeffsAtNode(1,inode),ui,Ei)
#endif
    end do
    
  end subroutine hydro_calcMatDiag1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGalMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 1D
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

    ! local variable
    real(DP) :: ui,uj
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),uj)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),ui)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),uj)

      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),ui)
#endif
    end do

  end subroutine hydro_calcMatGalMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGal1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 1D
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

    ! local variable
    real(DP) :: ui,uj,Ei,Ej
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      Ei = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      Ej = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),ui,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),ui,Ei)
#endif
    end do

  end subroutine hydro_calcMatGal1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDissMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies scalar artificial viscosities in 1D
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

    ! local variable
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: anorm,aux,H_ij,q_ij,ui,uj,u_ij
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),uj)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),ui)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),uj)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),ui)
#endif
      
      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ji}-C_{ij})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))
      anorm = abs(a(1))
      
      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)
        
        ! Compute auxiliary values
        q_ij = 0.5*(u_ij*u_ij)
        
        ! Compute scalar dissipation
#ifdef THERMALLY_IDEAL_GAS
        DcoefficientsAtEdge(:,1,idx) = dscale * (abs(a(1)*u_ij) +&
            anorm*sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)))
#else
#error "Speed of sound must be implemented!"
#endif
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0
        
      end if
    end do

  end subroutine hydro_calcMatScDissMatD1d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDiss1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies scalar artificial viscosities in 1D
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

    ! local variable
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: anorm,aux,Ei,Ej,H_ij,q_ij,ui,uj,u_ij
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
    
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      Ei = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      Ej = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),ui,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),ui,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))
      anorm = abs(a(1))
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)
        
        ! Compute auxiliary values
        q_ij = 0.5*(u_ij*u_ij)
        
        ! Compute scalar dissipation
#ifdef THERMALLY_IDEAL_GAS
        aux = dscale * (abs(a(1)*u_ij) +&
            anorm*sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)))
#else
#error "Speed of sound must be implemented!"
#endif
        
        DcoefficientsAtEdge(1,1,idx) = aux
        DcoefficientsAtEdge(5,1,idx) = aux
        DcoefficientsAtEdge(9,1,idx) = aux
      end if
    end do
    
  end subroutine hydro_calcMatScDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDissMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 1D
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

    ! local variable
    real(DP), dimension(NVAR1D,NVAR1D) :: R_ij,L_ij
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: aux,H_ij,q_ij,ui,uj,u_ij
    real(DP) :: l1,l2,l3,anorm,c_ij,cPow2_ij,b1,b2
    integer :: idx
    
    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),uj)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),ui)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),uj)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),ui)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(1,1,idx)-&
                  DmatrixCoeffsAtEdge(1,2,idx))
      anorm = abs(a(1))
      
      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)
                
        ! Compute auxiliary values
        q_ij = 0.5*(u_ij*u_ij)
#ifdef THERMALLY_IDEAL_GAS
        cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(cPow2_ij)
        
        b2   = (GAMMA-1.0)/cPow2_ij
        b1   = b2*q_ij

        ! Diagonal matrix of eigenvalues
        l1 = abs(u_ij-c_ij)
        l2 = abs(u_ij)
        l3 = abs(u_ij+c_ij)
        
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
        L_ij(1,1) = 0.5*(b1+u_ij/c_ij)
        L_ij(2,1) = 1.0-b1
        L_ij(3,1) = 0.5*(b1-u_ij/c_ij)
        
        L_ij(1,2) = -0.5*(b2*u_ij+1/c_ij)
        L_ij(2,2) =  b2*u_ij
        L_ij(3,2) = -0.5*(b2*u_ij-1/c_ij)
        
        L_ij(1,3) = 0.5*b2
        L_ij(2,3) = -b2
        L_ij(3,3) = 0.5*b2

        ! Include scaling parameter
        anorm = dscale*anorm

        ! Compute tensorial dissipation D_ij = diag(R_ij*|Lbd_ij|*L_ij)*I
        DcoefficientsAtEdge(1,1,idx) = anorm*( R_ij(1,1)*L_ij(1,1)+&
                                               R_ij(1,2)*L_ij(2,1)+&
                                               R_ij(1,3)*L_ij(3,1)  )
        DcoefficientsAtEdge(2,1,idx) = anorm*( R_ij(2,1)*L_ij(1,2)+&
                                               R_ij(2,2)*L_ij(2,2)+&
                                               R_ij(2,3)*L_ij(3,2)  )
        DcoefficientsAtEdge(3,1,idx) = anorm*( R_ij(3,1)*L_ij(1,3)+&
                                               R_ij(3,2)*L_ij(2,3)+&
                                               R_ij(3,3)*L_ij(3,3)  )
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0

      end if
    end do

  end subroutine hydro_calcMatRoeDissMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDiss1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 1D
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

    ! local variable
    real(DP), dimension(NVAR1D,NVAR1D) :: R_ij,L_ij
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: aux,Ei,Ej,H_ij,q_ij,ui,uj,u_ij
    real(DP) :: l1,l2,l3,anorm,c_ij,cPow2_ij,b1,b2
    integer :: idx,i,j,k

    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      Ei = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      Ej = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),ui,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),ui,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))
      anorm = abs(a(1))

      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)
        
        ! Compute auxiliary values
        q_ij = 0.5*u_ij*u_ij
#ifdef THERMALLY_IDEAL_GAS
        cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(cPow2_ij)

        b2   = (GAMMA-1.0)/cPow2_ij
        b1   = b2*q_ij
        
        ! Diagonal matrix of eigenvalues
        l1 = abs(u_ij-c_ij)
        l2 = abs(u_ij)
        l3 = abs(u_ij+c_ij)
        
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
        L_ij(1,1) = 0.5*(b1+u_ij/c_ij)
        L_ij(2,1) = 1.0-b1
        L_ij(3,1) = 0.5*(b1-u_ij/c_ij)
        
        L_ij(1,2) = -0.5*(b2*u_ij+1/c_ij)
        L_ij(2,2) =  b2*u_ij
        L_ij(3,2) = -0.5*(b2*u_ij-1/c_ij)
        
        L_ij(1,3) = 0.5*b2
        L_ij(2,3) = -b2
        L_ij(3,3) = 0.5*b2
        
        ! Include scaling parameter
        anorm = dscale*anorm

        ! Compute tensorial dissipation D_ij = R_ij*|Lbd_ij|*L_ij
        do i = 1, NVAR1D
          do j = 1, NVAR1D
            aux = 0.0
            do k = 1, NVAR1D
              aux = aux + R_ij(i,k)*L_ij(k,j)
            end do
            DcoefficientsAtEdge(NVAR1D*(j-1)+i,1,idx) = anorm*aux
          end do
        end do
        
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0
        
      end if
    end do

  end subroutine hydro_calcMatRoeDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDissMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 1D
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

    ! local variable
    real(DP) :: ui,uj,ci,cj,Ei,Ej
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      Ei = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      Ej = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),uj)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),ui)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),uj)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),ui)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*ui*ui), SYS_EPSREAL))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*uj*uj), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Compute dissipation tensor D_ij
      DcoefficientsAtEdge(:,1,idx) = dscale *&
          max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj)+&
               abs(DmatrixCoeffsAtEdge(1,1,idx))*cj,&
               abs(DmatrixCoeffsAtEdge(1,2,idx)*ui)+&
               abs(DmatrixCoeffsAtEdge(1,2,idx))*ci )
    end do
    
  end subroutine hydro_calcMatRusDissMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDiss1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 1D
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

    ! local variable
    real(DP) :: ui,uj,ci,cj,Ei,Ej,aux
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      Ei = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      Ej = TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),ui,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),uj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
       MATRIX_HYDRO_2T_1D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),ui,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*ui*ui), SYS_EPSREAL))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*uj*uj), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif

      ! Compute dissipation tensor D_ij
      aux = dscale * max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj)+&
                          abs(DmatrixCoeffsAtEdge(1,1,idx))*cj,&
                          abs(DmatrixCoeffsAtEdge(1,2,idx)*ui)+&
                          abs(DmatrixCoeffsAtEdge(1,2,idx))*ci )

      DcoefficientsAtEdge(:,1,idx) = 0.0
      DcoefficientsAtEdge(1,1,idx) = aux
      DcoefficientsAtEdge(5,1,idx) = aux
      DcoefficientsAtEdge(9,1,idx) = aux
    end do
    
  end subroutine hydro_calcMatRusDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcCharacteristics1d_sim(Dweight, DdataAtEdge,&
      DcharVariablesAtEdge, DeigenvaluesAtEdge,&
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

    ! local variables
    real(DP), dimension(NVAR1D) :: Diff
    real(DP) :: u_ij,H_ij,q_ij,c_ij,aux,ui,uj,cPow2_ij,anorm,b1_ij,b2_ij
    integer :: idx

    ! Compute norm of weighting coefficient
    anorm = abs(Dweight(1))

    ! Check if weighting coefficient is zero
    if (anorm .le. SYS_EPSREAL) then
      if (present(DcharVariablesAtEdge))     DcharVariablesAtEdge     = 0.0
      if (present(DeigenvaluesAtEdge))       DeigenvaluesAtEdge       = 0.0
      if (present(DrightEigenvectorsAtEdge)) DrightEigenvectorsAtEdge = 0.0
      if (present(DleftEigenvectorsAtEdge))  DleftEigenvectorsAtEdge  = 0.0

      ! That's it
      return
    end if
    
    
    ! Do we have to compute characteristic variables
    if (present(DcharVariablesAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)
        
        ! Compute velocities
        ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
        uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)
               
        ! Compute auxiliary variables
        q_ij = 0.5*(u_ij*u_ij)
#ifdef THERMALLY_IDEAL_GAS
        cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(cPow2_ij)

        b2_ij = (GAMMA-1.0)/cPow2_ij
        b1_ij  = b2_ij*q_ij
        
        ! Compute solution difference U_j-U_i
        Diff = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
          
        ! Compute characteristic variables
        DcharVariablesAtEdge(1,idx) = anorm * 0.5 *&
            (       (b1_ij+u_ij/c_ij)*Diff(1)-&
                (b2_ij*u_ij+1.0/c_ij)*Diff(2)+&
                                b2_ij*Diff(3))
        DcharVariablesAtEdge(2,idx) = anorm *&
            (             (1.0-b1_ij)*Diff(1)+&
                           b2_ij*u_ij*Diff(2)-&
                                b2_ij*Diff(3) )
        DcharVariablesAtEdge(3,idx) = anorm * 0.5 *&
            (       (b1_ij-u_ij/c_ij)*Diff(1)-&
                (b2_ij*u_ij-1.0/c_ij)*Diff(2)+&
                                b2_ij*Diff(3) )
      end do
    end if


    ! Do we have to compute eigenvalues
    if (present(DeigenvaluesAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)
        
        ! Compute velocities
        ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
        uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)
        
        ! Compute auxiliary variable
#ifdef THERMALLY_IDEAL_GAS
        c_ij = sqrt(max((GAMMA-1.0)*(H_ij-0.5*(u_ij*u_ij)), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif

        ! Compute eigenvalues
        DeigenvaluesAtEdge(1,idx) = u_ij-c_ij
        DeigenvaluesAtEdge(2,idx) = u_ij
        DeigenvaluesAtEdge(3,idx) = u_ij+c_ij
      end do
    end if


    ! Do we have to compute right eigenvectors
    if (present(DrightEigenvectorsAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)
        
        ! Compute velocities
        ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
        uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)

        ! Compute auxiliary variables
        q_ij = 0.5*(u_ij*u_ij)
#ifdef THERMALLY_IDEAL_GAS
        c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif

        ! Compute right eigenvectors
        DrightEigenvectorsAtEdge(1,idx) =  1.0
        DrightEigenvectorsAtEdge(2,idx) =  u_ij-c_ij
        DrightEigenvectorsAtEdge(3,idx) =  H_ij-u_ij*c_ij

        DrightEigenvectorsAtEdge(4,idx) =  1.0
        DrightEigenvectorsAtEdge(5,idx) =  u_ij
        DrightEigenvectorsAtEdge(6,idx) =  q_ij

        DrightEigenvectorsAtEdge(7,idx) =  1.0
        DrightEigenvectorsAtEdge(8,idx) =  u_ij+c_ij
        DrightEigenvectorsAtEdge(9,idx) =  H_ij+u_ij*c_ij
      end do
    end if


    ! Do we have to compute left eigenvectors
    if (present(DleftEigenvectorsAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)
        
        ! Compute velocities
        ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
        uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+\
               PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)

        ! Compute auxiliary variables
        q_ij = 0.5*(u_ij*u_ij)
#ifdef THERMALLY_IDEAL_GAS
        cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(cPow2_ij)

        b2_ij = (GAMMA-1.0)/cPow2_ij
        b1_ij = b2_ij*q_ij

        ! Compute left eigenvectors
        DleftEigenvectorsAtEdge(1,idx) = 0.5 * (b1_ij+u_ij/c_ij)
        DleftEigenvectorsAtEdge(2,idx) =       (1.0-b1_ij)
        DleftEigenvectorsAtEdge(3,idx) = 0.5 * (b1_ij-u_ij/c_ij)

        DleftEigenvectorsAtEdge(4,idx) =-0.5 * (b2_ij*u_ij+1.0/c_ij)
        DleftEigenvectorsAtEdge(5,idx) =       (b2_ij*u_ij)
        DleftEigenvectorsAtEdge(6,idx) =-0.5 * (b2_ij*u_ij-1.0/c_ij)

        DleftEigenvectorsAtEdge(7,idx) = 0.5*b2_ij
        DleftEigenvectorsAtEdge(8,idx) =    -b2_ij
        DleftEigenvectorsAtEdge(9,idx) = 0.5*b2_ij
      end do
    end if
    
  end subroutine hydro_calcCharacteristics1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTScalarDiss1d(&
      U1_i, U1_j, U2_i, U2_j, Coeff_ij, Coeff_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 1D using scalar dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: Coeff_ij,Coeff_ji

    ! scaling coefficients
    real(DP), intent(in) :: dscale1,dscale2

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! raw antidiffusive flux
    real(DP), dimension(:), intent(out) :: F_ij
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: ui,uj,d_ij,H_ij,q_ij,u_ij,aux,vel_ij,c_ij,anorm

    ! Compute velocities
    ui = X_VELOCITY_FROM_CONSVAR(U2_i,NVAR1D)
    uj = X_VELOCITY_FROM_CONSVAR(U2_j,NVAR1D)

    ! Compute skew-symmetric coefficient
    a = 0.5*(Coeff_ij-Coeff_ji); anorm = abs(a(1))

    ! Compute Roe mean values
    aux  = ROE_MEAN_RATIO(\
           DENSITY_FROM_CONSVAR(U2_i,NVAR1D),\
           DENSITY_FROM_CONSVAR(U2_j,NVAR1D))
    u_ij = ROE_MEAN_VALUE(ui,uj,aux)
    H_ij = ROE_MEAN_VALUE(\
           (TOTAL_ENERGY_FROM_CONSVAR(U2_i,NVAR1D)+\
            PRESSURE_FROM_CONSVAR_1D(U2_i,NVAR1D))/DENSITY_FROM_CONSVAR(U2_i,NVAR1D),\
           (TOTAL_ENERGY_FROM_CONSVAR(U2_j,NVAR1D)+\
            PRESSURE_FROM_CONSVAR_1D(U2_j,NVAR1D))/DENSITY_FROM_CONSVAR(U2_j,NVAR1D),aux)

    ! Compute auxiliary variables
    vel_ij = u_ij*a(1)
    q_ij   = 0.5*(u_ij*u_ij)
#ifdef THERMALLY_IDEAL_GAS
    c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented!"
#endif

    ! Scalar dissipation
    d_ij = abs(vel_ij) + anorm*c_ij

    ! Compute conservative fluxes
    F_ij = dscale1*(U1_i-U1_j) + dscale2*d_ij*(U2_i-U2_j)

  end subroutine hydro_calcFluxFCTScalarDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRoeDiss1d(&
      U1_i, U1_j, U2_i, U2_j, Coeff_ij, Coeff_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 1D using tensorial dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: Coeff_ij,Coeff_ji

    ! scaling coefficients
    real(DP), intent(in) :: dscale1,dscale2

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! raw antidiffusive flux
    real(DP), dimension(:), intent(out) :: F_ij
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: ui,uj,b1,b2,aux,H_ij,q_ij,u_ij
    real(DP) :: anorm,l1,l2,l3,w1,w2,w3,cPow2_ij,c_ij

    ! Compute velocities
    ui = X_VELOCITY_FROM_CONSVAR(U2_i,NVAR1D)
    uj = X_VELOCITY_FROM_CONSVAR(U2_j,NVAR1D)

    ! Compute skew-symmetric coefficient
    a = 0.5*(Coeff_ij-Coeff_ji); anorm = abs(a(1))

    ! Compute Roe mean values
    aux  = ROE_MEAN_RATIO(\
           DENSITY_FROM_CONSVAR(U2_i,NVAR1D),\
           DENSITY_FROM_CONSVAR(U2_j,NVAR1D))
    u_ij = ROE_MEAN_VALUE(ui,uj,aux)
    H_ij = ROE_MEAN_VALUE(\
           (TOTAL_ENERGY_FROM_CONSVAR(U2_i,NVAR1D)+\
            PRESSURE_FROM_CONSVAR_1D(U2_i,NVAR1D))/DENSITY_FROM_CONSVAR(U2_i,NVAR1D),\
           (TOTAL_ENERGY_FROM_CONSVAR(U2_j,NVAR1D)+\
            PRESSURE_FROM_CONSVAR_1D(U2_j,NVAR1D))/DENSITY_FROM_CONSVAR(U2_j,NVAR1D),aux)

    ! Compute auxiliary variables
    q_ij = 0.5*(u_ij*u_ij)
#ifdef THERMALLY_IDEAL_GAS
    cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL)
#else
#error "Speed of sound must be implemented!"
#endif
    c_ij = sqrt(cPow2_ij)

    ! Compute eigenvalues
    l1 = abs(u_ij-c_ij)
    l2 = abs(u_ij)
    l3 = abs(u_ij+c_ij)

    ! Compute solution difference U2_i-U2_j
    Diff = U2_i-U2_j

    ! Compute auxiliary quantities for characteristic variables
    b2 = (GAMMA-1.0)/cPow2_ij; b1 = b2*q_ij

    ! Compute characteristic variables multiplied by the
    ! corresponding eigenvalue
    w1 = l1 * 0.5 * (     (b1+u_ij/c_ij)*Diff(1) -&
                      (b2*u_ij+1.0/c_ij)*Diff(2) +&
                                      b2*Diff(3) )
    w2 = l2 *          (        (1.0-b1)*Diff(1) +&
                                 b2*u_ij*Diff(2) -&
                                      b2*Diff(3) )
    w3 = l3 * 0.5 * (     (b1-u_ij/c_ij)*Diff(1) -&
                      (b2*u_ij-1.0/c_ij)*Diff(2) +&
                                      b2*Diff(3) )

    ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
    Diff(1) = w1 + w2 + w3
    Diff(2) = (u_ij-c_ij)*w1 + u_ij*w2 + (u_ij+c_ij)*w3
    Diff(3) = (H_ij-u_ij*c_ij)*w1 + q_ij*w2 + (H_ij+u_ij*c_ij)*w3

    ! Compute conservative fluxes
    F_ij = dscale1*(U1_i-U1_j) + dscale2*anorm*Diff

  end subroutine hydro_calcFluxFCTRoeDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRusanov1d(&
      U1_i, U1_j, U2_i, U2_j, Coeff_ij, Coeff_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 1D using the Rusanov dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: Coeff_ij,Coeff_ji

    ! scaling coefficients
    real(DP), intent(in) :: dscale1,dscale2

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! raw antidiffusive flux
    real(DP), dimension(:), intent(out) :: F_ij
!</output>
!</subroutine>

    ! local variables
    real(DP) :: d_ij,ui,uj,ci,cj,Ei,Ej

    ! Compute velocities and energy
    ui = X_VELOCITY_FROM_CONSVAR(U2_i,NVAR1D)
    Ei = TOTAL_ENERGY_FROM_CONSVAR(U2_i,NVAR1D)
    uj = X_VELOCITY_FROM_CONSVAR(U2_j,NVAR1D)
    Ej = TOTAL_ENERGY_FROM_CONSVAR(U2_j,NVAR1D)

    ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
    ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*ui*ui), SYS_EPSREAL))
    cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*uj*uj), SYS_EPSREAL))
#else
#error "Speed of sound must be implemented"
#endif

    ! Scalar dissipation
    d_ij = max( abs(Coeff_ij(1)*uj) + abs(Coeff_ij(1))*cj,&
                abs(Coeff_ji(1)*ui) + abs(Coeff_ji(1))*ci )

    ! Compute conservative fluxes
    F_ij = dscale1*(U1_i-U1_j) + dscale2*d_ij*(U2_i-U2_j)

  end subroutine hydro_calcFluxFCTRusanov1d

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDensity1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
    end do

  end subroutine hydro_trafoFluxDensity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDensity1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density in 1D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
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
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed density difference
      DtransformedDataAtEdge(1,idx) =&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine hydro_trafoDiffDensity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxEnergy1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed total energy fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
    end do

  end subroutine hydro_trafoFluxEnergy1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffEnergy1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the energy in 1D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
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
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed total density difference
      DtransformedDataAtEdge(1,idx) =&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine hydro_trafoDiffEnergy1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxPressure1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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
    real(DP) :: ui,uj
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(1,1,idx) = (GAMMA-1.0)*&
          (0.5*ui*ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
      DtransformedFluxesAtEdge(1,2,idx) =-(GAMMA-1.0)*&
          (0.5*uj*uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do
    
  end subroutine hydro_trafoFluxPressure1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffPressure1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the pressure in 1D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
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

    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(1,idx) =&
          PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine hydro_trafoDiffPressure1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxVelocity1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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
    real(DP) :: ui,uj

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(1,1,idx) =&
          (X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -(X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
    end do
    
  end subroutine hydro_trafoFluxVelocity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffVelocity1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the x-velocity
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
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
    
    do idx = 1, size(DdataAtEdge,3)

      ! Transformed velocity difference in x-direction
      DtransformedDataAtEdge(1,idx) =&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine hydro_trafoDiffVelocity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxMomentum1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed momentum fluxes in x-direction
      DtransformedFluxesAtEdge(1,1,idx) =&
          X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
    end do
    
  end subroutine hydro_trafoFluxMomentum1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffMomentum1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the x-momentum
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
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
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed momentum difference in x-direction
      DtransformedDataAtEdge(1,idx) =&
          X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do
    
  end subroutine hydro_trafoDiffMomentum1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenEng1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)

      ! Transformed total energy fluxes
      DtransformedFluxesAtEdge(2,1,idx) =&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
    end do

  end subroutine hydro_trafoFluxDenEng1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenEng1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density and energy in 1D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
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
    
    do idx = 1, size(DdataAtEdge,3)

      ! Transformed density difference
      DtransformedDataAtEdge(1,idx) =&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)

      ! Transformed total energy difference
      DtransformedDataAtEdge(2,idx) =&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine hydro_trafoDiffDenEng1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenPre1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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
    real(DP) :: ui,uj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
           
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)

      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(2,1,idx) = (GAMMA-1.0)*&
          (0.5*ui*ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
      DtransformedFluxesAtEdge(2,2,idx) =-(GAMMA-1.0)*&
          (0.5*uj*uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do

  end subroutine hydro_trafoFluxDenPre1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenPre1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density and energy in 1D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
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
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed density difference
      DtransformedDataAtEdge(1,idx) = &
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(2,idx) =&
          PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine hydro_trafoDiffDenPre1d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoFluxDenPreVel1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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
    real(DP) :: ui,uj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)

      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(2,1,idx) =&
          (X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -(X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(3,1,idx) = (GAMMA-1.0)*&
          (0.5*ui*ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
      DtransformedFluxesAtEdge(3,2,idx) =-(GAMMA-1.0)*&
          (0.5*uj*uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do
      
  end subroutine hydro_trafoFluxDenPreVel1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_trafoDiffDenPreVel1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density, pressure and velocity in 1D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
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

    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed density difference
      DtransformedDataAtEdge(1,idx) =&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      
      ! Transformed velocity difference in x-direction
      DtransformedDataAtEdge(2,idx) =&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(3,idx) =&
          PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine hydro_trafoDiffDenPreVel1d_sim

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
      p   = (GAMMA-1.0)*rho*(E-0.5*v1*v1)
      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))

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
      if (4.0/(GAMMA-1.0)*c .le. vn_b-vn) then
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
      ppv  = p+0.5*(vn-vn_b)*cup
      ppv  = max(0.0, ppv)

      if (ppv .eq. p) then

        ! Select guessed pressure from PVRS Riemann solver
        pstar = ppv
      else
        if (ppv .lt. p) then

          ! Guess pressure from the Two-Rarefaction Riemann solver
          vm    = 0.5*(vn+vn_b)
          ptl   = 1.0+(GAMMA-1.0)/2.0*(vn-vm)/c
          ptr   = 1.0+(GAMMA-1.0)/2.0*(vm-vn_b)/c
          pstar = 0.5*(p*ptl + p*ptr)**(2.0*GAMMA/(GAMMA-1.0))
        else

          ! Guess pressure from the Two-Shock Riemann solver
          ! with PVRS as estimated pressure value
          ge    = sqrt((2.0/(GAMMA+1.0)/rho)/((GAMMA-1.0)/(GAMMA+1.0)*p+ppv))
          pstar = p - 0.5*(vn_b-vn)/ge
        end if
      end if

      ! Initialize solution difference and pressure
      vdiff = (vn_b-vn)/2.0
      pold  = pstar

      newton: do ite = 1, 100

        ! Compute pressure function f(pold) and its derivative f1(pold)
        if (pold .le. p) then

          ! Rarefaction wave
          prat = pold/p

          f  = 2.0/(GAMMA-1.0)*c*(prat**((GAMMA-1.0)/(2.0*GAMMA))-1.0)
          fd = (1.0/(rho*c))*prat**(-(GAMMA+1.0)/(2.0*GAMMA))
        else

          ! Shock wave
          auxA = 2.0/(GAMMA+1.0)/rho
          auxB = (GAMMA-1.0)/(GAMMA+1.0)*p
          qrt  = sqrt(auxA/(auxB + pold))

          f  = (pold-p)*qrt
          fd = (1.0 - 0.5*(pold - p)/(auxB + pold))*qrt
        end if

        pstar = pold - (f+vdiff)/fd
        if (pstar .lt. 0.0) then
          pold = 1.0E-6
          cycle newton
        end if

        aux = 2.0*abs((pstar-pold)/(pstar+pold))
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
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcBoundaryvalues1d')
          call sys_halt()
        end if
      end if

      !-------------------------------------------------------------------------
      ! Calculate the velocity in the star region
      !-------------------------------------------------------------------------

      ! Note that the contribution fR-fL vanishes due to constant states
      vn = 0.5*(vn+vn_b)


      !-------------------------------------------------------------------------
      ! Calculate the density in the star region
      !-------------------------------------------------------------------------

      if (pstar .le. p) then

        ! Rarefaction wave
        rho = rho*(pstar/p)**(1.0/GAMMA)
      else

        ! Shock wave
        rho = rho*(pstar/p+(GAMMA-1.0)/(GAMMA+1.0))/&
                  ((GAMMA-1.0)/(GAMMA+1.0)*(pstar/p)+1.0)
      end if

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*DpointNormal(1)*vn
      Du(3) = pstar/(GAMMA-1.0)+0.5*rho*(vn*vn)

    case(BDRC_VISCOUSWALL)
      !-------------------------------------------------------------------------

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      E   = Du(3)/rho
      p   = (GAMMA-1.0)*rho*(E-0.5*v1*v1)

      ! Update the solution vector and let vn:=0
      Du(2) = 0.0
      Du(3) = p/(GAMMA-1.0)

    case(BDRC_SUPERINLET)
      !-------------------------------------------------------------------------

      ! The free stream primitive variables are Deval=[rho,v1,p]
      rho = DbdrValue(1)
      p   = DbdrValue(3)
      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*DbdrValue(2)

      ! Compute Riemann invariants based on the free stream values
      W(1) = vn-2*c/(GAMMA-1.0)
      W(2) = vn+2*c/(GAMMA-1.0)
      W(3) = p/(rho**GAMMA)

      ! Transform back into conservative variables
      vn   = 0.5*(W(1)+W(2))
      c    = 0.25*(GAMMA-1.0)*(W(2)-W(1))
      rho  = (c*c/GAMMA/W(3))**(1.0/(GAMMA-1.0))
      p    = rho*c*c/GAMMA

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*DbdrNormal(1)*vn
      Du(3) = p/(GAMMA-1.0)+0.5*rho*(vn*vn)


    case(BDRC_FREESTREAM)
      !-------------------------------------------------------------------------

      ! The free stream primitive variables are Deval=[rho,v1,p]
      rho = DbdrValue(1)
      p   = DbdrValue(3)
      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*DbdrValue(2)

      ! Compute Riemann invariants based on the free stream values
      Winf(1) = vn-2.0*c/(GAMMA-1.0)
      Winf(2) = vn+2.0*c/(GAMMA-1.0)
      Winf(3) = p/(rho**GAMMA)

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      E   = Du(3)/rho
      p   = (GAMMA-1.0)*rho*(E-0.5*v1*v1)

      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*v1

      ! Compute Riemann invariants based on the solution values
      Wu(1) = vn-2.0*c/(GAMMA-1.0)
      Wu(2) = vn+2.0*c/(GAMMA-1.0)
      Wu(3) = p/(rho**GAMMA)

      ! Adopt free stream/computed values depending on the sign of the eigenvalue
      W(1) = merge(Winf(1), Wu(1), vn <  c)
      W(2) = merge(Winf(2), Wu(2), vn < -c)
      W(3) = merge(Winf(3), Wu(3), vn <  SYS_EPSREAL)

      ! Transform back into conservative variables
      vn   = 0.5*(W(1)+W(2))
      c    = 0.25*(GAMMA-1.0)*(W(2)-W(1))
      rho  = (c*c/GAMMA/W(3))**(1.0/(GAMMA-1.0))
      p    = rho*c*c/GAMMA

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*DbdrNormal(1)*vn
      Du(3) = p/(GAMMA-1.0)+0.5*rho*(vn*vn)


    case(BDRC_SUBINLET)
      !-------------------------------------------------------------------------

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      E   = Du(3)/rho
      p   = (GAMMA-1.0)*rho*(E-0.5*v1*v1)

      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*v1

      ! The specified density and pressure is Deval=[rho,p]
      rho = DbdrValue(1)
      p   = DbdrValue(2)

      ! Compute Riemann invariants
      W(1) = vn-2.0*c/(GAMMA-1.0)
      W(2) = vn+2.0*c/(GAMMA-1.0)
      W(3) = p/(rho**GAMMA)

      ! Transform back into conservative variables
      vn   = 0.5*(W(1)+W(2))
      c    = 0.25*(GAMMA-1.0)*(W(2)-W(1))
      rho  = (c*c/GAMMA/W(3))**(1.0/(GAMMA-1.0))
      p    = rho*c*c/GAMMA

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*DbdrNormal(1)*vn
      Du(3) = p/(GAMMA-1.0)+0.5*rho*(vn*vn)


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
      p   = (GAMMA-1.0)*rho*(E-0.5*v1*v1)

      vn  = DbdrNormal(1)*v1
      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))

      ! Compute Riemann invariants based on the solution values and prescribed exit pressure
      W(2) = 2*c/(GAMMA-1.0)-vn
      W(3) = p/(rho**GAMMA)
      W(1) = 4/(GAMMA-1.0)*sqrt(max(GAMMA*ps/rho*(p/ps)**(1.0/GAMMA), SYS_EPSREAL))-W(2)

      ! Transform back into conservative variables
      vn  = 0.5*(W(1)-W(2))
      c   = 0.25*(GAMMA-1.0)*(W(1)+W(2))
      rho = (c*c/GAMMA/W(3))**(1.0/(GAMMA-1.0))
      p   = rho*c*c/GAMMA

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*DbdrNormal(1)*vn
      Du(3) = p/(GAMMA-1.0)+0.5*rho*(vn*vn)


    case DEFAULT
      call output_line('Unsupported type of boundary condition!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcBoundaryvalues1d')
      call sys_halt()
    end select

  end subroutine hydro_calcBoundaryvalues1d

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
            0.5*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR1D+ivar)+&
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
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR1D+ivar) = 0.0
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
            0.5*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
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
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = 0.0
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
