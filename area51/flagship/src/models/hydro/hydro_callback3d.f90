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
!# 1.) hydro_calcFluxGal3d_sim
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
!# 12.) hydro_calcMatGal3d_sim
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
!# 25.) hydro_trafoFluxEnergy3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 26.) hydro_trafoDiffEnergy3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 27.) hydro_trafoFluxPressure3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 28.) hydro_trafoDiffPressure3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the pressure
!#
!# 29.) hydro_trafoFluxVelocity3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 30.) hydro_trafoDiffVelocity3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 31.) hydro_trafoFluxMomentum3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 32.) hydro_trafoDiffMomentum3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 33.) hydro_trafoFluxDenEng3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 34.) hydro_trafoDiffDenEng3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 35.) hydro_trafoFluxDenPre3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 36.) hydro_trafoDiffDenPre3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 37.) hydro_trafoFluxDenPreVel3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 38.) hydro_trafoDiffDenPreVel3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure 
!#         and the velocity
!#
!# 39.) hydro_calcBoundaryvalues3d
!#      -> Computes the boundary values for a given node
!#
!# 40.) hydro_hadaptCallbackScalar3d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 3D, whereby the vector is stored in interleave format
!#
!# 41.) hydro_hadaptCallbackBlock3d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 3D, whereby the vector is stored in block format
!#
!# </purpose>
!##############################################################################

module hydro_callback3d

#include "hydro.h"

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
  public :: hydro_calcFluxGal3d_sim
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
  public :: hydro_calcMatGal3d_sim
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
  public :: hydro_calcBoundaryvalues3d
  public :: hydro_hadaptCallbackScalar3d
  public :: hydro_hadaptCallbackBlock3d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGal3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)

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
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_HYDRO_2T_XDIR_3D(Fxi,DdataAtEdge,1,idx,ui,pi)
      FLUX_HYDRO_2T_XDIR_3D(Fxj,DdataAtEdge,2,idx,uj,pj)
      
      ! Compute fluxes for y-direction
      FLUX_HYDRO_2T_YDIR_3D(Fyi,DdataAtEdge,1,idx,vi,pi)
      FLUX_HYDRO_2T_YDIR_3D(Fyj,DdataAtEdge,2,idx,vj,pj)

      ! Compute fluxes for z-direction
      FLUX_HYDRO_2T_ZDIR_3D(Fzi,DdataAtEdge,1,idx,wi,pi)
      FLUX_HYDRO_2T_ZDIR_3D(Fzj,DdataAtEdge,2,idx,wj,pj)

      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                         DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                         DmatrixCoeffsAtEdge(3,1,idx)*Fzi )
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_HYDRO_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,pi,pj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_HYDRO_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,vi,vj,pi,pj)

      ! Compute flux difference for z-direction
      FLUXDIFF_HYDRO_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,wi,wj,pi,pj)

      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij)
#endif
    end do

  end subroutine hydro_calcFluxGal3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxGalNoBdr3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)

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
  real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
  real(DP) :: pi,pj,ui,uj,vi,vj,wi,wj
  integer :: idx
  
    
    do idx = 1, nedges

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

      ! Compute flux difference for x-direction
      FLUXDIFF_HYDRO_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,pi,pj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_HYDRO_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,vi,vj,pi,pj)

      ! Compute flux difference for z-direction
      FLUXDIFF_HYDRO_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,wi,wj,pi,pj)

      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale *&
          (0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))*Fx_ij+&
           0.5*(DmatrixCoeffsAtEdge(2,1,idx)-DmatrixCoeffsAtEdge(2,2,idx))*Fy_ij+&
           0.5*(DmatrixCoeffsAtEdge(3,1,idx)-DmatrixCoeffsAtEdge(3,2,idx))*Fz_ij)
      DfluxesAtEdge(:,2,idx) = DfluxesAtEdge(:,1,idx)
    end do

  end subroutine hydro_calcFluxGalNoBdr3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxScDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)
    
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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: pi,pj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,anorm,aux,c_ij,d_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    integer :: idx

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_HYDRO_2T_XDIR_3D(Fxi,DdataAtEdge,1,idx,ui,pi)
      FLUX_HYDRO_2T_XDIR_3D(Fxj,DdataAtEdge,2,idx,uj,pj)
      
      ! Compute fluxes for y-direction
      FLUX_HYDRO_2T_YDIR_3D(Fyi,DdataAtEdge,1,idx,vi,pi)
      FLUX_HYDRO_2T_YDIR_3D(Fyj,DdataAtEdge,2,idx,vj,pj)

      ! Compute fluxes for z-direction
      FLUX_HYDRO_2T_ZDIR_3D(Fzi,DdataAtEdge,1,idx,wi,pi)
      FLUX_HYDRO_2T_ZDIR_3D(Fzj,DdataAtEdge,2,idx,wj,pj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_HYDRO_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,pi,pj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_HYDRO_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,vi,vj,pi,pj)

      ! Compute flux difference for z-direction
      FLUXDIFF_HYDRO_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,wi,wj,pi,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      w_ij = ROE_MEAN_VALUE(wi,wj,aux)
      H_ij = ROE_MEAN_VALUE(\
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)+pi)/\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)+pj)/\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx),aux)

      ! Compute auxiliary variables
      vel_ij = u_ij*a(1) + v_ij*a(2) + w_ij*a(3)
      q_ij   = 0.5*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL_DP))
#else
#error "Speed of sound must be implemented!"
#endif

      ! Compute scalar dissipation
      d_ij = abs(vel_ij) + anorm*c_ij

      ! Multiply the solution difference by the scalar dissipation
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                         DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                         DmatrixCoeffsAtEdge(3,1,idx)*Fzi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij+ Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij+ Diff)
#endif
    end do

  end subroutine hydro_calcFluxScDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxScDissDiSp3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)
    

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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: pi,pj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,aux,c_ij,d_ij,q_ij,u_ij,v_ij,w_ij
    integer :: idx

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_HYDRO_2T_XDIR_3D(Fxi,DdataAtEdge,1,idx,ui,pi)
      FLUX_HYDRO_2T_XDIR_3D(Fxj,DdataAtEdge,2,idx,uj,pj)
      
      ! Compute fluxes for y-direction
      FLUX_HYDRO_2T_YDIR_3D(Fyi,DdataAtEdge,1,idx,vi,pi)
      FLUX_HYDRO_2T_YDIR_3D(Fyj,DdataAtEdge,2,idx,vj,pj)

      ! Compute fluxes for z-direction
      FLUX_HYDRO_2T_ZDIR_3D(Fzi,DdataAtEdge,1,idx,wi,pi)
      FLUX_HYDRO_2T_ZDIR_3D(Fzj,DdataAtEdge,2,idx,wj,pj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_HYDRO_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,pi,pj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_HYDRO_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,vi,vj,pi,pj)

      ! Compute flux difference for z-direction
      FLUXDIFF_HYDRO_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,wi,wj,pi,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !-------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      
      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      w_ij = ROE_MEAN_VALUE(vi,vj,aux)
      H_ij = ROE_MEAN_VALUE(\
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)+pi)/\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)+pj)/\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx),aux)

      ! Compute auxiliary variable
      q_ij = 0.5*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL_DP))
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Compute scalar dissipation with dimensional splitting
      d_ij = ( abs(a(1)*u_ij) + abs(a(1))*c_ij +&
               abs(a(2)*v_ij) + abs(a(2))*c_ij +&
               abs(a(3)*w_ij) + abs(a(3))*c_ij )

      ! Multiply the solution difference by the scalar dissipation
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                         DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                         DmatrixCoeffsAtEdge(3,1,idx)*Fzi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxScDissDiSp3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRoeDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)

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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: Diff
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: pi,pj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,c2_ij,c_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    real(DP) :: anorm,aux,aux1,aux2
    real(DP) :: l1,l2,l3,l4,l5,w1,w2,w3,w4,w5
    integer :: idx

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_HYDRO_2T_XDIR_3D(Fxi,DdataAtEdge,1,idx,ui,pi)
      FLUX_HYDRO_2T_XDIR_3D(Fxj,DdataAtEdge,2,idx,uj,pj)
      
      ! Compute fluxes for y-direction
      FLUX_HYDRO_2T_YDIR_3D(Fyi,DdataAtEdge,1,idx,vi,pi)
      FLUX_HYDRO_2T_YDIR_3D(Fyj,DdataAtEdge,2,idx,vj,pj)

      ! Compute fluxes for z-direction
      FLUX_HYDRO_2T_ZDIR_3D(Fzi,DdataAtEdge,1,idx,wi,pi)
      FLUX_HYDRO_2T_ZDIR_3D(Fzj,DdataAtEdge,2,idx,wj,pj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_HYDRO_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,pi,pj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_HYDRO_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,vi,vj,pi,pj)

      ! Compute flux difference for z-direction
      FLUXDIFF_HYDRO_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,wi,wj,pi,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Normalize the skew-symmetric coefficient
        a = a/anorm
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)+pi)/\
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)+pj)/\
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx),aux)

        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2) + w_ij*a(3)
        q_ij   = 0.5*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

        ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
        c2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL_DP)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(c2_ij)

        ! Compute eigenvalues
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)
        l5 = abs(vel_ij)

        ! Compute solution difference U_j-U_i
        Diff = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = (GAMMA-1.0)/2.0/c2_ij*(q_ij*Diff(1)&
                                     -u_ij*Diff(2)&
                                     -v_ij*Diff(3)&
                                     -w_ij*Diff(4)&
                                          +Diff(5))
        aux2 = (vel_ij*Diff(1)&
                 -a(1)*Diff(2)&
                 -a(2)*Diff(3)&
                 -a(3)*Diff(4))/2.0/c_ij

        ! Get the dimension with largest coefficient
        select case(maxloc(a,1))
        case(1)
          ! Compute characteristic variables multiplied by the corresponding eigenvalue
          w1 = l1 * (aux1 + aux2)
          w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/c2_ij)*Diff(1)&
                              +(GAMMA-1.0)*( u_ij*Diff(2)&
                                            +v_ij*Diff(3)&
                                            +w_ij*Diff(4)&
                                                 -Diff(5))/c2_ij)
          w3 = l3 * (aux1 - aux2)
          w4 = l4 * ( (v_ij-vel_ij*a(2))/a(1)*Diff(1)&
                                        +a(2)*Diff(2)&
                        +(a(2)*a(2)-1.0)/a(1)*Diff(3)&
                              +a(2)*a(3)/a(1)*Diff(4))
          w5 = l5 * ( (vel_ij*a(3)-w_ij)/a(1)*Diff(1)&
                                        -a(3)*Diff(2)&
                              -a(2)*a(3)/a(1)*Diff(3)&
                        +(1.0-a(3)*a(3))/a(1)*Diff(4))

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
          w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/c2_ij)*Diff(1)&
                               +(GAMMA-1.0)*(u_ij*Diff(2)&
                                            +v_ij*Diff(3)&
                                            +w_ij*Diff(4)&
                                                 -Diff(5))/c2_ij)
          w3 = l3 * (aux1 - aux2)
          w4 = l4 * ( (vel_ij*a(1)-u_ij)/a(2)*Diff(1)&
                        +(1.0-a(1)*a(1))/a(2)*Diff(2)&
                                        -a(1)*Diff(3)&
                              -a(1)*a(3)/a(2)*Diff(4))
          w5 = l5 * ( (w_ij-vel_ij*a(3))/a(2)*Diff(1)&
                              +a(1)*a(3)/a(2)*Diff(2)&
                                        +a(3)*Diff(3)&
                        +(a(3)*a(3)-1.0)/a(2)*Diff(4))

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
          w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/c2_ij)*Diff(1)&
                               +(GAMMA-1.0)*(u_ij*Diff(2)&
                                            +v_ij*Diff(3)&
                                            +w_ij*Diff(4)&
                                                 -Diff(5))/c2_ij)
          w3 = l3 * (aux1 - aux2)
          w4 = l4 * ( (u_ij-vel_ij*a(1))/a(3)*Diff(1)&
                        +(a(1)*a(1)-1.0)/a(3)*Diff(2)&
                              +a(1)*a(2)/a(3)*Diff(3)&
                                        +a(1)*Diff(4) )
          w5 = l5 * ( (vel_ij*a(2)-v_ij)/a(3)*Diff(1)&
                              -a(1)*a(2)/a(3)*Diff(2)&
                        +(1.0-a(2)*a(2))/a(3)*Diff(3)&
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
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                           DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                           DmatrixCoeffsAtEdge(3,1,idx)*Fzi + Diff)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij + Diff)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij + Diff)
#endif
      else
        
#ifdef HYDRO_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                           DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                           DmatrixCoeffsAtEdge(3,1,idx)*Fzi)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij)
#endif
      end if
    end do

  end subroutine hydro_calcFluxRoeDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRoeDissDiSp3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)
    

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
    real(DP), dimension(NVAR3D) :: Fxi,Fxj,Fyi,Fyj,Fzi,Fzj
#else
    real(DP), dimension(NVAR3D) :: Fx_ij,Fy_ij,Fz_ij
#endif
    real(DP), dimension(NVAR3D) :: DiffX,DiffY,DiffZ
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: pi,pj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,c2_ij,c_ij,q_ij,u_ij,v_ij,w_ij
    real(DP) :: anorm,aux,aux1,aux2
    real(DP) :: l1,l2,l3,l4,l5,w1,w2,w3,w4,w5
    integer :: idx

    
    do idx = 1, nedges
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_HYDRO_2T_XDIR_3D(Fxi,DdataAtEdge,1,idx,ui,pi)
      FLUX_HYDRO_2T_XDIR_3D(Fxj,DdataAtEdge,2,idx,uj,pj)
      
      ! Compute fluxes for y-direction
      FLUX_HYDRO_2T_YDIR_3D(Fyi,DdataAtEdge,1,idx,vi,pi)
      FLUX_HYDRO_2T_YDIR_3D(Fyj,DdataAtEdge,2,idx,vj,pj)

      ! Compute fluxes for z-direction
      FLUX_HYDRO_2T_ZDIR_3D(Fzi,DdataAtEdge,1,idx,wi,pi)
      FLUX_HYDRO_2T_ZDIR_3D(Fzj,DdataAtEdge,2,idx,wj,pj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_HYDRO_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,pi,pj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_HYDRO_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,vi,vj,pi,pj)

      ! Compute flux difference for z-direction
      FLUXDIFF_HYDRO_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,wi,wj,pi,pj)
#endif
      
      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Normalize the skew-symmetric coefficient
        a = abs(a)
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)+pi)/\
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)+pj)/\
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx),aux)

        ! Compute auxiliary variable
        q_ij = 0.5*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

        ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
        c2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL_DP)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij = sqrt(c2_ij)

        !-----------------------------------------------------------------------
        ! Dimensional splitting: x-direction
        !-----------------------------------------------------------------------

        ! Compute eigenvalues
        l1 = abs(u_ij-c_ij)
        l2 = abs(u_ij)
        l3 = abs(u_ij+c_ij)
        l4 = abs(u_ij)
        l5 = abs(u_ij)
        
        ! Compute solution difference U_j-U_i
        DiffX = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = (GAMMA-1.0)/2.0/c2_ij*(q_ij*DiffX(1)&
                                     -u_ij*DiffX(2)&
                                     -v_ij*DiffX(3)&
                                     -w_ij*DiffX(4)&
                                          +DiffX(5))
        aux2 = (u_ij*DiffX(1)-DiffX(2))/2.0/c_ij

        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/c2_ij)*DiffX(1)&
                            +(GAMMA-1.0)*( u_ij*DiffX(2)&
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

        ! Compute solution difference U_j-U_i
        DiffY = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)

        ! Compute auxiliary quantities for characteristic variables
        aux1 = (GAMMA-1.0)/2.0/c2_ij*(q_ij*DiffY(1)&
                                     -u_ij*DiffY(2)&
                                     -v_ij*DiffY(3)&
                                     -w_ij*DiffY(4)&
                                          +DiffY(5))
        aux2 = (v_ij*DiffY(1)-DiffY(3))/2.0/c_ij

        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/c2_ij)*DiffY(1)&
                             +(GAMMA-1.0)*(u_ij*DiffY(2)&
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

        ! Compute auxiliary quantities for characteristic variables
        aux1 = (GAMMA-1.0)/2.0/c2_ij*(q_ij*DiffZ(1)&
                                     -u_ij*DiffZ(2)&
                                     -v_ij*DiffZ(3)&
                                     -w_ij*DiffZ(4)&
                                          +DiffZ(5))
        aux2 = (w_ij*DiffZ(1)-DiffZ(3))/2.0/c_ij

        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/c2_ij)*DiffZ(1)&
                             +(GAMMA-1.0)*(u_ij*DiffZ(2)&
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
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                           DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                           DmatrixCoeffsAtEdge(3,1,idx)*Fzi+&
                                           DiffX+DiffY+DiffZ)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij+&
                                            DiffX+DiffY+DiffZ)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij+&
                                            DiffX+DiffY+DiffZ)
#endif
      else
        
#ifdef HYDRO_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                           DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                           DmatrixCoeffsAtEdge(3,1,idx)*Fzi)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                            DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij)
#endif
      end if
    end do

  end subroutine hydro_calcFluxRoeDissDiSp3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRusDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)
    
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
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

      ! Compute specific energies
      Ei = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      Ej = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_HYDRO_2T_XDIR_3D(Fxi,DdataAtEdge,1,idx,ui,pi)
      FLUX_HYDRO_2T_XDIR_3D(Fxj,DdataAtEdge,2,idx,uj,pj)
      
      ! Compute fluxes for y-direction
      FLUX_HYDRO_2T_YDIR_3D(Fyi,DdataAtEdge,1,idx,vi,pi)
      FLUX_HYDRO_2T_YDIR_3D(Fyj,DdataAtEdge,2,idx,vj,pj)

      ! Compute fluxes for z-direction
      FLUX_HYDRO_2T_ZDIR_3D(Fzi,DdataAtEdge,1,idx,wi,pi)
      FLUX_HYDRO_2T_ZDIR_3D(Fzj,DdataAtEdge,2,idx,wj,pj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_HYDRO_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,pi,pj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_HYDRO_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,vi,vj,pi,pj)

      ! Compute flux difference for z-direction
      FLUXDIFF_HYDRO_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,wi,wj,pi,pj)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*(ui*ui+vi*vi+wi*wi)), SYS_EPSREAL_DP))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*(uj*uj+vj*vj+wj*wj)), SYS_EPSREAL_DP))
#else
#error "Speed of sound must be implemented!"
#endif

#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(0.5_DP*(DmatrixCoeffsAtEdge(1,1,idx)-&
                              DmatrixCoeffsAtEdge(1,2,idx))*uj+&
                      0.5_DP*(DmatrixCoeffsAtEdge(2,1,idx)-&
                              DmatrixCoeffsAtEdge(2,2,idx))*vj+&
                      0.5_DP*(DmatrixCoeffsAtEdge(3,1,idx)-&
                              DmatrixCoeffsAtEdge(3,2,idx))*wj)+&
                 0.5_DP*sqrt((DmatrixCoeffsAtEdge(1,1,idx)-&
                              DmatrixCoeffsAtEdge(1,2,idx))**2+&
                             (DmatrixCoeffsAtEdge(2,1,idx)-&
                              DmatrixCoeffsAtEdge(2,2,idx))**2+&
                             (DmatrixCoeffsAtEdge(3,1,idx)-&
                              DmatrixCoeffsAtEdge(3,2,idx))**2)*cj,&
                  abs(0.5_DP*(DmatrixCoeffsAtEdge(1,2,idx)-&
                              DmatrixCoeffsAtEdge(1,1,idx))*ui+&
                      0.5_DP*(DmatrixCoeffsAtEdge(2,2,idx)-&
                              DmatrixCoeffsAtEdge(2,1,idx))*vi+&
                      0.5_DP*(DmatrixCoeffsAtEdge(3,2,idx)-&
                              DmatrixCoeffsAtEdge(3,1,idx))*wi)+&
                 0.5_DP*sqrt((DmatrixCoeffsAtEdge(1,2,idx)-&
                              DmatrixCoeffsAtEdge(1,1,idx))**2+&
                             (DmatrixCoeffsAtEdge(2,2,idx)-&
                              DmatrixCoeffsAtEdge(2,1,idx))**2+&
                             (DmatrixCoeffsAtEdge(3,2,idx)-&
                              DmatrixCoeffsAtEdge(3,1,idx))**2)*ci )
#else      
      ! Compute scalar dissipation
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj+&
                      DmatrixCoeffsAtEdge(2,1,idx)*vj+&
                      DmatrixCoeffsAtEdge(3,1,idx)*wj)+&
                 sqrt(DmatrixCoeffsAtEdge(1,1,idx)**2+&
                      DmatrixCoeffsAtEdge(2,1,idx)**2+&
                      DmatrixCoeffsAtEdge(3,1,idx)**2)*cj,&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui+&
                      DmatrixCoeffsAtEdge(2,2,idx)*vi+&
                      DmatrixCoeffsAtEdge(3,2,idx)*wi)+&
                 sqrt(DmatrixCoeffsAtEdge(1,2,idx)**2+&
                      DmatrixCoeffsAtEdge(2,2,idx)**2+&
                      DmatrixCoeffsAtEdge(3,2,idx)**2)*ci )
#endif

      ! Multiply the solution difference by the scalar dissipation
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------
      
#ifdef HYDRO_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                         DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                         DmatrixCoeffsAtEdge(3,1,idx)*Fzi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxRusDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxRusDissDiSp3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)

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
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

      ! Compute specific energies
      Ei = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      Ej = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

#ifdef HYDRO_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_HYDRO_2T_XDIR_3D(Fxi,DdataAtEdge,1,idx,ui,pi)
      FLUX_HYDRO_2T_XDIR_3D(Fxj,DdataAtEdge,2,idx,uj,pj)
      
      ! Compute fluxes for y-direction
      FLUX_HYDRO_2T_YDIR_3D(Fyi,DdataAtEdge,1,idx,vi,pi)
      FLUX_HYDRO_2T_YDIR_3D(Fyj,DdataAtEdge,2,idx,vj,pj)

      ! Compute fluxes for z-direction
      FLUX_HYDRO_2T_ZDIR_3D(Fzi,DdataAtEdge,1,idx,wi,pi)
      FLUX_HYDRO_2T_ZDIR_3D(Fzj,DdataAtEdge,2,idx,wj,pj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_HYDRO_2T_XDIR_3D(Fx_ij,DdataAtEdge,1,2,idx,ui,uj,pi,pj)
                        
      ! Compute flux difference for y-direction
      FLUXDIFF_HYDRO_2T_YDIR_3D(Fy_ij,DdataAtEdge,1,2,idx,vi,vj,pi,pj)

      ! Compute flux difference for z-direction
      FLUXDIFF_HYDRO_2T_ZDIR_3D(Fz_ij,DdataAtEdge,1,2,idx,wi,wj,pi,pj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !-------------------------------------------------------------------------

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*(ui*ui+vi*vi+wi*wi)), SYS_EPSREAL_DP))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*(uj*uj+vj*vj+wj*wj)), SYS_EPSREAL_DP))
#else
#error "Speed of sound must be implemented!"
#endif
      
#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation with dimensional splitting based on
      ! the skew-symmetric part which does not include the symmetric
      ! boundary contribution
      d_ij = max( abs(0.5_DP*(DmatrixCoeffsAtEdge(1,1,idx)-&
                              DmatrixCoeffsAtEdge(1,2,idx))*uj)+&
                  abs(0.5_DP*(DmatrixCoeffsAtEdge(1,1,idx)-&
                              DmatrixCoeffsAtEdge(1,2,idx)))*cj,&
                  abs(0.5_DP*(DmatrixCoeffsAtEdge(1,2,idx)-&
                              DmatrixCoeffsAtEdge(1,1,idx))*ui)+&
                  abs(0.5_DP*(DmatrixCoeffsAtEdge(1,2,idx)-&
                              DmatrixCoeffsAtEdge(1,1,idx)))*ci )&
           + max( abs(0.5_DP*(DmatrixCoeffsAtEdge(2,1,idx)-&
                              DmatrixCoeffsAtEdge(2,2,idx))*vj)+&
                  abs(0.5_DP*(DmatrixCoeffsAtEdge(2,1,idx)-&
                              DmatrixCoeffsAtEdge(2,2,idx)))*cj,&
                  abs(0.5_DP*(DmatrixCoeffsAtEdge(2,2,idx)-&
                              DmatrixCoeffsAtEdge(2,1,idx))*vi)+&
                  abs(0.5_DP*(DmatrixCoeffsAtEdge(2,2,idx)-&
                              DmatrixCoeffsAtEdge(2,1,idx)))*ci )&
           + max( abs(0.5_DP*(DmatrixCoeffsAtEdge(3,1,idx)-&
                              DmatrixCoeffsAtEdge(3,2,idx))*wj)+&
                  abs(0.5_DP*(DmatrixCoeffsAtEdge(3,1,idx)-&
                              DmatrixCoeffsAtEdge(3,2,idx)))*cj,&
                  abs(0.5_DP*(DmatrixCoeffsAtEdge(3,2,idx)-&
                              DmatrixCoeffsAtEdge(3,1,idx))*wi)+&
                  abs(0.5_DP*(DmatrixCoeffsAtEdge(3,2,idx)-&
                              DmatrixCoeffsAtEdge(3,1,idx)))*ci )
#else
      ! Compute scalar dissipation with dimensional splitting
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj)+&
                  abs(DmatrixCoeffsAtEdge(1,1,idx))*cj,&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui)+&
                  abs(DmatrixCoeffsAtEdge(1,2,idx))*ci )&
           + max( abs(DmatrixCoeffsAtEdge(2,1,idx)*vj)+&
                  abs(DmatrixCoeffsAtEdge(2,1,idx))*cj,&
                  abs(DmatrixCoeffsAtEdge(2,2,idx)*vi)+&
                  abs(DmatrixCoeffsAtEdge(2,2,idx))*ci )&
           + max( abs(DmatrixCoeffsAtEdge(3,1,idx)*wj)+&
                  abs(DmatrixCoeffsAtEdge(3,1,idx))*cj,&
                  abs(DmatrixCoeffsAtEdge(3,2,idx)*wi)+&
                  abs(DmatrixCoeffsAtEdge(3,2,idx))*ci )
#endif

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))
      
      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef HYDRO_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fxj+&
                                         DmatrixCoeffsAtEdge(2,2,idx)*Fyj+&
                                         DmatrixCoeffsAtEdge(3,2,idx)*Fzj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fxi-&
                                         DmatrixCoeffsAtEdge(2,1,idx)*Fyi-&
                                         DmatrixCoeffsAtEdge(3,1,idx)*Fzi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,1,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,1,idx)*Fz_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fx_ij+&
                                          DmatrixCoeffsAtEdge(2,2,idx)*Fy_ij+&
                                          DmatrixCoeffsAtEdge(3,2,idx)*Fz_ij + Diff)
#endif
    end do

  end subroutine hydro_calcFluxRusDissDiSp3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiagMatD3d_sim(DdataAtNode, DmatrixCoeffsAtNode,&
      IverticesAtNode, dscale, nnodes, DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 3D.
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
    real(DP) :: ui,vi,wi
    integer :: inode


    do inode = 1, nnodes
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_1T_FROM_CONSVAR(DdataAtNode,NVAR3D,inode)
      vi = Y_VELOCITY_1T_FROM_CONSVAR(DdataAtNode,NVAR3D,inode)
      wi = Z_VELOCITY_1T_FROM_CONSVAR(DdataAtNode,NVAR3D,inode)
      
#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ii = diag(A_i)*C_{ii}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtNode,1,inode,\
        dscale,DmatrixCoeffsAtNode(1,inode),DmatrixCoeffsAtNode(2,inode),\
        DmatrixCoeffsAtNode(3,inode),ui,vi,wi)
#else
      ! Compute Galerkin coefficient $K_ii = -diag(A_i)*C_{ii}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtNode,1,inode,\
        -dscale,DmatrixCoeffsAtNode(1,inode),DmatrixCoeffsAtNode(2,inode),\
        DmatrixCoeffsAtNode(3,inode),ui,vi,wi)
#endif
    end do

  end subroutine hydro_calcMatDiagMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatDiag3d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 3D.
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
    real(DP) :: Ei,ui,vi,wi
    integer :: inode


    do inode = 1, nnodes
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_1T_FROM_CONSVAR(DdataAtNode,NVAR3D,inode)
      vi = Y_VELOCITY_1T_FROM_CONSVAR(DdataAtNode,NVAR3D,inode)
      wi = Z_VELOCITY_1T_FROM_CONSVAR(DdataAtNode,NVAR3D,inode)
      Ei = SPECIFIC_TOTAL_ENERGY_1T_FROM_CONSVAR(DdataAtNode,NVAR3D,inode)

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ii = A_i*C_{ii}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtNode,1,inode,\
      dscale,DmatrixCoeffsAtNode(1,inode),DmatrixCoeffsAtNode(2,inode),\
      DmatrixCoeffsAtNode(3,inode),ui,vi,wi,Ei)
#else
      ! Compute Galerkin coefficient $K_ii = -A_i*C_{ii}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtNode,1,inode,\
        -dscale,DmatrixCoeffsAtNode(1,inode),DmatrixCoeffsAtNode(2,inode),\
        DmatrixCoeffsAtNode(3,inode),ui,vi,wi,Ei)
#endif
    end do

  end subroutine hydro_calcMatDiag3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGalMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 3D.
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
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),uj,vj,wj)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),ui,vi,wi)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),uj,vj,wj)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),ui,vi,wi)
#endif
    end do

  end subroutine hydro_calcMatGalMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatGal3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D.
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
    real(DP) :: Ei,Ej,ui,uj,vi,vj,wi,wj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      Ei = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      Ej = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),ui,vi,wi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),ui,vi,wi,Ei)
#endif
    end do

  end subroutine hydro_calcMatGal3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDissMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 3D and applies scalar artificial viscosities proportional to
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
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,anorm,aux,c_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    integer :: idx

    do idx = 1, nedges
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),uj,vj,wj)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),ui,vi,wi)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),uj,vj,wj)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),ui,vi,wi)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      
      if (anorm .gt. SYS_EPSREAL_DP) then

        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)+\
               PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)+\
               PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx),aux)

        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2) + w_ij*a(3)
        q_ij   = 0.5*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

        ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
        c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL_DP))
#else
#error "Speed of sound must be implemented!"
#endif

        ! Compute scalar dissipation
        DcoefficientsAtEdge(:,1,idx) = dscale * (abs(vel_ij) + anorm*c_ij)

      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0

      end if
    end do

  end subroutine hydro_calcMatScDissMatD3d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatScDiss3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D and applies
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
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: Ei,Ej,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,anorm,aux,c_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    integer :: idx


    do idx = 1, nedges
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      Ei = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      Ej = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),ui,vi,wi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),ui,vi,wi,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      
      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)+\
               PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)+\
               PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx),aux)
        
        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2) + w_ij*a(3)
        q_ij   = 0.5*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

        ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
        c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL_DP))
#else
#error "Speed of sound must be implemented!"
#endif

        ! Compute scalar dissipation
        aux = dscale * (abs(vel_ij) + anorm*c_ij)
              
        DcoefficientsAtEdge( 1,1,idx) = aux
        DcoefficientsAtEdge( 6,1,idx) = aux
        DcoefficientsAtEdge(11,1,idx) = aux
        DcoefficientsAtEdge(16,1,idx) = aux
      end if
    end do

  end subroutine hydro_calcMatScDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDissMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D and applies
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
    real(DP), dimension(NVAR3D,NVAR3D) :: R_ij,L_ij
    real(DP), dimension(NVAR3D) :: a
    real(DP) :: ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,anorm,aux,cPow2_ij,c_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    real(DP) :: l1,l2,l3,l4,l5
    integer :: idx


    do idx = 1, nedges
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),uj,vj,wj)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),ui,vi,wi)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),uj,vj,wj)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),ui,vi,wi)
#endif
        
       ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      if (anorm .gt. SYS_EPSREAL_DP) then

        ! Normalize the skew-symmetric coefficient
        a = a/anorm

        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)+\
               PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)+\
               PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx),aux)

        ! Compute auxiliary values
        vel_ij = u_ij*a(1)+v_ij*a(2)+w_ij*a(3)
        q_ij   = 0.5*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)
        
        ! Compute speed of sound
#ifdef THERMALLY_IDEAL_GAS
        cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL_DP)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij= sqrt(cPow2_ij)

        ! Diagonal matrix of eigenvectors
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)
        l5 = abs(vel_ij)

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
        
          R_ij(1,4) =  0.0
          R_ij(2,4) =  l4*a(2)
          R_ij(3,4) = -l4*a(1)
          R_ij(4,4) =  0.0
          R_ij(5,4) =  l4*(u_ij*a(2)-v_ij*a(1))

          R_ij(1,5) =  0.0
          R_ij(2,5) = -l5*a(3)
          R_ij(3,5) =  0.0
          R_ij(4,5) =  l5*a(1)
          R_ij(5,5) =  l5*(w_ij*a(1)-u_ij*a(3))

          ! Matrix of left eigenvectors
          L_ij(1,1) =  0.5*((GAMMA-1.0)*q_ij+c_ij*vel_ij)/cPow2_ij
          L_ij(2,1) =  (cPow2_ij-(GAMMA-1.0)*q_ij)/cPow2_ij
          L_ij(3,1) =  0.5*((GAMMA-1.0)*q_ij-c_ij*vel_ij)/cPow2_ij
          L_ij(4,1) =  (v_ij-vel_ij*a(2))/a(1)
          L_ij(5,1) =  (vel_ij*a(3)-w_ij)/a(1)
          
          L_ij(1,2) =  0.5*(-(GAMMA-1.0)*u_ij-c_ij*a(1))/cPow2_ij
          L_ij(2,2) =  (GAMMA-1.0)*u_ij/cPow2_ij
          L_ij(3,2) =  0.5*(-(GAMMA-1.0)*u_ij+c_ij*a(1))/cPow2_ij
          L_ij(4,2) =  a(2)
          L_ij(5,2) = -a(3)
          
          L_ij(1,3) =  0.5*(-(GAMMA-1.0)*v_ij-c_ij*a(2))/cPow2_ij
          L_ij(2,3) =  (GAMMA-1.0)*v_ij/cPow2_ij
          L_ij(3,3) =  0.5*(-(GAMMA-1.0)*v_ij+c_ij*a(2))/cPow2_ij
          L_ij(4,3) =  (a(2)*a(2)-1.0)/a(1)
          L_ij(5,3) = -a(2)*a(3)/a(1)
          
          L_ij(1,4) =  0.5*(-(GAMMA-1.0)*w_ij-c_ij*a(3))/cPow2_ij
          L_ij(2,4) =  (GAMMA-1.0)*w_ij/cPow2_ij
          L_ij(3,4) =  0.5*(-(GAMMA-1.0)*w_ij+c_ij*a(3))/cPow2_ij
          L_ij(4,4) =  a(2)*a(3)/a(1)
          L_ij(5,4) = -(1-a(3)*a(3))/a(1)

          L_ij(1,5) =  0.5*(GAMMA-1.0)/cPow2_ij
          L_ij(2,5) = -(GAMMA-1.0)/cPow2_ij
          L_ij(3,5) =  0.5*(GAMMA-1.0)/cPow2_ij
          L_ij(4,5) =  0.0
          L_ij(5,5) =  0.0

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
        
          R_ij(1,4) =  0.0
          R_ij(2,4) =  l4*a(2)
          R_ij(3,4) = -l4*a(1)
          R_ij(4,4) =  0.0
          R_ij(5,4) =  l4*(u_ij*a(2)-v_ij*a(1))

          R_ij(1,5) =  0.0
          R_ij(2,5) =  0.0
          R_ij(3,5) =  l5*a(3)
          R_ij(4,5) = -l5*a(2)
          R_ij(5,5) =  l5*(v_ij*a(3)-w_ij*a(2))

          ! Matrix of left eigenvectors
          L_ij(1,1) =  0.5*((GAMMA-1.0)*q_ij+c_ij*vel_ij)/cPow2_ij
          L_ij(2,1) =  (cPow2_ij-(GAMMA-1.0)*q_ij)/cPow2_ij
          L_ij(3,1) =  0.5*((GAMMA-1.0)*q_ij-c_ij*vel_ij)/cPow2_ij
          L_ij(4,1) =  (vel_ij*a(1)-u_ij)/a(2)
          L_ij(5,1) =  (w_ij-vel_ij*a(3))/a(2)
          
          L_ij(1,2) =  0.5*(-(GAMMA-1.0)*u_ij-c_ij*a(1))/cPow2_ij
          L_ij(2,2) =  (GAMMA-1.0)*u_ij/cPow2_ij
          L_ij(3,2) =  0.5*(-(GAMMA-1.0)*u_ij+c_ij*a(1))/cPow2_ij
          L_ij(4,2) =  (1.0-a(1)*a(1))/a(2)
          L_ij(5,2) =  a(1)*a(3)/a(2)
          
          L_ij(1,3) =  0.5*(-(GAMMA-1.0)*v_ij-c_ij*a(2))/cPow2_ij
          L_ij(2,3) =  (GAMMA-1.0)*v_ij/cPow2_ij
          L_ij(3,3) =  0.5*(-(GAMMA-1.0)*v_ij+c_ij*a(2))/cPow2_ij
          L_ij(4,3) = -a(1)
          L_ij(5,3) =  a(3)
          
          L_ij(1,4) =  0.5*(-(GAMMA-1.0)*w_ij-c_ij*a(3))/cPow2_ij
          L_ij(2,4) =  (GAMMA-1.0)*w_ij/cPow2_ij
          L_ij(3,4) =  0.5*(-(GAMMA-1.0)*w_ij+c_ij*a(3))/cPow2_ij
          L_ij(4,4) = -a(1)*a(3)/a(2)
          L_ij(5,4) =  (a(3)*a(3)-1.0)/a(2)

          L_ij(1,5) =  0.5*(GAMMA-1.0)/cPow2_ij
          L_ij(2,5) = -(GAMMA-1.0)/cPow2_ij
          L_ij(3,5) =  0.5*(GAMMA-1.0)/cPow2_ij
          L_ij(4,5) =  0.0
          L_ij(5,5) =  0.0
          
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
        
          R_ij(1,4) =  0.0
          R_ij(2,4) = -l4*a(3)
          R_ij(3,4) =  0.0
          R_ij(4,4) =  l4*a(1)
          R_ij(5,4) =  l4*(w_ij*a(1)-u_ij*a(3))

          R_ij(1,5) =  0.0
          R_ij(2,5) =  0.0
          R_ij(3,5) =  l5*a(3)
          R_ij(4,5) = -l5*a(2)
          R_ij(5,5) =  l5*(v_ij*a(3)-w_ij*a(2))

          ! Matrix of left eigenvectors
          L_ij(1,1) =  0.5*((GAMMA-1.0)*q_ij+c_ij*vel_ij)/cPow2_ij
          L_ij(2,1) =  (cPow2_ij-(GAMMA-1.0)*q_ij)/cPow2_ij
          L_ij(3,1) =  0.5*((GAMMA-1.0)*q_ij-c_ij*vel_ij)/cPow2_ij
          L_ij(4,1) =  (u_ij-vel_ij*a(1))/a(3)
          L_ij(5,1) =  (vel_ij*a(2)-v_ij)/a(3)
          
          L_ij(1,2) =  0.5*(-(GAMMA-1.0)*u_ij-c_ij*a(1))/cPow2_ij
          L_ij(2,2) =  (GAMMA-1.0)*u_ij/cPow2_ij
          L_ij(3,2) =  0.5*(-(GAMMA-1.0)*u_ij+c_ij*a(1))/cPow2_ij
          L_ij(4,2) =  (a(1)*a(1)-1.0)/a(3)
          L_ij(5,2) = -a(1)*a(2)/a(3)
          
          L_ij(1,3) =  0.5*(-(GAMMA-1.0)*v_ij-c_ij*a(2))/cPow2_ij
          L_ij(2,3) =  (GAMMA-1.0)*v_ij/cPow2_ij
          L_ij(3,3) =  0.5*(-(GAMMA-1.0)*v_ij+c_ij*a(2))/cPow2_ij
          L_ij(4,3) =  a(1)*a(2)/a(3)
          L_ij(5,3) =  (1.0-a(2)*a(2))/a(3)
          
          L_ij(1,4) =  0.5*(-(GAMMA-1.0)*w_ij-c_ij*a(3))/cPow2_ij
          L_ij(2,4) =  (GAMMA-1.0)*w_ij/cPow2_ij
          L_ij(3,4) =  0.5*(-(GAMMA-1.0)*w_ij+c_ij*a(3))/cPow2_ij
          L_ij(4,4) =  a(1)
          L_ij(5,4) = -a(2)

          L_ij(1,5) =  0.5*(GAMMA-1.0)/cPow2_ij
          L_ij(2,5) = -(GAMMA-1.0)/cPow2_ij
          L_ij(3,5) =  0.5*(GAMMA-1.0)/cPow2_ij
          L_ij(4,5) =  0.0
          L_ij(5,5) =  0.0

        end select

        ! Include scaling parameter
        anorm = dscale*anorm
        
        ! Compute tensorial dissipation D_ij = diag(R_ij*|Lbd_ij|*L_ij)*I
        DcoefficientsAtEdge(1,1,idx) = anorm*( R_ij(1,1)*L_ij(1,1)+&
                                               R_ij(1,2)*L_ij(2,1)+&
                                               R_ij(1,3)*L_ij(3,1)+&
                                               R_ij(1,4)*L_ij(4,1)+&
                                               R_ij(1,5)*L_ij(5,1)  )
        DcoefficientsAtEdge(2,1,idx) = anorm*( R_ij(2,1)*L_ij(1,2)+&
                                               R_ij(2,2)*L_ij(2,2)+&
                                               R_ij(2,3)*L_ij(3,2)+&
                                               R_ij(2,4)*L_ij(4,2)+&
                                               R_ij(2,5)*L_ij(5,2)  )
        DcoefficientsAtEdge(3,1,idx) = anorm*( R_ij(3,1)*L_ij(1,3)+&
                                               R_ij(3,2)*L_ij(2,3)+&
                                               R_ij(3,3)*L_ij(3,3)+&
                                               R_ij(3,4)*L_ij(4,3)+&
                                               R_ij(3,5)*L_ij(5,3)  )
        DcoefficientsAtEdge(4,1,idx) = anorm*( R_ij(4,1)*L_ij(1,4)+&
                                               R_ij(4,2)*L_ij(2,4)+&
                                               R_ij(4,3)*L_ij(3,4)+&
                                               R_ij(4,4)*L_ij(4,4)+&
                                               R_ij(4,5)*L_ij(5,4)  )
        DcoefficientsAtEdge(5,1,idx) = anorm*( R_ij(5,1)*L_ij(1,5)+&
                                               R_ij(5,2)*L_ij(2,5)+&
                                               R_ij(5,3)*L_ij(3,5)+&
                                               R_ij(5,4)*L_ij(4,5)+&
                                               R_ij(5,5)*L_ij(5,5)  )
      else

        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0
        
      end if
    end do

  end subroutine hydro_calcMatRoeDissMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRoeDiss3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D and applies
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
    real(DP), dimension(NVAR3D,NVAR3D) :: R_ij,L_ij
    real(DP), dimension(NVAR3D) :: a
    real(DP) :: Ei,Ej,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,anorm,aux,cPow2_ij,c_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    real(DP) :: l1,l2,l3,l4,l5
    integer :: idx,i,j,k


    do idx = 1, nedges
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      Ei = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      Ej = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),ui,vi,wi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),ui,vi,wi,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient $0.5*(C_{ij}-C_{ji})$ and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      
      if (anorm .gt. SYS_EPSREAL_DP) then

        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)+\
               PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)+\
               PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx),aux)

        ! Compute auxiliary values
        vel_ij = u_ij*a(1)+v_ij*a(2)+w_ij*a(3)
        q_ij   = 0.5*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)
        
        ! Compute speed of sound
#ifdef THERMALLY_IDEAL_GAS
        cPow2_ij = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL_DP)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij= sqrt(cPow2_ij)

        ! Diagonal matrix of eigenvectors
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)
        l5 = abs(vel_ij)

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
        
          R_ij(1,4) =  0.0
          R_ij(2,4) =  l4*a(2)
          R_ij(3,4) = -l4*a(1)
          R_ij(4,4) =  0.0
          R_ij(5,4) =  l4*(u_ij*a(2)-v_ij*a(1))

          R_ij(1,5) =  0.0
          R_ij(2,5) = -l5*a(3)
          R_ij(3,5) =  0.0
          R_ij(4,5) =  l5*a(1)
          R_ij(5,5) =  l5*(w_ij*a(1)-u_ij*a(3))

          ! Matrix of left eigenvectors
          L_ij(1,1) =  0.5*((GAMMA-1.0)*q_ij+c_ij*vel_ij)/cPow2_ij
          L_ij(2,1) =  (cPow2_ij-(GAMMA-1.0)*q_ij)/cPow2_ij
          L_ij(3,1) =  0.5*((GAMMA-1.0)*q_ij-c_ij*vel_ij)/cPow2_ij
          L_ij(4,1) =  (v_ij-vel_ij*a(2))/a(1)
          L_ij(5,1) =  (vel_ij*a(3)-w_ij)/a(1)
          
          L_ij(1,2) =  0.5*(-(GAMMA-1.0)*u_ij-c_ij*a(1))/cPow2_ij
          L_ij(2,2) =  (GAMMA-1.0)*u_ij/cPow2_ij
          L_ij(3,2) =  0.5*(-(GAMMA-1.0)*u_ij+c_ij*a(1))/cPow2_ij
          L_ij(4,2) =  a(2)
          L_ij(5,2) = -a(3)
          
          L_ij(1,3) =  0.5*(-(GAMMA-1.0)*v_ij-c_ij*a(2))/cPow2_ij
          L_ij(2,3) =  (GAMMA-1.0)*v_ij/cPow2_ij
          L_ij(3,3) =  0.5*(-(GAMMA-1.0)*v_ij+c_ij*a(2))/cPow2_ij
          L_ij(4,3) =  (a(2)*a(2)-1.0)/a(1)
          L_ij(5,3) = -a(2)*a(3)/a(1)
          
          L_ij(1,4) =  0.5*(-(GAMMA-1.0)*w_ij-c_ij*a(3))/cPow2_ij
          L_ij(2,4) =  (GAMMA-1.0)*w_ij/cPow2_ij
          L_ij(3,4) =  0.5*(-(GAMMA-1.0)*w_ij+c_ij*a(3))/cPow2_ij
          L_ij(4,4) =  a(2)*a(3)/a(1)
          L_ij(5,4) = -(1-a(3)*a(3))/a(1)

          L_ij(1,5) =  0.5*(GAMMA-1.0)/cPow2_ij
          L_ij(2,5) = -(GAMMA-1.0)/cPow2_ij
          L_ij(3,5) =  0.5*(GAMMA-1.0)/cPow2_ij
          L_ij(4,5) =  0.0
          L_ij(5,5) =  0.0

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
        
          R_ij(1,4) =  0.0
          R_ij(2,4) =  l4*a(2)
          R_ij(3,4) = -l4*a(1)
          R_ij(4,4) =  0.0
          R_ij(5,4) =  l4*(u_ij*a(2)-v_ij*a(1))

          R_ij(1,5) =  0.0
          R_ij(2,5) =  0.0
          R_ij(3,5) =  l5*a(3)
          R_ij(4,5) = -l5*a(2)
          R_ij(5,5) =  l5*(v_ij*a(3)-w_ij*a(2))

          ! Matrix of left eigenvectors
          L_ij(1,1) =  0.5*((GAMMA-1.0)*q_ij+c_ij*vel_ij)/cPow2_ij
          L_ij(2,1) =  (cPow2_ij-(GAMMA-1.0)*q_ij)/cPow2_ij
          L_ij(3,1) =  0.5*((GAMMA-1.0)*q_ij-c_ij*vel_ij)/cPow2_ij
          L_ij(4,1) =  (vel_ij*a(1)-u_ij)/a(2)
          L_ij(5,1) =  (w_ij-vel_ij*a(3))/a(2)
          
          L_ij(1,2) =  0.5*(-(GAMMA-1.0)*u_ij-c_ij*a(1))/cPow2_ij
          L_ij(2,2) =  (GAMMA-1.0)*u_ij/cPow2_ij
          L_ij(3,2) =  0.5*(-(GAMMA-1.0)*u_ij+c_ij*a(1))/cPow2_ij
          L_ij(4,2) =  (1.0-a(1)*a(1))/a(2)
          L_ij(5,2) =  a(1)*a(3)/a(2)
          
          L_ij(1,3) =  0.5*(-(GAMMA-1.0)*v_ij-c_ij*a(2))/cPow2_ij
          L_ij(2,3) =  (GAMMA-1.0)*v_ij/cPow2_ij
          L_ij(3,3) =  0.5*(-(GAMMA-1.0)*v_ij+c_ij*a(2))/cPow2_ij
          L_ij(4,3) = -a(1)
          L_ij(5,3) =  a(3)
          
          L_ij(1,4) =  0.5*(-(GAMMA-1.0)*w_ij-c_ij*a(3))/cPow2_ij
          L_ij(2,4) =  (GAMMA-1.0)*w_ij/cPow2_ij
          L_ij(3,4) =  0.5*(-(GAMMA-1.0)*w_ij+c_ij*a(3))/cPow2_ij
          L_ij(4,4) = -a(1)*a(3)/a(2)
          L_ij(5,4) =  (a(3)*a(3)-1.0)/a(2)

          L_ij(1,5) =  0.5*(GAMMA-1.0)/cPow2_ij
          L_ij(2,5) = -(GAMMA-1.0)/cPow2_ij
          L_ij(3,5) =  0.5*(GAMMA-1.0)/cPow2_ij
          L_ij(4,5) =  0.0
          L_ij(5,5) =  0.0
          
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
        
          R_ij(1,4) =  0.0
          R_ij(2,4) = -l4*a(3)
          R_ij(3,4) =  0.0
          R_ij(4,4) =  l4*a(1)
          R_ij(5,4) =  l4*(w_ij*a(1)-u_ij*a(3))

          R_ij(1,5) =  0.0
          R_ij(2,5) =  0.0
          R_ij(3,5) =  l5*a(3)
          R_ij(4,5) = -l5*a(2)
          R_ij(5,5) =  l5*(v_ij*a(3)-w_ij*a(2))

          ! Matrix of left eigenvectors
          L_ij(1,1) =  0.5*((GAMMA-1.0)*q_ij+c_ij*vel_ij)/cPow2_ij
          L_ij(2,1) =  (cPow2_ij-(GAMMA-1.0)*q_ij)/cPow2_ij
          L_ij(3,1) =  0.5*((GAMMA-1.0)*q_ij-c_ij*vel_ij)/cPow2_ij
          L_ij(4,1) =  (u_ij-vel_ij*a(1))/a(3)
          L_ij(5,1) =  (vel_ij*a(2)-v_ij)/a(3)
          
          L_ij(1,2) =  0.5*(-(GAMMA-1.0)*u_ij-c_ij*a(1))/cPow2_ij
          L_ij(2,2) =  (GAMMA-1.0)*u_ij/cPow2_ij
          L_ij(3,2) =  0.5*(-(GAMMA-1.0)*u_ij+c_ij*a(1))/cPow2_ij
          L_ij(4,2) =  (a(1)*a(1)-1.0)/a(3)
          L_ij(5,2) = -a(1)*a(2)/a(3)
          
          L_ij(1,3) =  0.5*(-(GAMMA-1.0)*v_ij-c_ij*a(2))/cPow2_ij
          L_ij(2,3) =  (GAMMA-1.0)*v_ij/cPow2_ij
          L_ij(3,3) =  0.5*(-(GAMMA-1.0)*v_ij+c_ij*a(2))/cPow2_ij
          L_ij(4,3) =  a(1)*a(2)/a(3)
          L_ij(5,3) =  (1.0-a(2)*a(2))/a(3)
          
          L_ij(1,4) =  0.5*(-(GAMMA-1.0)*w_ij-c_ij*a(3))/cPow2_ij
          L_ij(2,4) =  (GAMMA-1.0)*w_ij/cPow2_ij
          L_ij(3,4) =  0.5*(-(GAMMA-1.0)*w_ij+c_ij*a(3))/cPow2_ij
          L_ij(4,4) =  a(1)
          L_ij(5,4) = -a(2)

          L_ij(1,5) =  0.5*(GAMMA-1.0)/cPow2_ij
          L_ij(2,5) = -(GAMMA-1.0)/cPow2_ij
          L_ij(3,5) =  0.5*(GAMMA-1.0)/cPow2_ij
          L_ij(4,5) =  0.0
          L_ij(5,5) =  0.0

        end select

        ! Include scaling parameter
        anorm = dscale*anorm

        ! Compute tensorial dissipation D_ij = R_ij*|Lbd_ij|*L_ij
        do i = 1, NVAR3D
          do j = 1, NVAR3D
            aux = 0.0
            do k = 1, NVAR3D
              aux = aux + R_ij(i,k)*L_ij(k,j)
            end do
            DcoefficientsAtEdge(NVAR3D*(j-1)+i,1,idx) = anorm*aux
          end do
        end do
        
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0
        
      end if
    end do

  end subroutine hydro_calcMatRoeDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDissMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 3D and applies the scalar artificial viscosities of
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
    real(DP) :: Ei,Ej,ci,cj,ui,uj,vi,vj,wi,wj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      Ei = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      Ej = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP      
      ! Compute Galerkin coefficient $K_ij = diag(A_j)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),uj,vj,wj)
      
      ! Compute Galerkin coefficient $K_ji = diag(A_i)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),ui,vi,wi)
#else
      ! Compute Galerkin coefficient $K_ij = -diag(A_j)*C_{ij}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),uj,vj,wj)
      
      ! Compute Galerkin coefficient $K_ji = -diag(A_i)*C_{ji}$
      MATRIXDIAG_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),ui,vi,wi)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate scalar artificial dissipation of Rusanov-type
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*(ui*ui+vi*vi+wi*wi)), SYS_EPSREAL_DP))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*(uj*uj+vj*vj+wj*wj)), SYS_EPSREAL_DP))
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Compute dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = dscale *&
          max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj+&
                   DmatrixCoeffsAtEdge(2,1,idx)*vj+&
                   DmatrixCoeffsAtEdge(3,1,idx)*wj) +&
                   sqrt(DmatrixCoeffsAtEdge(1,1,idx)**2+&
                        DmatrixCoeffsAtEdge(2,1,idx)**2+&
                        DmatrixCoeffsAtEdge(3,1,idx)**2)*cj,&
               abs(DmatrixCoeffsAtEdge(1,2,idx)*ui+&
                   DmatrixCoeffsAtEdge(2,2,idx)*vi+&
                   DmatrixCoeffsAtEdge(3,2,idx)*wi) +&
                   sqrt(DmatrixCoeffsAtEdge(1,2,idx)**2+&
                        DmatrixCoeffsAtEdge(2,2,idx)**2+&
                        DmatrixCoeffsAtEdge(3,2,idx)**2)*ci )
    end do

  end subroutine hydro_calcMatRusDissMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcMatRusDiss3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D applies
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
    real(DP) :: Ei,Ej,aux,ci,cj,ui,uj,vi,vj,wi,wj
    integer :: idx


    do idx = 1, nedges
      
      ! Compute auxiliary variables
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      Ei = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      Ej = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0

#ifdef HYDRO_USE_IBP
      ! Compute Galerkin coefficient $K_ij = A_j*C_{ji}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = A_i*C_{ij}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),ui,vi,wi,Ei)
#else
      ! Compute Galerkin coefficient $K_ij = -A_j*C_{ij}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,2,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,1,idx),DmatrixCoeffsAtEdge(2,1,idx),\
        DmatrixCoeffsAtEdge(3,1,idx),uj,vj,wj,Ej)
      
      ! Compute Galerkin coefficient $K_ji = -A_i*C_{ji}$
      MATRIX_HYDRO_2T_3D(DcoefficientsAtEdge,3,idx,\
        -dscale,DmatrixCoeffsAtEdge(1,2,idx),DmatrixCoeffsAtEdge(2,2,idx),\
        DmatrixCoeffsAtEdge(3,2,idx),ui,vi,wi,Ei)
#endif

      !---------------------------------------------------------------------------
      ! Evaluate scalar artificial dissipation of Rusanov-type
      !---------------------------------------------------------------------------

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*(ui*ui+vi*vi+wi*wi)), SYS_EPSREAL_DP))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*(uj*uj+vj*vj+wj*wj)), SYS_EPSREAL_DP))
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Compute dissipation tensor
      aux = dscale *&
          max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj+&
                   DmatrixCoeffsAtEdge(2,1,idx)*vj+&
                   DmatrixCoeffsAtEdge(3,1,idx)*wj) +&
                   sqrt(DmatrixCoeffsAtEdge(1,1,idx)**2+&
                        DmatrixCoeffsAtEdge(2,1,idx)**2+&
                        DmatrixCoeffsAtEdge(3,1,idx)**2)*cj,&
               abs(DmatrixCoeffsAtEdge(1,2,idx)*ui+&
                   DmatrixCoeffsAtEdge(2,2,idx)*vi+&
                   DmatrixCoeffsAtEdge(3,2,idx)*wi) +&
                   sqrt(DmatrixCoeffsAtEdge(1,2,idx)**2+&
                        DmatrixCoeffsAtEdge(2,2,idx)**2+&
                        DmatrixCoeffsAtEdge(3,2,idx)**2)*ci )

      DcoefficientsAtEdge( :,1,idx) = 0.0
      DcoefficientsAtEdge( 1,1,idx) = aux
      DcoefficientsAtEdge( 6,1,idx) = aux
      DcoefficientsAtEdge(11,1,idx) = aux
      DcoefficientsAtEdge(16,1,idx) = aux
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

  pure subroutine hydro_calcFluxFCTScDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)


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
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: pi,pj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,anorm,aux,c_ij,d_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    integer :: idx


    do idx = 1, nedges

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute pressures
      pi = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
      pj = PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)

      ! Compute skew-symmetric coefficient
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      v_ij = ROE_MEAN_VALUE(vi,vj,aux)
      w_ij = ROE_MEAN_VALUE(wi,wj,aux)
      H_ij = ROE_MEAN_VALUE(\
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)+pi)/\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)+pj)/\
             DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx),aux)

      ! Compute auxiliary variables
      vel_ij = u_ij*a(1) + v_ij*a(2) + w_ij*a(3)
      q_ij   = 0.5*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      c_ij = sqrt(max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL_DP))
#else
#error "Speed of sound must be implemented!"
#endif

      ! Compute scalar dissipation
      d_ij = abs(vel_ij) + anorm*c_ij

      ! Compute antidiffusive flux
       DfluxesAtEdge(:,idx) = dscale*d_ij*(DdataAtEdge(:,1,idx)-DdataAtEdge(:,2,idx))
     end do

   end subroutine hydro_calcFluxFCTScDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRoeDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)

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
    real(DP), dimension(NVAR3D) :: Diff
    real(DP), dimension(NDIM3D) :: a
    real(DP) :: ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,c2_ij,c_ij,q_ij,u_ij,v_ij,vel_ij,w_ij
    real(DP) :: anorm,aux,aux1,aux2
    real(DP) :: l1,l2,l3,l4,l5,w1,w2,w3,w4,w5
    integer :: idx


    do idx = 1, nedges

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute the skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(:,1,idx)-DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      if (anorm .gt. SYS_EPSREAL_DP) then
        
        ! Normalize the skew-symmetric coefficient
        a = a/anorm
        
        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)+\
               PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)+\
               PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx))/\
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx),aux)

        ! Compute auxiliary variables
        vel_ij = u_ij*a(1) + v_ij*a(2) + w_ij*a(3)
        q_ij   = 0.5*(u_ij*u_ij+v_ij*v_ij+w_ij*w_ij)

        ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
        c2_ij  = max((GAMMA-1.0)*(H_ij-q_ij), SYS_EPSREAL_DP)
#else
#error "Speed of sound must be implemented!"
#endif
        c_ij   = sqrt(c2_ij)

        ! Compute eigenvalues
        l1 = abs(vel_ij-c_ij)
        l2 = abs(vel_ij)
        l3 = abs(vel_ij+c_ij)
        l4 = abs(vel_ij)
        l5 = abs(vel_ij)

        ! Compute solution difference U_i-U_j
        Diff = DdataAtEdge(:,1,idx)-DdataAtEdge(:,2,idx)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = (GAMMA-1.0)/2.0/c2_ij*(q_ij*Diff(1)&
                                     -u_ij*Diff(2)&
                                     -v_ij*Diff(3)&
                                     -w_ij*Diff(4)&
                                          +Diff(5))
        aux2 = (vel_ij*Diff(1)&
                 -a(1)*Diff(2)&
                 -a(2)*Diff(3)&
                 -a(3)*Diff(4))/2.0/c_ij

        ! Get the dimension with largest coefficient
        select case(maxloc(a,1))
        case(1)
          ! Compute characteristic variables multiplied by the corresponding eigenvalue
          w1 = l1 * (aux1 + aux2)
          w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/c2_ij)*Diff(1)&
                              +(GAMMA-1.0)*( u_ij*Diff(2)&
                                            +v_ij*Diff(3)&
                                            +w_ij*Diff(4)&
                                                 -Diff(5))/c2_ij)
          w3 = l3 * (aux1 - aux2)
          w4 = l4 * ( (v_ij-vel_ij*a(2))/a(1)*Diff(1)&
                                        +a(2)*Diff(2)&
                        +(a(2)*a(2)-1.0)/a(1)*Diff(3)&
                              +a(2)*a(3)/a(1)*Diff(4))
          w5 = l5 * ( (vel_ij*a(3)-w_ij)/a(1)*Diff(1)&
                                        -a(3)*Diff(2)&
                              -a(2)*a(3)/a(1)*Diff(3)&
                        +(1.0-a(3)*a(3))/a(1)*Diff(4))

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
          w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/c2_ij)*Diff(1)&
                               +(GAMMA-1.0)*(u_ij*Diff(2)&
                                            +v_ij*Diff(3)&
                                            +w_ij*Diff(4)&
                                                 -Diff(5))/c2_ij)
          w3 = l3 * (aux1 - aux2)
          w4 = l4 * ( (vel_ij*a(1)-u_ij)/a(2)*Diff(1)&
                        +(1.0-a(1)*a(1))/a(2)*Diff(2)&
                                        -a(1)*Diff(3)&
                              -a(1)*a(3)/a(2)*Diff(4))
          w5 = l5 * ( (w_ij-vel_ij*a(3))/a(2)*Diff(1)&
                              +a(1)*a(3)/a(2)*Diff(2)&
                                        +a(3)*Diff(3)&
                        +(a(3)*a(3)-1.0)/a(2)*Diff(4))

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
          w2 = l2 * ((1.0-(GAMMA-1.0)*q_ij/c2_ij)*Diff(1)&
                               +(GAMMA-1.0)*(u_ij*Diff(2)&
                                            +v_ij*Diff(3)&
                                            +w_ij*Diff(4)&
                                                 -Diff(5))/c2_ij)
          w3 = l3 * (aux1 - aux2)
          w4 = l4 * ( (u_ij-vel_ij*a(1))/a(3)*Diff(1)&
                        +(a(1)*a(1)-1.0)/a(3)*Diff(2)&
                              +a(1)*a(2)/a(3)*Diff(3)&
                                        +a(1)*Diff(4) )
          w5 = l5 * ( (vel_ij*a(2)-v_ij)/a(3)*Diff(1)&
                              -a(1)*a(2)/a(3)*Diff(2)&
                        +(1.0-a(2)*a(2))/a(3)*Diff(3)&
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
        DfluxesAtEdge(:,idx) = dscale*Diff
      else
        ! Cancel antidiffusive flux
        DfluxesAtEdge(:,idx) = 0
      end if
    end do

  end subroutine hydro_calcFluxFCTRoeDiss3d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine hydro_calcFluxFCTRusDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DfluxesAtEdge, rcollection)

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
    real(DP), dimension(NVAR3D) :: Diff
    real(DP) :: Ei,Ej,ci,cj,ui,uj,vi,vj,wi,wj
    real(DP) :: d_ij
    integer :: idx

    
    do idx = 1, nedges

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute specific energies
      Ei = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      Ej = SPECIFIC_TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      ci = sqrt(max((GAMMA-1.0)*GAMMA*(Ei-0.5*(ui*ui+vi*vi+wi*wi)), SYS_EPSREAL_DP))
      cj = sqrt(max((GAMMA-1.0)*GAMMA*(Ej-0.5*(uj*uj+vj*vj+wj*wj)), SYS_EPSREAL_DP))
#else
#error "Speed of sound must be implemented!"
#endif
      
#ifdef HYDRO_USE_IBP
      ! Compute scalar dissipation based on the skew-symmetric part
      ! which does not include the symmetric boundary contribution
      d_ij = max( abs(0.5_DP*(DmatrixCoeffsAtEdge(1,1,idx)-&
                              DmatrixCoeffsAtEdge(1,2,idx))*uj+&
                      0.5_DP*(DmatrixCoeffsAtEdge(2,1,idx)-&
                              DmatrixCoeffsAtEdge(2,2,idx))*vj+&
                      0.5_DP*(DmatrixCoeffsAtEdge(3,1,idx)-&
                              DmatrixCoeffsAtEdge(3,2,idx))*wj)+&
                 0.5_DP*sqrt((DmatrixCoeffsAtEdge(1,1,idx)-&
                              DmatrixCoeffsAtEdge(1,2,idx))**2+&
                             (DmatrixCoeffsAtEdge(2,1,idx)-&
                              DmatrixCoeffsAtEdge(2,2,idx))**2+&
                             (DmatrixCoeffsAtEdge(3,1,idx)-&
                              DmatrixCoeffsAtEdge(3,2,idx))**2)*cj,&
                  abs(0.5_DP*(DmatrixCoeffsAtEdge(1,2,idx)-&
                              DmatrixCoeffsAtEdge(1,1,idx))*ui+&
                      0.5_DP*(DmatrixCoeffsAtEdge(2,2,idx)-&
                              DmatrixCoeffsAtEdge(2,1,idx))*vi+&
                      0.5_DP*(DmatrixCoeffsAtEdge(3,2,idx)-&
                              DmatrixCoeffsAtEdge(3,1,idx))*wi)+&
                 0.5_DP*sqrt((DmatrixCoeffsAtEdge(1,2,idx)-&
                              DmatrixCoeffsAtEdge(1,1,idx))**2+&
                             (DmatrixCoeffsAtEdge(2,2,idx)-&
                              DmatrixCoeffsAtEdge(2,1,idx))**2+&
                             (DmatrixCoeffsAtEdge(3,2,idx)-&
                              DmatrixCoeffsAtEdge(3,1,idx))**2)*ci )
#else
      ! Compute scalar dissipation
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj+&
                      DmatrixCoeffsAtEdge(2,1,idx)*vj+&
                      DmatrixCoeffsAtEdge(3,1,idx)*wj)+&
                 sqrt(DmatrixCoeffsAtEdge(1,1,idx)**2+&
                      DmatrixCoeffsAtEdge(2,1,idx)**2+&
                      DmatrixCoeffsAtEdge(3,1,idx)**2)*cj,&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui+&
                      DmatrixCoeffsAtEdge(2,2,idx)*vi+&
                      DmatrixCoeffsAtEdge(3,2,idx)*wi)+&
                 sqrt(DmatrixCoeffsAtEdge(1,2,idx)**2+&
                      DmatrixCoeffsAtEdge(2,2,idx)**2+&
                      DmatrixCoeffsAtEdge(3,2,idx)**2)*ci )
#endif

      ! Compute antidiffusive flux
      DfluxesAtEdge(:,idx) = dscale* d_ij*(DdataAtEdge(:,1,idx)-DdataAtEdge(:,2,idx))
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
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
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
      DtransformedDataAtEdge(1,idx) =&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
    end do

  end subroutine hydro_trafoDiffDensity3d_sim

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
      DtransformedFluxesAtEdge(1,1,idx) =&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
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
      DtransformedDataAtEdge(1,idx) =&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
    end do

  end subroutine hydro_trafoDiffEnergy3d_sim

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
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(1,1,idx) = (GAMMA-1.0)*&
          (0.5*(ui*ui+vi*vi)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vi*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wi*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))
      DtransformedFluxesAtEdge(1,2,idx) =-(GAMMA-1.0)*&
          (0.5*(uj*uj+vj*vj)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vj*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wj*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
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
      DtransformedDataAtEdge(1,idx) =&
          PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
    end do

  end subroutine hydro_trafoDiffPressure3d_sim

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
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      
      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(1,1,idx) =&
          (X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -(X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Transformed velocity fluxes in y-direction
      DtransformedFluxesAtEdge(2,1,idx) =&
          (Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -(Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Transformed velocity fluxes in z-direction
      DtransformedFluxesAtEdge(3,1,idx) =&
          (Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      DtransformedFluxesAtEdge(3,2,idx) =&
         -(Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
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
      DtransformedDataAtEdge(1,idx) =&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed velocity difference in y-direction
      DtransformedDataAtEdge(2,idx) =&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed velocity difference in z-direction
      DtransformedDataAtEdge(3,idx) =&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
    end do
    
  end subroutine hydro_trafoDiffVelocity3d_sim

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
      DtransformedFluxesAtEdge(1,1,idx) =&
          X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)

      ! Transformed momentum fluxes in y-direction
      DtransformedFluxesAtEdge(2,1,idx) =&
          Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)

      ! Transformed momentum fluxes in z-direction
      DtransformedFluxesAtEdge(3,1,idx) =&
          Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(3,2,idx) =&
         -Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
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
      DtransformedDataAtEdge(1,idx) =&
          X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          X_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed momentum difference in y-direction
      DtransformedDataAtEdge(2,idx) =&
          Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed momentum difference in z-direction
      DtransformedDataAtEdge(3,idx) =&
          Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
    end do
    
  end subroutine hydro_trafoDiffMomentum3d_sim

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
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)

      ! Transformed total energy fluxes
      DtransformedFluxesAtEdge(2,1,idx) =&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
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
      DtransformedDataAtEdge(1,idx) =&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed total energy difference
      DtransformedDataAtEdge(2,idx) =&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
    end do

  end subroutine hydro_trafoDiffDenEng3d_sim

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
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      
      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(2,1,idx) = (GAMMA-1.0)*&
          (0.5*(ui*ui+vi*vi)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vi*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wi*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))
      DtransformedFluxesAtEdge(2,2,idx) =-(GAMMA-1.0)*&
          (0.5*(uj*uj+vj*vj)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vj*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wj*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
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
      DtransformedDataAtEdge(1,idx) =&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(2,idx) =&
          PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
    end do
    
  end subroutine hydro_trafoDiffDenPre3d_sim

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
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)
      
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)

      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(2,1,idx) =&
          (X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -(X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Transformed velocity fluxes in y-direction
      DtransformedFluxesAtEdge(3,1,idx) =&
          (Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      DtransformedFluxesAtEdge(3,2,idx) =&
         -(Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Transformed velocity fluxes in z-direction
      DtransformedFluxesAtEdge(4,1,idx) =&
          (Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      DtransformedFluxesAtEdge(4,2,idx) =&
         -(Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))/&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)

      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(5,1,idx) =(GAMMA-1.0)*&
          (0.5*(ui*ui+vi*vi)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vi*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wi*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))
      DtransformedFluxesAtEdge(5,2,idx) =-(GAMMA-1.0)*&
          (0.5*(uj*uj+vj*vj)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          vj*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)-&
          wj*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx)+&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR3D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
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
      DtransformedDataAtEdge(2,idx) =&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed velocity difference in x-direction
      DtransformedDataAtEdge(2,idx) =&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)
      
      ! Transformed velocity difference in y-direction
      DtransformedDataAtEdge(3,idx) =&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed velocity difference in z-direction
      DtransformedDataAtEdge(4,idx) =&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,2,idx)-&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR3D,1,idx)

      ! Transformed pressure difference
      DtransformedDataAtEdge(5,idx) =&
          PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR_3D(DdataAtEdge,NVAR3D,1,idx)
    end do

  end subroutine hydro_trafoDiffDenPreVel3d_sim

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

  !*****************************************************************************

!<subroutine>

  subroutine hydro_hadaptCallbackScalar3d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 3D. The solution vector is assumed
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
            OU_CLASS_WARNING,OU_MODE_STD,'hydro_hadaptCallbackScalar3d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      do ivar = 1, NVAR3D
        p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR3D+ivar) = &
            0.5_DP*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR3D+ivar)+&
                    p_Dsolution((rcollection%IquickAccess(3)-1)*NVAR3D+ivar))
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      do ivar = 1, NVAR3D
        p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR3D+ivar) = &
            0.25_DP*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR3D+ivar)+&
                     p_Dsolution((rcollection%IquickAccess(3)-1)*NVAR3D+ivar)+&
                     p_Dsolution((rcollection%IquickAccess(4)-1)*NVAR3D+ivar)+&
                     p_Dsolution((rcollection%IquickAccess(5)-1)*NVAR3D+ivar))
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        do ivar = 1, NVAR3D
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR3D+ivar) = &
              p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR3D+ivar)
        end do
      else
        do ivar = 1, NVAR3D
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR3D+ivar) = 0.0_DP
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)

    end select

  end subroutine hydro_hadaptCallbackScalar3d

  !*****************************************************************************

!<subroutine>

  subroutine hydro_hadaptCallbackBlock3d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 3D. The solution vector is assumed
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
      if (rsolution%nblocks .ne. NVAR3D) then
        call output_line('Vector is not in block format!',&
            OU_CLASS_WARNING,OU_MODE_STD,'hydro_hadaptCallbackBlock3d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      neq = rsolution%NEQ/NVAR3D
      do ivar = 1, NVAR3D
        p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = &
            0.5_DP*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
                    p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(3)) )
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      neq = rsolution%NEQ/NVAR3D
      do ivar = 1, NVAR3D
        p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) =&
            0.25_DP*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
                     p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(3))+&
                     p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(4))+&
                     p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(5)) )
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        neq = rsolution%NEQ/NVAR3D
        do ivar = 1, NVAR3D
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = &
              p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))
        end do
      else
        neq = rsolution%NEQ/NVAR3D
        do ivar = 1, NVAR3D
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = 0.0_DP
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)

    end select

  end subroutine hydro_hadaptCallbackBlock3d

end module hydro_callback3d
