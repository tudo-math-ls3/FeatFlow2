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
!#         adopting tensorial artificial viscosities
!#
!# 13.) mhd_calcMatRoeDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 14.) mhd_calcMatRusDissMatD1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov artificial viscosities
!#
!# 15.) mhd_calcMatRusDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov flux artificial viscosities
!#
!# 16.) mhd_calcCharacteristics1d_sim
!#      -> Computes characteristic variables
!#
!# 17.) mhd_calcFluxFCTScDiss1d_sim
!#      -> Computes fluxes for FCT algorithm
!#         adopting scalar artificial viscosities
!#
!# 18.) mhd_calcFluxFCTRoeDiss1d_sim
!#      -> Computes fluxes for FCT algorithm
!#         adopting tensorial artificial viscosities
!#
!# 19.) mhd_calcFluxFCTRusDiss1d_sim
!#      -> Computes fluxes for FCT algorithm
!#         adopting the Rusanov artificial viscosities
!#
!# 20.) mhd_trafoFluxDensity1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 21.) mhd_trafoDiffDensity1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 22.) mhd_trafoFluxEnergy1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 23.) mhd_trafoDiffEnergy1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 24.) mhd_trafoFluxPressure1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 25.) mhd_trafoDiffPressure1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the pressure
!#
!# 26.) mhd_trafoFluxVelocity1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 27.) mhd_trafoDiffVelocity1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 28.) mhd_trafoFluxMomentum1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 29.) mhd_trafoDiffMomentum1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 30.) mhd_trafoFluxDenEng1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 31.) mhd_trafoDiffDenEng1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 32.) mhd_trafoFluxDenPre1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 33.) mhd_trafoDiffDenPre1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 34.) mhd_trafoFluxDenPreVel1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 35.) mhd_trafoDiffDenPreVel1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure 
!#         and the velocity
!#
!# 36.) mhd_trafoDiffMagfield1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the magnetic field
!#
!# 37.) mhd_trafoFluxMagfield1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the magnetic field
!#
!# 38.) mhd_calcBoundaryvalues1d
!#      -> Computes the boundary values for a given node
!#
!# 39.) mhd_hadaptCallbackScalar1d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 1D, whereby the vector is stored in interleave format
!#
!# 40.) mhd_hadaptCallbackBlock1d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 1D, whereby the vector is stored in block format
!#
!# </purpose>
!##############################################################################

module mhd_callback1d

#include "mhd.h"

  use boundarycondaux
  use collection
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
  public :: mhd_calcBoundaryvalues1d
  public :: mhd_hadaptCallbackScalar1d
  public :: mhd_hadaptCallbackBlock1d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGal1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR1D) :: Fi,Fj
#else
    real(DP), dimension(NVAR1D) :: F_ij
#endif
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx


    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj
      
#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_MHD_2T_XDIR_1D(Fi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_1D(Fj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
      
      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fi )
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_1D(F_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
                 
      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * DmatrixCoeffsAtEdge(1,1,idx)*F_ij
      DfluxesAtEdge(:,2,idx) = -dscale * DmatrixCoeffsAtEdge(1,2,idx)*F_ij      
#endif
    end do

  end subroutine mhd_calcFluxGal1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxGalNoBdr1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

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
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    integer :: idx


    do idx = 1, size(DfluxesAtEdge,3)

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !-------------------------------------------------------------------------
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj
      
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_1D(F_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
      
      ! Assemble symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale *&
          0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))*F_ij
      DfluxesAtEdge(:,2,idx) = DfluxesAtEdge(:,1,idx)
    end do

  end subroutine mhd_calcFluxGalNoBdr1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxScDiss1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR1D) :: Fi,Fj
#else
    real(DP), dimension(NVAR1D) :: F_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,X_ij,anorm,aux,cf_ij,d_ij,q_ij,u_ij
    real(DP) :: aPow2_ij,astPow2_ij,bPow2_ij,bxPow2_ij
    integer :: idx


    do idx = 1, size(DfluxesAtEdge,3)
          
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_MHD_2T_XDIR_1D(Fi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_1D(Fj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_1D(F_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar artificial dissipation proportional to the spectral
      ! radius (largest eigenvalue) of the Roe-matrix
      !
      ! There are seven eigenvalues
      !
      !   u-cf, u-ca, u-cs, u, u+cs, u+ca, u+cf,
      !
      ! where u us the x-velocity component and ca, cs and cf are the
      ! velocities of th Alfveen waves, the slow and fast waves.
      !
      ! The largest in magnitude eigenvalue is |u|+cf
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
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+\
              TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))/\
              DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+\
              TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))/\
              DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx), aux)
        
      ! Compute the square of the Roe-averaged speed of the Alfven waves.
      ! Note that left and right states are interchanged!
      bxPow2_ij = (ROE_MEAN_VALUE(\
                    X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                    X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux))**2/\
                  ROE_MEAN_VALUE(\
                   DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),\
                   DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),aux)
    
      ! Compute the density-averaged magnetic field
      X_ij = ((X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-&
               X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2+&
              (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-&
               Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2+&
              (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-&
               Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2)/&
              (2.0*(sqrt(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))+&
                    sqrt(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))))

      ! Compute the square of the Roe-averaged magnetic field.
      ! Note that left and right states are interchanged!
      bPow2_ij = (ROE_MEAN_VALUE(\
                   X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                   X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2+\
                  ROE_MEAN_VALUE(\
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2+\
                  ROE_MEAN_VALUE(\
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2)/\
                  ROE_MEAN_VALUE(\
                   DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),\
                   DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),aux)

      ! Compute the magnitude of the Roe-averaged velocity
      q_ij = ROE_MEAN_VALUE(\
              X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
              X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2+\
             ROE_MEAN_VALUE(\
              Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
              Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2+\
             ROE_MEAN_VALUE(\
              Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
              Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2

      ! Compute the Roe-averaged speed of sound
#ifdef THERMALLY_IDEAL_GAS
      aPow2_ij = (2.0-GAMMA)*X_ij + (GAMMA-1.0)*(H_ij-0.5*q_ij-bPow2_ij)
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Compute auxiliary variables
      astPow2_ij = aPow2_ij+bPow2_ij
      aux = sqrt(astPow2_ij**2-4.0*aPow2_ij*bxPow2_ij)
            
      ! Compute the Roe-averagred speed of the fast waves
      cf_ij = sqrt(0.5*(astPow2_ij+aux))

      ! Scalar dissipation
      d_ij = abs(u_ij*a(1)) + anorm*cf_ij
      
      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------
      
#ifdef MHD_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*F_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*F_ij + Diff)
#endif
    end do

  end subroutine mhd_calcFluxScDiss1d_sim

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcFluxRoeDiss1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

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
#ifdef MHD_USE_IBP
    real(DP), dimension(NVAR1D) :: Fi,Fj
#else
    real(DP), dimension(NVAR1D) :: F_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP), dimension(7,7) :: Reig
    real(DP) :: pi,pj,qi,qj,ui,uj,vi,vj,wi,wj
    real(DP) :: caPow2_ij,ca_ij,cfPow2_ij,cf_ij,csPow2_ij,cs_ij
    real(DP) :: S,aux,auxf,auxs,auxy,auxz
    real(DP) :: H_ij,X_ij,q_ij,rho_ij,u_ij,v_ij,w_ij
    real(DP) :: aPow2_ij,astPow2_ij,bPow2_ij,bxPow2_ij
    real(DP) :: l1,l2,l3,l4,l5,l6,l7,w1,w2,w3,w4,w5,w6,w7
    real(DP) :: anorm
    integer, dimension(7) :: Ipiv
    integer :: idx,info

    
    do idx = 1, size(DfluxesAtEdge,3)
          
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_MHD_2T_XDIR_1D(Fi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_1D(Fj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_1D(F_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor of Roe-type
      !
      ! There are seven eigenvalues
      !
      !   u-cf, u-ca, u-cs, u, u+cs, u+ca, u+cf,
      !
      ! where u us the x-velocity component and ca, cs and cf are the
      ! velocities of th Alfveen waves, the slow and fast waves.
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
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+\
                TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))/\
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+\
                TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))/\
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx), aux)

        ! Compute the Roe-averaged density with left and right states interchanged!
        rho_ij = ROE_MEAN_VALUE(\
                  DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),\
                  DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),aux)
        
        ! Compute the square of the Roe-averaged speed of the Alfven waves.
        ! Note that left and right states are interchanged!
        bxPow2_ij = (ROE_MEAN_VALUE(\
                      X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                      X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux))**2/rho_ij
        ca_ij = sqrt(bxPow2_ij)
    
        ! Compute the density-averaged magnetic field
        X_ij = ((X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2+&
                (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2+&
                (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2)/&
                (2.0*(sqrt(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))+&
                      sqrt(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))))

        ! Compute the square of the Roe-averaged magnetic field.
        ! Note that left and right states are interchanged!
        bPow2_ij = (ROE_MEAN_VALUE(\
                     X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                     X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2+\
                    ROE_MEAN_VALUE(\
                     Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                     Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2+\
                    ROE_MEAN_VALUE(\
                     Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                     Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2)/rho_ij

        ! Compute the magnitude of the Roe-averaged velocity
        q_ij = ROE_MEAN_VALUE(\
                X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
                X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2+\
               ROE_MEAN_VALUE(\
                Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
                Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2+\
               ROE_MEAN_VALUE(\
                Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
                Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2

        ! Compute the Roe-averaged speed of sound
#ifdef THERMALLY_IDEAL_GAS
        aPow2_ij = (2.0-GAMMA)*X_ij+(GAMMA-1.0)*(H_ij-0.5*q_ij-bPow2_ij)
#else
#error "Speed of sound must be implemented!"
#endif

        ! Compute auxiliary quantities
        astPow2_ij = aPow2_ij+bPow2_ij
        aux = sqrt(astPow2_ij**2-4.0*aPow2_ij*bxPow2_ij)

        ! Compute the Roe-averagred speed of the slow and fast waves
        cfPow2_ij = 0.5*(astPow2_ij+aux); cf_ij=sqrt(cfPow2_ij)
        csPow2_ij = 0.5*(astPow2_ij-aux); cs_ij=sqrt(csPow2_ij)

        ! Compute eigenvalues
        l1 = abs(u_ij-cf_ij)
        l2 = abs(u_ij-ca_ij)
        l3 = abs(u_ij-cs_ij)
        l4 = abs(u_ij)
        l5 = abs(u_ij+cs_ij)
        l6 = abs(u_ij+ca_ij)
        l7 = abs(u_ij+cf_ij)
        
        ! Compute solution difference U_j-U_i
        Diff = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)

        ! The variable names are adopted from J.M. Stone, T.A. Gardiner,
        ! P. Teuben, J.F. Hawlay and J.B. Simon ATHENA: A New Code for
        ! Astrophysical MHD. The Astrophysocal Journal Supplement Series
        ! (ISSN 0067-0049), vol. 178, September 2008, p. 137-177.

        ! Compute the "alpha_f,s" values
        auxf = sqrt((aPow2_ij-csPow2_ij)/(cfPow2_ij-csPow2_ij))
        auxs = sqrt((cfPow2_ij-aPow2_ij)/(cfPow2_ij-csPow2_ij))

        ! Compute the "beta_"y,z" values (with left and right states interchanged!)
        auxy = ROE_MEAN_VALUE(\
                Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)
        auxz = ROE_MEAN_VALUE(\
                Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)
        aux  = sqrt(auxy**2+auxz**2)
        auxy = auxy/aux
        auxz = auxz/aux
        
        ! Compute the sign if the magnetic field
        S = sign(1.0, X_MAGNETICFIELD_CONSTANT_1D)

        ! Compute matrix of right eigenvectors       
        Reig(1,1) = auxf/aPow2_ij
        Reig(1,2) = 0.0
        Reig(1,3) = auxs/aPow2_ij
        Reig(1,4) = 1.0/aPow2_ij
        Reig(1,5) = auxs/aPow2_ij
        Reig(1,6) = 0.0
        Reig(1,7) = auxf/aPow2_ij

        Reig(2,1) = auxf*(u_ij-cf_ij)/aPow2_ij
        Reig(2,2) = 0.0
        Reig(2,3) = auxs*(u_ij-cs_ij)/aPow2_ij
        Reig(2,4) = u_ij/aPow2_ij
        Reig(2,5) = auxs*(u_ij+cs_ij)/aPow2_ij
        Reig(2,6) = 0.0
        Reig(2,7) = auxf*(u_ij+cf_ij)/aPow2_ij

        Reig(3,1) = (auxf*v_ij+auxs*cs_ij*auxy*S)/aPow2_ij
        Reig(3,2) = -rho_ij*auxz
        Reig(3,3) = (auxs*v_ij-auxf*cf_ij*auxy*S)/aPow2_ij
        Reig(3,4) = v_ij/aPow2_ij
        Reig(3,5) = (auxs*v_ij+auxf*cf_ij*auxy*S)/aPow2_ij
        Reig(3,6) = rho_ij*auxz
        Reig(3,7) = (auxf*v_ij-auxs*cs_ij*auxy*S)/aPow2_ij

        Reig(4,1) = (auxf*w_ij+auxs*cs_ij*auxz*S)/aPow2_ij
        Reig(4,2) = rho_ij*auxy
        Reig(4,3) = (auxs*w_ij-auxf*cf_ij*auxz*S)/aPow2_ij
        Reig(4,4) = w_ij/aPow2_ij
        Reig(4,5) = (auxs*w_ij+auxf*cf_ij*auxz*S)/aPow2_ij
        Reig(4,6) = -rho_ij*auxy
        Reig(4,7) = (auxf*w_ij-auxs*cs_ij*auxz*S)/aPow2_ij

        Reig(5,1) =  auxs*auxy/sqrt(rho_ij*aPow2_ij)
        Reig(5,2) = -S*sqrt(rho_ij)*auxz
        Reig(5,3) = -auxf*auxy/sqrt(rho_ij*aPow2_ij)
        Reig(5,4) =  0.0
        Reig(5,5) = -auxf*auxy/sqrt(rho_ij*aPow2_ij)
        Reig(5,6) = -S*sqrt(rho_ij)*auxz
        Reig(5,7) =  auxs*auxy/sqrt(rho_ij*aPow2_ij)

        Reig(6,1) = auxs*auxz/sqrt(rho_ij*aPow2_ij)
        Reig(6,2) = S*sqrt(rho_ij)*auxy
        Reig(6,3) = -auxf*auxz/sqrt(rho_ij*aPow2_ij)
        Reig(6,4) = 0.0
        Reig(6,5) = -auxf*auxz/sqrt(rho_ij*aPow2_ij)
        Reig(6,6) = S*sqrt(rho_ij)*auxy
        Reig(6,7) = auxs*auxz/sqrt(rho_ij*aPow2_ij)

        Reig(7,1) = (auxf*(H_ij-bPow2_ij-u_ij*cf_ij)+&
                     auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-&
                    auxs*aux/sqrt(rho_ij*aPow2_ij)
        Reig(7,2) = -rho_ij*(v_ij*auxz-w_ij*auxy)
        Reig(7,3) = (auxs*(H_ij-bPow2_ij-u_ij*cs_ij)-&
                     auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-&
                    auxf*aux/sqrt(rho_ij*aPow2_ij)
        Reig(7,4) = (0.5*q_ij+(GAMMA-2.0)/(GAMMA-1.0)*X_ij)/aPow2_ij
        Reig(7,5) = (auxs*(H_ij-bPow2_ij+u_ij*cs_ij)+&
                     auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-&
                    auxf*aux/sqrt(rho_ij*aPow2_ij)
        Reig(7,6) = rho_ij*(v_ij*auxz-w_ij*auxy)
        Reig(7,7) = (auxf*(H_ij-bPow2_ij+u_ij*cf_ij)-&
                     auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-&
                    auxs*aux/sqrt(rho_ij*aPow2_ij)

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
        Diff(7) = anorm * ( ((auxf*(H_ij-bPow2_ij-u_ij*cf_ij)+auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))*w1 +&
                             (auxf*(H_ij-bPow2_ij+u_ij*cf_ij)-auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))*w7)/aPow2_ij -&
                             auxs*aux*(w1+w7)/sqrt(rho_ij*aPow2_ij) +&
                             rho_ij*(v_ij*auxz-w_ij*auxy)*(-w2+w6) +&
                            ((auxs*(H_ij-bPow2_ij-u_ij*cs_ij)-auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))*w3 +&
                             (auxs*(H_ij-bPow2_ij+u_ij*cs_ij)+auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))*w5)/aPow2_ij -&
                             auxf*aux*(w3+w5)/sqrt(rho_ij*aPow2_ij) +&
                             (0.5*q_ij+(GAMMA-2.0)/(GAMMA-1.0)*X_ij)*w4/aPow2_ij )

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fi + Diff)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else        
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*F_ij + Diff)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*F_ij + Diff)
#endif
      else

#ifdef MHD_USE_IBP
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fj-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*Fi )
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*F_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*F_ij)
#endif
      end if
    end do

  end subroutine mhd_calcFluxRoeDiss1d_sim

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcFluxRusDiss1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)


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

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
      ! Compute total pressures
      pi = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
      pj = TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)

      ! Compute auxiliary quantities
      qi = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*ui +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*vi +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*wi
      qj = X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*uj +&
           Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*vj +&
           Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*wj

#ifdef MHD_USE_IBP
      ! Compute fluxes for x-direction
      FLUX_MHD_2T_XDIR_1D(Fi,DdataAtEdge,1,idx,ui,vi,wi,pi,qi)
      FLUX_MHD_2T_XDIR_1D(Fj,DdataAtEdge,2,idx,uj,vj,wj,pj,qj)
#else
      ! Compute flux difference for x-direction
      FLUXDIFF_MHD_2T_XDIR_1D(F_ij,DdataAtEdge,1,2,idx,ui,uj,vi,vj,wi,wj,pi,pj,qi,qj)
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
      
      ! Compute the speed of the Alfven waves
      cai = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))
      caj = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      aPow2i = GAMMA*PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)/&
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      aPow2j = GAMMA*PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)/&
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
#else
#error "Speed of sound must be implemented!"
#endif

      ! Compute auxiliary quantities
      astPow2i = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)/&
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx) + aPow2i
      astPow2j = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)/&
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx) + aPow2j

      ! Compute the speed of the fast waves
      cfi = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*cai**2)))
      cfj = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*caj**2)))
            
      ! Scalar dissipation for the Rusanov flux
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj)+&
                  abs(DmatrixCoeffsAtEdge(1,1,idx))*cfj,&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui)+&
                  abs(DmatrixCoeffsAtEdge(1,2,idx))*cfi )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef MHD_USE_IBP
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*Fj-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*Fi + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*F_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*F_ij + Diff)
#endif
    end do

  end subroutine mhd_calcFluxRusDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiagMatD1d_sim(DdataAtNode, DmatrixCoeffsAtNode,&
      IverticesAtNode, dscale, DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 1D.
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
    integer :: inode


    do inode = 1, size(DcoefficientsAtNode,3)
      
      ! Set coefficient to zero
      DcoefficientsAtNode(:,:,inode) = 0.0
    end do

  end subroutine mhd_calcMatDiagMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatDiag1d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 1D.
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

  
  end subroutine mhd_calcMatDiag1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGalMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 1D.
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


  end subroutine mhd_calcMatGalMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatGal1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 1D.
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

   
  end subroutine mhd_calcMatGal1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDissMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 1D and applies scalar artificial viscosities proportional to
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


  end subroutine mhd_calcMatScDissMatD1d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatScDiss1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 1D and applies
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

   
  end subroutine mhd_calcMatScDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDissMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
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


  end subroutine mhd_calcMatRoeDissMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRoeDiss1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 1D and applies
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

  
  end subroutine mhd_calcMatRoeDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDissMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! in 1D and applies the scalar artificial viscosities of Rusanov-type.
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

    
  end subroutine mhd_calcMatRusDissMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcMatRusDiss1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 1D applies
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

  
  end subroutine mhd_calcMatRusDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcCharacteristics1d_sim(Dweight, DdataAtEdge,&
      DcharVariablesAtEdge, DeigenvaluesAtEdge,&
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

    
  end subroutine mhd_calcCharacteristics1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTScDiss1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    
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
  !   DIMENSION(nvar,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: ui,uj,vi,vj,wi,wj
    real(DP) :: H_ij,X_ij,aux,cf_ij,d_ij,q_ij,u_ij
    real(DP) :: aPow2_ij,anorm,astPow2_ij,bPow2_ij,bxPow2_ij
    integer :: idx


    do idx = 1, size(DfluxesAtEdge,2)

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

      ! Compute skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))
      anorm = abs(a(1)) ! = sqrt(a(1)*a(1))

      ! Compute Roe mean values
      aux  = ROE_MEAN_RATIO(\
              DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
              DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
      u_ij = ROE_MEAN_VALUE(ui,uj,aux)
      H_ij = ROE_MEAN_VALUE(\
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+\
              TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))/\
              DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
             (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+\
              TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))/\
              DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx), aux)
    
      ! Compute the square of the Roe-averaged speed of the Alfven waves.
      ! Note that left and right states are interchanged!
      bxPow2_ij = (ROE_MEAN_VALUE(\
                    X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                    X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux))**2/\
                  ROE_MEAN_VALUE(\
                   DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),\
                   DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),aux)
    
      ! Compute the density-averaged magnetic field
      X_ij = ((X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-&
               X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2+&
              (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-&
               Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2+&
              (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-&
               Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2)/&
              (2.0*(sqrt(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))+&
                    sqrt(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))))

      ! Compute the square of the Roe-averaged magnetic field.
      ! Note that left and right states are interchanged!
      bPow2_ij = (ROE_MEAN_VALUE(\
                   X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                   X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2+\
                  ROE_MEAN_VALUE(\
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                   Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2+\
                  ROE_MEAN_VALUE(\
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                   Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2)/\
                  ROE_MEAN_VALUE(\
                   DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),\
                   DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),aux)

      ! Compute the magnitude of the Roe-averaged velocity
      q_ij = ROE_MEAN_VALUE(\
              X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
              X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2+\
             ROE_MEAN_VALUE(\
              Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
              Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2+\
             ROE_MEAN_VALUE(\
              Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
              Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2


      ! Compute the Roe-averaged speed of sound
#ifdef THERMALLY_IDEAL_GAS
      aPow2_ij = (2.0-GAMMA)*X_ij + (GAMMA-1.0)*(H_ij-0.5*q_ij-bPow2_ij)
#else
#error "Speed of sound must be implemented!"
#endif
      
      ! Compute auxiliary variables
      astPow2_ij = aPow2_ij+bPow2_ij
      aux = sqrt(astPow2_ij**2-4.0*aPow2_ij*bxPow2_ij)
            
      ! Compute the Roe-averagred speed of the fast waves
      cf_ij = sqrt(0.5*(astPow2_ij+aux))

      ! Scalar dissipation
      d_ij = abs(u_ij*a(1)) + anorm*cf_ij
      
      ! Compute conservative fluxes
      DfluxesAtEdge(:,idx) = dscale*d_ij*(DdataAtEdge(:,1,idx)-DdataAtEdge(:,2,idx))
    end do

  end subroutine mhd_calcFluxFCTScDiss1d_sim

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcFluxFCTRoeDiss1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

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
  !   DIMENSION(nvar,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP), dimension(7,7) :: Reig
    real(DP) :: ui,uj,vi,vj,wi,wj
    real(DP) :: caPow2_ij,ca_ij,cfPow2_ij,cf_ij,csPow2_ij,cs_ij
    real(DP) :: S,aux,auxf,auxs,auxy,auxz
    real(DP) :: H_ij,X_ij,q_ij,rho_ij,u_ij,v_ij,w_ij
    real(DP) :: aPow2_ij,astPow2_ij,bPow2_ij,bxPow2_ij
    real(DP) :: l1,l2,l3,l4,l5,l6,l7,w1,w2,w3,w4,w5,w6,w7
    real(DP) :: anorm
    integer, dimension(7) :: Ipiv
    integer :: idx,info

    
    do idx = 1, size(DfluxesAtEdge,2)

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

      ! Compute skew-symmetric coefficient and its norm
      a = 0.5*(DmatrixCoeffsAtEdge(1,1,idx)-DmatrixCoeffsAtEdge(1,2,idx))
      anorm = abs(a(1)) ! = sqrt(a(1)*a(1))

      if (anorm .gt. SYS_EPSREAL) then

        ! Compute Roe mean values
        aux  = ROE_MEAN_RATIO(\
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx))
        u_ij = ROE_MEAN_VALUE(ui,uj,aux)
        v_ij = ROE_MEAN_VALUE(vi,vj,aux)
        w_ij = ROE_MEAN_VALUE(wi,wj,aux)
        H_ij = ROE_MEAN_VALUE(\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)+\
                TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))/\
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
               (TOTAL_ENERGY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)+\
                TOTAL_PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))/\
                DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx), aux)

        ! Compute the Roe-averaged density with left and right states interchanged!
        rho_ij = ROE_MEAN_VALUE(\
                  DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),\
                  DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),aux)
        
        ! Compute the square of the Roe-averaged speed of the Alfven waves.
        ! Note that left and right states are interchanged!
        bxPow2_ij = (ROE_MEAN_VALUE(\
                      X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                      X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux))**2/rho_ij
        ca_ij = sqrt(bxPow2_ij)
    
        ! Compute the density-averaged magnetic field
        X_ij = ((X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-&
                 X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2+&
                (Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-&
                 Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2+&
                (Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)-&
                 Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))**2)/&
                (2.0*(sqrt(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))+&
                      sqrt(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))))

        ! Compute the square of the Roe-averaged magnetic field.
        ! Note that left and right states are interchanged!
        bPow2_ij = (ROE_MEAN_VALUE(\
                     X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                     X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2+\
                    ROE_MEAN_VALUE(\
                     Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                     Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2+\
                    ROE_MEAN_VALUE(\
                     Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                     Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)**2)/rho_ij

        ! Compute the magnitude of the Roe-averaged velocity
        q_ij = ROE_MEAN_VALUE(\
                X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
                X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2+\
               ROE_MEAN_VALUE(\
                Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
                Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2+\
               ROE_MEAN_VALUE(\
                Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx),\
                Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx),aux)**2

        ! Compute the Roe-averaged speed of sound
#ifdef THERMALLY_IDEAL_GAS
        aPow2_ij = (2.0-GAMMA)*X_ij+(GAMMA-1.0)*(H_ij-0.5*q_ij-bPow2_ij)
#else
#error "Speed of sound must be implemented!"
#endif

        ! Compute auxiliary quantities
        astPow2_ij = aPow2_ij+bPow2_ij
        aux = sqrt(astPow2_ij**2-4.0*aPow2_ij*bxPow2_ij)

        ! Compute the Roe-averagred speed of the slow and fast waves
        cfPow2_ij = 0.5*(astPow2_ij+aux); cf_ij=sqrt(cfPow2_ij)
        csPow2_ij = 0.5*(astPow2_ij-aux); cs_ij=sqrt(csPow2_ij)

        ! Compute eigenvalues
        l1 = abs(u_ij-cf_ij)
        l2 = abs(u_ij-ca_ij)
        l3 = abs(u_ij-cs_ij)
        l4 = abs(u_ij)
        l5 = abs(u_ij+cs_ij)
        l6 = abs(u_ij+ca_ij)
        l7 = abs(u_ij+cf_ij)
        
        ! Compute solution difference U_i-U_j
        Diff = DdataAtEdge(:,1,idx)-DdataAtEdge(:,2,idx)

        ! The variable names are adopted from J.M. Stone, T.A. Gardiner,
        ! P. Teuben, J.F. Hawlay and J.B. Simon ATHENA: A New Code for
        ! Astrophysical MHD. The Astrophysocal Journal Supplement Series
        ! (ISSN 0067-0049), vol. 178, September 2008, p. 137-177.

        ! Compute the "alpha_f,s" values
        auxf = sqrt((aPow2_ij-csPow2_ij)/(cfPow2_ij-csPow2_ij))
        auxs = sqrt((cfPow2_ij-aPow2_ij)/(cfPow2_ij-csPow2_ij))

        ! Compute the "beta_"y,z" values (with left and right states interchanged!)
        auxy = ROE_MEAN_VALUE(\
                Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)
        auxz = ROE_MEAN_VALUE(\
                Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx),\
                Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx),aux)
        aux  = sqrt(auxy**2+auxz**2)
        auxy = auxy/aux
        auxz = auxz/aux
        
        ! Compute the sign if the magnetic field
        S = sign(1.0, X_MAGNETICFIELD_CONSTANT_1D)

        ! Compute matrix of right eigenvectors       
        Reig(1,1) = auxf/aPow2_ij
        Reig(1,2) = 0.0
        Reig(1,3) = auxs/aPow2_ij
        Reig(1,4) = 1.0/aPow2_ij
        Reig(1,5) = auxs/aPow2_ij
        Reig(1,6) = 0.0
        Reig(1,7) = auxf/aPow2_ij

        Reig(2,1) = auxf*(u_ij-cf_ij)/aPow2_ij
        Reig(2,2) = 0.0
        Reig(2,3) = auxs*(u_ij-cs_ij)/aPow2_ij
        Reig(2,4) = u_ij/aPow2_ij
        Reig(2,5) = auxs*(u_ij+cs_ij)/aPow2_ij
        Reig(2,6) = 0.0
        Reig(2,7) = auxf*(u_ij+cf_ij)/aPow2_ij

        Reig(3,1) = (auxf*v_ij+auxs*cs_ij*auxy*S)/aPow2_ij
        Reig(3,2) = -rho_ij*auxz
        Reig(3,3) = (auxs*v_ij-auxf*cf_ij*auxy*S)/aPow2_ij
        Reig(3,4) = v_ij/aPow2_ij
        Reig(3,5) = (auxs*v_ij+auxf*cf_ij*auxy*S)/aPow2_ij
        Reig(3,6) = rho_ij*auxz
        Reig(3,7) = (auxf*v_ij-auxs*cs_ij*auxy*S)/aPow2_ij

        Reig(4,1) = (auxf*w_ij+auxs*cs_ij*auxz*S)/aPow2_ij
        Reig(4,2) = rho_ij*auxy
        Reig(4,3) = (auxs*w_ij-auxf*cf_ij*auxz*S)/aPow2_ij
        Reig(4,4) = w_ij/aPow2_ij
        Reig(4,5) = (auxs*w_ij+auxf*cf_ij*auxz*S)/aPow2_ij
        Reig(4,6) = -rho_ij*auxy
        Reig(4,7) = (auxf*w_ij-auxs*cs_ij*auxz*S)/aPow2_ij

        Reig(5,1) =  auxs*auxy/sqrt(rho_ij*aPow2_ij)
        Reig(5,2) = -S*sqrt(rho_ij)*auxz
        Reig(5,3) = -auxf*auxy/sqrt(rho_ij*aPow2_ij)
        Reig(5,4) =  0.0
        Reig(5,5) = -auxf*auxy/sqrt(rho_ij*aPow2_ij)
        Reig(5,6) = -S*sqrt(rho_ij)*auxz
        Reig(5,7) =  auxs*auxy/sqrt(rho_ij*aPow2_ij)

        Reig(6,1) = auxs*auxz/sqrt(rho_ij*aPow2_ij)
        Reig(6,2) = S*sqrt(rho_ij)*auxy
        Reig(6,3) = -auxf*auxz/sqrt(rho_ij*aPow2_ij)
        Reig(6,4) = 0.0
        Reig(6,5) = -auxf*auxz/sqrt(rho_ij*aPow2_ij)
        Reig(6,6) = S*sqrt(rho_ij)*auxy
        Reig(6,7) = auxs*auxz/sqrt(rho_ij*aPow2_ij)

        Reig(7,1) = (auxf*(H_ij-bPow2_ij-u_ij*cf_ij)+&
                     auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-&
                    auxs*aux/sqrt(rho_ij*aPow2_ij)
        Reig(7,2) = -rho_ij*(v_ij*auxz-w_ij*auxy)
        Reig(7,3) = (auxs*(H_ij-bPow2_ij-u_ij*cs_ij)-&
                     auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-&
                    auxf*aux/sqrt(rho_ij*aPow2_ij)
        Reig(7,4) = (0.5*q_ij+(GAMMA-2.0)/(GAMMA-1.0)*X_ij)/aPow2_ij
        Reig(7,5) = (auxs*(H_ij-bPow2_ij+u_ij*cs_ij)+&
                     auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-&
                    auxf*aux/sqrt(rho_ij*aPow2_ij)
        Reig(7,6) = rho_ij*(v_ij*auxz-w_ij*auxy)
        Reig(7,7) = (auxf*(H_ij-bPow2_ij+u_ij*cf_ij)-&
                     auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))/aPow2_ij-&
                    auxs*aux/sqrt(rho_ij*aPow2_ij)

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
        Diff(7) = anorm * ( ((auxf*(H_ij-bPow2_ij-u_ij*cf_ij)+auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))*w1 +&
                             (auxf*(H_ij-bPow2_ij+u_ij*cf_ij)-auxs*cs_ij*S*(v_ij*auxy+w_ij*auxz))*w7)/aPow2_ij -&
                             auxs*aux*(w1+w7)/sqrt(rho_ij*aPow2_ij) +&
                             rho_ij*(v_ij*auxz-w_ij*auxy)*(-w2+w6) +&
                            ((auxs*(H_ij-bPow2_ij-u_ij*cs_ij)-auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))*w3 +&
                             (auxs*(H_ij-bPow2_ij+u_ij*cs_ij)+auxf*cf_ij*S*(v_ij*auxy+w_ij*auxz))*w5)/aPow2_ij -&
                             auxf*aux*(w3+w5)/sqrt(rho_ij*aPow2_ij) +&
                             (0.5*q_ij+(GAMMA-2.0)/(GAMMA-1.0)*X_ij)*w4/aPow2_ij )

        ! Compute antidiffusive flux
        DfluxesAtEdge(:,idx) = dscale*Diff
      else
        ! Clear antidiffusive flux
        DfluxesAtEdge(:,idx) = 0
      end if
    end do

  end subroutine mhd_calcFluxFCTRoeDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_calcFluxFCTRusDiss1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

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
    

    do idx = 1, size(DfluxesAtEdge,2)

      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

      ! Compute the speed of the Alfven waves
      cai = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx))
      caj = abs(X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx))

      ! Compute the speed of sound
#ifdef THERMALLY_IDEAL_GAS
      aPow2i = GAMMA*PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)/&
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      aPow2j = GAMMA*PRESSURE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)/&
               DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
#else
#error "Speed of sound must be implemented!"
#endif

      ! Compute auxiliary quantities
      astPow2i = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)/&
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx) + aPow2i
      astPow2j = MAGNETICFIELD_MAGNITUDE_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)/&
                 DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx) + aPow2j

      ! Compute the speed of the fast waves
      cfi = sqrt(0.5*(astPow2i + sqrt(astPow2i**2 - 4.0*aPow2i*cai**2)))
      cfj = sqrt(0.5*(astPow2j + sqrt(astPow2j**2 - 4.0*aPow2j*caj**2)))
            
      ! Scalar dissipation for the Rusanov flux
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj)+&
                  abs(DmatrixCoeffsAtEdge(1,1,idx))*cfj,&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui)+&
                  abs(DmatrixCoeffsAtEdge(1,2,idx))*cfi )

      ! Compute conservative fluxes
      DfluxesAtEdge(:,idx) = dscale*d_ij*(DdataAtEdge(:,1,idx)-DdataAtEdge(:,2,idx))
    end do

  end subroutine mhd_calcFluxFCTRusDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDensity1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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

  end subroutine mhd_trafoFluxDensity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDensity1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density in 1D.
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

  end subroutine mhd_trafoDiffDensity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxEnergy1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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

  end subroutine mhd_trafoFluxEnergy1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffEnergy1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the energy in 1D.
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

  end subroutine mhd_trafoDiffEnergy1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxPressure1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

      
      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(1,1,idx) = (GAMMA-1.0)*&
          (0.5*(ui*ui+vi*vi+wi*wi)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          vi*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          wi*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
          X_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
          Y_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
          Z_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
      DtransformedFluxesAtEdge(1,2,idx) =-(GAMMA-1.0)*&
          (0.5*(uj*uj+vj*vj+wj*wj)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          vj*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          wj*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
          X_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
          Y_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
          Z_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do
    
  end subroutine mhd_trafoFluxPressure1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffPressure1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the pressure in 1D.
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
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine mhd_trafoDiffPressure1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxVelocity1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)


      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(1,1,idx) =&
          (X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           ui*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -(X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           uj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

      ! Transformed velocity fluxes in y-direction
      DtransformedFluxesAtEdge(2,1,idx) =&
          (Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           vi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -(Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           vj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

      ! Transformed velocity fluxes in z-direction
      DtransformedFluxesAtEdge(3,1,idx) =&
          (Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           wi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      DtransformedFluxesAtEdge(3,2,idx) =&
         -(Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           wj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
    end do
    
  end subroutine mhd_trafoFluxVelocity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffVelocity1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the velocity in 1D.
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

      ! Transformed velocity difference in y-direction
      DtransformedDataAtEdge(2,idx) =&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)

      ! Transformed velocity difference in z-direction
      DtransformedDataAtEdge(3,idx) =&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine mhd_trafoDiffVelocity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxMomentum1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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

      ! Transformed momentum fluxes in y-direction
      DtransformedFluxesAtEdge(2,1,idx) =&
          Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)

      ! Transformed momentum fluxes in z-direction
      DtransformedFluxesAtEdge(3,1,idx) =&
          Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(3,2,idx) =&
         -Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
    end do
    
  end subroutine mhd_trafoFluxMomentum1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffMomentum1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the momentum in 1D.
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

      ! Transformed momentum difference in y-direction
      DtransformedDataAtEdge(2,idx) =&
          Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          Y_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)

      ! Transformed momentum difference in z-direction
      DtransformedDataAtEdge(3,idx) =&
          Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          Z_MOMENTUM_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do
    
  end subroutine mhd_trafoDiffMomentum1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDenEng1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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

  end subroutine mhd_trafoFluxDenEng1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDenEng1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density and energy in 1D.
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

  end subroutine mhd_trafoDiffDenEng1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDenPre1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
           
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) =&
          DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)

      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(2,1,idx) = (GAMMA-1.0)*&
          (0.5*(ui*ui+vi*vi+wi*wi)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          vi*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          wi*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
          X_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
          Y_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
          Z_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
      DtransformedFluxesAtEdge(2,2,idx) =-(GAMMA-1.0)*&
          (0.5*(uj*uj+vj*vj+wj*wj)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          vj*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          wj*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
          X_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
          Y_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
          Z_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do

  end subroutine mhd_trafoFluxDenPre1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDenPre1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of the given
    ! conservative variables to differences for the density and energy in 1D.
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
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine mhd_trafoDiffDenPre1d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxDenPreVel1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      uj = X_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      vi = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      vj = Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      wi = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      wj = Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)
      
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

      ! Transformed velocity fluxes in y-direction
      DtransformedFluxesAtEdge(3,1,idx) =&
          (Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           vi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      DtransformedFluxesAtEdge(3,2,idx) =&
         -(Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           vj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

      ! Transformed velocity fluxes in z-direction
      DtransformedFluxesAtEdge(4,1,idx) =&
          (Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           wi*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      DtransformedFluxesAtEdge(4,2,idx) =&
         -(Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
           wj*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))/&
           DENSITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)

      ! Transformed pressure fluxes
#ifdef PERFECT_GAS
      DtransformedFluxesAtEdge(5,1,idx) = (GAMMA-1.0)*&
          (0.5*(ui*ui+vi*vi+wi*wi)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          ui*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          vi*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          wi*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
          X_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
          Y_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)*&
          Z_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
      DtransformedFluxesAtEdge(5,2,idx) =-(GAMMA-1.0)*&
          (0.5*(uj*uj+vj*vj+wj*wj)*DENSITY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          uj*X_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          vj*Y_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          wj*Z_MOMENTUM_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx)-&
          X_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
          X_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
          Y_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)*&
          Z_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)-&
          TOTAL_ENERGY_1T_FROM_CONSVAR(DfluxesAtEdge,NVAR1D,idx))
#else
#error "Pressure for nonperfect gas must be implemented!"
#endif
    end do
      
  end subroutine mhd_trafoFluxDenPreVel1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffDenPreVel1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

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

      ! Transformed velocity difference in y-direction
      DtransformedDataAtEdge(3,idx) =&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          Y_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)

      ! Transformed velocity difference in z-direction
      DtransformedDataAtEdge(4,idx) =&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          Z_VELOCITY_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(5,idx) =&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,2,idx)-&
          PRESSURE_2T_FROM_CONSVAR(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine mhd_trafoDiffDenPreVel1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoFluxMagfield1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

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

      ! Transformed magnetic field fluxes in y-direction
      DtransformedFluxesAtEdge(1,1,idx) =&
          Y_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(1,2,idx) =&
         -Y_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)

      ! Transformed magnetic field fluxes in z-direction
      DtransformedFluxesAtEdge(2,1,idx) =&
          Z_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)
      DtransformedFluxesAtEdge(2,2,idx) =&
         -Z_MAGNETICFIELD_1T_FROM_CONSVAR_1D(DfluxesAtEdge,NVAR1D,idx)
    end do

  end subroutine mhd_trafoFluxMagfield1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine mhd_trafoDiffMagfield1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

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
      
      ! Transformed magnetic field difference in y-direction
      DtransformedDataAtEdge(1,idx) =&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)-&
          Y_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)

      ! Transformed magnetic field difference in z-direction
      DtransformedDataAtEdge(2,idx) =&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,2,idx)-&
          Z_MAGNETICFIELD_2T_FROM_CONSVAR_1D(DdataAtEdge,NVAR1D,1,idx)
    end do

  end subroutine mhd_trafoDiffMagfield1d_sim

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

  subroutine mhd_hadaptCallbackScalar1d(iOperation, rcollection)

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
            OU_CLASS_WARNING,OU_MODE_STD,'mhd_hadaptCallbackScalar1d')
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

  end subroutine mhd_hadaptCallbackScalar1d

  !*****************************************************************************

!<subroutine>

  subroutine mhd_hadaptCallbackBlock1d(iOperation, rcollection)

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
            OU_CLASS_WARNING,OU_MODE_STD,'mhd_hadaptCallbackBlock1d')
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

  end subroutine mhd_hadaptCallbackBlock1d

end module mhd_callback1d
