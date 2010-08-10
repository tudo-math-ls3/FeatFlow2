!##############################################################################
!# ****************************************************************************
!# <name> euler_callback2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible Euler/Navier-Stokes equations in 2D.
!#
!# The following callback functions are available:
!#
!# 1.) euler_calcFluxGal2d_sim
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#
!# 2.) euler_calcFluxGalNoBdr2d_sim
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#        without assembling the symmetric boundary contribution
!#
!# 3.) euler_calcFluxScDiss2d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting scalar artificial viscosities
!#
!# 4.) euler_calcFluxScDissDiSp2d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting scalar artificial viscosities based on
!#        dimensional splitting approach
!#
!# 5.) euler_calcFluxRoeDiss2d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities
!#
!# 6.) euler_calcFluxRoeDissDiSp2d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities based on
!#        dimensional splitting approach
!#
!# 7.) euler_calcFluxRusDiss2d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion based on
!#        dimensional splitting approach
!#
!# 8.) euler_calcFluxRusDiSp2d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion
!#
!# 9.) euler_calcMatDiagMatD2d_sim
!#     -> Computes local matrix for diagonal entry
!#
!# 10.) euler_calcMatDiag2d_sim
!#      -> Computes local matrix for diagonal entry
!#
!# 11.) euler_calcMatGalMatD2d_sim
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 12.) euler_calcMatGal2d_sim
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 13.) euler_calcMatScDissMatD2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 14.) euler_calcMatScDiss2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 15.) euler_calcMatRoeDissMatD2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 16.) euler_calcMatRoeDiss2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 17.) euler_calcMatRusDissMatD2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov artificial viscosities
!#
!# 18.) euler_calcMatRusDiss2d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov flux artificial viscosities
!#
!# 19.) euler_calcCharacteristics2d_sim
!#      -> Computes characteristic variables
!#
!# 20.) euler_calcFluxFCTScalarDiss2d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting scalar artificial viscosities
!#
!# 21.) euler_calcFluxFCTTensorDiss2d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting tensorial artificial viscosities
!#
!# 22.) euler_calcFluxFCTRusanov2d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting the Rusanov artificial viscosities
!#
!# 23.) euler_trafoFluxDensity2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 24.) euler_trafoDiffDensity2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 25.) euler_trafoFluxEnergy2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 26.) euler_trafoDiffEnergy2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 27.) euler_trafoFluxPressure2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 28.) euler_trafoFluxVelocity2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 29.) euler_trafoDiffVelocity2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 30.) euler_trafoFluxMomentum2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 31.) euler_trafoDiffMomentum2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 32.) euler_trafoFluxDenEng2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 33.) euler_trafoDiffDenEng2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 34.) euler_trafoFluxDenPre2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 35.) euler_trafoDiffDenPre2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 36.) euler_trafoFluxDenPreVel2d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 37.) euler_trafoDiffDenPreVel2d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure 
!#         and the velocity
!#
!# 38.) euler_calcBoundaryvalues2d
!#      -> Computes the boundary values for a given node
!#
!# 39.) euler_hadaptCallbackScalar2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in interleave format
!#
!# 40.) euler_hadaptCallbackBlock2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in block format
!#
!# 41.) euler_coeffVectorBdr2d_sim
!#      -> Calculates the coefficients for the linear form in 2D
!#
!# 42.) euler_coeffMatrixBdr2d_sim
!#      -> Calculates the coefficients for the bilinear form in 2D
!# </purpose>
!##############################################################################

module euler_callback2d

  use boundaryfilter
  use collection
  use derivatives
  use euler_basic
  use flagship_callback
  use fsystem
  use genoutput
  use graph
  use groupfemsystem
  use hadaptaux
  use linearsystemblock
  use linearsystemscalar
  use problem
  use solveraux
  use storage
  use thermodynamics

  implicit none

  private
  public :: euler_calcFluxGal2d_sim
  public :: euler_calcFluxGalNoBdr2d_sim
  public :: euler_calcFluxScDiss2d_sim
  public :: euler_calcFluxScDissDiSp2d_sim
  public :: euler_calcFluxRoeDiss2d_sim
  public :: euler_calcFluxRoeDissDiSp2d_sim
  public :: euler_calcFluxRusDiss2d_sim
  public :: euler_calcFluxRusDissDiSp2d_sim
  public :: euler_calcMatDiagMatD2d_sim
  public :: euler_calcMatDiag2d_sim
  public :: euler_calcMatGalMatD2d_sim
  public :: euler_calcMatGal2d_sim
  public :: euler_calcMatScDissMatD2d_sim
  public :: euler_calcMatScDiss2d_sim
  public :: euler_calcMatRoeDissMatD2d_sim
  public :: euler_calcMatRoeDiss2d_sim
  public :: euler_calcMatRusDissMatD2d_sim
  public :: euler_calcMatRusDiss2d_sim
  public :: euler_calcCharacteristics2d_sim
  public :: euler_calcFluxFCTScalarDiss2d
  public :: euler_calcFluxFCTTensorDiss2d
  public :: euler_calcFluxFCTRusanov2d
  public :: euler_trafoFluxDensity2d_sim
  public :: euler_trafoFluxEnergy2d_sim
  public :: euler_trafoFluxPressure2d_sim
  public :: euler_trafoFluxVelocity2d_sim
  public :: euler_trafoFluxMomentum2d_sim
  public :: euler_trafoFluxDenEng2d_sim
  public :: euler_trafoFluxDenPre2d_sim
  public :: euler_trafoFluxDenPreVel2d_sim
  public :: euler_trafoDiffDensity2d_sim
  public :: euler_trafoDiffEnergy2d_sim
  public :: euler_trafoDiffPressure2d_sim
  public :: euler_trafoDiffVelocity2d_sim
  public :: euler_trafoDiffMomentum2d_sim
  public :: euler_trafoDiffDenEng2d_sim
  public :: euler_trafoDiffDenPre2d_sim
  public :: euler_trafoDiffDenPreVel2d_sim
  public :: euler_calcBoundaryvalues2d
  public :: euler_coeffVectorBdr2d_sim
  public :: euler_hadaptCallbackScalar2d
  public :: euler_hadaptCallbackBlock2d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxGal2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the standard
    ! Galerkin discretisation in 2D.
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
#ifdef USE_EULER_INTEGRATEBYPARTS
    real(DP), dimension(NVAR2D) :: dF1_i, dF2_i, dF1_j, dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
#endif
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      !
      !      / rho*u       \              / rho*v     \
      ! F1 = | rho*u*u + p |   and   F2 = | rho*u*v   |
      !      | rho*u*v     |              | rho*v*v+p |
      !      \ rho*u*H     /              \ rho*v*H   /
      !
      ! Here, we do not compute the pressure p and the enthalpy H but we
      ! calculate the fluxes from the conservative variables as follows:
      !
      !      / U2                           \
      ! F1 = | G1*U4-G14*ru2i-G2*rv2i       |
      !      | U3*ui                        |
      !      \ (gamma*U4-G2*(ru2i+rv2i))*ui /
      !
      !      / U3                           \
      ! F2 = | U3*ui = U2*vi                |
      !      | G1*U4-G14*rv2i-G2*ru2i       |
      !      \ (gamma*U4-G2*(ru2i+rv2i))*vi /
      !
      ! where the auxiliary values for node i are defined as follows:
      !
      ! ru2i = U2*U2/U1 = ui*U2
      ! rv2i = U3*U3/U1 = vi*U3
      ! ui = U2/U1
      ! vi = U3/U1
      !
      ! and the predefined constants are given by:
      !
      ! G14 = (gamma-3)/2   and   G2 = (gamma-1)/2   and   G1 = gamma-1
      !
      ! The auxiliary values for node j are defined accordingly.
      ! -------------------------------------------------------------------------
      
      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)

      ! Compute auxiliary variables
      ru2i = ui*DdataAtEdge(2,1,idx)
      rv2i = vi*DdataAtEdge(3,1,idx)
      ru2j = uj*DdataAtEdge(2,2,idx)
      rv2j = vj*DdataAtEdge(3,2,idx)
      
#ifdef USE_EULER_INTEGRATEBYPARTS
      ! Compute fluxes for x-direction
      dF1_i(1) = DdataAtEdge(2,1,idx)
      dF1_i(2) = G1*DdataAtEdge(4,1,idx)-G14*ru2i-G2*rv2i
      dF1_i(3) = DdataAtEdge(3,1,idx)*ui
      dF1_i(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*ui
      
      dF1_j(1) = DdataAtEdge(2,2,idx)
      dF1_j(2) = G1*DdataAtEdge(4,2,idx)-G14*ru2j-G2*rv2j
      dF1_j(3) = DdataAtEdge(3,2,idx)*uj
      dF1_j(4) = (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*uj
      
      ! Compute fluxes for y-direction
      dF2_i(1) = DdataAtEdge(3,1,idx)
      dF2_i(2) = DdataAtEdge(3,1,idx)*ui
      dF2_i(3) = G1*DdataAtEdge(4,1,idx)-G14*rv2i-G2*ru2i
      dF2_i(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*vi
      
      dF2_j(1) = DdataAtEdge(3,2,idx)
      dF2_j(2) = DdataAtEdge(3,2,idx)*uj
      dF2_j(3) = (G1*DdataAtEdge(4,2,idx)-G14*rv2j-G2*ru2j)
      dF2_j(4) = (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*vj

      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*dF2_i )
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = DdataAtEdge(2,1,idx) - DdataAtEdge(2,2,idx)
      dF1_ij(2) = (G1*DdataAtEdge(4,1,idx)-G14*ru2i-G2*rv2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*ru2j-G2*rv2j)
      dF1_ij(3) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF1_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*ui-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*uj
      
      ! Compute flux difference for y-direction
      dF2_ij(1) = DdataAtEdge(3,1,idx) - DdataAtEdge(3,2,idx)
      dF2_ij(2) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF2_ij(3) = (G1*DdataAtEdge(4,1,idx)-G14*rv2i-G2*ru2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*rv2j-G2*ru2j)
      dF2_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*vi-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*vj
      
      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij)
#endif

    end do

  end subroutine euler_calcFluxGal2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxGalNoBdr2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the TVD
    ! discretisation in 2D. The symmetric boundary contributions
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
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j
    integer :: idx

    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      ! For a detailed description of algorithm and the definition of auxiliary
      ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
      !-------------------------------------------------------------------------
      
      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)

      ! Compute auxiliary variables
      ru2i = ui*DdataAtEdge(2,1,idx)
      rv2i = vi*DdataAtEdge(3,1,idx)
      ru2j = uj*DdataAtEdge(2,2,idx)
      rv2j = vj*DdataAtEdge(3,2,idx)

      ! Compute flux difference for x-direction
      dF1_ij(1) = DdataAtEdge(2,1,idx) - DdataAtEdge(2,2,idx)
      dF1_ij(2) = (G1*DdataAtEdge(4,1,idx)-G14*ru2i-G2*rv2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*ru2j-G2*rv2j)
      dF1_ij(3) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF1_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*ui-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*uj
      
      ! Compute flux difference for y-direction
      dF2_ij(1) = DdataAtEdge(3,1,idx) - DdataAtEdge(3,2,idx)
      dF2_ij(2) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF2_ij(3) = (G1*DdataAtEdge(4,1,idx)-G14*rv2i-G2*ru2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*rv2j-G2*ru2j)
      dF2_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*vi-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*vj
      
      ! Assemble symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * ((DmatrixCoeffsAtEdge(1,1,idx)-&
                                            DmatrixCoeffsAtEdge(1,2,idx))/2._DP*dF1_ij+&
                                           (DmatrixCoeffsAtEdge(2,1,idx)-&
                                            DmatrixCoeffsAtEdge(2,2,idx))/2._DP*dF2_ij)
      DfluxesAtEdge(:,2,idx) = DfluxesAtEdge(:,1,idx)

    end do

  end subroutine euler_calcFluxGalNoBdr2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxScDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using scalar dissipation.
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
#ifdef USE_EULER_INTEGRATEBYPARTS
    real(DP), dimension(NVAR2D) :: dF1_i, dF2_i, dF1_j, dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j
    real(DP) :: d_ij,hi,hj,H_ij,q_ij,u_ij,v_ij,aux,vel,c_ij
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      ! For a detailed description of algorithm and the definition of auxiliary
      ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)

      ! Compute auxiliary variables
      ru2i = ui*DdataAtEdge(2,1,idx)
      rv2i = vi*DdataAtEdge(3,1,idx)
      ru2j = uj*DdataAtEdge(2,2,idx)
      rv2j = vj*DdataAtEdge(3,2,idx)

#ifdef USE_EULER_INTEGRATEBYPARTS
      ! Compute fluxes for x-direction
      dF1_i(1) = DdataAtEdge(2,1,idx)
      dF1_i(2) = G1*DdataAtEdge(4,1,idx)-G14*ru2i-G2*rv2i
      dF1_i(3) = DdataAtEdge(3,1,idx)*ui
      dF1_i(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*ui
      
      dF1_j(1) = DdataAtEdge(2,2,idx)
      dF1_j(2) = G1*DdataAtEdge(4,2,idx)-G14*ru2j-G2*rv2j
      dF1_j(3) = DdataAtEdge(3,2,idx)*uj
      dF1_j(4) = (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*uj
      
      ! Compute fluxes for y-direction
      dF2_i(1) = DdataAtEdge(3,1,idx)
      dF2_i(2) = DdataAtEdge(3,1,idx)*ui
      dF2_i(3) = G1*DdataAtEdge(4,1,idx)-G14*rv2i-G2*ru2i
      dF2_i(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*vi
      
      dF2_j(1) = DdataAtEdge(3,2,idx)
      dF2_j(2) = DdataAtEdge(3,2,idx)*uj
      dF2_j(3) = (G1*DdataAtEdge(4,2,idx)-G14*rv2j-G2*ru2j)
      dF2_j(4) = (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*vj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = DdataAtEdge(2,1,idx) - DdataAtEdge(2,2,idx)
      dF1_ij(2) = (G1*DdataAtEdge(4,1,idx)-G14*ru2i-G2*rv2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*ru2j-G2*rv2j)
      dF1_ij(3) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF1_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*ui-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*uj
      
      ! Compute flux difference for y-direction
      dF2_ij(1) = DdataAtEdge(3,1,idx) - DdataAtEdge(3,2,idx)
      dF2_ij(2) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF2_ij(3) = (G1*DdataAtEdge(4,1,idx)-G14*rv2i-G2*ru2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*rv2j-G2*ru2j)
      dF2_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*vi-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*vj
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral radius
      !-------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient
      a = 0.5_DP*(DmatrixCoeffsAtEdge(:,1,idx)-&
                  DmatrixCoeffsAtEdge(:,2,idx))
      
      ! Compute Roe mean values
      aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
      u_ij = (aux*ui+uj)/(aux+1.0_DP)
      v_ij = (aux*vi+vj)/(aux+1.0_DP)
      hi   = GAMMA*DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui+vi*vi)
      hj   = GAMMA*DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj+vj*vj)
      H_ij = (aux*hi+hj)/(aux+1.0_DP)
      
      ! Compute auxiliary variables
      aux  = sqrt(a(1)*a(1)+a(2)*a(2))
      vel  = u_ij*a(1) + v_ij*a(2)
      q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
      c_ij = sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL))
      
      ! Scalar dissipation
      d_ij = abs(vel) + aux*c_ij

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef USE_EULER_INTEGRATEBYPARTS
      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*dF2_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij + Diff)
#endif

    end do

  end subroutine euler_calcFluxScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxScDissDiSp2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using scalar dissipation,
    ! whereby dimensional splitting is employed.
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
#ifdef USE_EULER_INTEGRATEBYPARTS
    real(DP), dimension(NVAR2D) :: dF1_i, dF2_i, dF1_j, dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j
    real(DP) :: d_ij,hi,hj,H_ij,q_ij,u_ij,v_ij,aux
    integer :: idx
    
    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      ! For a detailed description of algorithm and the definition of auxiliary
      ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
      !-------------------------------------------------------------------------

   ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)

      ! Compute auxiliary variables
      ru2i = ui*DdataAtEdge(2,1,idx)
      rv2i = vi*DdataAtEdge(3,1,idx)
      ru2j = uj*DdataAtEdge(2,2,idx)
      rv2j = vj*DdataAtEdge(3,2,idx)

#ifdef USE_EULER_INTEGRATEBYPARTS
      ! Compute fluxes for x-direction
      dF1_i(1) = DdataAtEdge(2,1,idx)
      dF1_i(2) = G1*DdataAtEdge(4,1,idx)-G14*ru2i-G2*rv2i
      dF1_i(3) = DdataAtEdge(3,1,idx)*ui
      dF1_i(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*ui
      
      dF1_j(1) = DdataAtEdge(2,2,idx)
      dF1_j(2) = G1*DdataAtEdge(4,2,idx)-G14*ru2j-G2*rv2j
      dF1_j(3) = DdataAtEdge(3,2,idx)*uj
      dF1_j(4) = (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*uj
      
      ! Compute fluxes for y-direction
      dF2_i(1) = DdataAtEdge(3,1,idx)
      dF2_i(2) = DdataAtEdge(3,1,idx)*ui
      dF2_i(3) = G1*DdataAtEdge(4,1,idx)-G14*rv2i-G2*ru2i
      dF2_i(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*vi
      
      dF2_j(1) = DdataAtEdge(3,2,idx)
      dF2_j(2) = DdataAtEdge(3,2,idx)*uj
      dF2_j(3) = (G1*DdataAtEdge(4,2,idx)-G14*rv2j-G2*ru2j)
      dF2_j(4) = (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*vj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = DdataAtEdge(2,1,idx) - DdataAtEdge(2,2,idx)
      dF1_ij(2) = (G1*DdataAtEdge(4,1,idx)-G14*ru2i-G2*rv2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*ru2j-G2*rv2j)
      dF1_ij(3) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF1_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*ui-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*uj
      
      ! Compute flux difference for y-direction
      dF2_ij(1) = DdataAtEdge(3,1,idx) - DdataAtEdge(3,2,idx)
      dF2_ij(2) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF2_ij(3) = (G1*DdataAtEdge(4,1,idx)-G14*rv2i-G2*ru2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*rv2j-G2*ru2j)
      dF2_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*vi-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*vj
#endif
      
      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spectral radius
      !-------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient
      a = 0.5_DP*(DmatrixCoeffsAtEdge(:,1,idx)-&
                  DmatrixCoeffsAtEdge(:,2,idx))
      
      ! Compute Roe mean values
      aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
      u_ij = (aux*ui+uj)/(aux+1.0_DP)
      v_ij = (aux*vi+vj)/(aux+1.0_DP)
      hi   = GAMMA*DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui+vi*vi)
      hj   = GAMMA*DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj+vj*vj)
      H_ij = (aux*hi+hj)/(aux+1.0_DP)
      
      ! Compute auxiliary variable
      q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
      aux  = sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL))
      
      ! Scalar dissipation for x- and y-direction
      d_ij = ( abs(a(1)*u_ij) + abs(a(1))*aux +&
               abs(a(2)*v_ij) + abs(a(2))*aux )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef USE_EULER_INTEGRATEBYPARTS
      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*dF2_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij + Diff)
#endif

    end do

  end subroutine euler_calcFluxScDissDiSp2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxRoeDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using tensorial dissipation.
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
#ifdef USE_EULER_INTEGRATEBYPARTS
    real(DP), dimension(NVAR2D) :: dF1_i, dF2_i, dF1_j, dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j
    real(DP) :: aux,aux1,aux2,uPow2,vPow2,hi,hj,H_ij,q_ij,u_ij,v_ij
    real(DP) :: anorm,l1,l2,l3,l4,w1,w2,w3,w4,cPow2,c_ij
    integer :: idx


    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      ! For a detailed description of algorithm and the definition of auxiliary
      ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)

      ! Compute auxiliary variables
      ru2i = ui*DdataAtEdge(2,1,idx)
      rv2i = vi*DdataAtEdge(3,1,idx)
      ru2j = uj*DdataAtEdge(2,2,idx)
      rv2j = vj*DdataAtEdge(3,2,idx)

#ifdef USE_EULER_INTEGRATEBYPARTS
      ! Compute fluxes for x-direction
      dF1_i(1) = DdataAtEdge(2,1,idx)
      dF1_i(2) = G1*DdataAtEdge(4,1,idx)-G14*ru2i-G2*rv2i
      dF1_i(3) = DdataAtEdge(3,1,idx)*ui
      dF1_i(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*ui
      
      dF1_j(1) = DdataAtEdge(2,2,idx)
      dF1_j(2) = G1*DdataAtEdge(4,2,idx)-G14*ru2j-G2*rv2j
      dF1_j(3) = DdataAtEdge(3,2,idx)*uj
      dF1_j(4) = (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*uj
      
      ! Compute fluxes for y-direction
      dF2_i(1) = DdataAtEdge(3,1,idx)
      dF2_i(2) = DdataAtEdge(3,1,idx)*ui
      dF2_i(3) = G1*DdataAtEdge(4,1,idx)-G14*rv2i-G2*ru2i
      dF2_i(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*vi
      
      dF2_j(1) = DdataAtEdge(3,2,idx)
      dF2_j(2) = DdataAtEdge(3,2,idx)*uj
      dF2_j(3) = (G1*DdataAtEdge(4,2,idx)-G14*rv2j-G2*ru2j)
      dF2_j(4) = (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*vj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = DdataAtEdge(2,1,idx) - DdataAtEdge(2,2,idx)
      dF1_ij(2) = (G1*DdataAtEdge(4,1,idx)-G14*ru2i-G2*rv2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*ru2j-G2*rv2j)
      dF1_ij(3) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF1_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*ui-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*uj
      
      ! Compute flux difference for y-direction
      dF2_ij(1) = DdataAtEdge(3,1,idx) - DdataAtEdge(3,2,idx)
      dF2_ij(2) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF2_ij(3) = (G1*DdataAtEdge(4,1,idx)-G14*rv2i-G2*ru2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*rv2j-G2*ru2j)
      dF2_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*vi-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*vj
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor by Roe
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient
      a = 0.5_DP*(DmatrixCoeffsAtEdge(:,1,idx)-&
                  DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      if (anorm .gt. SYS_EPSREAL) then
        
        ! Normalize the skew-symmetric coefficient
        a = a/anorm
        
        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        v_ij = (aux*vi+vj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui+vi*vi)
        hj   = GAMMA*DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj+vj*vj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)

        ! Compute auxiliary variables
        aux   = u_ij*a(1) + v_ij*a(2)
        uPow2 = u_ij*u_ij
        vPow2 = v_ij*v_ij
        q_ij  = 0.5_DP*(uPow2+vPow2)
        cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
        c_ij= sqrt(cPow2)
        
        ! Compute eigenvalues
        l1 = abs(aux-c_ij)
        l2 = abs(aux)
        l3 = abs(aux+c_ij)
        l4 = abs(aux)
        
        ! Compute solution difference U_j-U_i
        Diff = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = G2/cPow2*(q_ij*Diff(1)-u_ij*Diff(2)-v_ij*Diff(3)+Diff(4))
        aux2 = 0.5_DP*(aux*Diff(1)-a(1)*Diff(2)-a(2)*Diff(3))/c_ij
        
        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((1.0_DP-G1*q_ij/cPow2)*Diff(1)+G1*(u_ij*Diff(2)+v_ij*Diff(3)-Diff(4))/cPow2)
        w3 = l3 * (aux1 - aux2)
        w4 = l4 * ((a(1)*v_ij-a(2)*u_ij)*Diff(1)+a(2)*Diff(2)-a(1)*Diff(3))
        
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        Diff(1) = anorm * ( w1 + w2 + w3 )
        Diff(2) = anorm * ( (u_ij-c_ij*a(1))*w1 + u_ij*w2 + (u_ij+c_ij*a(1))*w3 + a(2)*w4 )
        Diff(3) = anorm * ( (v_ij-c_ij*a(2))*w1 + v_ij*w2 + (v_ij+c_ij*a(2))*w3 - a(1)*w4 )
        Diff(4) = anorm * ( (H_ij-c_ij*aux)*w1  + q_ij*w2 +&
                            (H_ij+c_ij*aux)*w3  + (u_ij*a(2)-v_ij*a(1))*w4 )

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef USE_EULER_INTEGRATEBYPARTS
        ! Assemble skew-symmetric fluxes
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                             DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                             DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                             DmatrixCoeffsAtEdge(2,1,idx)*dF2_i + Diff)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        ! Assemble fluxes
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                              DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij + Diff)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                              DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij + Diff)
#endif

      else

#ifdef USE_EULER_INTEGRATEBYPARTS
        ! Assemble skew-symmetric fluxes
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                             DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                             DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                             DmatrixCoeffsAtEdge(2,1,idx)*dF2_i)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        ! Assemble fluxes
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                              DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                              DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij)
#endif

      end if
      
    end do

  end subroutine euler_calcFluxRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxRoeDissDiSp2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using tensorial dissipation,
    ! whereby dimensional splitting is employed.
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
#ifdef USE_EULER_INTEGRATEBYPARTS
    real(DP), dimension(NVAR2D) :: dF1_i, dF2_i, dF1_j, dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff1, Diff2
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j
    real(DP) :: aux,aux1,aux2,uPow2,vPow2,hi,hj,H_ij,q_ij,u_ij,v_ij
    real(DP) :: anorm,l1,l2,l3,l4,w1,w2,w3,w4,cPow2,c_ij
    integer :: idx

    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      ! For a detailed description of algorithm and the definition of auxiliary
      ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)

      ! Compute auxiliary variables
      ru2i = ui*DdataAtEdge(2,1,idx)
      rv2i = vi*DdataAtEdge(3,1,idx)
      ru2j = uj*DdataAtEdge(2,2,idx)
      rv2j = vj*DdataAtEdge(3,2,idx)

#ifdef USE_EULER_INTEGRATEBYPARTS
      ! Compute fluxes for x-direction
      dF1_i(1) = DdataAtEdge(2,1,idx)
      dF1_i(2) = G1*DdataAtEdge(4,1,idx)-G14*ru2i-G2*rv2i
      dF1_i(3) = DdataAtEdge(3,1,idx)*ui
      dF1_i(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*ui
      
      dF1_j(1) = DdataAtEdge(2,2,idx)
      dF1_j(2) = G1*DdataAtEdge(4,2,idx)-G14*ru2j-G2*rv2j
      dF1_j(3) = DdataAtEdge(3,2,idx)*uj
      dF1_j(4) = (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*uj
      
      ! Compute fluxes for y-direction
      dF2_i(1) = DdataAtEdge(3,1,idx)
      dF2_i(2) = DdataAtEdge(3,1,idx)*ui
      dF2_i(3) = G1*DdataAtEdge(4,1,idx)-G14*rv2i-G2*ru2i
      dF2_i(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*vi
      
      dF2_j(1) = DdataAtEdge(3,2,idx)
      dF2_j(2) = DdataAtEdge(3,2,idx)*uj
      dF2_j(3) = (G1*DdataAtEdge(4,2,idx)-G14*rv2j-G2*ru2j)
      dF2_j(4) = (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*vj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = DdataAtEdge(2,1,idx) - DdataAtEdge(2,2,idx)
      dF1_ij(2) = (G1*DdataAtEdge(4,1,idx)-G14*ru2i-G2*rv2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*ru2j-G2*rv2j)
      dF1_ij(3) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF1_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*ui-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*uj
      
      ! Compute flux difference for y-direction
      dF2_ij(1) = DdataAtEdge(3,1,idx) - DdataAtEdge(3,2,idx)
      dF2_ij(2) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF2_ij(3) = (G1*DdataAtEdge(4,1,idx)-G14*rv2i-G2*ru2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*rv2j-G2*ru2j)
      dF2_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*vi-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*vj
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor by Roe
      !-------------------------------------------------------------------------

      ! Compute the skew-symmetric coefficient
      a = 0.5_DP*(DmatrixCoeffsAtEdge(:,1,idx)-&
                  DmatrixCoeffsAtEdge(:,2,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))

      if (anorm .gt. SYS_EPSREAL) then

        ! Compute the absolute value
        a = abs(a)
        
        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        v_ij = (aux*vi+vj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui+vi*vi)
        hj   = GAMMA*DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj+vj*vj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)

        ! Compute auxiliary variables
        aux   = u_ij*a(1) + v_ij*a(2)
        uPow2 = u_ij*u_ij
        vPow2 = v_ij*v_ij
        q_ij  = 0.5_DP*(uPow2+vPow2)
        cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
        c_ij  = sqrt(cPow2)

        !-----------------------------------------------------------------------
        ! Dimensional splitting: x-direction
        !-----------------------------------------------------------------------
        
        ! Compute eigenvalues
        l1 = abs(u_ij-c_ij)
        l2 = abs(u_ij)
        l3 = abs(u_ij+c_ij)
        l4 = abs(u_ij)
        
        ! Compute solution difference U_j-U_i
        Diff1 = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = G2/cPow2*(q_ij*Diff1(1)-u_ij*Diff1(2)-v_ij*Diff1(3)+Diff1(4))
        aux2 = 0.5_DP*(u_ij*Diff1(1)-Diff1(2))/c_ij
        
        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((1.0_DP-G1*q_ij/cPow2)*Diff1(1)+G1*(u_ij*Diff1(2)+v_ij*Diff1(3)-Diff1(4))/cPow2)
        w3 = l3 * (aux1 - aux2)
        w4 = l4 * (v_ij*Diff1(1)-Diff1(3))
        
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        Diff1(1) = a(1) * ( w1 + w2 + w3 )
        Diff1(2) = a(1) * ( (u_ij-c_ij)*w1 + u_ij*w2 + (u_ij+c_ij)*w3 )
        Diff1(3) = a(1) * ( v_ij*w1 + v_ij*w2 + v_ij*w3 - w4 )
        Diff1(4) = a(1) * ( (H_ij-c_ij*u_ij)*w1  + q_ij*w2 + (H_ij+c_ij*u_ij)*w3  -v_ij*w4 )
        
        !-----------------------------------------------------------------------
        ! Dimensional splitting: y-direction
        !-----------------------------------------------------------------------

        ! Compute eigenvalues
        l1 = abs(v_ij-c_ij)
        l2 = abs(v_ij)
        l3 = abs(v_ij+c_ij)
        l4 = abs(v_ij)
        
        ! Compute solution difference U_j-U_i
        Diff2 = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = G2/cPow2*(q_ij*Diff2(1)-u_ij*Diff2(2)-v_ij*Diff2(3)+Diff2(4))
        aux2 = 0.5_DP*(v_ij*Diff2(1)-Diff2(3))/c_ij

        ! Compute characteristic variables multiplied by the corresponding eigenvalue
        w1 = l1 * (aux1 + aux2)
        w2 = l2 * ((1.0_DP-G1*q_ij/cPow2)*Diff2(1)+G1*(u_ij*Diff2(2)+v_ij*Diff2(3)-Diff2(4))/cPow2)
        w3 = l3 * (aux1 - aux2)
        w4 = l4 * (-u_ij*Diff2(1)+Diff2(2))
        
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        Diff2(1) = a(2) * ( w1 + w2 + w3 )
        Diff2(2) = a(2) * ( u_ij*w1 + u_ij*w2 + u_ij*w3 + w4 )
        Diff2(3) = a(2) * ( (v_ij-c_ij)*w1 + v_ij*w2 + (v_ij+c_ij)*w3 )
        Diff2(4) = a(2) * ( (H_ij-c_ij*v_ij)*w1  + q_ij*w2 + (H_ij+c_ij*v_ij)*w3  + u_ij*w4 )

        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-----------------------------------------------------------------------

#ifdef USE_EULER_INTEGRATEBYPARTS
        ! Assemble skew-symmetric fluxes
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                             DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                             DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                             DmatrixCoeffsAtEdge(2,1,idx)*dF2_i+&
                                             Diff1+Diff2)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        ! Assemble fluxes
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                              DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij+&
                                              Diff1+Diff2)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                              DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij+&
                                              Diff1+Diff2)
#endif

      else

#ifdef USE_EULER_INTEGRATEBYPARTS
        ! Assemble skew-symmetric fluxes
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                             DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                             DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                             DmatrixCoeffsAtEdge(2,1,idx)*dF2_i)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
        ! Assemble fluxes
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                              DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                              DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij)
#endif

      end if

    end do

  end subroutine euler_calcFluxRoeDissDiSp2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxRusDiss2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using the Rusanov dissipation.
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
#ifdef USE_EULER_INTEGRATEBYPARTS
    real(DP), dimension(NVAR2D) :: dF1_i, dF2_i, dF1_j, dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j
    real(DP) :: d_ij,ci,cj,Ei,Ej
    integer :: idx

    do idx = 1, size(DfluxesAtEdge,3)

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      ! For a detailed description of algorithm and the definition of auxiliary
      ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      Ei = DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      Ej = DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)

      ! Compute auxiliary variables
      ru2i = ui*DdataAtEdge(2,1,idx)
      rv2i = vi*DdataAtEdge(3,1,idx)
      ru2j = uj*DdataAtEdge(2,2,idx)
      rv2j = vj*DdataAtEdge(3,2,idx)

#ifdef USE_EULER_INTEGRATEBYPARTS
      ! Compute fluxes for x-direction
      dF1_i(1) = DdataAtEdge(2,1,idx)
      dF1_i(2) = G1*DdataAtEdge(4,1,idx)-G14*ru2i-G2*rv2i
      dF1_i(3) = DdataAtEdge(3,1,idx)*ui
      dF1_i(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*ui
      
      dF1_j(1) = DdataAtEdge(2,2,idx)
      dF1_j(2) = G1*DdataAtEdge(4,2,idx)-G14*ru2j-G2*rv2j
      dF1_j(3) = DdataAtEdge(3,2,idx)*uj
      dF1_j(4) = (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*uj
      
      ! Compute fluxes for y-direction
      dF2_i(1) = DdataAtEdge(3,1,idx)
      dF2_i(2) = DdataAtEdge(3,1,idx)*ui
      dF2_i(3) = G1*DdataAtEdge(4,1,idx)-G14*rv2i-G2*ru2i
      dF2_i(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*vi
      
      dF2_j(1) = DdataAtEdge(3,2,idx)
      dF2_j(2) = DdataAtEdge(3,2,idx)*uj
      dF2_j(3) = (G1*DdataAtEdge(4,2,idx)-G14*rv2j-G2*ru2j)
      dF2_j(4) = (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*vj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = DdataAtEdge(2,1,idx) - DdataAtEdge(2,2,idx)
      dF1_ij(2) = (G1*DdataAtEdge(4,1,idx)-G14*ru2i-G2*rv2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*ru2j-G2*rv2j)
      dF1_ij(3) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF1_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*ui-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*uj
      
      ! Compute flux difference for y-direction
      dF2_ij(1) = DdataAtEdge(3,1,idx) - DdataAtEdge(3,2,idx)
      dF2_ij(2) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF2_ij(3) = (G1*DdataAtEdge(4,1,idx)-G14*rv2i-G2*ru2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*rv2j-G2*ru2j)
      dF2_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*vi-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*vj
#endif

      !---------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov type
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
      ci = sqrt(max(G15*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
      cj = sqrt(max(G15*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))
      
      ! Scalar dissipation for the Rusanov flux
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj+&
                      DmatrixCoeffsAtEdge(2,1,idx)*vj)+&
                 sqrt(DmatrixCoeffsAtEdge(1,1,idx)**2+&
                      DmatrixCoeffsAtEdge(2,1,idx)**2)*cj,&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui+&
                      DmatrixCoeffsAtEdge(2,2,idx)*vi)+&
                 sqrt(DmatrixCoeffsAtEdge(1,2,idx)**2+&
                      DmatrixCoeffsAtEdge(2,2,idx)**2)*ci )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef USE_EULER_INTEGRATEBYPARTS
      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*dF2_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij + Diff)
#endif

    end do

  end subroutine euler_calcFluxRusDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxRusDissDiSp2d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using the Rusanov dissipation,
    ! whereby dimensional splitting is employed.
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
#ifdef USE_EULER_INTEGRATEBYPARTS
    real(DP), dimension(NVAR2D) :: dF1_i, dF2_i, dF1_j, dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j
    real(DP) :: d_ij,ci,cj,Ei,Ej
    integer :: idx
    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      ! For a detailed description of algorithm and the definition of auxiliary
      ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      Ei = DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      Ej = DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)

      ! Compute auxiliary variables
      ru2i = ui*DdataAtEdge(2,1,idx)
      rv2i = vi*DdataAtEdge(3,1,idx)
      ru2j = uj*DdataAtEdge(2,2,idx)
      rv2j = vj*DdataAtEdge(3,2,idx)

#ifdef USE_EULER_INTEGRATEBYPARTS
      ! Compute fluxes for x-direction
      dF1_i(1) = DdataAtEdge(2,1,idx)
      dF1_i(2) = G1*DdataAtEdge(4,1,idx)-G14*ru2i-G2*rv2i
      dF1_i(3) = DdataAtEdge(3,1,idx)*ui
      dF1_i(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*ui
      
      dF1_j(1) = DdataAtEdge(2,2,idx)
      dF1_j(2) = G1*DdataAtEdge(4,2,idx)-G14*ru2j-G2*rv2j
      dF1_j(3) = DdataAtEdge(3,2,idx)*uj
      dF1_j(4) = (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*uj
      
      ! Compute fluxes for y-direction
      dF2_i(1) = DdataAtEdge(3,1,idx)
      dF2_i(2) = DdataAtEdge(3,1,idx)*ui
      dF2_i(3) = G1*DdataAtEdge(4,1,idx)-G14*rv2i-G2*ru2i
      dF2_i(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*vi
      
      dF2_j(1) = DdataAtEdge(3,2,idx)
      dF2_j(2) = DdataAtEdge(3,2,idx)*uj
      dF2_j(3) = (G1*DdataAtEdge(4,2,idx)-G14*rv2j-G2*ru2j)
      dF2_j(4) = (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*vj
#else
      ! Compute flux difference for x-direction
      dF1_ij(1) = DdataAtEdge(2,1,idx) - DdataAtEdge(2,2,idx)
      dF1_ij(2) = (G1*DdataAtEdge(4,1,idx)-G14*ru2i-G2*rv2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*ru2j-G2*rv2j)
      dF1_ij(3) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF1_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*ui-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*uj
      
      ! Compute flux difference for y-direction
      dF2_ij(1) = DdataAtEdge(3,1,idx) - DdataAtEdge(3,2,idx)
      dF2_ij(2) = DdataAtEdge(3,1,idx)*ui - DdataAtEdge(3,2,idx)*uj
      dF2_ij(3) = (G1*DdataAtEdge(4,1,idx)-G14*rv2i-G2*ru2i)-&
                  (G1*DdataAtEdge(4,2,idx)-G14*rv2j-G2*ru2j)
      dF2_ij(4) = (GAMMA*DdataAtEdge(4,1,idx)-G2*(ru2i+rv2i))*vi-&
                  (GAMMA*DdataAtEdge(4,2,idx)-G2*(ru2j+rv2j))*vj
#endif
      
      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov type
      !-------------------------------------------------------------------------

      ! Compute the speed of sound
      ci = sqrt(max(G15*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
      cj = sqrt(max(G15*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))
      
      ! Scalar dissipation
      d_ij = max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj)+&
                  abs(DmatrixCoeffsAtEdge(1,1,idx))*cj,&
                  abs(DmatrixCoeffsAtEdge(1,2,idx)*ui)+&
                  abs(DmatrixCoeffsAtEdge(1,2,idx))*ci )&
           + max( abs(DmatrixCoeffsAtEdge(2,1,idx)*vj)+&
                  abs(DmatrixCoeffsAtEdge(2,1,idx))*cj,&
                  abs(DmatrixCoeffsAtEdge(2,2,idx)*vi)+&
                  abs(DmatrixCoeffsAtEdge(2,2,idx))*ci )

      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))
      
      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef USE_EULER_INTEGRATEBYPARTS
      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_j+&
                                           DmatrixCoeffsAtEdge(2,2,idx)*dF2_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF1_i-&
                                           DmatrixCoeffsAtEdge(2,1,idx)*dF2_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,1,idx)*dF2_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF1_ij+&
                                            DmatrixCoeffsAtEdge(2,2,idx)*dF2_ij + Diff)
#endif

    end do

  end subroutine euler_calcFluxRusDissDiSp2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatDiagMatD2d_sim(DdataAtNode, DmatrixCoeffsAtNode,&
      IverticesAtNode, dscale, DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 2D
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
    real(DP) :: ui,vi
    integer :: inode

    do inode = 1, size(DcoefficientsAtNode,3)
      
      ! Compute auxiliary variables
      ui = DdataAtNode(2,inode)/DdataAtNode(1,inode)
      vi = DdataAtNode(3,inode)/DdataAtNode(1,inode)
      
      ! Compute Galerkin coefficient K_ii
      DcoefficientsAtNode(1,1,inode) = 0.0_DP
      DcoefficientsAtNode(2,1,inode) = dscale * (G13*ui*DmatrixCoeffsAtNode(1,inode)+&
                                                     vi*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode(3,1,inode) = dscale * (ui*DmatrixCoeffsAtNode(1,inode)+&
                                                 G13*vi*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode(4,1,inode) = dscale * (GAMMA*(ui*DmatrixCoeffsAtNode(1,inode)+&
                                                        vi*DmatrixCoeffsAtNode(2,inode)))
    end do

  end subroutine euler_calcMatDiagMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatDiag2d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 2D
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
    real(DP) :: ui,vi,qi,Ei,uvi,uPow2i,vPow2i,aux
    integer :: inode

    do inode = 1, size(DcoefficientsAtNode,3)
      
      ! Compute auxiliary variables
      ui = DdataAtNode(2,inode)/DdataAtNode(1,inode)
      vi = DdataAtNode(3,inode)/DdataAtNode(1,inode)
      Ei = DdataAtNode(4,inode)/DdataAtNode(1,inode)
      uvi = ui*vi; qi = ui*ui+vi*vi; uPow2i = ui*ui; vPow2i = vi*vi
      aux = ui*DmatrixCoeffsAtNode(1,inode)+vi*DmatrixCoeffsAtNode(2,inode)
      
      ! Compute Galerkin coefficient K_ii
      DcoefficientsAtNode( 1,1,inode) = 0.0_DP
      DcoefficientsAtNode( 2,1,inode) = dscale * ((G2*qi-uPow2i)*DmatrixCoeffsAtNode(1,inode)-&
                                                  uvi*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode( 3,1,inode) = dscale * ((G2*qi-vPow2i)*DmatrixCoeffsAtNode(2,inode)-&
                                                  uvi*DmatrixCoeffsAtNode(1,inode))
      DcoefficientsAtNode( 4,1,inode) = dscale * (G1*qi-GAMMA*Ei)*aux
      
      DcoefficientsAtNode( 5,1,inode) = dscale * DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode( 6,1,inode) = dscale * (G13*ui*DmatrixCoeffsAtNode(1,inode)+&
                                                  vi*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode( 7,1,inode) = dscale * (vi*DmatrixCoeffsAtNode(1,inode)-&
                                                  G1*ui*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode( 8,1,inode) = dscale * ((GAMMA*Ei-G2*qi)*DmatrixCoeffsAtNode(1,inode)-G1*ui*aux)
      
      DcoefficientsAtNode( 9,1,inode) = dscale * DmatrixCoeffsAtNode(2,inode)
      DcoefficientsAtNode(10,1,inode) = dscale * (ui*DmatrixCoeffsAtNode(2,inode)-&
                                                  G1*vi*DmatrixCoeffsAtNode(1,inode))
      DcoefficientsAtNode(11,1,inode) = dscale * (ui*DmatrixCoeffsAtNode(1,inode)+&
                                                  G13*vi*DmatrixCoeffsAtNode(2,inode))
      DcoefficientsAtNode(12,1,inode) = dscale * ((GAMMA*Ei-G2*qi)*DmatrixCoeffsAtNode(2,inode)-G1*vi*aux)
      
      DcoefficientsAtNode(13,1,inode) = 0.0_DP
      DcoefficientsAtNode(14,1,inode) = dscale * G1*DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(15,1,inode) = dscale * G1*DmatrixCoeffsAtNode(2,inode)
      DcoefficientsAtNode(16,1,inode) = dscale * (GAMMA*(ui*DmatrixCoeffsAtNode(1,inode)+&
                                                  vi*DmatrixCoeffsAtNode(2,inode)))
    end do

  end subroutine euler_calcMatDiag2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatGalMatD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 2D
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
    real(DP) :: ui,uj,vi,vj
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)

      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0_DP
      
      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge(1,2,idx) = 0.0_DP
      DcoefficientsAtEdge(2,2,idx) = dscale * (G13*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                 vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(3,2,idx) = dscale * (uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                 G13*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(4,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                 vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge(1,3,idx) = 0.0_DP
      DcoefficientsAtEdge(2,3,idx) = dscale * (G13*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                 vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(3,3,idx) = dscale * (ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                 G13*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(4,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                 vi*DmatrixCoeffsAtEdge(2,2,idx)))
    end do

  end subroutine euler_calcMatGalMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatGal2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 2D
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
    real(DP) :: Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj,uPow2i,uPow2j,vPow2i,vPow2j,aux1,aux2
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      Ei = DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)
      uvi = ui*vi; qi = ui*ui+vi*vi; uPow2i = ui*ui; vPow2i = vi*vi
      
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      Ej = DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)
      uvj = uj*vj; qj = uj*uj+vj*vj; uPow2j = uj*uj; vPow2j = vj*vj
      
      aux1 = uj*DmatrixCoeffsAtEdge(1,1,idx)+vj*DmatrixCoeffsAtEdge(2,1,idx)
      aux2 = ui*DmatrixCoeffsAtEdge(1,2,idx)+vi*DmatrixCoeffsAtEdge(2,2,idx)
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0_DP

      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge( 1,2,idx) = 0.0_DP
      DcoefficientsAtEdge( 2,2,idx) = dscale * ((G2*qj-uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                  uvj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 3,2,idx) = dscale * ((G2*qj-vPow2j)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                  uvj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge( 4,2,idx) = dscale * (G1*qj-GAMMA*Ej)*aux1
      
      DcoefficientsAtEdge( 5,2,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge( 6,2,idx) = dscale * (G13*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                  vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 7,2,idx) = dscale * (vj*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                  G1*uj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 8,2,idx) = dscale * ((GAMMA*Ej-G2*qj)*DmatrixCoeffsAtEdge(1,1,idx)-G1*uj*aux1)
      
      DcoefficientsAtEdge( 9,2,idx) = dscale * DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(10,2,idx) = dscale * (uj*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                  G1*vj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge(11,2,idx) = dscale * (uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                  G13*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(12,2,idx) = dscale * ((GAMMA*Ej-G2*qj)*DmatrixCoeffsAtEdge(2,1,idx)-G1*vj*aux1)
      
      DcoefficientsAtEdge(13,2,idx) = 0.0_DP
      DcoefficientsAtEdge(14,2,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(15,2,idx) = dscale * G1*DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(16,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                  vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge( 1,3,idx) = 0.0_DP
      DcoefficientsAtEdge( 2,3,idx) = dscale * ((G1*qi-uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                  uvi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 3,3,idx) = dscale * ((G1*qi-vPow2i)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                  uvi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge( 4,3,idx) = dscale * (G1*qi-GAMMA*Ei)*aux2
      
      DcoefficientsAtEdge( 5,3,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge( 6,3,idx) = dscale * (G13*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                  vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 7,3,idx) = dscale * (vi*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                  G1*ui*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 8,3,idx) = dscale * ((GAMMA*Ei-G2*qi)*DmatrixCoeffsAtEdge(1,2,idx)-G1*ui*aux2)
      
      DcoefficientsAtEdge( 9,3,idx) = dscale * DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(10,3,idx) = dscale * (ui*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                  G1*vi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge(11,3,idx) = dscale * (ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                  G13*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(12,3,idx) = dscale * ((GAMMA*Ei-G2*qi)*DmatrixCoeffsAtEdge(2,2,idx)-G1*vi*aux2)
      
      DcoefficientsAtEdge(13,3,idx) = 0.0_DP
      DcoefficientsAtEdge(14,3,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(15,3,idx) = dscale * G1*DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(16,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                  vi*DmatrixCoeffsAtEdge(2,2,idx)))
    end do
      
  end subroutine euler_calcMatGal2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatScDissMatD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies scalar artificial viscosities in 2D
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
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: anorm,aux,hi,hj,H_ij,q_ij,ui,uj,u_ij,vi,vj,v_ij
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      
      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge(1,2,idx) = 0.0_DP
      DcoefficientsAtEdge(2,2,idx) = dscale * (G13*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                 vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(3,2,idx) = dscale * (uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                 G13*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(4,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                 vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge(1,3,idx) = 0.0_DP
      DcoefficientsAtEdge(2,3,idx) = dscale * (G13*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                 vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(3,3,idx) = dscale * (ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                 G13*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(4,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                 vi*DmatrixCoeffsAtEdge(2,2,idx)))

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient and its norm
      a = 0.5_DP*(DmatrixCoeffsAtEdge(:,2,idx)-&
                  DmatrixCoeffsAtEdge(:,1,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        v_ij = (aux*vi+vj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui+vi*vi)
        hj   = GAMMA*DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj+vj*vj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)
        
        ! Compute auxiliary variables
        q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
        
        ! Compute scalar dissipation
        DcoefficientsAtEdge(:,1,idx) = -dscale * (abs(a(1)*u_ij+a(2)*v_ij) +&
            anorm*sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL)))
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0_DP

      end if
    end do

  end subroutine euler_calcMatScDissMatD2d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatScDiss2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies scalar artificial viscosities in 2D
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
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: anorm,aux,aux1,aux2,hi,hj,H_ij,q_ij,u_ij,v_ij
    real(DP) :: Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj,uPow2i,uPow2j,vPow2i,vPow2j
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      Ei = DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)
      uvi = ui*vi; qi = ui*ui+vi*vi; uPow2i = ui*ui; vPow2i = vi*vi
      
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      Ej = DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)
      uvj = uj*vj; qj = uj*uj+vj*vj; uPow2j = uj*uj; vPow2j = vj*vj
      
      aux1 = uj*DmatrixCoeffsAtEdge(1,1,idx)+vj*DmatrixCoeffsAtEdge(2,1,idx)
      aux2 = ui*DmatrixCoeffsAtEdge(1,2,idx)+vi*DmatrixCoeffsAtEdge(2,2,idx)
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0_DP
      
      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge( 1,2,idx) = 0.0_DP
      DcoefficientsAtEdge( 2,2,idx) = dscale * ((G2*qj-uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                  uvj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 3,2,idx) = dscale * ((G2*qj-vPow2j)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                  uvj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge( 4,2,idx) = dscale * (G1*qj-GAMMA*Ej)*aux1
      
      DcoefficientsAtEdge( 5,2,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge( 6,2,idx) = dscale * (G13*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                  vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 7,2,idx) = dscale * (vj*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                  G1*uj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 8,2,idx) = dscale * ((GAMMA*Ej-G2*qj)*DmatrixCoeffsAtEdge(1,1,idx)-G1*uj*aux1)
      
      DcoefficientsAtEdge( 9,2,idx) = dscale * DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(10,2,idx) = dscale * (uj*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                  G1*vj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge(11,2,idx) = dscale * (uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                  G13*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(12,2,idx) = dscale * ((GAMMA*Ej-G2*qj)*DmatrixCoeffsAtEdge(2,1,idx)-G1*vj*aux1)
      
      DcoefficientsAtEdge(13,2,idx) = 0.0_DP
      DcoefficientsAtEdge(14,2,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(15,2,idx) = dscale * G1*DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(16,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                  vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge( 1,3,idx) = 0.0_DP
      DcoefficientsAtEdge( 2,3,idx) = dscale * ((G1*qi-uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                  uvi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 3,3,idx) = dscale * ((G1*qi-vPow2i)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                  uvi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge( 4,3,idx) = dscale * (G1*qi-GAMMA*Ei)*aux2
      
      DcoefficientsAtEdge( 5,3,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge( 6,3,idx) = dscale * (G13*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                  vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 7,3,idx) = dscale * (vi*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                  G1*ui*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 8,3,idx) = dscale * ((GAMMA*Ei-G2*qi)*DmatrixCoeffsAtEdge(1,2,idx)-G1*ui*aux2)
      
      DcoefficientsAtEdge( 9,3,idx) = dscale * DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(10,3,idx) = dscale * (ui*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                  G1*vi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge(11,3,idx) = dscale * (ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                  G13*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(12,3,idx) = dscale * ((GAMMA*Ei-G2*qi)*DmatrixCoeffsAtEdge(2,2,idx)-G1*vi*aux2)
      
      DcoefficientsAtEdge(13,3,idx) = 0.0_DP
      DcoefficientsAtEdge(14,3,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(15,3,idx) = dscale * G1*DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(16,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                  vi*DmatrixCoeffsAtEdge(2,2,idx)))

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient and its norm
      a = 0.5_DP*(DmatrixCoeffsAtEdge(:,2,idx)-&
                  DmatrixCoeffsAtEdge(:,1,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        v_ij = (aux*vi+vj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui+vi*vi)
        hj   = GAMMA*DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj+vj*vj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)
        
        ! Compute auxiliary variables
        q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)

        ! Compute scalar dissipation
        aux = -dscale * (abs(a(1)*u_ij+a(2)*v_ij) +&
            anorm*sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL)))
              
        DcoefficientsAtEdge( 1,1,idx) = aux
        DcoefficientsAtEdge( 6,1,idx) = aux
        DcoefficientsAtEdge(11,1,idx) = aux
        DcoefficientsAtEdge(16,1,idx) = aux

      end if
    end do

  end subroutine euler_calcMatScDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatRoeDissMatD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 2D
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
    real(DP), dimension(NVAR2D,NVAR2D) :: R_ij,L_ij
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: aux,hi,hj,H_ij,q_ij,ui,uj,u_ij,vi,vj,v_ij
    real(DP) :: l1,l2,l3,l4,anorm,c1,c2,c_ij,cPow2,vel
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)

      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0_DP
      
      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge(1,2,idx) = 0.0_DP
      DcoefficientsAtEdge(2,2,idx) = dscale * (G13*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                 vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(3,2,idx) = dscale * (uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                 G13*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(4,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                 vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge(1,3,idx) = 0.0_DP
      DcoefficientsAtEdge(2,3,idx) = dscale * (G13*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                 vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(3,3,idx) = dscale * (ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                 G13*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(4,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                 vi*DmatrixCoeffsAtEdge(2,2,idx)))

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient and its norm
      a = 0.5_DP*(DmatrixCoeffsAtEdge(:,2,idx)-&
                  DmatrixCoeffsAtEdge(:,1,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL) then

        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        v_ij = (aux*vi+vj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui+vi*vi)
        hj   = GAMMA*DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj+vj*vj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)

        ! Compute auxiliary values
        c1    = a(1)/anorm
        c2    = a(2)/anorm
        vel   = c1*u_ij+c2*v_ij
        q_ij  = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
        cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
        c_ij= sqrt(cPow2)
        
        ! Diagonal matrix of eigenvalues
        l1 = abs(vel-c_ij)
        l2 = abs(vel)
        l3 = abs(vel+c_ij)
        l4 = abs(vel)

        ! Matrix of right eigenvectors
        R_ij(1,1) =  l1
        R_ij(2,1) =  l1*(u_ij-c_ij*c1)
        R_ij(3,1) =  l1*(v_ij-c_ij*c2)
        R_ij(4,1) =  l1*(H_ij-c_ij*vel)
        
        R_ij(1,2) =  l2
        R_ij(2,2) =  l2*u_ij
        R_ij(3,2) =  l2*v_ij
        R_ij(4,2) =  l2*q_ij
        
        R_ij(1,3) =  l3
        R_ij(2,3) =  l3*(u_ij+c_ij*c1)
        R_ij(3,3) =  l3*(v_ij+c_ij*c2)
        R_ij(4,3) =  l3*(H_ij+c_ij*vel)
        
        R_ij(1,4) =  0.0_DP
        R_ij(2,4) =  l4*c2
        R_ij(3,4) = -l4*c1
        R_ij(4,4) =  l4*(u_ij*c2-v_ij*c1)
        
        ! Matrix of left eigenvectors
        L_ij(1,1) = 0.5_DP*(G1*q_ij+c_ij*vel)/cPow2
        L_ij(2,1) = (cPow2-G1*q_ij)/cPow2
        L_ij(3,1) = 0.5_DP*(G1*q_ij-c_ij*vel)/cPow2
        L_ij(4,1) = v_ij*c1-u_ij*c2
        
        L_ij(1,2) = 0.5_DP*(-G1*u_ij-c_ij*c1)/cPow2
        L_ij(2,2) = G1*u_ij/cPow2
        L_ij(3,2) = 0.5_DP*(-G1*u_ij+c_ij*c1)/cPow2
        L_ij(4,2) = c2

        L_ij(1,3) = 0.5_DP*(-G1*v_ij-c_ij*c2)/cPow2
        L_ij(2,3) = G1*v_ij/cPow2
        L_ij(3,3) = 0.5_DP*(-G1*v_ij+c_ij*c2)/cPow2
        L_ij(4,3) = -c1
        
        L_ij(1,4) =  G2/cPow2
        L_ij(2,4) = -G1/cPow2
        L_ij(3,4) =  G2/cPow2
        L_ij(4,4) =  0.0_DP
        
        ! Include scaling parameter
        anorm = -dscale*anorm
        
        ! Compute tensorial dissipation D_ij = diag(R_ij*|Lbd_ij|*L_ij)*I
        DcoefficientsAtEdge(1,1,idx) = anorm*( R_ij(1,1)*L_ij(1,1)+&
            R_ij(1,2)*L_ij(2,1)+R_ij(1,3)*L_ij(3,1)+R_ij(1,4)*L_ij(4,1)  )
        DcoefficientsAtEdge(2,1,idx) = anorm*( R_ij(2,1)*L_ij(1,2)+&
            R_ij(2,2)*L_ij(2,2)+R_ij(2,3)*L_ij(3,2)+R_ij(2,4)*L_ij(4,2)  )
        DcoefficientsAtEdge(3,1,idx) = anorm*( R_ij(3,1)*L_ij(1,3)+&
            R_ij(3,2)*L_ij(2,3)+R_ij(3,3)*L_ij(3,3)+R_ij(3,4)*L_ij(4,3)  )
        DcoefficientsAtEdge(4,1,idx) = anorm*( R_ij(4,1)*L_ij(1,4)+&
            R_ij(4,2)*L_ij(2,4)+R_ij(4,3)*L_ij(3,4)+R_ij(4,4)*L_ij(4,4)  )
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0_DP
        
      end if
    end do

  end subroutine euler_calcMatRoeDissMatD2d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatRoeDiss2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 2D
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
    real(DP), dimension(NVAR2D,NVAR2D) :: R_ij,L_ij
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: anorm,aux,hi,hj,H_ij,q_ij,u_ij,v_ij,vel,c1,c2,cPow2,c_ij,l1,l2,l3,l4
    real(DP) :: Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj,uPow2i,uPow2j,vPow2i,vPow2j,aux1,aux2
    integer :: idx,i,j,k

    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      Ei = DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)
      uvi = ui*vi; qi = ui*ui+vi*vi; uPow2i = ui*ui; vPow2i = vi*vi
      
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      Ej = DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)
      uvj = uj*vj; qj = uj*uj+vj*vj; uPow2j = uj*uj; vPow2j = vj*vj
      
      aux1 = uj*DmatrixCoeffsAtEdge(1,1,idx)+vj*DmatrixCoeffsAtEdge(2,1,idx)
      aux2 = ui*DmatrixCoeffsAtEdge(1,2,idx)+vi*DmatrixCoeffsAtEdge(2,2,idx)

      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge( 1,2,idx) = 0.0_DP
      DcoefficientsAtEdge( 2,2,idx) = dscale * ((G2*qj-uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                 uvj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 3,2,idx) = dscale * ((G2*qj-vPow2j)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                 uvj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge( 4,2,idx) = dscale * (G1*qj-GAMMA*Ej)*aux1
      
      DcoefficientsAtEdge( 5,2,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge( 6,2,idx) = dscale * (G13*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 7,2,idx) = dscale * (vj*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                G1*uj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 8,2,idx) = dscale * ((GAMMA*Ej-G2*qj)*DmatrixCoeffsAtEdge(1,1,idx)-G1*uj*aux1)
      
      DcoefficientsAtEdge( 9,2,idx) = dscale * DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(10,2,idx) = dscale * (uj*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                G1*vj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge(11,2,idx) = dscale * (uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                G13*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(12,2,idx) = dscale * ((GAMMA*Ej-G2*qj)*DmatrixCoeffsAtEdge(2,1,idx)-G1*vj*aux1)
      
      DcoefficientsAtEdge(13,2,idx) = 0.0_DP
      DcoefficientsAtEdge(14,2,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(15,2,idx) = dscale * G1*DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(16,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge( 1,3,idx) = 0.0_DP
      DcoefficientsAtEdge( 2,3,idx) = dscale * ((G1*qi-uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                uvi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 3,3,idx) = dscale * ((G1*qi-vPow2i)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                uvi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge( 4,3,idx) = dscale * (G1*qi-GAMMA*Ei)*aux2
      
      DcoefficientsAtEdge( 5,3,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge( 6,3,idx) = dscale * (G13*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 7,3,idx) = dscale * (vi*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                G1*ui*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 8,3,idx) = dscale * ((GAMMA*Ei-G2*qi)*DmatrixCoeffsAtEdge(1,2,idx)-G1*ui*aux2)
      
      DcoefficientsAtEdge( 9,3,idx) = dscale * DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(10,3,idx) = dscale * (ui*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                G1*vi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge(11,3,idx) = dscale * (ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                G13*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(12,3,idx) = dscale * ((GAMMA*Ei-G2*qi)*DmatrixCoeffsAtEdge(2,2,idx)-G1*vi*aux2)
      
      DcoefficientsAtEdge(13,3,idx) = 0.0_DP
      DcoefficientsAtEdge(14,3,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(15,3,idx) = dscale * G1*DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(16,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                vi*DmatrixCoeffsAtEdge(2,2,idx)))

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient and its norm
      a = 0.5_DP*(DmatrixCoeffsAtEdge(:,2,idx)-&
                  DmatrixCoeffsAtEdge(:,1,idx))
      anorm = sqrt(a(1)*a(1)+a(2)*a(2))
      
      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        v_ij = (aux*vi+vj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui+vi*vi)
        hj   = GAMMA*DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj+vj*vj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)
        
        ! Compute auxiliary variables
        c1    = a(1)/anorm
        c2    = a(2)/anorm
        vel   = c1*u_ij+c2*v_ij
        q_ij  = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
        cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
        c_ij  = sqrt(cPow2)
        
        ! Diagonal matrix of eigenvalues
        l1 = abs(vel-c_ij)
        l2 = abs(vel)
        l3 = abs(vel+c_ij)
        l4 = abs(vel)
        
        ! Matrix of right eigenvectors
        R_ij(1,1) =  l1
        R_ij(2,1) =  l1*(u_ij-c_ij*c1)
        R_ij(3,1) =  l1*(v_ij-c_ij*c2)
        R_ij(4,1) =  l1*(H_ij-c_ij*vel)
        
        R_ij(1,2) =  l2
        R_ij(2,2) =  l2*u_ij
        R_ij(3,2) =  l2*v_ij
        R_ij(4,2) =  l2*q_ij
        
        R_ij(1,3) =  l3
        R_ij(2,3) =  l3*(u_ij+c_ij*c1)
        R_ij(3,3) =  l3*(v_ij+c_ij*c2)
        R_ij(4,3) =  l3*(H_ij+c_ij*vel)
        
        R_ij(1,4) =  0.0_DP
        R_ij(2,4) =  l4*c2
        R_ij(3,4) = -l4*c1
        R_ij(4,4) =  l4*(u_ij*c2-v_ij*c1)
        
        ! Matrix of left eigenvectors
        L_ij(1,1) = 0.5_DP*(G1*q_ij+c_ij*vel)/cPow2
        L_ij(2,1) = (cPow2-G1*q_ij)/cPow2
        L_ij(3,1) = 0.5_DP*(G1*q_ij-c_ij*vel)/cPow2
        L_ij(4,1) = v_ij*c1-u_ij*c2
        
        L_ij(1,2) = 0.5_DP*(-G1*u_ij-c_ij*c1)/cPow2
        L_ij(2,2) = G1*u_ij/cPow2
        L_ij(3,2) = 0.5_DP*(-G1*u_ij+c_ij*c1)/cPow2
        L_ij(4,2) = c2
        
        L_ij(1,3) = 0.5_DP*(-G1*v_ij-c_ij*c2)/cPow2
        L_ij(2,3) = G1*v_ij/cPow2
        L_ij(3,3) = 0.5_DP*(-G1*v_ij+c_ij*c2)/cPow2
        L_ij(4,3) = -c1
        
        L_ij(1,4) =  G2/cPow2
        L_ij(2,4) = -G1/cPow2
        L_ij(3,4) =  G2/cPow2
        L_ij(4,4) =  0.0_DP
        
        ! Include scaling parameter
        anorm = -dscale*anorm

        ! Compute tensorial dissipation D_ij = R_ij*|Lbd_ij|*L_ij
        do i = 1, NVAR2D
          do j = 1, NVAR2D
            aux = 0.0_DP
            do k = 1, NVAR2D
              aux = aux + R_ij(i,k)*L_ij(k,j)
            end do
            DcoefficientsAtEdge(NVAR2D*(j-1)+i,1,idx) = anorm*aux
          end do
        end do
        
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0_DP
        
      end if
    end do

  end subroutine euler_calcMatRoeDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatRusDissMatD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 2D
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
    real(DP) :: ui,uj,vi,vj,ci,cj,Ei,Ej
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      Ei = DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      Ej = DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)

      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge(1,2,idx) = 0.0_DP
      DcoefficientsAtEdge(2,2,idx) = dscale * (G13*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                               vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(3,2,idx) = dscale * (uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                               G13*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(4,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                               vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge(1,3,idx) = 0.0_DP
      DcoefficientsAtEdge(2,3,idx) = dscale * (G13*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                               vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(3,3,idx) = dscale * (ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                               G13*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(4,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                               vi*DmatrixCoeffsAtEdge(2,2,idx)))

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
      ci = sqrt(max(G15*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
      cj = sqrt(max(G15*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))
      
      ! Compute dissipation tensor D_ij
      DcoefficientsAtEdge(:,1,idx) = -dscale *&
          max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj+&
                   DmatrixCoeffsAtEdge(2,1,idx)*vj) +&
                   sqrt(DmatrixCoeffsAtEdge(1,1,idx)**2+&
                        DmatrixCoeffsAtEdge(2,1,idx)**2)*cj,&
               abs(DmatrixCoeffsAtEdge(1,2,idx)*ui+&
                   DmatrixCoeffsAtEdge(2,2,idx)*vi) +&
                   sqrt(DmatrixCoeffsAtEdge(1,2,idx)**2+&
                        DmatrixCoeffsAtEdge(2,2,idx)**2)*ci )
    end do

  end subroutine euler_calcMatRusDissMatD2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatRusDiss2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 2D
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
    real(DP) :: ci,cj,Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj
    real(DP) :: uPow2i,uPow2j,vPow2i,vPow2j,aux1,aux2
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      Ei = DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)
      uvi = ui*vi; qi = ui*ui+vi*vi; uPow2i = ui*ui; vPow2i = vi*vi
      
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      Ej = DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)
      uvj = uj*vj; qj = uj*uj+vj*vj; uPow2j = uj*uj; vPow2j = vj*vj
      
      aux1 = uj*DmatrixCoeffsAtEdge(1,1,idx)+vj*DmatrixCoeffsAtEdge(2,1,idx)
      aux2 = ui*DmatrixCoeffsAtEdge(1,2,idx)+vi*DmatrixCoeffsAtEdge(2,2,idx)

      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge( 1,2,idx) = 0.0_DP
      DcoefficientsAtEdge( 2,2,idx) = dscale * ((G2*qj-uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                uvj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 3,2,idx) = dscale * ((G2*qj-vPow2j)*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                uvj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge( 4,2,idx) = dscale * (G1*qj-GAMMA*Ej)*aux1
      
      DcoefficientsAtEdge( 5,2,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge( 6,2,idx) = dscale * (G13*uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 7,2,idx) = dscale * (vj*DmatrixCoeffsAtEdge(1,1,idx)-&
                                                G1*uj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge( 8,2,idx) = dscale * ((GAMMA*Ej-G2*qj)*DmatrixCoeffsAtEdge(1,1,idx)-G1*uj*aux1)
      
      DcoefficientsAtEdge( 9,2,idx) = dscale * DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(10,2,idx) = dscale * (uj*DmatrixCoeffsAtEdge(2,1,idx)-&
                                                G1*vj*DmatrixCoeffsAtEdge(1,1,idx))
      DcoefficientsAtEdge(11,2,idx) = dscale * (uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                G13*vj*DmatrixCoeffsAtEdge(2,1,idx))
      DcoefficientsAtEdge(12,2,idx) = dscale * ((GAMMA*Ej-G2*qj)*DmatrixCoeffsAtEdge(2,1,idx)-G1*vj*aux1)
      
      DcoefficientsAtEdge(13,2,idx) = 0.0_DP
      DcoefficientsAtEdge(14,2,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(15,2,idx) = dscale * G1*DmatrixCoeffsAtEdge(2,1,idx)
      DcoefficientsAtEdge(16,2,idx) = dscale * (GAMMA*(uj*DmatrixCoeffsAtEdge(1,1,idx)+&
                                                vj*DmatrixCoeffsAtEdge(2,1,idx)))
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge( 1,3,idx) = 0.0_DP
      DcoefficientsAtEdge( 2,3,idx) = dscale * ((G1*qi-uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                uvi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 3,3,idx) = dscale * ((G1*qi-vPow2i)*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                uvi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge( 4,3,idx) = dscale * (G1*qi-GAMMA*Ei)*aux2
      
      DcoefficientsAtEdge( 5,3,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge( 6,3,idx) = dscale * (G13*ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 7,3,idx) = dscale * (vi*DmatrixCoeffsAtEdge(1,2,idx)-&
                                                G1*ui*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge( 8,3,idx) = dscale * ((GAMMA*Ei-G2*qi)*DmatrixCoeffsAtEdge(1,2,idx)-G1*ui*aux2)
      
      DcoefficientsAtEdge( 9,3,idx) = dscale * DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(10,3,idx) = dscale * (ui*DmatrixCoeffsAtEdge(2,2,idx)-&
                                                G1*vi*DmatrixCoeffsAtEdge(1,2,idx))
      DcoefficientsAtEdge(11,3,idx) = dscale * (ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                G13*vi*DmatrixCoeffsAtEdge(2,2,idx))
      DcoefficientsAtEdge(12,3,idx) = dscale * ((GAMMA*Ei-G2*qi)*DmatrixCoeffsAtEdge(2,2,idx)-G1*vi*aux2)
      
      DcoefficientsAtEdge(13,3,idx) = 0.0_DP
      DcoefficientsAtEdge(14,3,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(15,3,idx) = dscale * G1*DmatrixCoeffsAtEdge(2,2,idx)
      DcoefficientsAtEdge(16,3,idx) = dscale * (GAMMA*(ui*DmatrixCoeffsAtEdge(1,2,idx)+&
                                                vi*DmatrixCoeffsAtEdge(2,2,idx)))

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------

      ! Compute the speed of sound
      ci = sqrt(max(G15*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
      cj = sqrt(max(G15*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))
      
      ! Compute dissipation tensor D_ij
      aux1 = -dscale *&
          max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj+&
                   DmatrixCoeffsAtEdge(2,1,idx)*vj) +&
                   sqrt(DmatrixCoeffsAtEdge(1,1,idx)**2+&
                        DmatrixCoeffsAtEdge(2,1,idx)**2)*cj,&
               abs(DmatrixCoeffsAtEdge(1,2,idx)*ui+&
                   DmatrixCoeffsAtEdge(2,2,idx)*vi) +&
                   sqrt(DmatrixCoeffsAtEdge(1,2,idx)**2+&
                        DmatrixCoeffsAtEdge(2,2,idx)**2)*ci )

      DcoefficientsAtEdge(:,1,idx) = 0.0_DP
      DcoefficientsAtEdge( 1,1,idx) = aux1
      DcoefficientsAtEdge( 6,1,idx) = aux1
      DcoefficientsAtEdge(11,1,idx) = aux1
      DcoefficientsAtEdge(16,1,idx) = aux1
    end do

  end subroutine euler_calcMatRusDiss2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcCharacteristics2d_sim(Dweight, DdataAtEdge,&
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
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: u_ij,v_ij,H_ij,q_ij,c_ij,aux,aux1,aux2,ui,uj,vi,vj,hi,hj,cPow2,a1,a2,anorm
    integer :: idx

    ! Compute norm of weighting coefficient
    anorm = sqrt(Dweight(1)*Dweight(1)+Dweight(2)*Dweight(2))

    ! Check if weighting coefficient is zero
    if (anorm .le. SYS_EPSREAL) then
      if (present(DcharVariablesAtEdge))     DcharVariablesAtEdge     = 0.0_DP
      if (present(DeigenvaluesAtEdge))       DeigenvaluesAtEdge       = 0.0_DP
      if (present(DrightEigenvectorsAtEdge)) DrightEigenvectorsAtEdge = 0.0_DP
      if (present(DleftEigenvectorsAtEdge))  DleftEigenvectorsAtEdge  = 0.0_DP

      ! That`s it
      return
    end if

    ! Compute normalised weighting coefficient
    a1  = Dweight(1)/anorm
    a2  = Dweight(2)/anorm

    ! Do we have to compute characteristic variables
    if (present(DcharVariablesAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)
        
        ! Compute velocities
        ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
        vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
        uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
        vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
        
        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        v_ij = (aux*vi+vj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui+vi*vi)
        hj   = GAMMA*DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj+vj*vj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)
        
        ! Compute auxiliary variables
        q_ij  = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
        cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
        c_ij  = sqrt(cPow2)
        aux   = a1*u_ij+a2*v_ij

        ! Compute solution difference U_j-U_i
        Diff = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
        
        ! Compute auxiliary quantities for characteristic variables
        aux1 = G2/cPow2*(q_ij*Diff(1)-u_ij*Diff(2)-v_ij*Diff(3)+Diff(4) )
        aux2 = 0.5_DP*(aux*Diff(1)-a1*Diff(2)-a2*Diff(3) )/c_ij

        ! Compute characteristic variables
        DcharVariablesAtEdge(1,idx) = anorm * (aux1 + aux2)
        DcharVariablesAtEdge(2,idx) = anorm * ((1.0_DP-G1*q_ij/cPow2)*Diff(1)+&
                                                             G1*(u_ij*Diff(2)+&
                                                                 v_ij*Diff(3)-&
                                                                      Diff(4))/cPow2)
        DcharVariablesAtEdge(3,idx) = anorm * (aux1 - aux2)
        DcharVariablesAtEdge(4,idx) = anorm * ((a1*v_ij-a2*u_ij)*Diff(1)+&
                                                              a2*Diff(2)-&
                                                              a1*Diff(3))
      end do
    end if

    ! Do we have to compute eigenvalues
    if (present(DeigenvaluesAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)

        ! Compute velocities
        ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
        vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
        uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
        vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
        
        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        v_ij = (aux*vi+vj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui+vi*vi)
        hj   = GAMMA*DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj+vj*vj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)
        
        ! Compute auxiliary variables
        c_ij= sqrt(max(G1*(H_ij-0.5_DP*(u_ij*u_ij+v_ij*v_ij)), SYS_EPSREAL))
        aux   = a1*u_ij+a2*v_ij

        ! Compute eigenvalues
        DeigenvaluesAtEdge(1,idx) = aux-c_ij
        DeigenvaluesAtEdge(2,idx) = aux
        DeigenvaluesAtEdge(3,idx) = aux+c_ij
        DeigenvaluesAtEdge(4,idx) = aux
      end do
    end if

    ! Do we have to compute right eigenvectors
    if (present(DrightEigenvectorsAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)

        ! Compute velocities
        ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
        vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
        uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
        vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
        
        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        v_ij = (aux*vi+vj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui+vi*vi)
        hj   = GAMMA*DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj+vj*vj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)
        
        ! Compute auxiliary variables
        q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
        c_ij  = sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL))
        aux  = a1*u_ij+a2*v_ij

        ! Compute right eigenvectors
        DrightEigenvectorsAtEdge( 1,idx) =  1.0_DP
        DrightEigenvectorsAtEdge( 2,idx) =  u_ij-c_ij*a1
        DrightEigenvectorsAtEdge( 3,idx) =  v_ij-c_ij*a2
        DrightEigenvectorsAtEdge( 4,idx) =  H_ij-c_ij*aux

        DrightEigenvectorsAtEdge( 5,idx) =  1.0_DP
        DrightEigenvectorsAtEdge( 6,idx) =  u_ij
        DrightEigenvectorsAtEdge( 7,idx) =  v_ij
        DrightEigenvectorsAtEdge( 8,idx) =  q_ij

        DrightEigenvectorsAtEdge( 9,idx) =  1.0_DP
        DrightEigenvectorsAtEdge(10,idx) =  u_ij+c_ij*a1
        DrightEigenvectorsAtEdge(11,idx) =  v_ij+c_ij*a2
        DrightEigenvectorsAtEdge(12,idx) =  H_ij+c_ij*aux

        DrightEigenvectorsAtEdge(13,idx) =  0.0_DP
        DrightEigenvectorsAtEdge(14,idx) =  a2
        DrightEigenvectorsAtEdge(15,idx) = -a1
        DrightEigenvectorsAtEdge(16,idx) =  u_ij*a2-v_ij*a1
      end do
    end if

    ! Do we have to compute left eigenvectors
    if (present(DleftEigenvectorsAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)

        ! Compute velocities
        ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
        vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
        uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
        vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
        
        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        v_ij = (aux*vi+vj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui+vi*vi)
        hj   = GAMMA*DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj+vj*vj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)

        ! Compute auxiliary variables
        q_ij  = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
        cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
        c_ij   = sqrt(cPow2)
        aux   = a1*u_ij+a2*v_ij

        ! Compute left eigenvectors
        DleftEigenvectorsAtEdge( 1,idx) =  0.5_DP*(G1*q_ij+c_ij*aux)/cPow2
        DleftEigenvectorsAtEdge( 2,idx) = (cPow2-G1*q_ij)/cPow2
        DleftEigenvectorsAtEdge( 3,idx) =  0.5_DP*(G1*q_ij-c_ij*aux)/cPow2
        DleftEigenvectorsAtEdge( 4,idx) =  v_ij*a1-u_ij*a2

        DleftEigenvectorsAtEdge( 5,idx) =  0.5_DP*(-G1*u_ij-c_ij*a1)/cPow2
        DleftEigenvectorsAtEdge( 6,idx) =  G1*u_ij/cPow2
        DleftEigenvectorsAtEdge( 7,idx) =  0.5_DP*(-G1*u_ij+c_ij*a1)/cPow2
        DleftEigenvectorsAtEdge( 8,idx) =  a2

        DleftEigenvectorsAtEdge( 9,idx) =  0.5_DP*(-G1*v_ij-c_ij*a2)/cPow2
        DleftEigenvectorsAtEdge(10,idx) =  G1*v_ij/cPow2
        DleftEigenvectorsAtEdge(11,idx) =  0.5_DP*(-G1*v_ij+c_ij*a2)/cPow2
        DleftEigenvectorsAtEdge(12,idx) = -a1

        DleftEigenvectorsAtEdge(13,idx) =  G2/cPow2
        DleftEigenvectorsAtEdge(14,idx) = -G1/cPow2
        DleftEigenvectorsAtEdge(15,idx) =  G2/cPow2
        DleftEigenvectorsAtEdge(16,idx) =  0.0_DP
      end do
    end if

  end subroutine euler_calcCharacteristics2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxFCTScalarDiss2d(&
      U1_i, U1_j, U2_i, U2_j, C_ij, C_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 2D using scalar dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

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
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: ui,vi,uj,vj
    real(DP) :: d_ij,hi,hj,H_ij,q_ij,u_ij,v_ij,aux,vel,cs

    ! Compute velocities
    ui = U2_i(2)/U2_i(1); vi = U2_i(3)/U2_i(1)
    uj = U2_j(2)/U2_j(1); vj = U2_j(3)/U2_j(1)

    ! Compute skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji)

    ! Compute Roe mean values
    aux  = sqrt(max(U2_i(1)/U2_j(1), SYS_EPSREAL))
    u_ij = (aux*ui+uj)/(aux+1.0_DP)
    v_ij = (aux*vi+vj)/(aux+1.0_DP)
    hi   = GAMMA*U2_i(4)/U2_i(1)-G2*(U2_i(2)*U2_i(2)+U2_i(3)*U2_i(3))/(U2_i(1)*U2_i(1))
    hj   = GAMMA*U2_j(4)/U2_j(1)-G2*(U2_j(2)*U2_j(2)+U2_j(3)*U2_j(3))/(U2_j(1)*U2_j(1))
    H_ij = (aux*hi+hj)/(aux+1.0_DP)

    ! Compute auxiliary variables
    aux  = sqrt(a(1)*a(1)+a(2)*a(2))
    vel  = u_ij*a(1) + v_ij*a(2)
    q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
    cs   = sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL))

    ! Scalar dissipation
    d_ij = abs(vel) + aux*cs

    ! Compute conservative fluxes
    F_ij = dscale1*(U1_i-U1_j) + dscale2*d_ij*(U2_i-U2_j)

  end subroutine euler_calcFluxFCTScalarDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxFCTTensorDiss2d(&
      U1_i, U1_j, U2_i, U2_j, C_ij, C_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 2D using tensorial dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

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
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: ui,vi,uj,vj
    real(DP) :: aux,aux1,aux2,uPow2,vPow2,hi,hj,H_ij,q_ij,u_ij,v_ij
    real(DP) :: anorm,l1,l2,l3,l4,w1,w2,w3,w4,cPow2,cs

    ! Compute velocities
    ui = U2_i(2)/U2_i(1); vi = U2_i(3)/U2_i(1)
    uj = U2_j(2)/U2_j(1); vj = U2_j(3)/U2_j(1)

    ! Compute the skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji); anorm = sqrt(a(1)*a(1)+a(2)*a(2))

    if (anorm .gt. SYS_EPSREAL) then

      ! Normalize the skew-symmetric coefficient
      a = a/anorm

      ! Compute Roe mean values
      aux  = sqrt(max(U2_i(1)/U2_j(1), SYS_EPSREAL))
      u_ij = (aux*ui+uj)/(aux+1.0_DP)
      v_ij = (aux*vi+vj)/(aux+1.0_DP)
      hi   = GAMMA*U2_i(4)/U2_i(1)-G2*(U2_i(2)*U2_i(2)+U2_i(3)*U2_i(3))/(U2_i(1)*U2_i(1))
      hj   = GAMMA*U2_j(4)/U2_j(1)-G2*(U2_j(2)*U2_j(2)+U2_j(3)*U2_j(3))/(U2_j(1)*U2_j(1))
      H_ij = (aux*hi+hj)/(aux+1.0_DP)

      ! Compute auxiliary variables
      aux   = u_ij*a(1) + v_ij*a(2)
      uPow2 = u_ij*u_ij
      vPow2 = v_ij*v_ij
      q_ij  = 0.5_DP*(uPow2+vPow2)
      cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
      cs  = sqrt(cPow2)

      ! Compute eigenvalues
      l1 = abs(aux-cs)
      l2 = abs(aux)
      l3 = abs(aux+cs)
      l4 = abs(aux)

      ! Compute solution difference U2_i-U2_j
      Diff = U2_i-U2_j

      ! Compute auxiliary quantities for characteristic variables
      aux1 = G2/cPow2*(q_ij*Diff(1)-u_ij*Diff(2)-v_ij*Diff(3)+Diff(4))
      aux2 = 0.5_DP*(aux*Diff(1)-a(1)*Diff(2)-a(2)*Diff(3))/cs

      ! Compute characteristic variables multiplied by the corresponding eigenvalue
      w1 = l1 * (aux1 + aux2)
      w2 = l2 * ((1.0_DP-G1*q_ij/cPow2)*Diff(1)+G1*(u_ij*Diff(2)+v_ij*Diff(3)-Diff(4))/cPow2)
      w3 = l3 * (aux1 - aux2)
      w4 = l4 * ((a(1)*v_ij-a(2)*u_ij)*Diff(1)+a(2)*Diff(2)-a(1)*Diff(3))

      ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
      Diff(1) = w1 + w2 + w3
      Diff(2) = (u_ij-cs*a(1))*w1 + u_ij*w2 + (u_ij+cs*a(1))*w3 + a(2)*w4
      Diff(3) = (v_ij-cs*a(2))*w1 + v_ij*w2 + (v_ij+cs*a(2))*w3 - a(1)*w4
      Diff(4) = (H_ij-cs*aux)*w1  + q_ij*w2 + (H_ij+cs*aux)*w3  + (u_ij*a(2)-v_ij*a(1))*w4

      ! Compute conservative flux
      F_ij = dscale1*(U1_i-U1_j) + dscale2*anorm*Diff

    else

      ! Compute conservative flux without spatial contribution
      F_ij = dscale1*(U1_i-U1_j)

    end if
  end subroutine euler_calcFluxFCTTensorDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxFCTRusanov2d(&
      U1_i, U1_j, U2_i, U2_j, C_ij, C_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 2D using the Rusanov dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

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
    real(DP) :: ui,vi,uj,vj
    real(DP) :: d_ij,ci,cj,Ei,Ej

    ! Compute velocities and energy
    ui = U2_i(2)/U2_i(1); vi = U2_i(3)/U2_i(1); Ei = U2_i(4)/U2_i(1)
    uj = U2_j(2)/U2_j(1); vj = U2_j(3)/U2_j(1); Ej = U2_j(4)/U2_j(1)

    ! Compute the speed of sound
    ci = sqrt(max(G15*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
    cj = sqrt(max(G15*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))

    ! Scalar dissipation for the Rusanov flux
    d_ij = max( abs(C_ij(1)*uj+C_ij(2)*vj) +&
                sqrt(C_ij(1)*C_ij(1)+C_ij(2)*C_ij(2))*cj,&
                abs(C_ji(1)*ui+C_ji(2)*vi) +&
                sqrt(C_ji(1)*C_ji(1)+C_ji(2)*C_ji(2))*ci )

    ! Compute conservative fluxes
    F_ij = dscale1*(U1_i-U1_j) + dscale2*d_ij*(U2_i-U2_j)

  end subroutine euler_calcFluxFCTRusanov2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxDensity2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density in 2D
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
      DtransformedFluxesAtEdge(1,1,idx) = DfluxesAtEdge(1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-DfluxesAtEdge(1,idx)
    end do

  end subroutine euler_trafoFluxDensity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffDensity2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density in 2D
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
      DtransformedDataAtEdge(1,idx) = DdataAtEdge(1,2,idx)-DdataAtEdge(1,1,idx)
    end do

  end subroutine euler_trafoDiffDensity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxEnergy2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the energy in 2D
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
      DtransformedFluxesAtEdge(1,1,idx) = DfluxesAtEdge(4,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-DfluxesAtEdge(4,idx)
    end do

  end subroutine euler_trafoFluxEnergy2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffEnergy2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the energy in 2D
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
      DtransformedDataAtEdge(1,idx) = DdataAtEdge(4,2,idx)-DdataAtEdge(4,1,idx)
    end do

  end subroutine euler_trafoDiffEnergy2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxPressure2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the pressure in 2D
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
    real(DP) :: ui,uj,vi,vj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      
      ! Transformed pressure fluxes
      DtransformedFluxesAtEdge(1,1,idx) = G1*(0.5_DP*(ui*ui+vi*vi)*DfluxesAtEdge(1,idx)-&
                                              ui*DfluxesAtEdge(2,idx)-&
                                              vi*DfluxesAtEdge(3,idx)+DfluxesAtEdge(4,idx))
      DtransformedFluxesAtEdge(1,2,idx) =-G1*(0.5_DP*(uj*uj+vj*vj)*DfluxesAtEdge(1,idx)-&
                                              uj*DfluxesAtEdge(2,idx)-&
                                              vj*DfluxesAtEdge(3,idx)+DfluxesAtEdge(4,idx))
    end do

  end subroutine euler_trafoFluxPressure2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffPressure2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the pressure in 2D
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
    real(DP) :: pi,pj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute pressures
      pi = G1*(DdataAtEdge(4,1,idx)-0.5_DP*(&
          DdataAtEdge(2,1,idx)*DdataAtEdge(2,1,idx)+&
          DdataAtEdge(3,1,idx)*DdataAtEdge(3,1,idx))/DdataAtEdge(1,1,idx))
      pj = G1*(DdataAtEdge(4,2,idx)-0.5_DP*(&
          DdataAtEdge(2,2,idx)*DdataAtEdge(2,2,idx)+&
          DdataAtEdge(3,2,idx)*DdataAtEdge(3,2,idx))/DdataAtEdge(1,2,idx))
     
      ! Transformed pressure difference
      DtransformedDataAtEdge(1,idx) = pj-pi
    end do

  end subroutine euler_trafoDiffPressure2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxVelocity2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the x- and y-velocity
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
    real(DP) :: ui,uj,vi,vj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)

      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(1,1,idx) = (DfluxesAtEdge(2,idx)-&
                                           ui*DfluxesAtEdge(1,idx))/DdataAtEdge(1,1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-(DfluxesAtEdge(2,idx)-&
                                           uj*DfluxesAtEdge(1,idx))/DdataAtEdge(1,2,idx)

      ! Transformed velocity fluxes in y-direction
      DtransformedFluxesAtEdge(2,1,idx) = (DfluxesAtEdge(3,idx)-&
                                           vi*DfluxesAtEdge(1,idx))/DdataAtEdge(1,1,idx)
      DtransformedFluxesAtEdge(2,2,idx) =-(DfluxesAtEdge(3,idx)-&
                                           vj*DfluxesAtEdge(1,idx))/DdataAtEdge(1,2,idx)
    end do
    
  end subroutine euler_trafoFluxVelocity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffVelocity2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the x- and y-velocity
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
      DtransformedDataAtEdge(1,idx) = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)-&
                                      DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)

      ! Transformed velocity difference in y-direction
      DtransformedDataAtEdge(2,idx) = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)-&
                                      DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
    end do
    
  end subroutine euler_trafoDiffVelocity2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxMomentum2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the x- and y-momentum
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
      DtransformedFluxesAtEdge(1,1,idx) = DfluxesAtEdge(2,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-DfluxesAtEdge(2,idx)

      ! Transformed momentum fluxes in y-direction
      DtransformedFluxesAtEdge(2,1,idx) = DfluxesAtEdge(3,idx)
      DtransformedFluxesAtEdge(2,2,idx) =-DfluxesAtEdge(3,idx)
    end do
    
  end subroutine euler_trafoFluxMomentum2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffMomentum2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the x- and y-momentum
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
      DtransformedDataAtEdge(1,idx) = DdataAtEdge(2,2,idx)-DdataAtEdge(2,1,idx)

      ! Transformed momentum difference in y-direction
      DtransformedDataAtEdge(2,idx) = DdataAtEdge(3,2,idx)-DdataAtEdge(3,1,idx)
    end do
    
  end subroutine euler_trafoDiffMomentum2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxDenEng2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 2D
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
      DtransformedFluxesAtEdge(1,1,idx) = DfluxesAtEdge(1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-DfluxesAtEdge(1,idx)

      ! Transformed total energy fluxes
      DtransformedFluxesAtEdge(2,1,idx) = DfluxesAtEdge(4,idx)
      DtransformedFluxesAtEdge(2,2,idx) =-DfluxesAtEdge(4,idx)
    end do

  end subroutine euler_trafoFluxDenEng2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffDenEng2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 2D
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
      DtransformedDataAtEdge(1,idx) = DdataAtEdge(1,2,idx)-DdataAtEdge(1,1,idx)

      ! Transformed total energy difference
      DtransformedDataAtEdge(2,idx) = DdataAtEdge(4,2,idx)-DdataAtEdge(4,1,idx)
    end do

  end subroutine euler_trafoDiffDenEng2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxDenPre2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 2D
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
    real(DP) :: ui,uj,vi,vj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)

      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) = DfluxesAtEdge(1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-DfluxesAtEdge(1,idx)
      
      ! Transformed pressure fluxes
      DtransformedFluxesAtEdge(2,1,idx) = G1*(0.5_DP*(ui*ui+vi*vi)*DfluxesAtEdge(1,idx)-&
                                          ui*DfluxesAtEdge(2,idx)-&
                                          vi*DfluxesAtEdge(3,idx)+DfluxesAtEdge(4,idx))
      DtransformedFluxesAtEdge(2,2,idx) =-G1*(0.5_DP*(uj*uj+vj*vj)*DfluxesAtEdge(1,idx)-&
                                          uj*DfluxesAtEdge(2,idx)-&
                                          vj*DfluxesAtEdge(3,idx)+DfluxesAtEdge(4,idx))
    end do

  end subroutine euler_trafoFluxDenPre2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffDenPre2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 2D
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
    real(DP) :: pi,pj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)

      ! Compute pressures
      pi = G1*(DdataAtEdge(4,1,idx)-0.5_DP*(&
          DdataAtEdge(2,1,idx)*DdataAtEdge(2,1,idx)+&
          DdataAtEdge(3,1,idx)*DdataAtEdge(3,1,idx))/DdataAtEdge(1,1,idx))
      pj = G1*(DdataAtEdge(4,2,idx)-0.5_DP*(&
          DdataAtEdge(2,2,idx)*DdataAtEdge(2,2,idx)+&
          DdataAtEdge(3,2,idx)*DdataAtEdge(3,2,idx))/DdataAtEdge(1,2,idx))
      
      ! Transformed density difference
      DtransformedDataAtEdge(1,idx) = DdataAtEdge(1,2,idx)-DdataAtEdge(1,1,idx)
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(2,idx) = pj-pi
    end do
    
  end subroutine euler_trafoDiffDenPre2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxDenPreVel2d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation
    ! of the given flux into primitive variables in 2D
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
    real(DP) :: ui,uj,vi,vj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) = DfluxesAtEdge(1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-DfluxesAtEdge(1,idx)

      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(2,1,idx) = (DfluxesAtEdge(2,idx)-&
                                           ui*DfluxesAtEdge(1,idx))/DdataAtEdge(1,1,idx)
      DtransformedFluxesAtEdge(2,2,idx) =-(DfluxesAtEdge(2,idx)-&
                                           uj*DfluxesAtEdge(1,idx))/DdataAtEdge(1,2,idx)

      ! Transformed velocity fluxes in y-direction
      DtransformedFluxesAtEdge(3,1,idx) = (DfluxesAtEdge(3,idx)-&
                                           vi*DfluxesAtEdge(1,idx))/DdataAtEdge(1,1,idx)
      DtransformedFluxesAtEdge(3,2,idx) =-(DfluxesAtEdge(3,idx)-&
                                           vj*DfluxesAtEdge(1,idx))/DdataAtEdge(1,2,idx)

      ! Transformed pressure fluxes
      DtransformedFluxesAtEdge(4,1,idx) = G1*(0.5_DP*(ui*ui+vi*vi)*DfluxesAtEdge(1,idx)-&
                                              ui*DfluxesAtEdge(2,idx)-&
                                              vi*DfluxesAtEdge(3,idx)+DfluxesAtEdge(4,idx))
      DtransformedFluxesAtEdge(4,2,idx) =-G1*(0.5_DP*(uj*uj+vj*vj)*DfluxesAtEdge(1,idx)-&
                                              uj*DfluxesAtEdge(2,idx)-&
                                              vj*DfluxesAtEdge(3,idx)+DfluxesAtEdge(4,idx))
    end do

  end subroutine euler_trafoFluxDenPreVel2d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffDenPreVel2d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density, pressure and velocity in 2D
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
    real(DP) :: pi,pj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)

      ! Compute pressures
      pi = G1*(DdataAtEdge(4,1,idx)-0.5_DP*(&
          DdataAtEdge(2,1,idx)*DdataAtEdge(2,1,idx)+&
          DdataAtEdge(3,1,idx)*DdataAtEdge(3,1,idx))/DdataAtEdge(1,1,idx))
      pj = G1*(DdataAtEdge(4,2,idx)-0.5_DP*(&
          DdataAtEdge(2,2,idx)*DdataAtEdge(2,2,idx)+&
          DdataAtEdge(3,2,idx)*DdataAtEdge(3,2,idx))/DdataAtEdge(1,2,idx))
      
      ! Transformed density difference
      DtransformedDataAtEdge(2,idx) = DdataAtEdge(1,2,idx)-DdataAtEdge(1,1,idx)

      ! Transformed velocity difference in x-direction
      DtransformedDataAtEdge(2,idx) = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)-&
                                      DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      
      ! Transformed velocity difference in y-direction
      DtransformedDataAtEdge(3,idx) = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)-&
                                      DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)

      ! Transformed pressure difference
      DtransformedDataAtEdge(4,idx) = pj-pi
    end do

  end subroutine euler_trafoDiffDenPreVel2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcBoundaryvalues2d(DbdrNormal, DpointNormal,&
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

    ! local variables
    real(DP), dimension(NVAR2D) :: W,Wu,Winf    ! Riemann invariants, eigenvalues, etc.
    real(DP) :: rho,v1,v2,p,E,c,v1_0,v2_0       ! primitive variables
    real(DP) :: v1_b,v2_b,vn_b,vn,vt,pstar,ps   ! velocities and boundary values
    real(DP) :: cup,f,fd,ge,qrt                 ! auxiliary variables ...
    real(DP) :: pold,ppv,prat,ptl,ptr,vdiff,vm  ! ... for the Riemann solver
    real(DP) :: auxA,auxB,aux,dnx2,dny2,dnxy
    integer:: ite

    ! What type of boundary condition is given?
    select case(ibdrCondType)
    case(BDRC_EULERWALL,&
         BDRC_RLXEULERWALL)
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
      p   = thdyn_pressure(GAMMA, E, rho, v1, v2)
      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))

      ! Precompute auxiliary data
      dnxy = DbdrNormal(1)*DbdrNormal(2)
      dnx2 = DbdrNormal(1)*DbdrNormal(1)
      dny2 = DbdrNormal(2)*DbdrNormal(2)

      if (ibdrCondType .eq. BDRC_EULERWALL) then
        ! Compute reflected velocities at the boundary
        v1_b = (dny2-dnx2)*v1 - 2.0_DP*dnxy*v2
        v2_b = (dnx2-dny2)*v2 - 2.0_DP*dnxy*v1
      else
        ! Compute initial velocity from previous time step
        v1_0 = Du0(2)/Du0(1)
        v2_0 = Du0(3)/Du0(1)

        ! Compute semi-reflected velocities at the boundary
        v1_b = (dny2-dnx2)*v1-2.0_DP*dnxy*v2 + 2*DbdrValue(1)*(v1_0*dnx2+v2_0*dnxy)
        v2_b = (dnx2-dny2)*v2-2.0_DP*dnxy*v1 + 2*DbdrValue(1)*(v2_0*dny2+v1_0*dnxy)
      end if

      ! Compute normal velocities at the boundary and the ghost state
      ! w.r.t. the numerical/approximate  outward unit normal vector
      vn   = DpointNormal(1)*v1   + DpointNormal(2)*v2
      vn_b = DpointNormal(1)*v1_b + DpointNormal(2)*v2_b

      ! Compute the tangential velocity depending on the sign of N*v
      if (vn .gt. 0.0_DP) then
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
      if (2.0_DP*G9*c .le. vn_b-vn) then
        if (present(istatus)) then
          istatus = -ibdrCondType
          return
        else
          call output_line('Riemann solver failed due to vacuum',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcBoundaryvalues2d')
          call sys_halt()
        end if
      end if

      ! Provide a guess value for pressure in the "star region"
      ! by using the PVRS Riemann solver as suggested by Toro

      cup  = rho*c
      ppv  = p+0.5_DP*(vn-vn_b)*cup
      ppv  = max(0.0_DP, ppv)

      if (ppv .eq. p) then

        ! Select guessed pressure from PVRS Riemann solver
        pstar = ppv
      else
        if (ppv .lt. p) then

          ! Guess pressure from the Two-Rarefaction Riemann solver
          vm    = 0.5_DP*(vn+vn_b)
          ptl   = 1.0_DP + G2*(vn-vm)/c
          ptr   = 1.0_DP + G2*(vm-vn_b)/c
          pstar = 0.5_DP*(p*ptl + p*ptr)**G11
        else

          ! Guess pressure from the Two-Shock Riemann solver
          ! with PVRS as estimated pressure value
          ge    = sqrt((G10/rho)/(G12*p+ppv))
          pstar = p - 0.5_DP*(vn_b-vn)/ge
        end if
      end if

      ! Initialize solution difference and pressure
      vdiff = (vn_b-vn)/2.0_DP
      pold  = pstar

      newton: do ite = 1, 100

        ! Compute pressure function f(pold) and its derivative f1(pold)
        if (pold .le. p) then

          ! Rarefaction wave
          prat = pold/p

          f  = G9*c*(prat**G7 - 1.0_DP)
          fd = (1.0_DP/(rho*c))*prat**(-G8)
        else

          ! Shock wave
          auxA = G10/rho
          auxB = G12*p
          qrt  = sqrt(auxA/(auxB + pold))

          f  = (pold-p)*qrt
          fd = (1.0_DP - 0.5_DP*(pold - p)/(auxB + pold))*qrt
        end if

        pstar = pold - (f+vdiff)/fd
        if (pstar .lt. 0.0_DP) then
          pold = 1.0E-6
          cycle newton
        end if

        aux = 2.0_DP*abs((pstar-pold)/(pstar+pold))
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
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcBoundaryvalues2d')
          call sys_halt()
        end if
      end if

      !-------------------------------------------------------------------------
      ! Calculate the velocity in the star region
      !-------------------------------------------------------------------------

      ! Note that the contribution fR-fL vanishes due to constant states
      vn = 0.5_DP*(vn+vn_b)


      !-------------------------------------------------------------------------
      ! Calculate the density in the star region
      !-------------------------------------------------------------------------

      if (pstar .le. p) then

        ! Rarefaction wave
        rho = rho*(pstar/p)**G4
      else

        ! Shock wave
        rho = rho*(pstar/p+G12)/(G12*(pstar/p)+1.0_DP)
      end if

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DpointNormal(2)*vt+DpointNormal(1)*vn)
      Du(3) = rho*(-DpointNormal(1)*vt+DpointNormal(2)*vn)
      Du(4) = pstar/G1+0.5_DP*rho*(vn*vn+vt*vt)


    case(BDRC_VISCOUSWALL)
      !-------------------------------------------------------------------------

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      v2  = Du(3)/rho
      E   = Du(4)/rho
      p   = thdyn_pressure(GAMMA, E, rho, v1, v2)

      ! Update the solution vector and let vn:=0 and vt:=0
      Du(2) = 0.0_DP
      Du(3) = 0.0_DP
      Du(4) = p/G1


    case(BDRC_SUPERINLET)
      !-------------------------------------------------------------------------

      ! The free stream primitive variables are Deval=[rho,v1,v2,p]
      rho = DbdrValue(1)
      p   = DbdrValue(4)
      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*DbdrValue(2)+DbdrNormal(2)*DbdrValue(3)
      vt  = DbdrNormal(2)*DbdrValue(2)-DbdrNormal(1)*DbdrValue(3)

      ! Compute Riemann invariants based on the free stream values
      W(1) = vn-2*c/G1
      W(2) = vn+2*c/G1
      W(3) = p/(rho**GAMMA)
      W(4) = vt

      ! Transform back into conservative variables
      vn   = 0.5_DP*(W(1)+W(2))
      c    = 0.25_DP*G1*(W(2)-W(1))
      rho  = (c*c/GAMMA/W(3))**G3
      p    = rho*c*c/GAMMA
      vt   = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/G1+0.5_DP*rho*(vn*vn+vt*vt)


    case(BDRC_FREESTREAM)
      !-------------------------------------------------------------------------

      ! The free stream primitive variables are Deval=[rho,v1,v2,p]
      rho = DbdrValue(1)
      p   = DbdrValue(4)
      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*DbdrValue(2)+DbdrNormal(2)*DbdrValue(3)
      vt  = DbdrNormal(2)*DbdrValue(2)-DbdrNormal(1)*DbdrValue(3)

      ! Compute Riemann invariants based on the free stream values
      Winf(1) = vn-2.0_DP*c/G1
      Winf(2) = vn+2.0_DP*c/G1
      Winf(3) = p/(rho**GAMMA)
      Winf(4) = vt

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      v2  = Du(3)/rho
      E   = Du(4)/rho
      p   = thdyn_pressure(GAMMA, E, rho, v1, v2)

      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*v1+DbdrNormal(2)*v2
      vt  = DbdrNormal(2)*v1-DbdrNormal(1)*v2

      ! Compute Riemann invariants based on the solution values
      Wu(1) = vn-2.0_DP*c/G1
      Wu(2) = vn+2.0_DP*c/G1
      Wu(3) = p/(rho**GAMMA)
      Wu(4) = vt

      ! Adopt free stream/computed values depending on the sign of the eigenvalue
      W(1) = merge(Winf(1), Wu(1), vn <  c)
      W(2) = merge(Winf(2), Wu(2), vn < -c)
      W(3) = merge(Winf(3), Wu(3), vn <  SYS_EPSREAL)
      W(4) = merge(Winf(4), Wu(4), vn <  SYS_EPSREAL)

      ! Transform back into conservative variables
      vn   = 0.5_DP*(W(1)+W(2))
      c    = 0.25_DP*G1*(W(2)-W(1))
      rho  = (c*c/GAMMA/W(3))**G3
      p    = rho*c*c/GAMMA
      vt   = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/G1+0.5_DP*rho*(vn*vn+vt*vt)


    case(BDRC_SUBINLET)
      !-------------------------------------------------------------------------

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      v2  = Du(3)/rho
      E   = Du(4)/rho
      p   = thdyn_pressure(GAMMA, E, rho, v1, v2)

      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*v1+DbdrNormal(2)*v2
      vt  = DbdrNormal(2)*v1-DbdrNormal(1)*v2

      ! The specified density and pressure is Deval=[rho,p]
      rho = DbdrValue(1)
      p   = DbdrValue(2)

      ! Compute Riemann invariants
      W(1) = vn-2.0_DP*c/G1
      W(2) = vn+2.0_DP*c/G1
      W(3) = p/(rho**GAMMA)
      W(4) = vt

      ! Transform back into conservative variables
      vn   = 0.5_DP*(W(1)+W(2))
      c    = 0.25_DP*G1*(W(2)-W(1))
      rho  = (c*c/GAMMA/W(3))**G3
      p    = rho*c*c/GAMMA
      vt   = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/G1+0.5_DP*rho*(vn*vn+vt*vt)


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
      p   = thdyn_pressure(GAMMA, E, rho, v1, v2)

      vn  = DbdrNormal(1)*v1+DbdrNormal(2)*v2
      vt  = DbdrNormal(2)*v1-DbdrNormal(1)*v2
      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))

      ! Compute Riemann invariants based on the solution values and prescribed exit pressure
      W(2) = 2*c/G1-vn
      W(3) = p/(rho**GAMMA)
      W(4) = vt
      W(1) = 4/G1*sqrt(max(GAMMA*ps/rho*(p/ps)**G4, SYS_EPSREAL))-W(2)

      ! Transform back into conservative variables
      vn  = 0.5_DP*(W(1)-W(2))
      c   = 0.25_DP*G1*(W(1)+W(2))
      rho = (c*c/GAMMA/W(3))**G3
      p   = rho*c*c/GAMMA
      vt  = W(4)

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*( DbdrNormal(2)*vt+DbdrNormal(1)*vn)
      Du(3) = rho*(-DbdrNormal(1)*vt+DbdrNormal(2)*vn)
      Du(4) = p/G1+0.5_DP*rho*(vn*vn+vt*vt)


    case DEFAULT
      call output_line('Unsupported type of boundary condition!',&
          OU_CLASS_ERROR,OU_MODE_STD,'euler_calcBoundaryvalues2d')
      call sys_halt()
    end select

  end subroutine euler_calcBoundaryvalues2d

  ! ***************************************************************************

!<subroutine>

  subroutine euler_coeffVectorBdr2d_sim(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

    use basicgeometry
    use boundary
    use boundaryfilter
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

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional
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

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution
    real(DP), dimension(:,:), pointer :: Daux1
    real(DP), dimension(:,:,:), pointer :: Daux2
    real(DP), dimension(NVAR2D) :: DstateI,DstateM,Dflux,Diff
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dnx,dny,dtime,dscale,pI,cI,rM,pM,cM,dvnI,dvtI,dvnM,dvtM
    integer :: ibdrtype,isegment,iel,ipoint,ndim,ivar,nvar,iexpr,nmaxExpr

    ! This subroutine assumes that the first quick access string
    ! value holds the name of the function parser in the collection.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(1)))

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
      allocate(Daux1(ubound(Dpoints,2)*nvar, ubound(Dpoints,3)))
      
      ! Evaluate the solution in the cubature points on the boundary
      call fevl_evaluate_sim(DER_FUNC, Daux1, p_rsolution%RvectorBlock(1),&
          Dpoints, rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! What type of boundary conditions are we?
      select case(ibdrtype)
        
      case (BDRC_FREESTREAM_WEAK)
        !-----------------------------------------------------------------------
        ! Free-stream boundary conditions:

        ! Initialize values for function parser
        Dvalue = 0.0_DP
        Dvalue(NDIM3D+1) = dtime

        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)

        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)

            ! Get the normal vector in the point from the boundary.
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)
            
            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute free stream values from function parser given in
            ! term of the primitive variables [rho,v1,v2,p]
            do iexpr = 1, 4
              call fparser_evalFunction(p_rfparser,&
                  nmaxExpr*(isegment-1)+iexpr, Dvalue, DstateM(iexpr))
            end do

            ! Compute auxiliary quantities based on free stream state vector
            rM = DstateM(1)
            pM = DstateM(4)
            cM = sqrt(max(GAMMA*pM/rM, SYS_EPSREAL))
            dvnM =  dnx*DstateM(2)+dny*DstateM(3)
            dvtM = -dny*DstateM(2)+dnx*DstateM(3)
            
            ! Compute auxiliary quantities based on internal state vector
            pI = G1*(Daux1((ipoint-1)*NVAR2D+4,iel)-0.5_DP*&
                     (Daux1((ipoint-1)*NVAR2D+2,iel)**2+&
                      Daux1((ipoint-1)*NVAR2D+3,iel)**2))/&
                 Daux1((ipoint-1)*NVAR2D+1,iel)
            cI = sqrt(max(GAMMA*pI/Daux1((ipoint-1)*NVAR2D+1,iel), SYS_EPSREAL))

            ! Compute the normal and tangential velocities based
            ! on internal state vector
            dvnI = ( dnx*Daux1((ipoint-1)*NVAR2D+2,iel)+&
                     dny*Daux1((ipoint-1)*NVAR2D+3,iel) )/&
                     Daux1((ipoint-1)*NVAR2D+1,iel)
            dvtI = (-dny*Daux1((ipoint-1)*NVAR2D+2,iel)+&
                     dnx*Daux1((ipoint-1)*NVAR2D+3,iel) )/&
                     Daux1((ipoint-1)*NVAR2D+1,iel)
            
            ! Select free stream or computed Riemann invariant depending
            ! on the sign of the corresponding eigenvalue
            if (dvnI .lt. cI) then
              DstateM(1) = dvnM-G9*cM
            else
              DstateM(1) = dvnI-G9*cI
            end if

            if (dvnI .lt. -cI) then
              DstateM(2) = dvnM+G9*cM
            else
              DstateM(2) = dvnI+G9*cI
            end if

            if (dvnI .lt. SYS_EPSREAL) then
              DstateM(3) = pM/(rM**GAMMA)
              DstateM(4) = dvtM
            else
              DstateM(3) = pI/(Daux1((ipoint-1)*NVAR2D+1,iel)**GAMMA)
              DstateM(4) = dvtI
            end if

            ! Convert Riemann invariants into conservative state variables
            cM = 0.25_DP*G1*(DstateM(2)-DstateM(1))
            rM = (cM*cM/(GAMMA*DstateM(3)))**G3
            pM = rM*cM*cM/GAMMA
            dvnM = 0.5_DP*(DstateM(1)+DstateM(2))
            dvtM = DstateM(4)

            ! Setup the state vector based on Riemann invariants
            DstateM(1) = rM
            DstateM(2) = rM*( dnx*dvnM+dny*dvtM)
            DstateM(3) = rM*(-dny*dvnM+dnx*dvtM)
            DstateM(4) = pM/G1+0.5_DP*(dvnM*dvnM+dvtM*dvtM)
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux1((ipoint-1)*NVAR2D+1,iel)
            DstateI(2) = Daux1((ipoint-1)*NVAR2D+2,iel)
            DstateI(3) = Daux1((ipoint-1)*NVAR2D+3,iel)
            DstateI(4) = Daux1((ipoint-1)*NVAR2D+4,iel)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*0.5_DP*(Dflux-Diff)
          end do
        end do

        
      case (BDRC_EULERWALL_WEAK)
        !-----------------------------------------------------------------------
        ! Euler wall boundary condition:
        
        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)
            
            ! Get the normal vector in the point from the boundary.
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux1((ipoint-1)*NVAR2D+1,iel)
            DstateI(2) = Daux1((ipoint-1)*NVAR2D+2,iel)
            DstateI(3) = Daux1((ipoint-1)*NVAR2D+3,iel)
            DstateI(4) = Daux1((ipoint-1)*NVAR2D+4,iel)
            
            ! Compute the normal and tangential velocities based
            ! on the internal state vector
            dvnI = ( dnx*Daux1((ipoint-1)*NVAR2D+2,iel)+&
                     dny*Daux1((ipoint-1)*NVAR2D+3,iel) )/&
                     Daux1((ipoint-1)*NVAR2D+1,iel)
            dvtI = (-dny*Daux1((ipoint-1)*NVAR2D+2,iel)+&
                     dnx*Daux1((ipoint-1)*NVAR2D+3,iel) )/&
                     Daux1((ipoint-1)*NVAR2D+1,iel)

            ! Compute the mirrored state vector
            DstateM(1) = DstateI(1)
            DstateM(2) = DstateM(1)*(-dvnI*dnx - dvtI*dny)
            DstateM(3) = DstateM(1)*(-dvnI*dny + dvtI*dnx)
            DstateM(4) = DstateI(4)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)

            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*0.5_DP*(Dflux-Diff)
          end do
        end do

      case default

        print *, ibdrtype
        call output_line('Invalid type of boundary conditions!',&
            OU_CLASS_ERROR,OU_MODE_STD,'euler_coeffVectorBdr2d_sim')
        call sys_halt()
        
      end select

      ! Deallocate temporal memory
      deallocate(Daux1)
        
    else

      !-------------------------------------------------------------------------
      ! Solution is stored in block format
      !-------------------------------------------------------------------------
      
      ! Allocate temporal memory
      allocate(Daux2(ubound(Dpoints,2), ubound(Dpoints,3), nvar))
      
      ! Evaluate the solution in the cubature points on the boundary
      do ivar = 1, nvar
        call fevl_evaluate_sim(DER_FUNC, Daux2(:,:,ivar),&
            p_rsolution%RvectorBlock(ivar), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      end do

      ! What type of boundary conditions are we?
      select case(ibdrtype)
        
      case (BDRC_FREESTREAM_WEAK)
        !-----------------------------------------------------------------------
        ! Free-stream boundary conditions:

        ! Initialize values for function parser
        Dvalue = 0.0_DP
        Dvalue(NDIM3D+1) = dtime

        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)

        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)

            ! Get the normal vector in the point from the boundary.
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)
            
            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

            ! Compute free stream values from function parser given in
            ! term of the primitive variables [rho,v1,v2,p]
            do iexpr = 1, 4
              call fparser_evalFunction(p_rfparser,&
                  nmaxExpr*(isegment-1)+iexpr, Dvalue, DstateM(iexpr))
            end do

            ! Compute auxiliary quantities based on free stream state vector
            rM = DstateM(1)
            pM = DstateM(4)
            cM = sqrt(max(GAMMA*pM/rM, SYS_EPSREAL))
            dvnM =  dnx*DstateM(2)+dny*DstateM(3)
            dvtM = -dny*DstateM(2)+dnx*DstateM(3)

            ! Compute auxiliary quantities based on internal state vector
            pI = G1*(Daux2(ipoint,iel,4)-0.5_DP*&
                     (Daux2(ipoint,iel,2)**2+&
                      Daux2(ipoint,iel,3)**2))/Daux2(ipoint,iel,1)
            cI = sqrt(max(GAMMA*pI/Daux2(ipoint,iel,1), SYS_EPSREAL))

            ! Compute the normal and tangential velocities based
            ! on internal state vector
            dvnI = ( dnx*Daux2(ipoint,iel,2)+&
                     dny*Daux2(ipoint,iel,3))/Daux2(ipoint,iel,1)
            dvtI = (-dny*Daux2(ipoint,iel,2)+&
                     dnx*Daux2(ipoint,iel,3))/Daux2(ipoint,iel,1)

            ! Select free stream or computed Riemann invariant depending
            ! on the sign of the corresponding eigenvalue
            if (dvnI .lt. cI) then
              DstateM(1) = dvnM-G9*cM
            else
              DstateM(1) = dvnI-G9*cI
            end if

            if (dvnI .lt. -cI) then
              DstateM(2) = dvnM+G9*cM
            else
              DstateM(2) = dvnI+G9*cI
            end if

            if (dvnI .lt. SYS_EPSREAL) then
              DstateM(3) = pM/(rM**GAMMA)
              DstateM(4) = dvtM
            else
              DstateM(3) = pI/(Daux2(ipoint,iel,1)**GAMMA)
              DstateM(4) = dvtI
            end if

            ! Convert Riemann invariants into conservative state variables
            cM = 0.25_DP*G1*(DstateM(2)-DstateM(1))
            rM = (cM*cM/(GAMMA*DstateM(3)))**G3
            pM = rM*cM*cM/GAMMA
            dvnM = 0.5_DP*(DstateM(1)+DstateM(2))
            dvtM = DstateM(4)

            ! Setup the state vector based on Riemann invariants
            DstateM(1) = rM
            DstateM(2) = rM*( dnx*dvnM+dny*dvtM)
            DstateM(3) = rM*(-dny*dvnM+dnx*dvtM)
            DstateM(4) = pM/G1+0.5_DP*(dvnM*dvnM+dvtM*dvtM)
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux2(ipoint,iel,1)
            DstateI(2) = Daux2(ipoint,iel,2)
            DstateI(3) = Daux2(ipoint,iel,3)
            DstateI(4) = Daux2(ipoint,iel,4)
            
            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*0.5_DP*(Dflux-Diff)
          end do
        end do


      case (BDRC_EULERWALL_WEAK)
       
 !-----------------------------------------------------------------------
        ! Euler wall boundary condition:
        
        do iel = 1, size(rdomainIntSubset%p_Ielements)
          do ipoint = 1, ubound(Dpoints,2)
            
            ! Get the normal vector in the point from the boundary
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, DpointPar(ipoint,iel), dnx, dny, cparType=BDR_PAR_LENGTH)
            
            ! Setup the computed internal state vector
            DstateI(1) = Daux2(ipoint,iel,1)
            DstateI(2) = Daux2(ipoint,iel,2)
            DstateI(3) = Daux2(ipoint,iel,3)
            DstateI(4) = Daux2(ipoint,iel,4)

            ! Compute the normal and tangential velocities based
            ! on the internal state vector
            dvnI = ( dnx*Daux2(ipoint,iel,2)+&
                     dny*Daux2(ipoint,iel,3) )/Daux2(ipoint,iel,1)
            dvtI = (-dny*Daux2(ipoint,iel,2)+&
                      dnx*Daux2(ipoint,iel,3) )/Daux2(ipoint,iel,1)

            ! Compute the mirrored state vector
            DstateM(1) = DstateI(1)
            DstateM(2) = DstateM(1)*(-dvnI*dnx - dvtI*dny)
            DstateM(3) = DstateM(1)*(-dvnI*dny + dvtI*dnx)
            DstateM(4) = DstateI(4)

            ! Invoke Riemann solver
            call doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)
            
            ! Store flux in the cubature points
            Dcoefficients(:,1,ipoint,iel) = dscale*0.5_DP*(Dflux-Diff)
          end do
        end do

      case default
        call output_line('Invalid type of boundary conditions!',&
            OU_CLASS_ERROR,OU_MODE_STD,'euler_coeffVectorBdr2d_sim')
        call sys_halt()
        
      end select

      ! Deallocate temporal memory
      deallocate(Daux2)
      
    end if

  contains

    ! Here comes the working routine

    !***************************************************************************
    ! Approximate Riemann solver along the outward unit normal
    !***************************************************************************

    subroutine doRiemannSolver(DstateI, DstateM, dnx, dny, Dflux, Diff)

      ! input parameters
      real(DP), dimension(NVAR2D), intent(in) :: DstateI, DstateM
      real(DP), intent(in) :: dnx, dny

      ! output parameters
      real(DP), dimension(NVAR2D), intent(out) :: Dflux, Diff

      ! local variables
      real(DP) :: uI,vI,ru2I,rv2I,uM,vM,ru2M,rv2M,hI,hM
      real(DP) :: cPow2,uPow2,vPow2,c_IM,h_IM,q_IM,u_IM,v_IM
      real(DP) :: l1,l2,l3,l4,w1,w2,w3,w4,aux,aux1,aux2,dveln
      
      ! Compute auxiliary quantities
      uI = DstateI(2)/DstateI(1)
      vI = DstateI(3)/DstateI(1)
      ru2I = uI*DstateI(2)
      rv2I = vI*DstateI(3)
      
      ! Compute auxiliary quantities
      uM = DstateM(2)/DstateM(1)
      vM = DstateM(3)/DstateM(1)
      ru2M = uM*DstateM(2)
      rv2M = vM*DstateM(3)

      ! Calculate $\frac12{\bf n}\cdot[{\bf F}(U_I)+{\bf F}(U_M)]$
      Dflux(1) = dnx*(DstateI(2)+DstateM(2))+&
                 dny*(DstateI(3)+DstateM(3))
      Dflux(2) = dnx*(G1*(DstateI(4)+DstateM(4))-&
                      G14*(ru2I+ru2M)-&
                       G2*(rv2I+rv2M))+&
                 dny*(DstateI(2)*vI+DstateM(2)*vM)
      Dflux(3) = dnx*(DstateI(3)*uI+DstateM(3)*uM)+&
                 dny*(G1*(DstateI(4)+DstateM(4))-&
                      G14*(rv2I+rv2M)-&
                       G2*(ru2I+ru2M))
      Dflux(4) = dnx*((GAMMA*DstateI(4)-G2*(ru2I+rv2I))*uI+&
                      (GAMMA*DstateM(4)-G2*(ru2M+rv2M))*uM)+&
                 dny*((GAMMA*DstateI(4)-G2*(ru2I+rv2I))*vI+&
                      (GAMMA*DstateM(4)-G2*(ru2M+rv2M))*vM)
      
      ! Compute Roe mean values
      aux  = sqrt(max(DstateI(1)/DstateM(1), SYS_EPSREAL))
      u_IM = (aux*uI+uM)/(aux+1.0_DP)
      v_IM = (aux*vI+vM)/(aux+1.0_DP)
      hI   = GAMMA*DstateI(4)/DstateI(1)-G2*(uI**2+vI**2)
      hM   = GAMMA*DstateM(4)/DstateM(1)-G2*(uM**2+vM**2)
      h_IM =(aux*hI+hM)/(aux+1.0_DP)
      
      ! Compute auxiliary variables
      uPow2 = u_IM*u_IM
      vPow2 = v_IM*v_IM
      q_IM  = 0.5_DP*(uPow2+vPow2)
      cPow2 = max(G1*(H_IM-q_IM), SYS_EPSREAL)
      c_IM  = sqrt(cPow2)

      ! Compute normal velocity
      dveln =  dnx*uI+dny*vI

      ! Compute eigenvalues
      l1 = abs(dveln-c_IM)
      l2 = abs(dveln)
      l3 = abs(dveln+c_IM)
      l4 = abs(dveln)
      
      ! Compute solution difference U_M-U_I
      Diff = DstateM-DstateI
      
      ! Compute auxiliary quantities for characteristic variables
      aux1 = G2/cPow2*(q_IM*Diff(1)-u_IM*Diff(2)-v_IM*Diff(3)+Diff(4))
      aux2 = 0.5_DP*(aux*Diff(1)-dnx*Diff(2)-dny*Diff(3))/c_IM
      
      ! Compute characteristic variables multiplied by the corresponding eigenvalue
      w1 = l1 * (aux1 + aux2)
      w2 = l2 * ((1.0_DP-G1*q_IM/cPow2)*Diff(1)+G1*(u_IM*Diff(2)+v_IM*Diff(3)-Diff(4))/cPow2)
      w3 = l3 * (aux1 - aux2)
      w4 = l4 * ((dnx*v_IM-dny*u_IM)*Diff(1)+dny*Diff(2)-dnx*Diff(3))
      
      ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
      Diff(1) = w1 + w2 + w3
      Diff(2) = (u_IM-c_IM*dnx)*w1 + u_IM*w2 + (u_IM+c_IM*dnx)*w3 + dny*w4
      Diff(3) = (v_IM-c_IM*dny)*w1 + v_IM*w2 + (v_IM+c_IM*dny)*w3 - dnx*w4
      Diff(4) = (H_IM-c_IM*dveln)*w1  + q_IM*w2 +&
                (H_IM+c_IM*dveln)*w3  + (u_IM*dny-v_IM*dnx)*w4

    end subroutine doRiemannSolver
    
  end subroutine euler_coeffVectorBdr2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine euler_hadaptCallbackScalar2d(iOperation, rcollection)

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
            OU_CLASS_WARNING,OU_MODE_STD,'euler_hadaptCallbackScalar2d')
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
            0.5_DP*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR2D+ivar)+&
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
            0.25_DP*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR2D+ivar)+&
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
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR2D+ivar) = 0.0_DP
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)

    end select

  end subroutine euler_hadaptCallbackScalar2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_hadaptCallbackBlock2d(iOperation, rcollection)

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
      if (rsolution%nblocks .ne. NVAR2D) then
        call output_line('Vector is not in block format!',&
            OU_CLASS_WARNING,OU_MODE_STD,'euler_hadaptCallbackBlock2d')
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
            0.5_DP*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
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
            0.25_DP*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
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
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = 0.0_DP
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)

    end select

  end subroutine euler_hadaptCallbackBlock2d

end module euler_callback2d
