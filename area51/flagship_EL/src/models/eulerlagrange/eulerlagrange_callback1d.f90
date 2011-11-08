!##############################################################################
!# ****************************************************************************
!# <name> eulerlagrange_callback1d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible Euler/Navier-Stokes equations in 1D.
!#
!# The following callback functions are available:
!#
!# 1.) eulerlagrange_calcFluxGal1d_sim
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#
!# 2.) eulerlagrange_calcFluxGalNoBdr1d_sim
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#        without assembling the symmetric boundary contribution
!#
!# 3.) eulerlagrange_calcFluxScDiss1d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting scalar artificial viscosities
!#
!# 4.) eulerlagrange_calcFluxRoeDiss1d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities
!#
!# 5.) eulerlagrange_calcFluxRusDiss1d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion
!#
!# 6.) eulerlagrange_calcMatDiagMatD1d_sim
!#     -> Computes local matrix for diagonal entry
!#
!# 7.) eulerlagrange_calcMatDiag1d_sim
!#     -> Computes local matrix for diagonal entry
!#
!# 8.) eulerlagrange_calcMatGalMatD1d_sim
!#     -> Computes local matrices for standard Galerkin scheme
!#
!# 9.) eulerlagrange_calcMatGal1d_sim
!#     -> Computes local matrices for standard Galerkin scheme
!#
!# 10.) eulerlagrange_calcMatScDissMatD1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 11.) eulerlagrange_calcMatScDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 12.) eulerlagrange_calcMatRoeDissMatD1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 13.) eulerlagrange_calcMatRoeDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 14.) eulerlagrange_calcMatRusDissMatD1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov artificial viscosities
!#
!# 15.) eulerlagrange_calcMatRusDiss1d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov flux artificial viscosities
!#
!# 16.) eulerlagrange_calcCharacteristics1d_sim
!#      -> Computes characteristic variables
!#
!# 17.) eulerlagrange_calcFluxFCTScalarDiss1d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting scalar artificial viscosities
!#
!# 18.) eulerlagrange_calcFluxFCTTensorDiss1d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting tensorial artificial viscosities
!#
!# 19.) eulerlagrange_calcFluxFCTRusanov1d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting the Rusanov artificial viscosities
!#
!# 20.) eulerlagrange_trafoFluxDensity1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 21.) eulerlagrange_trafoDiffDensity1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 22.) eulerlagrange_trafoFluxEnergy1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 23.) eulerlagrange_trafoDiffEnergy1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 24.) eulerlagrange_trafoFluxPressure1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 25.) eulerlagrange_trafoDiffPressure1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the pressure
!#
!# 26.) eulerlagrange_trafoFluxVelocity1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 27.) eulerlagrange_trafoDiffVelocity1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 28.) eulerlagrange_trafoFluxMomentum1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 29.) eulerlagrange_trafoDiffMomentum1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 30.) eulerlagrange_trafoFluxDenEng1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 31.) eulerlagrange_trafoDiffDenEng1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 32.) eulerlagrange_trafoFluxDenPre1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 33.) eulerlagrange_trafoDiffDenPre1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 34.) eulerlagrange_trafoFluxDenPreVel1d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 35.) eulerlagrange_trafoDiffDenPreVel1d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure
!#         and the velocity
!#
!# 36.) eulerlagrange_calcBoundaryvalues1d
!#      -> Computes the boundary values for a given node
!#
!# 37.) eulerlagrange_hadaptCallbackScalar1d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 1D, whereby the vector is stored in interleave format
!#
!# 38.) eulerlagrange_hadaptCallbackBlock1d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 1D, whereby the vector is stored in block format
!#
!# </purpose>
!##############################################################################

module eulerlagrange_callback1d

  use boundaryfilter
  use collection
  use eulerlagrange_basic
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
  public :: eulerlagrange_calcFluxGal1d_sim
  public :: eulerlagrange_calcFluxGalNoBdr1d_sim
  public :: eulerlagrange_calcFluxScDiss1d_sim
  public :: eulerlagrange_calcFluxRoeDiss1d_sim
  public :: eulerlagrange_calcFluxRusDiss1d_sim
  public :: eulerlagrange_calcMatDiagMatD1d_sim
  public :: eulerlagrange_calcMatDiag1d_sim
  public :: eulerlagrange_calcMatGalMatD1d_sim
  public :: eulerlagrange_calcMatGal1d_sim
  public :: eulerlagrange_calcMatScDissMatD1d_sim
  public :: eulerlagrange_calcMatScDiss1d_sim
  public :: eulerlagrange_calcMatRoeDissMatD1d_sim
  public :: eulerlagrange_calcMatRoeDiss1d_sim
  public :: eulerlagrange_calcMatRusDissMatD1d_sim
  public :: eulerlagrange_calcMatRusDiss1d_sim
  public :: eulerlagrange_calcCharacteristics1d_sim
  public :: eulerlagrange_calcFluxFCTScalarDiss1d
  public :: eulerlagrange_calcFluxFCTTensorDiss1d
  public :: eulerlagrange_calcFluxFCTRusanov1d
  public :: eulerlagrange_trafoFluxDensity1d_sim
  public :: eulerlagrange_trafoFluxEnergy1d_sim
  public :: eulerlagrange_trafoFluxPressure1d_sim
  public :: eulerlagrange_trafoFluxVelocity1d_sim
  public :: eulerlagrange_trafoFluxMomentum1d_sim
  public :: eulerlagrange_trafoFluxDenEng1d_sim
  public :: eulerlagrange_trafoFluxDenPre1d_sim
  public :: eulerlagrange_trafoFluxDenPreVel1d_sim
  public :: eulerlagrange_trafoDiffDensity1d_sim
  public :: eulerlagrange_trafoDiffEnergy1d_sim
  public :: eulerlagrange_trafoDiffPressure1d_sim
  public :: eulerlagrange_trafoDiffVelocity1d_sim
  public :: eulerlagrange_trafoDiffMomentum1d_sim
  public :: eulerlagrange_trafoDiffDenEng1d_sim
  public :: eulerlagrange_trafoDiffDenPre1d_sim
  public :: eulerlagrange_trafoDiffDenPreVel1d_sim
  public :: eulerlagrange_calcBoundaryvalues1d
  public :: eulerlagrange_hadaptCallbackScalar1d
  public :: eulerlagrange_hadaptCallbackBlock1d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxGal1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IdofsAtEdge, dscale, DfluxesAtEdge, rcollection)

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
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

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
#ifdef USE_EULER_IBP
    real(DP), dimension(NVAR1D) :: dF_i, dF_j
#else
    real(DP), dimension(NVAR1D) :: dF_ij
#endif
    real(DP) :: ui,uj,ru2i,ru2j
    integer :: idx


    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin flux
      !
      !     / rho*u       \
      ! F = | rho*u*u + p |
      !     \ rho*u*H     /
      !
      ! Here, we do not compute the pressure p and the enthalpy H but we
      ! calculate the flux from the conservative variables as follows:
      !
      !     / U2                      \
      ! F = | G1*U3-G14*ru2i          |
      !     \ (gamma*U3-G2*ru2i)*ui /
      !
      ! where the auxiliary values for node i are defined as follows:
      !
      ! ru2i = U2*U2/U1 = ui*U2
      ! ui = U2/U1
      !
      ! and the predefined constants are given by:
      !
      ! G14 = (gamma-3)/2   and   G2 = (gamma-1)/2   and   G1 = gamma-1
      !
      ! The auxiliary values for node j are defined accordingly.
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      
      ! Compute auxiliary variables
      ru2i = ui*DdataAtEdge(2,1,idx)
      ru2j = uj*DdataAtEdge(2,2,idx)

#ifdef USE_EULER_IBP
      ! Compute fluxes for x-direction
      dF_i(1) = DdataAtEdge(2,1,idx)
      dF_i(2) = G1*DdataAtEdge(3,1,idx)-G14*ru2i
      dF_i(3) = (GAMMA*DdataAtEdge(3,1,idx)-G2*ru2i)*ui
      
      dF_j(1) = DdataAtEdge(2,2,idx)
      dF_j(2) = G1*DdataAtEdge(3,2,idx)-G14*ru2j
      dF_j(3) = (GAMMA*DdataAtEdge(3,2,idx)-G2*ru2j)*uj
      
      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF_i )
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
#else
      ! Compute flux difference for x-direction
      dF_ij(1) = DdataAtEdge(2,1,idx) - DdataAtEdge(2,2,idx)
      dF_ij(2) = (G1*DdataAtEdge(3,1,idx)-G14*ru2i)-&
                 (G1*DdataAtEdge(3,2,idx)-G14*ru2j)
      dF_ij(3) = (GAMMA*DdataAtEdge(3,1,idx)-G2*ru2i)*ui-&
                 (GAMMA*DdataAtEdge(3,2,idx)-G2*ru2j)*uj

      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * DmatrixCoeffsAtEdge(1,1,idx)*dF_ij
      DfluxesAtEdge(:,2,idx) = -dscale * DmatrixCoeffsAtEdge(1,2,idx)*dF_ij
#endif

    end do

  end subroutine eulerlagrange_calcFluxGal1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxGalNoBdr1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IdofsAtEdge, dscale, DfluxesAtEdge, rcollection)

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
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

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
    real(DP), dimension(NVAR1D) :: dF_ij
    real(DP) :: ui,uj,ru2i,ru2j
    integer :: idx


    do idx = 1, size(DfluxesAtEdge,3)

      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      ! For a detailed description of algorithm and the definition of auxiliary
      ! quantities have a look at the subroutine "eulerlagrange_calcFluxGalerkin1d".
      !-------------------------------------------------------------------------
      
      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      
      ! Compute auxiliary variables
      ru2i = ui*DdataAtEdge(2,1,idx)
      ru2j = uj*DdataAtEdge(2,2,idx)
      
      ! Compute flux difference for x-direction
      dF_ij(1) = DdataAtEdge(2,1,idx) - DdataAtEdge(2,2,idx)
      dF_ij(2) = (G1*DdataAtEdge(3,1,idx)-G14*ru2i)-&
                 (G1*DdataAtEdge(3,2,idx)-G14*ru2j)
      dF_ij(3) = (GAMMA*DdataAtEdge(3,1,idx)-G2*ru2i)*ui-&
                 (GAMMA*DdataAtEdge(3,2,idx)-G2*ru2j)*uj

      ! Assemble symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,1,idx)-&
                                         DmatrixCoeffsAtEdge(1,2,idx))/2._DP*dF_ij
      DfluxesAtEdge(:,2,idx) = DfluxesAtEdge(:,1,idx)

    end do

  end subroutine eulerlagrange_calcFluxGalNoBdr1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxScDiss1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IdofsAtEdge, dscale, DfluxesAtEdge, rcollection)

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
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

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
#ifdef USE_EULER_IBP
    real(DP), dimension(NVAR1D) :: dF_i, dF_j
#else
    real(DP), dimension(NVAR1D) :: dF_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: ui,uj,ru2i,ru2j
    real(DP) :: d_ij,hi,hj,H_ij,q_ij,u_ij,aux,vel,cs
    integer :: idx

    do idx = 1, size(DfluxesAtEdge,3)
          
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      ! For a detailed description of algorithm and the definition of auxiliary
      ! quantities have a look at the subroutine "eulerlagrange_calcFluxGalerkin1d".
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      
      ! Compute auxiliary variables
      ru2i = ui*DdataAtEdge(2,1,idx)
      ru2j = uj*DdataAtEdge(2,2,idx)

#ifdef USE_EULER_IBP
      ! Compute fluxes for x-direction
      dF_i(1) = DdataAtEdge(2,1,idx)
      dF_i(2) = G1*DdataAtEdge(3,1,idx)-G14*ru2i
      dF_i(3) = (GAMMA*DdataAtEdge(3,1,idx)-G2*ru2i)*ui
      
      dF_j(1) = DdataAtEdge(2,2,idx)
      dF_j(2) = G1*DdataAtEdge(3,2,idx)-G14*ru2j
      dF_j(3) = (GAMMA*DdataAtEdge(3,2,idx)-G2*ru2j)*uj
#else
      ! Compute flux difference for x-direction
      dF_ij(1) = DdataAtEdge(2,1,idx) - DdataAtEdge(2,2,idx)
      dF_ij(2) = (G1*DdataAtEdge(3,1,idx)-G14*ru2i)-&
                 (G1*DdataAtEdge(3,2,idx)-G14*ru2j)
      dF_ij(3) = (GAMMA*DdataAtEdge(3,1,idx)-G2*ru2i)*ui-&
                 (GAMMA*DdataAtEdge(3,2,idx)-G2*ru2j)*uj
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation proportional to the spaectral radius
      !-------------------------------------------------------------------------

      ! Compute skew-symmetric coefficient
      a = 0.5_DP*(DmatrixCoeffsAtEdge(1,1,idx)-&
                  DmatrixCoeffsAtEdge(1,2,idx))

      ! Compute Roe mean values
      aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
      u_ij = (aux*ui+uj)/(aux+1.0_DP)
      hi   = GAMMA*DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui)
      hj   = GAMMA*DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj)
      H_ij = (aux*hi+hj)/(aux+1.0_DP)
      
      ! Compute auxiliary variables
      aux  = abs(a(1)) ! = sqrt(a(1)*a(1))
      vel  = u_ij*a(1)
      q_ij = 0.5_DP*(u_ij*u_ij)
      cs   = sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL))
      
      ! Scalar dissipation
      d_ij = abs(vel) + aux*cs
      
      ! Multiply the solution difference by the artificial diffusion factor
      Diff = d_ij*(DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx))

      !-------------------------------------------------------------------------
      ! Build both contributions into the fluxes
      !-------------------------------------------------------------------------

#ifdef USE_EULER_IBP
      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)

#else
      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_ij + Diff)
#endif

    end do

  end subroutine eulerlagrange_calcFluxScDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxRoeDiss1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IdofsAtEdge, dscale, DfluxesAtEdge, rcollection)

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
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

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
#ifdef USE_EULER_IBP
    real(DP), dimension(NVAR1D) :: dF_i, dF_j
#else
    real(DP), dimension(NVAR1D) :: dF_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: ui,uj,ru2i,ru2j,b1,b2
    real(DP) :: aux,uPow2,hi,hj,H_ij,q_ij,u_ij
    real(DP) :: anorm,l1,l2,l3,w1,w2,w3,cPow2,cs
    integer :: idx
    
    
    do idx = 1, size(DfluxesAtEdge,3)
          
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      ! For a detailed description of algorithm and the definition of auxiliary
      ! quantities have a look at the subroutine "eulerlagrange_calcFluxGalerkin1d".
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      
      ! Compute auxiliary variables
      ru2i = ui*DdataAtEdge(2,1,idx)
      ru2j = uj*DdataAtEdge(2,2,idx)

#ifdef USE_EULER_IBP
      ! Compute fluxes for x-direction
      dF_i(1) = DdataAtEdge(2,1,idx)
      dF_i(2) = G1*DdataAtEdge(3,1,idx)-G14*ru2i
      dF_i(3) = (GAMMA*DdataAtEdge(3,1,idx)-G2*ru2i)*ui
      
      dF_j(1) = DdataAtEdge(2,2,idx)
      dF_j(2) = G1*DdataAtEdge(3,2,idx)-G14*ru2j
      dF_j(3) = (GAMMA*DdataAtEdge(3,2,idx)-G2*ru2j)*uj
#else
      ! Compute flux difference for x-direction
      dF_ij(1) = DdataAtEdge(2,1,idx) - DdataAtEdge(2,2,idx)
      dF_ij(2) = (G1*DdataAtEdge(3,1,idx)-G14*ru2i)-&
                 (G1*DdataAtEdge(3,2,idx)-G14*ru2j)
      dF_ij(3) = (GAMMA*DdataAtEdge(3,1,idx)-G2*ru2i)*ui-&
                 (GAMMA*DdataAtEdge(3,2,idx)-G2*ru2j)*uj
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the dissipation tensor by Roe
      !-------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient
      a = 0.5_DP*(DmatrixCoeffsAtEdge(1,1,idx)-&
                  DmatrixCoeffsAtEdge(1,2,idx))
      anorm = abs(a(1))
      
      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui)
        hj   = GAMMA*DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)
        
        ! Compute auxiliary variables
        uPow2 = u_ij*u_ij
        q_ij  = 0.5_DP*uPow2
        cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
        cs    = sqrt(cPow2)
        
        ! Compute eigenvalues
        l1 = abs(u_ij-cs)
        l2 = abs(u_ij)
        l3 = abs(u_ij+cs)
        
        ! Compute solution difference U_j-U_i
        Diff = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
        
        ! Compute auxiliary quantities for characteristic variables
        b2 = G1/cPow2; b1 = b2*q_ij
        
        ! Compute characteristic variables multiplied by the
        ! corresponding eigenvalue
        w1 = l1 * 0.5_DP * (       (b1+u_ij/cs)*Diff(1) -&
                            (b2*u_ij+1.0_DP/cs)*Diff(2) +&
                                             b2*Diff(3) )
        w2 = l2 *          (             (1-b1)*Diff(1) +&
                                        b2*u_ij*Diff(2) -&
                                             b2*Diff(3) )
        w3 = l3 * 0.5_DP * (       (b1-u_ij/cs)*Diff(1) -&
                            (b2*u_ij-1.0_DP/cs)*Diff(2) +&
                                             b2*Diff(3) )
        
        ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
        Diff(1) = anorm * ( w1 + w2 + w3 )
        Diff(2) = anorm * ( (u_ij-cs)*w1 + u_ij*w2 + (u_ij+cs)*w3 )
        Diff(3) = anorm * ( (H_ij-u_ij*cs)*w1 + q_ij*w2 + (H_ij+u_ij*cs)*w3 )
        
        !-----------------------------------------------------------------------
        ! Build both contributions into the fluxes
        !-------------------------------------------------------------------------

#ifdef USE_EULER_IBP
        ! Assemble skew-symmetric fluxes
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF_i + Diff)
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
        
#else
        ! Assemble fluxes
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF_ij + Diff)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_ij + Diff)
#endif

      else

#ifdef USE_EULER_IBP
        ! Assemble skew-symmetric fluxes
        DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_j-&
                                           DmatrixCoeffsAtEdge(1,1,idx)*dF_i )
        DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)
        
#else
        ! Assemble fluxes
        DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF_ij)
        DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_ij)
#endif

      end if

    end do

  end subroutine eulerlagrange_calcFluxRoeDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxRusDiss1d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IdofsAtEdge, dscale, DfluxesAtEdge, rcollection)


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
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

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
#ifdef USE_EULER_IBP
    real(DP), dimension(NVAR1D) :: dF_i, dF_j
#else
    real(DP), dimension(NVAR1D) :: dF_ij
#endif
    real(DP), dimension(NVAR1D) :: Diff
    real(DP) :: ui,uj,ru2i,ru2j
    real(DP) :: d_ij,ci,cj,Ei,Ej
    integer :: idx

    
    do idx = 1, size(DfluxesAtEdge,3)
      
      !-------------------------------------------------------------------------
      ! Evaluate the Galerkin fluxes
      ! For a detailed description of algorithm and the definition of auxiliary
      ! quantities have a look at the subroutine "eulerlagrange_calcFluxGalerkin1d".
      !-------------------------------------------------------------------------

      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      Ei = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      Ej = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)

      ! Compute auxiliary variables
      ru2i = ui*DdataAtEdge(2,1,idx)
      ru2j = uj*DdataAtEdge(2,2,idx)
      
#ifdef USE_EULER_IBP
      ! Compute fluxes for x-direction
      dF_i(1) = DdataAtEdge(2,1,idx)
      dF_i(2) = G1*DdataAtEdge(3,1,idx)-G14*ru2i
      dF_i(3) = (GAMMA*DdataAtEdge(3,1,idx)-G2*ru2i)*ui
      
      dF_j(1) = DdataAtEdge(2,2,idx)
      dF_j(2) = G1*DdataAtEdge(3,2,idx)-G14*ru2j
      dF_j(3) = (GAMMA*DdataAtEdge(3,2,idx)-G2*ru2j)*uj
#else
      ! Compute flux difference for x-direction
      dF_ij(1) = DdataAtEdge(2,1,idx) - DdataAtEdge(2,2,idx)
      dF_ij(2) = (G1*DdataAtEdge(3,1,idx)-G14*ru2i)-&
                 (G1*DdataAtEdge(3,2,idx)-G14*ru2j)
      dF_ij(3) = (GAMMA*DdataAtEdge(3,1,idx)-G2*ru2i)*ui-&
                 (GAMMA*DdataAtEdge(3,2,idx)-G2*ru2j)*uj
#endif

      !-------------------------------------------------------------------------
      ! Evaluate the scalar dissipation of Rusanov-type
      !-------------------------------------------------------------------------
      
      ! Compute the speed of sound
      ci = sqrt(max(G15*(Ei-0.5_DP*ui*ui), SYS_EPSREAL))
      cj = sqrt(max(G15*(Ej-0.5_DP*uj*uj), SYS_EPSREAL))
      
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

#ifdef USE_EULER_IBP
      ! Assemble skew-symmetric fluxes
      DfluxesAtEdge(:,1,idx) = dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_j-&
                                         DmatrixCoeffsAtEdge(1,1,idx)*dF_i + Diff)
      DfluxesAtEdge(:,2,idx) = -DfluxesAtEdge(:,1,idx)

#else
      ! Assemble fluxes
      DfluxesAtEdge(:,1,idx) =  dscale * (DmatrixCoeffsAtEdge(1,1,idx)*dF_ij + Diff)
      DfluxesAtEdge(:,2,idx) = -dscale * (DmatrixCoeffsAtEdge(1,2,idx)*dF_ij + Diff)
#endif

    end do

  end subroutine eulerlagrange_calcFluxRusDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatDiagMatD1d_sim(DdataAtNode, DmatrixCoeffsAtNode,&
      IdofsAtNode, dscale, DcoefficientsAtNode, rcollection)

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
  integer, dimension(:,:), intent(in) :: IdofsAtNode

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
      
      ! Compute auxiliary variables
      ui = DdataAtNode(2,inode)/DdataAtNode(1,inode)
      
      ! Compute Galerkin coefficient K_ii
      DcoefficientsAtNode(1,1,inode) = 0.0_DP
      DcoefficientsAtNode(2,1,inode) = dscale * G13*ui*DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(3,1,inode) = dscale * GAMMA*ui*DmatrixCoeffsAtNode(1,inode)
    end do

  end subroutine eulerlagrange_calcMatDiagMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatDiag1d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IdofsAtNode, dscale,&
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
  integer, dimension(:,:), intent(in) :: IdofsAtNode

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
    real(DP) :: ui,Ei,uPow2i
    integer :: inode

    do inode = 1, size(DcoefficientsAtNode,3)
      
      ! Compute auxiliary variables
      ui = DdataAtNode(2,inode)/DdataAtNode(1,inode)
      Ei = DdataAtNode(3,inode)/DdataAtNode(1,inode)
      uPow2i = ui*ui
      
      ! Compute Galerkin coefficient K_ii
      DcoefficientsAtNode(1,1,inode) = 0.0_DP
      DcoefficientsAtNode(2,1,inode) = dscale * G14*uPow2i*DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(3,1,inode) = dscale * (G1*uPow2i-GAMMA*Ei)*ui*DmatrixCoeffsAtNode(1,inode)
      
      DcoefficientsAtNode(4,1,inode) = dscale * DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(5,1,inode) = dscale * G13*ui*DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(6,1,inode) = dscale * (GAMMA*Ei-G16*uPow2i)*DmatrixCoeffsAtNode(1,inode)
      
      DcoefficientsAtNode(7,1,inode) = 0.0_DP
      DcoefficientsAtNode(8,1,inode) = dscale * G1*DmatrixCoeffsAtNode(1,inode)
      DcoefficientsAtNode(9,1,inode) = dscale * GAMMA*ui*DmatrixCoeffsAtNode(1,inode)
    end do
    
  end subroutine eulerlagrange_calcMatDiag1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatGalMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
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
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

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
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0_DP

      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge(1,2,idx) = 0.0_DP
      DcoefficientsAtEdge(2,2,idx) = dscale * G13*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge(1,3,idx) = 0.0_DP
      DcoefficientsAtEdge(2,3,idx) = dscale * G13*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)
    end do

  end subroutine eulerlagrange_calcMatGalMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatGal1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
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
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

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
    real(DP) :: ui,uj,Ei,Ej,uPow2i,uPow2j
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      Ei = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uPow2i = ui*ui
      
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      Ej = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      uPow2j = uj*uj
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0_DP

      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge(1,2,idx) = 0.0_DP
      DcoefficientsAtEdge(2,2,idx) = dscale * G14*uPow2j*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * (G1*uPow2j-GAMMA*Ej)*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(4,2,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(5,2,idx) = dscale * G13*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(6,2,idx) = dscale * (GAMMA*Ej-G16*uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(7,2,idx) = 0.0_DP
      DcoefficientsAtEdge(8,2,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(9,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge(1,3,idx) = 0.0_DP
      DcoefficientsAtEdge(2,3,idx) = dscale * G14*uPow2i*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * (G1*uPow2i-GAMMA*Ei)*ui*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(4,3,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(5,3,idx) = dscale * G13*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(6,3,idx) = dscale * (GAMMA*Ei-G16*uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(7,3,idx) = 0.0_DP
      DcoefficientsAtEdge(8,3,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(9,3,idx) = dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)
    end do

  end subroutine eulerlagrange_calcMatGal1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatScDissMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
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
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

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
    real(DP) :: anorm,aux,hi,hj,H_ij,q_ij,ui,uj,u_ij
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
      
      ! Compute auxiliary variables
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)

      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge(1,2,idx) = 0.0_DP
      DcoefficientsAtEdge(2,2,idx) = dscale * G13*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge(1,3,idx) = 0.0_DP
      DcoefficientsAtEdge(2,3,idx) = dscale * G13*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)
      
      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient and its norm
      a = 0.5_DP*(DmatrixCoeffsAtEdge(1,2,idx)-&
                  DmatrixCoeffsAtEdge(1,1,idx))
      anorm = abs(a(1))
      
      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui)
        hj   = GAMMA*DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)
        
        ! Compute auxiliary values
        q_ij = 0.5_DP*u_ij*u_ij
        
        ! Compute scalar dissipation
        DcoefficientsAtEdge(:,1,idx) = -dscale * (abs(a(1)*u_ij) +&
            anorm*sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL)))
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0_DP
        
      end if
    end do

  end subroutine eulerlagrange_calcMatScDissMatD1d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatScDiss1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
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
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

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
    real(DP) :: anorm,aux,hi,hj,Ei,Ej,H_ij,q_ij,ui,uj,u_ij,uPow2i,uPow2j
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)
    
      ! Compute auxiliary variables
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      Ei = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uPow2i = ui*ui
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      Ej = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      uPow2j = uj*uj
      
      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge(1,2,idx) = 0.0_DP
      DcoefficientsAtEdge(2,2,idx) = dscale * G14*uPow2j*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * (G1*uPow2j-GAMMA*Ej)*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(4,2,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(5,2,idx) = dscale * G13*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(6,2,idx) = dscale * (GAMMA*Ej-G16*uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(7,2,idx) = 0.0_DP
      DcoefficientsAtEdge(8,2,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(9,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge(1,3,idx) = 0.0_DP
      DcoefficientsAtEdge(2,3,idx) = dscale * G14*uPow2i*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * (G1*uPow2i-GAMMA*Ei)*ui*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(4,3,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(5,3,idx) = dscale * G13*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(6,3,idx) = dscale * (GAMMA*Ei-G16*uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(7,3,idx) = 0.0_DP
      DcoefficientsAtEdge(8,3,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(9,3,idx) = dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient and its norm
      a = 0.5_DP*(DmatrixCoeffsAtEdge(1,2,idx)-&
                  DmatrixCoeffsAtEdge(1,1,idx))
      anorm = abs(a(1))
      
      ! Nullify dissipation tensor
      DcoefficientsAtEdge(:,1,idx) = 0.0_DP

      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui)
        hj   = GAMMA*DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)

        ! Compute auxiliary values
        q_ij = 0.5_DP*u_ij*u_ij
        
        ! Compute scalar dissipation
        aux = -dscale * (abs(a(1)*u_ij) +&
            anorm*sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL)))
        
        DcoefficientsAtEdge(1,1,idx) = aux
        DcoefficientsAtEdge(5,1,idx) = aux
        DcoefficientsAtEdge(9,1,idx) = aux
      end if
    end do
    
  end subroutine eulerlagrange_calcMatScDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatRoeDissMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
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
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

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
    real(DP) :: aux,hi,hj,H_ij,q_ij,ui,uj,u_ij
    real(DP) :: l1,l2,l3,anorm,cs,cPow2,b1,b2
    integer :: idx
    
    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)

      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge(1,2,idx) = 0.0_DP
      DcoefficientsAtEdge(2,2,idx) = dscale * G13*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge(1,3,idx) = 0.0_DP
      DcoefficientsAtEdge(2,3,idx) = dscale * G13*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient and its norm
      a = 0.5_DP*(DmatrixCoeffsAtEdge(1,2,idx)-&
                  DmatrixCoeffsAtEdge(1,1,idx))
      anorm = abs(a(1))
      
      if (anorm .gt. SYS_EPSREAL) then
        
        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui)
        hj   = GAMMA*DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)
        
        ! Compute auxiliary values
        q_ij  = 0.5_DP*u_ij*u_ij
        cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
        cs    = sqrt(cPow2)
        
        b2    = G1/cPow2
        b1    = b2*q_ij

        ! Diagonal matrix of eigenvalues
        l1 = abs(u_ij-cs)
        l2 = abs(u_ij)
        l3 = abs(u_ij+cs)
        
        ! Matrix of right eigenvectors
        R_ij(1,1) =  l1
        R_ij(2,1) =  l1*(u_ij-cs)
        R_ij(3,1) =  l1*(H_ij-cs*u_ij)
        
        R_ij(1,2) =  l2
        R_ij(2,2) =  l2*u_ij
        R_ij(3,2) =  l2*q_ij
        
        R_ij(1,3) =  l3
        R_ij(2,3) =  l3*(u_ij+cs)
        R_ij(3,3) =  l3*(H_ij+cs*u_ij)
        
        ! Matrix of left eigenvectors
        L_ij(1,1) = 0.5_DP*(b1+u_ij/cs)
        L_ij(2,1) = 1.0_DP-b1
        L_ij(3,1) = 0.5_DP*(b1-u_ij/cs)
        
        L_ij(1,2) = -0.5_DP*(b2*u_ij+1/cs)
        L_ij(2,2) =  b2*u_ij
        L_ij(3,2) = -0.5_DP*(b2*u_ij-1/cs)
        
        L_ij(1,3) = 0.5_DP*b2
        L_ij(2,3) = -b2
        L_ij(3,3) = 0.5_DP*b2

        ! Include scaling parameter
        anorm = -dscale*anorm

        ! Compute tensorial dissipation D_ij = diag(R_ij*|Lbd_ij|*L_ij)*I
        DcoefficientsAtEdge(1,1,idx) = anorm*( R_ij(1,1)*L_ij(1,1)+&
            R_ij(1,2)*L_ij(2,1)+R_ij(1,3)*L_ij(3,1)  )
        DcoefficientsAtEdge(2,1,idx) = anorm*( R_ij(2,1)*L_ij(1,2)+&
            R_ij(2,2)*L_ij(2,2)+R_ij(2,3)*L_ij(3,2)  )
        DcoefficientsAtEdge(3,1,idx) = anorm*( R_ij(3,1)*L_ij(1,3)+&
            R_ij(3,2)*L_ij(2,3)+R_ij(3,3)*L_ij(3,3)  )
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0_DP

      end if
    end do

  end subroutine eulerlagrange_calcMatRoeDissMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatRoeDiss1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
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
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

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
    real(DP) :: aux,Ei,Ej,hi,hj,H_ij,q_ij,ui,uj,u_ij
    real(DP) :: l1,l2,l3,anorm,cs,cPow2,b1,b2,uPow2i,uPow2j
    integer :: idx,i,j,k

    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      Ei = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uPow2i = ui*ui
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      Ej = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      uPow2j = uj*uj
      
      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge(1,2,idx) = 0.0_DP
      DcoefficientsAtEdge(2,2,idx) = dscale * G14*uPow2j*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * (G1*uPow2j-GAMMA*Ej)*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(4,2,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(5,2,idx) = dscale * G13*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(6,2,idx) = dscale * (GAMMA*Ej-G16*uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(7,2,idx) = 0.0_DP
      DcoefficientsAtEdge(8,2,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(9,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge(1,3,idx) = 0.0_DP
      DcoefficientsAtEdge(2,3,idx) = dscale * G14*uPow2i*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * (G1*uPow2i-GAMMA*Ei)*ui*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(4,3,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(5,3,idx) = dscale * G13*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(6,3,idx) = dscale * (GAMMA*Ei-G16*uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(7,3,idx) = 0.0_DP
      DcoefficientsAtEdge(8,3,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(9,3,idx) = dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute skew-symmetric coefficient and its norm
      a = 0.5_DP*(DmatrixCoeffsAtEdge(1,2,idx)-&
                  DmatrixCoeffsAtEdge(1,1,idx))
      anorm = abs(a(1))

      if (anorm .gt. SYS_EPSREAL) then

        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui)
        hj   = GAMMA*DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)

        ! Compute auxiliary values
        q_ij  = 0.5_DP*u_ij*u_ij
        cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
        cs    = sqrt(cPow2)
        b2    = G1/cPow2
        b1    = b2*q_ij
        
        ! Diagonal matrix of eigenvalues
        l1 = abs(u_ij-cs)
        l2 = abs(u_ij)
        l3 = abs(u_ij+cs)
        
        ! Matrix of right eigenvectors
        R_ij(1,1) =  l1
        R_ij(2,1) =  l1*(u_ij-cs)
        R_ij(3,1) =  l1*(H_ij-cs*u_ij)
        
        R_ij(1,2) =  l2
        R_ij(2,2) =  l2*u_ij
        R_ij(3,2) =  l2*q_ij
        
        R_ij(1,3) =  l3
        R_ij(2,3) =  l3*(u_ij+cs)
        R_ij(3,3) =  l3*(H_ij+cs*u_ij)
        
        ! Matrix of left eigenvectors
        L_ij(1,1) = 0.5_DP*(b1+u_ij/cs)
        L_ij(2,1) = 1.0_DP-b1
        L_ij(3,1) = 0.5_DP*(b1-u_ij/cs)
        
        L_ij(1,2) = -0.5_DP*(b2*u_ij+1/cs)
        L_ij(2,2) =  b2*u_ij
        L_ij(3,2) = -0.5_DP*(b2*u_ij-1/cs)
        
        L_ij(1,3) = 0.5_DP*b2
        L_ij(2,3) = -b2
        L_ij(3,3) = 0.5_DP*b2
        
        ! Include scaling parameter
        anorm = -dscale*anorm

        ! Compute tensorial dissipation D_ij = R_ij*|Lbd_ij|*L_ij
        do i = 1, NVAR1D
          do j = 1, NVAR1D
            aux = 0.0_DP
            do k = 1, NVAR1D
              aux = aux + R_ij(i,k)*L_ij(k,j)
            end do
            DcoefficientsAtEdge(NVAR1D*(j-1)+i,1,idx) = anorm*aux
          end do
        end do
        
      else
        
        ! Nullify dissipation tensor
        DcoefficientsAtEdge(:,1,idx) = 0.0_DP
        
      end if
    end do

  end subroutine eulerlagrange_calcMatRoeDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatRusDissMatD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
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
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

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
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      Ei = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      Ej = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      
      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge(1,2,idx) = 0.0_DP
      DcoefficientsAtEdge(2,2,idx) = dscale * G13*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge(1,3,idx) = 0.0_DP
      DcoefficientsAtEdge(2,3,idx) = dscale * G13*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
      ci = sqrt(max(G15*(Ei-0.5_DP*ui*ui), SYS_EPSREAL))
      cj = sqrt(max(G15*(Ej-0.5_DP*uj*uj), SYS_EPSREAL))
      
      ! Compute dissipation tensor D_ij
      DcoefficientsAtEdge(:,1,idx) = -dscale *&
          max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj)+&
               abs(DmatrixCoeffsAtEdge(1,1,idx))*cj,&
               abs(DmatrixCoeffsAtEdge(1,2,idx)*ui)+&
               abs(DmatrixCoeffsAtEdge(1,2,idx))*ci )
    end do
    
  end subroutine eulerlagrange_calcMatRusDissMatD1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatRusDiss1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
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
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

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
    real(DP) :: ui,uj,ci,cj,Ei,Ej,uPow2i,uPow2j,aux
    integer :: idx

    do idx = 1, size(DcoefficientsAtEdge,3)

      ! Compute auxiliary variables
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      Ei = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      uPow2i = ui*ui
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      Ej = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      uPow2j = uj*uj
      
      ! Compute Galerkin coefficient K_ij
      DcoefficientsAtEdge(1,2,idx) = 0.0_DP
      DcoefficientsAtEdge(2,2,idx) = dscale * G14*uPow2j*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(3,2,idx) = dscale * (G1*uPow2j-GAMMA*Ej)*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(4,2,idx) = dscale * DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(5,2,idx) = dscale * G13*uj*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(6,2,idx) = dscale * (GAMMA*Ej-G16*uPow2j)*DmatrixCoeffsAtEdge(1,1,idx)
      
      DcoefficientsAtEdge(7,2,idx) = 0.0_DP
      DcoefficientsAtEdge(8,2,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,1,idx)
      DcoefficientsAtEdge(9,2,idx) = dscale * GAMMA*uj*DmatrixCoeffsAtEdge(1,1,idx)
      
      ! Compute Galerkin coefficient K_ji
      DcoefficientsAtEdge(1,3,idx) = 0.0_DP
      DcoefficientsAtEdge(2,3,idx) = dscale * G14*uPow2i*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(3,3,idx) = dscale * (G1*uPow2i-GAMMA*Ei)*ui*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(4,3,idx) = dscale * DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(5,3,idx) = dscale * G13*ui*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(6,3,idx) = dscale * (GAMMA*Ei-G16*uPow2i)*DmatrixCoeffsAtEdge(1,2,idx)
      
      DcoefficientsAtEdge(7,3,idx) = 0.0_DP
      DcoefficientsAtEdge(8,3,idx) = dscale * G1*DmatrixCoeffsAtEdge(1,2,idx)
      DcoefficientsAtEdge(9,3,idx) = dscale * GAMMA*ui*DmatrixCoeffsAtEdge(1,2,idx)

      !---------------------------------------------------------------------------
      ! Evaluate the dissipation
      !---------------------------------------------------------------------------
      
      ! Compute the speed of sound
      ci = sqrt(max(G15*(Ei-0.5_DP*ui*ui), SYS_EPSREAL))
      cj = sqrt(max(G15*(Ej-0.5_DP*uj*uj), SYS_EPSREAL))

      ! Compute dissipation tensor D_ij
      aux = -dscale * max( abs(DmatrixCoeffsAtEdge(1,1,idx)*uj)+&
                           abs(DmatrixCoeffsAtEdge(1,1,idx))*cj,&
                           abs(DmatrixCoeffsAtEdge(1,2,idx)*ui)+&
                           abs(DmatrixCoeffsAtEdge(1,2,idx))*ci )

      DcoefficientsAtEdge(:,1,idx) = 0.0_DP
      DcoefficientsAtEdge(1,1,idx) = aux
      DcoefficientsAtEdge(5,1,idx) = aux
      DcoefficientsAtEdge(9,1,idx) = aux
    end do
    
  end subroutine eulerlagrange_calcMatRusDiss1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcCharacteristics1d_sim(Dweight, DdataAtEdge,&
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
    real(DP) :: u_ij,H_ij,q_ij,cs,aux,aux1,aux2,ui,uj,hi,hj,cPow2,anorm,b1,b2
    integer :: idx

    ! Compute norm of weighting coefficient
    anorm = abs(Dweight(1))

    ! Check if weighting coefficient is zero
    if (anorm .le. SYS_EPSREAL) then
      if (present(DcharVariablesAtEdge))     DcharVariablesAtEdge     = 0.0_DP
      if (present(DeigenvaluesAtEdge))       DeigenvaluesAtEdge       = 0.0_DP
      if (present(DrightEigenvectorsAtEdge)) DrightEigenvectorsAtEdge = 0.0_DP
      if (present(DleftEigenvectorsAtEdge))  DleftEigenvectorsAtEdge  = 0.0_DP

      ! That's it
      return
    end if
    
    
    ! Do we have to compute characteristic variables
    if (present(DcharVariablesAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)
        
        ! Compute velocities
        ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
        uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)

        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui)
        hj   = GAMMA*DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)
        
        ! Compute auxiliary variables
        q_ij  = 0.5_DP*(u_ij*u_ij)
        cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
        cs    = sqrt(cPow2)
        b2    = G1/cPow2
        b1    = b2*q_ij
        
        ! Compute solution difference U_j-U_i
        Diff = DdataAtEdge(:,2,idx)-DdataAtEdge(:,1,idx)
          
        ! Compute characteristic variables
        DcharVariablesAtEdge(1,idx) = anorm * 0.5_DP *&
                                     (       (b1+u_ij/cs)*Diff(1)-&
                                      (b2*u_ij+1.0_DP/cs)*Diff(2)+&
                                                       b2*Diff(3))
        DcharVariablesAtEdge(2,idx) = anorm *&
                                     (             (1-b1)*Diff(1)+&
                                                  b2*u_ij*Diff(2)-&
                                                        b2*Diff(3) )
        DcharVariablesAtEdge(3,idx) = anorm * 0.5_DP *&
                                     (       (b1-u_ij/cs)*Diff(1)-&
                                      (b2*u_ij-1.0_DP/cs)*Diff(2)+&
                                                       b2*Diff(3) )
      end do
    end if


    ! Do we have to compute eigenvalues
    if (present(DeigenvaluesAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)
        
        ! Compute velocities
        ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
        uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)

        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui)
        hj   = GAMMA*DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)
        
        ! Compute auxiliary variable
        cs = sqrt(max(G1*(H_ij-0.5_DP*(u_ij*u_ij)), SYS_EPSREAL))
       
        ! Compute eigenvalues
        DeigenvaluesAtEdge(1,idx) = u_ij-cs
        DeigenvaluesAtEdge(2,idx) = u_ij
        DeigenvaluesAtEdge(3,idx) = u_ij+cs
      end do
    end if


    ! Do we have to compute right eigenvectors
    if (present(DrightEigenvectorsAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)
        
        ! Compute velocities
        ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
        uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)

        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui)
        hj   = GAMMA*DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)

        ! Compute auxiliary variables
        q_ij  = 0.5_DP*(u_ij*u_ij)
        cs    = sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL))

        ! Compute right eigenvectors
        DrightEigenvectorsAtEdge(1,idx) =  1.0_DP
        DrightEigenvectorsAtEdge(2,idx) =  u_ij-cs
        DrightEigenvectorsAtEdge(3,idx) =  H_ij-u_ij*cs

        DrightEigenvectorsAtEdge(4,idx) =  1.0_DP
        DrightEigenvectorsAtEdge(5,idx) =  u_ij
        DrightEigenvectorsAtEdge(6,idx) =  q_ij

        DrightEigenvectorsAtEdge(7,idx) =  1.0_DP
        DrightEigenvectorsAtEdge(8,idx) =  u_ij+cs
        DrightEigenvectorsAtEdge(9,idx) =  H_ij+u_ij*cs
      end do
    end if


    ! Do we have to compute left eigenvectors
    if (present(DleftEigenvectorsAtEdge)) then
      do idx = 1, size(DdataAtEdge,3)
        
        ! Compute velocities
        ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
        uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)

        ! Compute Roe mean values
        aux  = sqrt(max(DdataAtEdge(1,1,idx)/DdataAtEdge(1,2,idx), SYS_EPSREAL))
        u_ij = (aux*ui+uj)/(aux+1.0_DP)
        hi   = GAMMA*DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)-G2*(ui*ui)
        hj   = GAMMA*DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)-G2*(uj*uj)
        H_ij = (aux*hi+hj)/(aux+1.0_DP)

        ! Compute auxiliary variables
        q_ij  = 0.5_DP*(u_ij*u_ij)
        cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
        cs    = sqrt(cPow2)
        b2    = G1/cPow2
        b1    = b2*q_ij

        ! Compute left eigenvectors
        DleftEigenvectorsAtEdge(1,idx) = 0.5_DP * (b1+u_ij/cs)
        DleftEigenvectorsAtEdge(2,idx) =          (1-b1)
        DleftEigenvectorsAtEdge(3,idx) = 0.5_DP * (b1-u_ij/cs)

        DleftEigenvectorsAtEdge(4,idx) =-0.5_DP * (b2*u_ij+1.0_DP/cs)
        DleftEigenvectorsAtEdge(5,idx) =          (b2*u_ij)
        DleftEigenvectorsAtEdge(6,idx) =-0.5_DP * (b2*u_ij-1.0_DP/cs)

        DleftEigenvectorsAtEdge(7,idx) = 0.5_DP*b2
        DleftEigenvectorsAtEdge(8,idx) =       -b2
        DleftEigenvectorsAtEdge(9,idx) = 0.5_DP*b2
      end do
    end if
    
  end subroutine eulerlagrange_calcCharacteristics1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxFCTScalarDiss1d(&
      U1_i, U1_j, U2_i, U2_j, C_ij, C_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 1D using scalar dissipation.
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
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: ui,uj,d_ij,hi,hj,H_ij,q_ij,aux,vel,cs

    ! Compute velocities
    ui = U2_i(2)/U2_i(1); uj = U2_j(2)/U2_j(1)

    ! Compute skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji)

    ! Compute Roe mean values
    aux  = sqrt(max(U2_i(1)/U2_j(1), SYS_EPSREAL))
    hi   = GAMMA*U2_i(3)/U2_i(1)-G2*(U2_i(2)*U2_i(2))/(U2_i(1)*U2_i(1))
    hj   = GAMMA*U2_j(3)/U2_j(1)-G2*(U2_j(2)*U2_j(2))/(U2_j(1)*U2_j(1))
    H_ij = (aux*hi+hj)/(aux+1.0_DP)
    aux  = (aux*ui+uj)/(aux+1.0_DP)

    ! Compute auxiliary variables
    vel  = aux*a(1)
    q_ij = 0.5_DP*(aux*aux)
    cs   = sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL))

    ! Scalar dissipation
    d_ij = abs(vel) + abs(a(1))*cs

    ! Compute conservative fluxes
    F_ij = dscale1*(U1_i-U1_j) + dscale2*d_ij*(U2_i-U2_j)

  end subroutine eulerlagrange_calcFluxFCTScalarDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxFCTTensorDiss1d(&
      U1_i, U1_j, U2_i, U2_j, C_ij, C_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 1D using tensorial dissipation.
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
    real(DP), dimension(NVAR1D) :: Diff
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: ui,uj,b1,b2
    real(DP) :: aux,uPow2,hi,hj,H_ij,q_ij,u_ij
    real(DP) :: anorm,l1,l2,l3,w1,w2,w3,cPow2,cs

    ! Compute velocities
    ui = U2_i(2)/U2_i(1); uj = U2_j(2)/U2_j(1)

    ! Compute skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji); anorm = abs(a(1))

    ! Compute Roe mean values
    aux  = sqrt(max(U2_i(1)/U2_j(1), SYS_EPSREAL))
    u_ij = (aux*ui+uj)/(aux+1.0_DP)
    hi   = GAMMA*U2_i(3)/U2_i(1)-G2*(U2_i(2)*U2_i(2))/(U2_i(1)*U2_i(1))
    hj   = GAMMA*U2_j(3)/U2_j(1)-G2*(U2_j(2)*U2_j(2))/(U2_j(1)*U2_j(1))
    H_ij = (aux*hi+hj)/(aux+1.0_DP)

    ! Compute auxiliary variables
    uPow2 = u_ij*u_ij
    q_ij  = 0.5_DP*uPow2
    cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
    cs    = sqrt(cPow2)

    ! Compute eigenvalues
    l1 = abs(u_ij-cs)
    l2 = abs(u_ij)
    l3 = abs(u_ij+cs)

    ! Compute solution difference U2_i-U2_j
    Diff = U2_i-U2_j

    ! Compute auxiliary quantities for characteristic variables
    b2 = G1/cPow2; b1 = b2*q_ij

    ! Compute characteristic variables multiplied by the
    ! corresponding eigenvalue
    w1 = l1 * 0.5_DP * (       (b1+u_ij/cs)*Diff(1) -&
                        (b2*u_ij+1.0_DP/cs)*Diff(2) +&
                                         b2*Diff(3) )
    w2 = l2 *          (             (1-b1)*Diff(1) +&
                                    b2*u_ij*Diff(2) -&
                                         b2*Diff(3) )
    w3 = l3 * 0.5_DP * (       (b1-u_ij/cs)*Diff(1) -&
                        (b2*u_ij-1.0_DP/cs)*Diff(2) +&
                                         b2*Diff(3) )

    ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
    Diff(1) = w1 + w2 + w3
    Diff(2) = (u_ij-cs)*w1 + u_ij*w2 + (u_ij+cs)*w3
    Diff(3) = (H_ij-u_ij*cs)*w1 + q_ij*w2 + (H_ij+u_ij*cs)*w3

    ! Compute conservative fluxes
    F_ij = dscale1*(U1_i-U1_j) + dscale2*anorm*Diff

  end subroutine eulerlagrange_calcFluxFCTTensorDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxFCTRusanov1d(&
      U1_i, U1_j, U2_i, U2_j, C_ij, C_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 1D using the Rusanov dissipation.
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
    real(DP) :: d_ij,ui,uj,ci,cj,Ei,Ej

    ! Compute velocities and energy
    ui = U2_i(2)/U2_i(1); uj = U2_j(2)/U2_j(1)
    Ei = U2_i(3)/U2_i(1); Ej = U2_j(3)/U2_j(1)

    ! Compute the speed of sound
    ci = sqrt(max(G15*(Ei-0.5_DP*ui*ui), SYS_EPSREAL))
    cj = sqrt(max(G15*(Ej-0.5_DP*uj*uj), SYS_EPSREAL))

    ! Scalar dissipation
    d_ij = max( abs(C_ij(1)*uj) + abs(C_ij(1))*cj,&
                abs(C_ji(1)*ui) + abs(C_ji(1))*ci )

    ! Compute conservative fluxes
    F_ij = dscale1*(U1_i-U1_j) + dscale2*d_ij*(U2_i-U2_j)

  end subroutine eulerlagrange_calcFluxFCTRusanov1d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxDensity1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density in 1D
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

  end subroutine eulerlagrange_trafoFluxDensity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffDensity1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density in 1D
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

  end subroutine eulerlagrange_trafoDiffDensity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxEnergy1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the energy in 1D
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
      DtransformedFluxesAtEdge(1,1,idx) = DfluxesAtEdge(3,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-DfluxesAtEdge(3,idx)
    end do

  end subroutine eulerlagrange_trafoFluxEnergy1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffEnergy1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the energy in 1D
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
      DtransformedDataAtEdge(1,idx) = DdataAtEdge(3,2,idx)-DdataAtEdge(3,1,idx)
    end do

  end subroutine eulerlagrange_trafoDiffEnergy1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxPressure1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the pressure in 1D
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
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      
      ! Transformed pressure fluxes
      DtransformedFluxesAtEdge(1,1,idx) = G1*(0.5_DP*ui*ui*DfluxesAtEdge(1,idx)-&
                                          ui*DfluxesAtEdge(2,idx)+DfluxesAtEdge(3,idx))
      DtransformedFluxesAtEdge(1,2,idx) =-G1*(0.5_DP*uj*uj*DfluxesAtEdge(1,idx)-&
                                          uj*DfluxesAtEdge(2,idx)+DfluxesAtEdge(3,idx))
    end do
    
  end subroutine eulerlagrange_trafoFluxPressure1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffPressure1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the pressure in 1D
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
      pi = G1*(DdataAtEdge(3,1,idx)-&
          0.5_DP*DdataAtEdge(2,1,idx)*DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx))
      pj = G1*(DdataAtEdge(3,2,idx)-&
          0.5_DP*DdataAtEdge(2,2,idx)*DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx))

      ! Transformed pressure difference
      DtransformedDataAtEdge(1,idx) = pj-pi
    end do

  end subroutine eulerlagrange_trafoDiffPressure1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxVelocity1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the x-velocity
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
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)

      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(1,1,idx) = (DfluxesAtEdge(2,idx)-&
                                           ui*DfluxesAtEdge(1,idx))/DdataAtEdge(1,1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-(DfluxesAtEdge(2,idx)-&
                                           uj*DfluxesAtEdge(1,idx))/DdataAtEdge(1,2,idx)
    end do
    
  end subroutine eulerlagrange_trafoFluxVelocity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffVelocity1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the x-velocity
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
    end do

  end subroutine eulerlagrange_trafoDiffVelocity1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxMomentum1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the x-momentum
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
    end do
    
  end subroutine eulerlagrange_trafoFluxMomentum1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffMomentum1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the x-momentum
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
    end do
    
  end subroutine eulerlagrange_trafoDiffMomentum1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxDenEng1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 1D
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
      DtransformedFluxesAtEdge(2,1,idx) = DfluxesAtEdge(3,idx)
      DtransformedFluxesAtEdge(2,2,idx) =-DfluxesAtEdge(3,idx)
    end do

  end subroutine eulerlagrange_trafoFluxDenEng1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffDenEng1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 1D
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
      DtransformedDataAtEdge(2,idx) = DdataAtEdge(3,2,idx)-DdataAtEdge(3,1,idx)
    end do

  end subroutine eulerlagrange_trafoDiffDenEng1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxDenPre1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 1D
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
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) = DfluxesAtEdge(1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-DfluxesAtEdge(1,idx)

      ! Transformed pressure fluxes
      DtransformedFluxesAtEdge(2,1,idx) = G1*(0.5_DP*ui*ui*DfluxesAtEdge(1,idx)-&
                                          ui*DfluxesAtEdge(2,idx)+DfluxesAtEdge(3,idx))
      DtransformedFluxesAtEdge(2,2,idx) =-G1*(0.5_DP*uj*uj*DfluxesAtEdge(1,idx)-&
                                          uj*DfluxesAtEdge(2,idx)+DfluxesAtEdge(3,idx))
    end do

  end subroutine eulerlagrange_trafoFluxDenPre1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffDenPre1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 1D
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
      pi = G1*(DdataAtEdge(3,1,idx)-0.5_DP*&
          DdataAtEdge(2,1,idx)*DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx))
      pj = G1*(DdataAtEdge(3,2,idx)-0.5_DP*&
          DdataAtEdge(2,2,idx)*DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx))
      
      ! Transformed density difference
      DtransformedDataAtEdge(1,idx) = DdataAtEdge(1,2,idx)-DdataAtEdge(1,1,idx)
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(2,idx) = pj-pi
    end do

  end subroutine eulerlagrange_trafoDiffDenPre1d_sim
  
  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxDenPreVel1d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density, pressure and velocity in 1D
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
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) = DfluxesAtEdge(1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-DfluxesAtEdge(1,idx)

      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(2,1,idx) = (DfluxesAtEdge(2,idx)-&
                                           ui*DfluxesAtEdge(1,idx))/DdataAtEdge(1,1,idx)
      DtransformedFluxesAtEdge(2,2,idx) =-(DfluxesAtEdge(2,idx)-&
                                           uj*DfluxesAtEdge(1,idx))/DdataAtEdge(1,2,idx)

      ! Transformed pressure fluxes
      DtransformedFluxesAtEdge(3,1,idx) = G1*(0.5_DP*ui*ui*DfluxesAtEdge(1,idx)-&
                                          ui*DfluxesAtEdge(2,idx)+DfluxesAtEdge(3,idx))
      DtransformedFluxesAtEdge(3,2,idx) =-G1*(0.5_DP*uj*uj*DfluxesAtEdge(1,idx)-&
                                          uj*DfluxesAtEdge(2,idx)+DfluxesAtEdge(3,idx))
    end do
      
  end subroutine eulerlagrange_trafoFluxDenPreVel1d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffDenPreVel1d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density, pressure and velocity in 1D
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
      pi = G1*(DdataAtEdge(3,1,idx)-0.5_DP*&
               DdataAtEdge(2,1,idx)*DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx))
      pj = G1*(DdataAtEdge(3,2,idx)-0.5_DP*&
               DdataAtEdge(2,2,idx)*DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx))
      
      ! Transformed density difference
      DtransformedDataAtEdge(1,idx) = DdataAtEdge(1,2,idx)-DdataAtEdge(1,1,idx)
      
      ! Transformed velocity difference in x-direction
      DtransformedDataAtEdge(2,idx) = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)-&
                                      DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(3,idx) = pj-pi
    end do

  end subroutine eulerlagrange_trafoDiffDenPreVel1d_sim

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_calcBoundaryvalues1d(DbdrNormal, DpointNormal,&
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
    select case(ibdrCondType)
    case(BDR_EULERWALL,&
         BDR_RLXEULERWALL)
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
      p   = thdyn_pressure(GAMMA, E, rho, v1)
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
      if (2.0_DP*G9*c .le. vn_b-vn) then
        if (present(istatus)) then
          istatus = -ibdrCondType
          return
        else
          call output_line('Riemann solver failed due to vacuum',&
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcBoundaryvalues1d')
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
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcBoundaryvalues1d')
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
      Du(2) = rho*DpointNormal(1)*vn
      Du(3) = pstar/G1+0.5_DP*rho*(vn*vn)

    case(BDR_VISCOUSWALL)
      !-------------------------------------------------------------------------

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      E   = Du(3)/rho
      p   = thdyn_pressure(GAMMA, E, rho, v1)

      ! Update the solution vector and let vn:=0
      Du(2) = 0.0_DP
      Du(3) = p/G1

    case(BDR_SUPERINLET)
      !-------------------------------------------------------------------------

      ! The free stream primitive variables are Deval=[rho,v1,p]
      rho = DbdrValue(1)
      p   = DbdrValue(3)
      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*DbdrValue(2)

      ! Compute Riemann invariants based on the free stream values
      W(1) = vn-2*c/G1
      W(2) = vn+2*c/G1
      W(3) = p/(rho**GAMMA)

      ! Transform back into conservative variables
      vn   = 0.5_DP*(W(1)+W(2))
      c    = 0.25_DP*G1*(W(2)-W(1))
      rho  = (c*c/GAMMA/W(3))**G3
      p    = rho*c*c/GAMMA

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*DbdrNormal(1)*vn
      Du(3) = p/G1+0.5_DP*rho*(vn*vn)


    case(BDR_FARFIELD)
      !-------------------------------------------------------------------------

      ! The free stream primitive variables are Deval=[rho,v1,p]
      rho = DbdrValue(1)
      p   = DbdrValue(3)
      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*DbdrValue(2)

      ! Compute Riemann invariants based on the free stream values
      Winf(1) = vn-2.0_DP*c/G1
      Winf(2) = vn+2.0_DP*c/G1
      Winf(3) = p/(rho**GAMMA)

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      E   = Du(3)/rho
      p   = thdyn_pressure(GAMMA, E, rho, v1)

      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*v1

      ! Compute Riemann invariants based on the solution values
      Wu(1) = vn-2.0_DP*c/G1
      Wu(2) = vn+2.0_DP*c/G1
      Wu(3) = p/(rho**GAMMA)

      ! Adopt free stream/computed values depending on the sign of the eigenvalue
      W(1) = merge(Winf(1), Wu(1), vn <  c)
      W(2) = merge(Winf(2), Wu(2), vn < -c)
      W(3) = merge(Winf(3), Wu(3), vn <  SYS_EPSREAL)

      ! Transform back into conservative variables
      vn   = 0.5_DP*(W(1)+W(2))
      c    = 0.25_DP*G1*(W(2)-W(1))
      rho  = (c*c/GAMMA/W(3))**G3
      p    = rho*c*c/GAMMA

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*DbdrNormal(1)*vn
      Du(3) = p/G1+0.5_DP*rho*(vn*vn)


    case(BDR_SUBINLET)
      !-------------------------------------------------------------------------

      ! Compute primitive variables
      rho = Du(1)
      v1  = Du(2)/rho
      E   = Du(3)/rho
      p   = thdyn_pressure(GAMMA, E, rho, v1)

      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))
      vn  = DbdrNormal(1)*v1

      ! The specified density and pressure is Deval=[rho,p]
      rho = DbdrValue(1)
      p   = DbdrValue(2)

      ! Compute Riemann invariants
      W(1) = vn-2.0_DP*c/G1
      W(2) = vn+2.0_DP*c/G1
      W(3) = p/(rho**GAMMA)

      ! Transform back into conservative variables
      vn   = 0.5_DP*(W(1)+W(2))
      c    = 0.25_DP*G1*(W(2)-W(1))
      rho  = (c*c/GAMMA/W(3))**G3
      p    = rho*c*c/GAMMA

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*DbdrNormal(1)*vn
      Du(3) = p/G1+0.5_DP*rho*(vn*vn)


    case(BDR_SUBOUTLET)
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
      p   = thdyn_pressure(GAMMA, E, rho, v1)

      vn  = DbdrNormal(1)*v1
      c   = sqrt(max(GAMMA*p/rho, SYS_EPSREAL))

      ! Compute Riemann invariants based on the solution values and prescribed exit pressure
      W(2) = 2*c/G1-vn
      W(3) = p/(rho**GAMMA)
      W(1) = 4/G1*sqrt(max(GAMMA*ps/rho*(p/ps)**G4, SYS_EPSREAL))-W(2)

      ! Transform back into conservative variables
      vn  = 0.5_DP*(W(1)-W(2))
      c   = 0.25_DP*G1*(W(1)+W(2))
      rho = (c*c/GAMMA/W(3))**G3
      p   = rho*c*c/GAMMA

      ! Update the solution vector
      Du(1) = rho
      Du(2) = rho*DbdrNormal(1)*vn
      Du(3) = p/G1+0.5_DP*rho*(vn*vn)


    case DEFAULT
      call output_line('Unsupported type of boundary condition!',&
          OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcBoundaryvalues1d')
      call sys_halt()
    end select

  end subroutine eulerlagrange_calcBoundaryvalues1d

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_hadaptCallbackScalar1d(iOperation, rcollection)

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
            OU_CLASS_WARNING,OU_MODE_STD,'eulerlagrange_hadaptCallbackScalar1d')
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
            0.5_DP*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR1D+ivar)+&
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
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR1D+ivar) = 0.0_DP
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)

    end select

  end subroutine eulerlagrange_hadaptCallbackScalar1d

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_hadaptCallbackBlock1d(iOperation, rcollection)

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
            OU_CLASS_WARNING,OU_MODE_STD,'eulerlagrange_hadaptCallbackBlock1d')
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
            0.5_DP*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
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
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = 0.0_DP
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)

    end select

  end subroutine eulerlagrange_hadaptCallbackBlock1d

end module eulerlagrange_callback1d
