!##############################################################################
!# ****************************************************************************
!# <name> euler_callback1d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible Euler/Navier-Stokes equations in 1D.
!#
!# The following callback functions are available:
!#
!# 1.) euler_calcFluxGalerkin1d
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#
!# 2.) euler_calcFluxGalerkinNoBdr1d
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#        without assembling the symmetric boundary contribution
!#
!# 3.) euler_calcFluxScalarDiss1d
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting scalar artificial viscosities
!#
!# 4.) euler_calcFluxTensorDiss1d
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities
!#
!# 5.) euler_calcFluxRusanov1d
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion
!#
!# 6.) euler_calcMatrixDiagonalDiag1d
!#     -> Computes local matrix for diagonal entry
!#
!# 7.) euler_calcMatrixDiagonal1d
!#     -> Computes local matrix for diagonal entry
!#
!# 8.) euler_calcMatrixGalerkinDiag1d
!#     -> Computes local matrices for standard Galerkin scheme
!#
!# 9.) euler_calcMatrixGalerkin1d
!#     -> Computes local matrices for standard Galerkin scheme
!#
!# 10.) euler_calcMatrixScalarDissDiag1d
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 11.) euler_calcMatrixScalarDiss1d
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 12.) euler_calcMatrixTensorDissDiag1d
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 13.) euler_calcMatrixTensorDiss1d
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 14.) euler_calcMatrixRusanovDiag1d
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov artificial viscosities
!#
!# 15.) euler_calcMatrixRusanov1d
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov flux artificial viscosities
!#
!# 16.) euler_calcCharacteristics1d
!#      -> Computes characteristic variables
!#
!# 17.) euler_calcFluxFCTScalarDiss1d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting scalar artificial viscosities
!#
!# 18.) euler_calcFluxFCTTensorDiss1d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting tensorial artificial viscosities
!#
!# 19.) euler_calcFluxFCTRusanov1d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting the Rusanov artificial viscosities
!#
!# 20.) euler_trafoFluxDensity1d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 21.) euler_trafoDiffDensity1d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 22.) euler_trafoFluxEnergy1d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 23.) euler_trafoDiffEnergy1d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 24.) euler_trafoFluxPressure1d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 25.) euler_trafoDiffPressure1d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the pressure
!#
!# 26.) euler_trafoFluxVelocity1d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 27.) euler_trafoDiffVelocity1d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 28.) euler_trafoFluxMomentum1d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 29.) euler_trafoDiffMomentum1d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 30.) euler_trafoFluxDenEng1d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 31.) euler_trafoDiffDenEng1d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 32.) euler_trafoFluxDenPre1d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 33.) euler_trafoDiffDenPre1d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 34.) euler_trafoFluxDenPreVel1d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 35.) euler_trafoDiffDenPreVel1d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure 
!#         and the velocity
!#
!# 36.) euler_calcBoundaryvalues1d
!#      -> Computes the boundary values for a given node
!#
!# 37.) euler_hadaptCallbackScalar1d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 1D, whereby the vector is stored in interleave format
!#
!# 38.) euler_hadaptCallbackBlock1d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 1D, whereby the vector is stored in block format
!#
!# </purpose>
!##############################################################################

module euler_callback1d

  use boundaryfilter
  use collection
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
  public :: euler_calcFluxGalerkin1d
  public :: euler_calcFluxGalerkinNoBdr1d
  public :: euler_calcFluxScalarDiss1d
  public :: euler_calcFluxTensorDiss1d
  public :: euler_calcFluxRusanov1d
  public :: euler_calcMatrixDiagonalDiag1d
  public :: euler_calcMatrixDiagonal1d
  public :: euler_calcMatrixGalerkinDiag1d
  public :: euler_calcMatrixGalerkin1d
  public :: euler_calcMatrixScalarDissDiag1d
  public :: euler_calcMatrixScalarDiss1d
  public :: euler_calcMatrixTensorDissDiag1d
  public :: euler_calcMatrixTensorDiss1d
  public :: euler_calcMatrixRusanovDiag1d
  public :: euler_calcMatrixRusanov1d
  public :: euler_calcCharacteristics1d
  public :: euler_calcFluxFCTScalarDiss1d
  public :: euler_calcFluxFCTTensorDiss1d
  public :: euler_calcFluxFCTRusanov1d
  public :: euler_trafoFluxDensity1d
  public :: euler_trafoFluxEnergy1d
  public :: euler_trafoFluxPressure1d
  public :: euler_trafoFluxVelocity1d
  public :: euler_trafoFluxMomentum1d
  public :: euler_trafoFluxDenEng1d
  public :: euler_trafoFluxDenPre1d
  public :: euler_trafoFluxDenPreVel1d
  public :: euler_trafoDiffDensity1d  
  public :: euler_trafoDiffEnergy1d
  public :: euler_trafoDiffPressure1d
  public :: euler_trafoDiffVelocity1d
  public :: euler_trafoDiffMomentum1d
  public :: euler_trafoDiffDenEng1d
  public :: euler_trafoDiffDenPre1d
  public :: euler_trafoDiffDenPreVel1d
  public :: euler_calcBoundaryvalues1d
  public :: euler_hadaptCallbackScalar1d
  public :: euler_hadaptCallbackBlock1d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxGalerkin1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the standard
    ! Galerkin discretisation in 1D.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(out) :: F_ij, F_ji
!</output>
!</subroutine>

    ! local variables
#ifdef USE_EULER_IBP
    real(DP), dimension(NVAR1D) :: dF_i, dF_j
#else
    real(DP), dimension(NVAR1D) :: dF_ij
#endif
    real(DP) :: ui,uj,ru2i,ru2j

    !---------------------------------------------------------------------------
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
    ! ---------------------------------------------------------------------------

    ! Compute velocities
    ui = U_i(2)/U_i(1); uj = U_j(2)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); ru2j = uj*U_j(2)

#ifdef USE_EULER_IBP
    ! Compute fluxes for x-direction
    dF_i(1) = U_i(2)
    dF_i(2) = G1*U_i(3)-G14*ru2i
    dF_i(3) = (GAMMA*U_i(3)-G2*ru2i)*ui

    dF_j(1) = U_j(2)
    dF_j(2) = G1*U_j(3)-G14*ru2j
    dF_j(3) = (GAMMA*U_j(3)-G2*ru2j)*uj

    ! Assembly fluxes
    F_ij = dscale * ( C_ji(1)*dF_j - C_ij(1)*dF_i )
    F_ji = -F_ij
#else
    ! Compute flux difference for x-direction
    dF_ij(1) = U_i(2)                    - U_j(2)
    dF_ij(2) = G1*U_i(3)-G14*ru2i        - (G1*U_j(3)-G14*ru2j)
    dF_ij(3) = (GAMMA*U_i(3)-G2*ru2i)*ui - (GAMMA*U_j(3)-G2*ru2j)*uj

    ! Assembly fluxes
    F_ij =   dscale * C_ij(1)*dF_ij
    F_ji = - dscale * C_ji(1)*dF_ij
#endif

  end subroutine euler_calcFluxGalerkin1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxGalerkinNoBdr1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the TVD
    ! discretisation in 1D. The symmetric boundary contributions
    ! are neglected and incorporated in the antidiffusive flux.
    ! Hence, this is simply the standard Galerkin flux for the
    ! skew-symmetric internal contributions.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(out) :: F_ij, F_ji
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(NVAR1D) :: dF_ij
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: ui,uj,ru2i,ru2j

    !---------------------------------------------------------------------------
    ! Evaluate the Galerkin fluxes
    ! For a detailed description of algorithm and the definition of auxiliary
    ! quantities have a look at the subroutine "euler_calcFluxGalerkin1d".
    !---------------------------------------------------------------------------

    ! Compute velocities
    ui = U_i(2)/U_i(1); uj = U_j(2)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); ru2j = uj*U_j(2)

    ! Compute flux difference for x-direction
    dF_ij(1) = U_i(2)                    - U_j(2)
    dF_ij(2) = G1*U_i(3)-G14*ru2i        - (G1*U_j(3)-G14*ru2j)
    dF_ij(3) = (GAMMA*U_i(3)-G2*ru2i)*ui - (GAMMA*U_j(3)-G2*ru2j)*uj

    ! Compute skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji)

    ! Assembly fluxes
    F_ij = dscale * a(1)*dF_ij
    F_ji = F_ij

  end subroutine euler_calcFluxGalerkinNoBdr1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxScalarDiss1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 1D using scalar dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(out) :: F_ij, F_ji
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

    !---------------------------------------------------------------------------
    ! Evaluate the Galerkin fluxes
    ! For a detailed description of algorithm and the definition of auxiliary
    ! quantities have a look at the subroutine "euler_calcFluxGalerkin1d".
    !---------------------------------------------------------------------------

    ! Compute velocities
    ui = U_i(2)/U_i(1); uj = U_j(2)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); ru2j = uj*U_j(2)

#ifdef USE_EULER_IBP
    ! Compute fluxes for x-direction
    dF_i(1) = U_i(2)
    dF_i(2) = G1*U_i(3)-G14*ru2i
    dF_i(3) = (GAMMA*U_i(3)-G2*ru2i)*ui

    dF_j(1) = U_j(2)
    dF_j(2) = G1*U_j(3)-G14*ru2j
    dF_j(3) = (GAMMA*U_j(3)-G2*ru2j)*uj

    ! Assembly fluxes
    F_ij = dscale * ( C_ji(1)*dF_j - C_ij(1)*dF_i )
    F_ji = -F_ij
#else
    ! Compute flux difference for x-direction
    dF_ij(1) = U_i(2)                    - U_j(2)
    dF_ij(2) = G1*U_i(3)-G14*ru2i        - (G1*U_j(3)-G14*ru2j)
    dF_ij(3) = (GAMMA*U_i(3)-G2*ru2i)*ui - (GAMMA*U_j(3)-G2*ru2j)*uj

    ! Assembly fluxes
    F_ij =   dscale * C_ij(1)*dF_ij
    F_ji = - dscale * C_ji(1)*dF_ij
#endif

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji)

    ! Compute Roe mean values
    aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
    u_ij = (aux*ui+uj)/(aux+1.0_DP)
    hi   = GAMMA*U_i(3)/U_i(1)-G2*(U_i(2)*U_i(2))/(U_i(1)*U_i(1))
    hj   = GAMMA*U_j(3)/U_j(1)-G2*(U_j(2)*U_j(2))/(U_j(1)*U_j(1))
    H_ij = (aux*hi+hj)/(aux+1.0_DP)

    ! Compute auxiliary variables
    aux  = abs(a(1)) ! = sqrt(a(1)*a(1))
    vel  = u_ij*a(1)
    q_ij = 0.5_DP*(u_ij*u_ij)
    cs   = sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL))

    ! Scalar dissipation
    d_ij = dscale * (abs(vel) + aux*cs)

    ! Multiply the solution difference by the artificial diffusion factor
    Diff = d_ij*(U_j-U_i)

    ! Add the artificial diffusion to the fluxes
    F_ij = F_ij+Diff
    F_ji = F_ji-Diff

  end subroutine euler_calcFluxScalarDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxTensorDiss1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 1D using tensorial dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(out) :: F_ij, F_ji
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

    !---------------------------------------------------------------------------
    ! Evaluate the Galerkin fluxes
    ! For a detailed description of algorithm and the definition of auxiliary
    ! quantities have a look at the subroutine "euler_calcFluxGalerkin1d".
    !---------------------------------------------------------------------------

    ! Compute velocities
    ui = U_i(2)/U_i(1); uj = U_j(2)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); ru2j = uj*U_j(2)

#ifdef USE_EULER_IBP
    ! Compute fluxes for x-direction
    dF_i(1) = U_i(2)
    dF_i(2) = G1*U_i(3)-G14*ru2i
    dF_i(3) = (GAMMA*U_i(3)-G2*ru2i)*ui

    dF_j(1) = U_j(2)
    dF_j(2) = G1*U_j(3)-G14*ru2j
    dF_j(3) = (GAMMA*U_j(3)-G2*ru2j)*uj

    ! Assembly fluxes
    F_ij = dscale * ( C_ji(1)*dF_j - C_ij(1)*dF_i )
    F_ji = -F_ij
#else
    ! Compute flux difference for x-direction
    dF_ij(1) = U_i(2)                    - U_j(2)
    dF_ij(2) = G1*U_i(3)-G14*ru2i        - (G1*U_j(3)-G14*ru2j)
    dF_ij(3) = (GAMMA*U_i(3)-G2*ru2i)*ui - (GAMMA*U_j(3)-G2*ru2j)*uj

    ! Assembly fluxes
    F_ij =   dscale * C_ij(1)*dF_ij
    F_ji = - dscale * C_ji(1)*dF_ij
#endif

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute the skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji); anorm = abs(a(1))

    if (anorm .gt. SYS_EPSREAL) then

      ! Compute Roe mean values
      aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
      u_ij = (aux*ui+uj)/(aux+1.0_DP)
      hi   = GAMMA*U_i(3)/U_i(1)-G2*(U_i(2)*U_i(2))/(U_i(1)*U_i(1))
      hj   = GAMMA*U_j(3)/U_j(1)-G2*(U_j(2)*U_j(2))/(U_j(1)*U_j(1))
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
      Diff = U_j-U_i

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

      ! Compute "anorm * dscale"
      anorm = anorm*dscale

      ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
      Diff(1) = anorm * ( w1 + w2 + w3 )
      Diff(2) = anorm * ( (u_ij-cs)*w1 + u_ij*w2 + (u_ij+cs)*w3 )
      Diff(3) = anorm * ( (H_ij-u_ij*cs)*w1 + q_ij*w2 + (H_ij+u_ij*cs)*w3 )

      ! Add the artificial diffusion to the fluxes
      F_ij = F_ij+Diff
      F_ji = F_ji-Diff

    end if

  end subroutine euler_calcFluxTensorDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxRusanov1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 1D using the Rusanov dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(out) :: F_ij, F_ji
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

    !---------------------------------------------------------------------------
    ! Evaluate the Galerkin fluxes
    ! For a detailed description of algorithm and the definition of auxiliary
    ! quantities have a look at the subroutine "euler_calcFluxGalerkin1d".
    !---------------------------------------------------------------------------

    ! Compute velocities and energy
    ui = U_i(2)/U_i(1); uj = U_j(2)/U_j(1)
    Ei = U_i(3)/U_i(1); Ej = U_j(3)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); ru2j = uj*U_j(2)

#ifdef USE_EULER_IBP
    ! Compute fluxes for x-direction
    dF_i(1) = U_i(2)
    dF_i(2) = G1*U_i(3)-G14*ru2i
    dF_i(3) = (GAMMA*U_i(3)-G2*ru2i)*ui

    dF_j(1) = U_j(2)
    dF_j(2) = G1*U_j(3)-G14*ru2j
    dF_j(3) = (GAMMA*U_j(3)-G2*ru2j)*uj

    ! Assembly fluxes
    F_ij = dscale * ( C_ji(1)*dF_j - C_ij(1)*dF_i )
    F_ji = -F_ij
#else
    ! Compute flux difference for x-direction
    dF_ij(1) = U_i(2)                    - U_j(2)
    dF_ij(2) = G1*U_i(3)-G14*ru2i        - (G1*U_j(3)-G14*ru2j)
    dF_ij(3) = (GAMMA*U_i(3)-G2*ru2i)*ui - (GAMMA*U_j(3)-G2*ru2j)*uj

    ! Assembly fluxes
    F_ij =   dscale * C_ij(1)*dF_ij
    F_ji = - dscale * C_ji(1)*dF_ij
#endif

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute the speed of sound
    ci = sqrt(max(G15*(Ei-0.5_DP*ui*ui), SYS_EPSREAL))
    cj = sqrt(max(G15*(Ej-0.5_DP*uj*uj), SYS_EPSREAL))

    ! Scalar dissipation
    d_ij = max( abs(C_ij(1)*uj) + abs(C_ij(1))*cj,&
                abs(C_ji(1)*ui) + abs(C_ji(1))*ci )

    ! Multiply the solution difference by the artificial diffusion factor
    Diff = dscale * d_ij*(U_j-U_i)

    ! Add the artificial diffusion to the fluxes
    F_ij = F_ij+Diff
    F_ji = F_ji-Diff

  end subroutine euler_calcFluxRusanov1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixDiagonalDiag1d(U_i, C_ii, i, dscale, K_ii)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 1D
!</description>

!<input>
    ! local solution at node I
    real(DP), dimension(:), intent(in) :: U_i

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ii

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! row number
    integer, intent(in) :: i
!</input>

!<output>
    ! local Jacobian matrix
    real(DP), dimension(:), intent(out) :: K_ii
!</output>
!</subroutine>

    ! local variable
    real(DP) :: ui

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1)

    ! Compute Galerkin coefficient K_ii
    K_ii(1) = 0.0_DP
    K_ii(2) = dscale * G13*ui*C_ii(1)
    K_ii(3) = dscale * GAMMA*ui*C_ii(1)

  end subroutine euler_calcMatrixDiagonalDiag1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixDiagonal1d(U_i, C_ii, i, dscale, K_ii)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 1D
!</description>

!<input>
    ! local solution at node I
    real(DP), dimension(:), intent(in) :: U_i

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ii

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! row number
    integer, intent(in) :: i
!</input>

!<output>
    ! local Jacobian matrix
    real(DP), dimension(:), intent(out) :: K_ii
!</output>
!</subroutine>

    ! local variable
    real(DP) :: ui,Ei,uPow2i

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   Ei = U_i(3)/U_i(1);   uPow2i = ui*ui

    ! Compute Galerkin coefficient K_ii
    K_ii(1) = 0.0_DP
    K_ii(2) = dscale * G14*uPow2i*C_ii(1)
    K_ii(3) = dscale * (G1*uPow2i-GAMMA*Ei)*ui*C_ii(1)

    K_ii(4) = dscale * C_ii(1)
    K_ii(5) = dscale * G13*ui*C_ii(1)
    K_ii(6) = dscale * (GAMMA*Ei-G16*uPow2i)*C_ii(1)

    K_ii(7) = 0.0_DP
    K_ii(8) = dscale * G1*C_ii(1)
    K_ii(9) = dscale * GAMMA*ui*C_ii(1)

  end subroutine euler_calcMatrixDiagonal1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixGalerkinDiag1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

    ! local variable
    real(DP) :: ui,uj

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   uj = U_j(2)/U_j(1)

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale * G13*uj*C_ij(1)
    K_ij(3) = dscale * GAMMA*uj*C_ij(1)

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale * G13*ui*C_ji(1)
    K_ji(3) = dscale * GAMMA*ui*C_ji(1)

    ! Nullify dissipation tensor
    D_ij = 0.0_DP

  end subroutine euler_calcMatrixGalerkinDiag1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixGalerkin1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

    ! local variable
    real(DP) :: ui,uj,Ei,Ej,uPow2i,uPow2j

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   Ei = U_i(3)/U_i(1);   uPow2i = ui*ui
    uj = U_j(2)/U_j(1);   Ej = U_j(3)/U_j(1);   uPow2j = uj*uj

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale * G14*uPow2j*C_ij(1)
    K_ij(3) = dscale * (G1*uPow2j-GAMMA*Ej)*uj*C_ij(1)

    K_ij(4) = dscale * C_ij(1)
    K_ij(5) = dscale * G13*uj*C_ij(1)
    K_ij(6) = dscale * (GAMMA*Ej-G16*uPow2j)*C_ij(1)

    K_ij(7) = 0.0_DP
    K_ij(8) = dscale * G1*C_ij(1)
    K_ij(9) = dscale * GAMMA*uj*C_ij(1)

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale * G14*uPow2i*C_ji(1)
    K_ji(3) = dscale * (G1*uPow2i-GAMMA*Ei)*ui*C_ji(1)

    K_ji(4) = dscale * C_ji(1)
    K_ji(5) = dscale * G13*ui*C_ji(1)
    K_ji(6) = dscale * (GAMMA*Ei-G16*uPow2i)*C_ji(1)

    K_ji(7) = 0.0_DP
    K_ji(8) = dscale * G1*C_ji(1)
    K_ji(9) = dscale * GAMMA*ui*C_ji(1)

    ! Nullify dissipation tensor
    D_ij = 0.0_DP

  end subroutine euler_calcMatrixGalerkin1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixScalarDissDiag1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies scalar artificial viscosities in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: anorm,aux,hi,hj,H_ij,q_ij,ui,uj,u_ij

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   uj = U_j(2)/U_j(1)

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale * G13*uj*C_ij(1)
    K_ij(3) = dscale * GAMMA*uj*C_ij(1)

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale * G13*ui*C_ji(1)
    K_ji(3) = dscale * GAMMA*ui*C_ji(1)

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute skew-symmetric coefficient and its norm
    a = 0.5_DP*(C_ji-C_ij); anorm = abs(a(1))

    if (anorm .gt. SYS_EPSREAL) then

      ! Compute Roe mean values
      aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
      u_ij = (aux*ui+uj)/(aux+1.0_DP)
      hi   = GAMMA*U_i(3)/U_i(1)-G2*(ui*ui)
      hj   = GAMMA*U_j(3)/U_j(1)-G2*(uj*uj)
      H_ij = (aux*hi+hj)/(aux+1.0_DP)

      ! Compute auxiliary values
      q_ij = 0.5_DP*u_ij*u_ij

      ! Compute scalar dissipation
      D_ij = dscale * (abs(a(1)*u_ij) +&
                       anorm*sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL)))
    else

      ! Nullify dissipation tensor
      D_ij = 0.0_DP

    end if

  end subroutine euler_calcMatrixScalarDissDiag1d

!*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixScalarDiss1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies scalar artificial viscosities in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: anorm,aux,hi,hj,Ei,Ej,H_ij,q_ij,ui,uj,u_ij,uPow2i,uPow2j

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   Ei = U_i(3)/U_i(1);   uPow2i = ui*ui
    uj = U_j(2)/U_j(1);   Ej = U_j(3)/U_j(1);   uPow2j = uj*uj

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale * G14*uPow2j*C_ij(1)
    K_ij(3) = dscale * (G1*uPow2j-GAMMA*Ej)*uj*C_ij(1)

    K_ij(4) = dscale * C_ij(1)
    K_ij(5) = dscale * G13*uj*C_ij(1)
    K_ij(6) = dscale * (GAMMA*Ej-G16*uPow2j)*C_ij(1)

    K_ij(7) = 0.0_DP
    K_ij(8) = dscale * G1*C_ij(1)
    K_ij(9) = dscale * GAMMA*uj*C_ij(1)

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale * G14*uPow2i*C_ji(1)
    K_ji(3) = dscale * (G1*uPow2i-GAMMA*Ei)*ui*C_ji(1)

    K_ji(4) = dscale * C_ji(1)
    K_ji(5) = dscale * G13*ui*C_ji(1)
    K_ji(6) = dscale * (GAMMA*Ei-G16*uPow2i)*C_ji(1)

    K_ji(7) = 0.0_DP
    K_ji(8) = dscale * G1*C_ji(1)
    K_ji(9) = dscale * GAMMA*ui*C_ji(1)

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute coefficients
    a = 0.5_DP*(C_ji-C_ij); anorm = abs(a(1))

    if (anorm .gt. SYS_EPSREAL) then

      ! Compute Roe mean values
      aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
      u_ij = (aux*ui+uj)/(aux+1.0_DP)
      hi   = GAMMA*U_i(3)/U_i(1)-G2*(ui*ui)
      hj   = GAMMA*U_j(3)/U_j(1)-G2*(uj*uj)
      H_ij = (aux*hi+hj)/(aux+1.0_DP)

      ! Compute auxiliary values
      q_ij = 0.5_DP*u_ij*u_ij

      ! Compute scalar dissipation
      aux = dscale * (abs(a(1)*u_ij) +&
                      anorm*sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL)))

      D_ij    = 0.0_DP
      D_ij(1) = aux
      D_ij(5) = aux
      D_ij(9) = aux

    else

      ! Nullify dissipation tensor
      D_ij = 0.0_DP

    end if

  end subroutine euler_calcMatrixScalarDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixTensorDissDiag1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NVAR1D,NVAR1D) :: R_ij,L_ij
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: aux,hi,hj,H_ij,q_ij,ui,uj,u_ij
    real(DP) :: l1,l2,l3,anorm,cs,cPow2,b1,b2

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   uj = U_j(2)/U_j(1)

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale * G13*uj*C_ij(1)
    K_ij(3) = dscale * GAMMA*uj*C_ij(1)

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale * G13*ui*C_ji(1)
    K_ji(3) = dscale * GAMMA*ui*C_ji(1)

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute skew-symmetric coefficient and its norm
    a = 0.5_DP*(C_ji-C_ij); anorm = abs(a(1))

    if (anorm .gt. SYS_EPSREAL) then

      ! Compute Roe mean values
      aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
      u_ij = (aux*ui+uj)/(aux+1.0_DP)
      hi   = GAMMA*U_i(3)/U_i(1)-G2*(ui*ui)
      hj   = GAMMA*U_j(3)/U_j(1)-G2*(uj*uj)
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

      ! Compute tensorial dissipation D_ij = diag(R_ij*|Lbd_ij|*L_ij)*I
      D_ij    = 0.0_DP
      D_ij(1) = anorm*( R_ij(1,1)*L_ij(1,1)+&
                        R_ij(1,2)*L_ij(2,1)+&
                        R_ij(1,3)*L_ij(3,1)  )
      D_ij(2) = anorm*( R_ij(2,1)*L_ij(1,2)+&
                        R_ij(2,2)*L_ij(2,2)+&
                        R_ij(2,3)*L_ij(3,2)  )
      D_ij(3) = anorm*( R_ij(3,1)*L_ij(1,3)+&
                        R_ij(3,2)*L_ij(2,3)+&
                        R_ij(3,3)*L_ij(3,3)  )
    else

      ! Nullify dissipation tensor
      D_ij = 0.0_DP

    end if

  end subroutine euler_calcMatrixTensorDissDiag1d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixTensorDiss1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NVAR1D,NVAR1D) :: R_ij,L_ij
    real(DP), dimension(NDIM1D) :: a
    real(DP) :: aux,Ei,Ej,hi,hj,H_ij,q_ij,ui,uj,u_ij
    real(DP) :: l1,l2,l3,anorm,cs,cPow2,b1,b2,uPow2i,uPow2j

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   Ei = U_i(3)/U_i(1);   uPow2i = ui*ui
    uj = U_j(2)/U_j(1);   Ej = U_j(3)/U_j(1);   uPow2j = uj*uj

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale * G14*uPow2j*C_ij(1)
    K_ij(3) = dscale * (G1*uPow2j-GAMMA*Ej)*uj*C_ij(1)

    K_ij(4) = dscale * C_ij(1)
    K_ij(5) = dscale * G13*uj*C_ij(1)
    K_ij(6) = dscale * (GAMMA*Ej-G16*uPow2j)*C_ij(1)

    K_ij(7) = 0.0_DP
    K_ij(8) = dscale * G1*C_ij(1)
    K_ij(9) = dscale * GAMMA*uj*C_ij(1)

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale * G14*uPow2i*C_ji(1)
    K_ji(3) = dscale * (G1*uPow2i-GAMMA*Ei)*ui*C_ji(1)

    K_ji(4) = dscale * C_ji(1)
    K_ji(5) = dscale * G13*ui*C_ji(1)
    K_ji(6) = dscale * (GAMMA*Ei-G16*uPow2i)*C_ji(1)

    K_ji(7) = 0.0_DP
    K_ji(8) = dscale * G1*C_ji(1)
    K_ji(9) = dscale * GAMMA*ui*C_ji(1)

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute coefficients
    a = 0.5_DP*(C_ji-C_ij); anorm = abs(a(1))

    if (anorm .gt. SYS_EPSREAL) then

      ! Compute Roe mean values
      aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
      u_ij = (aux*ui+uj)/(aux+1.0_DP)
      hi   = GAMMA*U_i(3)/U_i(1)-G2*(ui*ui)
      hj   = GAMMA*U_j(3)/U_j(1)-G2*(uj*uj)
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

      ! Compute tensorial dissipation D_ij = R_ij*|Lbd_ij|*L_ij
      call DGEMM('n', 'n', NVAR1D, NVAR1D, NVAR1D, anorm,&
          R_ij, NVAR1D, L_ij, NVAR1D, 0.0_DP, D_ij, NVAR1D)

    else

      ! Nullify dissipation tensor
      D_ij = 0.0_DP

    end if

  end subroutine euler_calcMatrixTensorDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixRusanovDiag1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

    ! local variable
    real(DP) :: ui,uj,ci,cj,Ei,Ej

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   uj = U_j(2)/U_j(1)
    Ei = U_i(3)/U_i(1);   Ej = U_j(3)/U_j(1)

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale * G13*uj*C_ij(1)
    K_ij(3) = dscale * GAMMA*uj*C_ij(1)

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale * G13*ui*C_ji(1)
    K_ji(3) = dscale * GAMMA*ui*C_ji(1)

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute the speed of sound
    ci = sqrt(max(G15*(Ei-0.5_DP*ui*ui), SYS_EPSREAL))
    cj = sqrt(max(G15*(Ej-0.5_DP*uj*uj), SYS_EPSREAL))

    ! Compute dissipation tensor D_ij
    D_ij = dscale * max( abs(C_ij(1)*uj) + abs(C_ij(1))*cj,&
                         abs(C_ji(1)*ui) + abs(C_ji(1))*ci )

  end subroutine euler_calcMatrixRusanovDiag1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixRusanov1d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! local matrices
    real(DP), dimension(:), intent(out) :: K_ij,K_ji,D_ij
!</output>
!</subroutine>

    ! local variable
    real(DP) :: ui,uj,ci,cj,Ei,Ej,uPow2i,uPow2j,aux

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   Ei = U_i(3)/U_i(1);   uPow2i = ui*ui
    uj = U_j(2)/U_j(1);   Ej = U_j(3)/U_j(1);   uPow2j = uj*uj

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale * G14*uPow2j*C_ij(1)
    K_ij(3) = dscale * (G1*uPow2j-GAMMA*Ej)*uj*C_ij(1)

    K_ij(4) = dscale * C_ij(1)
    K_ij(5) = dscale * G13*uj*C_ij(1)
    K_ij(6) = dscale * (GAMMA*Ej-G16*uPow2j)*C_ij(1)

    K_ij(7) = 0.0_DP
    K_ij(8) = dscale * G1*C_ij(1)
    K_ij(9) = dscale * GAMMA*uj*C_ij(1)

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale * G14*uPow2i*C_ji(1)
    K_ji(3) = dscale * (G1*uPow2i-GAMMA*Ei)*ui*C_ji(1)

    K_ji(4) = dscale * C_ji(1)
    K_ji(5) = dscale * G13*ui*C_ji(1)
    K_ji(6) = dscale * (GAMMA*Ei-G16*uPow2i)*C_ji(1)

    K_ji(7) = 0.0_DP
    K_ji(8) = dscale * G1*C_ji(1)
    K_ji(9) = dscale * GAMMA*ui*C_ji(1)

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute the speed of sound
    ci = sqrt(max(G15*(Ei-0.5_DP*ui*ui), SYS_EPSREAL))
    cj = sqrt(max(G15*(Ej-0.5_DP*uj*uj), SYS_EPSREAL))

    ! Compute dissipation tensor D_ij
    aux = dscale * max( abs(C_ij(1)*uj) + abs(C_ij(1))*cj,&
                        abs(C_ji(1)*ui) + abs(C_ji(1))*ci )

    D_ij = 0.0_DP
    D_ij(1) = aux
    D_ij(5) = aux
    D_ij(9) = aux

  end subroutine euler_calcMatrixRusanov1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcCharacteristics1d(&
      U_i, U_j, Dweight, W_ij, Lbd_ij, R_ij, L_ij)

!<description>
    ! This subroutine computes the characteristic variables in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! weighting vector
    real(DP), dimension(:), intent(in) :: Dweight
!</input>

!<output>
    ! vector of characteristic variables
    real(DP), dimension(:), intent(out), optional :: W_ij

    ! OPTIONAL: diagonal matrix of eigenvalues
    real(DP), dimension(:), intent(out), optional :: Lbd_ij

    ! OPTIONAL: transformation matrix into conservative variables
    real(DP), dimension(:), intent(out), optional :: R_ij

    ! OPTIONAL: transformation matrix into characteristic variables
    real(DP), dimension(:), intent(out), optional :: L_ij
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(NVAR1D) :: Diff
    real(DP) :: u_ij,H_ij,q_ij,cs,aux,aux1,aux2,hi,hj
    real(DP) :: cPow2,uPow2,a1,anorm,b1,b2

    ! Compute norm of weighting coefficient
    anorm = abs(Dweight(1))

    ! Check if weighting coefficient is zero
    if (anorm .gt. SYS_EPSREAL) then

      ! Compute normalised weighting coefficient
      a1  = Dweight(1)/anorm

      ! Compute Roe mean values
      aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
      u_ij = (aux*U_i(2)/U_i(1)+U_j(2)/U_j(1))/(aux+1.0_DP)
      hi   = GAMMA*U_i(3)/U_i(1)-G2*(U_i(2)*U_i(2))/(U_i(1)*U_i(1))
      hj   = GAMMA*U_j(3)/U_j(1)-G2*(U_j(2)*U_j(2))/(U_j(1)*U_j(1))
      H_ij = (aux*hi+hj)/(aux+1.0_DP)

      ! Compute auxiliary variables
      uPow2 = u_ij*u_ij
      q_ij  = 0.5_DP*(uPow2)
      cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
      cs    = sqrt(cPow2)
      b2    = G1/cPow2
      b1    = b2*q_ij

      ! Compute diagonal matrix of eigenvalues (if present)
      if (present(Lbd_ij)) then
        Lbd_ij(1) = u_ij-cs
        Lbd_ij(2) = u_ij
        Lbd_ij(3) = u_ij+cs
      end if

      ! Compute matrix of right eigenvectors
      if (present(R_ij)) then
        R_ij(1) =  1.0_DP
        R_ij(2) =  u_ij-cs
        R_ij(3) =  H_ij-u_ij*cs

        R_ij(4) =  1.0_DP
        R_ij(5) =  u_ij
        R_ij(6) =  q_ij

        R_ij(7) =  1.0_DP
        R_ij(8) =  u_ij+cs
        R_ij(9) =  H_ij+u_ij*cs
      end if

      ! Compute matrix of left eigenvectors
      if (present(L_ij)) then
        L_ij(1) = 0.5_DP * (b1+u_ij/cs)
        L_ij(2) =          (1-b1)
        L_ij(3) = 0.5_DP * (b1-u_ij/cs)

        L_ij(4) =-0.5_DP * (b2*u_ij+1.0_DP/cs)
        L_ij(5) =          (b2*u_ij)
        L_ij(6) =-0.5_DP * (b2*u_ij-1.0_DP/cs)

        L_ij(7) = 0.5_DP*b2
        L_ij(8) =       -b2
        L_ij(9) = 0.5_DP*b2
      end if

      ! Compute characteristic solution difference
      if (present(W_ij)) then
        ! Compute solution difference U_j-U_i
        Diff = U_j-U_i

        ! Compute characteristic variables
        W_ij(1) = anorm * 0.5_DP * (       (b1+u_ij/cs)*Diff(1) -&
                                    (b2*u_ij+1.0_DP/cs)*Diff(2) +&
                                                     b2*Diff(3) )
        W_ij(2) = anorm * (                      (1-b1)*Diff(1) +&
                                                b2*u_ij*Diff(2) -&
                                                     b2*Diff(3) )
        W_ij(3) = anorm * 0.5_DP * (       (b1-u_ij/cs)*Diff(1) -&
                                    (b2*u_ij-1.0_DP/cs)*Diff(2) +&
                                                     b2*Diff(3) )
      end if

    else   ! |dweight| = 0

      if (present(Lbd_ij)) Lbd_ij = 0.0_DP
      if (present(R_ij))   R_ij   = 0.0_DP
      if (present(L_ij))   L_ij   = 0.0_DP
      if (present(W_ij))   W_ij   = 0.0_DP

    end if

  end subroutine euler_calcCharacteristics1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxFCTScalarDiss1d(&
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

  end subroutine euler_calcFluxFCTScalarDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxFCTTensorDiss1d(&
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

  end subroutine euler_calcFluxFCTTensorDiss1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxFCTRusanov1d(&
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

  end subroutine euler_calcFluxFCTRusanov1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxDensity1d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! flux
    real(DP), dimension(:), intent(in) :: F_ij
!</input>

!<output>
    ! transformed flux
    real(DP), dimension(:), intent(out) :: G_ij,G_ji
!</output>
!</subroutine>

    ! density fluxes
    G_ij(1) =  F_ij(1)
    G_ji(1) = -F_ij(1)

  end subroutine euler_trafoFluxDensity1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffDensity1d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j
!</input>

!<output>
    ! transformed difference
    real(DP), dimension(:), intent(out) :: U_ij
!</output>
!</subroutine>

    ! density difference
    U_ij(1) =  U_j(1)-U_i(1)

  end subroutine euler_trafoDiffDensity1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxEnergy1d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the energy in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! flux
    real(DP), dimension(:), intent(in) :: F_ij
!</input>

!<output>
    ! transformed flux
    real(DP), dimension(:), intent(out) :: G_ij,G_ji
!</output>
!</subroutine>
    
    ! energy fluxes
    G_ij(1) =  F_ij(3)
    G_ji(1) = -F_ij(3)

  end subroutine euler_trafoFluxEnergy1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffEnergy1d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the energy in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j
!</input>

!<output>
    ! transformed difference
    real(DP), dimension(:), intent(out) :: U_ij
!</output>
!</subroutine>
    
    ! energy difference
    U_ij(1) =  U_j(3)-U_i(3)

  end subroutine euler_trafoDiffEnergy1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxPressure1d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the pressure in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! flux
    real(DP), dimension(:), intent(in) :: F_ij
!</input>

!<output>
    ! transformed flux
    real(DP), dimension(:), intent(out) :: G_ij,G_ji
!</output>
!</subroutine>

    ! local variables
    real(DP) :: ui,uj

    ! velocities
    ui = U_i(2)/U_i(1)
    uj = U_j(2)/U_j(1)

    ! pressure fluxes
    G_ij(1) =  G1*(0.5_DP*ui*ui*F_ij(1)-ui*F_ij(2)+F_ij(3))
    G_ji(1) = -G1*(0.5_DP*uj*uj*F_ij(1)-uj*F_ij(2)+F_ij(3))
    
  end subroutine euler_trafoFluxPressure1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffPressure1d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the pressure in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j
!</input>

!<output>
    ! transformed difference
    real(DP), dimension(:), intent(out) :: U_ij
!</output>
!</subroutine>

    ! local variables
    real(DP) :: pi,pj

    ! pressures
    pi = G1*(U_i(3)-0.5_DP*U_i(2)*U_i(2)/U_i(1))
    pj = G1*(U_j(3)-0.5_DP*U_j(2)*U_j(2)/U_j(1))

    ! pressure difference
    U_ij(1) = pj-pi

  end subroutine euler_trafoDiffPressure1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxVelocity1d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the x-velocity
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! flux
    real(DP), dimension(:), intent(in) :: F_ij
!</input>

!<output>
    ! transformed flux
    real(DP), dimension(:), intent(out) :: G_ij,G_ji
!</output>
!</subroutine>

    ! local variables
    real(DP) :: ui,uj

    ! velocities
    ui = U_i(2)/U_i(1)
    uj = U_j(2)/U_j(1)

    ! velocity fluxes in x-direction
    G_ij(1) =  (F_ij(2)-ui*F_ij(1))/U_i(1)
    G_ji(1) = -(F_ij(2)-uj*F_ij(1))/U_j(1)
    
  end subroutine euler_trafoFluxVelocity1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffVelocity1d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the x-velocity
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j
!</input>

!<output>
    ! transformed differences
    real(DP), dimension(:), intent(out) :: U_ij
!</output>
!</subroutine>

    ! velocity difference in x-direction
    U_ij(1) =  U_j(2)/U_j(1)-U_i(2)/U_i(1)

  end subroutine euler_trafoDiffVelocity1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxMomentum1d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the x-momentum
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! flux
    real(DP), dimension(:), intent(in) :: F_ij
!</input>

!<output>
    ! transformed flux
    real(DP), dimension(:), intent(out) :: G_ij,G_ji
!</output>
!</subroutine>

    ! momentum fluxes in x-direction
    G_ij(1) =  F_ij(2)
    G_ji(1) = -F_ij(2)
    
  end subroutine euler_trafoFluxMomentum1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffMomentum1d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the x-momentum
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j
!</input>

!<output>
    ! transformed differences
    real(DP), dimension(:), intent(out) :: U_ij
!</output>
!</subroutine>

    ! momentum difference in x-direction
    U_ij(1) =  U_j(2)-U_i(2)
    
  end subroutine euler_trafoDiffMomentum1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxDenEng1d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! flux
    real(DP), dimension(:), intent(in) :: F_ij
!</input>

!<output>
    ! transformed flux
    real(DP), dimension(:), intent(out) :: G_ij,G_ji
!</output>
!</subroutine>

    ! density fluxes
    G_ij(1) =  F_ij(1)
    G_ji(1) = -F_ij(1)

    ! energy fluxes
    G_ij(2) =  F_ij(3)
    G_ji(2) = -F_ij(3)

  end subroutine euler_trafoFluxDenEng1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffDenEng1d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j
!</input>

!<output>
    ! transformed differences
    real(DP), dimension(:), intent(out) :: U_ij
!</output>
!</subroutine>

    ! density difference
    U_ij(1) =  U_j(1)-U_i(1)

    ! energy difference
    U_ij(2) =  U_j(3)-U_i(3)

  end subroutine euler_trafoDiffDenEng1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxDenPre1d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! flux
    real(DP), dimension(:), intent(in) :: F_ij
!</input>

!<output>
    ! transformed flux
    real(DP), dimension(:), intent(out) :: G_ij,G_ji
!</output>
!</subroutine>

    ! local variables
    real(DP) :: ui,uj

    ! velocities
    ui = U_i(2)/U_i(1)
    uj = U_j(2)/U_j(1)

    ! density fluxes
    G_ij(1) =  F_ij(1)
    G_ji(1) = -F_ij(1)

    ! pressure fluxes
    G_ij(2) =  G1*(0.5_DP*ui*ui*F_ij(1)-ui*F_ij(2)+F_ij(3))
    G_ji(2) = -G1*(0.5_DP*uj*uj*F_ij(1)-uj*F_ij(2)+F_ij(3))

  end subroutine euler_trafoFluxDenPre1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffDenPre1d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j
!</input>

!<output>
    ! transformed differences
    real(DP), dimension(:), intent(out) :: U_ij
!</output>
!</subroutine>

    ! local variables
    real(DP) :: pi,pj

    ! pressures
    pi = G1*(U_i(3)-0.5_DP*U_i(2)*U_i(2)/U_i(1))
    pj = G1*(U_j(3)-0.5_DP*U_j(2)*U_j(2)/U_j(1))

    ! density difference
    U_ij(1) = U_j(1)-U_i(1)

    ! pressure difference
    U_ij(2) = pj-pi

  end subroutine euler_trafoDiffDenPre1d
  
  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxDenPreVel1d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density, pressure and velocity in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! flux
    real(DP), dimension(:), intent(in) :: F_ij
!</input>

!<output>
    ! transformed flux
    real(DP), dimension(:), intent(out) :: G_ij,G_ji
!</output>
!</subroutine>

    ! local variables
    real(DP) :: ui,uj

    ! velocities
    ui = U_i(2)/U_i(1)
    uj = U_j(2)/U_j(1)

    ! density fluxes
    G_ij(1) =  F_ij(1)
    G_ji(1) = -F_ij(1)

    ! velocity fluxes in x-direction
    G_ij(2) =  (F_ij(2)-ui*F_ij(1))/U_i(1)
    G_ji(2) = -(F_ij(2)-uj*F_ij(1))/U_j(1)

    ! pressure fluxes
    G_ij(3) =  G1*(0.5_DP*ui*ui*F_ij(1)-ui*F_ij(2)+F_ij(3))
    G_ji(3) = -G1*(0.5_DP*uj*uj*F_ij(1)-uj*F_ij(2)+F_ij(3))

  end subroutine euler_trafoFluxDenPreVel1d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffDenPreVel1d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density, pressure and velocity in 1D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j
!</input>

!<output>
    ! transformed differences
    real(DP), dimension(:), intent(out) :: U_ij
!</output>
!</subroutine>

    ! local variables
    real(DP) :: pi,pj

    ! pressures
    pi = G1*(U_i(3)-0.5_DP*U_i(2)*U_i(2)/U_i(1))
    pj = G1*(U_j(3)-0.5_DP*U_j(2)*U_j(2)/U_j(1))

    ! density difference
    U_ij(1) = U_j(1)-U_i(1)

    ! velocity difference in x-direction
    U_ij(2) =  U_j(2)/U_j(1)-U_i(2)/U_i(1)
    
    ! pressure difference
    U_ij(3) = pj-pi

  end subroutine euler_trafoDiffDenPreVel1d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcBoundaryvalues1d(DbdrNormal, DpointNormal,&
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
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcBoundaryvalues1d')
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
          OU_CLASS_ERROR,OU_MODE_STD,'euler_calcBoundaryvalues1d')
      call sys_halt()
    end select

  end subroutine euler_calcBoundaryvalues1d

  !*****************************************************************************

!<subroutine>

  subroutine euler_hadaptCallbackScalar1d(iOperation, rcollection)

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
            OU_CLASS_WARNING,OU_MODE_STD,'euler_hadaptCallbackScalar1d')
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

  end subroutine euler_hadaptCallbackScalar1d

  !*****************************************************************************

!<subroutine>

  subroutine euler_hadaptCallbackBlock1d(iOperation, rcollection)

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
            OU_CLASS_WARNING,OU_MODE_STD,'euler_hadaptCallbackBlock1d')
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

  end subroutine euler_hadaptCallbackBlock1d

end module euler_callback1d
