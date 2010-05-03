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
!# 1.) euler_calcFluxGalerkin2d
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#
!# 2.) euler_calcFluxGalerkinNoBdr2d
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#        without assembling the symmetric boundary contribution
!#
!# 3.) euler_calcFluxScalarDiss2d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting scalar artificial viscosities
!#
!# 4.) euler_calcFluxScalarDissDiSp2d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting scalar artificial viscosities based on
!#        dimensional splitting approach
!#
!# 5.) euler_calcFluxTensorDiss2d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting tensorial artificial viscosities
!#
!# 6.) euler_calcFluxTensorDissDiSp2d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting tensorial artificial viscosities based on
!#        dimensional splitting approach
!#
!# 7.) euler_calcFluxRusanov2d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting the Rusanov artificial diffusion based on
!#        dimensional splitting approach
!#
!# 8.) euler_calcFluxRusanovDiSp2d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting the Rusanov artificial diffusion
!#
!# 9.) euler_calcMatrixDiagonalDiag2d
!#     -> Computes local matrix for diagonal entry
!#
!# 10.) euler_calcMatrixDiagonal2d
!#      -> Computes local matrix for diagonal entry
!#
!# 11.) euler_calcMatrixGalerkinDiag2d
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 12.) euler_calcMatrixGalerkin2d
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 13.) euler_calcMatrixScalarDissDiag2d
!#      -> Computes local matrices for low-order discretization
!#         adopting scalar artificial viscosities
!#
!# 14.) euler_calcMatrixScalarDiss2d
!#      -> Computes local matrices for low-order discretization
!#         adopting scalar artificial viscosities
!#
!# 15.) euler_calcMatrixTensorDissDiag2d
!#      -> Computes local matrices for low-order discretization
!#         adopting tensorial artificial viscosities
!#
!# 16.) euler_calcMatrixTensorDiss2d
!#      -> Computes local matrices for low-order discretization
!#         adopting tensorial artificial viscosities
!#
!# 17.) euler_calcMatrixRusanovDiag2d
!#      -> Computes local matrices for low-order discretization
!#         adopting the Rusanov artificial viscosities
!#
!# 18.) euler_calcMatrixRusanov2d
!#      -> Computes local matrices for low-order discretization
!#         adopting the Rusanov flux artificial viscosities
!#
!# 19.) euler_calcCharacteristics2d
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
!# 23.) euler_trafoFluxDensity2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 24.) euler_trafoDiffDensity2d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 25.) euler_trafoFluxEnergy2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 26.) euler_trafoDiffEnergy2d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 27.) euler_trafoFluxPressure2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 28.) euler_trafoFluxVelocity2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 29.) euler_trafoDiffVelocity2d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 30.) euler_trafoFluxMomentum2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 31.) euler_trafoDiffMomentum2d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 32.) euler_trafoFluxDenEng2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 33.) euler_trafoDiffDenEng2d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 34.) euler_trafoFluxDenPre2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 35.) euler_trafoDiffDenPre2d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 36.) euler_trafoFluxDenPreVel2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 37.) euler_trafoDiffDenPreVel2d
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
!# 41.) euler_limitEdgewiseVelocity2D
!#      -> Performs synchronised flux correction for the velocity in 2D
!#
!# </purpose>
!##############################################################################

module euler_callback2d

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
  public :: euler_calcFluxGalerkin2d
  public :: euler_calcFluxGalerkinNoBdr2d
  public :: euler_calcFluxScalarDiss2d
  public :: euler_calcFluxScalarDissDiSp2d
  public :: euler_calcFluxTensorDiss2d
  public :: euler_calcFluxTensorDissDiSp2d
  public :: euler_calcFluxRusanov2d
  public :: euler_calcFluxRusanovDiSp2d
  public :: euler_calcMatrixDiagonalDiag2d
  public :: euler_calcMatrixDiagonal2d
  public :: euler_calcMatrixGalerkinDiag2d
  public :: euler_calcMatrixGalerkin2d
  public :: euler_calcMatrixScalarDissDiag2d
  public :: euler_calcMatrixScalarDiss2d
  public :: euler_calcMatrixTensorDissDiag2d
  public :: euler_calcMatrixTensorDiss2d
  public :: euler_calcMatrixRusanovDiag2d
  public :: euler_calcMatrixRusanov2d
  public :: euler_calcCharacteristics2d
  public :: euler_calcFluxFCTScalarDiss2d
  public :: euler_calcFluxFCTTensorDiss2d
  public :: euler_calcFluxFCTRusanov2d
  public :: euler_trafoFluxDensity2d
  public :: euler_trafoFluxEnergy2d
  public :: euler_trafoFluxPressure2d
  public :: euler_trafoFluxVelocity2d
  public :: euler_trafoFluxMomentum2d
  public :: euler_trafoFluxDenEng2d
  public :: euler_trafoFluxDenPre2d
  public :: euler_trafoFluxDenPreVel2d
  public :: euler_trafoDiffDensity2d
  public :: euler_trafoDiffEnergy2d
  public :: euler_trafoDiffPressure2d
  public :: euler_trafoDiffVelocity2d
  public :: euler_trafoDiffMomentum2d
  public :: euler_trafoDiffDenEng2d
  public :: euler_trafoDiffDenPre2d
  public :: euler_trafoDiffDenPreVel2d
  public :: euler_calcBoundaryvalues2d
  public :: euler_hadaptCallbackScalar2d
  public :: euler_hadaptCallbackBlock2d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxGalerkin2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the standard
    ! Galerkin discretization in 2D.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP), dimension(NVAR2D) :: dF1_i, dF2_i, dF1_j, dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
#endif
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j

    !---------------------------------------------------------------------------
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
    ! ---------------------------------------------------------------------------

    ! Compute velocities
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); rv2i = vi*U_i(3)
    ru2j = uj*U_j(2); rv2j = vj*U_j(3)

#ifdef USE_EULER_IBP
    ! Compute fluxes for x-direction
    dF1_i(1) = U_i(2)
    dF1_i(2) = G1*U_i(4)-G14*ru2i-G2*rv2i
    dF1_i(3) = U_i(3)*ui
    dF1_i(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*ui

    dF1_j(1) = U_j(2)
    dF1_j(2) = G1*U_j(4)-G14*ru2j-G2*rv2j
    dF1_j(3) = U_j(3)*uj
    dF1_j(4) = (GAMMA*U_j(4)-G2*(ru2j+rv2j))*uj

    ! Compute fluxes for y-direction
    dF2_i(1) = U_i(3)
    dF2_i(2) = U_i(3)*ui
    dF2_i(3) = G1*U_i(4)-G14*rv2i-G2*ru2i
    dF2_i(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*vi

    dF2_j(1) = U_j(3)
    dF2_j(2) = U_j(3)*uj
    dF2_j(3) = (G1*U_j(4)-G14*rv2j-G2*ru2j)
    dF2_j(4) = (GAMMA*U_j(4)-G2*(ru2j+rv2j))*vj

    ! Assembly fluxes
    F_ij = dscale * ( C_ji(1)*dF1_j + C_ji(2)*dF2_j &
                     -C_ij(1)*dF1_i - C_ij(2)*dF2_i)
    F_ji = -F_ij
#else
    ! Compute flux difference for x-direction
    dF1_ij(1) = U_i(2)                           - U_j(2)
    dF1_ij(2) = G1*U_i(4)-G14*ru2i-G2*rv2i       - (G1*U_j(4)-G14*ru2j-G2*rv2j)
    dF1_ij(3) = U_i(3)*ui                        - U_j(3)*uj
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*ui - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*uj

    ! Compute flux difference for y-direction
    dF2_ij(1) = U_i(3)                           - U_j(3)
    dF2_ij(2) = U_i(3)*ui                        - U_j(3)*uj
    dF2_ij(3) = G1*U_i(4)-G14*rv2i-G2*ru2i       - (G1*U_j(4)-G14*rv2j-G2*ru2j)
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*vi - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*vj

    ! Assembly fluxes
    F_ij =   dscale * ( C_ij(1)*dF1_ij + C_ij(2)*dF2_ij)
    F_ji = - dscale * ( C_ji(1)*dF1_ij + C_ji(2)*dF2_ij)
#endif

  end subroutine euler_calcFluxGalerkin2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxGalerkinNoBdr2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the TVD
    ! discretization in 2D. The symmetric boundary contributions
    ! are neglected and incorporated in the antidiffusive flux.
    ! Hence, this is simply the standard Galerkin flux for the
    ! skew-symmetric internal contributions.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j

    !---------------------------------------------------------------------------
    ! Evaluate the Galerkin fluxes
    ! For a detailed description of algorithm and the definition of auxiliary
    ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------

    ! Compute velocities
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); rv2i = vi*U_i(3)
    ru2j = uj*U_j(2); rv2j = vj*U_j(3)

    ! Compute flux difference for x-direction
    dF1_ij(1) = U_i(2)                           - U_j(2)
    dF1_ij(2) = G1*U_i(4)-G14*ru2i-G2*rv2i       - (G1*U_j(4)-G14*ru2j-G2*rv2j)
    dF1_ij(3) = U_i(3)*ui                        - U_j(3)*uj
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*ui - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*uj

    ! Compute flux difference for y-direction
    dF2_ij(1) = U_i(3)                           - U_j(3)
    dF2_ij(2) = U_i(3)*ui                        - U_j(3)*uj
    dF2_ij(3) = G1*U_i(4)-G14*rv2i-G2*ru2i       - (G1*U_j(4)-G14*rv2j-G2*ru2j)
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*vi - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*vj

    ! Compute skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji)

    ! Assembly fluxes and exploit skew-symmetry of a_ij and F_ij
    F_ij = dscale * (a(1)*dF1_ij + a(2)*dF2_ij)
    F_ji = F_ij

  end subroutine euler_calcFluxGalerkinNoBdr2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxScalarDiss2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using scalar dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP), dimension(NVAR2D) :: dF1_i, dF2_i, dF1_j, dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j
    real(DP) :: d_ij,hi,hj,H_ij,q_ij,u_ij,v_ij,aux,vel,cs

    !---------------------------------------------------------------------------
    ! Evaluate the Galerkin fluxes
    ! For a detailed description of algorithm and the definition of auxiliary
    ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------

    ! Compute velocities
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); rv2i = vi*U_i(3)
    ru2j = uj*U_j(2); rv2j = vj*U_j(3)

#ifdef USE_EULER_IBP
    ! Compute fluxes for x-direction
    dF1_i(1) = U_i(2)
    dF1_i(2) = G1*U_i(4)-G14*ru2i-G2*rv2i
    dF1_i(3) = U_i(3)*ui
    dF1_i(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*ui

    dF1_j(1) = U_j(2)
    dF1_j(2) = G1*U_j(4)-G14*ru2j-G2*rv2j
    dF1_j(3) = U_j(3)*uj
    dF1_j(4) = (GAMMA*U_j(4)-G2*(ru2j+rv2j))*uj

    ! Compute fluxes for y-direction
    dF2_i(1) = U_i(3)
    dF2_i(2) = U_i(3)*ui
    dF2_i(3) = G1*U_i(4)-G14*rv2i-G2*ru2i
    dF2_i(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*vi

    dF2_j(1) = U_j(3)
    dF2_j(2) = U_j(3)*uj
    dF2_j(3) = (G1*U_j(4)-G14*rv2j-G2*ru2j)
    dF2_j(4) = (GAMMA*U_j(4)-G2*(ru2j+rv2j))*vj

    ! Assembly fluxes
    F_ij = dscale * ( C_ji(1)*dF1_j + C_ji(2)*dF2_j &
                     -C_ij(1)*dF1_i - C_ij(2)*dF2_i)
    F_ji = -F_ij
#else
    ! Compute flux difference for x-direction
    dF1_ij(1) = U_i(2)                           - U_j(2)
    dF1_ij(2) = G1*U_i(4)-G14*ru2i-G2*rv2i       - (G1*U_j(4)-G14*ru2j-G2*rv2j)
    dF1_ij(3) = U_i(3)*ui                        - U_j(3)*uj
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*ui - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*uj

    ! Compute flux difference for y-direction
    dF2_ij(1) = U_i(3)                           - U_j(3)
    dF2_ij(2) = U_i(3)*ui                        - U_j(3)*uj
    dF2_ij(3) = G1*U_i(4)-G14*rv2i-G2*ru2i       - (G1*U_j(4)-G14*rv2j-G2*ru2j)
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*vi - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*vj

    ! Assembly fluxes
    F_ij =   dscale * ( C_ij(1)*dF1_ij + C_ij(2)*dF2_ij)
    F_ji = - dscale * ( C_ji(1)*dF1_ij + C_ji(2)*dF2_ij)
#endif

    !---------------------------------------------------------------------------
    ! Evaluate the issipation
    !---------------------------------------------------------------------------

    ! Compute skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji)

    ! Compute Roe mean values
    aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
    u_ij = (aux*ui+uj)/(aux+1.0_DP)
    v_ij = (aux*vi+vj)/(aux+1.0_DP)
    hi   = GAMMA*U_i(4)/U_i(1)-G2*(U_i(2)*U_i(2)+U_i(3)*U_i(3))/(U_i(1)*U_i(1))
    hj   = GAMMA*U_j(4)/U_j(1)-G2*(U_j(2)*U_j(2)+U_j(3)*U_j(3))/(U_j(1)*U_j(1))
    H_ij = (aux*hi+hj)/(aux+1.0_DP)

    ! Compute auxiliary variables
    aux  = sqrt(a(1)*a(1)+a(2)*a(2))
    vel  = u_ij*a(1) + v_ij*a(2)
    q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
    cs   = sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL))

    ! Scalar dissipation
    d_ij = dscale * (abs(vel) + aux*cs)

    ! Multiply the solution difference by the artificial diffusion factor
    Diff = d_ij*(U_j-U_i)

    ! Add the artificial diffusion to the fluxes
    F_ij = F_ij+Diff
    F_ji = F_ji-Diff

  end subroutine euler_calcFluxScalarDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxScalarDissDiSp2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using scalar dissipation,
    ! whereby dimensional splitting is employed.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP), dimension(NVAR2D) :: dF1_i, dF2_i, dF1_j, dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j
    real(DP) :: d_ij,hi,hj,H_ij,q_ij,u_ij,v_ij,aux

    !---------------------------------------------------------------------------
    ! Evaluate the Galerkin fluxes
    ! For a detailed description of algorithm and the definition of auxiliary
    ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------

    ! Compute velocities
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); rv2i = vi*U_i(3)
    ru2j = uj*U_j(2); rv2j = vj*U_j(3)

#ifdef USE_EULER_IBP
    ! Compute fluxes for x-direction
    dF1_i(1) = U_i(2)
    dF1_i(2) = G1*U_i(4)-G14*ru2i-G2*rv2i
    dF1_i(3) = U_i(3)*ui
    dF1_i(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*ui

    dF1_j(1) = U_j(2)
    dF1_j(2) = G1*U_j(4)-G14*ru2j-G2*rv2j
    dF1_j(3) = U_j(3)*uj
    dF1_j(4) = (GAMMA*U_j(4)-G2*(ru2j+rv2j))*uj

    ! Compute fluxes for y-direction
    dF2_i(1) = U_i(3)
    dF2_i(2) = U_i(3)*ui
    dF2_i(3) = G1*U_i(4)-G14*rv2i-G2*ru2i
    dF2_i(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*vi

    dF2_j(1) = U_j(3)
    dF2_j(2) = U_j(3)*uj
    dF2_j(3) = (G1*U_j(4)-G14*rv2j-G2*ru2j)
    dF2_j(4) = (GAMMA*U_j(4)-G2*(ru2j+rv2j))*vj

    ! Assembly fluxes
    F_ij = dscale * ( C_ji(1)*dF1_j + C_ji(2)*dF2_j &
                     -C_ij(1)*dF1_i - C_ij(2)*dF2_i)
    F_ji = -F_ij
#else
    ! Compute flux difference for x-direction
    dF1_ij(1) = U_i(2)                           - U_j(2)
    dF1_ij(2) = G1*U_i(4)-G14*ru2i-G2*rv2i       - (G1*U_j(4)-G14*ru2j-G2*rv2j)
    dF1_ij(3) = U_i(3)*ui                        - U_j(3)*uj
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*ui - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*uj

    ! Compute flux difference for y-direction
    dF2_ij(1) = U_i(3)                           - U_j(3)
    dF2_ij(2) = U_i(3)*ui                        - U_j(3)*uj
    dF2_ij(3) = G1*U_i(4)-G14*rv2i-G2*ru2i       - (G1*U_j(4)-G14*rv2j-G2*ru2j)
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*vi - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*vj

    ! Assembly fluxes
    F_ij =   dscale * ( C_ij(1)*dF1_ij + C_ij(2)*dF2_ij)
    F_ji = - dscale * ( C_ji(1)*dF1_ij + C_ji(2)*dF2_ij)
#endif

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji)

    ! Compute Roe mean values
    aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
    u_ij = (aux*ui+uj)/(aux+1.0_DP)
    v_ij = (aux*vi+vj)/(aux+1.0_DP)
    hi   = GAMMA*U_i(4)/U_i(1)-G2*(U_i(2)*U_i(2)+U_i(3)*U_i(3))/(U_i(1)*U_i(1))
    hj   = GAMMA*U_j(4)/U_j(1)-G2*(U_j(2)*U_j(2)+U_j(3)*U_j(3))/(U_j(1)*U_j(1))
    H_ij = (aux*hi+hj)/(aux+1.0_DP)

    ! Compute auxiliary variable
    q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
    aux  = sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL))

    ! Scalar dissipation for x- and y-direction
    d_ij = dscale * ( abs(a(1)*u_ij) + abs(a(1))*aux +&
                      abs(a(2)*v_ij) + abs(a(2))*aux )

    ! Multiply the solution difference by the artificial diffusion factor
    Diff = d_ij*(U_j-U_i)

    ! Add the artificial diffusion to the fluxes
    F_ij = F_ij+Diff
    F_ji = F_ji-Diff

  end subroutine euler_calcFluxScalarDissDiSp2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxTensorDiss2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using tensorial dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP), dimension(NVAR2D) :: dF1_i, dF2_i, dF1_j, dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j
    real(DP) :: aux,aux1,aux2,uPow2,vPow2,hi,hj,H_ij,q_ij,u_ij,v_ij
    real(DP) :: anorm,l1,l2,l3,l4,w1,w2,w3,w4,cPow2,cs


    !---------------------------------------------------------------------------
    ! Evaluate the Galerkin fluxes
    ! For a detailed description of algorithm and the definition of auxiliary
    ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------

    ! Compute velocities
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); rv2i = vi*U_i(3)
    ru2j = uj*U_j(2); rv2j = vj*U_j(3)

#ifdef USE_EULER_IBP
    ! Compute fluxes for x-direction
    dF1_i(1) = U_i(2)
    dF1_i(2) = G1*U_i(4)-G14*ru2i-G2*rv2i
    dF1_i(3) = U_i(3)*ui
    dF1_i(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*ui

    dF1_j(1) = U_j(2)
    dF1_j(2) = G1*U_j(4)-G14*ru2j-G2*rv2j
    dF1_j(3) = U_j(3)*uj
    dF1_j(4) = (GAMMA*U_j(4)-G2*(ru2j+rv2j))*uj

    ! Compute fluxes for y-direction
    dF2_i(1) = U_i(3)
    dF2_i(2) = U_i(3)*ui
    dF2_i(3) = G1*U_i(4)-G14*rv2i-G2*ru2i
    dF2_i(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*vi

    dF2_j(1) = U_j(3)
    dF2_j(2) = U_j(3)*uj
    dF2_j(3) = (G1*U_j(4)-G14*rv2j-G2*ru2j)
    dF2_j(4) = (GAMMA*U_j(4)-G2*(ru2j+rv2j))*vj

    ! Assembly fluxes
    F_ij = dscale * ( C_ji(1)*dF1_j + C_ji(2)*dF2_j &
                     -C_ij(1)*dF1_i - C_ij(2)*dF2_i)
    F_ji = -F_ij
#else
    ! Compute flux difference for x-direction
    dF1_ij(1) = U_i(2)                           - U_j(2)
    dF1_ij(2) = G1*U_i(4)-G14*ru2i-G2*rv2i       - (G1*U_j(4)-G14*ru2j-G2*rv2j)
    dF1_ij(3) = U_i(3)*ui                        - U_j(3)*uj
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*ui - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*uj

    ! Compute flux difference for y-direction
    dF2_ij(1) = U_i(3)                           - U_j(3)
    dF2_ij(2) = U_i(3)*ui                        - U_j(3)*uj
    dF2_ij(3) = G1*U_i(4)-G14*rv2i-G2*ru2i       - (G1*U_j(4)-G14*rv2j-G2*ru2j)
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*vi - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*vj

    ! Assembly fluxes
    F_ij =   dscale * ( C_ij(1)*dF1_ij + C_ij(2)*dF2_ij)
    F_ji = - dscale * ( C_ji(1)*dF1_ij + C_ji(2)*dF2_ij)
#endif

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute the skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji); anorm = sqrt(a(1)*a(1)+a(2)*a(2))

    if (anorm .gt. SYS_EPSREAL) then

      ! Normalize the skew-symmetric coefficient
      a = a/anorm

      ! Compute Roe mean values
      aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
      u_ij = (aux*ui+uj)/(aux+1.0_DP)
      v_ij = (aux*vi+vj)/(aux+1.0_DP)
      hi   = GAMMA*U_i(4)/U_i(1)-G2*(U_i(2)*U_i(2)+U_i(3)*U_i(3))/(U_i(1)*U_i(1))
      hj   = GAMMA*U_j(4)/U_j(1)-G2*(U_j(2)*U_j(2)+U_j(3)*U_j(3))/(U_j(1)*U_j(1))
      H_ij = (aux*hi+hj)/(aux+1.0_DP)

      ! Compute auxiliary variables
      aux   = u_ij*a(1) + v_ij*a(2)
      uPow2 = u_ij*u_ij
      vPow2 = v_ij*v_ij
      q_ij  = 0.5_DP*(uPow2+vPow2)
      cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
      cs = sqrt(cPow2)

      ! Compute eigenvalues
      l1 = abs(aux-cs)
      l2 = abs(aux)
      l3 = abs(aux+cs)
      l4 = abs(aux)

      ! Compute solution difference U_j-U_i
      Diff = U_j-U_i

      ! Compute auxiliary quantities for characteristic variables
      aux1 = G2/cPow2*(q_ij*Diff(1)-u_ij*Diff(2)-v_ij*Diff(3)+Diff(4))
      aux2 = 0.5_DP*(aux*Diff(1)-a(1)*Diff(2)-a(2)*Diff(3))/cs

      ! Compute characteristic variables multiplied by the corresponding eigenvalue
      w1 = l1 * (aux1 + aux2)
      w2 = l2 * ((1.0_DP-G1*q_ij/cPow2)*Diff(1)+G1*(u_ij*Diff(2)+v_ij*Diff(3)-Diff(4))/cPow2)
      w3 = l3 * (aux1 - aux2)
      w4 = l4 * ((a(1)*v_ij-a(2)*u_ij)*Diff(1)+a(2)*Diff(2)-a(1)*Diff(3))

      ! Compute "anorm * dscale"
      anorm = anorm*dscale

      ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
      Diff(1) = anorm * ( w1 + w2 + w3 )
      Diff(2) = anorm * ( (u_ij-cs*a(1))*w1 + u_ij*w2 + (u_ij+cs*a(1))*w3 + a(2)*w4 )
      Diff(3) = anorm * ( (v_ij-cs*a(2))*w1 + v_ij*w2 + (v_ij+cs*a(2))*w3 - a(1)*w4 )
      Diff(4) = anorm * ( (H_ij-cs*aux)*w1  + q_ij*w2 + (H_ij+cs*aux)*w3  + (u_ij*a(2)-v_ij*a(1))*w4 )

      ! Add the artificial diffusion to the fluxes
      F_ij = F_ij+Diff
      F_ji = F_ji-Diff

    end if

  end subroutine euler_calcFluxTensorDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxTensorDissDiSp2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using tensorial dissipation,
    ! whereby dimensional splitting is employed.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP), dimension(NVAR2D) :: dF1_i, dF2_i, dF1_j, dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j
    real(DP) :: aux,aux1,aux2,uPow2,vPow2,hi,hj,H_ij,q_ij,u_ij,v_ij
    real(DP) :: anorm,l1,l2,l3,l4,w1,w2,w3,w4,cPow2,cs


    !---------------------------------------------------------------------------
    ! Evaluate the Galerkin fluxes
    ! For a detailed description of algorithm and the definition of auxiliary
    ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------

    ! Compute velocities
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); rv2i = vi*U_i(3)
    ru2j = uj*U_j(2); rv2j = vj*U_j(3)

#ifdef USE_EULER_IBP
    ! Compute fluxes for x-direction
    dF1_i(1) = U_i(2)
    dF1_i(2) = G1*U_i(4)-G14*ru2i-G2*rv2i
    dF1_i(3) = U_i(3)*ui
    dF1_i(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*ui

    dF1_j(1) = U_j(2)
    dF1_j(2) = G1*U_j(4)-G14*ru2j-G2*rv2j
    dF1_j(3) = U_j(3)*uj
    dF1_j(4) = (GAMMA*U_j(4)-G2*(ru2j+rv2j))*uj

    ! Compute fluxes for y-direction
    dF2_i(1) = U_i(3)
    dF2_i(2) = U_i(3)*ui
    dF2_i(3) = G1*U_i(4)-G14*rv2i-G2*ru2i
    dF2_i(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*vi

    dF2_j(1) = U_j(3)
    dF2_j(2) = U_j(3)*uj
    dF2_j(3) = (G1*U_j(4)-G14*rv2j-G2*ru2j)
    dF2_j(4) = (GAMMA*U_j(4)-G2*(ru2j+rv2j))*vj

    ! Assembly fluxes
    F_ij = dscale * ( C_ji(1)*dF1_j + C_ji(2)*dF2_j &
                     -C_ij(1)*dF1_i - C_ij(2)*dF2_i)
    F_ji = -F_ij
#else
    ! Compute flux difference for x-direction
    dF1_ij(1) = U_i(2)                           - U_j(2)
    dF1_ij(2) = G1*U_i(4)-G14*ru2i-G2*rv2i       - (G1*U_j(4)-G14*ru2j-G2*rv2j)
    dF1_ij(3) = U_i(3)*ui                        - U_j(3)*uj
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*ui - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*uj

    ! Compute flux difference for y-direction
    dF2_ij(1) = U_i(3)                           - U_j(3)
    dF2_ij(2) = U_i(3)*ui                        - U_j(3)*uj
    dF2_ij(3) = G1*U_i(4)-G14*rv2i-G2*ru2i       - (G1*U_j(4)-G14*rv2j-G2*ru2j)
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*vi - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*vj

    ! Assembly fluxes
    F_ij =   dscale * ( C_ij(1)*dF1_ij + C_ij(2)*dF2_ij)
    F_ji = - dscale * ( C_ji(1)*dF1_ij + C_ji(2)*dF2_ij)
#endif

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute the skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji); anorm = sqrt(a(1)*a(1)+a(2)*a(2))

    if (anorm .gt. SYS_EPSREAL) then

      ! Compute the absolute value
      a = abs(a)

      ! Compute Roe mean values
      aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
      u_ij = (aux*ui+uj)/(aux+1.0_DP)
      v_ij = (aux*vi+vj)/(aux+1.0_DP)
      hi   = GAMMA*U_i(4)/U_i(1)-G2*(U_i(2)*U_i(2)+U_i(3)*U_i(3))/(U_i(1)*U_i(1))
      hj   = GAMMA*U_j(4)/U_j(1)-G2*(U_j(2)*U_j(2)+U_j(3)*U_j(3))/(U_j(1)*U_j(1))
      H_ij = (aux*hi+hj)/(aux+1.0_DP)

      ! Compute auxiliary variables
      uPow2 = u_ij*u_ij
      vPow2 = v_ij*v_ij
      q_ij  = 0.5_DP*(uPow2+vPow2)
      cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
      cs = sqrt(cPow2)

      !-------------------------------------------------------------------------
      ! Dimensional splitting: x-direction
      !-------------------------------------------------------------------------

      ! Compute eigenvalues
      l1 = abs(u_ij-cs)
      l2 = abs(u_ij)
      l3 = abs(u_ij+cs)
      l4 = abs(u_ij)

      ! Compute solution difference U_j-U_i
      Diff = U_j-U_i

      ! Compute auxiliary quantities for characteristic variables
      aux1 = G2/cPow2*(q_ij*Diff(1)-u_ij*Diff(2)-v_ij*Diff(3)+Diff(4))
      aux2 = 0.5_DP*(u_ij*Diff(1)-Diff(2))/cs

      ! Compute characteristic variables multiplied by the corresponding eigenvalue
      w1 = l1 * (aux1 + aux2)
      w2 = l2 * ((1.0_DP-G1*q_ij/cPow2)*Diff(1)+G1*(u_ij*Diff(2)+v_ij*Diff(3)-Diff(4))/cPow2)
      w3 = l3 * (aux1 - aux2)
      w4 = l4 * (v_ij*Diff(1)-Diff(3))

      ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
      Diff(1) = dscale * a(1) * ( w1 + w2 + w3 )
      Diff(2) = dscale * a(1) * ( (u_ij-cs)*w1 + u_ij*w2 + (u_ij+cs)*w3 )
      Diff(3) = dscale * a(1) * ( v_ij*w1 + v_ij*w2 + v_ij*w3 - w4 )
      Diff(4) = dscale * a(1) * ( (H_ij-cs*u_ij)*w1  + q_ij*w2 + (H_ij+cs*u_ij)*w3  -v_ij*w4 )

      ! Add the artificial diffusion to the fluxes
      F_ij = F_ij+Diff
      F_ji = F_ji-Diff


      !-------------------------------------------------------------------------
      ! Dimensional splitting: y-direction
      !-------------------------------------------------------------------------
      ! Compute eigenvalues
      l1 = abs(v_ij-cs)
      l2 = abs(v_ij)
      l3 = abs(v_ij+cs)
      l4 = abs(v_ij)

      ! Compute solution difference U_j-U_i
      Diff = U_j-U_i

      ! Compute auxiliary quantities for characteristic variables
      aux1 = G2/cPow2*(q_ij*Diff(1)-u_ij*Diff(2)-v_ij*Diff(3)+Diff(4))
      aux2 = 0.5_DP*(v_ij*Diff(1)-Diff(3))/cs

      ! Compute characteristic variables multiplied by the corresponding eigenvalue
      w1 = l1 * (aux1 + aux2)
      w2 = l2 * ((1.0_DP-G1*q_ij/cPow2)*Diff(1)+G1*(u_ij*Diff(2)+v_ij*Diff(3)-Diff(4))/cPow2)
      w3 = l3 * (aux1 - aux2)
      w4 = l4 * (-u_ij*Diff(1)+Diff(2))

      ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
      Diff(1) = dscale * a(2) * ( w1 + w2 + w3 )
      Diff(2) = dscale * a(2) * ( u_ij*w1 + u_ij*w2 + u_ij*w3 + w4 )
      Diff(3) = dscale * a(2) * ( (v_ij-cs)*w1 + v_ij*w2 + (v_ij+cs)*w3 )
      Diff(4) = dscale * a(2) * ( (H_ij-cs*v_ij)*w1  + q_ij*w2 + (H_ij+cs*v_ij)*w3  + u_ij*w4 )

      ! Add the artificial diffusion to the fluxes
      F_ij = F_ij+Diff
      F_ji = F_ji-Diff

    end if

  end subroutine euler_calcFluxTensorDissDiSp2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxRusanov2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using the Rusanov dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP), dimension(NVAR2D) :: dF1_i, dF2_i, dF1_j, dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j
    real(DP) :: d_ij,ci,cj,Ei,Ej

    !---------------------------------------------------------------------------
    ! Evaluate the Galerkin fluxes
    ! For a detailed description of algorithm and the definition of auxiliary
    ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------

    ! Compute velocities and energy
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1); Ei = U_i(4)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1); Ej = U_j(4)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); rv2i = vi*U_i(3)
    ru2j = uj*U_j(2); rv2j = vj*U_j(3)

#ifdef USE_EULER_IBP
    ! Compute fluxes for x-direction
    dF1_i(1) = U_i(2)
    dF1_i(2) = G1*U_i(4)-G14*ru2i-G2*rv2i
    dF1_i(3) = U_i(3)*ui
    dF1_i(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*ui

    dF1_j(1) = U_j(2)
    dF1_j(2) = G1*U_j(4)-G14*ru2j-G2*rv2j
    dF1_j(3) = U_j(3)*uj
    dF1_j(4) = (GAMMA*U_j(4)-G2*(ru2j+rv2j))*uj

    ! Compute fluxes for y-direction
    dF2_i(1) = U_i(3)
    dF2_i(2) = U_i(3)*ui
    dF2_i(3) = G1*U_i(4)-G14*rv2i-G2*ru2i
    dF2_i(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*vi

    dF2_j(1) = U_j(3)
    dF2_j(2) = U_j(3)*uj
    dF2_j(3) = (G1*U_j(4)-G14*rv2j-G2*ru2j)
    dF2_j(4) = (GAMMA*U_j(4)-G2*(ru2j+rv2j))*vj

    ! Assembly fluxes
    F_ij = dscale * ( C_ji(1)*dF1_j + C_ji(2)*dF2_j &
                     -C_ij(1)*dF1_i - C_ij(2)*dF2_i)
    F_ji = -F_ij
#else
    ! Compute flux difference for x-direction
    dF1_ij(1) = U_i(2)                           - U_j(2)
    dF1_ij(2) = G1*U_i(4)-G14*ru2i-G2*rv2i       - (G1*U_j(4)-G14*ru2j-G2*rv2j)
    dF1_ij(3) = U_i(3)*ui                        - U_j(3)*uj
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*ui - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*uj

    ! Compute flux difference for y-direction
    dF2_ij(1) = U_i(3)                           - U_j(3)
    dF2_ij(2) = U_i(3)*ui                        - U_j(3)*uj
    dF2_ij(3) = G1*U_i(4)-G14*rv2i-G2*ru2i       - (G1*U_j(4)-G14*rv2j-G2*ru2j)
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*vi - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*vj

    ! Assembly fluxes
    F_ij =   dscale * ( C_ij(1)*dF1_ij + C_ij(2)*dF2_ij)
    F_ji = - dscale * ( C_ji(1)*dF1_ij + C_ji(2)*dF2_ij)
#endif

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

!!$    ! Compute enthalpy
!!$    hi = GAMMA*Ei-G2*(ui*ui+vi*vi)
!!$    hj = GAMMA*Ej-G2*(uj*uj+vj*vj)
!!$
!!$    ! Compute the speed of sound
!!$    ci = sqrt(max((GAMMA-1)*(hi-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
!!$    cj = sqrt(max((GAMMA-1)*(hj-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))

    ! Compute the speed of sound
    ci = sqrt(max(G15*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
    cj = sqrt(max(G15*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))

    ! Scalar dissipation for the Rusanov flux
    d_ij = max( abs(C_ij(1)*uj+C_ij(2)*vj) +&
                sqrt(C_ij(1)*C_ij(1)+C_ij(2)*C_ij(2))*cj,&
                abs(C_ji(1)*ui+C_ji(2)*vi) +&
                sqrt(C_ji(1)*C_ji(1)+C_ji(2)*C_ji(2))*ci )

    ! Multiply the solution difference by the artificial diffusion factor
    Diff = dscale * d_ij*(U_j-U_i)

    ! Add the artificial diffusion to the fluxes
    F_ij = F_ij+Diff
    F_ji = F_ji-Diff

  end subroutine euler_calcFluxRusanov2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcFluxRusanovDiSp2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using the Rusanov dissipation,
    ! whereby dimensional splitting is employed.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP), dimension(NVAR2D) :: dF1_i, dF2_i, dF1_j, dF2_j
#else
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
#endif
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: ui,vi,uj,vj,ru2i,ru2j,rv2i,rv2j
    real(DP) :: d_ij,ci,cj,Ei,Ej

    !---------------------------------------------------------------------------
    ! Evaluate the Galerkin fluxes
    ! For a detailed description of algorithm and the definition of auxiliary
    ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------

    ! Compute velocities and energy
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1); Ei = U_i(4)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1); Ej = U_j(4)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); rv2i = vi*U_i(3)
    ru2j = uj*U_j(2); rv2j = vj*U_j(3)

#ifdef USE_EULER_IBP
    ! Compute fluxes for x-direction
    dF1_i(1) = U_i(2)
    dF1_i(2) = G1*U_i(4)-G14*ru2i-G2*rv2i
    dF1_i(3) = U_i(3)*ui
    dF1_i(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*ui

    dF1_j(1) = U_j(2)
    dF1_j(2) = G1*U_j(4)-G14*ru2j-G2*rv2j
    dF1_j(3) = U_j(3)*uj
    dF1_j(4) = (GAMMA*U_j(4)-G2*(ru2j+rv2j))*uj

    ! Compute fluxes for y-direction
    dF2_i(1) = U_i(3)
    dF2_i(2) = U_i(3)*ui
    dF2_i(3) = G1*U_i(4)-G14*rv2i-G2*ru2i
    dF2_i(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*vi

    dF2_j(1) = U_j(3)
    dF2_j(2) = U_j(3)*uj
    dF2_j(3) = (G1*U_j(4)-G14*rv2j-G2*ru2j)
    dF2_j(4) = (GAMMA*U_j(4)-G2*(ru2j+rv2j))*vj

    ! Assembly fluxes
    F_ij = dscale * ( C_ji(1)*dF1_j + C_ji(2)*dF2_j &
                     -C_ij(1)*dF1_i - C_ij(2)*dF2_i)
    F_ji = -F_ij
#else
    ! Compute flux difference for x-direction
    dF1_ij(1) = U_i(2)                           - U_j(2)
    dF1_ij(2) = G1*U_i(4)-G14*ru2i-G2*rv2i       - (G1*U_j(4)-G14*ru2j-G2*rv2j)
    dF1_ij(3) = U_i(3)*ui                        - U_j(3)*uj
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*ui - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*uj

    ! Compute flux difference for y-direction
    dF2_ij(1) = U_i(3)                           - U_j(3)
    dF2_ij(2) = U_i(3)*ui                        - U_j(3)*uj
    dF2_ij(3) = G1*U_i(4)-G14*rv2i-G2*ru2i       - (G1*U_j(4)-G14*rv2j-G2*ru2j)
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(ru2i+rv2i))*vi - (GAMMA*U_j(4)-G2*(ru2j+rv2j))*vj

    ! Assembly fluxes
    F_ij =   dscale * ( C_ij(1)*dF1_ij + C_ij(2)*dF2_ij)
    F_ji = - dscale * ( C_ji(1)*dF1_ij + C_ji(2)*dF2_ij)
#endif

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

!!$    ! Compute enthalpy
!!$    hi = GAMMA*Ei+(1-GAMMA)*0.5*(ui*ui+vi*vi)
!!$    hj = GAMMA*Ej+(1-GAMMA)*0.5*(uj*uj+vj*vj)
!!$
!!$    ! Compute speed of sound
!!$    ci = sqrt(max((GAMMA-1)*(hi-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
!!$    cj = sqrt(max((GAMMA-1)*(hj-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))

    ! Compute the speed of sound
    ci = sqrt(max(G15*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
    cj = sqrt(max(G15*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))

    ! Scalar dissipation
    d_ij = max( abs(C_ij(1)*uj) + abs(C_ij(1))*cj,&
                abs(C_ji(1)*ui) + abs(C_ji(1))*ci )&
         + max( abs(C_ij(2)*vj) + abs(C_ij(2))*cj,&
                abs(C_ji(2)*vi) + abs(C_ji(2))*ci )

    ! Multiply the solution difference by the artificial diffusion factor
    Diff = dscale * d_ij*(U_j-U_i)

    ! Add the artificial diffusion to the fluxes
    F_ij = F_ij+Diff
    F_ji = F_ji-Diff

  end subroutine euler_calcFluxRusanovDiSp2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixDiagonalDiag2d(U_i, C_ii, i, dscale, K_ii)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 2D
!</description>

!<input>
    ! local solution at node I
    real(DP), dimension(:), intent(in) :: U_i

    ! coefficients from spatial discretization
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
    real(DP) :: ui,vi

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1)

    ! Compute Galerkin coefficient K_ii
    K_ii(1) = 0.0_DP
    K_ii(2) = dscale * (G13*ui*C_ii(1)+vi*C_ii(2))
    K_ii(3) = dscale * (ui*C_ii(1)+G13*vi*C_ii(2))
    K_ii(4) = dscale * (GAMMA*(ui*C_ii(1)+vi*C_ii(2)))

  end subroutine euler_calcMatrixDiagonalDiag2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixDiagonal2d(U_i, C_ii, i, dscale, K_ii)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 2D
!</description>

!<input>
    ! local solution at node I
    real(DP), dimension(:), intent(in) :: U_i

    ! coefficients from spatial discretization
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
    real(DP) :: ui,vi,qi,Ei,uvi,uPow2i,vPow2i,aux

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1);   Ei = U_i(4)/U_i(1)
    uvi = ui*vi;   qi = ui*ui+vi*vi;   uPow2i = ui*ui;   vPow2i = vi*vi
    aux = ui*C_ii(1)+vi*C_ii(2)

    ! Compute Galerkin coefficient K_ii
    K_ii( 1) = 0.0_DP
    K_ii( 2) = dscale * ((G2*qi-uPow2i)*C_ii(1)-uvi*C_ii(2))
    K_ii( 3) = dscale * ((G2*qi-vPow2i)*C_ii(2)-uvi*C_ii(1))
    K_ii( 4) = dscale * (G1*qi-GAMMA*Ei)*aux

    K_ii( 5) = dscale * C_ii(1)
    K_ii( 6) = dscale * (G13*ui*C_ii(1)+vi*C_ii(2))
    K_ii( 7) = dscale * (vi*C_ii(1)-G1*ui*C_ii(2))
    K_ii( 8) = dscale * ((GAMMA*Ei-G2*qi)*C_ii(1)-G1*ui*aux)

    K_ii( 9) = dscale * C_ii(2)
    K_ii(10) = dscale * (ui*C_ii(2)-G1*vi*C_ii(1))
    K_ii(11) = dscale * (ui*C_ii(1)+G13*vi*C_ii(2))
    K_ii(12) = dscale * ((GAMMA*Ei-G2*qi)*C_ii(2)-G1*vi*aux)

    K_ii(13) = 0.0_DP
    K_ii(14) = dscale * G1*C_ii(1)
    K_ii(15) = dscale * G1*C_ii(2)
    K_ii(16) = dscale * (GAMMA*(ui*C_ii(1)+vi*C_ii(2)))

  end subroutine euler_calcMatrixDiagonal2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixGalerkinDiag2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP) :: ui,uj,vi,vj

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1)

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale * (G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij(3) = dscale * (uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(4) = dscale * (GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale * (G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji(3) = dscale * (ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(4) = dscale * (GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

    ! Nullify dissipation tensor
    D_ij = 0.0_DP

  end subroutine euler_calcMatrixGalerkinDiag2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixGalerkin2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP) :: Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj,uPow2i,uPow2j,vPow2i,vPow2j,aux1,aux2

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1);   Ei = U_i(4)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1);   Ej = U_j(4)/U_j(1)

    uvi = ui*vi;   qi = ui*ui+vi*vi;   uPow2i = ui*ui;   vPow2i = vi*vi
    uvj = uj*vj;   qj = uj*uj+vj*vj;   uPow2j = uj*uj;   vPow2j = vj*vj

    aux1 = uj*C_ij(1)+vj*C_ij(2)
    aux2 = ui*C_ji(1)+vi*C_ji(2)

    ! Compute Galerkin coefficient K_ij
    K_ij( 1) = 0.0_DP
    K_ij( 2) = dscale * ((G2*qj-uPow2j)*C_ij(1)-uvj*C_ij(2))
    K_ij( 3) = dscale * ((G2*qj-vPow2j)*C_ij(2)-uvj*C_ij(1))
    K_ij( 4) = dscale * (G1*qj-GAMMA*Ej)*aux1

    K_ij( 5) = dscale * C_ij(1)
    K_ij( 6) = dscale * (G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij( 7) = dscale * (vj*C_ij(1)-G1*uj*C_ij(2))
    K_ij( 8) = dscale * ((GAMMA*Ej-G2*qj)*C_ij(1)-G1*uj*aux1)

    K_ij( 9) = dscale * C_ij(2)
    K_ij(10) = dscale * (uj*C_ij(2)-G1*vj*C_ij(1))
    K_ij(11) = dscale * (uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(12) = dscale * ((GAMMA*Ej-G2*qj)*C_ij(2)-G1*vj*aux1)

    K_ij(13) = 0.0_DP
    K_ij(14) = dscale * G1*C_ij(1)
    K_ij(15) = dscale * G1*C_ij(2)
    K_ij(16) = dscale * (GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji( 1) = 0.0_DP
    K_ji( 2) = dscale * ((G1*qi-uPow2i)*C_ji(1)-uvi*C_ji(2))
    K_ji( 3) = dscale * ((G1*qi-vPow2i)*C_ji(2)-uvi*C_ji(1))
    K_ji( 4) = dscale * (G1*qi-GAMMA*Ei)*aux2

    K_ji( 5) = dscale * C_ji(1)
    K_ji( 6) = dscale * (G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji( 7) = dscale * (vi*C_ji(1)-G1*ui*C_ji(2))
    K_ji( 8) = dscale * ((GAMMA*Ei-G2*qi)*C_ji(1)-G1*ui*aux2)

    K_ji( 9) = dscale * C_ji(2)
    K_ji(10) = dscale * (ui*C_ji(2)-G1*vi*C_ji(1))
    K_ji(11) = dscale * (ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(12) = dscale * ((GAMMA*Ei-G2*qi)*C_ji(2)-G1*vi*aux2)

    K_ji(13) = 0.0_DP
    K_ji(14) = dscale * G1*C_ji(1)
    K_ji(15) = dscale * G1*C_ji(2)
    K_ji(16) = dscale * (GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

    ! Nullify dissipation tensor
    D_ij = 0.0_DP

  end subroutine euler_calcMatrixGalerkin2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixScalarDissDiag2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies scalar artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: anorm,aux,hi,hj,H_ij,q_ij,ui,uj,u_ij,vi,vj,v_ij

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1)

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale * (G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij(3) = dscale * (uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(4) = dscale * (GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale * (G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji(3) = dscale * (ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(4) = dscale * (GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute skew-symmetric coefficient and its norm
    a = 0.5_DP*(C_ji-C_ij); anorm = sqrt(a(1)*a(1)+a(2)*a(2))

    if (anorm .gt. SYS_EPSREAL) then

      ! Compute Roe mean values
      aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
      u_ij = (aux*ui+uj)/(aux+1.0_DP)
      v_ij = (aux*vi+vj)/(aux+1.0_DP)
      hi   = GAMMA*U_i(4)/U_i(1)-G2*(ui*ui+vi*vi)
      hj   = GAMMA*U_j(4)/U_j(1)-G2*(uj*uj+vj*vj)
      H_ij = (aux*hi+hj)/(aux+1.0_DP)

      ! Compute auxiliary variables
      q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)

      ! Compute scalar dissipation
      D_ij = dscale * (abs(a(1)*u_ij+a(2)*v_ij) +&
                       anorm*sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL)))
    else

      ! Nullify dissipation tensor
      D_ij = 0.0_DP

    end if

  end subroutine euler_calcMatrixScalarDissDiag2d

!*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixScalarDiss2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies scalar artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: anorm,aux,aux1,aux2,hi,hj,H_ij,q_ij,u_ij,v_ij
    real(DP) :: Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj,uPow2i,uPow2j,vPow2i,vPow2j

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1);   Ei = U_i(4)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1);   Ej = U_j(4)/U_j(1)

    uvi = ui*vi;   qi = ui*ui+vi*vi;   uPow2i = ui*ui;   vPow2i = vi*vi
    uvj = uj*vj;   qj = uj*uj+vj*vj;   uPow2j = uj*uj;   vPow2j = vj*vj

    aux1 = uj*C_ij(1)+vj*C_ij(2)
    aux2 = ui*C_ji(1)+vi*C_ji(2)

    ! Compute Galerkin coefficient K_ij
    K_ij( 1) = 0.0_DP
    K_ij( 2) = dscale * ((G2*qj-uPow2j)*C_ij(1)-uvj*C_ij(2))
    K_ij( 3) = dscale * ((G2*qj-vPow2j)*C_ij(2)-uvj*C_ij(1))
    K_ij( 4) = dscale * (G1*qj-GAMMA*Ej)*aux1

    K_ij( 5) = dscale * C_ij(1)
    K_ij( 6) = dscale * (G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij( 7) = dscale * (vj*C_ij(1)-G1*uj*C_ij(2))
    K_ij( 8) = dscale * ((GAMMA*Ej-G2*qj)*C_ij(1)-G1*uj*aux1)

    K_ij( 9) = dscale * C_ij(2)
    K_ij(10) = dscale * (uj*C_ij(2)-G1*vj*C_ij(1))
    K_ij(11) = dscale * (uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(12) = dscale * ((GAMMA*Ej-G2*qj)*C_ij(2)-G1*vj*aux1)

    K_ij(13) = 0.0_DP
    K_ij(14) = dscale * G1*C_ij(1)
    K_ij(15) = dscale * G1*C_ij(2)
    K_ij(16) = dscale * (GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji( 1) = 0.0_DP
    K_ji( 2) = dscale * ((G1*qi-uPow2i)*C_ji(1)-uvi*C_ji(2))
    K_ji( 3) = dscale * ((G1*qi-vPow2i)*C_ji(2)-uvi*C_ji(1))
    K_ji( 4) = dscale * (G1*qi-GAMMA*Ei)*aux2

    K_ji( 5) = dscale * C_ji(1)
    K_ji( 6) = dscale * (G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji( 7) = dscale * (vi*C_ji(1)-G1*ui*C_ji(2))
    K_ji( 8) = dscale * ((GAMMA*Ei-G2*qi)*C_ji(1)-G1*ui*aux2)

    K_ji( 9) = dscale * C_ji(2)
    K_ji(10) = dscale * (ui*C_ji(2)-G1*vi*C_ji(1))
    K_ji(11) = dscale * (ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(12) = dscale * ((GAMMA*Ei-G2*qi)*C_ji(2)-G1*vi*aux2)

    K_ji(13) = 0.0_DP
    K_ji(14) = dscale * G1*C_ji(1)
    K_ji(15) = dscale * G1*C_ji(2)
    K_ji(16) = dscale * (GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute coefficients
    a = 0.5_DP*(C_ji-C_ij); anorm = sqrt(a(1)*a(1)+a(2)*a(2))

    if (anorm .gt. SYS_EPSREAL) then

      ! Compute Roe mean values
      aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
      u_ij = (aux*ui+uj)/(aux+1.0_DP)
      v_ij = (aux*vi+vj)/(aux+1.0_DP)
      hi   = GAMMA*U_i(4)/U_i(1)-G2*(ui*ui+vi*vi)
      hj   = GAMMA*U_j(4)/U_j(1)-G2*(uj*uj+vj*vj)
      H_ij = (aux*hi+hj)/(aux+1.0_DP)

      ! Compute auxiliary variables
      q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)

      ! Compute scalar dissipation
      aux = dscale * (abs(a(1)*u_ij+a(2)*v_ij) +&
                      anorm*sqrt(max(G1*(H_ij-q_ij), SYS_EPSREAL)))

      D_ij     = 0.0_DP
      D_ij( 1) = aux
      D_ij( 6) = aux
      D_ij(11) = aux
      D_ij(16) = aux

    else

      D_ij = 0.0_DP

    end if

  end subroutine euler_calcMatrixScalarDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixTensorDissDiag2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP), dimension(NVAR2D,NVAR2D) :: R_ij,L_ij
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: aux,hi,hj,H_ij,q_ij,ui,uj,u_ij,vi,vj,v_ij
    real(DP) :: l1,l2,l3,l4,anorm,c1,c2,cs,cPow2,vel

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1)

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale * (G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij(3) = dscale * (uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(4) = dscale * (GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale * (G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji(3) = dscale * (ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(4) = dscale * (GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute skew-symmetric coefficient and its norm
    a = 0.5_DP*(C_ji-C_ij); anorm = sqrt(a(1)*a(1)+a(2)*a(2))

    if (anorm .gt. SYS_EPSREAL) then

      ! Compute Roe mean values
      aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
      u_ij = (aux*ui+uj)/(aux+1.0_DP)
      v_ij = (aux*vi+vj)/(aux+1.0_DP)
      hi   = GAMMA*U_i(4)/U_i(1)-G2*(ui*ui+vi*vi)
      hj   = GAMMA*U_j(4)/U_j(1)-G2*(uj*uj+vj*vj)
      H_ij = (aux*hi+hj)/(aux+1.0_DP)

      ! Compute auxiliary values
      c1    = a(1)/anorm
      c2    = a(2)/anorm
      vel   = c1*u_ij+c2*v_ij
      q_ij  = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
      cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
      cs = sqrt(cPow2)

      ! Diagonal matrix of eigenvalues
      l1 = abs(vel-cs)
      l2 = abs(vel)
      l3 = abs(vel+cs)
      l4 = abs(vel)

      ! Matrix of right eigenvectors
      R_ij(1,1) =  l1
      R_ij(2,1) =  l1*(u_ij-cs*c1)
      R_ij(3,1) =  l1*(v_ij-cs*c2)
      R_ij(4,1) =  l1*(H_ij-cs*vel)

      R_ij(1,2) =  l2
      R_ij(2,2) =  l2*u_ij
      R_ij(3,2) =  l2*v_ij
      R_ij(4,2) =  l2*q_ij

      R_ij(1,3) =  l3
      R_ij(2,3) =  l3*(u_ij+cs*c1)
      R_ij(3,3) =  l3*(v_ij+cs*c2)
      R_ij(4,3) =  l3*(H_ij+cs*vel)

      R_ij(1,4) =  0.0_DP
      R_ij(2,4) =  l4*c2
      R_ij(3,4) = -l4*c1
      R_ij(4,4) =  l4*(u_ij*c2-v_ij*c1)

      ! Matrix of left eigenvectors
      L_ij(1,1) = 0.5_DP*(G1*q_ij+cs*vel)/cPow2
      L_ij(2,1) = (cPow2-G1*q_ij)/cPow2
      L_ij(3,1) = 0.5_DP*(G1*q_ij-cs*vel)/cPow2
      L_ij(4,1) = v_ij*c1-u_ij*c2

      L_ij(1,2) = 0.5_DP*(-G1*u_ij-cs*c1)/cPow2
      L_ij(2,2) = G1*u_ij/cPow2
      L_ij(3,2) = 0.5_DP*(-G1*u_ij+cs*c1)/cPow2
      L_ij(4,2) = c2

      L_ij(1,3) = 0.5_DP*(-G1*v_ij-cs*c2)/cPow2
      L_ij(2,3) = G1*v_ij/cPow2
      L_ij(3,3) = 0.5_DP*(-G1*v_ij+cs*c2)/cPow2
      L_ij(4,3) = -c1

      L_ij(1,4) =  G2/cPow2
      L_ij(2,4) = -G1/cPow2
      L_ij(3,4) =  G2/cPow2
      L_ij(4,4) =  0.0_DP


      ! Compute tensorial dissipation D_ij = diag(R_ij*|Lbd_ij|*L_ij)*I
      D_ij    = 0.0_DP
      D_ij(1) = anorm*( R_ij(1,1)*L_ij(1,1)+&
                        R_ij(1,2)*L_ij(2,1)+&
                        R_ij(1,3)*L_ij(3,1)+&
                        R_ij(1,4)*L_ij(4,1)  )
      D_ij(2) = anorm*( R_ij(2,1)*L_ij(1,2)+&
                        R_ij(2,2)*L_ij(2,2)+&
                        R_ij(2,3)*L_ij(3,2)+&
                        R_ij(2,4)*L_ij(4,2)  )
      D_ij(3) = anorm*( R_ij(3,1)*L_ij(1,3)+&
                        R_ij(3,2)*L_ij(2,3)+&
                        R_ij(3,3)*L_ij(3,3)+&
                        R_ij(3,4)*L_ij(4,3)  )
      D_ij(4) = anorm*( R_ij(4,1)*L_ij(1,4)+&
                        R_ij(4,2)*L_ij(2,4)+&
                        R_ij(4,3)*L_ij(3,4)+&
                        R_ij(4,4)*L_ij(4,4)  )

    else

      D_ij = 0.0_DP

    end if

  end subroutine euler_calcMatrixTensorDissDiag2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixTensorDiss2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP), dimension(NVAR2D,NVAR2D) :: R_ij,L_ij
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: anorm,aux,hi,hj,H_ij,q_ij,u_ij,v_ij,vel,c1,c2,cPow2,cs,l1,l2,l3,l4
    real(DP) :: Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj,uPow2i,uPow2j,vPow2i,vPow2j,aux1,aux2

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1);   Ei = U_i(4)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1);   Ej = U_j(4)/U_j(1)

    uvi = ui*vi;   qi = ui*ui+vi*vi;   uPow2i = ui*ui;   vPow2i = vi*vi
    uvj = uj*vj;   qj = uj*uj+vj*vj;   uPow2j = uj*uj;   vPow2j = vj*vj

    aux1 = uj*C_ij(1)+vj*C_ij(2)
    aux2 = ui*C_ji(1)+vi*C_ji(2)

    ! Compute Galerkin coefficient K_ij
    K_ij( 1) = 0.0_DP
    K_ij( 2) = dscale * ((G2*qj-uPow2j)*C_ij(1)-uvj*C_ij(2))
    K_ij( 3) = dscale * ((G2*qj-vPow2j)*C_ij(2)-uvj*C_ij(1))
    K_ij( 4) = dscale * (G1*qj-GAMMA*Ej)*aux1

    K_ij( 5) = dscale * C_ij(1)
    K_ij( 6) = dscale * (G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij( 7) = dscale * (vj*C_ij(1)-G1*uj*C_ij(2))
    K_ij( 8) = dscale * ((GAMMA*Ej-G2*qj)*C_ij(1)-G1*uj*aux1)

    K_ij( 9) = dscale * C_ij(2)
    K_ij(10) = dscale * (uj*C_ij(2)-G1*vj*C_ij(1))
    K_ij(11) = dscale * (uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(12) = dscale * ((GAMMA*Ej-G2*qj)*C_ij(2)-G1*vj*aux1)

    K_ij(13) = 0.0_DP
    K_ij(14) = dscale * G1*C_ij(1)
    K_ij(15) = dscale * G1*C_ij(2)
    K_ij(16) = dscale * (GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji( 1) = 0.0_DP
    K_ji( 2) = dscale * ((G1*qi-uPow2i)*C_ji(1)-uvi*C_ji(2))
    K_ji( 3) = dscale * ((G1*qi-vPow2i)*C_ji(2)-uvi*C_ji(1))
    K_ji( 4) = dscale * (G1*qi-GAMMA*Ei)*aux2

    K_ji( 5) = dscale * C_ji(1)
    K_ji( 6) = dscale * (G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji( 7) = dscale * (vi*C_ji(1)-G1*ui*C_ji(2))
    K_ji( 8) = dscale * ((GAMMA*Ei-G2*qi)*C_ji(1)-G1*ui*aux2)

    K_ji( 9) = dscale * C_ji(2)
    K_ji(10) = dscale * (ui*C_ji(2)-G1*vi*C_ji(1))
    K_ji(11) = dscale * (ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(12) = dscale * ((GAMMA*Ei-G2*qi)*C_ji(2)-G1*vi*aux2)

    K_ji(13) = 0.0_DP
    K_ji(14) = dscale * G1*C_ji(1)
    K_ji(15) = dscale * G1*C_ji(2)
    K_ji(16) = dscale * (GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

    ! Compute coefficients
    a = 0.5_DP*(C_ji-C_ij); anorm = sqrt(a(1)*a(1)+a(2)*a(2))

    if (anorm .gt. SYS_EPSREAL) then

      ! Compute Roe mean values
      aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
      u_ij = (aux*ui+uj)/(aux+1.0_DP)
      v_ij = (aux*vi+vj)/(aux+1.0_DP)
      hi   = GAMMA*U_i(4)/U_i(1)-G2*(ui*ui+vi*vi)
      hj   = GAMMA*U_j(4)/U_j(1)-G2*(uj*uj+vj*vj)
      H_ij = (aux*hi+hj)/(aux+1.0_DP)

      ! Compute auxiliary variables
      c1    = a(1)/anorm
      c2    = a(2)/anorm
      vel   = c1*u_ij+c2*v_ij
      q_ij  = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
      cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
      cs    = sqrt(cPow2)

      ! Diagonal matrix of eigenvalues
      l1 = abs(vel-cs)
      l2 = abs(vel)
      l3 = abs(vel+cs)
      l4 = abs(vel)

      ! Matrix of right eigenvectors
      R_ij(1,1) =  l1
      R_ij(2,1) =  l1*(u_ij-cs*c1)
      R_ij(3,1) =  l1*(v_ij-cs*c2)
      R_ij(4,1) =  l1*(H_ij-cs*vel)

      R_ij(1,2) =  l2
      R_ij(2,2) =  l2*u_ij
      R_ij(3,2) =  l2*v_ij
      R_ij(4,2) =  l2*q_ij

      R_ij(1,3) =  l3
      R_ij(2,3) =  l3*(u_ij+cs*c1)
      R_ij(3,3) =  l3*(v_ij+cs*c2)
      R_ij(4,3) =  l3*(H_ij+cs*vel)

      R_ij(1,4) =  0.0_DP
      R_ij(2,4) =  l4*c2
      R_ij(3,4) = -l4*c1
      R_ij(4,4) =  l4*(u_ij*c2-v_ij*c1)

      ! Matrix of left eigenvectors
      L_ij(1,1) = 0.5_DP*(G1*q_ij+cs*vel)/cPow2
      L_ij(2,1) = (cPow2-G1*q_ij)/cPow2
      L_ij(3,1) = 0.5_DP*(G1*q_ij-cs*vel)/cPow2
      L_ij(4,1) = v_ij*c1-u_ij*c2

      L_ij(1,2) = 0.5_DP*(-G1*u_ij-cs*c1)/cPow2
      L_ij(2,2) = G1*u_ij/cPow2
      L_ij(3,2) = 0.5_DP*(-G1*u_ij+cs*c1)/cPow2
      L_ij(4,2) = c2

      L_ij(1,3) = 0.5_DP*(-G1*v_ij-cs*c2)/cPow2
      L_ij(2,3) = G1*v_ij/cPow2
      L_ij(3,3) = 0.5_DP*(-G1*v_ij+cs*c2)/cPow2
      L_ij(4,3) = -c1

      L_ij(1,4) =  G2/cPow2
      L_ij(2,4) = -G1/cPow2
      L_ij(3,4) =  G2/cPow2
      L_ij(4,4) =  0.0_DP

      ! Compute tensorial dissipation D_ij = R_ij*|Lbd_ij|*L_ij
      call DGEMM('n', 'n', NVAR2D, NVAR2D, NVAR2D, anorm,&
                 R_ij, NVAR2D, L_ij, NVAR2D, 0.0_DP, D_ij, NVAR2D)

    else

      D_ij = 0.0_DP

    end if

  end subroutine euler_calcMatrixTensorDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixRusanovDiag2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP) :: ui,uj,vi,vj,ci,cj,Ei,Ej


    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1);   Ei = U_i(4)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1);   Ej = U_j(4)/U_j(1)

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale * (G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij(3) = dscale * (uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(4) = dscale * (GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale * (G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji(3) = dscale * (ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(4) = dscale * (GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

!!$    ! Compute auxiliary quantities
!!$    hi = GAMMA*Ei+(1-GAMMA)*0.5*(ui*ui+vi*vi)
!!$    hj = GAMMA*Ej+(1-GAMMA)*0.5*(uj*uj+vj*vj)
!!$
!!$    ci = sqrt(max((GAMMA-1)*(hi-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
!!$    cj = sqrt(max((GAMMA-1)*(hj-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))

    ! Compute the speed of sound
    ci = sqrt(max(G15*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
    cj = sqrt(max(G15*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))

    ! Compute dissipation tensor D_ij
    D_ij = dscale * max( abs(C_ij(1)*uj+C_ij(2)*vj) +&
                         sqrt(C_ij(1)**2+C_ij(2)**2)*cj,&
                         abs(C_ji(1)*ui+C_ji(2)*vi) +&
                         sqrt(C_ji(1)**2+C_ji(2)**2)*ci )

  end subroutine euler_calcMatrixRusanovDiag2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcMatrixRusanov2d(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U_i,U_j

    ! coefficients from spatial discretization
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
    real(DP) :: ci,cj,Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj
    real(DP) :: uPow2i,uPow2j,vPow2i,vPow2j,aux1,aux2

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1);   Ei = U_i(4)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1);   Ej = U_j(4)/U_j(1)

    uvi = ui*vi;   qi = ui*ui+vi*vi;   uPow2i = ui*ui;   vPow2i = vi*vi
    uvj = uj*vj;   qj = uj*uj+vj*vj;   uPow2j = uj*uj;   vPow2j = vj*vj

    aux1 = uj*C_ij(1)+vj*C_ij(2)
    aux2 = ui*C_ji(1)+vi*C_ji(2)

    ! Compute Galerkin coefficient K_ij
    K_ij( 1) = 0.0_DP
    K_ij( 2) = dscale * ((G2*qj-uPow2j)*C_ij(1)-uvj*C_ij(2))
    K_ij( 3) = dscale * ((G2*qj-vPow2j)*C_ij(2)-uvj*C_ij(1))
    K_ij( 4) = dscale * (G1*qj-GAMMA*Ej)*aux1

    K_ij( 5) = dscale * C_ij(1)
    K_ij( 6) = dscale * (G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij( 7) = dscale * (vj*C_ij(1)-G1*uj*C_ij(2))
    K_ij( 8) = dscale * ((GAMMA*Ej-G2*qj)*C_ij(1)-G1*uj*aux1)

    K_ij( 9) = dscale * C_ij(2)
    K_ij(10) = dscale * (uj*C_ij(2)-G1*vj*C_ij(1))
    K_ij(11) = dscale * (uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(12) = dscale * ((GAMMA*Ej-G2*qj)*C_ij(2)-G1*vj*aux1)

    K_ij(13) = 0.0_DP
    K_ij(14) = dscale * G1*C_ij(1)
    K_ij(15) = dscale * G1*C_ij(2)
    K_ij(16) = dscale * (GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji( 1) = 0.0_DP
    K_ji( 2) = dscale * ((G1*qi-uPow2i)*C_ji(1)-uvi*C_ji(2))
    K_ji( 3) = dscale * ((G1*qi-vPow2i)*C_ji(2)-uvi*C_ji(1))
    K_ji( 4) = dscale * (G1*qi-GAMMA*Ei)*aux2

    K_ji( 5) = dscale * C_ji(1)
    K_ji( 6) = dscale * (G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji( 7) = dscale * (vi*C_ji(1)-G1*ui*C_ji(2))
    K_ji( 8) = dscale * ((GAMMA*Ei-G2*qi)*C_ji(1)-G1*ui*aux2)

    K_ji( 9) = dscale * C_ji(2)
    K_ji(10) = dscale * (ui*C_ji(2)-G1*vi*C_ji(1))
    K_ji(11) = dscale * (ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(12) = dscale * ((GAMMA*Ei-G2*qi)*C_ji(2)-G1*vi*aux2)

    K_ji(13) = 0.0_DP
    K_ji(14) = dscale * G1*C_ji(1)
    K_ji(15) = dscale * G1*C_ji(2)
    K_ji(16) = dscale * (GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

    !---------------------------------------------------------------------------
    ! Evaluate the dissipation
    !---------------------------------------------------------------------------

!!$    ! Compute auxiliary quantities
!!$    hi = GAMMA*Ei+(1-GAMMA)*0.5*(ui*ui+vi*vi)
!!$    hj = GAMMA*Ej+(1-GAMMA)*0.5*(uj*uj+vj*vj)
!!$
!!$    ci = sqrt(max((GAMMA-1)*(hi-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
!!$    cj = sqrt(max((GAMMA-1)*(hj-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))

    ! Compute the speed of sound
    ci = sqrt(max(G15*(Ei-0.5_DP*(ui*ui+vi*vi)), SYS_EPSREAL))
    cj = sqrt(max(G15*(Ej-0.5_DP*(uj*uj+vj*vj)), SYS_EPSREAL))

    ! Compute dissipation tensor D_ij
    aux1 = dscale * max( abs(C_ij(1)*uj+C_ij(2)*vj) +&
                         sqrt(C_ij(1)**2+C_ij(2)**2)*cj,&
                         abs(C_ji(1)*ui+C_ji(2)*vi) +&
                         sqrt(C_ji(1)**2+C_ji(2)**2)*ci )

    D_ij = 0.0_DP
    D_ij( 1) = aux1
    D_ij( 6) = aux1
    D_ij(11) = aux1
    D_ij(16) = aux1

  end subroutine euler_calcMatrixRusanov2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_calcCharacteristics2d(&
      U_i, U_j, Dweight, W_ij, Lbd_ij, R_ij, L_ij)

!<description>
    ! This subroutine computes the characteristic variables in 2D
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
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: u_ij,v_ij,H_ij,q_ij,cs,aux,aux1,aux2,hi,hj
    real(DP) :: cPow2,uPow2,vPow2,a1,a2,anorm

    ! Compute norm of weighting coefficient
    anorm = sqrt(Dweight(1)*Dweight(1)+Dweight(2)*Dweight(2))

    ! Check if weighting coefficient is zero
    if (anorm .gt. SYS_EPSREAL) then

      ! Compute normalised weighting coefficient
      a1  = Dweight(1)/anorm
      a2  = Dweight(2)/anorm

      ! Compute Roe mean values
      aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
      u_ij = (aux*U_i(2)/U_i(1)+U_j(2)/U_j(1))/(aux+1.0_DP)
      v_ij = (aux*U_i(3)/U_i(1)+U_j(3)/U_j(1))/(aux+1.0_DP)
      hi   = GAMMA*U_i(4)/U_i(1)-G2*(U_i(2)*U_i(2)+U_i(3)*U_i(3))/(U_i(1)*U_i(1))
      hj   = GAMMA*U_j(4)/U_j(1)-G2*(U_j(2)*U_j(2)+U_j(3)*U_j(3))/(U_j(1)*U_j(1))
      H_ij = (aux*hi+hj)/(aux+1.0_DP)

      ! Compute auxiliary variables
      uPow2 = u_ij*u_ij
      vPow2 = v_ij*v_ij
      q_ij  = 0.5_DP*(uPow2+vPow2)
      cPow2 = max(G1*(H_ij-q_ij), SYS_EPSREAL)
      cs = sqrt(cPow2)
      aux   = a1*u_ij+a2*v_ij

      ! Compute diagonal matrix of eigenvalues (if present)
      if (present(Lbd_ij)) then
        Lbd_ij(1) = aux-cs
        Lbd_ij(2) = aux
        Lbd_ij(3) = aux+cs
        Lbd_ij(4) = aux
      end if

      ! Compute matrix of right eigenvectors
      if (present(R_ij)) then
        R_ij( 1) =  1.0_DP
        R_ij( 2) =  u_ij-cs*a1
        R_ij( 3) =  v_ij-cs*a2
        R_ij( 4) =  H_ij-cs*aux

        R_ij( 5) =  1.0_DP
        R_ij( 6) =  u_ij
        R_ij( 7) =  v_ij
        R_ij( 8) =  q_ij

        R_ij( 9) =  1.0_DP
        R_ij(10) =  u_ij+cs*a1
        R_ij(11) =  v_ij+cs*a2
        R_ij(12) =  H_ij+cs*aux

        R_ij(13) =  0.0_DP
        R_ij(14) =  a2
        R_ij(15) = -a1
        R_ij(16) =  u_ij*a2-v_ij*a1
      end if

      ! Compute matrix of left eigenvectors
      if (present(L_ij)) then
        L_ij( 1) =  0.5_DP*(G1*q_ij+cs*aux)/cPow2
        L_ij( 2) = (cPow2-G1*q_ij)/cPow2
        L_ij( 3) =  0.5_DP*(G1*q_ij-cs*aux)/cPow2
        L_ij( 4) =  v_ij*a1-u_ij*a2

        L_ij( 5) =  0.5_DP*(-G1*u_ij-cs*a1)/cPow2
        L_ij( 6) =  G1*u_ij/cPow2
        L_ij( 7) =  0.5_DP*(-G1*u_ij+cs*a1)/cPow2
        L_ij( 8) =  a2

        L_ij( 9) =  0.5_DP*(-G1*v_ij-cs*a2)/cPow2
        L_ij(10) =  G1*v_ij/cPow2
        L_ij(11) =  0.5_DP*(-G1*v_ij+cs*a2)/cPow2
        L_ij(12) = -a1

        L_ij(13) =  G2/cPow2
        L_ij(14) = -G1/cPow2
        L_ij(15) =  G2/cPow2
        L_ij(16) =  0.0_DP
      end if

      ! Compute characteristic solution difference
      if (present(W_ij)) then
        ! Compute solution difference U_j-U_i
        Diff = U_j-U_i

        ! Compute auxiliary quantities for characteristic variables
        aux1 = G2/cPow2*(q_ij*Diff(1)-&
                         u_ij*Diff(2)-&
                         v_ij*Diff(3)+&
                              Diff(4) )
        aux2 = 0.5_DP*(aux*Diff(1)-&
                        a1*Diff(2)-&
                        a2*Diff(3) )/cs

        ! Compute characteristic variables
        W_ij(1) = anorm * (aux1 + aux2)
        W_ij(2) = anorm * ((1.0_DP-G1*q_ij/cPow2)*Diff(1)+G1*(u_ij*Diff(2)+v_ij*Diff(3)-Diff(4))/cPow2)
        W_ij(3) = anorm * (aux1 - aux2)
        W_ij(4) = anorm * ((a1*v_ij-a2*u_ij)*Diff(1)+a2*Diff(2)-a1*Diff(3))
      end if

    else   ! |dweight| = 0

      if (present(Lbd_ij)) Lbd_ij = 0.0_DP
      if (present(R_ij))   R_ij   = 0.0_DP
      if (present(L_ij))   L_ij   = 0.0_DP
      if (present(W_ij))   W_ij   = 0.0_DP

    end if

  end subroutine euler_calcCharacteristics2d

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

    ! coefficients from spatial discretization
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

    ! coefficients from spatial discretization
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
      cs = sqrt(cPow2)

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

    ! coefficients from spatial discretization
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

  pure subroutine euler_trafoFluxDensity2d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density in 2D
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

  end subroutine euler_trafoFluxDensity2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffDensity2d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density in 2D
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

  end subroutine euler_trafoDiffDensity2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxEnergy2d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the energy in 2D
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
    G_ij(1) =  F_ij(4)
    G_ji(1) = -F_ij(4)

  end subroutine euler_trafoFluxEnergy2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffEnergy2d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the energy in 2D
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
    U_ij(1) =  U_j(4)-U_i(4)

  end subroutine euler_trafoDiffEnergy2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxPressure2d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the pressure in 2D
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
    real(DP) :: ui,uj,vi,vj

    ! velocities
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1)

    ! pressure fluxes
    G_ij(1) =  G1*(0.5_DP*(ui*ui+vi*vi)*F_ij(1)-ui*F_ij(2)-vi*F_ij(3)+F_ij(4))
    G_ji(1) = -G1*(0.5_DP*(uj*uj+vj*vj)*F_ij(1)-uj*F_ij(2)-vj*F_ij(3)+F_ij(4))

  end subroutine euler_trafoFluxPressure2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffPressure2d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the pressure in 2D
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
    pi = G1*(U_i(4)-0.5_DP*(U_i(2)*U_i(2)+U_i(3)*U_i(3))/U_i(1))
    pj = G1*(U_j(4)-0.5_DP*(U_j(2)*U_j(2)+U_j(3)*U_j(3))/U_j(1))

    ! pressure difference
    U_ij(1) = pj-pi

  end subroutine euler_trafoDiffPressure2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxVelocity2d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the y-velocity
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
    real(DP) :: ui,uj,vi,vj

    ! velocities
    ui = U_i(2)/U_i(1);   uj = U_j(2)/U_j(1)
    vi = U_i(3)/U_i(1);   vj = U_j(3)/U_j(1)

    ! velocity fluxes in x-direction
    G_ij(1) =  (F_ij(2)-ui*F_ij(1))/U_i(1)
    G_ji(1) = -(F_ij(2)-uj*F_ij(1))/U_j(1)

    ! velocity fluxes in y-direction
    G_ij(2) =  (F_ij(3)-vi*F_ij(1))/U_i(1)
    G_ji(2) = -(F_ij(3)-vj*F_ij(1))/U_j(1)
    
  end subroutine euler_trafoFluxVelocity2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffVelocity2d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the y-velocity
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

    ! velocity difference in y-direction
    U_ij(2) =  U_j(3)/U_j(1)-U_i(3)/U_i(1)
    
  end subroutine euler_trafoDiffVelocity2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxMomentum2d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the y-momentum
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

    ! momentum fluxes in y-direction
    G_ij(2) =  F_ij(3)
    G_ji(2) = -F_ij(3)
    
  end subroutine euler_trafoFluxMomentum2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffMomentum2d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the y-momentum
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

    ! momentum difference in y-direction
    U_ij(2) =  U_j(3)-U_i(3)
    
  end subroutine euler_trafoDiffMomentum2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxDenEng2d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 2D
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
    G_ij(2) =  F_ij(4)
    G_ji(2) = -F_ij(4)

  end subroutine euler_trafoFluxDenEng2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffDenEng2d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 2D
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
    U_ij(2) =  U_j(4)-U_i(4)

  end subroutine euler_trafoDiffDenEng2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxDenPre2d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 2D
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
    real(DP) :: ui,uj,vi,vj

    ! velocities
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1)

    ! density fluxes
    G_ij(1) =  F_ij(1)
    G_ji(1) = -F_ij(1)

    ! pressure fluxes
    G_ij(2) =  G1*(0.5_DP*(ui*ui+vi*vi)*F_ij(1)-ui*F_ij(2)-vi*F_ij(3)+F_ij(4))
    G_ji(2) = -G1*(0.5_DP*(uj*uj+vj*vj)*F_ij(1)-uj*F_ij(2)-vj*F_ij(3)+F_ij(4))

  end subroutine euler_trafoFluxDenPre2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffDenPre2d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 2D
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
    pi = G1*(U_i(4)-0.5_DP*(U_i(2)*U_i(2)+U_i(3)*U_i(3))/U_i(1))
    pj = G1*(U_j(4)-0.5_DP*(U_j(2)*U_j(2)+U_j(3)*U_j(3))/U_j(1))

    ! density difference
    U_ij(1) = U_j(1)-U_i(1)

    ! pressure difference
    U_ij(2) = pj-pi

  end subroutine euler_trafoDiffDenPre2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoFluxDenPreVel2d(U_i, U_j, F_ij, G_ij, G_ji)

!<description>
    ! This subroutine computes the transformation
    ! of the given flux into primitive variables in 2D
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
    real(DP) :: ui,uj,vi,vj

    ! velocities
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1)

    ! density fluxes
    G_ij(1) =  F_ij(1)
    G_ji(1) = -F_ij(1)

    ! velocity fluxes in x-direction
    G_ij(2) =  (F_ij(2)-ui*F_ij(1))/U_i(1)
    G_ji(2) = -(F_ij(2)-uj*F_ij(1))/U_j(1)

    ! velocity fluxes in y-direction
    G_ij(3) =  (F_ij(3)-vi*F_ij(1))/U_i(1)
    G_ji(3) = -(F_ij(3)-vj*F_ij(1))/U_j(1)

    ! pressure fluxes
    G_ij(4) =  G1*(0.5_DP*(ui*ui+vi*vi)*F_ij(1)-ui*F_ij(2)-vi*F_ij(3)+F_ij(4))
    G_ji(4) = -G1*(0.5_DP*(uj*uj+vj*vj)*F_ij(1)-uj*F_ij(2)-vj*F_ij(3)+F_ij(4))

  end subroutine euler_trafoFluxDenPreVel2d

  !*****************************************************************************

!<subroutine>

  pure subroutine euler_trafoDiffDenPreVel2d(U_i, U_j, U_ij)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density, pressure and velocity in 2D
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
    pi = G1*(U_i(4)-0.5_DP*(U_i(2)*U_i(2)+U_i(3)*U_i(3))/U_i(1))
    pj = G1*(U_j(4)-0.5_DP*(U_j(2)*U_j(2)+U_j(3)*U_j(3))/U_j(1))

    ! density difference
    U_ij(1) = U_j(1)-U_i(1)

    ! velocity difference in x-direction
    U_ij(2) =  U_j(2)/U_j(1)-U_i(2)/U_i(1)

    ! velocity difference in y-direction
    U_ij(3) =  U_j(3)/U_j(1)-U_i(3)/U_i(1)
    
    ! pressure difference
    U_ij(4) = pj-pi

  end subroutine euler_trafoDiffDenPreVel2d

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
    case(BDR_EULERWALL,&
         BDR_RLXEULERWALL)
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

      if (ibdrCondType .eq. BDR_EULERWALL) then
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


    case(BDR_VISCOUSWALL)
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


    case(BDR_SUPERINLET)
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


    case(BDR_FARFIELD)
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


    case(BDR_SUBINLET)
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
