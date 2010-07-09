!##############################################################################
!# ****************************************************************************
!# <name> eulerlagrange_callback2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible Euler/Navier-Stokes equations in 2D.
!#
!# The following callback functions are available:
!#
!# 1.) eulerlagrange_calcFluxGalerkin2d
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#
!# 2.) eulerlagrange_calcFluxGalerkinNoBdr2d
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#        without assembling the symmetric boundary contribution
!#
!# 3.) eulerlagrange_calcFluxScalarDiss2d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting scalar artificial viscosities
!#
!# 4.) eulerlagrange_calcFluxScalarDissDiSp2d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting scalar artificial viscosities based on
!#        dimensional splitting approach
!#
!# 5.) eulerlagrange_calcFluxTensorDiss2d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting tensorial artificial viscosities
!#
!# 6.) eulerlagrange_calcFluxTensorDissDiSp2d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting tensorial artificial viscosities based on
!#        dimensional splitting approach
!#
!# 7.) eulerlagrange_calcFluxRusanov2d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting the Rusanov artificial diffusion based on
!#        dimensional splitting approach
!#
!# 8.) eulerlagrange_calcFluxRusanovDiSp2d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting the Rusanov artificial diffusion
!#
!# 9.) eulerlagrange_calcMatrixDiagonalDiag2d
!#     -> Computes local matrix for diagonal entry
!#
!# 10.) eulerlagrange_calcMatrixDiagonal2d
!#      -> Computes local matrix for diagonal entry
!#
!# 11.) eulerlagrange_calcMatrixGalerkinDiag2d
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 12.) eulerlagrange_calcMatrixGalerkin2d
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 13.) eulerlagrange_calcMatrixScalarDissDiag2d
!#      -> Computes local matrices for low-order discretization
!#         adopting scalar artificial viscosities
!#
!# 14.) eulerlagrange_calcMatrixScalarDiss2d
!#      -> Computes local matrices for low-order discretization
!#         adopting scalar artificial viscosities
!#
!# 15.) eulerlagrange_calcMatrixTensorDissDiag2d
!#      -> Computes local matrices for low-order discretization
!#         adopting tensorial artificial viscosities
!#
!# 16.) eulerlagrange_calcMatrixTensorDiss2d
!#      -> Computes local matrices for low-order discretization
!#         adopting tensorial artificial viscosities
!#
!# 17.) eulerlagrange_calcMatrixRusanovDiag2d
!#      -> Computes local matrices for low-order discretization
!#         adopting the Rusanov artificial viscosities
!#
!# 18.) eulerlagrange_calcMatrixRusanov2d
!#      -> Computes local matrices for low-order discretization
!#         adopting the Rusanov flux artificial viscosities
!#
!# 19.) eulerlagrange_calcCharacteristics2d
!#      -> Computes characteristic variables
!#
!# 20.) eulerlagrange_calcFluxFCTScalarDiss2d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting scalar artificial viscosities
!#
!# 21.) eulerlagrange_calcFluxFCTTensorDiss2d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting tensorial artificial viscosities
!#
!# 22.) eulerlagrange_calcFluxFCTRusanov2d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting the Rusanov artificial viscosities
!#
!# 23.) eulerlagrange_trafoFluxDensity2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 24.) eulerlagrange_trafoDiffDensity2d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 25.) eulerlagrange_trafoFluxEnergy2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 26.) eulerlagrange_trafoDiffEnergy2d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 27.) eulerlagrange_trafoFluxPressure2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 28.) eulerlagrange_trafoFluxVelocity2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 29.) eulerlagrange_trafoDiffVelocity2d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 30.) eulerlagrange_trafoFluxMomentum2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 31.) eulerlagrange_trafoDiffMomentum2d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 32.) eulerlagrange_trafoFluxDenEng2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 33.) eulerlagrange_trafoDiffDenEng2d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 34.) eulerlagrange_trafoFluxDenPre2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 35.) eulerlagrange_trafoDiffDenPre2d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 36.) eulerlagrange_trafoFluxDenPreVel2d
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 37.) eulerlagrange_trafoDiffDenPreVel2d
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure 
!#         and the velocity
!#
!# 38.) eulerlagrange_calcBoundaryvalues2d
!#      -> Computes the boundary values for a given node
!#
!# 39.) eulerlagrange_hadaptCallbackScalar2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in interleave format
!#
!# 40.) eulerlagrange_hadaptCallbackBlock2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in block format
!#
!# 41.) eulerlagrange_calcbarycoords
!#      -> Calculates the barycentric coordinates of the particle position
!#
!# 42.) eulerlagrange_findelement
!#      -> Finds the element in wich the particle is
!#
!# 43.) eulerlagrange_wrongelement
!#      -> Finds the right element with the findalgorithm finds the wrong element
!#
!# 44.) eulerlagrange_moveparticlesoneway
!#      -> Calculates the new position of the particles
!#
!# 45.) eulerlagrange_moveparticlestwoway
!#      -> Calculates the new position of the particles with two way coppling
!#
!# 46.) eulerlagrange_calcvolpart
!#      -> Calculates the volume fraction of the particles in the gridpoints
!#
!# 47.) eulerlagrange_calcvelopart
!#      -> Calculates the velocity of the particles in the gridpoints
!#
!# </purpose>
!##############################################################################

module eulerlagrange_callback2d

  use basicgeometry
  use boundaryfilter
  use collection
  use eulerlagrange_basic
  use flagship_callback
  use fsystem
  use genoutput
  use geometryaux
  use graph
  use groupfemsystem
  use hadaptaux
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use problem
  use pprocsolution
  use solveraux
  use storage
  use thermodynamics
  use triasearch

  implicit none

  private
  public :: eulerlagrange_calcFluxGalerkin2d
  public :: eulerlagrange_calcFluxGalerkinNoBdr2d
  public :: eulerlagrange_calcFluxScalarDiss2d
  public :: eulerlagrange_calcFluxScalarDissDiSp2d
  public :: eulerlagrange_calcFluxTensorDiss2d
  public :: eulerlagrange_calcFluxTensorDissDiSp2d
  public :: eulerlagrange_calcFluxRusanov2d
  public :: eulerlagrange_calcFluxRusanovDiSp2d
  public :: eulerlagrange_calcMatrixDiagonalDiag2d
  public :: eulerlagrange_calcMatrixDiagonal2d
  public :: eulerlagrange_calcMatrixGalerkinDiag2d
  public :: eulerlagrange_calcMatrixGalerkin2d
  public :: eulerlagrange_calcMatrixScalarDissDiag2d
  public :: eulerlagrange_calcMatrixScalarDiss2d
  public :: eulerlagrange_calcMatrixTensorDissDiag2d
  public :: eulerlagrange_calcMatrixTensorDiss2d
  public :: eulerlagrange_calcMatrixRusanovDiag2d
  public :: eulerlagrange_calcMatrixRusanov2d
  public :: eulerlagrange_calcCharacteristics2d
  public :: eulerlagrange_calcFluxFCTScalarDiss2d
  public :: eulerlagrange_calcFluxFCTTensorDiss2d
  public :: eulerlagrange_calcFluxFCTRusanov2d
  public :: eulerlagrange_trafoFluxDensity2d
  public :: eulerlagrange_trafoFluxEnergy2d
  public :: eulerlagrange_trafoFluxPressure2d
  public :: eulerlagrange_trafoFluxVelocity2d
  public :: eulerlagrange_trafoFluxMomentum2d
  public :: eulerlagrange_trafoFluxDenEng2d
  public :: eulerlagrange_trafoFluxDenPre2d
  public :: eulerlagrange_trafoFluxDenPreVel2d
  public :: eulerlagrange_trafoDiffDensity2d
  public :: eulerlagrange_trafoDiffEnergy2d
  public :: eulerlagrange_trafoDiffPressure2d
  public :: eulerlagrange_trafoDiffVelocity2d
  public :: eulerlagrange_trafoDiffMomentum2d
  public :: eulerlagrange_trafoDiffDenEng2d
  public :: eulerlagrange_trafoDiffDenPre2d
  public :: eulerlagrange_trafoDiffDenPreVel2d
  public :: eulerlagrange_calcBoundaryvalues2d
  public :: eulerlagrange_hadaptCallbackScalar2d
  public :: eulerlagrange_hadaptCallbackBlock2d
  public :: eulerlagrange_calcbarycoords
  public :: eulerlagrange_findelement
  public :: eulerlagrange_wrongelement
  public :: eulerlagrange_moveparticlestwoway
  public :: eulerlagrange_moveparticlesoneway
  public :: eulerlagrange_getbarycoordelm
  public :: eulerlagrange_calcvolpart
  public :: eulerlagrange_calcvelopart

  type , public :: t_Particles
      !number of particles
      integer :: nPart 
      ! element
      integer(I32) :: h_element
      integer, dimension(:), pointer :: p_element
      ! position
      integer(I32) :: h_xpos, h_ypos, h_zpos
      real(DP), dimension(:), pointer :: p_xpos, p_ypos, p_zpos
      ! old position
      integer(I32) :: h_xpos_old, h_ypos_old, h_zpos_old
      real(DP), dimension(:), pointer :: p_xpos_old, p_ypos_old, p_zpos_old
      ! velocity
      integer(I32) :: h_xvelo, h_yvelo, h_zvelo
      real(DP), dimension(:), pointer :: p_xvelo, p_yvelo, p_zvelo
      ! old velocity
      integer(I32) :: h_xvelo_old, h_yvelo_old, h_zvelo_old
      real(DP), dimension(:), pointer :: p_xvelo_old, p_yvelo_old, p_zvelo_old
      ! velocity of the gas
      integer(I32) :: h_xvelo_gas, h_yvelo_gas, h_zvelo_gas
      real(DP), dimension(:), pointer :: p_xvelo_gas, p_yvelo_gas, p_zvelo_gas
      ! old velocity of the gas
      integer(I32) :: h_xvelo_gas_old, h_yvelo_gas_old, h_zvelo_gas_old
      real(DP), dimension(:), pointer :: p_xvelo_gas_old, p_yvelo_gas_old, p_zvelo_gas_old
      ! barycentric coordinates
      integer(I32) :: h_lambda1, h_lambda2, h_lambda3, h_lambda4
      real(DP), dimension(:), pointer :: p_lambda1, p_lambda2, p_lambda3, p_lambda4
      ! diameter, mass and alpha_n 
      integer(I32) :: h_diam, h_mass, h_alpha_n, h_temp, h_density
      real(DP), dimension(:), pointer :: p_diam, p_mass, p_alpha_n, p_temp, p_density
      ! midpoints of the element 
      integer(I32) :: h_midpoints_el
      real(DP), dimension(:,:), pointer :: p_midpoints_el
      ! volume fraction of the particles 
      integer(I32) :: h_PartVol, h_PartVolAver
      real(DP), dimension(:), pointer :: p_PartVol, p_PartVolAver
      integer :: iPartVolCount
      ! volumepart of the particles 
      integer(I32) :: h_PartVelox, h_PartVeloy
      real(DP), dimension(:), pointer :: p_PartVelox, p_PartVeloy
      ! gravity
      real(DP), dimension(2)  :: gravity
      ! viscosity of the gas
      real(DP) :: nu_g
      ! parameter for particle-wall collisions
      real(DP)  :: tang_val, norm_val
      ! variable for particle-wall collisions
      integer(I32) :: h_bdy_time
      real(DP), dimension(:), pointer :: p_bdy_time
      ! timestep for video
      integer :: iTimestep
      ! maximum value in x-direction
      real(DP):: maxvalx
      
   end type t_Particles

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxGalerkin2d(&
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
#ifdef USE_EULERLAGRANGE_IBP
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

#ifdef USE_EULERLAGRANGE_IBP
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

  end subroutine eulerlagrange_calcFluxGalerkin2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxGalerkinNoBdr2d(&
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
    ! quantities have a look at the subroutine "eulerlagrange_calcFluxGalerkin2d".
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

  end subroutine eulerlagrange_calcFluxGalerkinNoBdr2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxScalarDiss2d(&
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
#ifdef USE_EULERLAGRANGE_IBP
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
    ! quantities have a look at the subroutine "eulerlagrange_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------

    ! Compute velocities
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); rv2i = vi*U_i(3)
    ru2j = uj*U_j(2); rv2j = vj*U_j(3)

#ifdef USE_EULERLAGRANGE_IBP
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

  end subroutine eulerlagrange_calcFluxScalarDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxScalarDissDiSp2d(&
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
#ifdef USE_EULERLAGRANGE_IBP
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
    ! quantities have a look at the subroutine "eulerlagrange_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------

    ! Compute velocities
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); rv2i = vi*U_i(3)
    ru2j = uj*U_j(2); rv2j = vj*U_j(3)

#ifdef USE_EULERLAGRANGE_IBP
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

  end subroutine eulerlagrange_calcFluxScalarDissDiSp2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxTensorDiss2d(&
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
#ifdef USE_EULERLAGRANGE_IBP
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
    ! quantities have a look at the subroutine "eulerlagrange_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------

    ! Compute velocities
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); rv2i = vi*U_i(3)
    ru2j = uj*U_j(2); rv2j = vj*U_j(3)

#ifdef USE_EULERLAGRANGE_IBP
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

  end subroutine eulerlagrange_calcFluxTensorDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxTensorDissDiSp2d(&
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
#ifdef USE_EULERLAGRANGE_IBP
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
    ! quantities have a look at the subroutine "eulerlagrange_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------

    ! Compute velocities
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); rv2i = vi*U_i(3)
    ru2j = uj*U_j(2); rv2j = vj*U_j(3)

#ifdef USE_EULERLAGRANGE_IBP
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

  end subroutine eulerlagrange_calcFluxTensorDissDiSp2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxRusanov2d(&
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
#ifdef USE_EULERLAGRANGE_IBP
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
    ! quantities have a look at the subroutine "eulerlagrange_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------

    ! Compute velocities and energy
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1); Ei = U_i(4)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1); Ej = U_j(4)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); rv2i = vi*U_i(3)
    ru2j = uj*U_j(2); rv2j = vj*U_j(3)

#ifdef USE_EULERLAGRANGE_IBP
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

  end subroutine eulerlagrange_calcFluxRusanov2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxRusanovDiSp2d(&
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
#ifdef USE_EULERLAGRANGE_IBP
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
    ! quantities have a look at the subroutine "eulerlagrange_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------

    ! Compute velocities and energy
    ui = U_i(2)/U_i(1); vi = U_i(3)/U_i(1); Ei = U_i(4)/U_i(1)
    uj = U_j(2)/U_j(1); vj = U_j(3)/U_j(1); Ej = U_j(4)/U_j(1)

    ! Compute auxiliary variables
    ru2i = ui*U_i(2); rv2i = vi*U_i(3)
    ru2j = uj*U_j(2); rv2j = vj*U_j(3)

#ifdef USE_EULERLAGRANGE_IBP
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

  end subroutine eulerlagrange_calcFluxRusanovDiSp2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatrixDiagonalDiag2d(U_i, C_ii, i, dscale, K_ii)

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

  end subroutine eulerlagrange_calcMatrixDiagonalDiag2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatrixDiagonal2d(U_i, C_ii, i, dscale, K_ii)

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

  end subroutine eulerlagrange_calcMatrixDiagonal2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatrixGalerkinDiag2d(&
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

  end subroutine eulerlagrange_calcMatrixGalerkinDiag2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatrixGalerkin2d(&
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

  end subroutine eulerlagrange_calcMatrixGalerkin2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatrixScalarDissDiag2d(&
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

  end subroutine eulerlagrange_calcMatrixScalarDissDiag2d

!*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatrixScalarDiss2d(&
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

  end subroutine eulerlagrange_calcMatrixScalarDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatrixTensorDissDiag2d(&
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

  end subroutine eulerlagrange_calcMatrixTensorDissDiag2d

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_calcMatrixTensorDiss2d(&
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

  end subroutine eulerlagrange_calcMatrixTensorDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatrixRusanovDiag2d(&
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

  end subroutine eulerlagrange_calcMatrixRusanovDiag2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatrixRusanov2d(&
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

  end subroutine eulerlagrange_calcMatrixRusanov2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcCharacteristics2d(&
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

      ! Compute normalized weighting coefficient
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

  end subroutine eulerlagrange_calcCharacteristics2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxFCTScalarDiss2d(&
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

  end subroutine eulerlagrange_calcFluxFCTScalarDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxFCTTensorDiss2d(&
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
  end subroutine eulerlagrange_calcFluxFCTTensorDiss2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxFCTRusanov2d(&
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

  end subroutine eulerlagrange_calcFluxFCTRusanov2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxDensity2d(U_i, U_j, F_ij, G_ij, G_ji)

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

  end subroutine eulerlagrange_trafoFluxDensity2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffDensity2d(U_i, U_j, U_ij)

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

  end subroutine eulerlagrange_trafoDiffDensity2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxEnergy2d(U_i, U_j, F_ij, G_ij, G_ji)

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

  end subroutine eulerlagrange_trafoFluxEnergy2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffEnergy2d(U_i, U_j, U_ij)

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

  end subroutine eulerlagrange_trafoDiffEnergy2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxPressure2d(U_i, U_j, F_ij, G_ij, G_ji)

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

  end subroutine eulerlagrange_trafoFluxPressure2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffPressure2d(U_i, U_j, U_ij)

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

  end subroutine eulerlagrange_trafoDiffPressure2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxVelocity2d(U_i, U_j, F_ij, G_ij, G_ji)

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
    
  end subroutine eulerlagrange_trafoFluxVelocity2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffVelocity2d(U_i, U_j, U_ij)

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
    
  end subroutine eulerlagrange_trafoDiffVelocity2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxMomentum2d(U_i, U_j, F_ij, G_ij, G_ji)

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
    
  end subroutine eulerlagrange_trafoFluxMomentum2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffMomentum2d(U_i, U_j, U_ij)

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
    
  end subroutine eulerlagrange_trafoDiffMomentum2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxDenEng2d(U_i, U_j, F_ij, G_ij, G_ji)

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

  end subroutine eulerlagrange_trafoFluxDenEng2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffDenEng2d(U_i, U_j, U_ij)

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

  end subroutine eulerlagrange_trafoDiffDenEng2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxDenPre2d(U_i, U_j, F_ij, G_ij, G_ji)

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

  end subroutine eulerlagrange_trafoFluxDenPre2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffDenPre2d(U_i, U_j, U_ij)

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

  end subroutine eulerlagrange_trafoDiffDenPre2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxDenPreVel2d(U_i, U_j, F_ij, G_ij, G_ji)

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

  end subroutine eulerlagrange_trafoFluxDenPreVel2d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffDenPreVel2d(U_i, U_j, U_ij)

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

  end subroutine eulerlagrange_trafoDiffDenPreVel2d

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_calcBoundaryvalues2d(DbdrNormal, DpointNormal,&
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
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcBoundaryvalues2d')
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
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcBoundaryvalues2d')
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
          OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcBoundaryvalues2d')
      call sys_halt()
    end select

  end subroutine eulerlagrange_calcBoundaryvalues2d

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_hadaptCallbackScalar2d(iOperation, rcollection)

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
            OU_CLASS_WARNING,OU_MODE_STD,'eulerlagrange_hadaptCallbackScalar2d')
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

  end subroutine eulerlagrange_hadaptCallbackScalar2d

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_hadaptCallbackBlock2d(iOperation, rcollection)

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
            OU_CLASS_WARNING,OU_MODE_STD,'eulerlagrange_hadaptCallbackBlock2d')
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

  end subroutine eulerlagrange_hadaptCallbackBlock2d

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_calcbarycoords(p_rproblemLevel,rParticles,iPart)

!<description>
    ! This subroutine computes the barycentric coordinates of the position of the particle in the element.

!<input>
    ! Particles
    type(t_Particles), intent(inout) :: rParticles
    
    ! Current number of particle
    integer, intent(inout) :: iPart

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! Pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! Pointer to vertices at each element
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Pointer to the vertex coordinates
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! Current element
    integer :: currentElement

    ! Coordinates of the vertices of the actual element
    real(DP), dimension(2,3) :: vert_coord

    ! Determinates
    real(DP) :: det_A, det_A1, det_A2, det_A3

    ! Local variables
    integer :: ivt

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation
  
    ! Get vertices at element
    call storage_getbase_int2D(&
        p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    ! Get coordinates of the vertices
    call storage_getbase_double2D(&
        p_rtriangulation%h_DvertexCoords, p_DvertexCoords)

    ! store current element
    currentElement = rParticles%p_element(iPart)

    ! store coordinates of the vertices
    do ivt=1,3

      vert_coord(1,ivt)= p_DvertexCoords(1,p_IverticesAtElement(ivt,currentElement))
      vert_coord(2,ivt)= p_DvertexCoords(2,p_IverticesAtElement(ivt,currentElement))

    end do

    ! compute determinate
    !	    1 x1 y1
    !    A=	1 x2 y2
    !		1 x3 y3
    ! detA=x2y3+x1y2+x3y1-x2y1-x3y2-x1y3
    det_A=  vert_coord(1,2) * vert_coord(2,3)+&
            vert_coord(1,1) * vert_coord(2,2)+&
            vert_coord(1,3) * vert_coord(2,1)-&
            vert_coord(1,2) * vert_coord(2,1)-&
            vert_coord(1,3) * vert_coord(2,2)-&
            vert_coord(1,1) * vert_coord(2,3)

    ! calculate barycentric coorinates (lambda1,lambda2,lambda3)
    ! lambda1
    !		1 x  y
    !   A1=	1 x2 y2
    !		1 x3 y3
    ! detA1=x2y3+xy2+x3y-x2y-x3y2-xy3
    det_A1= vert_coord(1,2)         * vert_coord(2,3)+&
            rParticles%p_xpos(iPart)* vert_coord(2,2)+&
            vert_coord(1,3)         * rParticles%p_ypos(iPart)-&
            vert_coord(1,2)         * rParticles%p_ypos(iPart)-&
            vert_coord(1,3)         * vert_coord(2,2)-&
            rParticles%p_xpos(iPart)* vert_coord(2,3)
 
    !lambda1=|det_A1/det_A|
    rParticles%p_lambda1(iPart)= abs(det_A1/det_A)

    ! lambda2
    !		1 x1 y1
    !   A2=	1 x  y
    !		1 x3 y3
    ! detA2=xy3+x1y+x3y1-xy1-x3y-x1y3
    det_A2= rParticles%p_xpos(iPart)* vert_coord(2,3)+&
            vert_coord(1,1)         * rParticles%p_ypos(iPart)+&
            vert_coord(1,3)         * vert_coord(2,1)-&
            rParticles%p_xpos(iPart)* vert_coord(2,1)-&
            vert_coord(1,3)         * rParticles%p_ypos(iPart)-&
            vert_coord(1,1)         * vert_coord(2,3)

    !lambda2=|det_A2/det_A|
    rParticles%p_lambda2(iPart)= abs(det_A2/det_A)

    ! lambda3
    !		1 x1 y1
    !   A3=	1 x2 y2
    !		1 x  y
    ! detA3=x2y+x1y2+xy1-x2y1-xy2-x1y
    det_A3= vert_coord(1,2)         * rParticles%p_ypos(iPart)+&
            vert_coord(1,1)         * vert_coord(2,2)+&
            rParticles%p_xpos(iPart)* vert_coord(2,1)-&
            vert_coord(1,2)         * vert_coord(2,1)-&
            rParticles%p_xpos(iPart)* vert_coord(2,2)-&
            vert_coord(1,1)         * rParticles%p_ypos(iPart)

    ! lambda3=|det_A3/det_A|
    rParticles%p_lambda3(iPart)= abs(det_A3/det_A)


  end subroutine eulerlagrange_calcbarycoords

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_findelement(rparlist,p_rproblemLevel,rParticles,iPart)

!<description>
    ! This subroutine searchs in which element the particle is.

!<input>
    ! parameterlist
    type(t_parlist), intent(inout) :: rparlist
    
    ! particles
    type(t_Particles), intent(inout) :: rParticles

    ! pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation


    ! pointer to the neighbour elements adjacent to an element
    !
    ! Handle to 
    !       p_IneighboursAtElement = array [1..TRIA_MAXNME2D,1..NEL] of integer
    ! For each element, the numbers of adjacent elements
    ! in mathematically positive sense, meeting the element in an edge.
    ! p_RneighbourElement(IEL)\%Ineighbours(.) describes the elements adjacent 
    ! to IEL along the edges (p_RedgesOnElement(IEL)\%Iedges(.)-NVT).
    ! This is the old KADJ array.
    !
    ! Note:  For meshes with hanging vertices, this array is slightly
    ! modified. For 'big' elements this array contains the element
    ! numbers of the 'first' adjacent element via an edge/face.
    ! Note: To access all elements adjacent to an element via a
    ! hanging vertex, calculate the vertex number of the hanging
    ! vertex via InodalProperty and access all adjacent elements.
    ! For small 'hanging' elements, this array contains as usual
    ! the number of the 'big' adjacent element(s).
    integer, dimension(:,:), pointer :: p_IneighboursAtElement


    ! pointer to the coordinates of the vertices
    !
    ! A list of all corner(!)-vertices of the elements in the triangulation.
    ! Handle to 
    !       p_RcornerCoordinates = array [1..ndim,1..NVT] of double
    ! with
    !   p_DvertexCoords(1,.) = X-coordinate.
    ! for 1D meshes,
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    ! for 2D meshes and
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    !   p_DvertexCoords(3,.) = Z-coordinate.
    ! for 3D meshes.
    ! This is a handle to the old DCORVG-array.
    !
    ! Note that the array may be longer than NVT in general!
    ! (May happen in case of a mesh hierarchy generated by a 2-level
    ! refinement, where the coordinates of the points on the
    ! coarser levels are contained in te coordinates of the
    ! finer levels.)
    ! In such a case, only the first NVT n-tuples in this array are valid!
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement

   
    ! Current number of particle
    integer, intent(inout) :: iPart
  
    ! Position of the element
    real(DP), dimension(2) :: particlepos

    ! Element number
    integer :: iel
 
    ! Particle is in element 
    logical :: binside
    
    ! Search mode
    character(LEN=15) :: searchmode

    ! Variables for midpoints_el
	real(DP), dimension(1:4) :: distances
	integer :: i, adj, minl, ite
	integer, parameter :: itemax = 100000
	real(DP) :: distToMid

    ! Variables for raytrace2D
	integer :: iresult
	    ! =1 : The element was found successfully.
        ! =0 : The raytracing search broke down inside of the domain. 
        ! =-1: The search broke down because the domain was left.

    integer :: ilastElement
        !Last analysed element.
        ! If iresult= 1: ilastElement = iel
        ! If iresult= 0: Number of the last analysed element before the search
        !                was stopped.
        ! If iresult=-1: Number of the element through which the
        !                domain was left. 
    integer :: ilastEdge
        ! Number of the last analysed edge. Range 1..NMT.
        ! If iresult= 1: ilastEdge=0
        ! If iresult= 0: Number of the last analysed edge before the search
        !                was stopped.
        ! If iresult=-1: Number of the edge through which the domain was left. 


    real(DP), dimension(2,4) :: DcornerCoords


    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation
   
    ! Get vertices at element
    call storage_getbase_int2D(&
        p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    ! Get coordinates for elements
    call storage_getbase_double2d (&
        p_rtriangulation%h_DvertexCoords,p_DvertexCoords)
        
    ! Get neighboured elements
    call storage_getbase_int2D(&
        p_rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)
    
    ! Set particles position
    particlepos(1)= rParticles%p_xpos(iPart)
    particlepos(2)= rParticles%p_ypos(iPart) 
   
    ! Get searchmode (brute force, raytrace, etc.)  
    call parlst_getvalue_string(rparlist, 'Eulerlagrange', "search", searchmode)

    select case(trim(searchmode))
    case('bruteforce')
        call tsrch_getElem_BruteForce(particlepos,p_DvertexCoords,p_IverticesAtElement,iel)
    case('raytrace2D')
        call tsrch_getElem_raytrace2D(&
                particlepos,p_rtriangulation,rParticles%p_element(iPart),iresult,ilastElement,ilastEdge,itemax)
    case('midpoint')

	    distToMid = 10000.0_dp

	    gotoNextElm: do ite = 1, itemax
	    
		    distances = 10000.0_dp
   
            if (rParticles%p_element(iPart)==0) rParticles%p_element(iPart)=1
   
		    ! Calculate the distances to the midpoints
		    do i = 1, 3 
		      if (p_IneighboursAtElement(i,rParticles%p_element(iPart)) > 0) then
			    adj = p_IneighboursAtElement(i,rParticles%p_element(iPart))
			    distances(i+1) = (rParticles%p_xpos(iPart)-rParticles%p_midpoints_el(1,adj))**2.0_dp +&
								 (rParticles%p_ypos(iPart)-rParticles%p_midpoints_el(2,adj))**2.0_dp
			  end if
		    end do
		    
		    ! Distance of the current element to the new position
		    distances(1) = (rParticles%p_xpos(iPart)-&
		                        rParticles%p_midpoints_el(1,rParticles%p_element(iPart)))**2.0_dp+&
		                   (rParticles%p_ypos(iPart)-&
		                        rParticles%p_midpoints_el(2,rParticles%p_element(iPart)))**2.0_dp
		
		    ! Position with the smallest distance
		    minl = minloc(distances,1)

		    ! Check if the distance didn't change
		    if (minl .eq. 1) exit gotoNextElm

		    ! Store element with the lowest distance
		    rParticles%p_element(iPart) = p_IneighboursAtElement(minl-1,rParticles%p_element(iPart))
		    distToMid = distances(minl)

	    end do gotoNextElm
	    
    case default
        call output_line('Invalid search mode!',&
                     OU_CLASS_WARNING,OU_MODE_STD,'flagship')
        call sys_halt()
    end select

  end subroutine eulerlagrange_findelement


  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_wrongelement(rparlist,p_rproblemLevel,rParticles,iPart)

!<description>
    ! If the search algorithim finds the wrong element, the subroutine searchh in the neighboured elements.
    
!<input>
    ! Parameterlist
    type(t_parlist), intent(inout) :: rparlist

    ! Particles
    type(t_Particles), intent(inout) :: rParticles

    ! Current number of particle
    integer, intent(inout) :: iPart

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! pointer to array containing the elements adjacent to a vertex.
    !
    ! Handle to 
    !       p_IelementsAtVertex = array(1..*) of integer
    ! p_IelementsAtVertex ( p_IelementsAtVertexIdx(IVT)..p_IelementsAtVertexIdx(IVT+1)-1 )
    ! contains the number of the adjacent element in a vertex.
    ! This replaces the old KVEL array.
    !
    ! Note: For hanging vertices, this array contains only those
    ! elements which are 'corner adjacent' to a vertex (i.e. the 'smaller' elements).
    ! The 'big' elements adjacent to the edge which the hanging vertex
    ! is a midpoint of are not part of the vertex neighbourhood
    ! in this array.
    integer, dimension(:), pointer :: p_IelementsAtVertex

    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Handle to 
    !       p_IelementsAtVertexIdx=array [1..NVT+1] of integer.
    ! Index array for p_IelementsAtVertex of length NVT+1 for describing the
    ! elements adjacent to a corner vertex. for vertex IVT, the array
    ! p_IelementsAtVertex contains the numbers of the elements around this
    ! vertex at indices 
    !     p_IelementsAtVertexIdx(IVT)..p_IelementsAtVertexIdx(IVT+1)-1.
    ! By subtracting
    !     p_IelementsAtVertexIdx(IVT+1)-p_IelementsAtVertexIdx(IVT)
    ! One can get the number of elements adjacent to a vertex IVT.
    integer, dimension(:), pointer ::  p_IelementsAtVertexIdx

    ! pointer to the coordinates of the vertices
    !
    ! A list of all corner(!)-vertices of the elements in the triangulation.
    ! Handle to 
    !       p_RcornerCoordinates = array [1..ndim,1..NVT] of double
    ! with
    !   p_DvertexCoords(1,.) = X-coordinate.
    ! for 1D meshes,
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    ! for 2D meshes and
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    !   p_DvertexCoords(3,.) = Z-coordinate.
    ! for 3D meshes.
    ! This is a handle to the old DCORVG-array.
    !
    ! Note that the array may be longer than NVT in general!
    ! (May happen in case of a mesh hierarchy generated by a 2-level
    ! refinement, where the coordinates of the points on the
    ! coarser levels are contained in te coordinates of the
    ! finer levels.)
    ! In such a case, only the first NVT n-tuples in this array are valid!
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! Pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

 	! Local variables
 	integer :: Vert, Elm, currentelm, nVertex, Neighbour, i_NVBD, i
    real(DP) :: dxi1,dxi2,dxi3
    real(DP) :: dx,dy
    ! Coordinates of the cornervertices
    real(DP), dimension(2,4) :: DcornerCoords

    logical :: binside

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation
  
    call storage_getbase_int2D(&
         p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
         
    call storage_getbase_int(p_rtriangulation%h_IelementsAtVertex,&
        p_IelementsAtVertex)
        
    call storage_getbase_int(p_rtriangulation%h_IelementsAtVertexIdx,&
        p_IelementsAtVertexIdx)

    call storage_getbase_double2d (&
        p_rtriangulation%h_DvertexCoords,p_DvertexCoords)


    ! Store the current element
	currentelm = rParticles%p_element(iPart)
	dxi1=0.0_dp
	dxi2=0.0_dp
	dxi3=0.0_dp
	dx= rParticles%p_xpos(iPart)
	dy= rParticles%p_ypos(iPart)
	
    ! Loop over the vertices of the element
	SearchVertex: do Vert = 1, 3												
		
		!Current vertex
		nVertex = p_IverticesAtElement(Vert, currentelm)
			
		! Loop over the element containing to the vertex
		SearchElement: do Elm = 1, (p_IelementsAtVertexIdx(nVertex+1)-p_IelementsAtVertexIdx(nVertex))		
											
	    	if (p_IelementsAtVertex(p_IelementsAtVertexIdx(nVertex)+Elm-1) == 0) then
				exit SearchElement
			end if

			rParticles%p_element(iPart) = p_IelementsAtVertex(p_IelementsAtVertexIdx(nVertex)+Elm-1)

            ! Store coordinates of cornervertices
            DcornerCoords(1,1)= p_DvertexCoords(1,p_IverticesAtElement(1,rParticles%p_element(iPart)))
            DcornerCoords(1,2)= p_DvertexCoords(1,p_IverticesAtElement(2,rParticles%p_element(iPart)))
            DcornerCoords(1,3)= p_DvertexCoords(1,p_IverticesAtElement(3,rParticles%p_element(iPart)))
            DcornerCoords(2,1)= p_DvertexCoords(2,p_IverticesAtElement(1,rParticles%p_element(iPart)))
            DcornerCoords(2,2)= p_DvertexCoords(2,p_IverticesAtElement(2,rParticles%p_element(iPart)))
            DcornerCoords(2,3)= p_DvertexCoords(2,p_IverticesAtElement(3,rParticles%p_element(iPart)))

            ! Check if the particle is in the element
            call gaux_isInElement_tri2D(dx,dy,DcornerCoords,binside)
          
            ! If the particle is in the element, then exit loop
            if (binside) then
                exit SearchVertex
            end if 
            
		end do SearchElement

	end do SearchVertex

    ! Store coordinates of cornervertices
    DcornerCoords(1,1)= p_DvertexCoords(1,p_IverticesAtElement(1,rParticles%p_element(iPart)))
    DcornerCoords(1,2)= p_DvertexCoords(1,p_IverticesAtElement(2,rParticles%p_element(iPart)))
    DcornerCoords(1,3)= p_DvertexCoords(1,p_IverticesAtElement(3,rParticles%p_element(iPart)))
    DcornerCoords(2,1)= p_DvertexCoords(2,p_IverticesAtElement(1,rParticles%p_element(iPart)))
    DcornerCoords(2,2)= p_DvertexCoords(2,p_IverticesAtElement(2,rParticles%p_element(iPart)))
    DcornerCoords(2,3)= p_DvertexCoords(2,p_IverticesAtElement(3,rParticles%p_element(iPart)))

    ! Check if the particle is in the element
    call gaux_isInElement_tri2D(dx,dy,DcornerCoords,binside)

    ! If the particle is still outside the element
    if (binside == .FALSE.) then
      ! Take the old elementnumber
	  rParticles%p_element(iPart) = currentelm

	  ! Check if the is/was a particle-wall-collision
	  call eulerlagrange_partwallcollision(rparlist,p_rproblemLevel,rParticles,iPart)
    end if

  end subroutine eulerlagrange_wrongelement

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_moveparticlestwoway(rparlist,p_rproblemLevel,rsolutionPrimal,rParticles)

!<description>
    ! This subroutine calculates the movement of the particles.

!<input>
    ! Parameterlist
    type(t_parlist), intent(inout) :: rparlist

    ! Particles
    type(t_Particles), intent(inout) :: rParticles

    ! Primal solution vector
    type(t_vectorBlock), intent(inout) :: rsolutionPrimal

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel
    
    ! pointer to elements adjacent to the boundary. 
    !
    ! Handle to 
    !       p_IelementsAtBoundary = array [1..NVBD] of integer.
    ! This array contains a list of all elements on the (real) boundary
    ! in mathematically positive sense.
    ! p_IelementsAtBoundary(i) is the element adjacent to edge
    ! h_IedgesAtBoundary - therefore one element number might appear
    ! more than once in this array!
    ! The boundary elements of boundary component i are saved at
    !        p_IboundaryCpIdx(i)..p_IboundaryCpIdx(i+1)-1.
    ! This is the old KEBD array.
    integer, dimension(:), pointer :: p_IelementsAtBoundary


    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! pointer to the coordinates of the vertices
    !
    ! A list of all corner(!)-vertices of the elements in the triangulation.
    ! Handle to 
    !       p_RcornerCoordinates = array [1..ndim,1..NVT] of double
    ! with
    !   p_DvertexCoords(1,.) = X-coordinate.
    ! for 1D meshes,
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    ! for 2D meshes and
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    !   p_DvertexCoords(3,.) = Z-coordinate.
    ! for 3D meshes.
    ! This is a handle to the old DCORVG-array.
    !
    ! Note that the array may be longer than NVT in general!
    ! (May happen in case of a mesh hierarchy generated by a 2-level
    ! refinement, where the coordinates of the points on the
    ! coarser levels are contained in te coordinates of the
    ! finer levels.)
    ! In such a case, only the first NVT n-tuples in this array are valid!
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! Pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! Coordinates of the three cornervertices
	real(DP) :: ux1_part, uy1_part, ux2_part, uy2_part, ux3_part, uy3_part
	! Density of the gas in the Cornervertices
	real(DP), dimension(3) :: rho_gas
	
    real(DP) :: rho_g, C_W, Re_p, Velo_rel, dt, c_pi
    
    type(t_vectorScalar) :: rvector1, rvector2, rvector3
    real(DP), dimension(:), pointer :: p_Ddata1, p_Ddata2, p_Ddata3
    integer :: iPart, istartpos

    ! Startingpostions of the particles
    real(DP) :: partxmin, partxmax, partymin, partymax

    ! Variables for particle-diameter, -mass and -temperature
    real(DP) :: particledensity, particledensitymin, particledensitymax
    real(DP) :: particlediam, particlediammin, particlediammax
    real(DP) :: parttemp, parttempmin, parttempmax
    integer :: idensitypart, idiampart, itemppart
    
    ! Velocity of the particles
    real(DP) :: velopartx, veloparty, random1, random2, random3

    ! Scalars for the velocity of the gas phase
    real(DP) :: velogasx, velogasy

    ! Volume fraction of the particles in the current position of the particle
    real(DP) :: dVolFrac
    
    ! Heat transfer, Nusselt number, Prandtl number and gas temperatue in the position of the particle 
    real(DP) :: HeatTransfer_T, Nusselt, Prandtl, Temp_gas
    
    ! Thermal conductivity and heat capacity at constant pressure
    real(DP) :: ThermConductivity_g, HeatCapa_g

    ! Energies of the gas and the particles
    real(DP) :: E_intern_gas, E_intern_part
    real(DP) :: E_total_gas, E_total_part

    ! Specific heats at constant volume
    real(DP) :: c_v_gas, c_v_part
    
    ! Forces for the movement
    real(DP), dimension(2) :: F_ges, F_D, F_G
    
    ! Volumefraction of the position of the particle
    real(DP) ::  alpha_p

    ! Variables for starting position from PGM-file
    integer, dimension(:,:), pointer :: p_Idata
    real(DP) :: x,y,xmin,ymin,xmax,ymax
    integer :: nvt,ix,iy,ivt
    type(t_pgm) :: rpgm
    real(DP), dimension(:), pointer :: p_Ddata
    character(LEN=SYS_STRLEN) :: ssolutionname

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation
 
    ! Get vertices of the elements
    call storage_getbase_int2D(&
         p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    ! Get coordinates of the vertices
    call storage_getbase_double2d (&
        p_rtriangulation%h_DvertexCoords,p_DvertexCoords)

    ! Initialize local variables
    rho_g= 0.0_dp 
    C_W=0.0_dp 
    Re_p=0.0_dp 
    Velo_rel=0.0_dp
    dt=0.0_dp
    c_pi=0.0_dp
    parttemp= 0.0_dp
    HeatTransfer_T= 0.0_dp
    Nusselt= 0.0_dp 
    Prandtl= 0.0_dp
    Temp_gas= 0.0_dp
    ThermConductivity_g= 0.0_dp
    HeatCapa_g= 0.0_dp

    ! Get values for the startingpositions of the particles
    call parlst_getvalue_double(rparlist, 'Timestepping', "dinitialStep", dt)

    ! Get scalars for gasvelocity
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "velogasx", velogasx)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "velogasy", velogasy)

    ! Get data from solution
    call eulerlagrange_getVariable(rsolutionPrimal, 'velocity_x', rvector1)
    call eulerlagrange_getVariable(rsolutionPrimal, 'velocity_y', rvector2)
    call eulerlagrange_getVariable(rsolutionPrimal, 'density', rvector3)
    call lsyssc_getbase_double(rvector1, p_Ddata1)
    call lsyssc_getbase_double(rvector2, p_Ddata2)
    call lsyssc_getbase_double(rvector3, p_Ddata3)
 
    c_pi= 3.14159265358979323846264338327950288_dp


    ! Loop over the particles
    do iPart = 1, rParticles%nPart

	! Store old data
	rParticles%p_xpos_old(iPart)=	   rParticles%p_xpos(iPart)
	rParticles%p_ypos_old(iPart)=	   rParticles%p_ypos(iPart)
	rParticles%p_xvelo_old(iPart)=	   rParticles%p_xvelo(iPart)
	rParticles%p_yvelo_old(iPart)=	   rParticles%p_yvelo(iPart)
	rParticles%p_xvelo_gas_old(iPart)= rParticles%p_xvelo_gas(iPart)
	rParticles%p_yvelo_gas_old(iPart)= rParticles%p_xvelo_gas(iPart)

	! Velocity and density of the gas in the first corner (in mathematically positive sense)
	ux1_part= p_Ddata1(p_IverticesAtElement(1,rParticles%p_element(iPart)))
	uy1_part= p_Ddata2(p_IverticesAtElement(1,rParticles%p_element(iPart)))
	rho_gas(1)= p_Ddata3(p_IverticesAtElement(1,rParticles%p_element(iPart)))

	! Velocity and density of the gas in the second corner (in mathematically positive sense)
	ux2_part= p_Ddata1(p_IverticesAtElement(2,rParticles%p_element(iPart)))
	uy2_part= p_Ddata2(p_IverticesAtElement(2,rParticles%p_element(iPart)))
	rho_gas(2)= p_Ddata3(p_IverticesAtElement(2,rParticles%p_element(iPart)))

	! Velocity and density of the gas in the third corner (in mathematically positive sense)
	ux3_part= p_Ddata1(p_IverticesAtElement(3,rParticles%p_element(iPart)))
	uy3_part= p_Ddata2(p_IverticesAtElement(3,rParticles%p_element(iPart)))
	rho_gas(3)= p_Ddata3(p_IverticesAtElement(3,rParticles%p_element(iPart)))

	! Calculate velocity of the gas
	rParticles%p_xvelo_gas(iPart)= 	rParticles%p_lambda1(iPart)*ux1_part + &
									rParticles%p_lambda2(iPart)*ux2_part + &
									rParticles%p_lambda3(iPart)*ux3_part 
	rParticles%p_yvelo_gas(iPart)= 	rParticles%p_lambda1(iPart)*uy1_part + &
									rParticles%p_lambda2(iPart)*uy2_part + &
									rParticles%p_lambda3(iPart)*uy3_part

    ! Scaling the velocity of the gasphase
    rParticles%p_xvelo_gas(iPart)=rParticles%p_xvelo_gas(iPart)*velogasx
    rParticles%p_yvelo_gas(iPart)=rParticles%p_yvelo_gas(iPart)*velogasy

	! Calculate the density of the gas in the position of the particle
	rho_g= 	rParticles%p_lambda1(iPart)*rho_gas(1) + rParticles%p_lambda2(iPart)*&
	        rho_gas(2) + rParticles%p_lambda3(iPart)*rho_gas(3) 


	! Calculate the relative velocity
	Velo_rel= sqrt((rParticles%p_xvelo_old(iPart)-rParticles%p_xvelo_gas(iPart))**2.0_dp +&
	               (rParticles%p_yvelo_old(iPart)-rParticles%p_yvelo_gas(iPart))**2.0_dp)


	! Calculate particle Reynoldsnumber
	!
	! Re_p= \frac{d_p\ \left|\textbf{u}_g-\textbf{u}_p\right|}{\nu_g}
	! with \nu_g=\frac{\eta_g}{\rho_g}
	!
    Re_p= rho_g*0.5_dp*rParticles%p_diam(iPart)*Velo_rel/rParticles%nu_g

    ! get values for the energie calculation
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "temp_gas", Temp_gas)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "c_v_gas", c_v_gas)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "c_v_part", c_v_part)

    ! Calculate internal energies for phase k (gas and particles)
    !
    ! e_k = c_vk * T_k
    !
    E_intern_gas= c_v_gas*Temp_gas
    E_intern_part= c_v_part*rParticles%p_temp(iPart)

    ! Calculate total energies for phase k (gas and particles)
    ! 
    ! E_k= e_k + \frac{1}{2}|u_k|^2
    !
    E_total_gas= E_intern_gas + (rParticles%p_xvelo_gas(iPart)**2.0_dp + &
                    rParticles%p_yvelo_gas(iPart)**2.0_dp)/2.0_dp
    E_total_part= E_intern_part + (rParticles%p_xvelo(iPart)**2.0_dp + &
                    rParticles%p_yvelo(iPart)**2.0_dp)/2.0_dp


	! Calculate the drag force coefficient
	if (Re_p<1000) then
		C_W= 24.0_dp/Re_p*(1.0_dp+0.15_dp*Re_p**0.687_dp)
	else
		C_W= 0.44_dp   !24.0_dp/Re_p
	end if

	! Calculate alpha_n
	rParticles%p_alpha_n(iPart)= C_W*c_pi*rho_g/8.0_dp 

    alpha_p= 0.00002_dp
    if (Velo_rel.ge.500) pause
    ! Compute the dragforce
    ! F_D= \frac{3}{4} * C_W * \frac{\alpha_p \rho_g}{d_p} |u_g - u_p| (u_g - u_p) 
    F_D(1)= 3*C_W*alpha_p*rho_g*Velo_rel*(rParticles%p_xvelo_gas(iPart)-&
            rParticles%p_xvelo_old(iPart))/(rParticles%p_diam(iPart)*4)
    F_D(2)= 3*C_W*alpha_p*rho_g*Velo_rel*(rParticles%p_yvelo_gas(iPart)-&
            rParticles%p_yvelo_old(iPart))/(rParticles%p_diam(iPart)*4)

    ! Compute the gravity
    F_G(1)= rParticles%gravity(1)*rParticles%p_mass(iPart)
    F_G(2)= rParticles%gravity(2)*rParticles%p_mass(iPart)
    
    ! Set force
    F_ges(1)= F_D(1) + F_G(1)
    F_ges(2)= F_D(2) + F_G(2)

    ! Compute new velocity of the particles
    ! a = F_ges/m_p
    rParticles%p_xvelo(iPart)=  rParticles%p_xvelo_old(iPart)+dt*F_ges(1)/rParticles%p_mass(iPart)
    rParticles%p_yvelo(iPart)=  rParticles%p_yvelo_old(iPart)+dt*F_ges(2)/rParticles%p_mass(iPart)

	!---------------------------------------------------------------------------------
	! Calculate the new position of the particle
    !
	! x_new= x_old + delta t * ux_old
	!---------------------------------------------------------------------------------
	
	rParticles%p_xpos(iPart) = rParticles%p_xpos_old(iPart) + dt * rParticles%p_xvelo(iPart)
	rParticles%p_ypos(iPart) = rParticles%p_ypos_old(iPart) + dt * rParticles%p_yvelo(iPart)
	
	! If the particle comes to the outlet
	if (rParticles%p_xpos(iPart).ge.rParticles%maxvalx) then
	
	    ! get values for the startingpositions of the particles
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "xmin", partxmin)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "xmax", partxmax)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "ymin", partymin)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "ymax", partymax)
 
        ! get particlevelocity
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "velopartx", velopartx)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "veloparty", veloparty)
		
        ! Get particle-mass, -temp and -diameter
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particledensity", particledensity)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particlediam", particlediam)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "parttemp", parttemp)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particledensitymin", particledensitymin)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particlediammin", particlediammin)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "parttempmin", parttempmin)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particledensitymax", particledensitymax)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particlediammax", particlediammax)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "parttempmax", parttempmax)
	
	    ! get variable for startingposition
        call parlst_getvalue_int(rparlist, 'Eulerlagrange', "startpos", istartpos)

        ! Get variable for mass of the particles
        call parlst_getvalue_int(rparlist, 'Eulerlagrange', "idensitypart", idensitypart)

        ! Get variable for diameter of the particles
        call parlst_getvalue_int(rparlist, 'Eulerlagrange', "idiampart", idiampart)

        ! Get variable for temperature of the particles
        call parlst_getvalue_int(rparlist, 'Eulerlagrange', "itemppart", itemppart)

        ! Initialisation for starting position from PGM-file
        if (istartpos == 2) then
            ! Get global configuration from parameter list
            call parlst_getvalue_string(rparlist,&
                  'Eulerlagrange', 'filestartpoints', ssolutionName)

            ! Initialize solution from portable graymap image
            call ppsol_readPGM(0, ssolutionName, rpgm)

            ! Set pointer for image data
            call storage_getbase_int2D(rpgm%h_Idata, p_Idata)
            
            ! Determine minimum/maximum values of array
            xmin = huge(DP); xmax = -huge(DP)
            ymin = huge(DP); ymax = -huge(DP)

            do ivt = 1, p_rtriangulation%nvt
                xmin = min(xmin, partxmin)
                xmax = max(xmax, partxmax)
                ymin = min(ymin, partymin)
                ymax = max(ymax, partymax)
            end do

        end if

	    select case(istartpos)
        case(0)
  		  ! Get randomnumber
		  call random_number(random1)
		  call random_number(random2)
		  
          partxmin= minval(p_DvertexCoords(1,:))
          partxmax= minval(p_DvertexCoords(1,:))+&
                    0.2_dp*(maxval(p_DvertexCoords(1,:))-minval(p_DvertexCoords(1,:)))
          partymin= minval(p_DvertexCoords(2,:))
          partymax= maxval(p_DvertexCoords(2,:))

          ! Set startingpositions of the particle
          rParticles%p_xpos(iPart)= partxmin + random1*(partxmax - partxmin)
          rParticles%p_ypos(iPart)= partymin + random2*0.999_dp*(partymax - partymin)
          rParticles%p_xpos_old(iPart)= partxmin + random1*(partxmax - partxmin)
          rParticles%p_ypos_old(iPart)= partymin + random2*0.999_dp*(partymax - partymin)


        case(1)
  		  ! Get randomnumber
		  call random_number(random1)
		  call random_number(random2)
		  
          ! Set startingpositions of the particle
          rParticles%p_xpos(iPart)= partxmin + random1*(partxmax - partxmin)
          rParticles%p_ypos(iPart)= partymin + random2*0.999_dp*(partymax - partymin)
          rParticles%p_xpos_old(iPart)= partxmin + random1*(partxmax - partxmin)
          rParticles%p_ypos_old(iPart)= partymin + random2*0.999_dp*(partymax - partymin)
        
        case(2)
         call random_number(random3)
  
         do
            ! Get random numbers
            call random_number(random1)
		    call random_number(random2)
	
	        partxmin= minval(p_DvertexCoords(1,:))
            partxmax= minval(p_DvertexCoords(1,:))+&
                     (maxval(p_DvertexCoords(1,:))-minval(p_DvertexCoords(1,:)))
            partymin= minval(p_DvertexCoords(2,:))
            partymax= maxval(p_DvertexCoords(2,:))
	    
		    ! Get point in the array
            rParticles%p_xpos(iPart)= partxmin + random1*(partxmax - partxmin)
            rParticles%p_ypos(iPart)= partymin + random2*(partymax - partymin)

            ix = 1+(rpgm%width-1)*(rParticles%p_xpos(iPart)-xmin)/(xmax-xmin)
            if (ix .lt. 1 .or. ix .gt. rpgm%width) cycle

            iy = rpgm%height-(rpgm%height-1)*(rParticles%p_ypos(iPart)-ymin)/(ymax-ymin)
            if (iy .lt. 1 .or. iy .gt. rpgm%height) cycle

            ! If there can be particles, exit loop
            if (random3 .le. real(p_Idata(ix,iy),DP)/real(rpgm%maxgray,DP)) exit
          end do
        
          ! Set particle positions
          rParticles%p_xpos_old(iPart)= rParticles%p_xpos(iPart)
          rParticles%p_ypos_old(iPart)= rParticles%p_ypos(iPart)
         
        case(3)
         
          ! Get random numbers
          call random_number(random1)
          call random_number(random2)
               
          partymin= minval(p_DvertexCoords(2,:))
          partymax= maxval(p_DvertexCoords(2,:))
          partxmin= minval(p_DvertexCoords(1,:))
        
          ! Set startingpositions of the particle
          rParticles%p_xpos(iPart)= partxmin + random1*dt*rParticles%p_xvelo(iPart)
          rParticles%p_ypos(iPart)= partymin + random2*0.999_dp*(partymax - partymin)
          rParticles%p_xpos_old(iPart)= rParticles%p_xpos(iPart)
          rParticles%p_ypos_old(iPart)= rParticles%p_ypos(iPart)
          
         case default
            call output_line('Invalid starting position!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'flagship_startpos')
            call sys_halt()
             
        end select
  
         ! Set diameter of the particles
        select case(idiampart)
        case (0)
            rParticles%p_diam(iPart)= particlediam
            
        case (1)
            ! Get random number
            call random_number(random1)
            
            rParticles%p_diam(iPart)= random1*particlediam

        case (2)
            ! Get random number
            call random_number(random1)
            
            rParticles%p_diam(iPart)= particlediammin+random1*(particlediammax-particlediammin)
            
        case default
          call output_line('Invalid mass type mode!', &
                           OU_CLASS_ERROR,OU_MODE_STD,'flagship_diamtype')
          call sys_halt()
        end select

        ! Set mass of the particles
        select case(idensitypart)
        case (0)
            !Set particle density
            rParticles%p_density(iPart)= particledensity
            
        case (1)
            ! Get random number
            call random_number(random2)
            
            !Set particle density
            rParticles%p_density(iPart)= random2*particledensity
            
        case (2)
            ! Get random number
            call random_number(random2)
            
            !Set particle density
            rParticles%p_density(iPart)= particledensitymin+random2*(particledensitymax-particledensitymin)

        case default
          call output_line('Invalid diameter type mode!', &
                           OU_CLASS_ERROR,OU_MODE_STD,'flagship_masstype')
          call sys_halt()
        end select
 
        ! Set temperature of the particles
        select case(itemppart)
        case (0)
            ! Set particle temperature (all particles have the same temperatur)
            rParticles%p_temp(iPart)= parttemp
            
        case (1)
            ! Get random number
            call random_number(random2)
            
            ! Set particle temperature
            rParticles%p_temp(iPart)= random2*parttemp
            
        case (2)
            ! Get random number
            call random_number(random2)
            
            ! Set particle temperature (between tempmin an tempmax)
            rParticles%p_temp(iPart)= parttempmin+random2*(parttempmax-parttempmin)

        case default
          call output_line('Invalid temp type mode!', &
                           OU_CLASS_ERROR,OU_MODE_STD,'flagship_temptype')
          call sys_halt()
        end select
 
        ! Set initial values for the particle
        rParticles%p_xvelo(iPart)= velopartx
        rParticles%p_yvelo(iPart)= veloparty
        rParticles%p_xvelo_old(iPart)= velopartx
        rParticles%p_yvelo_old(iPart)= veloparty
        rParticles%p_xvelo_gas(iPart)= 0.0_dp
        rParticles%p_yvelo_gas(iPart)= 0.0_dp
        rParticles%p_xvelo_gas_old(iPart)= 0.0_dp
        rParticles%p_yvelo_gas_old(iPart)= 0.0_dp
        rParticles%p_alpha_n(iPart)= 0.0_dp
        rParticles%p_element(iPart)= 1
        rParticles%p_bdy_time(iPart)= 0.0_dp

 
 
         ! Set particle mass
        rParticles%p_mass(iPart)= &
               rParticles%p_density(iPart)*(rParticles%p_diam(iPart)**2 * 3.14159265358_dp /4.0_dp)

 
        if (istartpos == 2) then
            ! Release portable graymap image
            call ppsol_releasePGM(rpgm)
        end if

        ! Find the new start element for the particle
        call eulerlagrange_findelement(rparlist,p_rproblemLevel,rParticles,iPart)

        ! Calculate barycentric coordinates
        call eulerlagrange_calcbarycoords(p_rproblemLevel,rParticles,iPart)

        ! Wrong element
        if ((abs(rParticles%p_lambda1(iPart))+abs(rParticles%p_lambda2(iPart))+&
                  abs(rParticles%p_lambda3(iPart))-1) .ge. 0.00001) then
            call eulerlagrange_wrongelement(rparlist,p_rproblemLevel,rParticles,iPart)
        end if

	end if	

    if (rParticles%p_xpos(iPart) .le. minval(p_DvertexCoords(1,:))) rParticles%p_xpos(iPart) = minval(p_DvertexCoords(1,:))
    if (rParticles%p_xpos(iPart) .ge. maxval(p_DvertexCoords(1,:))) rParticles%p_xpos(iPart) = maxval(p_DvertexCoords(1,:))
    if (rParticles%p_ypos(iPart) .le. minval(p_DvertexCoords(2,:))) rParticles%p_ypos(iPart) = minval(p_DvertexCoords(2,:))
    if (rParticles%p_ypos(iPart) .ge. maxval(p_DvertexCoords(2,:))) rParticles%p_ypos(iPart) = maxval(p_DvertexCoords(2,:))

	end do
	
    ! Release temporal data
    call lsyssc_releasevector(rvector1)
    call lsyssc_releasevector(rvector2)
    call lsyssc_releasevector(rvector3)

  end subroutine eulerlagrange_moveparticlestwoway

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_moveparticlesoneway(rparlist,p_rproblemLevel,rsolutionPrimal,rParticles)

!<description>
    ! This subroutine calculates the movement of the particles.

!<input>
    ! Parameterlist
    type(t_parlist), intent(inout) :: rparlist

    ! Particles
    type(t_Particles), intent(inout) :: rParticles

    ! Primal solution vector
    type(t_vectorBlock), intent(inout) :: rsolutionPrimal

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel
    
    ! pointer to elements adjacent to the boundary. 
    !
    ! Handle to 
    !       p_IelementsAtBoundary = array [1..NVBD] of integer.
    ! This array contains a list of all elements on the (real) boundary
    ! in mathematically positive sense.
    ! p_IelementsAtBoundary(i) is the element adjacent to edge
    ! h_IedgesAtBoundary - therefore one element number might appear
    ! more than once in this array!
    ! The boundary elements of boundary component i are saved at
    !        p_IboundaryCpIdx(i)..p_IboundaryCpIdx(i+1)-1.
    ! This is the old KEBD array.
    integer, dimension(:), pointer :: p_IelementsAtBoundary


    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! pointer to the coordinates of the vertices
    !
    ! A list of all corner(!)-vertices of the elements in the triangulation.
    ! Handle to 
    !       p_RcornerCoordinates = array [1..ndim,1..NVT] of double
    ! with
    !   p_DvertexCoords(1,.) = X-coordinate.
    ! for 1D meshes,
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    ! for 2D meshes and
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    !   p_DvertexCoords(3,.) = Z-coordinate.
    ! for 3D meshes.
    ! This is a handle to the old DCORVG-array.
    !
    ! Note that the array may be longer than NVT in general!
    ! (May happen in case of a mesh hierarchy generated by a 2-level
    ! refinement, where the coordinates of the points on the
    ! coarser levels are contained in te coordinates of the
    ! finer levels.)
    ! In such a case, only the first NVT n-tuples in this array are valid!
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! Pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! Coordinates of the three cornervertices
	real(DP) :: ux1_part, uy1_part, ux2_part, uy2_part, ux3_part, uy3_part
	! Density of the gas in the Cornervertices
	real(DP), dimension(3) :: rho_gas
	
	! Current particle number
	integer :: iPart
	
	! type for solution
	integer :: isolutionpart
	
	! current element
	integer :: currentElement
	
    real(DP) :: rho_g, C_W, Re_p, Velo_rel, dt, c_pi
 
    ! Startingpostions of the particles
    real(DP) :: partxmin, partxmax, partymin, partymax

    ! Variables for particle-diameter, -mass and -temperature
    real(DP) :: particledensity, particledensitymin, particledensitymax
    real(DP) :: particlediam, particlediammin, particlediammax
    real(DP) :: parttemp, parttempmin, parttempmax
    integer :: idensitypart, idiampart, itemppart, istartpos

         
    ! Constants for Runge-Kutta
    real(DP), dimension(2) :: rk_k1, rk_k2, rk_k3, rk_k4
    
    ! Velocity vectors for Runge-Kutta
    real(DP), dimension(2) :: v1_i1, v2_i1, v3_i1
    
    ! Position vectors for Runge-Kutta
    real(DP), dimension(2) :: p1_i1, p2_i1, p3_i1
    
    ! Forces for particle movement
    real(DP), dimension(2) :: F_G, F_D, F_ges
    
    ! Scalar for velocity of the particles and the gas
    real(DP) :: velopartx, veloparty
    real(DP) :: velogasx, velogasy
     
    ! Variables for the velocity of the gas phase in the corner vertices
    type(t_vectorScalar) :: rvector1, rvector2, rvector3
    real(DP), dimension(:), pointer :: p_Ddata1, p_Ddata2, p_Ddata3
    
    ! Variables for randome numbers
    real(DP) :: random1, random2, random3

    ! Coordinates and velocity of the gas phase in the current position
    real(DP), dimension(2) :: Dcurrpos 
    real(DP), dimension(3) :: Dbarycoords

    ! Variables for starting position from PGM-file
    integer, dimension(:,:), pointer :: p_Idata
    real(DP) :: x,y,xmin,ymin,xmax,ymax
    integer :: nvt,ix,iy,ivt
    type(t_pgm) :: rpgm
    real(DP), dimension(:), pointer :: p_Ddata
    character(LEN=SYS_STRLEN) :: ssolutionname

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation
 
    ! Get vertices of the elements
    call storage_getbase_int2D(&
         p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    ! Get coordinates of the vertices
    call storage_getbase_double2d (&
        p_rtriangulation%h_DvertexCoords,p_DvertexCoords)

    ! Get values for the startingpositions of the particles
    call parlst_getvalue_double(rparlist, 'Timestepping', "dinitialStep", dt)

    ! Get values for the startingpositions of the particles
    call parlst_getvalue_int(rparlist, 'Eulerlagrange', "isolutionpart", isolutionpart)

    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "velogasx", velogasx)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "velogasy", velogasy)


    ! Get data from solution
    call eulerlagrange_getVariable(rsolutionPrimal, 'velocity_x', rvector1)
    call eulerlagrange_getVariable(rsolutionPrimal, 'velocity_y', rvector2)
    call eulerlagrange_getVariable(rsolutionPrimal, 'density', rvector3)
    call lsyssc_getbase_double(rvector1, p_Ddata1)
    call lsyssc_getbase_double(rvector2, p_Ddata2)
    call lsyssc_getbase_double(rvector3, p_Ddata3)


    ! Loop over the particles
    do iPart = 1, rParticles%nPart
   
  	! Store old data
	!rParticles%p_xpos_old(iPart)=	   rParticles%p_xpos(iPart)
	!rParticles%p_ypos_old(iPart)=	   rParticles%p_ypos(iPart)
	rParticles%p_xvelo_old(iPart)=	   rParticles%p_xvelo(iPart)
	rParticles%p_yvelo_old(iPart)=	   rParticles%p_yvelo(iPart)
	rParticles%p_xvelo_gas_old(iPart)= rParticles%p_xvelo_gas(iPart)
	rParticles%p_yvelo_gas_old(iPart)= rParticles%p_xvelo_gas(iPart)

	! Velocity and density of the gas in the first corner (in mathematically positive sense)
	ux1_part= p_Ddata1(p_IverticesAtElement(1,rParticles%p_element(iPart)))
	uy1_part= p_Ddata2(p_IverticesAtElement(1,rParticles%p_element(iPart)))
	rho_gas(1)= p_Ddata3(p_IverticesAtElement(1,rParticles%p_element(iPart)))

	! Velocity and density of the gas in the second corner (in mathematically positive sense)
	ux2_part= p_Ddata1(p_IverticesAtElement(2,rParticles%p_element(iPart)))
	uy2_part= p_Ddata2(p_IverticesAtElement(2,rParticles%p_element(iPart)))
	rho_gas(2)= p_Ddata3(p_IverticesAtElement(2,rParticles%p_element(iPart)))

	! Velocity and density of the gas in the third corner (in mathematically positive sense)
	ux3_part= p_Ddata1(p_IverticesAtElement(3,rParticles%p_element(iPart)))
	uy3_part= p_Ddata2(p_IverticesAtElement(3,rParticles%p_element(iPart)))
	rho_gas(3)= p_Ddata3(p_IverticesAtElement(3,rParticles%p_element(iPart)))

	! Calculate velocity of the gas
	rParticles%p_xvelo_gas(iPart)= 	rParticles%p_lambda1(iPart)*ux1_part + &
									rParticles%p_lambda2(iPart)*ux2_part + &
									rParticles%p_lambda3(iPart)*ux3_part 
	rParticles%p_yvelo_gas(iPart)= 	rParticles%p_lambda1(iPart)*uy1_part + &
									rParticles%p_lambda2(iPart)*uy2_part + &
									rParticles%p_lambda3(iPart)*uy3_part
   
	! Calculate density of the gas
	rho_g= 	rParticles%p_lambda1(iPart)*rho_gas(1) + &
			rParticles%p_lambda2(iPart)*rho_gas(2) + &
			rParticles%p_lambda3(iPart)*rho_gas(3) 

    ! Set velocity of the particle in the current psoition
    rParticles%p_xvelo(iPart)= rParticles%p_xvelo_gas(iPart)
    rParticles%p_yvelo(iPart)= rParticles%p_yvelo_gas(iPart)

    ! Set current position as the particle position
    Dcurrpos(1) = rParticles%p_xpos(iPart)
    Dcurrpos(2) = rParticles%p_ypos(iPart)

    call eulerlagrange_getbarycoordelm(rparlist,p_rproblemLevel,Dcurrpos,Dbarycoords,&
            rParticles,currentElement,iPart)

    ! Set temperature of the particles
    select case(isolutionpart)
    case (1)  
        !*******************************************************************
        ! explicit Euler method
        ! x_i+1 = x_i + \Delta t v_i
        !*******************************************************************

        ! Calculate the relative velocity
        Velo_rel= sqrt((rParticles%p_xvelo(iPart)-rParticles%p_xvelo_old(iPart))**2.0_dp +&
                       (rParticles%p_yvelo(iPart)-rParticles%p_yvelo_old(iPart))**2.0_dp)


        ! Calculate particle Reynoldsnumber
        !
        ! Re_p= \frac{d_p\ \left|\textbf{u}_g-\textbf{u}_p\right|}{\nu_g}
        ! with \nu_g=\frac{\eta_g}{\rho_g}
        !
        Re_p= rho_g*0.5_dp*rParticles%p_diam(iPart)*Velo_rel/rParticles%nu_g


        ! Calculate the drag force coefficient
        if (Re_p<1000) then
	        C_W= 24.0_dp/Re_p*(1.0_dp+0.15_dp*Re_p**0.687_dp)
        else
	        C_W= 24.0_dp/Re_p
        end if

        if (Velo_rel.ge.500) pause

        ! Compute the dragforce
        F_D(1)= C_W*3.14*rho_g*rParticles%p_diam(iPart)**2*Velo_rel*&
                (rParticles%p_yvelo(iPart)-rParticles%p_yvelo_old(iPart))/8.0_dp
        F_D(2)= C_W*3.14*rho_g*rParticles%p_diam(iPart)**2*Velo_rel*&
                (rParticles%p_yvelo(iPart)-rParticles%p_yvelo_old(iPart))/8.0_dp

        ! Compute the gravity
        F_G(1)= rParticles%gravity(1)*rParticles%p_mass(iPart)
        F_G(2)= rParticles%gravity(2)*rParticles%p_mass(iPart)

        ! Set force
        F_ges(1)= F_D(1) + F_G(1)
        F_ges(2)= F_D(2) + F_G(2)

        ! Compute new velocity of the particles
        rParticles%p_xvelo(iPart)= rParticles%p_xvelo_old(iPart)+dt*F_ges(1)/rParticles%p_mass(iPart)
        rParticles%p_yvelo(iPart)= rParticles%p_yvelo_old(iPart)+dt*F_ges(2)/rParticles%p_mass(iPart)
        
        ! Compute the new position
        rParticles%p_xpos(iPart)= rParticles%p_xpos_old(iPart) + dt* rParticles%p_xvelo(iPart)/velogasx
        rParticles%p_ypos(iPart)= rParticles%p_ypos_old(iPart) + dt* rParticles%p_yvelo(iPart)/velogasy

                             
    case (2)
        !*******************************************************************
        ! improved Euler method
        ! x_i+1 = x_i + \frac12 \Delta t (v_i + v_i+1)
        !*******************************************************************
                
        ! Calculate the relative velocity
        Velo_rel= sqrt((rParticles%p_xvelo(iPart)-rParticles%p_xvelo_old(iPart))**2.0_dp +&
                       (rParticles%p_yvelo(iPart)-rParticles%p_yvelo_old(iPart))**2.0_dp)


        ! Calculate particle Reynoldsnumber
        !
        ! Re_p= \frac{d_p\ \left|\textbf{u}_g-\textbf{u}_p\right|}{\nu_g}
        ! with \nu_g=\frac{\eta_g}{\rho_g}
        !
        Re_p= rho_g*0.5_dp*rParticles%p_diam(iPart)*Velo_rel/rParticles%nu_g


        ! Calculate the drag force coefficient
        if (Re_p<1000) then
	        C_W= 24.0_dp/Re_p*(1.0_dp+0.15_dp*Re_p**0.687_dp)
        else
	        C_W= 24.0_dp/Re_p
        end if

        if (Velo_rel.ge.500) pause

        ! Compute the dragforce
        F_D(1)= C_W*3.14*rho_g*rParticles%p_diam(iPart)**2*Velo_rel*&
                (rParticles%p_yvelo(iPart)-rParticles%p_yvelo_old(iPart))/8.0_dp
        F_D(2)= C_W*3.14*rho_g*rParticles%p_diam(iPart)**2*Velo_rel*&
                (rParticles%p_yvelo(iPart)-rParticles%p_yvelo_old(iPart))/8.0_dp

        ! Compute the gravity
        F_G(1)= rParticles%gravity(1)*rParticles%p_mass(iPart)
        F_G(2)= rParticles%gravity(2)*rParticles%p_mass(iPart)

        ! Set force
        F_ges(1)= F_D(1) + F_G(1)
        F_ges(2)= F_D(2) + F_G(2)

        ! Compute new velocity of the particles
        rParticles%p_xvelo(iPart)= rParticles%p_xvelo_old(iPart)+dt*F_ges(1)/rParticles%p_mass(iPart)
        rParticles%p_yvelo(iPart)= rParticles%p_yvelo_old(iPart)+dt*F_ges(2)/rParticles%p_mass(iPart)
                
        ! Compute first constant for Runge-Kutta
        rk_k1(1)= rParticles%p_xvelo(iPart)/velogasx
        rk_k1(2)= rParticles%p_yvelo(iPart)/velogasy
        
        ! first position
        p1_i1(1)= rParticles%p_xpos(iPart) + dt*rk_k1(1)
        p1_i1(2)= rParticles%p_ypos(iPart) + dt*rk_k1(2)
 
        ! Get element an barycentric coordinates for the position p1
        call eulerlagrange_getbarycoordelm(rparlist,p_rproblemLevel,p1_i1,Dbarycoords,&
            rParticles,currentElement,iPart)

        ! Velocity and density of the gas in the first corner (in mathematically positive sense)
        ux1_part= p_Ddata1(p_IverticesAtElement(1,currentElement))
        uy1_part= p_Ddata2(p_IverticesAtElement(1,currentElement))
        rho_gas(1)= p_Ddata3(p_IverticesAtElement(1,currentElement))

        ! Velocity and density of the gas in the second corner (in mathematically positive sense)
        ux2_part= p_Ddata1(p_IverticesAtElement(2,currentElement))
        uy2_part= p_Ddata2(p_IverticesAtElement(2,currentElement))
        rho_gas(2)= p_Ddata3(p_IverticesAtElement(2,currentElement))

        ! Velocity and density of the gas in the third corner (in mathematically positive sense)
        ux3_part= p_Ddata1(p_IverticesAtElement(3,currentElement))
        uy3_part= p_Ddata2(p_IverticesAtElement(3,currentElement))
        rho_gas(3)= p_Ddata3(p_IverticesAtElement(3,currentElement))

        ! first velocity
        v1_i1(1)= 	Dbarycoords(1) * ux1_part / rho_gas(1) + &
					Dbarycoords(2) * ux2_part / rho_gas(2) + &
					Dbarycoords(3) * ux3_part / rho_gas(3)
        v1_i1(2)= 	Dbarycoords(1) * uy1_part / rho_gas(1) + &
					Dbarycoords(2) * uy2_part / rho_gas(2) + &
					Dbarycoords(3) * uy3_part / rho_gas(3)
   
        ! Calculate the relative velocity
        Velo_rel= sqrt((rk_k1(1)-v1_i1(1))**2.0_dp +&
                       (rk_k1(2)-v1_i1(2))**2.0_dp)

        ! Calculate particle Reynoldsnumber
        !
        ! Re_p= \frac{d_p\ \left|\textbf{u}_g-\textbf{u}_p\right|}{\nu_g}
        ! with \nu_g=\frac{\eta_g}{\rho_g}
        !
        Re_p= rho_g*0.5_dp*rParticles%p_diam(iPart)*Velo_rel/rParticles%nu_g


        ! Calculate the drag force coefficient
        if (Re_p<1000) then
	        C_W= 24.0_dp/Re_p*(1.0_dp+0.15_dp*Re_p**0.687_dp)
        else
	        C_W= 24.0_dp/Re_p
        end if

        if (Velo_rel.ge.500) pause

        ! Compute the dragforce
        F_D(1)= C_W*3.14*rho_g*rParticles%p_diam(iPart)**2*Velo_rel*&
                (v1_i1(1)-rk_k1(1))/8.0_dp
        F_D(2)= C_W*3.14*rho_g*rParticles%p_diam(iPart)**2*Velo_rel*&
                (v1_i1(2)-rk_k1(2))/8.0_dp


        ! Compute the gravity
        F_G(1)= rParticles%gravity(1)*rParticles%p_mass(iPart)
        F_G(2)= rParticles%gravity(2)*rParticles%p_mass(iPart)

        ! Set force
        F_ges(1)= F_D(1) + F_G(1)
        F_ges(2)= F_D(2) + F_G(2)

        ! Compute new velocity of the particles
        rk_k2(1)=  rk_k1(1)+dt*F_ges(1)/rParticles%p_mass(iPart)/velogasx
        rk_k2(2)=  rk_k1(2)+dt*F_ges(2)/rParticles%p_mass(iPart)/velogasy
        
        ! Compute the new position
        rParticles%p_xpos(iPart)= rParticles%p_xpos_old(iPart) + dt* &
            (rk_k1(1) + v1_i1(1))/2.0_dp
        rParticles%p_ypos(iPart)= rParticles%p_ypos_old(iPart) + dt* &
            (rk_k1(2) + v1_i1(2))/2.0_dp
            
        ! Scale the new velocity
        rParticles%p_xvelo(iPart)=rk_k2(1)*velogasx
        rParticles%p_yvelo(iPart)=rk_k2(2)*velogasy
          
    case (3)
        !*******************************************************************
        ! classic Runge-Kutta algorithmn
        ! x_i+1= x_i + \frac16 \Delta t ( v_i + 2v_i+1^1 + 2v_i+1^2 + v_i+1^3)
        !*******************************************************************
          
        ! Calculate the relative velocity
        Velo_rel= sqrt((rParticles%p_xvelo(iPart)-rParticles%p_xvelo_old(iPart))**2.0_dp +&
                       (rParticles%p_yvelo(iPart)-rParticles%p_yvelo_old(iPart))**2.0_dp)


        ! Calculate particle Reynoldsnumber
        !
        ! Re_p= \frac{d_p\ \left|\textbf{u}_g-\textbf{u}_p\right|}{\nu_g}
        ! with \nu_g=\frac{\eta_g}{\rho_g}
        !
        Re_p= rho_g*0.5_dp*rParticles%p_diam(iPart)*Velo_rel/rParticles%nu_g


        ! Calculate the drag force coefficient
        if (Re_p<1000) then
	        C_W= 24.0_dp/Re_p*(1.0_dp+0.15_dp*Re_p**0.687_dp)
        else
	        C_W= 24.0_dp/Re_p
        end if

        if (Velo_rel.ge.500) pause

        ! Compute the dragforce
        F_D(1)= C_W*3.14*rho_g*rParticles%p_diam(iPart)**2*Velo_rel*&
                (rParticles%p_yvelo(iPart)-rParticles%p_yvelo_old(iPart))/8.0_dp
        F_D(2)= C_W*3.14*rho_g*rParticles%p_diam(iPart)**2*Velo_rel*&
                (rParticles%p_yvelo(iPart)-rParticles%p_yvelo_old(iPart))/8.0_dp

        ! Compute the gravity
        F_G(1)= rParticles%gravity(1)*rParticles%p_mass(iPart)
        F_G(2)= rParticles%gravity(2)*rParticles%p_mass(iPart)

        ! Set force
        F_ges(1)= F_D(1) + F_G(1)
        F_ges(2)= F_D(2) + F_G(2)

        ! Compute new velocity of the particles
        rParticles%p_xvelo(iPart)= rParticles%p_xvelo_old(iPart)+dt*F_ges(1)/rParticles%p_mass(iPart)
        rParticles%p_yvelo(iPart)= rParticles%p_yvelo_old(iPart)+dt*F_ges(2)/rParticles%p_mass(iPart)
                
        ! Compute first constant for Runge-Kutta
        rk_k1(1)= rParticles%p_xvelo(iPart)/velogasx
        rk_k1(2)= rParticles%p_yvelo(iPart)/velogasy
        
        ! first position
        p1_i1(1)= rParticles%p_xpos(iPart) + dt*rk_k1(1)/2.0_dp
        p1_i1(2)= rParticles%p_ypos(iPart) + dt*rk_k1(2)/2.0_dp
 
        ! Get element an barycentric coordinates for the position p1
        call eulerlagrange_getbarycoordelm(rparlist,p_rproblemLevel,p1_i1,Dbarycoords,&
            rParticles,currentElement,iPart)

        ! Velocity and density of the gas in the first corner (in mathematically positive sense)
        ux1_part= p_Ddata1(p_IverticesAtElement(1,currentElement))
        uy1_part= p_Ddata2(p_IverticesAtElement(1,currentElement))
        rho_gas(1)= p_Ddata3(p_IverticesAtElement(1,currentElement))

        ! Velocity and density of the gas in the second corner (in mathematically positive sense)
        ux2_part= p_Ddata1(p_IverticesAtElement(2,currentElement))
        uy2_part= p_Ddata2(p_IverticesAtElement(2,currentElement))
        rho_gas(2)= p_Ddata3(p_IverticesAtElement(2,currentElement))

        ! Velocity and density of the gas in the third corner (in mathematically positive sense)
        ux3_part= p_Ddata1(p_IverticesAtElement(3,currentElement))
        uy3_part= p_Ddata2(p_IverticesAtElement(3,currentElement))
        rho_gas(3)= p_Ddata3(p_IverticesAtElement(3,currentElement))

        ! first velocity
        v1_i1(1)= 	Dbarycoords(1) * ux1_part / rho_gas(1) + &
					Dbarycoords(2) * ux2_part / rho_gas(2) + &
					Dbarycoords(3) * ux3_part / rho_gas(3)
        v1_i1(2)= 	Dbarycoords(1) * uy1_part / rho_gas(1) + &
					Dbarycoords(2) * uy2_part / rho_gas(2) + &
					Dbarycoords(3) * uy3_part / rho_gas(3)
   
        ! Calculate the relative velocity
        Velo_rel= sqrt((rk_k1(1)-v1_i1(1))**2.0_dp +&
                       (rk_k1(2)-v1_i1(2))**2.0_dp)

        ! Calculate particle Reynoldsnumber
        !
        ! Re_p= \frac{d_p\ \left|\textbf{u}_g-\textbf{u}_p\right|}{\nu_g}
        ! with \nu_g=\frac{\eta_g}{\rho_g}
        !
        Re_p= rho_g*0.5_dp*rParticles%p_diam(iPart)*Velo_rel/rParticles%nu_g


        ! Calculate the drag force coefficient
        if (Re_p<1000) then
	        C_W= 24.0_dp/Re_p*(1.0_dp+0.15_dp*Re_p**0.687_dp)
        else
	        C_W= 24.0_dp/Re_p
        end if

        if (Velo_rel.ge.500) pause

        ! Compute the dragforce
        F_D(1)= C_W*3.14*rho_g*rParticles%p_diam(iPart)**2*Velo_rel*&
                (v1_i1(1)-rk_k1(1))/8.0_dp
        F_D(2)= C_W*3.14*rho_g*rParticles%p_diam(iPart)**2*Velo_rel*&
                (v1_i1(2)-rk_k1(2))/8.0_dp


        ! Compute the gravity
        F_G(1)= rParticles%gravity(1)*rParticles%p_mass(iPart)
        F_G(2)= rParticles%gravity(2)*rParticles%p_mass(iPart)

        ! Set force
        F_ges(1)= F_D(1) + F_G(1)
        F_ges(2)= F_D(2) + F_G(2)

        ! Compute new velocity of the particles
        rk_k2(1)=  rk_k1(1)+dt*F_ges(1)/rParticles%p_mass(iPart)/velogasx
        rk_k2(2)=  rk_k1(2)+dt*F_ges(2)/rParticles%p_mass(iPart)/velogasy
        
        ! second position
        p2_i1(1)= rParticles%p_xpos(iPart) + dt*rk_k2(1)/2.0_dp
        p2_i1(2)= rParticles%p_ypos(iPart) + dt*rk_k2(2)/2.0_dp
        
        ! Get element an barycentric coordinates for the position p2
        call eulerlagrange_getbarycoordelm(rparlist,p_rproblemLevel,p2_i1,Dbarycoords,&
            rParticles,currentElement,iPart)
 
         ! Velocity and density of the gas in the first corner (in mathematically positive sense)
        ux1_part= p_Ddata1(p_IverticesAtElement(1,currentElement))
        uy1_part= p_Ddata2(p_IverticesAtElement(1,currentElement))
        rho_gas(1)= p_Ddata3(p_IverticesAtElement(1,currentElement))

        ! Velocity and density of the gas in the second corner (in mathematically positive sense)
        ux2_part= p_Ddata1(p_IverticesAtElement(2,currentElement))
        uy2_part= p_Ddata2(p_IverticesAtElement(2,currentElement))
        rho_gas(2)= p_Ddata3(p_IverticesAtElement(2,currentElement))

        ! Velocity and density of the gas in the third corner (in mathematically positive sense)
        ux3_part= p_Ddata1(p_IverticesAtElement(3,currentElement))
        uy3_part= p_Ddata2(p_IverticesAtElement(3,currentElement))
        rho_gas(3)= p_Ddata3(p_IverticesAtElement(3,currentElement))

        ! second velocity
        v2_i1(1)= 	Dbarycoords(1) * ux1_part / rho_gas(1) + &
					Dbarycoords(2) * ux2_part / rho_gas(2) + &
					Dbarycoords(3) * ux3_part / rho_gas(3)
        v2_i1(2)= 	Dbarycoords(1) * uy1_part / rho_gas(1) + &
					Dbarycoords(2) * uy2_part / rho_gas(2) + &
					Dbarycoords(3) * uy3_part / rho_gas(3)
 
 
         ! Calculate the relative velocity
        Velo_rel= sqrt((rk_k2(1)-v2_i1(1))**2.0_dp +&
                       (rk_k2(2)-v2_i1(2))**2.0_dp)


        ! Calculate particle Reynoldsnumber
        !
        ! Re_p= \frac{d_p\ \left|\textbf{u}_g-\textbf{u}_p\right|}{\nu_g}
        ! with \nu_g=\frac{\eta_g}{\rho_g}
        !
        Re_p= rho_g*0.5_dp*rParticles%p_diam(iPart)*Velo_rel/rParticles%nu_g


        ! Calculate the drag force coefficient
        if (Re_p<1000) then
	        C_W= 24.0_dp/Re_p*(1.0_dp+0.15_dp*Re_p**0.687_dp)
        else
	        C_W= 24.0_dp/Re_p
        end if

        if (Velo_rel.ge.500) pause

        ! Compute the dragforce
        F_D(1)= C_W*3.14*rho_g*rParticles%p_diam(iPart)**2*Velo_rel*&
                (v2_i1(1)-rk_k2(1))/8.0_dp
        F_D(2)= C_W*3.14*rho_g*rParticles%p_diam(iPart)**2*Velo_rel*&
                (v2_i1(2)-rk_k2(2))/8.0_dp


        ! Compute the gravity
        F_G(1)= rParticles%gravity(1)*rParticles%p_mass(iPart)
        F_G(2)= rParticles%gravity(2)*rParticles%p_mass(iPart)

        ! Set force
        F_ges(1)= F_D(1) + F_G(1)
        F_ges(2)= F_D(2) + F_G(2)

        ! Compute new velocity of the particles
        rk_k3(1)=  rk_k2(1)+dt*F_ges(1)/rParticles%p_mass(iPart)/velogasx
        rk_k3(2)=  rk_k2(2)+dt*F_ges(2)/rParticles%p_mass(iPart)/velogasy
        
        ! third position
        p3_i1(1)= rParticles%p_xpos(iPart) + dt*rk_k3(1)
        p3_i1(2)= rParticles%p_ypos(iPart) + dt*rk_k3(2)
        
        ! Get element an barycentric coordinates for the position p3
        call eulerlagrange_getbarycoordelm(rparlist,p_rproblemLevel,p3_i1,Dbarycoords,&
            rParticles,currentElement,iPart)

        ! Velocity and density of the gas in the first corner (in mathematically positive sense)
        ux1_part= p_Ddata1(p_IverticesAtElement(1,currentElement))
        uy1_part= p_Ddata2(p_IverticesAtElement(1,currentElement))
        rho_gas(1)= p_Ddata3(p_IverticesAtElement(1,currentElement))

        ! Velocity and density of the gas in the second corner (in mathematically positive sense)
        ux2_part= p_Ddata1(p_IverticesAtElement(2,currentElement))
        uy2_part= p_Ddata2(p_IverticesAtElement(2,currentElement))
        rho_gas(2)= p_Ddata3(p_IverticesAtElement(2,currentElement))

        ! Velocity and density of the gas in the third corner (in mathematically positive sense)
        ux3_part= p_Ddata1(p_IverticesAtElement(3,currentElement))
        uy3_part= p_Ddata2(p_IverticesAtElement(3,currentElement))
        rho_gas(3)= p_Ddata3(p_IverticesAtElement(3,currentElement))

        ! third velocity
        v3_i1(1)= 	Dbarycoords(1) * ux1_part / rho_gas(1) + &
					Dbarycoords(2) * ux2_part / rho_gas(2) + &
					Dbarycoords(3) * ux3_part / rho_gas(3)
        v3_i1(2)= 	Dbarycoords(1) * uy1_part / rho_gas(1) + &
					Dbarycoords(2) * uy2_part / rho_gas(2) + &
					Dbarycoords(3) * uy3_part / rho_gas(3)

 
         ! Calculate the relative velocity
        Velo_rel= sqrt((rk_k3(1)-v3_i1(1))**2.0_dp +&
                       (rk_k3(2)-v3_i1(2))**2.0_dp)


        ! Calculate particle Reynoldsnumber
        !
        ! Re_p= \frac{d_p\ \left|\textbf{u}_g-\textbf{u}_p\right|}{\nu_g}
        ! with \nu_g=\frac{\eta_g}{\rho_g}
        !
        Re_p= rho_g*0.5_dp*rParticles%p_diam(iPart)*Velo_rel/rParticles%nu_g


        ! Calculate the drag force coefficient
        if (Re_p<1000) then
	        C_W= 24.0_dp/Re_p*(1.0_dp+0.15_dp*Re_p**0.687_dp)
        else
	        C_W= 24.0_dp/Re_p
        end if

        if (Velo_rel.ge.500) pause

        ! Compute the dragforce
        F_D(1)= C_W*3.14*rho_g*rParticles%p_diam(iPart)**2*Velo_rel*&
                (v3_i1(1)-rk_k3(1))/8.0_dp
        F_D(2)= C_W*3.14*rho_g*rParticles%p_diam(iPart)**2*Velo_rel*&
                (v3_i1(2)-rk_k3(2))/8.0_dp


        ! Compute the gravity
        F_G(1)= rParticles%gravity(1)*rParticles%p_mass(iPart)
        F_G(2)= rParticles%gravity(2)*rParticles%p_mass(iPart)

        ! Set force
        F_ges(1)= F_D(1) + F_G(1)
        F_ges(2)= F_D(2) + F_G(2)

        ! Compute new velocity of the particles
        rk_k4(1)=  rk_k3(1)+dt*F_ges(1)/rParticles%p_mass(iPart)/velogasx
        rk_k4(2)=  rk_k3(2)+dt*F_ges(2)/rParticles%p_mass(iPart)/velogasy
 
        ! Set new particle-position
        rParticles%p_xpos(iPart)= rParticles%p_xpos_old(iPart) + &
                dt*(rk_k1(1) + 2 * rk_k2(1) + 2 * rk_k3(1) + rk_k4(1))/6.0_dp
        rParticles%p_ypos(iPart)= rParticles%p_ypos_old(iPart) + &
                dt*(rk_k1(2) + 2 * rk_k2(2) + 2 * rk_k3(2) + rk_k4(2))/6.0_dp

        ! Set direction of the velocita of the particle
        rParticles%p_xvelo(iPart)= rk_k4(1)*velogasx
        rParticles%p_yvelo(iPart)= rk_k4(2)*velogasy


    case default
      call output_line('Invalid solution type!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'flagship_isolutionpart')
      call sys_halt()
    end select
 
 	! If the particle comes to the outlet
	if (rParticles%p_xpos(iPart).ge.rParticles%maxvalx) then
	
	    ! get values for the startingpositions of the particles
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "xmin", partxmin)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "xmax", partxmax)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "ymin", partymin)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "ymax", partymax)
 
        ! get particlevelocity
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "velopartx", velopartx)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "veloparty", veloparty)
		
        ! Get particle-mass, -temp and -diameter
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particledensity", particledensity)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particlediam", particlediam)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "parttemp", parttemp)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particledensitymin", particledensitymin)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particlediammin", particlediammin)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "parttempmin", parttempmin)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particledensitymax", particledensitymax)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particlediammax", particlediammax)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "parttempmax", parttempmax)
	
	    ! get variable for startingposition
        call parlst_getvalue_int(rparlist, 'Eulerlagrange', "startpos", istartpos)

        ! Get variable for mass of the particles
        call parlst_getvalue_int(rparlist, 'Eulerlagrange', "idensitypart", idensitypart)

        ! Get variable for diameter of the particles
        call parlst_getvalue_int(rparlist, 'Eulerlagrange', "idiampart", idiampart)

        ! Get variable for temperature of the particles
        call parlst_getvalue_int(rparlist, 'Eulerlagrange', "itemppart", itemppart)

        ! Initialisation for starting position from PGM-file
        if (istartpos == 2) then
            ! Get global configuration from parameter list
            call parlst_getvalue_string(rparlist,&
                  'Eulerlagrange', 'filestartpoints', ssolutionName)

            ! Initialize solution from portable graymap image
            call ppsol_readPGM(0, ssolutionName, rpgm)

            ! Set pointer for image data
            call storage_getbase_int2D(rpgm%h_Idata, p_Idata)
            
            ! Determine minimum/maximum values of array
            xmin = huge(DP); xmax = -huge(DP)
            ymin = huge(DP); ymax = -huge(DP)

            do ivt = 1, p_rtriangulation%nvt
                xmin = min(xmin, partxmin)
                xmax = max(xmax, partxmax)
                ymin = min(ymin, partymin)
                ymax = max(ymax, partymax)
            end do

        end if

	    select case(istartpos)
        case(0)
  		  ! Get randomnumber
		  call random_number(random1)
		  call random_number(random2)
		  
          partxmin= minval(p_DvertexCoords(1,:))
          partxmax= minval(p_DvertexCoords(1,:))+&
                    0.2_dp*(maxval(p_DvertexCoords(1,:))-minval(p_DvertexCoords(1,:)))
          partymin= minval(p_DvertexCoords(2,:))
          partymax= maxval(p_DvertexCoords(2,:))

          ! Set startingpositions of the particle
          rParticles%p_xpos(iPart)= partxmin + random1*(partxmax - partxmin)
          rParticles%p_ypos(iPart)= partymin + random2*0.999_dp*(partymax - partymin)
          rParticles%p_xpos_old(iPart)= partxmin + random1*(partxmax - partxmin)
          rParticles%p_ypos_old(iPart)= partymin + random2*0.999_dp*(partymax - partymin)


        case(1)
  		  ! Get randomnumber
		  call random_number(random1)
		  call random_number(random2)
		  
          ! Set startingpositions of the particle
          rParticles%p_xpos(iPart)= partxmin + random1*(partxmax - partxmin)
          rParticles%p_ypos(iPart)= partymin + random2*0.999_dp*(partymax - partymin)
          rParticles%p_xpos_old(iPart)= partxmin + random1*(partxmax - partxmin)
          rParticles%p_ypos_old(iPart)= partymin + random2*0.999_dp*(partymax - partymin)
        
        case(2)
         call random_number(random3)
  
         do
            ! Get random numbers
            call random_number(random1)
		    call random_number(random2)
	
	        partxmin= minval(p_DvertexCoords(1,:))
            partxmax= minval(p_DvertexCoords(1,:))+&
                     (maxval(p_DvertexCoords(1,:))-minval(p_DvertexCoords(1,:)))
            partymin= minval(p_DvertexCoords(2,:))
            partymax= maxval(p_DvertexCoords(2,:))
	    
		    ! Get point in the array
            rParticles%p_xpos(iPart)= partxmin + random1*(partxmax - partxmin)
            rParticles%p_ypos(iPart)= partymin + random2*(partymax - partymin)

            ix = 1+(rpgm%width-1)*(rParticles%p_xpos(iPart)-xmin)/(xmax-xmin)
            if (ix .lt. 1 .or. ix .gt. rpgm%width) cycle

            iy = rpgm%height-(rpgm%height-1)*(rParticles%p_ypos(iPart)-ymin)/(ymax-ymin)
            if (iy .lt. 1 .or. iy .gt. rpgm%height) cycle

            ! If there can be particles, exit loop
            if (random3 .le. real(p_Idata(ix,iy),DP)/real(rpgm%maxgray,DP)) exit
          end do
        
          ! Set particle positions
          rParticles%p_xpos_old(iPart)= rParticles%p_xpos(iPart)
          rParticles%p_ypos_old(iPart)= rParticles%p_ypos(iPart)
         
        case(3)
         
          ! Get random numbers
          call random_number(random1)
          call random_number(random2)
               
          partymin= minval(p_DvertexCoords(2,:))
          partymax= maxval(p_DvertexCoords(2,:))
          partxmin= minval(p_DvertexCoords(1,:))
        
          ! Set startingpositions of the particle
          rParticles%p_xpos(iPart)= partxmin + random1*dt*rParticles%p_xvelo(iPart)
          rParticles%p_ypos(iPart)= partymin + random2*0.999_dp*(partymax - partymin)
          rParticles%p_xpos_old(iPart)= rParticles%p_xpos(iPart)
          rParticles%p_ypos_old(iPart)= rParticles%p_ypos(iPart)
          
         case default
            call output_line('Invalid starting position!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'flagship_startpos')
            call sys_halt()
             
        end select
  
         ! Set diameter of the particles
        select case(idiampart)
        case (0)
            rParticles%p_diam(iPart)= particlediam
            
        case (1)
            ! Get random number
            call random_number(random1)
            
            rParticles%p_diam(iPart)= random1*particlediam

        case (2)
            ! Get random number
            call random_number(random1)
            
            rParticles%p_diam(iPart)= particlediammin+random1*(particlediammax-particlediammin)
            
        case default
          call output_line('Invalid mass type mode!', &
                           OU_CLASS_ERROR,OU_MODE_STD,'flagship_diamtype')
          call sys_halt()
        end select

        ! Set mass of the particles
        select case(idensitypart)
        case (0)
            !Set particle density
            rParticles%p_density(iPart)= particledensity
            
        case (1)
            ! Get random number
            call random_number(random2)
            
            !Set particle density
            rParticles%p_density(iPart)= random2*particledensity
            
        case (2)
            ! Get random number
            call random_number(random2)
            
            !Set particle density
            rParticles%p_density(iPart)= particledensitymin+random2*(particledensitymax-particledensitymin)

        case default
          call output_line('Invalid diameter type mode!', &
                           OU_CLASS_ERROR,OU_MODE_STD,'flagship_masstype')
          call sys_halt()
        end select
 
        ! Set temperature of the particles
        select case(itemppart)
        case (0)
            ! Set particle temperature (all particles have the same temperatur)
            rParticles%p_temp(iPart)= parttemp
            
        case (1)
            ! Get random number
            call random_number(random2)
            
            ! Set particle temperature
            rParticles%p_temp(iPart)= random2*parttemp
            
        case (2)
            ! Get random number
            call random_number(random2)
            
            ! Set particle temperature (between tempmin an tempmax)
            rParticles%p_temp(iPart)= parttempmin+random2*(parttempmax-parttempmin)

        case default
          call output_line('Invalid temp type mode!', &
                           OU_CLASS_ERROR,OU_MODE_STD,'flagship_temptype')
          call sys_halt()
        end select
 
        ! Set initial values for the particle
        rParticles%p_xvelo(iPart)= velopartx
        rParticles%p_yvelo(iPart)= veloparty
        rParticles%p_xvelo_old(iPart)= velopartx
        rParticles%p_yvelo_old(iPart)= veloparty
        rParticles%p_xvelo_gas(iPart)= 0.0_dp
        rParticles%p_yvelo_gas(iPart)= 0.0_dp
        rParticles%p_xvelo_gas_old(iPart)= 0.0_dp
        rParticles%p_yvelo_gas_old(iPart)= 0.0_dp
        rParticles%p_alpha_n(iPart)= 0.0_dp
        rParticles%p_element(iPart)= 1
        rParticles%p_bdy_time(iPart)= 0.0_dp

 
 
         ! Set particle mass
        rParticles%p_mass(iPart)= &
               rParticles%p_density(iPart)*(rParticles%p_diam(iPart)**2 * 3.14159265358_dp /4.0_dp)

 
        if (istartpos == 2) then
            ! Release portable graymap image
            call ppsol_releasePGM(rpgm)
        end if

        ! Find the new start element for the particle
        call eulerlagrange_findelement(rparlist,p_rproblemLevel,rParticles,iPart)

        ! Calculate barycentric coordinates
        call eulerlagrange_calcbarycoords(p_rproblemLevel,rParticles,iPart)

        ! Wrong element
        if ((abs(rParticles%p_lambda1(iPart))+abs(rParticles%p_lambda2(iPart))+&
                  abs(rParticles%p_lambda3(iPart))-1) .ge. 0.00001) then
            call eulerlagrange_wrongelement(rparlist,p_rproblemLevel,rParticles,iPart)
        end if

	end if	
     
   	rParticles%p_xpos_old(iPart)=	   rParticles%p_xpos(iPart)
	rParticles%p_ypos_old(iPart)=	   rParticles%p_ypos(iPart)

            
    if (rParticles%p_xpos(iPart) .le. minval(p_DvertexCoords(1,:))) rParticles%p_xpos(iPart) = minval(p_DvertexCoords(1,:))
    if (rParticles%p_xpos(iPart) .ge. maxval(p_DvertexCoords(1,:))) rParticles%p_xpos(iPart) = minval(p_DvertexCoords(1,:))
    if (rParticles%p_ypos(iPart) .le. minval(p_DvertexCoords(2,:))) rParticles%p_ypos(iPart) = minval(p_DvertexCoords(2,:))
    if (rParticles%p_ypos(iPart) .ge. maxval(p_DvertexCoords(2,:))) rParticles%p_ypos(iPart) = maxval(p_DvertexCoords(2,:))

    
    end do  !Loop over all particles


  end subroutine eulerlagrange_moveparticlesoneway

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_getbarycoordelm(rparlist,p_rproblemLevel,Dcurrpos,Dbarycoords,rParticles,currentElement,iPart)

!<description>
    ! This subroutine gets the velocity of the gas in a position.

!<input>
    ! Particles
    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! gas velocity at he position
    real(DP), dimension(2), intent(inout) :: Dcurrpos

    ! gas velocity at he position
    real(DP), dimension(3), intent(inout) :: Dbarycoords

    ! parameterlist
    type(t_parlist), intent(inout) :: rparlist

    ! Particles
    type(t_Particles), intent(inout) :: rParticles

    ! Number of the current particle
    integer, intent(inout) :: iPart
    
    ! Current element
    integer, intent(inout) :: currentElement

    ! pointer to array containing the elements adjacent to a vertex.
    !
    ! Handle to 
    !       p_IelementsAtVertex = array(1..*) of integer
    ! p_IelementsAtVertex ( p_IelementsAtVertexIdx(IVT)..p_IelementsAtVertexIdx(IVT+1)-1 )
    ! contains the number of the adjacent element in a vertex.
    ! This replaces the old KVEL array.
    !
    ! Note: For hanging vertices, this array contains only those
    ! elements which are 'corner adjacent' to a vertex (i.e. the 'smaller' elements).
    ! The 'big' elements adjacent to the edge which the hanging vertex
    ! is a midpoint of are not part of the vertex neighbourhood
    ! in this array.
    integer, dimension(:), pointer :: p_IelementsAtVertex


    ! Handle to 
    !       p_IelementsAtVertexIdx=array [1..NVT+1] of integer.
    ! Index array for p_IelementsAtVertex of length NVT+1 for describing the
    ! elements adjacent to a corner vertex. for vertex IVT, the array
    ! p_IelementsAtVertex contains the numbers of the elements around this
    ! vertex at indices 
    !     p_IelementsAtVertexIdx(IVT)..p_IelementsAtVertexIdx(IVT+1)-1.
    ! By subtracting
    !     p_IelementsAtVertexIdx(IVT+1)-p_IelementsAtVertexIdx(IVT)
    ! One can get the number of elements adjacent to a vertex IVT.
    integer, dimension(:), pointer ::  p_IelementsAtVertexIdx


    ! Pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! Pointer to vertices at each element
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Pointer to the vertex coordinates
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! Coordinates of the vertices of the actual element
    real(DP), dimension(2,4) :: Dvert_coord

    ! Local variables
    integer :: ivt, Vert, Elm, nVertex   
    logical :: binside

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation
  
    ! Get vertices at element
    call storage_getbase_int2D(&
        p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    ! Get coordinates of the vertices
    call storage_getbase_double2D(&
        p_rtriangulation%h_DvertexCoords, p_DvertexCoords)

    ! Get elements at vertex
    call storage_getbase_int(p_rtriangulation%h_IelementsAtVertex,&
        p_IelementsAtVertex)
    call storage_getbase_int(p_rtriangulation%h_IelementsAtVertexIdx,&
        p_IelementsAtVertexIdx)

    ! store current element
    currentElement = rParticles%p_element(iPart)

    ! store coordinates of the vertices
    do ivt=1,3

      Dvert_coord(1,ivt)= p_DvertexCoords(1,p_IverticesAtElement(ivt,currentElement))
      Dvert_coord(2,ivt)= p_DvertexCoords(2,p_IverticesAtElement(ivt,currentElement))

    end do

    ! Check if the particle is in the element
    call gaux_isInElement_tri2D(Dcurrpos(1),Dcurrpos(2),Dvert_coord,binside)
  
    ! If the position is not in the element, then search in neighoured elements
    if (binside == .false.) then
    
        ! Loop over the vertices of the element
        SearchVertex: do Vert = 1, 3												

            !Current vertex
            nVertex = p_IverticesAtElement(Vert, currentElement)

            ! Loop over the element containing to the vertex
            SearchElement: do Elm = 1, (p_IelementsAtVertexIdx(nVertex+1)-p_IelementsAtVertexIdx(nVertex))		
    							
                if (p_IelementsAtVertex(p_IelementsAtVertexIdx(nVertex)+Elm-1) == 0) then
                exit SearchElement
                end if

                currentElement = p_IelementsAtVertex(p_IelementsAtVertexIdx(nVertex)+Elm-1)

                ! Store coordinates of cornervertices
                Dvert_coord(1,1)= p_DvertexCoords(1,p_IverticesAtElement(1,currentElement))
                Dvert_coord(1,2)= p_DvertexCoords(1,p_IverticesAtElement(2,currentElement))
                Dvert_coord(1,3)= p_DvertexCoords(1,p_IverticesAtElement(3,currentElement))
                Dvert_coord(2,1)= p_DvertexCoords(2,p_IverticesAtElement(1,currentElement))
                Dvert_coord(2,2)= p_DvertexCoords(2,p_IverticesAtElement(2,currentElement))
                Dvert_coord(2,3)= p_DvertexCoords(2,p_IverticesAtElement(3,currentElement))

                ! Check if the particle is in the element
                call gaux_isInElement_tri2D(Dcurrpos(1),Dcurrpos(2),Dvert_coord,binside)

                ! If the particle is in the element, then exit loop
                if (binside==.true.) exit SearchVertex

            end do SearchElement

        end do SearchVertex
        
    end if 
    if (binside==.true.) then
        ! Store coordinates of cornervertices
        Dvert_coord(1,1)= p_DvertexCoords(1,p_IverticesAtElement(1,currentElement))
        Dvert_coord(1,2)= p_DvertexCoords(1,p_IverticesAtElement(2,currentElement))
        Dvert_coord(1,3)= p_DvertexCoords(1,p_IverticesAtElement(3,currentElement))
        Dvert_coord(2,1)= p_DvertexCoords(2,p_IverticesAtElement(1,currentElement))
        Dvert_coord(2,2)= p_DvertexCoords(2,p_IverticesAtElement(2,currentElement))
        Dvert_coord(2,3)= p_DvertexCoords(2,p_IverticesAtElement(3,currentElement))

        ! Get barycentric coordinates for the current position    
        call gaux_getBarycentricCoords_tri2D(Dvert_coord,Dcurrpos(1),Dcurrpos(2),&
                Dbarycoords(1),Dbarycoords(2),Dbarycoords(3))
 
     else
        call eulerlagrange_findelement(rparlist,p_rproblemLevel,rParticles,iPart)

        currentElement= rParticles%p_element(iPart)

        ! Store coordinates of cornervertices
        Dvert_coord(1,1)= p_DvertexCoords(1,p_IverticesAtElement(1,currentElement))
        Dvert_coord(1,2)= p_DvertexCoords(1,p_IverticesAtElement(2,currentElement))
        Dvert_coord(1,3)= p_DvertexCoords(1,p_IverticesAtElement(3,currentElement))
        Dvert_coord(2,1)= p_DvertexCoords(2,p_IverticesAtElement(1,currentElement))
        Dvert_coord(2,2)= p_DvertexCoords(2,p_IverticesAtElement(2,currentElement))
        Dvert_coord(2,3)= p_DvertexCoords(2,p_IverticesAtElement(3,currentElement))

        ! Check if the particle is in the element
        call gaux_isInElement_tri2D(Dcurrpos(1),Dcurrpos(2),Dvert_coord,binside)

        ! Check for particle-wall-collisions
        if (binside == .false.) then
            call eulerlagrange_partwallcollision(rparlist,p_rproblemLevel,rParticles,iPart)     
        end if
     
     end if
                 
  end subroutine eulerlagrange_getbarycoordelm


  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_partwallcollision(rparlist,p_rproblemLevel,rParticles,iPart)

!<description>
    ! This subroutine checks if there is/was a collision beetween the particle and the wall

!<input>
    ! parameterlist
    type(t_parlist), intent(inout) :: rparlist

    ! particles
    type(t_Particles), intent(inout) :: rParticles

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel
    
    ! current number of particle
    integer, intent(inout) :: iPart

    ! Pointer to nodal property array. 
    ! Handle to 
    !       p_InodalProperty=array [1..NVT+NMT+NAT] of integer.
    ! p_InodalProperty(i) defines for each vertex i=(1..NVT),
    ! each edge i=(NVT+1..NVT+NMT) and face i=NVT+NMT+1..NVT+NMT+NAT
    ! its function inside of the geometry.
    ! Generally said, the range of the p_InodalProperty-array 
    ! characterizes the type of the node (=vertex/edge):
    ! = 0    : The vertex/edge is an inner vertex/edge
    ! > 0    : The vertex/edge is a boundary vertex/edge on the real
    !           boundary. KNPR(.) defines the number of the boundary
    !           component.
    ! This is the old KNPR-array, slightly modified for edges and
    ! hanging nodes.
    !
    ! In case there are hanging nodes in the mesh, this array
    ! has a special meaning for all hanging vertices and all edges
    ! containing hanging vertices.  Values < 0 indicate hanging
    ! vertices at an edge. 
    ! Let iedgeC (NVT+1..NVT+NMT) be the number of
    ! a a 'full' edge containing the hanging vertex jvertex. 
    ! Let iedge be one of the sub-edges inside of edge iedgeC.
    ! Then there is:
    !   p_InodalProperty(jvertex) = -iedgeC
    !   p_InodalProperty(iedge)   = -iedgeC
    !   p_InodalProperty(iedgeC)  = -jvertex
    ! Let kfaceC (NVT+NMT+1..NVT+NMT+NAT) be the number of a 'full' face
    ! containing the hanging vertex jvertex. 
    ! Let kface be the number of a one of the subfaces inside
    ! the face kfaceC. Let iedge be the number of one of the sub-edges
    ! inside face kfaceC.
    ! Then there is:
    !   p_InodalProperty(jvertex) = -kfaceC
    !   p_InodalProperty(kface)   = -kfaceC
    !   p_InodalProperty(iedge)   = -kfaceC
    !   p_InodalProperty(kfaceC)  = -jvertex
    ! A hanging vertex is either the midpoint of a face or of an edge,
    ! therefore this assignment is unique due to the range of the number.
    ! 'Hanging edges' (only appear in 3D) without a hanging vertex
    ! in the center of an edge/face are not supported.
    integer, dimension(:), pointer :: p_InodalProperty
 
    ! Pointer to edges Adjacent to an Element.
    ! Handle to 
    !       p_IedgesAtElement = array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the edges following the
    ! corner vertices in mathematically positive sense.
    ! This is the old KMID array.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IedgesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! To be able to distinguish a number of an edge from a vertex number, 
    ! edges are numbered in the range NVT+1..NVT+NMT. 
    integer, dimension(:,:), pointer :: p_IedgesAtElement

    ! Pointer to vertices Adjacent to an Edge. 
    ! Handle to 
    !       p_IverticesAtEdge = array [1..2,1..NMT]
    ! The numbers of the two vertices adjacent to an edge IMT. 
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
   
    ! pointer to vertices on boundary. 
    !
    ! Handle to 
    !       p_IverticesAtBoundary = array [1..NVBD] of integer.
    ! This array contains a list of all vertices on the (real) boundary
    ! in mathematically positive sense.
    ! The boundary vertices of boundary component i are saved at
    !        p_IboundaryCpIdx(i)..p_IboundaryCpIdx(i+1)-1.
    ! This is the old KVBD array.
    integer, dimension(:), pointer :: p_IverticesAtBoundary

    ! pointer to edges adjacent to the boundary. 
    !
    ! Handle to 
    !       p_IedgesAtBoundary = array [1..NMBD] of integer.
    ! This array contains a list of all edges on the (real) boundary.
    ! 2D: in mathematically positive sense. 
    ! 3D: with increasing number.
    ! The boundary edges of boundary component i are saved at
    !        p_IboundaryCpEdgesIdx(i)..p_IboundaryCpEdgesIdx(i+1)-1.
    ! This is the old KMBD array.
    ! (Note: In 2D, the above index pointer coincides with
    !        p_IboundaryCpEdgesIdx(i)..p_IboundaryCpEdgesIdx(i+1)-1 ).
    integer, dimension(:), pointer :: p_IedgesAtBoundary

    ! pointer to elements adjacent to the boundary. 
    !
    ! Handle to 
    !       p_IelementsAtBoundary = array [1..NVBD] of integer.
    ! This array contains a list of all elements on the (real) boundary
    ! in mathematically positive sense.
    ! p_IelementsAtBoundary(i) is the element adjacent to edge
    ! h_IedgesAtBoundary - therefore one element number might appear
    ! more than once in this array!
    ! The boundary elements of boundary component i are saved at
    !        p_IboundaryCpIdx(i)..p_IboundaryCpIdx(i+1)-1.
    ! This is the old KEBD array.
    integer, dimension(:), pointer :: p_IelementsAtBoundary

    ! pointer to the coordinates of the vertices
    !
    ! A list of all corner(!)-vertices of the elements in the triangulation.
    ! Handle to 
    !       p_RcornerCoordinates = array [1..ndim,1..NVT] of double
    ! with
    !   p_DvertexCoords(1,.) = X-coordinate.
    ! for 1D meshes,
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    ! for 2D meshes and
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    !   p_DvertexCoords(3,.) = Z-coordinate.
    ! for 3D meshes.
    ! This is a handle to the old DCORVG-array.
    !
    ! Note that the array may be longer than NVT in general!
    ! (May happen in case of a mesh hierarchy generated by a 2-level
    ! refinement, where the coordinates of the points on the
    ! coarser levels are contained in te coordinates of the
    ! finer levels.)
    ! In such a case, only the first NVT n-tuples in this array are valid!
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! local variables	
	real(DP), dimension(2) :: velo_rest
	real(DP) :: proj_tang, proj_norm
	real(DP), dimension(2) :: bdy_move

    integer :: i, iedge

    ! Local varibales for the calculation of the intersection point
    ! Coordinates of the two rays
    ! (x1,y1)->(x2,y2) and (x3,y3)->(x4,y4).
    real(DP) :: dx0,dy0,dx1,dy1,dx2,dy2,dx3,dy3
    ! Coordinates of the intersection point
    real(DP) :: dx,dy
    ! Check if the two rays intersect
    ! =-1: The rays are the same
    ! = 0: The rays do not intersect.
    ! = 1: The rays intersect in exactly one point.
    integer :: iintersect
   ! Parameter value of the intersection.
   ! The intersection point (dx,dy) can be found at position
   ! (dx,dy) = (dx0,dy0) + da*(dx1-dx0,dy1-dy0).
    real(DP) :: da

    ! Tangent and normal for collision with the boundary
    real(DP), dimension(2) :: tang, norm
    
    ! Timestep
    real(DP) :: dt

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation

    call storage_getbase_int(&
         p_rtriangulation%h_IelementsAtBoundary, p_IelementsAtBoundary)
    call storage_getbase_double2D(&
         p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
    call storage_getbase_int2D(&
         p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int(&
         p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)
    call storage_getbase_int(&
         p_rtriangulation%h_Inodalproperty, p_InodalProperty)
    call storage_getbase_int2D(&
         p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call storage_getbase_int2D(&
         p_rtriangulation%h_IverticesAtEdge, p_IverticesAtEdge)


    ! Find the edge on the boundary
    edgesearch: do i=1, 3
        iedge= p_IedgesAtElement(i,rParticles%p_element(iPart))
        if (p_InodalProperty(p_rtriangulation%NVT+iedge)>0) exit edgesearch
    end do edgesearch

    ! Set ray for the vertices of the edge
    dx2= p_DvertexCoords(1,p_IverticesAtEdge(1,iedge))
    dy2= p_DvertexCoords(2,p_IverticesAtEdge(1,iedge))
    dx3= p_DvertexCoords(1,p_IverticesAtEdge(2,iedge))
    dy3= p_DvertexCoords(2,p_IverticesAtEdge(2,iedge))

    ! Set ray for the particle movement from the old and the new particle position
    dx0= rParticles%p_xpos_old(iPart)
    dy0= rParticles%p_ypos_old(iPart)
    dx1= rParticles%p_xpos(iPart)
    dy1= rParticles%p_ypos(iPart)

    ! Calculate the intersection point of the rays (collision point with the boundary)
    call gaux_getIntersection_ray2D(&
          dx0,dy0,dx1,dy1,dx2,dy2,dx3,dy3, dx,dy, iintersect, da)

    if (iintersect == 1) then

        ! Get values for the startingpositions of the particles
        call parlst_getvalue_double(rparlist, 'Timestepping', "dinitialStep", dt)

        ! Set the new particle position on the boundary
        rParticles%p_xpos(iPart) = dx
        rParticles%p_ypos(iPart) = dy

        ! Movementvector to the boundary
        velo_rest(1)= rParticles%p_xpos_old(iPart) - dx
        velo_rest(2)= rParticles%p_ypos_old(iPart) - dy
        ! Norming the vector
        if (sqrt(velo_rest(1)**2.0_dp+velo_rest(2)**2.0_dp) .ne. 0) then
		    velo_rest(1)= velo_rest(1)/sqrt(velo_rest(1)**2.0_dp+velo_rest(2)**2.0_dp)
		    velo_rest(2)= velo_rest(2)/sqrt(velo_rest(1)**2.0_dp+velo_rest(2)**2.0_dp)
        end if

		! Tangent in x-direction 
		tang(1)= p_DvertexCoords(1,p_IverticesAtEdge(2,iedge)) - p_DvertexCoords(1,p_IverticesAtEdge(2,iedge))
		! Tangent in y-direction
		tang(2)= p_DvertexCoords(2,p_IverticesAtEdge(2,iedge)) - p_DvertexCoords(1,p_IverticesAtEdge(2,iedge))
        ! Norming of the tangent
		tang(1)= tang(1)/sqrt(tang(1)**2.0_dp+tang(2)**2.0_dp)
		tang(2)= tang(2)/sqrt(tang(1)**2.0_dp+tang(2)**2.0_dp)
        ! Normal in x-direction
        norm(1)= tang(2)
        ! Normal in y-direction
        norm(2)= -tang(1)

        ! Projection to the boundary
        proj_tang = velo_rest(1)*tang(1)+velo_rest(2)*tang(2)
        proj_norm = velo_rest(1)*norm(1)+velo_rest(2)*norm(2)

	    ! Calculate the rest-vector of the movement
	    bdy_move(1)= rParticles%tang_val*proj_tang*tang(1)+&
					 rParticles%norm_val*proj_norm*norm(1)
	    bdy_move(2)= rParticles%tang_val*proj_tang*tang(2)+&
					 rParticles%norm_val*proj_norm*norm(2)

		! Compute the new coordinates after the collision
		rParticles%p_xvelo(iPart) = -bdy_move(1)*sqrt(rParticles%p_xvelo(iPart)**2+rParticles%p_yvelo(iPart)**2)
		rParticles%p_yvelo(iPart) = bdy_move(2)*sqrt(rParticles%p_xvelo(iPart)**2+rParticles%p_yvelo(iPart)**2)

        ! Set the new "old" velocity
        rParticles%p_xvelo_old(iPart)= rParticles%p_xvelo(iPart)
        rParticles%p_yvelo_old(iPart)= rParticles%p_yvelo(iPart)

        ! Set new "old position" and "old velocity" of the particle
        rParticles%p_xpos_old(iPart)= dx
        rParticles%p_ypos_old(iPart)= dy

        ! Set new position of the particle
        rParticles%p_xpos(iPart) = rParticles%p_xpos_old(iPart) + dt*(1-da) * rParticles%p_xvelo(iPart)
        rParticles%p_ypos(iPart) = rParticles%p_ypos_old(iPart) + dt*(1-da) * rParticles%p_yvelo(iPart)

    else
        ! Set new "old position" and "old velocity" of the particle
        rParticles%p_xpos(iPart)= dx
        rParticles%p_ypos(iPart)= dy
        rParticles%p_xvelo(iPart)= 0.0_dp
        rParticles%p_yvelo(iPart)= 0.0_dp
        rParticles%p_xvelo_gas(iPart)= 0.0_dp
        rParticles%p_yvelo_gas(iPart)= 0.0_dp
    end if
	
    if (rParticles%p_xpos(iPart) .le. minval(p_DvertexCoords(1,:))) rParticles%p_xpos(iPart) = minval(p_DvertexCoords(1,:))
    if (rParticles%p_xpos(iPart) .ge. maxval(p_DvertexCoords(1,:))) rParticles%p_xpos(iPart) = maxval(p_DvertexCoords(1,:))
    if (rParticles%p_ypos(iPart) .le. minval(p_DvertexCoords(2,:))) rParticles%p_ypos(iPart) = minval(p_DvertexCoords(2,:))
    if (rParticles%p_ypos(iPart) .ge. maxval(p_DvertexCoords(2,:))) rParticles%p_ypos(iPart) = maxval(p_DvertexCoords(2,:))

    ! Find the new start element for the particle
    call eulerlagrange_findelement(rparlist,p_rproblemLevel,rParticles,iPart)

    ! Calculate barycentric coordinates
    call eulerlagrange_calcbarycoords(p_rproblemLevel,rParticles,iPart)

    ! Wrong element
    if ((abs(rParticles%p_lambda1(iPart))+abs(rParticles%p_lambda2(iPart))+&
              abs(rParticles%p_lambda3(iPart))-1) .ge. 0.00001) then
        call eulerlagrange_wrongelement(rparlist,p_rproblemLevel,rParticles,iPart)
    end if


  end subroutine eulerlagrange_partwallcollision

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_calcvolpart(p_rproblemLevel,rParticles)

!<description>
    ! This subroutine computes volumefraction of the particles in the gridpoints

!<input>
    ! Particles
    type(t_Particles), intent(inout) :: rParticles

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! pointer to array of the volume for each element
    !
    ! 2D triangulation: Array with area of each element.
    ! 3D triangulation: Array with volume of each element.
    ! Handle to 
    !       p_DelementArea = array [1..NEL+1] of double.
    ! p_DelementArea [NEL+1] gives the total area/voloume of the domain.
    real(DP), dimension(:), pointer :: p_DelementVolume

    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! Local variables
	integer :: i, current, ivt
	real(DP) :: c_pi

    ! Clean array for new computation of the volume fraction 
    rParticles%p_PartVol= 0
    
    ! Rise the counter for new volume fraction computation
    rParticles%iPartVolCount= rParticles%iPartVolCount+1

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation
 
    ! Get vertices at element
    call storage_getbase_int2D(&
         p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    ! Get area of each element
    call storage_getbase_double(&
         p_rtriangulation%h_DelementVolume, p_DelementVolume)

    c_pi= 3.14159265358979323846264338327950288_dp

	do i= 1, rParticles%nPart

        ! Element of the current particle
		current= rParticles%p_element(i)

        ! Calculate barycentric coordinates
        call eulerlagrange_calcbarycoords(p_rproblemLevel,rParticles,i)

        ! Store the volumefraction of the particle in the gridpoints (with barycentric coordinates)
		rParticles%p_PartVol(p_IverticesAtElement(1,current))= &
		                rParticles%p_PartVol(p_IverticesAtElement(1,current)) + &
		                (abs(rParticles%p_lambda1(i))*c_pi*0.25_dp*rParticles%p_diam(i)**2.0_dp)/&
		                p_DelementVolume(current)
		rParticles%p_PartVol(p_IverticesAtElement(2,current))= &
		                rParticles%p_PartVol(p_IverticesAtElement(2,current)) + &
		                (abs(rParticles%p_lambda2(i))*c_pi*0.25_dp*rParticles%p_diam(i)**2.0_dp)/&
		                p_DelementVolume(current)
		rParticles%p_PartVol(p_IverticesAtElement(3,current))= &
		                rParticles%p_PartVol(p_IverticesAtElement(3,current)) + &
		                (abs(rParticles%p_lambda3(i))*c_pi*0.25_dp*rParticles%p_diam(i)**2.0_dp)/&
		                p_DelementVolume(current)

	end do

    ! Store the volume fraction for the average volume fraction
    do ivt= 1, p_rtriangulation%NVT
        rParticles%p_PartVolAver(ivt)= rParticles%p_PartVolAver(ivt) + rParticles%p_PartVol(ivt)
    end do


  end subroutine eulerlagrange_calcvolpart

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_calcvelopart(p_rproblemLevel,rParticles)

!<description>
    ! This subroutine computes velocity of the particles in the gridpoints

!<input>
    ! particles
    type(t_Particles), intent(inout) :: rParticles

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    
    ! pointer to array of the volume for each element
    !
    ! 2D triangulation: Array with area of each element.
    ! 3D triangulation: Array with volume of each element.
    ! Handle to 
    !       p_DelementArea = array [1..NEL+1] of double.
    ! p_DelementArea [NEL+1] gives the total area/voloume of the domain.
    real(DP), dimension(:), pointer :: p_DelementVolume

    ! Pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! Local variables
	integer :: i, current
	real(DP) :: c_pi

    rParticles%p_PartVelox=0.0_dp

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation
 
    ! Get vertices at element
    call storage_getbase_int2D(&
         p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    ! Get area of each element
    call storage_getbase_double(&
         p_rtriangulation%h_DelementVolume, p_DelementVolume)

	do i= 1, rParticles%nPart

        ! Element of the current particle
		current= rParticles%p_element(i)

        ! Store the velocity of the particle in the gridpoints (with barycentric coordinates)
		rParticles%p_PartVelox(p_IverticesAtElement(1,current))= &
		                rParticles%p_PartVelox(p_IverticesAtElement(1,current)) + &
		                (rParticles%p_xvelo(i))/&
		                p_DelementVolume(p_IverticesAtElement(1,current))
		rParticles%p_PartVelox(p_IverticesAtElement(2,current))= &
		                rParticles%p_PartVelox(p_IverticesAtElement(2,current)) + &
		                (rParticles%p_xvelo(i))/&
		                p_DelementVolume(current)
		rParticles%p_PartVelox(p_IverticesAtElement(3,current))= &
		                rParticles%p_PartVelox(p_IverticesAtElement(3,current)) + &
		                (rParticles%p_xvelo(i))/&
		                p_DelementVolume(p_IverticesAtElement(3,current))
		rParticles%p_PartVeloy(p_IverticesAtElement(1,current))= &
		                rParticles%p_PartVeloy(p_IverticesAtElement(1,current)) + &
		                (rParticles%p_yvelo(i))/&
		                p_DelementVolume(current)
		rParticles%p_PartVeloy(p_IverticesAtElement(2,current))= &
		                rParticles%p_PartVeloy(p_IverticesAtElement(2,current)) + &
		                (rParticles%p_yvelo(i))/&
		                p_DelementVolume(p_IverticesAtElement(2,current))
		rParticles%p_PartVeloy(p_IverticesAtElement(3,current))= &
		                rParticles%p_PartVeloy(p_IverticesAtElement(3,current)) + &
		                (rParticles%p_yvelo(i))/&
		                p_DelementVolume(current)
	end do

  end subroutine eulerlagrange_calcvelopart


end module eulerlagrange_callback2d
