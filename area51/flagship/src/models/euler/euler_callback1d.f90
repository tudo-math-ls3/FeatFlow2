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
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting scalar artificial viscosities
!#
!# 4.) euler_calcFluxTensorDiss1d
!#     -> Computes inviscid fluxes for low-order discretization
!#        adopting tensorial artificial viscosities
!#
!# 5.) euler_calcFluxRusanov1d
!#     -> Computes inviscid fluxes for low-order discretization
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
!#      -> Computes local matrices for low-order discretization
!#         adopting scalar artificial viscosities
!#
!# 11.) euler_calcMatrixScalarDiss1d
!#      -> Computes local matrices for low-order discretization
!#         adopting scalar artificial viscosities
!#
!# 12.) euler_calcMatrixTensorDissDiag1d
!#      -> Computes local matrices for low-order discretization
!#         adopting tensorial artificial viscosities
!#
!# 13.) euler_calcMatrixTensorDiss1d
!#      -> Computes local matrices for low-order discretization
!#         adopting tensorial artificial viscosities
!#
!# 14.) euler_calcMatrixRusanovDiag1d
!#      -> Computes local matrices for low-order discretization
!#         adopting the Rusanov artificial viscosities
!#
!# 15.) euler_calcMatrixRusanov1d
!#      -> Computes local matrices for low-order discretization
!#         adopting the Rusanov flux artificial viscosities
!#
!# 16.) euler_calcCharacteristics1d
!#      -> Computes characteristic variables in 1D
!#
!# 17.) euler_calcBoundaryvalues1d
!#      -> Computes the boundary values for a given node
!#
!# 18.) euler_hadaptCallbackScalar1d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 1D, whereby the vector is stored in interleave format
!#
!# 19.) euler_hadaptCallbackBlock1d
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
    ! Galerkin discretization in 1D.
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
    ! discretization in 1D. The symmetric boundary contributions
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
    a = dscale * 0.5_DP*(C_ij-C_ji)

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
    cs   = sqrt(max(-G1*(q_ij-H_ij), SYS_EPSREAL))

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
      cPow2 = max(-G1*(q_ij-H_ij), SYS_EPSREAL)
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
      !  corresponding eigenvalue
      w1 = l1 * 0.5_DP * ( (b1+u_ij/cs)*Diff(1) - (b2*u_ij+1.0_DP/cs)*Diff(2) + b2*Diff(3) )
      w2 = l2 *          ( (1-b1)*Diff(1)       +  b2*u_ij*Diff(2)            - b2*Diff(3) )
      w3 = l3 * 0.5_DP * ( (b1-u_ij/cs)*Diff(1) - (b2*u_ij-1.0_DP/cs)*Diff(2) + b2*Diff(3) )

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
                       anorm*sqrt(max(-G1*(q_ij-H_ij), SYS_EPSREAL)))
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
                      anorm*sqrt(max(-G1*(q_ij-H_ij), SYS_EPSREAL)))

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
      cPow2 = max(-G1*(q_ij-H_ij), SYS_EPSREAL)
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
      cPow2 = max(-G1*(q_ij-H_ij), SYS_EPSREAL)
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
    real(DP) :: ui,uj,ci,cj,Ei,Ej,uPow2i,uPow2j

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
    D_ij = dscale * max( abs(C_ij(1)*uj) + abs(C_ij(1))*cj,&
                         abs(C_ji(1)*ui) + abs(C_ji(1))*ci )

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
    
  end subroutine euler_calcCharacteristics1d

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
