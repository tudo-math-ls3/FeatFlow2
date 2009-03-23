!##############################################################################
!# ****************************************************************************
!# <name> euler_callback2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible Euler/Navier-Stokes equations in 2D
!#
!# The following callback functions are available:
!#
!# 1.) euler_calcFluxGalerkin2d
!#     -> compute inviscid fluxes for standard Galerkin scheme
!#
!# 2.) euler_calcFluxTVD2d
!#     -> compute inviscid fluxes for TVD scheme
!#
!# 3.) euler_calcFluxScalarDiss2d
!#     -> compute inviscid fluxes for low-order discretization
!#        adopting scalar artificial viscosities
!#
!# 4.) euler_calcFluxTensorDiss2d
!#     -> compute inviscid fluxes for low-order discretization
!#        adopting tensorial artificial viscosities
!#
!# 5.) euler_calcFluxRusanov2d
!#     -> compute inviscid fluxes for low-order discretization
!#        adopting the Rusanov artificial diffusion
!#
!# 6.) euler_calcMatrixGalerkinDiag2d
!#     -> compute local matrices for standard Galerkin scheme
!#
!# 7.) euler_calcMatrixGalerkin2d
!#     -> compute local matrices for standard Galerkin scheme
!#
!# 8.) euler_calcMatrixScalarDissDiag2d
!#     -> compute local matrices for low-order discretization
!#        adopting scalar artificial viscosities
!#
!# 9.) euler_calcMatrixScalarDiss2d
!#     -> compute local matrices for low-order discretization
!#        adopting scalar artificial viscosities
!#
!# 10.) euler_calcMatrixTensorDissDiag2d
!#      -> compute local matrices for low-order discretization
!#         adopting tensorial artificial viscosities
!#
!# 11.) euler_calcMatrixTensorDiss2d
!#      -> compute local matrices for low-order discretization
!#         adopting tensorial artificial viscosities
!#
!# 12.) euler_calcMatrixRusanovDiag2d
!#      -> compute local matrices for low-order discretization
!#         adopting the Rusanov artificial viscosities
!#
!# 13.) euler_calcMatrixRusanov2d
!#      -> compute local matrices for low-order discretization
!#         adopting the Rusanov flux artificial viscosities
!#
!# 14.) euler_calcCharacteristics2d
!#      -> compute characteristic variables in 2D
!#
!# 15.) euler_calcBoundaryvalues2d
!#      -> compute the boundary values for a given node
!#
!# 16.) euler_hadaptCallbackScalar2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in interleave format
!#
!# 17.) euler_hadaptCallbackBlock2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in block format
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
  public :: euler_calcFluxTVD2d
  public :: euler_calcFluxScalarDiss2d
  public :: euler_calcFluxTensorDiss2d
  public :: euler_calcFluxRusanov2d
  public :: euler_calcMatrixGalerkinDiag2d
  public :: euler_calcMatrixGalerkin2d
  public :: euler_calcMatrixScalarDissDiag2d
  public :: euler_calcMatrixScalarDiss2d
  public :: euler_calcMatrixTensorDissDiag2d
  public :: euler_calcMatrixTensorDiss2d
  public :: euler_calcMatrixRusanovDiag2d
  public :: euler_calcMatrixRusanov2d
  public :: euler_calcCharacteristics2d
  public :: euler_calcBoundaryvalues2d
  public :: euler_hadaptCallbackScalar2d
  public :: euler_hadaptCallbackBlock2d

contains
  
  !*****************************************************************************
  
!<subroutine>

  subroutine euler_calcFluxGalerkin2d(U_i, U_j, C_ij, C_ji, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the standard
    ! Galerkin discretization in 2D without resorting to matrix-vector.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(IN) :: dscale
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(OUT) :: F_ij, F_ji
!</output>
!</subroutine>
    
    ! local variables
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
    real(DP) :: aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8

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
    !      / U2                             \
    ! F1 = | G1*U4-G14*aux1-G2*aux2         |
    !      | U3*aux3                        |
    !      \ (gamma*U4-G2*(aux1+aux2))*aux3 /
    !
    !      / U3                             \
    ! F2 = | U3*aux3 = U2*aux4              |
    !      | G1*U4-G14*aux2-G2*aux1         |
    !      \ (gamma*U4-G2*(aux1+aux2))*aux4 /
    !
    ! where the auxiliary values are defined as follows:
    !
    ! aux1 = U2*U2/U1 = aux3*U2
    ! aux2 = U3*U3/U1 = aux4*U3
    ! aux3 = U2/U1
    ! aux4 = U3/U1
    !
    ! and the predefined constants are given by:
    !
    ! G14 = (gamma-3)/2   and   G2 = (gamma-1)/2   and   G1 = gamma-1
    !
    ! The auxiliary quantities aux[1-4] and aux[5-8] are used 
    ! for the values in nodes i and j, respectively.
    ! ---------------------------------------------------------------------------

    ! Compute auxiliary values
    aux3 = U_i(2)/U_i(1);  aux4 = U_i(3)/U_i(1)
    aux1 = aux3*U_i(2);    aux2 = aux4*U_i(3)
    aux7 = U_j(2)/U_j(1);  aux8 = U_j(3)/U_j(1)
    aux5 = aux7*U_j(2);    aux6 = aux8*U_j(3)
    
    ! Compute fluxes for x-direction
    dF1_ij(1) = U_i(2)                             - U_j(2)
    dF1_ij(2) = G1*U_i(4)-G14*aux1-G2*aux2         - (G1*U_j(4)-G14*aux5-G2*aux6)
    dF1_ij(3) = U_i(3)*aux3                        - U_j(3)*aux7
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux3 - ((GAMMA*U_j(4)-G2*(aux5+aux6))*aux7)

    ! Compute fluxes for y-direction
    dF2_ij(1) = U_i(3)                             - U_j(3)
    dF2_ij(2) = U_i(3)*aux3                        - U_j(3)*aux7
    dF2_ij(3) = G1*U_i(4)-G14*aux2-G2*aux1         - (G1*U_j(4)-G14*aux6-G2*aux5)
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux4 - (GAMMA*U_j(4)-G2*(aux5+aux6))*aux8

    ! Assembly fluxes
    F_ij = dscale * ( C_ij(1)*dF1_ij + C_ij(2)*dF2_ij)
    F_ji = dscale * (-C_ji(1)*dF1_ij - C_ji(2)*dF2_ij)

  end subroutine euler_calcFluxGalerkin2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcFluxTVD2d(U_i, U_j, C_ij, C_ji, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the TVD
    ! discretization in 2D. The symmetric boundary contributions
    ! are neglected and incorporated in the antidiffusive flux.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(IN) :: dscale
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(OUT) :: F_ij, F_ji
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8

    !---------------------------------------------------------------------------
    ! Evaluate the Galerkin fluxes
    ! For a detailed description of algorithm and the definition of auxiliary
    ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------
    
    ! Compute auxiliary values
    aux3 = U_i(2)/U_i(1);  aux4 = U_i(3)/U_i(1)
    aux1 = aux3*U_i(2);    aux2 = aux4*U_i(3)
    aux7 = U_j(2)/U_j(1);  aux8 = U_j(3)/U_j(1)
    aux5 = aux7*U_j(2);    aux6 = aux8*U_j(3)
    
    ! Compute fluxes for x-direction
    dF1_ij(1) = U_i(2)                             -  U_j(2)
    dF1_ij(2) = G1*U_i(4)-G14*aux1-G2*aux2         - (G1*U_j(4)-G14*aux5-G2*aux6)
    dF1_ij(3) = U_i(3)*aux3                        - U_j(3)*aux7
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux3 - ((GAMMA*U_j(4)-G2*(aux5+aux6))*aux7)

    ! Compute fluxes for y-direction
    dF2_ij(1) = U_i(3)                             - U_j(3)
    dF2_ij(2) = U_i(3)*aux3                        - U_j(3)*aux7
    dF2_ij(3) = G1*U_i(4)-G14*aux2-G2*aux1         - (G1*U_j(4)-G14*aux6-G2*aux5)
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux4 - (GAMMA*U_j(4)-G2*(aux5+aux6))*aux8

    ! Compute skew-symmetric coefficient
    a = dscale * 0.5_DP*(C_ij-C_ji)

    ! Assembly fluxes and exploit skew-symmetry of a_ij and F_ij
    F_ij = a(1)*dF1_ij + a(2)*dF2_ij
    F_ji = F_ij

  end subroutine euler_calcFluxTVD2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcFluxScalarDiss2d(U_i, U_j, C_ij, C_ji, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using scalar dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(IN) :: dscale
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(OUT) :: F_ij, F_ji
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
    real(DP) :: d_ij,hi,hj,H_ij,q_ij,u_ij,v_ij,aux

    !---------------------------------------------------------------------------
    ! Evaluate the Galerkin fluxes
    ! For a detailed description of algorithm and the definition of auxiliary
    ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------
    
    ! Compute auxiliary values
    aux3 = U_i(2)/U_i(1);  aux4 = U_i(3)/U_i(1)
    aux1 = aux3*U_i(2);    aux2 = aux4*U_i(3)
    aux7 = U_j(2)/U_j(1);  aux8 = U_j(3)/U_j(1)
    aux5 = aux7*U_j(2);    aux6 = aux8*U_j(3)
    
    ! Compute fluxes for x-direction
    dF1_ij(1) = U_i(2)                             - U_j(2)
    dF1_ij(2) = G1*U_i(4)-G14*aux1-G2*aux2         - (G1*U_j(4)-G14*aux5-G2*aux6)
    dF1_ij(3) = U_i(3)*aux3                        - U_j(3)*aux7
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux3 - ((GAMMA*U_j(4)-G2*(aux5+aux6))*aux7)

    ! Compute fluxes for y-direction
    dF2_ij(1) = U_i(3)                             - U_j(3)
    dF2_ij(2) = U_i(3)*aux3                        - U_j(3)*aux7
    dF2_ij(3) = G1*U_i(4)-G14*aux2-G2*aux1         - (G1*U_j(4)-G14*aux6-G2*aux5)
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux4 - (GAMMA*U_j(4)-G2*(aux5+aux6))*aux8

    ! Assembly fluxes
    F_ij = dscale * ( C_ij(1)*dF1_ij + C_ij(2)*dF2_ij)
    F_ji = dscale * (-C_ji(1)*dF1_ij - C_ji(2)*dF2_ij)


    !---------------------------------------------------------------------------
    ! Evaluate the scalar dissipation 
    !---------------------------------------------------------------------------

    ! Compute skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji)

    ! Compute Roe mean values
    aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
    u_ij = (aux*aux3+aux7)/(aux+1.0_DP)
    v_ij = (aux*aux4+aux8)/(aux+1.0_DP)
    q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
    hi   = GAMMA*U_i(4)/U_i(1)-G2*(U_i(2)*U_i(2)+U_i(3)*U_i(3))/(U_i(1)*U_i(1))
    hj   = GAMMA*U_j(4)/U_j(1)-G2*(U_j(2)*U_j(2)+U_j(3)*U_j(3))/(U_j(1)*U_j(1))
    H_ij = (aux*hi+hj)/(aux+1.0_DP)

    ! Compute auxiliary variables
    aux  = sqrt(a(1)*a(1)+a(2)*a(2))
    aux1 = u_ij*a(1) + v_ij*a(2)
    aux2 = sqrt(max(-G1*(q_ij-H_ij), SYS_EPSREAL))

    ! Scalar dissipation
    d_ij = dscale * (abs(aux1) + aux2*aux)

    ! Multiply the solution difference by the artificial diffusion factor
    Diff = d_ij*(U_j-U_i)

    ! Add the artificial diffusion to the fluxes
    F_ij = F_ij+Diff
    F_ji = F_ji-Diff

  end subroutine euler_calcFluxScalarDiss2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcFluxTensorDiss2d(U_i, U_j, C_ij, C_ji, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using tensorial dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(IN) :: dscale
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(OUT) :: F_ij, F_ji
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(NVAR2D) :: dF1_ij,dF2_ij,Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
    real(DP) :: aux,uPow2,vPow2,hi,hj,H_ij,q_ij,u_ij,v_ij
    real(DP) :: anorm,l1,l2,l3,l4,w1,w2,w3,w4,cPow2,cs_ij


    !---------------------------------------------------------------------------
    ! Evaluate the Galerkin fluxes
    ! For a detailed description of algorithm and the definition of auxiliary
    ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------
    
    ! Compute auxiliary values
    aux3 = U_i(2)/U_i(1);  aux4 = U_i(3)/U_i(1)
    aux1 = aux3*U_i(2);    aux2 = aux4*U_i(3)
    aux7 = U_j(2)/U_j(1);  aux8 = U_j(3)/U_j(1)
    aux5 = aux7*U_j(2);    aux6 = aux8*U_j(3)
    
    ! Compute fluxes for x-direction
    dF1_ij(1) = U_i(2)                             - U_j(2)
    dF1_ij(2) = G1*U_i(4)-G14*aux1-G2*aux2         - (G1*U_j(4)-G14*aux5-G2*aux6)
    dF1_ij(3) = U_i(3)*aux3                        - U_j(3)*aux7
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux3 - ((GAMMA*U_j(4)-G2*(aux5+aux6))*aux7)

    ! Compute fluxes for y-direction
    dF2_ij(1) = U_i(3)                             - U_j(3)
    dF2_ij(2) = U_i(3)*aux3                        - U_j(3)*aux7
    dF2_ij(3) = G1*U_i(4)-G14*aux2-G2*aux1         - (G1*U_j(4)-G14*aux6-G2*aux5)
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux4 - (GAMMA*U_j(4)-G2*(aux5+aux6))*aux8

    ! Assembly fluxes
    F_ij = dscale * ( C_ij(1)*dF1_ij + C_ij(2)*dF2_ij)
    F_ji = dscale * (-C_ji(1)*dF1_ij - C_ji(2)*dF2_ij)

    
    !---------------------------------------------------------------------------
    ! Evaluate the tensorial dissipation 
    !---------------------------------------------------------------------------
    
    ! Compute the skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji); anorm = sqrt(a(1)*a(1)+a(2)*a(2))

    if (anorm .gt. SYS_EPSREAL) then

      ! Normalize the skew-symmetric coefficient
      a = a/anorm
      
      ! Compute Roe mean values
      aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
      u_ij = (aux*aux3+aux7)/(aux+1.0_DP)
      v_ij = (aux*aux4+aux8)/(aux+1.0_DP)
      hi   = GAMMA*U_i(4)/U_i(1)-G2*(U_i(2)*U_i(2)+U_i(3)*U_i(3))/(U_i(1)*U_i(1))
      hj   = GAMMA*U_j(4)/U_j(1)-G2*(U_j(2)*U_j(2)+U_j(3)*U_j(3))/(U_j(1)*U_j(1))
      H_ij = (aux*hi+hj)/(aux+1.0_DP)
      
      ! Compute auxiliary variables
      aux   = u_ij*a(1) + v_ij*a(2)
      uPow2 = u_ij*u_ij
      vPow2 = v_ij*v_ij
      q_ij  = 0.5_DP*(uPow2+vPow2)
      cPow2 = max(-G1*(q_ij-H_ij), SYS_EPSREAL)
      cs_ij = sqrt(cPow2)
      
      ! Compute eigenvalues
      l1 = abs(aux-cs_ij)
      l2 = abs(aux)
      l3 = abs(aux+cs_ij)
      l4 = abs(aux)

      ! Compute solution difference U_i-U_j
      Diff = U_j-U_i
      
      ! Compute auxiliary quantities for characteristic variables
      aux1 = G2/cPow2*(q_ij*Diff(1)-u_ij*Diff(2)-v_ij*Diff(3)+Diff(4))
      aux2 = 0.5_DP*(aux*Diff(1)-a(1)*Diff(2)-a(2)*Diff(3))/cs_ij

      ! Compute characteristic variables multiplied by the corresponding eigenvalue
      w1 = l1 * (aux1 + aux2)
      w2 = l2 * ((1.0_DP-G1*q_ij/cPow2)*Diff(1)+G1*(u_ij*Diff(2)+v_ij*Diff(3)-Diff(4))/cPow2)
      w3 = l3 * (aux1 - aux2)
      w4 = l4 * ((a(1)*v_ij-a(2)*u_ij)*Diff(1)+a(2)*Diff(2)-a(1)*Diff(3))

      ! Compute "R_ij * |Lbd_ij| * L_ij * dU"
      Diff(1) = dscale * anorm * ( w1 + w2 + w3 )
      Diff(2) = dscale * anorm * ( (u_ij-cs_ij*a(1))*w1 + u_ij*w2 + (u_ij+cs_ij*a(1))*w3 + a(2)*w4 )
      Diff(3) = dscale * anorm * ( (v_ij-cs_ij*a(2))*w1 + v_ij*w2 + (v_ij+cs_ij*a(2))*w3 - a(1)*w4 )
      Diff(4) = dscale * anorm * ( (H_ij-cs_ij*aux)*w1  + q_ij*w2 + (H_ij+cs_ij*aux)*w3  + (u_ij*a(2)-v_ij*a(1))*w4 )
      
      ! Add the artificial diffusion to the fluxes
      F_ij = F_ij+Diff
      F_ji = F_ji-Diff

    end if

  end subroutine euler_calcFluxTensorDiss2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcFluxRusanov2d(U_i, U_j, C_ij, C_ji, dscale, F_ij, F_ji)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 2D using the Rusanov dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji

    ! scaling coefficient
    real(DP), intent(IN) :: dscale
!</input>

!<output>
    ! inviscid fluxes
    real(DP), dimension(:), intent(OUT) :: F_ij, F_ji
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij
    real(DP), dimension(NVAR2D) :: Diff
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
    real(DP) :: d_ij,hi,hj,H_ij,q_ij,u_ij,v_ij,aux

    !---------------------------------------------------------------------------
    ! Evaluate the Galerkin fluxes
    ! For a detailed description of algorithm and the definition of auxiliary
    ! quantities have a look at the subroutine "euler_calcFluxGalerkin2d".
    !---------------------------------------------------------------------------
    
    ! Compute auxiliary values
    aux3 = U_i(2)/U_i(1);  aux4 = U_i(3)/U_i(1)
    aux1 = aux3*U_i(2);    aux2 = aux4*U_i(3)
    aux7 = U_j(2)/U_j(1);  aux8 = U_j(3)/U_j(1)
    aux5 = aux7*U_j(2);    aux6 = aux8*U_j(3)
    
    ! Compute fluxes for x-direction
    dF1_ij(1) = U_i(2)                             - U_j(2)
    dF1_ij(2) = G1*U_i(4)-G14*aux1-G2*aux2         - (G1*U_j(4)-G14*aux5-G2*aux6)
    dF1_ij(3) = U_i(3)*aux3                        - U_j(3)*aux7
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux3 - ((GAMMA*U_j(4)-G2*(aux5+aux6))*aux7)

    ! Compute fluxes for y-direction
    dF2_ij(1) = U_i(3)                             - U_j(3)
    dF2_ij(2) = U_i(3)*aux3                        - U_j(3)*aux7
    dF2_ij(3) = G1*U_i(4)-G14*aux2-G2*aux1         - (G1*U_j(4)-G14*aux6-G2*aux5)
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux4 - (GAMMA*U_j(4)-G2*(aux5+aux6))*aux8

    ! Assembly fluxes
    F_ij = dscale * ( C_ij(1)*dF1_ij + C_ij(2)*dF2_ij)
    F_ji = dscale * (-C_ji(1)*dF1_ij - C_ji(2)*dF2_ij)


    !---------------------------------------------------------------------------
    ! Evaluate the scalar dissipation 
    !---------------------------------------------------------------------------

    ! Compute skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji)

    ! Compute enthalpy
    hi = GAMMA*U_i(4)/U_i(1)-G2*(U_i(2)*U_i(2)+U_i(3)*U_i(3))/(U_i(1)*U_i(1))
    hj = GAMMA*U_j(4)/U_j(1)-G2*(U_j(2)*U_j(2)+U_j(3)*U_j(3))/(U_j(1)*U_j(1))

    ! Compute auxiliary variables
    aux  = sqrt(a(1)*a(1)+a(2)*a(2))
    aux1 = aux3*a(1) + aux4*a(2)
    aux2 = aux7*a(1) + aux8*a(2)
    aux5 = sqrt(max(-G1*(0.5_DP*(aux3*aux3+aux4*aux4)-hi), SYS_EPSREAL))
    aux6 = sqrt(max(-G1*(0.5_DP*(aux7*aux7+aux8*aux8)-hj), SYS_EPSREAL))
    
    ! Scalar dissipation for the Rusanov flux
    d_ij = dscale*max( abs(aux1) + aux5*aux, abs(aux2) + aux6*aux)

    ! Multiply the solution difference by the artificial diffusion factor
    Diff = d_ij*(U_j-U_i)

    ! Add the artificial diffusion to the fluxes
    F_ij = F_ij+Diff
    F_ji = F_ji-Diff

  end subroutine euler_calcFluxRusanov2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixGalerkinDiag2d(U_i, U_j, C_ij, C_ji, dscale, K_ij, K_ji)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(IN) :: dscale
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(OUT) :: K_ij, K_ji
!</output>
!</subroutine>

    ! local variable
    real(DP) :: ui,uj,vi,vj

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1)

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale*(G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij(3) = dscale*(uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(4) = dscale*(GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale*(G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji(3) = dscale*(ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(4) = dscale*(GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

  end subroutine euler_calcMatrixGalerkinDiag2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixGalerkin2d(U_i, U_j, C_ij, C_ji, dscale, K_ij, K_ji)

!<description>
    ! This subroutine computes the Galerkin matrices in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(IN) :: dscale
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(OUT) :: K_ij,K_ji
!</output>
!</subroutine>

    ! local variable
    real(DP) :: Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj,uPow2i,uPow2j,vPow2i,vPow2j,aux1,aux2

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1);   Ei = U_i(4)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1);   Ej = U_j(4)/U_j(1)

    uvi = ui*vi;   qi = ui*ui+vi*vi;   uPow2i = ui*ui;   uPow2i = vi*vi
    uvj = uj*vj;   qj = uj*uj+vj*vj;   uPow2j = uj*uj;   vPow2j = vj*vj

    aux1 = uj*C_ij(1)+vj*C_ij(2)
    aux2 = ui*C_ji(1)+vi*C_ji(2)
    
    ! Compute Galerkin coefficient K_ij
    K_ij( 1) = 0.0_DP
    K_ij( 2) = dscale*((G2*qj-uPow2j)*C_ij(1)-uvj*C_ij(2))
    K_ij( 3) = dscale*((G2*qj-vPow2j)*C_ij(2)-uvj*C_ij(1))
    K_ij( 4) = dscale*(G1*qj-GAMMA*Ej)*aux1

    K_ij( 5) = dscale*C_ij(1)
    K_ij( 6) = dscale*(G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij( 7) = dscale*(vj*C_ij(1)-G1*uj*C_ij(2))
    K_ij( 8) = dscale*((GAMMA*Ej-G2*qj)*C_ij(1)-G1*uj*aux1)

    K_ij( 9) = dscale*C_ij(2)
    K_ij(10) = dscale*(uj*C_ij(2)-G1*vj*C_ij(1))
    K_ij(11) = dscale*(uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(12) = dscale*((GAMMA*Ej-G2*qj)*C_ij(2)-G1*vj*aux1)

    K_ij(13) = 0.0_DP
    K_ij(14) = dscale*G1*C_ij(1)
    K_ij(15) = dscale*G1*C_ij(2)
    K_ij(16) = dscale*(GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji( 1) = 0.0_DP
    K_ji( 2) = dscale*((G1*qi-uPow2i)*C_ji(1)-uvi*C_ji(2))
    K_ji( 3) = dscale*((G1*qi-vPow2i)*C_ji(2)-uvi*C_ji(1))
    K_ji( 4) = dscale*(G1*qi-GAMMA*Ei)*aux2

    K_ji( 5) = dscale*C_ji(1)
    K_ji( 6) = dscale*(G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji( 7) = dscale*(vi*C_ji(1)-G1*ui*C_ji(2))
    K_ji( 8) = dscale*((GAMMA*Ei-G2*qi)*C_ji(1)-G1*ui*aux2)

    K_ji( 9) = dscale*C_ji(2)
    K_ji(10) = dscale*(ui*C_ji(2)-G1*vi*C_ji(1))
    K_ji(11) = dscale*(ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(12) = dscale*((GAMMA*Ei-G2*qi)*C_ji(2)-G1*vi*aux2)

    K_ji(13) = 0.0_DP
    K_ji(14) = dscale*G1*C_ji(1)
    K_ji(15) = dscale*G1*C_ji(2)
    K_ji(16) = dscale*(GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

  end subroutine euler_calcMatrixGalerkin2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixScalarDissDiag2d(U_i, U_j, C_ij, C_ji, dscale, K_ij, K_ji)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies scalar artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(IN) :: dscale
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(OUT) :: K_ij,K_ji
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: anorm,aux,aux1,hi,hj,H_ij,q_ij,ui,uj,u_ij,vi,vj,v_ij,d_ij

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1)

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale*(G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij(3) = dscale*(uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(4) = dscale*(GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale*(G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji(3) = dscale*(ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(4) = dscale*(GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

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
      q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
         
      ! Compute scalar dissipation
      d_ij = dscale*(abs(a(1)*u_ij+a(2)*v_ij) + anorm*sqrt(max(-G1*(q_ij-H_ij), SYS_EPSREAL)))
      
      ! Apply scalar dissipation
      K_ij = K_ij - d_ij
      K_ji = K_ji - d_ij

    end if

  end subroutine euler_calcMatrixScalarDissDiag2d

!*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixScalarDiss2d(U_i, U_j, C_ij, C_ji, dscale, K_ij, K_ji)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies scalar artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(IN) :: dscale
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(OUT) :: K_ij, K_ji
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: anorm,aux,hi,hj,H_ij,q_ij,u_ij,v_ij,d_ij
    real(DP) :: Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj,uPow2i,uPow2j,vPow2i,vPow2j,aux1,aux2

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1);   Ei = U_i(4)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1);   Ej = U_j(4)/U_j(1)

    uvi = ui*vi;   qi = ui*ui+vi*vi;   uPow2i = ui*ui;   uPow2i = vi*vi
    uvj = uj*vj;   qj = uj*uj+vj*vj;   uPow2j = uj*uj;   vPow2j = vj*vj

    aux1 = uj*C_ij(1)+vj*C_ij(2)
    aux2 = ui*C_ji(1)+vi*C_ji(2)
    
    ! Compute Galerkin coefficient K_ij
    K_ij( 1) = 0.0_DP
    K_ij( 2) = dscale*((G2*qj-uPow2j)*C_ij(1)-uvj*C_ij(2))
    K_ij( 3) = dscale*((G2*qj-vPow2j)*C_ij(2)-uvj*C_ij(1))
    K_ij( 4) = dscale*(G1*qj-GAMMA*Ej)*aux1

    K_ij( 5) = dscale*C_ij(1)
    K_ij( 6) = dscale*(G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij( 7) = dscale*(vj*C_ij(1)-G1*uj*C_ij(2))
    K_ij( 8) = dscale*((GAMMA*Ej-G2*qj)*C_ij(1)-G1*uj*aux1)

    K_ij( 9) = dscale*C_ij(2)
    K_ij(10) = dscale*(uj*C_ij(2)-G1*vj*C_ij(1))
    K_ij(11) = dscale*(uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(12) = dscale*((GAMMA*Ej-G2*qj)*C_ij(2)-G1*vj*aux1)

    K_ij(13) = 0.0_DP
    K_ij(14) = dscale*G1*C_ij(1)
    K_ij(15) = dscale*G1*C_ij(2)
    K_ij(16) = dscale*(GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji( 1) = 0.0_DP
    K_ji( 2) = dscale*((G1*qi-uPow2i)*C_ji(1)-uvi*C_ji(2))
    K_ji( 3) = dscale*((G1*qi-vPow2i)*C_ji(2)-uvi*C_ji(1))
    K_ji( 4) = dscale*(G1*qi-GAMMA*Ei)*aux2

    K_ji( 5) = dscale*C_ji(1)
    K_ji( 6) = dscale*(G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji( 7) = dscale*(vi*C_ji(1)-G1*ui*C_ji(2))
    K_ji( 8) = dscale*((GAMMA*Ei-G2*qi)*C_ji(1)-G1*ui*aux2)

    K_ji( 9) = dscale*C_ji(2)
    K_ji(10) = dscale*(ui*C_ji(2)-G1*vi*C_ji(1))
    K_ji(11) = dscale*(ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(12) = dscale*((GAMMA*Ei-G2*qi)*C_ji(2)-G1*vi*aux2)

    K_ji(13) = 0.0_DP
    K_ji(14) = dscale*G1*C_ji(1)
    K_ji(15) = dscale*G1*C_ji(2)
    K_ji(16) = dscale*(GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

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
      q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
         
      ! Compute scalar dissipation
      d_ij = dscale*(abs(a(1)*u_ij+a(2)*v_ij) + anorm*sqrt(max(-G1*(q_ij-H_ij), SYS_EPSREAL)))
      
      ! Apply scalar dissipation
      K_ij( 1) = K_ij( 1) - d_ij
      K_ij( 6) = K_ij( 6) - d_ij
      K_ij(11) = K_ij(11) - d_ij
      K_ij(16) = K_ij(16) - d_ij

      K_ji( 1) = K_ji( 1) - d_ij
      K_ji( 6) = K_ji( 6) - d_ij
      K_ji(11) = K_ji(11) - d_ij
      K_ji(16) = K_ji(16) - d_ij

    end if

  end subroutine euler_calcMatrixScalarDiss2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixTensorDissDiag2d(U_i, U_j, C_ij, C_ji, dscale, K_ij, K_ji)

!<description>
    ! This subroutine computes the diagonal of the local Roe matrices 
    ! and applies tensorial artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(IN) :: dscale
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(OUT) :: K_ij, K_ji
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NVAR2D,NVAR2D) :: R_ij,L_ij
    real(DP), dimension(NVAR2D) :: D_ij
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: aux,hi,hj,H_ij,q_ij,ui,uj,u_ij,vi,vj,v_ij
    real(DP) :: l1,l2,l3,l4,anorm,c1,c2,cs_ij,cPow2,vel

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1)

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale*(G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij(3) = dscale*(uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(4) = dscale*(GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale*(G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji(3) = dscale*(ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(4) = dscale*(GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

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
      q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)

      ! Compute auxiliary values
      c1    = a(1)/anorm
      c2    = a(2)/anorm
      vel   = c1*u_ij+c2*v_ij
      cPow2 = max(-G1*(q_ij-H_ij), SYS_EPSREAL)
      cs_ij = sqrt(cPow2)
      
      ! Diagonal matrix of eigenvalues
      l1 = abs(vel-cs_ij)
      l2 = abs(vel)
      l3 = abs(vel+cs_ij)
      l4 = abs(vel)
      
      ! Matrix of right eigenvectors
      R_ij(1,1) =  l1
      R_ij(2,1) =  l1*(u_ij-cs_ij*c1)
      R_ij(3,1) =  l1*(v_ij-cs_ij*c2)
      R_ij(4,1) =  l1*(H_ij-cs_ij*vel)
      
      R_ij(1,2) =  l2
      R_ij(2,2) =  l2*u_ij
      R_ij(3,2) =  l2*v_ij
      R_ij(4,2) =  l2*q_ij
      
      R_ij(1,3) =  l3
      R_ij(2,3) =  l3*(u_ij+cs_ij*c1)
      R_ij(3,3) =  l3*(v_ij+cs_ij*c2)
      R_ij(4,3) =  l3*(H_ij+cs_ij*vel)
      
      R_ij(1,4) =  0.0_DP
      R_ij(2,4) =  l4*c2
      R_ij(3,4) = -l4*c1
      R_ij(4,4) =  l4*(u_ij*c2-v_ij*c1)
      
      ! Matrix of left eigenvectors
      L_ij(1,1) = 0.5_DP*(G1*q_ij+cs_ij*vel)/cPow2
      L_ij(2,1) = (cPow2-G1*q_ij)/cPow2
      L_ij(3,1) = 0.5_DP*(G1*q_ij-cs_ij*vel)/cPow2
      L_ij(4,1) = v_ij*c1-u_ij*c2
      
      L_ij(1,2) = 0.5_DP*(-G1*u_ij-cs_ij*c1)/cPow2
      L_ij(2,2) = G1*u_ij/cPow2
      L_ij(3,2) = 0.5_DP*(-G1*u_ij+cs_ij*c1)/cPow2
      L_ij(4,2) = c2
      
      L_ij(1,3) = 0.5_DP*(-G1*v_ij-cs_ij*c2)/cPow2
      L_ij(2,3) = G1*v_ij/cPow2
      L_ij(3,3) = 0.5_DP*(-G1*v_ij+cs_ij*c2)/cPow2
      L_ij(4,3) = -c1
      
      L_ij(1,4) =  G2/cPow2
      L_ij(2,4) = -G1/cPow2
      L_ij(3,4) =  G2/cPow2
      L_ij(4,4) =  0.0_DP

      
      ! Compute tensorial dissipation D_ij = diag(R_ij*|Lbd_ij|*L_ij)*I
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

      ! Apply tensorial dissipation
      K_ij = K_ij - D_ij
      K_ji = K_ji - D_ij

    end if

  end subroutine euler_calcMatrixTensorDissDiag2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixTensorDiss2d(U_i, U_j, C_ij, C_ji, dscale, K_ij, K_ji)

!<description>
    ! This subroutine computes the diagonal of the local Roe matrices 
    ! and applies tensorial artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(IN) :: dscale
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(OUT) :: K_ij,K_ji
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NVAR2D,NVAR2D) :: R_ij,L_ij
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: anorm,aux,hi,hj,H_ij,q_ij,u_ij,v_ij,d_ij,vel,c1,c2,cPow2,cs_ij,l1,l2,l3,l4
    real(DP) :: Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj,uPow2i,uPow2j,vPow2i,vPow2j,aux1,aux2

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1);   Ei = U_i(4)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1);   Ej = U_j(4)/U_j(1)

    uvi = ui*vi;   qi = ui*ui+vi*vi;   uPow2i = ui*ui;   uPow2i = vi*vi
    uvj = uj*vj;   qj = uj*uj+vj*vj;   uPow2j = uj*uj;   vPow2j = vj*vj

    aux1 = uj*C_ij(1)+vj*C_ij(2)
    aux2 = ui*C_ji(1)+vi*C_ji(2)
    
    ! Compute Galerkin coefficient K_ij
    K_ij( 1) = 0.0_DP
    K_ij( 2) = dscale*((G2*qj-uPow2j)*C_ij(1)-uvj*C_ij(2))
    K_ij( 3) = dscale*((G2*qj-vPow2j)*C_ij(2)-uvj*C_ij(1))
    K_ij( 4) = dscale*(G1*qj-GAMMA*Ej)*aux1

    K_ij( 5) = dscale*C_ij(1)
    K_ij( 6) = dscale*(G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij( 7) = dscale*(vj*C_ij(1)-G1*uj*C_ij(2))
    K_ij( 8) = dscale*((GAMMA*Ej-G2*qj)*C_ij(1)-G1*uj*aux1)

    K_ij( 9) = dscale*C_ij(2)
    K_ij(10) = dscale*(uj*C_ij(2)-G1*vj*C_ij(1))
    K_ij(11) = dscale*(uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(12) = dscale*((GAMMA*Ej-G2*qj)*C_ij(2)-G1*vj*aux1)

    K_ij(13) = 0.0_DP
    K_ij(14) = dscale*G1*C_ij(1)
    K_ij(15) = dscale*G1*C_ij(2)
    K_ij(16) = dscale*(GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji( 1) = 0.0_DP
    K_ji( 2) = dscale*((G1*qi-uPow2i)*C_ji(1)-uvi*C_ji(2))
    K_ji( 3) = dscale*((G1*qi-vPow2i)*C_ji(2)-uvi*C_ji(1))
    K_ji( 4) = dscale*(G1*qi-GAMMA*Ei)*aux2

    K_ji( 5) = dscale*C_ji(1)
    K_ji( 6) = dscale*(G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji( 7) = dscale*(vi*C_ji(1)-G1*ui*C_ji(2))
    K_ji( 8) = dscale*((GAMMA*Ei-G2*qi)*C_ji(1)-G1*ui*aux2)

    K_ji( 9) = dscale*C_ji(2)
    K_ji(10) = dscale*(ui*C_ji(2)-G1*vi*C_ji(1))
    K_ji(11) = dscale*(ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(12) = dscale*((GAMMA*Ei-G2*qi)*C_ji(2)-G1*vi*aux2)

    K_ji(13) = 0.0_DP
    K_ji(14) = dscale*G1*C_ji(1)
    K_ji(15) = dscale*G1*C_ji(2)
    K_ji(16) = dscale*(GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

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
      q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)

      ! Compute auxiliary variables
      c1   = a(1)/anorm
      c2   = a(2)/anorm
      vel  = c1*u_ij+c2*v_ij

      ! Compute speed of sound
      cPow2 = max(-G1*(q_ij-H_ij), SYS_EPSREAL)
      cs_ij = sqrt(cPow2)
      
      ! Diagonal matrix of eigenvalues
      l1 = abs(vel-cs_ij)
      l2 = abs(vel)
      l3 = abs(vel+cs_ij)
      l4 = abs(vel)
      
      ! Matrix of right eigenvectors
      R_ij(1,1) =  l1
      R_ij(2,1) =  l1*(u_ij-cs_ij*c1)
      R_ij(3,1) =  l1*(v_ij-cs_ij*c2)
      R_ij(4,1) =  l1*(H_ij-cs_ij*vel)
      
      R_ij(1,2) =  l2
      R_ij(2,2) =  l2*u_ij
      R_ij(3,2) =  l2*v_ij
      R_ij(4,2) =  l2*q_ij
      
      R_ij(1,3) =  l3
      R_ij(2,3) =  l3*(u_ij+cs_ij*c1)
      R_ij(3,3) =  l3*(v_ij+cs_ij*c2)
      R_ij(4,3) =  l3*(H_ij+cs_ij*vel)
      
      R_ij(1,4) =  0.0_DP
      R_ij(2,4) =  l4*c2
      R_ij(3,4) = -l4*c1
      R_ij(4,4) =  l4*(u_ij*c2-v_ij*c1)
      
      ! Matrix of left eigenvectors
      L_ij(1,1) = 0.5_DP*(G1*q_ij+cs_ij*vel)/cPow2
      L_ij(2,1) = (cPow2-G1*q_ij)/cPow2
      L_ij(3,1) = 0.5_DP*(G1*q_ij-cs_ij*vel)/cPow2
      L_ij(4,1) = v_ij*c1-u_ij*c2

      L_ij(1,2) = 0.5_DP*(-G1*u_ij-cs_ij*c1)/cPow2
      L_ij(2,2) = G1*u_ij/cPow2
      L_ij(3,2) = 0.5_DP*(-G1*u_ij+cs_ij*c1)/cPow2
      L_ij(4,2) = c2

      L_ij(1,3) = 0.5_DP*(-G1*v_ij-cs_ij*c2)/cPow2
      L_ij(2,3) = G1*v_ij/cPow2
      L_ij(3,3) = 0.5_DP*(-G1*v_ij+cs_ij*c2)/cPow2
      L_ij(4,3) = -c1

      L_ij(1,4) =  G2/cPow2
      L_ij(2,4) = -G1/cPow2
      L_ij(3,4) =  G2/cPow2
      L_ij(4,4) =  0.0_DP
      
      ! Compute tensorial dissipation D_ij = R_ij*|Lbd_ij|*L_ij
      call DGEMM('n', 'n', NVAR2D, NVAR2D, NVAR2D, anorm,&
                 R_ij, NVAR2D, L_ij, NVAR2D, 0.0_DP, D_ij, NVAR2D)

      ! Apply tensorial dissipation
      K_ij = K_ij - D_ij
      K_ji = K_ji - D_ij
      
    end if

  end subroutine euler_calcMatrixTensorDiss2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixRusanovDiag2d(U_i, U_j, C_ij, C_ji, dscale, K_ij, K_ji)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(IN) :: dscale
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(OUT) :: K_ij,K_ji
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: anorm,aux1,aux2,hi,hj,ui,uj,vi,vj,d_ij

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1)

    ! Compute Galerkin coefficient K_ij
    K_ij(1) = 0.0_DP
    K_ij(2) = dscale*(G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij(3) = dscale*(uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(4) = dscale*(GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji(1) = 0.0_DP
    K_ji(2) = dscale*(G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji(3) = dscale*(ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(4) = dscale*(GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

    ! Compute skew-symmetric coefficient and its norm
    a = 0.5_DP*(C_ji-C_ij); anorm = sqrt(a(1)*a(1)+a(2)*a(2))

    if (anorm .gt. SYS_EPSREAL) then
      
      ! Compute enthalpy
      hi   = GAMMA*U_i(4)/U_i(1)-G2*(ui*ui+vi*vi)
      hj   = GAMMA*U_j(4)/U_j(1)-G2*(uj*uj+vj*vj)
      
      ! Compute auxiliary values
      aux1 = abs(a(1)*ui+a(2)*vi) + anorm*sqrt(max(-G1*(ui*ui+vi*vi-hi), SYS_EPSREAL))
      aux2 = abs(a(1)*uj+a(2)*vj) + anorm*sqrt(max(-G1*(uj*uj+vj*vj-hj), SYS_EPSREAL))

      ! Compute scalar dissipation
      d_ij = dscale*max(aux1,aux2)
      
      ! Apply scalar dissipation
      K_ij = K_ij - d_ij
      K_ji = K_ji - d_ij

    end if

  end subroutine euler_calcMatrixRusanovDiag2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixRusanov2d(U_i, U_j, C_ij, C_ji, dscale, K_ij, K_ji)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji

    ! scaling parameter
    real(DP), intent(IN) :: dscale
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(OUT) :: K_ij, K_ji
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: anorm,aux,hi,hj,d_ij
    real(DP) :: Ei,Ej,ui,uj,vi,vj,qi,qj,uvi,uvj,uPow2i,uPow2j,vPow2i,vPow2j,aux1,aux2

    ! Compute auxiliary variables
    ui = U_i(2)/U_i(1);   vi = U_i(3)/U_i(1);   Ei = U_i(4)/U_i(1)
    uj = U_j(2)/U_j(1);   vj = U_j(3)/U_j(1);   Ej = U_j(4)/U_j(1)

    uvi = ui*vi;   qi = ui*ui+vi*vi;   uPow2i = ui*ui;   uPow2i = vi*vi
    uvj = uj*vj;   qj = uj*uj+vj*vj;   uPow2j = uj*uj;   vPow2j = vj*vj

    aux1 = uj*C_ij(1)+vj*C_ij(2)
    aux2 = ui*C_ji(1)+vi*C_ji(2)
    
    ! Compute Galerkin coefficient K_ij
    K_ij( 1) = 0.0_DP
    K_ij( 2) = dscale*((G2*qj-uPow2j)*C_ij(1)-uvj*C_ij(2))
    K_ij( 3) = dscale*((G2*qj-vPow2j)*C_ij(2)-uvj*C_ij(1))
    K_ij( 4) = dscale*(G1*qj-GAMMA*Ej)*aux1

    K_ij( 5) = dscale*C_ij(1)
    K_ij( 6) = dscale*(G13*uj*C_ij(1)+vj*C_ij(2))
    K_ij( 7) = dscale*(vj*C_ij(1)-G1*uj*C_ij(2))
    K_ij( 8) = dscale*((GAMMA*Ej-G2*qj)*C_ij(1)-G1*uj*aux1)

    K_ij( 9) = dscale*C_ij(2)
    K_ij(10) = dscale*(uj*C_ij(2)-G1*vj*C_ij(1))
    K_ij(11) = dscale*(uj*C_ij(1)+G13*vj*C_ij(2))
    K_ij(12) = dscale*((GAMMA*Ej-G2*qj)*C_ij(2)-G1*vj*aux1)

    K_ij(13) = 0.0_DP
    K_ij(14) = dscale*G1*C_ij(1)
    K_ij(15) = dscale*G1*C_ij(2)
    K_ij(16) = dscale*(GAMMA*(uj*C_ij(1)+vj*C_ij(2)))

    ! Compute Galerkin coefficient K_ji
    K_ji( 1) = 0.0_DP
    K_ji( 2) = dscale*((G1*qi-uPow2i)*C_ji(1)-uvi*C_ji(2))
    K_ji( 3) = dscale*((G1*qi-vPow2i)*C_ji(2)-uvi*C_ji(1))
    K_ji( 4) = dscale*(G1*qi-GAMMA*Ei)*aux2

    K_ji( 5) = dscale*C_ji(1)
    K_ji( 6) = dscale*(G13*ui*C_ji(1)+vi*C_ji(2))
    K_ji( 7) = dscale*(vi*C_ji(1)-G1*ui*C_ji(2))
    K_ji( 8) = dscale*((GAMMA*Ei-G2*qi)*C_ji(1)-G1*ui*aux2)

    K_ji( 9) = dscale*C_ji(2)
    K_ji(10) = dscale*(ui*C_ji(2)-G1*vi*C_ji(1))
    K_ji(11) = dscale*(ui*C_ji(1)+G13*vi*C_ji(2))
    K_ji(12) = dscale*((GAMMA*Ei-G2*qi)*C_ji(2)-G1*vi*aux2)

    K_ji(13) = 0.0_DP
    K_ji(14) = dscale*G1*C_ji(1)
    K_ji(15) = dscale*G1*C_ji(2)
    K_ji(16) = dscale*(GAMMA*(ui*C_ji(1)+vi*C_ji(2)))

    ! Compute coefficients
    a = 0.5_DP*(C_ji-C_ij); anorm = sqrt(a(1)*a(1)+a(2)*a(2))

    if (anorm .gt. SYS_EPSREAL) then
      
      ! Compute auxiliary values
      aux1 = abs(a(1)*ui+a(2)*vi) + anorm*sqrt(max(-G1*(qi-hi), SYS_EPSREAL))
      aux2 = abs(a(1)*uj+a(2)*vj) + anorm*sqrt(max(-G1*(qj-hj), SYS_EPSREAL))
      
      ! Compute scalar dissipation
      d_ij = dscale*max(aux1, aux2)
      
      ! Apply scalar dissipation
      K_ij( 1) = K_ij( 1) - d_ij
      K_ij( 6) = K_ij( 6) - d_ij
      K_ij(11) = K_ij(11) - d_ij
      K_ij(16) = K_ij(16) - d_ij

      K_ji( 1) = K_ji( 1) - d_ij
      K_ji( 6) = K_ji( 6) - d_ij
      K_ji(11) = K_ji(11) - d_ij
      K_ji(16) = K_ji(16) - d_ij

    end if

  end subroutine euler_calcMatrixRusanov2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcCharacteristics2d(U_i, U_j, Dweight, W_ij, Lbd_ij, R_ij, L_ij)

!<description>
    ! This subroutine computes the characteristic variables in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! weighting vector
    real(DP), dimension(:), intent(IN) :: Dweight
!</input>

!<output>
    ! vector of characteristic variables
    real(DP), dimension(:), intent(OUT) :: W_ij

    ! OPTIONAL: diagonal matrix of eigenvalues
    real(DP), dimension(:), intent(OUT), optional :: Lbd_ij

    ! OPTIONAL: transformation matrix into conservative variables
    real(DP), dimension(:), intent(OUT), optional :: R_ij

    ! OPTIONAL: transformation matrix into characteristic variables
    real(DP), dimension(:), intent(OUT), optional :: L_ij
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(NVAR2D) :: Diff
    real(DP) :: u_ij,v_ij,H_ij,q_ij,cs_ij,aux,aux1,aux2,hi,hj,cPow2,uPow2,vPow2,c1,c2,cnorm    
    
    ! Compute norm of weighting coefficient
    cnorm = sqrt(Dweight(1)*Dweight(1)+Dweight(2)*Dweight(2))
    
    ! Check if weighting coefficient is zero
    if (cnorm .gt. SYS_EPSREAL) then

      ! Compute normalized weighting coefficient
      c1  = Dweight(1)/cnorm
      c2  = Dweight(2)/cnorm

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
      cPow2 = max(-G1*(q_ij-H_ij), SYS_EPSREAL)
      cs_ij = sqrt(cPow2)  
      aux   = c1*u_ij+c2*v_ij

      ! Compute diagonal matrix of eigenvalues (if present)
      if (present(Lbd_ij)) then
        Lbd_ij(1) = aux-cs_ij
        Lbd_ij(2) = aux
        Lbd_ij(3) = aux+cs_ij
        Lbd_ij(4) = aux
      end if

      ! Compute matrix of right eigenvectors
      if (present(R_ij)) then
        ! Matrix of right eigenvectors
        R_ij( 1) =  1.0_DP
        R_ij( 2) =  u_ij-cs_ij*c1
        R_ij( 3) =  v_ij-cs_ij*c2
        R_ij( 4) =  H_ij-cs_ij*aux

        R_ij( 5) =  1.0_DP
        R_ij( 6) =  u_ij
        R_ij( 7) =  v_ij
        R_ij( 8) =  q_ij

        R_ij( 9) =  1.0_DP
        R_ij(10) =  u_ij+cs_ij*c1
        R_ij(11) =  v_ij+cs_ij*c2
        R_ij(12) =  H_ij+cs_ij*aux
        R_ij(13) =  0.0_DP

        R_ij(14) =  c2
        R_ij(15) = -c1
        R_ij(16) =  u_ij*c2-v_ij*c1
      end if

      ! Compute matrix of left eigenvectors
      if (present(L_ij)) then
        L_ij( 1) =  0.5_DP*(G1*q_ij+cs_ij*aux)/cPow2
        L_ij( 2) = (cPow2-G1*q_ij)/cPow2
        L_ij( 3) =  0.5_DP*(G1*q_ij-cs_ij*aux)/cPow2
        L_ij( 4) =  v_ij*c1-u_ij*c2

        L_ij( 5) =  0.5_DP*(-G1*u_ij-cs_ij*c1)/cPow2
        L_ij( 6) =  G1*u_ij/cPow2
        L_ij( 7) =  0.5_DP*(-G1*u_ij+cs_ij*c1)/cPow2
        L_ij( 8) =  c2

        L_ij( 9) =  0.5_DP*(-G1*v_ij-cs_ij*c2)/cPow2
        L_ij(10) =  G1*v_ij/cPow2
        L_ij(11) =  0.5_DP*(-G1*v_ij+cs_ij*c2)/cPow2
        L_ij(12) = -c1

        L_ij(13) =  G2/cPow2
        L_ij(14) = -G1/cPow2
        L_ij(15) =  G2/cPow2
        L_ij(16) =  0.0_DP
      end if

      ! Compute solution difference U_i-U_j
      Diff = U_j-U_i

      ! Compute auxiliary quantities for characteristic variables
      aux1 = G2/cPow2*(q_ij*Diff(1)-u_ij*Diff(2)-v_ij*Diff(3)+Diff(4))
      aux2 = 0.5_DP*(aux*Diff(1)-c1*Diff(2)-c2*Diff(3))/cs_ij
      
      ! Compute characteristic variables
      W_ij(1) = cnorm * (aux1 + aux2)
      W_ij(2) = cnorm * ((1.0_DP-G1*q_ij/cPow2)*Diff(1)+G1*(u_ij*Diff(2)+v_ij*Diff(3)-Diff(4))/cPow2)
      W_ij(3) = cnorm * (aux1 - aux2)
      W_ij(4) = cnorm * ((c1*v_ij-c2*u_ij)*Diff(1)+c2*Diff(2)-c1*Diff(3))
      
    else   ! |dweight| = 0
      
      W_ij = 0.0_DP
      if (present(Lbd_ij)) Lbd_ij = 0.0_DP
      if (present(R_ij))   R_ij   = 0.0_DP
      if (present(L_ij))   L_ij   = 0.0_DP
      
    end if
    
  end subroutine euler_calcCharacteristics2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcBoundaryvalues2d(DbdrNormal, DpointNormal, DbdrValue,&
                                        ibdrCondType, Du, Du0, istatus)

!<description>
    ! This subroutine computes the boundary values for a given node in 2D
!</description>

!<input>
    ! normal vector at the boundary
    real(DP), dimension(:), intent(IN) :: DbdrNormal

    ! normal vector at the point on the boundary
    real(DP), dimension(:), intent(IN) :: DpointNormal

    ! evaluated boundary values
    real(DP), dimension(:), intent(IN) :: DbdrValue

    ! initial solution from the previous time step
    real(DP), dimension(:), intent(IN) :: Du0

    ! type of boundary condition
    integer, intent(IN) :: ibdrCondType
!</input>

!<inputoutput>
    ! computed boundary values
    real(DP), dimension(:), intent(INOUT) :: Du

    ! OPTIONAL: status of the callback function
    integer, intent(INOUT), optional :: istatus
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
      !    "High-order accurate implementation of solid wall
      !     boundary conditions in curved geometries"
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
      !    "Unstructured Grid Adaptive Algorithms for 
      !     Fluid Dynamics and Heat Conduction"
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

      ! Compute normal elocities at the boundary and the ghost state
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

      ! Check if Newton's method converged
      if (ite .ge. 100) then
        if (present(istatus)) then
          istatus = -ibdrCondType
          return
        else
          call output_line('Riemann solver failed due to divergence in&
                           & Newton-Raphson iteration',&
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
      p   = DbdrValue(4); 
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
      p   = DbdrValue(4); 
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
      ! "Adaptive Finite Element Solution Algorithm
      !  for the Euler Equations", R.A. Shapiro

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

  subroutine euler_hadaptCallbackScalar2d(rcollection, iOperation, Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D. The solution vector is assumed
    ! to be store in scalar interleave format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(IN) :: iOperation

    ! Array of vertices involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ivertices

    ! Array of elements involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ielements
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolution
    real(DP), dimension(:), pointer, save :: p_Dsolution
    integer :: ivar


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! This subroutine assumes that the name of the first solution
      ! vector is stored in the second quick access string.

      ! Retrieve solution vector from colletion and set pointer
      rsolution => collct_getvalue_vec(rcollection,&
                                        trim(rcollection%SquickAccess(2)))

      ! Check if solution is stored in interleave format
      if (rsolution%nblocks .ne. 1) then
        call output_line('Vector is not in interleave format!',&
                         OU_CLASS_WARNING,OU_MODE_STD,'euler_hadaptCallbackScalar2d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)
      
      
    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)
      

    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, NVAR2D*Ivertices(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, NVAR2D*Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      do ivar = 1, NVAR2D
        p_Dsolution((Ivertices(1)-1)*NVAR2D+ivar) = &
            0.5_DP*(p_Dsolution((Ivertices(2)-1)*NVAR2D+ivar)+&
                    p_Dsolution((Ivertices(3)-1)*NVAR2D+ivar))
      end do

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, NVAR2D*Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      do ivar = 1, NVAR2D
        p_Dsolution((Ivertices(1)-1)*NVAR2D+ivar) = &
            0.25_DP*(p_Dsolution((Ivertices(2)-1)*NVAR2D+ivar)+&
                     p_Dsolution((Ivertices(3)-1)*NVAR2D+ivar)+&
                     p_Dsolution((Ivertices(4)-1)*NVAR2D+ivar)+&
                     p_Dsolution((Ivertices(5)-1)*NVAR2D+ivar))
      end do

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (Ivertices(2) .ne. 0) then
        do ivar = 1, NVAR2D
          p_Dsolution((Ivertices(1)-1)*NVAR2D+ivar) = &
              p_Dsolution((Ivertices(2)-1)*NVAR2D+ivar)
        end do
      else
        do ivar = 1, NVAR2D
          p_Dsolution((Ivertices(1)-1)*NVAR2D+ivar) = 0.0_DP
        end do
      end if
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)

    end select
    
  end subroutine euler_hadaptCallbackScalar2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_hadaptCallbackBlock2d(rcollection, iOperation, Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D. The solution vector is assumed
    ! to be store in block format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(IN) :: iOperation

    ! Array of vertices involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ivertices

    ! Array of elements involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ielements
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolution
    real(DP), dimension(:), pointer, save :: p_Dsolution
    integer :: ivar,neq


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! This subroutine assumes that the name of the first solution
      ! vector is stored in the second quick access string.

      ! Retrieve solution vector from colletion and set pointer
      rsolution => collct_getvalue_vec(rcollection,&
                                        trim(rcollection%SquickAccess(2)))

      ! Check if solution is stored in interleave format
      if (rsolution%nblocks .ne. NVAR2D) then
        call output_line('Vector is not in block format!',&
                         OU_CLASS_WARNING,OU_MODE_STD,'euler_hadaptCallbackBlock2d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)
      
      
    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)
      

    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, NVAR2D*Ivertices(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, NVAR2D*Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      neq = rsolution%NEQ/NVAR2D
      do ivar = 1, NVAR2D
        p_Dsolution((ivar-1)*neq+Ivertices(1)) = &
            0.5_DP*(p_Dsolution((ivar-1)*neq+Ivertices(2))+&
                    p_Dsolution((ivar-1)*neq+Ivertices(3)) )
      end do

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)

      
    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolution, NVAR2D*Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      neq = rsolution%NEQ/NVAR2D
      do ivar = 1, NVAR2D
        p_Dsolution((ivar-1)*neq+Ivertices(1)) =&
            0.25_DP*(p_Dsolution((ivar-1)*neq+Ivertices(2))+&
                     p_Dsolution((ivar-1)*neq+Ivertices(3))+&
                     p_Dsolution((ivar-1)*neq+Ivertices(4))+&
                     p_Dsolution((ivar-1)*neq+Ivertices(5)) )
      end do

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (Ivertices(2) .ne. 0) then
        neq = rsolution%NEQ/NVAR2D
        do ivar = 1, NVAR2D
          p_Dsolution((ivar-1)*neq+Ivertices(1)) = &
              p_Dsolution((ivar-1)*neq+Ivertices(2))
        end do
      else
        neq = rsolution%NEQ/NVAR2D
        do ivar = 1, NVAR2D
          p_Dsolution((ivar-1)*neq+Ivertices(1)) = 0.0_DP
        end do
      end if
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)

    end select
    
  end subroutine euler_hadaptCallbackBlock2d
 
end module euler_callback2d
