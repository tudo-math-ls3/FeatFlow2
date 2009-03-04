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
!# 5.) euler_calcMatrixGalerkinDiag2d
!#     -> compute local matrices for standard Galerkin scheme
!#
!# 6.) euler_calcMatrixGalerkin2d
!#     -> compute local matrices for standard Galerkin scheme
!#
!# 7.) euler_calcMatrixScalarDissDiag2d
!#     -> compute local matrices for low-order discretization
!#        adopting scalar artificial viscosities
!#
!# 8.) euler_calcMatrixScalarDiss2d
!#     -> compute local matrices for low-order discretization
!#        adopting scalar artificial viscosities
!#
!# 9.) euler_calcMatrixTensorDissDiag2d
!#     -> compute local matrices for low-order discretization
!#        adopting tensorial artificial viscosities
!#
!# 10.) euler_calcMatrixTensorDiss2d
!#      -> compute local matrices for low-order discretization
!#         adopting tensorial artificial viscosities
!#
!# 11.) euler_calcCharacteristics2d
!#      -> compute characteristic variables in 2D
!#
!# 12.) euler_calcBoundaryvalues2d
!#      -> compute the boundary values for a given node
!#
!# 13.) euler_hadaptCallbackScalar2d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 2D, whereby the vector is stored in interleave format
!#
!# 14.) euler_hadaptCallbackBlock2d
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
  public :: euler_calcRawFluxFCT2d
  public :: euler_calcMatrixGalerkinDiag2d
  public :: euler_calcMatrixGalerkin2d
  public :: euler_calcMatrixScalarDissDiag2d
  public :: euler_calcMatrixScalarDiss2d
  public :: euler_calcMatrixTensorDissDiag2d
  public :: euler_calcMatrixTensorDiss2d
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
    ! F1 = | -G14*aux1-G2*(aux2+U4)         |
    !      | U3*aux3                        |
    !      \ (gamma*U4-G2*(aux1+aux2))*aux3 /
    !
    !      / U3                             \
    ! F2 = | U3*aux3 = U2*aux4              |
    !      | -G14*aux2-G2*(aux1+U4)         |
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
    ! G14 = (gamma-3)/2   and   G2 = (gamma-1)/2
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
    dF1_ij(1) = U_i(2)                             -  U_j(2)
    dF1_ij(2) = -G14*aux1-G2*(aux2+U_i(4))         - (-G14*aux5-G2*(aux6+U_j(4)))
    dF1_ij(3) = U_i(3)*aux3                        - U_j(3)*aux7
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux3 - ((GAMMA*U_j(4)-G2*(aux5+aux6))*aux7)

    ! Compute fluxes for y-direction
    dF2_ij(1) = U_i(3)                             - U_j(3)
    dF2_ij(2) = U_i(3)*aux3                        - U_j(3)*aux7
    dF2_ij(3) = -G14*aux2-G2*(aux1+U_i(4))         - (-G14*aux6-G2*(aux5+U_j(4)))
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
    dF1_ij(2) = -G14*aux1-G2*(aux2+U_i(4))         - (-G14*aux5-G2*(aux6+U_j(4)))
    dF1_ij(3) = U_i(3)*aux3                        - U_j(3)*aux7
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux3 - ((GAMMA*U_j(4)-G2*(aux5+aux6))*aux7)

    ! Compute fluxes for y-direction
    dF2_ij(1) = U_i(3)                             - U_j(3)
    dF2_ij(2) = U_i(3)*aux3                        - U_j(3)*aux7
    dF2_ij(3) = -G14*aux2-G2*(aux1+U_i(4))         - (-G14*aux6-G2*(aux5+U_j(4)))
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux4 - (GAMMA*U_j(4)-G2*(aux5+aux6))*aux8

    ! Compute skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji)

    ! Assembly fluxes and exploit skew-symmetry of a_ij and F_ij
    F_ij = dscale * (a(1)*dF1_ij + a(2)*dF2_ij)
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
    dF1_ij(1) = U_i(2)                             -  U_j(2)
    dF1_ij(2) = -G14*aux1-G2*(aux2+U_i(4))         - (-G14*aux5-G2*(aux6+U_j(4)))
    dF1_ij(3) = U_i(3)*aux3                        - U_j(3)*aux7
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux3 - ((GAMMA*U_j(4)-G2*(aux5+aux6))*aux7)

    ! Compute fluxes for y-direction
    dF2_ij(1) = U_i(3)                             - U_j(3)
    dF2_ij(2) = U_i(3)*aux3                        - U_j(3)*aux7
    dF2_ij(3) = -G14*aux2-G2*(aux1+U_i(4))         - (-G14*aux6-G2*(aux5+U_j(4)))
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux4 - (GAMMA*U_j(4)-G2*(aux5+aux6))*aux8

    ! Assembly fluxes
    F_ij = dscale * ( C_ij(1)*dF1_ij + C_ij(2)*dF2_ij)
    F_ji = dscale * (-C_ji(1)*dF1_ij - C_ji(2)*dF2_ij)


    !---------------------------------------------------------------------------
    ! Evaluate the scalar dissipation 
    !---------------------------------------------------------------------------

    ! Compute skew-symmetric coefficient
    a = 0.5_DP*(C_ij-C_ji)

#ifdef RUSANOV_FLUX

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
#else

    ! Compute Roe mean values
    aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
    u_ij = (aux*U_i(2)/U_i(1)+U_j(2)/U_j(1))/(aux+1.0_DP)
    v_ij = (aux*U_i(3)/U_i(1)+U_j(3)/U_j(1))/(aux+1.0_DP)
    q_ij = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
    hi   = GAMMA*U_i(4)/U_i(1)-G2*(U_i(2)*U_i(2)+U_i(3)*U_i(3))/(U_i(1)*U_i(1))
    hj   = GAMMA*U_j(4)/U_j(1)-G2*(U_j(2)*U_j(2)+U_j(3)*U_j(3))/(U_j(1)*U_j(1))
    H_ij = (aux*hi+hj)/(aux+1.0_DP)

    ! Compute auxiliary variables
    aux  = a(1)*a(1)+a(2)*a(2)
    aux1 = u_ij*a(1) + v_ij*a(2)
    aux2 = sqrt(max(-G1*(q_ij-H_ij), SYS_EPSREAL))

    ! Scalar dissipation
    d_ij = dscale*(abs(aux1) + aux2*aux)

#endif

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
    real(DP), dimension(NVAR2D,NVAR2D) :: Roe,L_ij,R_ij,D_ij
    real(DP), dimension(NVAR2D) :: dF1_ij, dF2_ij,Diff,Daux
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
    real(DP) :: aux,u2,v2,uv,hi,hj,H_ij,q_ij,u_ij,v_ij
    real(DP) :: l1,l2,l3,l4,vel,cnrm,cx,cy,c2,cs_ij


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
    dF1_ij(2) = -G14*aux1-G2*(aux2+U_i(4))         - (-G14*aux5-G2*(aux6+U_j(4)))
    dF1_ij(3) = U_i(3)*aux3                        - U_j(3)*aux7
    dF1_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux3 - ((GAMMA*U_j(4)-G2*(aux5+aux6))*aux7)

    ! Compute fluxes for y-direction
    dF2_ij(1) = U_i(3)                             - U_j(3)
    dF2_ij(2) = U_i(3)*aux3                        - U_j(3)*aux7
    dF2_ij(3) = -G14*aux2-G2*(aux1+U_i(4))         - (-G14*aux6-G2*(aux5+U_j(4)))
    dF2_ij(4) = (GAMMA*U_i(4)-G2*(aux1+aux2))*aux4 - (GAMMA*U_j(4)-G2*(aux5+aux6))*aux8

    ! Assembly fluxes
    F_ij = dscale * ( C_ij(1)*dF1_ij + C_ij(2)*dF2_ij)
    F_ji = dscale * (-C_ji(1)*dF1_ij - C_ji(2)*dF2_ij)

    
    !---------------------------------------------------------------------------
    ! Evaluate the tensorial dissipation 
    !---------------------------------------------------------------------------

    ! Compute solution difference
    Diff = U_i-U_j

    ! Compute Roe mean values
    aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
    u_ij = (aux*U_i(2)/U_i(1)+U_j(2)/U_j(1))/(aux+1.0_DP)
    v_ij = (aux*U_i(3)/U_i(1)+U_j(3)/U_j(1))/(aux+1.0_DP)
    hi   = GAMMA*U_i(4)/U_i(1)-G2*(U_i(2)*U_i(2)+U_i(3)*U_i(3))/(U_i(1)*U_i(1))
    hj   = GAMMA*U_j(4)/U_j(1)-G2*(U_j(2)*U_j(2)+U_j(3)*U_j(3))/(U_j(1)*U_j(1))
    H_ij = (aux*hi+hj)/(aux+1.0_DP)

    ! Compute coefficient
    a = 0.5_DP*(C_ij-C_ji)
    
    ! Compute auxiliary variables
    aux1  = u_ij*a(1) + v_ij*a(2)
    u2    = u_ij*u_ij
    v2    = v_ij*v_ij
    uv    = u_ij*v_ij
    q_ij  = 0.5_DP*(u2+v2)
   
    ! Characteristic velocity
    cnrm = sqrt(a(1)*a(1)+a(2)*a(2))
    if (cnrm .ne. 0.0_DP) then
      
      ! Compute auxiliary variables
      cx  = a(1)/cnrm
      cy  = a(2)/cnrm
      vel = cx*u_ij+cy*v_ij
      
      ! Compute speed of sound
      c2    = max(-G1*(q_ij-H_ij), SYS_EPSREAL)
      cs_ij = sqrt(c2)
      
      ! Diagonal matrix of eigenvalues
      l1 = abs(vel-cs_ij)
      l2 = abs(vel)
      l3 = abs(vel+cs_ij)
      l4 = abs(vel)
      
      ! Matrix of right eigenvectors multiplied by the diagonal
      ! matrix of real eigenvalues from the left |Lbd_{ij}|*R_{ij}
      R_ij(1,1) =  l1
      R_ij(2,1) =  l1*(u_ij-cs_ij*cx)
      R_ij(3,1) =  l1*(v_ij-cs_ij*cy)
      R_ij(4,1) =  l1*(H_ij-cs_ij*vel)
      
      R_ij(1,2) =  l2
      R_ij(2,2) =  l2*u_ij
      R_ij(3,2) =  l2*v_ij
      R_ij(4,2) =  l2*q_ij
      
      R_ij(1,3) =  l3
      R_ij(2,3) =  l3*(u_ij+cs_ij*cx)
      R_ij(3,3) =  l3*(v_ij+cs_ij*cy)
      R_ij(4,3) =  l3*(H_ij+cs_ij*vel)
      
      R_ij(1,4) =  0.0_DP
      R_ij(2,4) =  l4*cy
      R_ij(3,4) = -l4*cx
      R_ij(4,4) =  l4*(u_ij*cy-v_ij*cx)
      
      if (abs(cx) > abs(cy)) then

        ! Matrix of left eigenvectors if C(x) is largest coefficient
        L_ij(1,1) = 0.5_DP*(G1*q_ij+cs_ij*vel)/c2
        L_ij(2,1) = (c2-G1*q_ij)/c2
        L_ij(3,1) = 0.5_DP*(G1*q_ij-cs_ij*vel)/c2
        L_ij(4,1) = (v_ij-vel*cy)/cx
        
        L_ij(1,2) = 0.5_DP*(-G1*u_ij-cs_ij*cx)/c2
        L_ij(2,2) = G1*u_ij/c2
        L_ij(3,2) = 0.5_DP*(-G1*u_ij+cs_ij*cx)/c2
        L_ij(4,2) = cy
        
        L_ij(1,3) = 0.5_DP*(-G1*v_ij-cs_ij*cy)/c2
        L_ij(2,3) = G1*v_ij/c2
        L_ij(3,3) = 0.5_DP*(-G1*v_ij+cs_ij*cy)/c2
        L_ij(4,3) = (cy*cy-1)/cx
        
        L_ij(1,4) =  G2/c2
        L_ij(2,4) = -G1/c2
        L_ij(3,4) =  G2/c2
        L_ij(4,4) =  0.0_DP

      else

        ! Matrix of left eigenvectors if C(Y) is largest coefficient
        L_ij(1,1) = 0.5_DP*(G1*q_ij+cs_ij*vel)/c2
        L_ij(2,1) = (c2-G1*q_ij)/c2
        L_ij(3,1) = 0.5_DP*(G1*q_ij-cs_ij*vel)/c2
        L_ij(4,1) = (vel*cx-u_ij)/cy
        
        L_ij(1,2) = 0.5_DP*(-G1*u_ij-cs_ij*cx)/c2
        L_ij(2,2) = G1*u_ij/c2
        L_ij(3,2) = 0.5_DP*(-G1*u_ij+cs_ij*cx)/c2
        L_ij(4,2) = (1-cx*cx)/cy
        
        L_ij(1,3) = 0.5_DP*(-G1*v_ij-cs_ij*cy)/c2
        L_ij(2,3) = G1*v_ij/c2
        L_ij(3,3) = 0.5_DP*(-G1*v_ij+cs_ij*cy)/c2
        L_ij(4,3) = -cx
        
        L_ij(1,4) =  G2/c2
        L_ij(2,4) = -G1/c2
        L_ij(3,4) =  G2/c2
        L_ij(4,4) =  0.0_DP

      end if
      
      ! Compute D_ij = R_ij*|Lbd_ij|*L_ij
      call DGEMM('n', 'n', NVAR2D, NVAR2D, NVAR2D, cnrm,&
                     R_ij, NVAR2D, L_ij, NVAR2D, 0.0_DP, D_ij, NVAR2D)
      
      ! Compute D_ij*(u_i-u_j)
      call DGEMV('n', NVAR2D, NVAR2D, dscale, D_ij,&
                      NVAR2D, Diff, 1, 0.0_DP, Daux, 1)

      ! Add the artificial diffusion to the fluxes
      F_ji = F_ji+Daux
      F_ij = F_ij-Daux

    end if

  end subroutine euler_calcFluxTensorDiss2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcRawFluxFCT2d(U_i, U_j, C_ij, C_ji, dscale, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive flux
    ! for the FCT discretization in 2D.
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
    real(DP), dimension(:), intent(OUT) :: F_ij
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(NVAR2D,NVAR2D) :: Roe
    real(DP), dimension(NVAR2D) :: diff
    real(DP), dimension(NDIM2D) :: a,s
    real(DP) :: aux,aux1,aux2,u2,v2,uv,hi,hj,H_ij,q_ij,u_ij,v_ij,cs_ij,d_ij

    ! Compute solution difference
    diff = u_i-u_j
    
    ! Compute Roe mean values
    aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
    u_ij = (aux*U_i(2)/U_i(1)+U_j(2)/U_j(1))/(aux+1.0_DP)
    v_ij = (aux*U_i(3)/U_i(1)+U_j(3)/U_j(1))/(aux+1.0_DP)
    
    hi   = GAMMA*U_i(4)/U_i(1)-G2*(U_i(2)*U_i(2)+U_i(3)*U_i(3))/(U_i(1)*U_i(1))
    hj   = GAMMA*U_j(4)/U_j(1)-G2*(U_j(2)*U_j(2)+U_j(3)*U_j(3))/(U_j(1)*U_j(1))
    H_ij = (aux*hi+hj)/(aux+1.0_DP)
    
    ! Compute coefficient
    a = 0.5_DP*(C_ji-C_ij)
    s = 0.5_DP*(C_ij+C_ji)
    
    ! Compute auxiliary variables
    aux1  = u_ij*a(1) + v_ij*a(2)
    aux2  = u_ij*s(1) + v_ij*s(2)
    u2    = u_ij*u_ij
    v2    = v_ij*v_ij
    uv    = u_ij*v_ij
    q_ij  = 0.5_DP*(u2+v2)
    cs_ij = sqrt(max(-G1*(q_ij-H_ij), SYS_EPSREAL))
   
    ! Compute Roe matrix for the skew-symmetric part
    Roe(1,1) =   0.0_DP
    Roe(2,1) =   (G1*q_ij-u2)*s(1)   - uv*s(2)
    Roe(3,1) = - uv*s(1)             + (G1*q_ij-v2)*s(2)
    Roe(4,1) =   (G1*q_ij-H_ij)*aux2
    
    Roe(1,2) =   s(1)
    Roe(2,2) =   aux2-G6*u_ij*s(1)
    Roe(3,2) =   v_ij*s(1)           - G1*u_ij*s(2)
    Roe(4,2) =   (H_ij-G1*u2)*s(1)   - G1*uv*s(2)
    
    Roe(1,3) =                         s(2)
    Roe(2,3) = - G1*v_ij*s(1)        + u_ij*s(2)
    Roe(3,3) =                         aux2-G6*v_ij*s(2)
    Roe(4,3) = - G1*uv*s(1)          + (H_ij-G1*v2)*(2)

    Roe(1,4) =   0.0_DP
    Roe(2,4) =   G1*s(1)
    Roe(3,4) =                         G1*s(2)
    Roe(4,4) =   GAMMA*aux2

    ! Compute Roe*(u_j-u_i) for skew-symmetric part
    call DGEMV('n', NVAR2D, NVAR2D, dscale, Roe,&
                    NVAR2D, diff, 1, 0.0_DP, F_ij, 1)

    ! Scalar dissipation
    d_ij = dscale*(abs(aux1) + cs_ij*sqrt(a(1)*a(1)+a(2)*a(2)))
    
    ! Assemble antidiffusive flux
    F_ij = F_ij+d_ij*diff

  end subroutine euler_calcRawFluxFCT2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixGalerkinDiag2d(U_i, U_j, C_ij, C_ji, A_ij, S_ij)

!<description>
    ! This subroutine computes the diagonal of the local Roe matrices in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(OUT) :: A_ij,S_ij
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NDIM2D) :: a,s
    real(DP) :: aux,aux1,aux2,u_ij,v_ij

    ! Compute Roe mean values
    aux  =   sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
    u_ij = ( aux*U_i(2)/U_i(1)+U_j(2)/U_j(1) )/(aux+1.0_DP)
    v_ij = ( aux*U_i(3)/U_i(1)+U_j(3)/U_j(1) )/(aux+1.0_DP)
    
    ! Compute coefficients
    a = 0.5_DP*(C_ji-C_ij)
    s = 0.5_DP*(C_ij+C_ji)

    ! Compute auxiliary variables
    aux1 = u_ij*a(1) + v_ij*a(2)
    aux2 = u_ij*s(1) + v_ij*s(2)

    ! Compute diagonal Roe matrix for skew-symmetric part
    A_ij(1) = 0.0_DP
    A_ij(2) = aux1-G6*u_ij*a(1)
    A_ij(3) = aux1-G6*v_ij*a(2)
    A_ij(4) = GAMMA*aux1

    ! Compute diagonal Roe matrix for symmetric part
    S_ij(1) = 0.0_DP
    S_ij(2) = aux2-G6*u_ij*s(1)
    S_ij(3) = aux2-G6*v_ij*s(2)
    S_ij(4) = GAMMA*aux2

  end subroutine euler_calcMatrixGalerkinDiag2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixGalerkin2d(U_i, U_j, C_ij, C_ji, A_ij, S_ij)

!<description>
    ! This subroutine computes the local Roe matrices in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(OUT) :: A_ij,S_ij
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(NDIM2D) :: a,s
    real(DP) :: aux,aux1,aux2,u2,v2,uv,hi,hj,H_ij,q_ij,u_ij,v_ij
    
    ! Compute Roe mean values
    aux  =   sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
    u_ij = ( aux*U_i(2)/U_i(1)+U_j(2)/U_j(1) )/(aux+1.0_DP)
    v_ij = ( aux*U_i(3)/U_i(1)+U_j(3)/U_j(1) )/(aux+1.0_DP)
    hi   =   GAMMA*U_i(4)/U_i(1)-&
             G2*( U_i(2)*U_i(2)+U_i(3)*U_i(3) )/(U_i(1)*U_i(1))
    hj   =   GAMMA*U_j(4)/U_j(1)-&
             G2*( U_j(2)*U_j(2)+U_j(3)*U_j(3) )/(U_j(1)*U_j(1))
    H_ij = ( aux*hi+hj )/(aux+1.0_DP)

    ! Compute coefficients
    a = 0.5_DP*(C_ji-C_ij)
    s = 0.5_DP*(C_ij+C_ji)

    ! Compute auxiliary variables
    aux1 = u_ij*a(1) + v_ij*a(2)
    aux2 = u_ij*s(1) + v_ij*s(2)
    u2   = u_ij*u_ij
    v2   = v_ij*v_ij
    uv   = u_ij*v_ij
    q_ij = 0.5_DP*(u2+v2)
    
    ! Compute Roe matrix for skew-symmetric part
    A_ij( 1) =   0.0_DP
    A_ij( 2) =   (G1*q_ij-u2)*a(1)   - uv*a(2)
    A_ij( 3) = - uv*a(1)             + (G1*q_ij-v2)*a(2)
    A_ij( 4) =   (G1*q_ij-H_ij)*aux1
    
    A_ij( 5) =   a(1)
    A_ij( 6) =   aux1-G6*u_ij*a(1)
    A_ij( 7) =   v_ij*a(1)           - G1*u_ij*a(2)
    A_ij( 8) =   (H_ij-G1*u2)*a(1)   - G1*uv*a(2)
    
    A_ij( 9) =                         a(2)
    A_ij(10) = - G1*v_ij*a(1)        + u_ij*a(2)
    A_ij(11) =                         aux1-G6*v_ij*a(2)
    A_ij(12) = - G1*uv*a(1)          + (H_ij-G1*v2)*a(2)

    A_ij(13) =   0.0_DP
    A_ij(14) =   G1*a(1)
    A_ij(15) =                         G1*a(2)
    A_ij(16) =   GAMMA*aux1

    ! Check if boundary integral vanishes
    if (abs(s(1))+abs(s(2)) .gt. SYS_EPSREAL) then
      
      ! Compute Roe matrix for symmetric part
      S_ij( 1) =   0.0_DP
      S_ij( 2) =   (G1*q_ij-u2)*s(1)   - uv*s(2)
      S_ij( 3) = - uv*s(1)             + (G1*q_ij-v2)*s(2)
      S_ij( 4) =   (G1*q_ij-H_ij)*aux2
      
      S_ij( 5) =   s(1)
      S_ij( 6) =   aux2-G6*u_ij*s(1)
      S_ij( 7) =   v_ij*s(1)           - G1*u_ij*s(2)
      S_ij( 8) =   (H_ij-G1*u2)*s(1)   - G1*uv*s(2)
      
      S_ij( 9) =                         s(2)
      S_ij(10) = - G1*v_ij*s(1)        + u_ij*s(2)
      S_ij(11) =                         aux2-G6*v_ij*s(2)
      S_ij(12) = - G1*uv*s(1)          + (H_ij-G1*v2)*s(2)
      
      S_ij(13) =   0.0_DP
      S_ij(14) =   G1*s(1)
      S_ij(15) =                         G1*s(2)
      S_ij(16) =   GAMMA*aux2

    else
      
      S_ij = 0.0_DP

    end if

  end subroutine euler_calcMatrixGalerkin2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixScalarDissDiag2d(U_i, U_j, C_ij, C_ji, A_ij, S_ij)

!<description>
    ! This subroutine computes the diagonal of the local Roe matrices 
    ! and applies scalar artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(OUT) :: A_ij,S_ij
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: aux,aux1,hi,hj,H_ij,q_ij,u_ij,v_ij
    real(DP) :: cnrm,cs_ij,vel
    
    ! Compute Roe mean values
    aux  =   sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
    u_ij = ( aux*U_i(2)/U_i(1)+U_j(2)/U_j(1) )/(aux+1.0_DP)
    v_ij = ( aux*U_i(3)/U_i(1)+U_j(3)/U_j(1) )/(aux+1.0_DP)
    hi   =   GAMMA*U_i(4)/U_i(1)-&
             G2*( U_i(2)*U_i(2)+U_i(3)*U_i(3) )/(U_i(1)*U_i(1))
    hj   =   GAMMA*U_j(4)/U_j(1)-&
             G2*( U_j(2)*U_j(2)+U_j(3)*U_j(3) )/(U_j(1)*U_j(1))
    H_ij = ( aux*hi+hj )/(aux+1.0_DP)

    ! Compute coefficients
    a = 0.5_DP*(C_ji-C_ij)

    ! Compute auxiliary variables
    aux1  = u_ij*a(1) + v_ij*a(2)
    q_ij  = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)

    ! Compute diagonal Roe matrix for skew-symmetric part
    A_ij(1) = 0.0_DP
    A_ij(2) = aux1-G6*u_ij*a(1)
    A_ij(3) = aux1-G6*v_ij*a(2)
    A_ij(4) = GAMMA*aux1

    ! Compute auxiliary variables
    cnrm = sqrt(a(1)*a(1)+a(2)*a(2))
    if (cnrm .eq. 0.0_DP) then
      S_ij = 0.0_DP
      return
    end if
    
    ! Compute velocities
    cs_ij = sqrt(max(-G1*(q_ij-H_ij), SYS_EPSREAL))
    vel   = a(1)*u_ij+a(2)*v_ij

    ! Compute scalar viscosities
    S_ij = -abs(vel)-cnrm*cs_ij

  end subroutine euler_calcMatrixScalarDissDiag2d

!*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixScalarDiss2d(U_i, U_j, C_ij, C_ji, A_ij, S_ij)

!<description>
    ! This subroutine computes the local Roe matrices
    ! and applies scalar artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(OUT) :: A_ij,S_ij
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: aux,aux1,u2,v2,uv,hi,hj,H_ij,q_ij,u_ij,v_ij
    real(DP) :: cnrm,cs_ij,vel
    
    ! Compute Roe mean values
    aux  =   sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
    u_ij = ( aux*U_i(2)/U_i(1)+U_j(2)/U_j(1) )/(aux+1.0_DP)
    v_ij = ( aux*U_i(3)/U_i(1)+U_j(3)/U_j(1) )/(aux+1.0_DP)
    hi   =   GAMMA*U_i(4)/U_i(1)-&
             G2*( U_i(2)*U_i(2)+U_i(3)*U_i(3) )/(U_i(1)*U_i(1))
    hj   =   GAMMA*U_j(4)/U_j(1)-&
             G2*( U_j(2)*U_j(2)+U_j(3)*U_j(3) )/(U_j(1)*U_j(1))
    H_ij = ( aux*hi+hj )/(aux+1.0_DP)

    ! Compute coefficients
    a = 0.5_DP*(C_ji-C_ij)

    ! Compute auxiliary variables
    aux1 = u_ij*a(1) + v_ij*a(2)
    u2   = u_ij*u_ij
    v2   = v_ij*v_ij
    uv   = u_ij*v_ij
    q_ij = 0.5_DP*(u2+v2)
    
    ! Compute Roe matrix for skew-symmetric part
    A_ij( 1) =   0.0_DP
    A_ij( 2) =   (G1*q_ij-u2)*a(1)   - uv*a(2)
    A_ij( 3) = - uv*a(1)             + (G1*q_ij-v2)*a(2)
    A_ij( 4) =   (G1*q_ij-H_ij)*aux1
    
    A_ij( 5) =   a(1)
    A_ij( 6) =   aux1-G6*u_ij*a(1)
    A_ij( 7) =   v_ij*a(1)           - G1*u_ij*a(2)
    A_ij( 8) =   (H_ij-G1*u2)*a(1)   - G1*uv*a(2)
    
    A_ij( 9) =                         a(2)
    A_ij(10) = - G1*v_ij*a(1)        + u_ij*a(2)
    A_ij(11) =                         aux1-G6*v_ij*a(2)
    A_ij(12) = - G1*uv*a(1)          + (H_ij-G1*v2)*a(2)

    A_ij(13) =   0.0_DP
    A_ij(14) =   G1*a(1)
    A_ij(15) =                         G1*a(2)
    A_ij(16) =   GAMMA*aux1

    ! Compute auxiliary variables
    cs_ij = sqrt(max(-G1*(q_ij-H_ij), SYS_EPSREAL))
    cnrm  = sqrt(a(1)*a(1)+a(2)*a(2))
    vel   = a(1)*u_ij+a(2)*v_ij

    ! Compute scalar viscosities
    S_ij     = 0.0_DP
    S_ij( 1) = -abs(vel)-cnrm*cs_ij
    S_ij( 6) = -abs(vel)-cnrm*cs_ij
    S_ij(11) = -abs(vel)-cnrm*cs_ij
    S_ij(16) = -abs(vel)-cnrm*cs_ij

  end subroutine euler_calcMatrixScalarDiss2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixTensorDissDiag2d(U_i, U_j, C_ij, C_ji, A_ij, S_ij)

!<description>
    ! This subroutine computes the diagonal of the local Roe matrices 
    ! and applies tensorial artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(OUT) :: A_ij,S_ij
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NVAR2D,NVAR2D) :: R_ij,L_ij
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: aux,aux1,hi,hj,H_ij,q_ij,u_ij,v_ij
    real(DP) :: l1,l2,l3,l4,cnrm,cx,cy,cs_ij,c2,vel

    ! Compute Roe mean values
    aux  =   sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
    u_ij = ( aux*U_i(2)/U_i(1)+U_j(2)/U_j(1) )/(aux+1.0_DP)
    v_ij = ( aux*U_i(3)/U_i(1)+U_j(3)/U_j(1) )/(aux+1.0_DP)
    hi   =   GAMMA*U_i(4)/U_i(1)-&
             G2*( U_i(2)*U_i(2)+U_i(3)*U_i(3) )/(U_i(1)*U_i(1))
    hj   =   GAMMA*U_j(4)/U_j(1)-&
             G2*( U_j(2)*U_j(2)+U_j(3)*U_j(3) )/(U_j(1)*U_j(1))
    H_ij = ( aux*hi+hj )/(aux+1.0_DP)

    ! Compute coefficients
    a = 0.5_DP*(C_ji-C_ij)

    ! Compute auxiliary variables
    aux1  = u_ij*a(1) + v_ij*a(2)
    q_ij  = 0.5_DP*(u_ij*u_ij+v_ij*v_ij)
    
    ! Compute diagonal Roe matrix for skew-symmetric part
    A_ij(1) = 0.0_DP
    A_ij(2) = aux1-G6*u_ij*a(1)
    A_ij(3) = aux1-G6*v_ij*a(2)
    A_ij(4) = GAMMA*aux1

    ! Compute auxiliary variables
    cnrm = sqrt(a(1)*a(1)+a(2)*a(2))
    if (cnrm .eq. 0.0_DP) then
      S_ij = 0.0_DP
      return
    end if

    cx    = a(1)/cnrm
    cy    = a(2)/cnrm
    vel   = cx*u_ij+cy*v_ij
    c2    = max(-G1*(q_ij-H_ij), SYS_EPSREAL)
    cs_ij = sqrt(c2)

    ! Diagonal matrix of eigenvalues
    l1 = abs(vel-cs_ij)
    l2 = abs(vel)
    l3 = abs(vel+cs_ij)
    l4 = abs(vel)
    
    ! Matrix of right eigenvectors
    R_ij(1,1) =  l1
    R_ij(2,1) =  l1*(u_ij-cs_ij*cx)
    R_ij(3,1) =  l1*(v_ij-cs_ij*cy)
    R_ij(4,1) =  l1*(H_ij-cs_ij*vel)

    R_ij(1,2) =  l2
    R_ij(2,2) =  l2*u_ij
    R_ij(3,2) =  l2*v_ij
    R_ij(4,2) =  l2*q_ij

    R_ij(1,3) =  l3
    R_ij(2,3) =  l3*(u_ij+cs_ij*cx)
    R_ij(3,3) =  l3*(v_ij+cs_ij*cy)
    R_ij(4,3) =  l3*(H_ij+cs_ij*vel)

    R_ij(1,4) =  0.0_DP
    R_ij(2,4) =  l4*cy
    R_ij(3,4) = -l4*cx
    R_ij(4,4) =  l4*(u_ij*cy-v_ij*cx)

    if (abs(cx) > abs(cy)) then

      ! Matrix of left eigenvectors if C(x) is largest coefficient
      L_ij(1,1) = 0.5_DP*(G1*q_ij+cs_ij*vel)/c2
      L_ij(2,1) = (c2-G1*q_ij)/c2
      L_ij(3,1) = 0.5_DP*(G1*q_ij-cs_ij*vel)/c2
      L_ij(4,1) = (v_ij-vel*cy)/cx

      L_ij(1,2) = 0.5_DP*(-G1*u_ij-cs_ij*cx)/c2
      L_ij(2,2) = G1*u_ij/c2
      L_ij(3,2) = 0.5_DP*(-G1*u_ij+cs_ij*cx)/c2
      L_ij(4,2) = cy

      L_ij(1,3) = 0.5_DP*(-G1*v_ij-cs_ij*cy)/c2
      L_ij(2,3) = G1*v_ij/c2
      L_ij(3,3) = 0.5_DP*(-G1*v_ij+cs_ij*cy)/c2
      L_ij(4,3) = (cy*cy-1)/cx

      L_ij(1,4) =  G2/c2
      L_ij(2,4) = -G1/c2
      L_ij(3,4) =  G2/c2
      L_ij(4,4) =  0.0_DP

    else

      ! Matrix of left eigenvectors if C(Y) is largest coefficient
      L_ij(1,1) = 0.5_DP*(G1*q_ij+cs_ij*vel)/c2
      L_ij(2,1) = (c2-G1*q_ij)/c2
      L_ij(3,1) = 0.5_DP*(G1*q_ij-cs_ij*vel)/c2
      L_ij(4,1) = (vel*cx-u_ij)/cy

      L_ij(1,2) = 0.5_DP*(-G1*u_ij-cs_ij*cx)/c2
      L_ij(2,2) = G1*u_ij/c2
      L_ij(3,2) = 0.5_DP*(-G1*u_ij+cs_ij*cx)/c2
      L_ij(4,2) = (1-cx*cx)/cy
      
      L_ij(1,3) = 0.5_DP*(-G1*v_ij-cs_ij*cy)/c2
      L_ij(2,3) = G1*v_ij/c2
      L_ij(3,3) = 0.5_DP*(-G1*v_ij+cs_ij*cy)/c2
      L_ij(4,3) = -cx
      
      L_ij(1,4) =  G2/c2
      L_ij(2,4) = -G1/c2
      L_ij(3,4) =  G2/c2
      L_ij(4,4) =  0.0_DP

    end if

    ! Compute D_ij = diag(R_ij*|Lbd_ij|*L_ij)*I
    S_ij(1) = -cnrm*( R_ij(1,1)*L_ij(1,1)+&
                      R_ij(1,2)*L_ij(2,1)+&
                      R_ij(1,3)*L_ij(3,1)+&
                      R_ij(1,4)*L_ij(4,1)  )
    S_ij(2) = -cnrm*( R_ij(2,1)*L_ij(1,2)+&
                      R_ij(2,2)*L_ij(2,2)+&
                      R_ij(2,3)*L_ij(3,2)+&
                      R_ij(2,4)*L_ij(4,2)  )
    S_ij(3) = -cnrm*( R_ij(3,1)*L_ij(1,3)+&
                      R_ij(3,2)*L_ij(2,3)+&
                      R_ij(3,3)*L_ij(3,3)+&
                      R_ij(3,4)*L_ij(4,3)  )
    S_ij(4) = -cnrm*( R_ij(4,1)*L_ij(1,4)+&
                      R_ij(4,2)*L_ij(2,4)+&
                      R_ij(4,3)*L_ij(3,4)+&
                      R_ij(4,4)*L_ij(4,4)  )

  end subroutine euler_calcMatrixTensorDissDiag2d

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcMatrixTensorDiss2d(U_i, U_j, C_ij, C_ji, A_ij, S_ij)

!<description>
    ! This subroutine computes the diagonal of the local Roe matrices 
    ! and applies tensorial artificial viscosities in 2D
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(IN) :: U_i,U_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!</input>

!<output>
    ! local Roe matrices
    real(DP), dimension(:), intent(OUT) :: A_ij,S_ij
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(NVAR2D,NVAR2D) :: R_ij,L_ij
    real(DP), dimension(NDIM2D) :: a
    real(DP) :: aux,aux1,u2,v2,uv,hi,hj,H_ij,q_ij,u_ij,v_ij
    real(DP) :: l1,l2,l3,l4,cnrm,cx,cy,c2,cs_ij,vel

    ! Compute Roe mean values
    aux  =   sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
    u_ij = ( aux*U_i(2)/U_i(1)+U_j(2)/U_j(1) )/(aux+1.0_DP)
    v_ij = ( aux*U_i(3)/U_i(1)+U_j(3)/U_j(1) )/(aux+1.0_DP)
    hi   =   GAMMA*U_i(4)/U_i(1)-&
             G2*( U_i(2)*U_i(2)+U_i(3)*U_i(3) )/(U_i(1)*U_i(1))
    hj   =   GAMMA*U_j(4)/U_j(1)-&
             G2*( U_j(2)*U_j(2)+U_j(3)*U_j(3) )/(U_j(1)*U_j(1))
    H_ij = ( aux*hi+hj )/(aux+1.0_DP)

    ! Compute coefficients
    a = 0.5_DP*(C_ji-C_ij)

    ! Compute auxiliary variables
    aux1  = u_ij*a(1) + v_ij*a(2)
    u2    = u_ij*u_ij
    v2    = v_ij*v_ij
    uv    = u_ij*v_ij
    q_ij  = 0.5_DP*(u2+v2)
    
    ! Compute Roe matrix for skew-symmetric part
    A_ij( 1) =   0.0_DP
    A_ij( 2) =   (G1*q_ij-u2)*a(1)   - uv*a(2)
    A_ij( 3) = - uv*a(1)             + (G1*q_ij-v2)*a(2)
    A_ij( 4) =   (G1*q_ij-H_ij)*aux1
    
    A_ij( 5) =   a(1)
    A_ij( 6) =   aux1-G6*u_ij*a(1)
    A_ij( 7) =   v_ij*a(1)           - G1*u_ij*a(2)
    A_ij( 8) =   (H_ij-G1*u2)*a(1)   - G1*uv*a(2)
    
    A_ij( 9) =                         a(2)
    A_ij(10) = - G1*v_ij*a(1)        + u_ij*a(2)
    A_ij(11) =                         aux1-G6*v_ij*a(2)
    A_ij(12) = - G1*uv*a(1)          + (H_ij-G1*v2)*a(2)

    A_ij(13) =   0.0_DP
    A_ij(14) =   G1*a(1)
    A_ij(15) =                         G1*a(2)
    A_ij(16) =   GAMMA*aux1

    ! Characteristic velocity
    cnrm = sqrt(a(1)*a(1)+a(2)*a(2))
    if (cnrm .eq. 0.0_DP) return
    cx    = a(1)/cnrm
    cy    = a(2)/cnrm
    vel   = cx*u_ij+cy*v_ij

    ! Compute speed of sound
    c2    = max(-G1*(q_ij-H_ij), SYS_EPSREAL)
    cs_ij = sqrt(c2)

    ! Diagonal matrix of eigenvalues
    l1 = abs(vel-cs_ij)
    l2 = abs(vel)
    l3 = abs(vel+cs_ij)
    l4 = abs(vel)
    
    ! Matrix of right eigenvectors
    R_ij(1,1) =  l1
    R_ij(2,1) =  l1*(u_ij-cs_ij*cx)
    R_ij(3,1) =  l1*(v_ij-cs_ij*cy)
    R_ij(4,1) =  l1*(H_ij-cs_ij*vel)

    R_ij(1,2) =  l2
    R_ij(2,2) =  l2*u_ij
    R_ij(3,2) =  l2*v_ij
    R_ij(4,2) =  l2*q_ij

    R_ij(1,3) =  l3
    R_ij(2,3) =  l3*(u_ij+cs_ij*cx)
    R_ij(3,3) =  l3*(v_ij+cs_ij*cy)
    R_ij(4,3) =  l3*(H_ij+cs_ij*vel)

    R_ij(1,4) =  0.0_DP
    R_ij(2,4) =  l4*cy
    R_ij(3,4) = -l4*cx
    R_ij(4,4) =  l4*(u_ij*cy-v_ij*cx)

    if (abs(cx) > abs(cy)) then
      ! Matrix of left eigenvectors if C(x) is largest coefficient
      L_ij(1,1) = 0.5_DP*(G1*q_ij+cs_ij*vel)/c2
      L_ij(2,1) = (c2-G1*q_ij)/c2
      L_ij(3,1) = 0.5_DP*(G1*q_ij-cs_ij*vel)/c2
      L_ij(4,1) = (v_ij-vel*cy)/cx

      L_ij(1,2) = 0.5_DP*(-G1*u_ij-cs_ij*cx)/c2
      L_ij(2,2) = G1*u_ij/c2
      L_ij(3,2) = 0.5_DP*(-G1*u_ij+cs_ij*cx)/c2
      L_ij(4,2) = cy

      L_ij(1,3) = 0.5_DP*(-G1*v_ij-cs_ij*cy)/c2
      L_ij(2,3) = G1*v_ij/c2
      L_ij(3,3) = 0.5_DP*(-G1*v_ij+cs_ij*cy)/c2
      L_ij(4,3) = (cy*cy-1)/cx

      L_ij(1,4) =  G2/c2
      L_ij(2,4) = -G1/c2
      L_ij(3,4) =  G2/c2
      L_ij(4,4) =  0.0_DP
    else
      ! Matrix of left eigenvectors if C(Y) is largest coefficient
      L_ij(1,1) = 0.5_DP*(G1*q_ij+cs_ij*vel)/c2
      L_ij(2,1) = (c2-G1*q_ij)/c2
      L_ij(3,1) = 0.5_DP*(G1*q_ij-cs_ij*vel)/c2
      L_ij(4,1) = (vel*cx-u_ij)/cy

      L_ij(1,2) = 0.5_DP*(-G1*u_ij-cs_ij*cx)/c2
      L_ij(2,2) = G1*u_ij/c2
      L_ij(3,2) = 0.5_DP*(-G1*u_ij+cs_ij*cx)/c2
      L_ij(4,2) = (1-cx*cx)/cy
      
      L_ij(1,3) = 0.5_DP*(-G1*v_ij-cs_ij*cy)/c2
      L_ij(2,3) = G1*v_ij/c2
      L_ij(3,3) = 0.5_DP*(-G1*v_ij+cs_ij*cy)/c2
      L_ij(4,3) = -cx
      
      L_ij(1,4) =  G2/c2
      L_ij(2,4) = -G1/c2
      L_ij(3,4) =  G2/c2
      L_ij(4,4) =  0.0_DP
    end if

    ! Compute D_ij = R_ij*|Lbd_ij|*L_ij
    call DGEMM('n', 'n', NVAR2D, NVAR2D, NVAR2D, -cnrm,&
               R_ij, NVAR2D, L_ij, NVAR2D, 0.0_DP, S_ij, NVAR2D)

  end subroutine euler_calcMatrixTensorDiss2d

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
    real(DP) :: u_ij,v_ij,H_ij,q_ij,c_ij,aux,hi,hj,c2,u2,v2,cx,cy,cnrm,vel
    real(DP), dimension(NVAR2D) :: diff

    ! Roe mean values
    aux  = sqrt(max(U_i(1)/U_j(1), SYS_EPSREAL))
    u_ij = (aux*U_i(2)/U_i(1)+U_j(2)/U_j(1))/(aux+1.0_DP)
    v_ij = (aux*U_i(3)/U_i(1)+U_j(3)/U_j(1))/(aux+1.0_DP)

    hi   = GAMMA*U_i(4)/U_i(1)-G2*(U_i(2)*U_i(2)+U_i(3)*U_i(3))/(U_i(1)*U_i(1))
    hj   = GAMMA*U_j(4)/U_j(1)-G2*(U_j(2)*U_j(2)+U_j(3)*U_j(3))/(U_j(1)*U_j(1))
    H_ij = (aux*hi+hj)/(aux+1.0_DP)

    ! auxiliary variables
    u2   = u_ij*u_ij
    v2   = v_ij*v_ij
    q_ij = 0.5_DP*(u2+v2)

    ! speed of sound
    c2   = max(-G1*(q_ij-H_ij), SYS_EPSREAL)
    c_ij = sqrt(c2)

    ! characteristic velocity
    cnrm = sqrt(Dweight(1)*Dweight(1)+Dweight(2)*Dweight(2))
    if (cnrm .eq. 0.0_DP) then
      if (present(Lbd_ij)) Lbd_ij = 0.0_DP
      if (present(R_ij))   R_ij   = 0.0_DP
      if (present(L_ij))   L_ij   = 0.0_DP
      W_ij = 0.0_DP
      return
    end if
    cx  = Dweight(1)/cnrm
    cy  = Dweight(2)/cnrm
    vel = cx*u_ij+cy*v_ij

    ! Compute solution difference
    diff = U_j-U_i

    ! Compute diagonal matrix of eigenvalues
    if (present(Lbd_ij)) then
      Lbd_ij(1) = vel-c_ij
      Lbd_ij(2) = vel
      Lbd_ij(3) = vel+c_ij
      Lbd_ij(4) = vel
    end if

    ! Compute matrix of right eigenvectors
    if (present(R_ij)) then
      ! Matrix of right eigenvectors
      R_ij( 1) =  1.0_DP
      R_ij( 2) =  u_ij-c_ij*cx
      R_ij( 3) =  v_ij-c_ij*cy
      R_ij( 4) =  H_ij-c_ij*vel
      R_ij( 5) =  1.0_DP
      R_ij( 6) =  u_ij
      R_ij( 7) =  v_ij
      R_ij( 8) =  q_ij
      R_ij( 9) =  1.0_DP
      R_ij(10) =  u_ij+c_ij*cx
      R_ij(11) =  v_ij+c_ij*cy
      R_ij(12) =  H_ij+c_ij*vel
      R_ij(13) =  0.0_DP
      R_ij(14) =  cy
      R_ij(15) = -cx
      R_ij(16) =  u_ij*cy-v_ij*cx
    end if

    ! Compute matrix of left eigenvectors
    if (present(L_ij)) then
      if (abs(cx) > abs(cy)) then
        ! Matrix of left eigenvectors if CX is largest coefficient
        L_ij( 1) =  0.5_DP*(G1*q_ij+c_ij*vel)/c2
        L_ij( 2) = (c2-G1*q_ij)/c2
        L_ij( 3) =  0.5_DP*(G1*q_ij-c_ij*vel)/c2
        L_ij( 4) = (v_ij-vel*cy)/cx
        L_ij( 5) =  0.5_DP*(-G1*u_ij-c_ij*cx)/c2
        L_ij( 6) =  G1*u_ij/c2
        L_ij( 7) =  0.5_DP*(-G1*u_ij+c_ij*cx)/c2
        L_ij( 8) =  cy
        L_ij( 9) =  0.5_DP*(-G1*v_ij-c_ij*cy)/c2
        L_ij(10) =  G1*v_ij/c2
        L_ij(11) =  0.5_DP*(-G1*v_ij+c_ij*cy)/c2
        L_ij(12) = (cy*cy-1)/cx
        L_ij(13) =  G2/c2
        L_ij(14) = -G1/c2
        L_ij(15) =  G2/c2
        L_ij(16) =  0.0_DP
      else
        ! Matrix of left eigenvectors if CY is largest coefficient
        L_ij( 1) =  0.5_DP*(G1*q_ij+c_ij*vel)/c2
        L_ij( 2) = (c2-G1*q_ij)/c2
        L_ij( 3) =  0.5_DP*(G1*q_ij-c_ij*vel)/c2
        L_ij( 4) = (vel*cx-u_ij)/cy
        L_ij( 5) =  0.5_DP*(-G1*u_ij-c_ij*cx)/c2
        L_ij( 6) =  G1*u_ij/c2
        L_ij( 7) =  0.5_DP*(-G1*u_ij+c_ij*cx)/c2
        L_ij( 8) = (1-cx*cx)/cy
        L_ij( 9) =  0.5_DP*(-G1*v_ij-c_ij*cy)/c2
        L_ij(10) =  G1*v_ij/c2
        L_ij(11) =  0.5_DP*(-G1*v_ij+c_ij*cy)/c2
        L_ij(12) = -cx
        L_ij(13) =  G2/c2
        L_ij(14) = -G1/c2
        L_ij(15) =  G2/c2
        L_ij(16) =  0.0_DP
      end if

      ! Compute W_ij = L_ij*(U_j-U_i)
      call DGEMV('n', NVAR2D, NVAR2D, 1.0_DP, L_ij,&
                      NVAR2D, Diff, 1, 0.0_DP, W_ij, 1)

    else

      ! Compute W_ij = L_ij*(U_j-U_i) directly
      if (abs(cx) > abs(cy)) then
        ! CX is largest coefficient
        W_ij(1) = 0.5_DP*(G1*q_ij+c_ij*vel)/c2*Diff(1) +&
                  0.5_DP*(-G1*u_ij-c_ij*cx)/c2*Diff(2)  +&
                  0.5_DP*(-G1*v_ij-c_ij*cy)/c2*Diff(3)  +&
                  G2/c2*Diff(4)
        W_ij(2) = (c2-G1*q_ij)/c2*Diff(1) +&
                  G1*u_ij/c2*Diff(2)      +&
                  G1*v_ij/c2*Diff(3)      -&
                  G1/c2*Diff(4)
        W_ij(3) = 0.5_DP*(G1*q_ij-c_ij*vel)/c2*Diff(1) +&
                  0.5_DP*(-G1*u_ij+c_ij*cx)/c2*Diff(2)  +&
                  0.5_DP*(-G1*v_ij+c_ij*cy)/c2*Diff(3)  +&
                  G2/c2*Diff(4)
        W_ij(4) = (v_ij-vel*cy)/cx*Diff(1) +&
                  cy*Diff(2)               +&
                  (cy*cy-1)/cx*Diff(3)
      else
        ! CY is largest coefficient
        W_ij(1) = 0.5_DP*(G1*q_ij+c_ij*vel)/c2*Diff(1)+&
                  0.5_DP*(-G1*u_ij-c_ij*cx)/c2*Diff(2)+&
                  0.5_DP*(-G1*v_ij-c_ij*cy)/c2*Diff(3)+&
                  G2/c2*Diff(4)
        W_ij(2) = (c2-G1*q_ij)/c2*Diff(1)+&
                  G1*u_ij/c2*Diff(2)+&
                  G1*v_ij/c2*Diff(3)-&
                  G1/c2*Diff(4)
        W_ij(3) = 0.5_DP*(G1*q_ij-c_ij*vel)/c2*Diff(1)+&
                  0.5_DP*(-G1*u_ij+c_ij*cx)/c2*Diff(2)+&
                  0.5_DP*(-G1*v_ij+c_ij*cy)/c2*Diff(3)+&
                  G2/c2*Diff(4)
        W_ij(4) = (vel*cx-u_ij)/cy*Diff(1)+&
                  (1-cx*cx)/cy*Diff(2)-&
                  cx*Diff(3)
      end if
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
