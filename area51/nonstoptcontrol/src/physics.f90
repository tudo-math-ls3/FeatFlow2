!##############################################################################
!# ****************************************************************************
!# <name> physics </name>
!# ****************************************************************************
!#
!# <purpose>
!# Encapsules the physics of the problem.
!# </purpose>
!##############################################################################

module physics

  use fsystem
  use storage
  
  implicit none

  ! Encapsules the physics of the problem.
  type t_physics
  
    ! Type of the underlying operator.
    ! optimal distributed control of the ...
    ! =0: 2D heat equation
    ! =1: 2D Stokes
    ! =2: 1D heat equation
    integer :: cequation = 0
    
    ! Type of reference problem. Depending on cequation.
    ! creferenceproblem = 0: no chosen reference problem.
    !
    ! cequation = 0/2: (1D/2D heat equation)
    !
    !   creferenceproblem = 1:
    !      y      = t^2 x
    !      lambda = -2 t x
    !
    !   creferenceproblem = 2:
    !      y      = t^2 (1-t)^2 x
    !      lambda = t (1-t) x
    !
    !   creferenceproblem = 3:
    !      y      = t^2 (1-t)^2 x
    !      lambda = t (1-t)^2 x
    !
    !   creferenceproblem = 4:
    !      y      = t^2 (1-t)^2 x
    !      lambda = passend, so dass f=0; incl. alpha.
    !
    ! cequation = 1: (2D Stokes equation)
    !
    !   creferenceproblem = 1:
    !      y      = t^2 (x1,-x2)
    !      lambda = -2 t (x1,-x2)
    !
    !   creferenceproblem = 2:
    !      y      = t^2 (1-t)^2 x2(1-x2)
    !      lambda = passend incl. alpha
    !
    !   creferenceproblem = 3:
    !      y      = t^2 (x1^2 x2,-x1 x2^2)
    !      lambda = passend incl. alpha
    !
    integer :: creferenceproblem = 0
    
    ! Viscosity parameter
    real(DP) :: dviscosity = 0.0_DP
    
    ! Relaxation parameter ALPHA of the optimal control
    real(DP) :: doptControlAlpha = 0.0_DP
    
    ! Relaxation parameter GAMMA for the terminal condition
    real(DP) :: doptControlGamma = 0.0_DP
    
    ! Couple the dual equation to the primal
    real(DP) :: dcoupleDualToPrimal = 1.0_DP

    ! Couple the primal equation to the dual
    real(DP) :: dcouplePrimalToDual = 1.0_DP
    
    ! Couple the terminal condition to the dual solution.
    real(DP) :: dcoupleTermCond = 1.0_DP

  end type

end module
