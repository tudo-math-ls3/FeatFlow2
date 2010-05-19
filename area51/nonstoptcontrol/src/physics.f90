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
    ! =0: optimal distributed control of the heat equation, constant viscosity
    integer :: cequation = 0
    
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
