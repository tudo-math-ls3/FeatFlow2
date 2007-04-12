!##############################################################################
!# ****************************************************************************
!# <name> cc2dmmediumm2optcanalysis </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to perform error analysis for the solution
!# of the optimal control problem. Here, routines can be found that calculate
!# the functional, which is to be minimised by the theory:
!#
!# $$ min J(y,u) = 1/2||y-z||_{L^2} + \gamma/2||y(T)-z(T)||_{L^2} + \alpha/2||u||^2 $$
!#
!# The following routines can be found here:
!#
!# 1.) c2d2_optc_stationaryFunctional
!#     -> Calculates the value of the stationary functional which is to be
!#        minimised:
!#     $$ J(y,u) = 1/2||y-z||_{L^2}  + \alpha/2||u||^2 $$
!#
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2optcanalysis

  USE fsystem
  USE linearsystemblock
  USE pprocerror
  USE cc2dmedium_callback
    
  IMPLICIT NONE
  
CONTAINS

!******************************************************************************

!<subroutine>

  SUBROUTINE c2d2_optc_stationaryFunctional (rsolution,dalpha,derror)

!<description>
  ! This function calculates the value of the functional which is to be
  ! minimised in the stationary optimal control problem.
  ! The functional is defined as
  !   $$ J(y,u) = 1/2||y-z||_{L^2}  + \alpha/2||u||^2 $$
  ! over a spatial domain $\Omega$.
  !
  ! For this purpose, the routine evaluates the user-defined callback functions 
  ! coeff_TARGET_x and coeff_TARGET_y.
!</description>
  
!<input>
  ! Solution vector to compute the norm/error from.
  TYPE(t_vectorBlock), INTENT(IN) :: rsolution
  
  ! Regularisation parameter $\alpha$.
  REAL(DP), INTENT(IN) :: dalpha
!</input>

!<output>
  ! Norm of the error functional.
  REAL(DP), INTENT(OUT) :: derror
!</output>
  
!</subroutine>
    
    ! local variables
    REAL(DP),DIMENSION(2) :: Derr
    
    ! Perform error analysis to calculate and add 1/2||y-z||_{L^2}.
    CALL pperr_scalar (rsolution%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
                       ffunction_TargetX)

    CALL pperr_scalar (rsolution%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
                       ffunction_TargetY)
                       
    derror = 0.5_DP * (0.5_DP*(Derr(1)**2+Derr(2)**2))
    
    ! Calculate \alpha/2||u||^2.
    IF (dalpha .NE. 0.0_DP) THEN
      CALL pperr_scalar (rsolution%RvectorBlock(4),PPERR_L2ERROR,Derr(1))

      CALL pperr_scalar (rsolution%RvectorBlock(5),PPERR_L2ERROR,Derr(2))
                         
      derror = derror + (0.5_DP/dalpha) * (0.5_DP*(Derr(1)**2+Derr(2)**2))
      
      ! Because of u=-lambda/alpha we have:
      !    alpha/2 ||u||^2 = alpha/2 ||lambda/alpha||^2 = 1/(2*alpha) ||lambda||^2
    END IF
    
  END SUBROUTINE

END MODULE
