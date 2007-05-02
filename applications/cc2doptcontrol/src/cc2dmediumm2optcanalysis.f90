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

  SUBROUTINE c2d2_optc_stationaryFunctional (rsolution,dalpha,Derror,rcollection)

!<description>
  ! This function calculates the value of the functional which is to be
  ! minimised in the stationary optimal control problem.
  ! The functional is defined as
  !   $$ J(y,u) = 1/2||y-z||_{L^2}^2  + \alpha/2||u||^2 $$
  ! over a spatial domain $\Omega$.
  !
  ! For this purpose, the routine evaluates the user-defined callback functions 
  ! ffunction_TargetX and ffunction_TargetY. The collection must be initialised
  ! for postprocessing before calling his routine!
!</description>
  
!<input>
  ! Solution vector to compute the norm/error from.
  TYPE(t_vectorBlock), INTENT(IN) :: rsolution
  
  ! Regularisation parameter $\alpha$.
  REAL(DP), INTENT(IN) :: dalpha
  
  ! Collection structure of the main application. Is passed to the callback routines.
  TYPE(t_collection), INTENT(IN) :: rcollection
!</input>

!<output>
  ! Returns information about the error.
  ! Derror(1) = ||y-z||_{L^2}.
  ! Derror(2) = ||u||.
  ! Derror(3) = J(y,u) = 1/2||y-z||_{L^2}^2  + \alpha/2||u||^2.
  ! Norm of the error functional.
  REAL(DP), DIMENSION(3), INTENT(OUT) :: derror
!</output>
  
!</subroutine>
    
    ! local variables
    REAL(DP),DIMENSION(2) :: Derr
    
    ! Perform error analysis to calculate and add 1/2||y-z||_{L^2}.
    CALL pperr_scalar (rsolution%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
                       ffunction_TargetX,rcollection)

    CALL pperr_scalar (rsolution%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
                       ffunction_TargetY,rcollection)
                       
    Derror(1) = (0.5_DP*(Derr(1)**2+Derr(2)**2))
    Derror(2) = 0.0_DP
    
    ! Calculate \alpha/2||u||^2.
    IF (dalpha .NE. 0.0_DP) THEN
      CALL pperr_scalar (rsolution%RvectorBlock(4),PPERR_L2ERROR,Derr(1))

      CALL pperr_scalar (rsolution%RvectorBlock(5),PPERR_L2ERROR,Derr(2))
                         
      Derror(2) = (0.5_DP*(Derr(1)**2+Derr(2)**2))
      
      ! Because of u=-lambda/alpha we have:
      !    alpha/2 ||u||^2 = alpha/2 ||lambda/alpha||^2 = 1/(2*alpha) ||lambda||^2
    END IF
    
    Derror(3) = 0.5_DP * Derror(1)+(0.5_DP/dalpha) * Derror(2)
    Derror(1) = SQRT(Derror(1))
    Derror(2) = SQRT(Derror(2))
    
  END SUBROUTINE

!******************************************************************************

!<subroutine>

  SUBROUTINE c2d2_optc_nonstatFunctional (rproblem,rsolution,dalpha,derror)

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
  
!<inputoutput>
  ! Problem structure. Defines the current point in time and provides a collection
  ! to be used for the callback routines
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>
  
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

    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    CALL c2d2_initCollectForAssembly (rproblem,rproblem%rcollection)

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
    
    ! Clean up the collection (as we are done with the assembly, that's it.
    CALL c2d2_doneCollectForAssembly (rproblem,rproblem%rcollection)
    
  END SUBROUTINE

END MODULE
