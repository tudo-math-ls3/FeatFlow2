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
!# 1.) cc_optc_stationaryFunctional
!#     -> Calculates the value of the stationary functional which is to be
!#        minimised:
!#     $$ J(y,u) = 1/2||y-z||_{L^2} + \alpha/2||u||^2 $$
!#
!# 1.) cc_optc_nonstatFunctional
!#     -> Calculates the value of the stationary functional which is to be
!#        minimised:
!#     $$ J(y,u) = 1/2||y-z||^2_{L^2} + \alpha/2||u||^2_{L^2} + gamma/2||y(T)-z(T)||^2_{L^2}$$
!#
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2optcanalysis

  USE fsystem
  USE linearsystemblock
  USE pprocerror
  USE spacetimevectors
  USE cc2dmedium_callback
    
  IMPLICIT NONE
  
CONTAINS

!******************************************************************************

!<subroutine>

  SUBROUTINE cc_optc_stationaryFunctional (rsolution,dalpha,Derror,rcollection)

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
  TYPE(t_collection), INTENT(INOUT) :: rcollection
!</input>

!<output>
  ! Returns information about the error.
  ! Derror(1) = ||y-z||_{L^2}.
  ! Derror(2) = ||u||.
  ! Derror(3) = J(y,u) = 1/2||y-z||_{L^2}^2  + \alpha/2||u||^2.
  ! Norm of the error functional.
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Derror
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

  SUBROUTINE cc_optc_nonstatFunctional (rproblem,rsolution,rtempVector,&
      dalpha,dgamma,Derror)

!<description>
  ! This function calculates the value of the functional which is to be
  ! minimised in the stationary optimal control problem.
  ! The functional is defined as
  !   $$ J(y,u) = 1/2||y-z||^2_{L^2} + \alpha/2||u||^2_{L^2} + gamma/2||y(T)-z(T)||^2_{L^2}$$
  ! over a spatial domain $\Omega$.
  !
  ! For this purpose, the routine evaluates the user-defined callback functions 
  ! ffunction_TargetX and ffunction_TargetY. The collection must be initialised
  ! for postprocessing before calling his routine!
  !
  ! The function z is given implicitely in the problem structure rproblem
  ! and evaluated in ffunction_TargetX and ffunction_TargetY!
!</description>
  
!<input>
  ! Solution vector to compute the norm/error from.
  TYPE(t_spacetimeVector), INTENT(IN) :: rsolution
  
  ! Regularisation parameter $\alpha$.
  REAL(DP), INTENT(IN) :: dalpha

  ! Regularisation parameter $\gamma$.
  REAL(DP), INTENT(IN) :: dgamma
!</input>

!<inputoutput>
  ! Problem structure defining z and the time discretisation.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

  ! A block temp vector with size and structure of the subvectors in rsolution.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector
!</inputoutput>

!<output>
  ! Returns information about the error.
  ! Derror(1) = ||y-z||_{L^2}.
  ! Derror(2) = ||u||_{L^2}.
  ! Derror(3) = ||y(T)-z(T)||_{L^2}.
  ! Derror(4) = J(y,u).
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Derror
!</output>
  
!</subroutine>
    
    ! local variables
    INTEGER :: isubstep
    REAL(DP) :: dtstep
    REAL(DP),DIMENSION(2) :: Derr
    
    Derror(1:4) = 0.0_DP
    dtstep = rsolution%p_rtimeDiscretisation%dtstep

    DO isubstep = 0,rsolution%NEQtime-1
      ! Current point in time
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + &
          isubstep*dtstep
      rproblem%rtimedependence%itimestep = isubstep

      ! Get the solution.
      ! Evaluate the space time function in rvector in the point
      ! in time dtime. Independent of the discretisation in time,
      ! this will give us a vector in space.
      !CALL sptivec_getTimestepData (rsolution, 1+isubstep, rtempVector)
      CALL tmevl_evaluate(rsolution,rproblem%rtimedependence%dtime,rtempVector)

      ! Initialise the collection for the assembly process with callback routines.
      ! This stores the simulation time in the collection and sets the
      ! current subvector z for the callback routines.
      CALL cc_initCollectForAssembly (rproblem,rproblem%rcollection)
      
      ! Compute:
      ! Derror(1) = ||y-z||^2_{L^2}.
      
      ! Perform error analysis to calculate and add 1/2||y-z||_{L^2}.
      CALL pperr_scalar (rtempVector%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
                        ffunction_TargetX,rproblem%rcollection)

      CALL pperr_scalar (rtempVector%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
                        ffunction_TargetY,rproblem%rcollection)

      Derror(1) = Derror(1) + 0.5_DP*(Derr(1)**2 + Derr(2)**2)

      ! Compute:
      ! Derror(3) = ||y(T)-z(T)||^2
      IF (isubstep .EQ. rsolution%NEQtime-1) THEN
        Derror(3) = 0.5_DP*(Derr(1)**2+Derr(2)**2)
      END IF
      
      ! Compute:
      ! Derror(2) = ||lambda||^2_{L^2}.
      ! (intermediate solution; ||u||=-1/alpha ||lambda|| !
      CALL pperr_scalar (rtempVector%RvectorBlock(4),PPERR_L2ERROR,Derr(1))
      CALL pperr_scalar (rtempVector%RvectorBlock(5),PPERR_L2ERROR,Derr(2))
                          
      Derror(2) = Derror(2) + 0.5_DP*(Derr(1)**2+Derr(2)**2)
      
      CALL cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
    END DO
    
    ! Normalise...
    Derror(1) = Derror(1) / REAL(rsolution%NEQtime,DP)
    Derror(2) = Derror(2) / REAL(rsolution%NEQtime,DP)
    
    ! Calculate J(.)
    Derror(4) = 0.5_DP * Derror(1)  +  0.5_DP * dgamma * Derror(3)
    
    IF (dalpha .NE. 0.0_DP) THEN
      ! Because of u=-lambda/alpha we have:
      !    alpha/2 ||u||^2 = alpha/2 ||lambda/alpha||^2 = 1/(2*alpha) ||lambda||^2
      Derror(4) = Derror(4) + 0.5_DP/dalpha * Derror(2)
      
      ! Calculate ||u|| = 1/alpha ||lambda||
      Derror(2) = 1.0_DP/dalpha * SQRT(Derror(2))
    ELSE
      Derror(2) = 0.0_DP
    END IF
    
    ! And the rest
    Derror(3) = SQRT(Derror(3))
    Derror(1) = SQRT(Derror(1))
    
  END SUBROUTINE

END MODULE
