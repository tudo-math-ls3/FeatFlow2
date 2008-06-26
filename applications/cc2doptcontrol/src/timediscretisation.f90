!##############################################################################
!# ****************************************************************************
!# <name> timediscretisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises a time discretisation structure, i.e. a discretisation
!# structure that realises different discretisations in time.
!# It currently supports implicit Euler as well as the dG(0) scheme.
!#
!# The following routines can be found here:
!#
!# 1.) tdiscr_initTheta
!#     -> Initialise a time discretisation structure according to a one-step
!#        theta scheme.
!#
!# 2.) tdiscr_initdG0
!#     -> Initialise a time discretisation structure according to the dG(0)
!#        scheme.
!#
!# 3.) tdiscr_igetNDofGlob
!#     -> Calculate the number of DOF's in time.
!# </purpose>
!##############################################################################

MODULE timediscretisation

  USE fsystem
  USE genoutput

  IMPLICIT NONE

!<constants>

!<constantblock description="Time discretisation type identifiers">

  ! Theta-scheme (implicit/explicit Euler, Crank Nicolson)
  INTEGER, PARAMETER :: TDISCR_THETA           = 0
  
  ! The dG(0) scheme.
  INTEGER, PARAMETER :: TDISCR_DG0             = 1

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! This block realises a time discretisation structure which contains all
  ! parameters that are necessary for the time discretisation of an ODE.
  TYPE t_timeDiscretisation
  
    ! A TDISCR_xxxx-Flag that identifies the time discretisation scheme.
    INTEGER :: ctype = TDISCR_THETA
    
    ! Number of time intervals
    INTEGER :: nintervals          = 0

    ! Absolute start time of the simulation
    REAL(DP) :: dtimeInit          = 0.0_DP     
    
    ! Maximum time of the simulation
    REAL(DP) :: dtimeMax           = 0.0_DP
    
    ! Time step length of the time discretisation
    REAL(DP) :: dtstep             = 0.0_DP

    ! If ctype=TDISCR_THETA: theta parameter that identifies the Theta scheme.
    !  =0.0: Forward Euler,
    !  =1.0: Backward Euler (Standard),
    !  =0.5: Crank Nicolson.
    REAL(DP) :: dtheta = 1.0_DP
  
  END TYPE

!</typeblock>

!</types>

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE tdiscr_initTheta (dstart, dend, nintervals, dtheta, rtimediscr)
  
!<description>
  ! Initialises a time discretisation structure for the discretisation with
  ! a one-step Theta scheme.
!</description>

!<input>
  ! Start point of the time interval.
  REAL(DP), INTENT(IN) :: dstart
  
  ! End point of the time interval
  REAL(DP), INTENT(IN) :: dend
  
  ! Number of subintervals
  INTEGER, INTENT(IN) :: nintervals
  
  ! Theta scheme parameter
  !  =0.0: Forward Euler,
  !  =1.0: Backward Euler (Standard),
  !  =0.5: Crank Nicolson.
  REAL(DP), INTENT(IN) :: dtheta
!</input>

!<output>
  ! The time discrtisation structure to be initialised.
  TYPE(t_timeDiscretisation), INTENT(OUT) :: rtimediscr
!</output>

    !IF (dtheta .NE. 1.0_DP) THEN
    !  CALL output_line ('Currently only implicit Euler is supported.', &
    !                    OU_CLASS_ERROR,OU_MODE_STD,'tdiscr_initTheta')
    !  CALL sys_halt()
    !END IF

    ! Initialise the parameters of the structure.
    rtimediscr%ctype = TDISCR_THETA
    rtimediscr%dtimeInit = dstart
    rtimediscr%dtimeMax = dend
    rtimediscr%nintervals = nintervals
    rtimediscr%dtstep = (dend-dstart)/REAL(nintervals,DP)
    rtimediscr%dtheta = dtheta

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE tdiscr_initdG0(dstart, dend, nintervals, rtimediscr)
  
!<description>
  ! Initialises a time discretisation structure for the discretisation with
  ! dG(0) scheme.
!</description>

!<input>
  ! Start point of the time interval.
  REAL(DP), INTENT(IN) :: dstart
  
  ! End point of the time interval
  REAL(DP), INTENT(IN) :: dend
  
  ! Number of subintervals
  INTEGER, INTENT(IN) :: nintervals
!</input>

!<output>
  ! The time discrtisation structure to be initialised.
  TYPE(t_timeDiscretisation), INTENT(OUT) :: rtimediscr
!</output>

    ! Initialise the parameters of the structure.
    rtimediscr%ctype = TDISCR_DG0
    rtimediscr%dtimeInit = dstart
    rtimediscr%dtimeMax = dend
    rtimediscr%nintervals = nintervals
    rtimediscr%dtstep = (dend-dstart)/REAL(nintervals,DP)

  END SUBROUTINE
  
  ! ***************************************************************************

!<function>  

  INTEGER FUNCTION tdiscr_igetNDofGlob(rtimediscr)

!<description>
  ! This function returns for a given time discretisation the number of 
  ! degrees of freedom in time.
!</description>

!<input>    
  ! The tiome discretisation structure that specifies the time discretisation.
  TYPE(t_timeDiscretisation), INTENT(IN) :: rtimediscr
!</input>

!<result>
  ! Global number of equations on current time grid
!</result>

!</function>

    tdiscr_igetNDofGlob = 0
    
    ! Cancel if the structure is not initialised.
    IF (rtimediscr%nintervals .EQ. 0) RETURN

    SELECT CASE(rtimediscr%ctype)
    CASE (TDISCR_THETA)
      ! One step scheme. We have as many DOF's as intervals + 1.
      tdiscr_igetNDofGlob = rtimediscr%nintervals + 1
      RETURN

    CASE (TDISCR_DG0)
      ! dG(0). The DOF's are in the interval midpoints.
      tdiscr_igetNDofGlob = rtimediscr%nintervals
      RETURN

    END SELECT
    
    CALL output_line ('Unsupported time discretisation.', &
                      OU_CLASS_ERROR,OU_MODE_STD,'tdiscr_igetNDofGlob')
    CALL sys_halt()
    
  END FUNCTION

END MODULE
