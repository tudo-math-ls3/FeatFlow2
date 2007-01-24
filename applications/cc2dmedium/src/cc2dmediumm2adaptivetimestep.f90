!##############################################################################
!# ****************************************************************************
!# <name> cc2dmmediumm2adaptivetimestep </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains time stepping routines to control the adaptive
!# time stepping in a nonstationary simulation.
!# </purpose>
!##############################################################################

MODULE cc2dmmediumm2adaptivetimestep

  USE fsystem
    
  IMPLICIT NONE
  
!<constants>

!<constantblock description="Identifiers for the type of the adaptive time stepping.">

  ! Fixed time step size, no adaptive time stepping.
  INTEGER, PARAMETER :: TADTS_FIXED = 0
  
  ! Adaptive time stepping with prediction,
  ! no repetition except for when the solver breaks down.
  INTEGER, PARAMETER :: TADTS_PREDICTION = 1
  
  ! Adaptive time stepping with prediction,
  ! repetition of the time step if nonlinear stopping
  ! criterion is too large or the solver breaks down.
  INTEGER, PARAMETER :: TADTS_PREDICTREPEAT = 2
  
  ! Adaptive time stepping with prediction,
  ! repetition of the time step if nonlinear stopping
  ! criterion is too large or the time error is too large
  ! or the solver breaks down.
  INTEGER, PARAMETER :: TADTS_PREDREPTIMECONTROL = 3

!</constantblock>

!<constantblock description="Identifiers for the type of the adaptive time stepping during the start phase.">

  ! Use depsadi for controlling the error in the start phase (Standard).
  INTEGER, PARAMETER :: TADTS_START_STANDARD    = 0

  ! Use linear combination of depsadi and depsadl for
  ! controlling the error in the start phase.
  INTEGER, PARAMETER :: TADTS_START_LINEAR      = 1
  
  ! Use logarithmic combination of depsadi and depsadl for
  ! controlling the error in the start phase.
  INTEGER, PARAMETER :: TADTS_START_LOGARITHMIC = 2

!</constantblock>

!</constants>  

  
!<types>

!<typeblock>

  ! Basic structure guiding the adaptive time stepping, in the start phase
  ! and during the simulation, as well as configuring how to estimate the time
  ! error.
  TYPE t_timeErrorControl
  
    ! Type of adaptive time stepping. One of the TADTS_xxxx constants.
    ! Standard is TADTS_FIXED = No adaptive time stepping.
    INTEGER                        :: ctype = TADTS_FIXED
    
    ! Maximum number of repetitions if ctype=TADTS_PREDICTREPEAT or
    ! ctype=TADTS_PREDREPTIMECONTROL.
    INTEGER                        :: nrepetitions = 3
    
    ! Initial time step.
    REAL(DP)                       :: dtinit = 0.1_DP
    
    ! Minimum time step.
    REAL(DP)                       :: dtmin = 0.001_DP
    
    ! Maximum time step.
    REAL(DP)                       :: dtmax = 1.0_DP
    
    ! Factor for modifying time step size. 
    ! Only if ctype != TADTS_FIXED.
    ! Time step size is reduced by sqrt(DTFACT) upon broke down
    ! of nonlinear solver, and by DTFACT upon broke down of linear solver.
    REAL(DP)                       :: dtfact = 9.0_DP
  
    ! Type of error control in space/time; former IEPSAD.
    ! =1: Calculate RELT=RELU2
    ! =2: Calculate RELT=RELUM
    ! =3: Calculate RELT=RELP2
    ! =4: Calculate RELT=RELPM
    ! =5: Calculate RELT=max(RELU2,RELP2)
    ! =6: Calculate RELT=max(RELUM,RELPM)
    ! =7: Calculate RELT=max(RELU2,RELP2,RELUM,RELPM)
    ! =8: Calculate RELT=min(RELU2,RELP2,RELUM,RELPM)  
    INTEGER                        :: ctimeErrorControl = 1
    
    ! Extrapolation in time. If .TRUE., the solution in the new time step
    ! is calculated by a weighted linear combination of the predictor and corrector
    ! step. Can only be used if ctype != TADTS_FIXED
    LOGICAL                        :: bextrapolationInTime = .FALSE.
    
    ! Upper limit for relation of calculated step sizes to accept
    ! a time step; if the relation 
    !                (new step size)/(old step size)
    ! is too large, repeat the step. Applies only to the
    ! simulation if ctype=TADTS_PREDREPTIMECONTROL!
    REAL(DP)                       :: depsadu = 0.5_DP
  
    ! Parameter for error control in the start phase. One of the
    ! TADTS_START_xxxx constants. Standard value is TADTS_START_STANDARD,
    ! which controls the time stepping in the start phase by depsadi only.
    ! This parameter affects the simulation only if 
    ! ctype=TADTS_PREDREPTIMECONTROL!
    INTEGER                        :: cadaptiveStartPhase = TADTS_START_STANDARD
    
    ! Parameter for time error limit in the start phase. 
    ! Applies only to the simulation if ctype=TADTS_PREDREPTIMECONTROL!
    REAL(DP)                       :: depsadi = 0.125_DP

    ! Parameter for time error limit after the start phase. 
    ! Applies only to the simulation if ctype=TADTS_PREDREPTIMECONTROL!
    REAL(DP)                       :: depsadl = 0.00125_DP
  
  END TYPE

!</typeblock>

!</types>
  
!CONTAINS
  
!******************************************************************************

END MODULE
