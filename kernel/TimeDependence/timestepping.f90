!##############################################################################
!# ****************************************************************************
!# <name> timestepping </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a realisation of 1D time-stepping schemes used 
!# for the time discretisation of PDE's. Examples for this are Explicit Euler,
!# Crank Nicolson or the Fractional Step Theta scheme.
!#
!# The basic time stepping is governed by the structure t_explicitTimeStepping,
!# which is maintained by the routines in this module. It contains information
!# about the length of the next time step, number of current time step,
!# current simulation time, etc.
!#
!# The following routines can be found here:
!#
!# 1.) timstp_init
!#     -> Initialise time stepping structure for advancing in time
!#
!# 2.) timstp_nextSubstep
!#     -> Go to next time (sub-)step
!#
!# 3.) timstp_getOrder
!#     -> Retrieves the order of the time error that is expected from a
!#        defined time stepping algorithm.
!#
!# How to use this module?
!#
!# Well, that's simple. It's like this:
!#
!#    TYPE(t_explicitTimeStepping) :: tstepping
!#
!#    ! Initialise time stepping  
!#    CALL timstp_init (tstepping,...)
!#
!#    DO istep = 1,#max number of time steps
!# 
!#      ... (do some work at the current time tstepping%dcurrentTime)
!#
!#      ! Proceed to next time step
!#      CALL timstp_nextSubstep (tstepping)
!#
!#    END DO
!#
!# First initialise the time stepping scheme by timstp_init, then call
!# timstp_nextSubstep in a loop to proceed from time dcurrentTime
!# to time dcurrentTime+dtstep. That's all.
!#
!# The time stepping structure contains definitions for various constants that
!# are to be used in front of the different terms in the PDE. See the 
!# documentation of t_explicitTimeStepping which constants are to be used
!# for what.
!# </purpose>
!##############################################################################

MODULE timestepping

  USE fsystem
  USE genoutput

  IMPLICIT NONE
  
!<constants>

!<constantblock description="Identifiers for time step schemes">

  ! Identifier for one step scheme (e.g. Explicit Euler)
  INTEGER, PARAMETER :: TSCHM_ONESTEP        = 0
  
  ! Identifier for Fractional Step scheme
  INTEGER, PARAMETER :: TSCHM_FRACTIONALSTEP = 1

!</constantnblock>

!</constants>
  
  
!<types>

!<typeblock>

  ! Time stepping structure used in explicit time stepping schemes like
  ! explicit Euler, Crank Nicolson etc.
  ! The time stepping scheme is designed to discretise a rather general
  ! equation in time for a solution $u$ and a time variable $t$:
  !
  !               $$ u_t(x,t) + N(u(x,t)) = f(x,t) $$
  !
  ! Given a time step $k$ and a time stepping scheme parameter $\theta$, 
  ! this is discretised in as:
  !
  !  $$ (u_{n+1} - u_n) / k  +  \theta N(u_{n+1})
  !
  !    = -(1-\theta) N(u_n)  +  \theta f_{n+1}  +  (1-\theta) f_n $$
  !
  ! Now, regroup this equation with $u_{n+1}$ on the LHS and $u_n$ on the
  ! RHS and give the coefficients in front of all terms an appropriate
  ! name. Then the time stepping scheme realised by this structure yields:
  !
  ! $$ u_{n+1} + dweightMatrixLHS*N(u_n+1) 
  !
  ! =  u_n + dweightMatrixRHS*N(u_n)  +  dweightNewRHS*f_{n+1}  +  dweightOldRHS*f_n $$
  !
  ! If the simulation has steady (in)homogenuous boundary conditions
  ! and a One-Step Theta-Scheme is used, the RHS build from f_{n+1}
  ! and f_n can be simplified. In this case:
  ! 
  !  $$  dweightNewRHS*f_n+1  +  dweightOldRHS*f_n
  !    = dweightNewRHS*f_n    +  dweightOldRHS*f_n
  !    = dweightStationaryRHS*f_n                    $$

  TYPE t_explicitTimeStepping
  
    ! Time step scheme identifier. One of the TSCHM_xxxx-constants.
    ! Usually TSCHM_ONESTEP for one step schemes
    INTEGER                  :: ctimestepType = -1
    
    ! Current substep in time step scheme. Normally =1.
    ! For time stepping schemes with multiple substeps (Fractional Step),
    ! this counts the current substep 1..nsubsteps.
    INTEGER                  :: isubstep
    
    ! Number of substeps of the time stepping scheme. Normally =1.
    ! For time stepping schemes with multiple substeps (Fractional Step),
    ! this is the number of substeps.
    INTEGER                  :: nsubsteps
    
    ! Current simulation time, this structure represents
    REAL(DP)                 :: dcurrentTime
    
    ! Current simulation time of the macrostep. For standard $\theta$-schemes,
    ! there is dcurrentTime=dtimeMacrostep. For schemes with multiple substeps
    ! like Fractional-Step, dcurrentTime is the current time of the substep, 
    ! while dtimeMacrostep saves the 'last' simulation time stamp with 
    ! "isubstep=1", i.e. the beginning of the 'macrostep'.
    REAL(DP)                 :: dtimeMacrostep
    
    ! Length of the next time step
    REAL(DP)                 :: dtstep
    
    ! Main configuration parameter of the Theta scheme
    REAL(DP)                 :: dtheta
    
    ! Theta-scheme identifier for the substeps of the step.
    ! For 1-step schemes:
    !  =0: Forward Euler,
    !  =1: Backward Euler,
    !  =0.5: Crank Nicolson.
    ! Special values for Fractional Step.
    REAL(DP)                 :: dthStep
    
    ! Weight in front of the matrix on the LHS.
    REAL(DP)                 :: dweightMatrixLHS

    ! Weight in front of the matrix on the RHS.
    REAL(DP)                 :: dweightMatrixRHS
    
    ! Weight to be used in the current time step for the "new" RHS $f_{n+1}$.
    REAL(DP)                 :: dweightNewRHS

    ! Weight to be used in the current time step for the "old" RHS $f_{n}$.
    REAL(DP)                 :: dweightOldRHS
    
    ! Weight to be used in the current time step if a stationary RHS is used
    ! (instead of a combination of dweightNewRHS and dweightOldRHS).
    REAL(DP)                 :: dweightStationaryRHS
    
    ! "Adjungated" parameter Theta' of the Theta scheme.
    ! Only used internally in the Fractional Step scheme.
    ! Standard = 0.0 = no fractional step used.
    REAL(DP)                 :: dthetaPrime = 0.0

    ! ALPHA-parameter for fractional step. Not used for standard
    ! time stepping scheme.
    ! Only used internally in the Fractional Step scheme.
    ! Standard = 0.0 = no fractional step used.
    REAL(DP)                 :: dalpha = 0.0
    
    ! BETA-parameter for fractional step. Not used for standard
    ! time stepping scheme.
    ! Only used internally in the Fractional Step scheme.
    ! Standard = 0.0 = no fractional step used.
    REAL(DP)                 :: dbeta = 0.0
    
    ! Length of the next time step when using an equidistant time
    ! discretisation.
    ! Only used internally.
    REAL(DP)                 :: dtstepFixed
    
  END TYPE

!</typeblock>

!</types>

CONTAINS

  !****************************************************************************

!<function>

  PURE INTEGER FUNCTION timstp_getOrder (rtstepScheme) RESULT(iorder)
  
!<description>
  ! This routine analyses the time stepping structure rtstepScheme and returns
  ! the expected order of the time stepping algorithm that is described by
  ! that structure (e.g. 1 for explicit Euler, 2 for Crank Nicolson).
!</description>

!<input>
  ! Time stepping structure that describes an explicit time stepping scheme.
  TYPE(t_explicitTimeStepping), INTENT(IN) :: rtstepScheme
!</input>

!<result>
  ! Expected order of the time stepping algorithm described by rtstepScheme.
!</result>
  
!</function>

    SELECT CASE (rtstepScheme%ctimestepType)
    CASE (TSCHM_FRACTIONALSTEP)
      iorder = 2
    CASE (TSCHM_ONESTEP)
      IF (rtstepScheme%dthStep .EQ. 0.5_DP) THEN
        iorder = 2
      ELSE
        iorder = 1
      END IF
    CASE DEFAULT
      ! Error case
      iorder = 0
    END SELECT

  END FUNCTION

  !****************************************************************************
  
!<subroutine>

  SUBROUTINE timstp_init (rtstepScheme, ctimestepType, dtime, dtstep, dtheta)
  
!<description>
  ! Initialisation of the time stepping scheme. This routine initialises
  ! the structure rtstepScheme for the discretisation in time.
!</description>
  
!<input>
  ! The type of time stepping to use. TSCHM_ONESTEP for a one-step scheme
  ! or TSCHM_FRACTIONALSTEP for Fractional Step.
  INTEGER, INTENT(IN)    :: ctimestepType
  
  ! The initial simulational time. 
  REAL(DP), INTENT(IN)   :: dtime
    
  ! (Theoretical) Time step length of the simulation.
  ! Might vary from the real time step size by the type of the
  ! scheme, which is indicated by rtstepScheme%dtstep
  REAL(DP), INTENT(IN)   :: dtstep

  ! OPTIONAL: Theta scheme identifier. 
  ! Only used for ctimestepType=TSCHM_ONESTEP.
  !  =0.0: Forward Euler
  !  =1.0: Backward Euler
  !  =0.5: Crank Nicolson
  ! Ignored for Fractional step. If not specified, dtheta=1.0 (Backward Euler)
  ! is used in one-step schemes.
  REAL(DP), INTENT(IN), OPTIONAL   :: dtheta
!</input>

!<output>
  ! The time stepping structure. This is initialised to perform the first
  ! time step.
  TYPE(t_explicitTimeStepping), INTENT(OUT) :: rtstepScheme
!</output>

!</subroutine>

    REAL(DP) :: dtheta1,dthetp1,dalpha,dbeta

    dtheta1 = 1.0
    IF (PRESENT(dtheta)) dtheta1 = dtheta

    ! Standard initialisation of the structure.
    rtstepScheme%isubstep         = 1
    rtstepScheme%dcurrentTime     = dtime
    rtstepScheme%dtimeMacrostep   = dtime
    rtstepScheme%dtstepFixed      = dtstep

    ! Do we have standard time stepping or Fractional Step?
    
    IF (ctimestepType .NE. TSCHM_FRACTIONALSTEP) THEN
      
      rtstepScheme%ctimestepType    = TSCHM_ONESTEP
      
      ! Standard time stepping. Here the parameters are a little bit  
      ! easier to initialize.
      rtstepScheme%nsubsteps        = 1
      
      rtstepScheme%dtstep           = dtstep
      rtstepScheme%dtheta           = dtheta1
      rtstepScheme%dthStep          = dtstep * dtheta1
      
      ! Initialize weights for matrices and RHS:   
      rtstepScheme%dweightMatrixLHS = dtstep * dtheta1
      rtstepScheme%dweightMatrixRHS = dtstep * (dtheta1 - 1.0_DP)
      rtstepScheme%dweightNewRHS    = dtstep * dtheta1
      rtstepScheme%dweightOldRHS    = dtstep * (dtheta1 - 1.0_DP)
      rtstepScheme%dweightStationaryRHS = dtstep
      
    ELSE
    
      rtstepScheme%ctimestepType    = TSCHM_FRACTIONALSTEP
      
      ! The FS Theta-Scheme uses by theory 4 parameters:           
      !                                                            
      !   Theta   = 1 - sqrt(2) / 2                                
      !   Theta'  = 1 - 2 * Theta                                  
      !   alpha   = ( 1 - 2 * Theta ) / ( 1 - Theta )              
      !   beta    = 1 - alpha                                      
      !                                                            
      ! The parameter THETA in the DAT-file is ignored and replaced
      ! by a hard-coded setting.                                   
      
      dtheta1 = 1.0_DP-sqrt(0.5_DP)      
      dthetp1 = 1.0_DP-2.0_DP*dtheta1       
      dalpha  = dthetp1 / (1.0_DP-dtheta1)
      dbeta   = dtheta1 / (1.0_DP-dtheta1)

      rtstepScheme%dtheta           = dtheta1
      rtstepScheme%dthetaPrime      = dthetp1
      rtstepScheme%dalpha           = dalpha
      rtstepScheme%dbeta            = dbeta
      
      ! Initialise for 1st substep

      rtstepScheme%dtstep           =  3.0_DP * dtstep * dtheta1
      rtstepScheme%dweightMatrixLHS =  3.0_DP * dtstep * dalpha * dtheta1
      rtstepScheme%dweightMatrixRHS = -3.0_DP * dtstep * dbeta * dtheta1
      rtstepScheme%dweightNewRHS    =  3.0_DP * dtstep * dalpha * dtheta1
      rtstepScheme%dweightOldRHS    =  3.0_DP * dtstep * dbeta * dtheta1
      rtstepScheme%dweightStationaryRHS = 3.0_DP * dtstep * dtheta1
      
    END IF
    
  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>

  SUBROUTINE timstp_nextSubstep (rtstepScheme)
  
!<description>
  ! Advances in time. In rtstepScheme, the current simulation time is increased
  ! to the next point in time. The weights for the terms in the differential
  ! equation are updated according to the next substep in the time stepping scheme.
!</description>
  
!<inputoutput>
  ! The time stepping structure. Is modified to represent the next time step.
  TYPE(t_explicitTimeStepping), INTENT(INOUT) :: rtstepScheme
!</inputoutput>

!</subroutine>

    REAL(DP) :: dtstep,dtheta1,dthetp1,dalpha ,dbeta

    IF (rtstepScheme%ctimestepType .LT. 0) THEN
      CALL output_line ('timstp_nextSubstep: Time stepping structure not initialised!')
    END IF

    ! Increase the simulation time
    rtstepScheme%dcurrentTime = rtstepScheme%dcurrentTime + rtstepScheme%dtstep
    
    ! Increase number of current substep
    rtstepScheme%isubstep = rtstepScheme%isubstep + 1
    
    IF (rtstepScheme%isubstep .GT. rtstepScheme%nsubsteps) THEN
      rtstepScheme%isubstep = 1
      rtstepScheme%dtimeMacrostep = rtstepScheme%dcurrentTime
    END IF

    IF (rtstepScheme%ctimestepType .EQ. TSCHM_FRACTIONALSTEP) THEN
    
      ! In case of Fractional step, we have to modify the length of the
      ! current time step according to the substep:
      !
      ! For fractional step the handling of the time step size dtstep
      ! is slightly different than for a 1-step scheme.                                   
      ! There we are orienting on the length of the macrostep of    
      ! step length 3*dtstep and break up that into three different  
      ! substeps at different points in time, not corresponding to  
      ! dtstep. Depending on the number of the current substep, we   
      ! have two settings for the weights and time length:          
      
      dtheta1 = rtstepScheme%dtheta           
      dthetp1 = rtstepScheme%dthetaPrime      
      dalpha  = rtstepScheme%dalpha           
      dbeta   = rtstepScheme%dbeta    
      
      dtstep  = rtstepScheme%dtstepFixed

      IF (rtstepScheme%isubstep .NE. 2) THEN     

        ! 1st and 3rd substep
        
        rtstepScheme%dtstep           =  3.0_DP * dtstep * dtheta1
        
        rtstepScheme%dweightMatrixLHS =  3.0_DP * dtstep * dalpha * dtheta1
        rtstepScheme%dweightMatrixRHS = -3.0_DP * dtstep * dbeta * dtheta1
        rtstepScheme%dweightNewRHS    =  3.0_DP * dtstep * dalpha * dtheta1
        rtstepScheme%dweightOldRHS    =  3.0_DP * dtstep * dbeta * dtheta1
        rtstepScheme%dweightStationaryRHS = 3.0_DP * dtstep * dtheta1

      ELSE
      
        ! 2nd substep

        rtstepScheme%dtstep           =  3.0_DP * dtstep * dthetp1
        
        rtstepScheme%dweightMatrixLHS =  3.0_DP * dtstep * dalpha * dtheta1
        rtstepScheme%dweightMatrixRHS = -3.0_DP * dtstep * dbeta * dthetp1
        rtstepScheme%dweightNewRHS    =  3.0_DP * dtstep * dalpha * dthetp1
        rtstepScheme%dweightOldRHS    =  3.0_DP * dtstep * dbeta * dthetp1
        rtstepScheme%dweightStationaryRHS = 3.0_DP * dtstep * dthetp1
      
      END IF
      
    END IF

  END SUBROUTINE

END MODULE
