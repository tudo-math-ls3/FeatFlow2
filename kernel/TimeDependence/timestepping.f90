!##############################################################################
!# ****************************************************************************
!# <name> timestepping </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a realisation of 1D time-stepping schemes used 
!# for the time discretisation of PDE`s. Examples for this are Explicit Euler,
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
!# 3.) timstp_nextSubstepTime
!#     -> Increase the simulation time without updating the weights
!#
!# 4.) timstp_nextSubstepWeights
!#     -> Set the weights for the next substep without updating the time
!#
!# 5.) timstp_getOrder
!#     -> Retrieves the order of the time error that is expected from a
!#        defined time stepping algorithm.
!#
!# 6.) timstp_setBaseSteplength
!#     Modify the base step length of a time stepping scheme
!#
!# How to use this module?
!#
!# Well, that is simple. It is like this:
!#
!# <code>
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
!# </code>
!#
!# First initialise the time stepping scheme by timstp_init, then call
!# timstp_nextSubstep in a loop to proceed from time dcurrentTime
!# to time dcurrentTime+dtstep. That is all.
!#
!# The time stepping structure contains definitions for various constants that
!# are to be used in front of the different terms in the PDE. See the 
!# documentation of t_explicitTimeStepping which constants are to be used
!# for what.
!#
!# Note: The timstp_nextSubstep routine increases the simulation time and
!#   updates the weights in one step. In some simulations, this has to be
!#   decoupled, e.g. when some things on the right hand side of a system
!#   must be assembled at timestep $t^n$ and some on timestep $t^{n+1}$.
!#   In this case, one can use timstp_nextSubstepTime and 
!#   timstp_nextSubstepWeights. The usage is as follows:
!#
!# <code>
!#    TYPE(t_explicitTimeStepping) :: tstepping
!#
!#    ! Initialise time stepping  
!#    CALL timstp_init (tstepping,...)
!#
!#    DO istep = 1,#max number of time steps
!# 
!#      ... (do some assembly at the timestep $t^n$)
!#
!#      ! Increase the simulation time
!#      CALL timstp_nextSubstepTime (tstepping)
!#
!#      ... (do some assembly at the timestep $t^{n+1}$
!#
!#      ! Proceed to next time step. Increase step number and update weights.
!#      CALL timstp_nextSubstepWeights (tstepping)
!#
!#    END DO
!# </code>
!#
!# </purpose>
!##############################################################################

module timestepping

  use fsystem
  use genoutput

  implicit none
  
  private
  
!<constants>

!<constantblock description="Identifiers for time step schemes">

  ! Identifier for one step scheme (e.g. Explicit Euler)
  integer, parameter, public :: TSCHM_ONESTEP        = 0
  
  ! Identifier for Fractional Step scheme
  integer, parameter, public :: TSCHM_FRACTIONALSTEP = 1

!</constantblock>

!</constants>
  
  
!<types>

!<typeblock>

  ! Time stepping structure used in explicit time stepping schemes like
  ! explicit Euler, Crank Nicolson etc.
  ! The time stepping scheme is designed to discretise a rather general
  ! equation in time for a solution $u$ and a time variable $t$:
  !
  !               <tex> $$ u_t(x,t) + N(u(x,t)) = f(x,t) $$ </tex>
  !
  ! Given a time step $k$ and a time stepping scheme parameter $\theta$, 
  ! this is discretised in as:
  ! <tex>
  !
  !  $$ (u_{n+1} - u_n) / k  +  \theta N(u_{n+1})
  !
  !    = -(1-\theta) N(u_n)  +  \theta f_{n+1}  +  (1-\theta) f_n $$
  !
  ! </tex>
  ! Now, regroup this equation with $u_{n+1}$ on the LHS and $u_n$ on the
  ! RHS and give the coefficients in front of all terms an appropriate
  ! name. Then the time stepping scheme realised by this structure yields:
  ! <tex>
  !
  ! $$ u_{n+1} + dweightMatrixLHS*N(u_n+1) 
  !
  ! =  u_n + dweightMatrixRHS*N(u_n)  +  dweightNewRHS*f_{n+1}  +  dweightOldRHS*f_n $$
  !
  ! </tex>
  ! If the simulation has steady (in)homogenuous boundary conditions
  ! and a One-Step Theta-Scheme is used, the RHS build from f_{n+1}
  ! and f_n can be simplified. In this case:
  ! <tex>
  !  $$  dweightNewRHS*f_n+1  +  dweightOldRHS*f_n
  !    = dweightNewRHS*f_n    +  dweightOldRHS*f_n
  !    = dweightStationaryRHS*f_n                    $$
  ! </tex>

  type t_explicitTimeStepping
  
    ! Time step scheme identifier. One of the TSCHM_xxxx-constants.
    ! Usually TSCHM_ONESTEP for one step schemes
    integer                  :: ctimestepType = -1
    
    ! Current substep in time step scheme. Normally =1.
    ! For time stepping schemes with multiple substeps (Fractional Step),
    ! this counts the current substep 1..nsubsteps.
    integer                  :: isubstep
    
    ! Number of substeps of the time stepping scheme. Normally =1.
    ! For time stepping schemes with multiple substeps (Fractional Step),
    ! this is the number of substeps.
    integer                  :: nsubsteps
    
    ! Current simulation time, this structure represents
    real(DP)                 :: dcurrentTime
    
    ! Current simulation time of the macrostep. For standard $\theta$-schemes,
    ! there is dcurrentTime=dtimeMacrostep. For schemes with multiple substeps
    ! like Fractional-Step, dcurrentTime is the current time of the substep, 
    ! while dtimeMacrostep saves the 'last' simulation time stamp with 
    ! "isubstep=1", i.e. the beginning of the 'macrostep'.
    real(DP)                 :: dtimeMacrostep
    
    ! Length of the previous timestep. =0 if there was no previous timestep.
    real(DP)                 :: dtlaststep
    
    ! Length of the next time step
    real(DP)                 :: dtstep
    
    ! Main configuration parameter of the Theta scheme
    real(DP)                 :: dtheta
    
    ! Theta-scheme identifier for the substeps of the step.
    ! For 1-step schemes:
    !  =0: Forward Euler,
    !  =1: Backward Euler,
    !  =0.5: Crank Nicolson.
    ! Special values for Fractional Step.
    real(DP)                 :: dthStep
    
    ! Weight in front of the matrix on the LHS.
    real(DP)                 :: dweightMatrixLHS

    ! Weight in front of the matrix on the RHS.
    real(DP)                 :: dweightMatrixRHS
    
    ! Weight to be used in the current time step for the "new" RHS $f_{n+1}$.
    real(DP)                 :: dweightNewRHS

    ! Weight to be used in the current time step for the "old" RHS $f_{n}$.
    real(DP)                 :: dweightOldRHS
    
    ! Weight to be used in the current time step if a stationary RHS is used
    ! (instead of a combination of dweightNewRHS and dweightOldRHS).
    real(DP)                 :: dweightStationaryRHS
    
    ! "Adjungated" parameter Theta` of the Theta scheme.
    ! Only used internally in the Fractional Step scheme.
    ! Standard = 0.0 = no fractional step used.
    real(DP)                 :: dthetaPrime = 0.0

    ! ALPHA-parameter for fractional step. Not used for standard
    ! time stepping scheme.
    ! Only used internally in the Fractional Step scheme.
    ! Standard = 0.0 = no fractional step used.
    real(DP)                 :: dalpha = 0.0
    
    ! BETA-parameter for fractional step. Not used for standard
    ! time stepping scheme.
    ! Only used internally in the Fractional Step scheme.
    ! Standard = 0.0 = no fractional step used.
    real(DP)                 :: dbeta = 0.0
    
    ! Length of the next time step when using an equidistant time
    ! discretisation.
    ! Only used internally.
    real(DP)                 :: dtstepFixed
    
  end type
  
  public :: t_explicitTimeStepping

!</typeblock>

!</types>

  public :: timstp_init
  public :: timstp_nextSubstep
  public :: timstp_nextSubstepTime
  public :: timstp_nextSubstepWeights
  public :: timstp_getOrder
  public :: timstp_setBaseSteplength
  
contains

  !****************************************************************************

!<function>

  pure integer function timstp_getOrder (rtstepScheme) result(iorder)
  
!<description>
  ! This routine analyses the time stepping structure rtstepScheme and returns
  ! the expected order of the time stepping algorithm that is described by
  ! that structure (e.g. 1 for explicit Euler, 2 for Crank Nicolson).
!</description>

!<input>
  ! Time stepping structure that describes an explicit time stepping scheme.
  type(t_explicitTimeStepping), intent(in) :: rtstepScheme
!</input>

!<result>
  ! Expected order of the time stepping algorithm described by rtstepScheme.
!</result>
  
!</function>

    select case (rtstepScheme%ctimestepType)
    case (TSCHM_FRACTIONALSTEP)
      iorder = 2
    case (TSCHM_ONESTEP)
      if (rtstepScheme%dthStep .eq. 0.5_DP) then
        iorder = 2
      else
        iorder = 1
      end if
    case DEFAULT
      ! Error case
      iorder = 0
    end select

  end function

  !****************************************************************************
  
!<subroutine>

  subroutine timstp_init (rtstepScheme, ctimestepType, dtime, dtstep, dtheta)
  
!<description>
  ! Initialisation of the time stepping scheme. This routine initialises
  ! the structure rtstepScheme for the discretisation in time.
!</description>
  
!<input>
  ! The type of time stepping to use. TSCHM_ONESTEP for a one-step scheme
  ! or TSCHM_FRACTIONALSTEP for Fractional Step.
  integer, intent(in)    :: ctimestepType
  
  ! The initial simulational time. 
  real(DP), intent(in)   :: dtime
    
  ! (Theoretical) Time step length of the simulation / 'base' step length of the
  ! time stepping scheme.
  ! Might vary from the real time step size by the type of the
  ! scheme. The real time step size is found in rtstepScheme%dtstep.
  real(DP), intent(in)   :: dtstep

  ! OPTIONAL: Theta scheme identifier. 
  ! Only used for ctimestepType=TSCHM_ONESTEP.
  !  =0.0: Forward Euler
  !  =1.0: Backward Euler
  !  =0.5: Crank Nicolson
  ! Ignored for Fractional step. If not specified, dtheta=1.0 (Backward Euler)
  ! is used in one-step schemes.
  real(DP), intent(in), optional   :: dtheta
!</input>

!<output>
  ! The time stepping structure. This is initialised to perform the first
  ! time step.
  type(t_explicitTimeStepping), intent(out) :: rtstepScheme
!</output>

!</subroutine>

    real(DP) :: dtheta1,dthetp1,dalpha,dbeta

    dtheta1 = 1.0
    if (present(dtheta)) dtheta1 = dtheta

    ! Standard initialisation of the structure.
    rtstepScheme%isubstep         = 1
    rtstepScheme%dcurrentTime     = dtime
    rtstepScheme%dtimeMacrostep   = dtime
    rtstepScheme%dtlaststep       = 0.0_DP
    
    if (ctimestepType .ne. TSCHM_FRACTIONALSTEP) then
      
      rtstepScheme%ctimestepType    = TSCHM_ONESTEP
      
      ! Standard time stepping. Here the parameters are a little bit  
      ! easier to initialise.
      rtstepScheme%nsubsteps        = 1
      
      rtstepScheme%dtstep           = dtstep
      rtstepScheme%dtheta           = dtheta1
      rtstepScheme%dthStep          = dtstep * dtheta1
      
    else
    
      rtstepScheme%ctimestepType    = TSCHM_FRACTIONALSTEP

      ! The FS-scheme has 3 substeps.
      rtstepScheme%nsubsteps        = 3
      
      ! The FS Theta-Scheme uses by theory 4 parameters:           
      !                                                            
      !   Theta   = 1 - sqrt(2) / 2                                
      !   Theta`  = 1 - 2 * Theta                                  
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
      
    end if

    ! Initialise the weights for the first time step by calling
    ! timstp_setBaseSteplength with time step length = dtstep
    call timstp_setBaseSteplength (rtstepScheme, dtstep)
    
  end subroutine

  !****************************************************************************
  
!<subroutine>

  subroutine timstp_setBaseSteplength (rtstepScheme, dtstep)
  
!<description>
  ! Modifies the base step length of the time stepping scheme in rtstepScheme
  ! to dtstep. The base step length is the theoretical time step length
  ! if an equidistant time discretisation is used. It coincides with the actual
  ! time step length only if a one step scheme is used.
  ! Reinitialises all weights according to the current time step configuration
  ! in rtstepScheme.
!</description>
  
!<input>
  ! (Theoretical) Time step length of the simulation.
  ! Might vary from the real time step size by the type of the
  ! scheme, which is indicated by rtstepScheme%dtstep
  real(DP), intent(in)   :: dtstep
!</input>

!<output>
  ! The time stepping structure. This is initialised according to the
  ! new time step length.
  type(t_explicitTimeStepping), intent(inout) :: rtstepScheme
!</output>

!</subroutine>

    ! local variables
    real(DP) :: dtheta1,dthetp1,dalpha,dbeta
    
    ! Set the new step length
    rtstepScheme%dtstepFixed        = dtstep

    if (rtstepScheme%ctimestepType .ne. TSCHM_FRACTIONALSTEP) then
      ! Standard time stepping scheme.
      rtstepScheme%nsubsteps        = 1
      
      dtheta1                       = rtstepScheme%dtheta
      
      rtstepScheme%dtstep           = dtstep
      rtstepScheme%dthStep          = dtstep * dtheta1
      
      ! Initialise weights for matrices and RHS:   
      rtstepScheme%dweightMatrixLHS = dtstep * dtheta1
      rtstepScheme%dweightMatrixRHS = - dtstep * (1.0_DP - dtheta1)
      rtstepScheme%dweightNewRHS    = dtstep * dtheta1
      rtstepScheme%dweightOldRHS    = dtstep * (1.0_DP - dtheta1)
      rtstepScheme%dweightStationaryRHS = dtstep
    
    else

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
      
      if (rtstepScheme%isubstep .ne. 2) then     

        ! 1st and 3rd substep
        
        rtstepScheme%dtstep           =  3.0_DP * dtstep * dtheta1
        
        rtstepScheme%dweightMatrixLHS =  3.0_DP * dtstep * dalpha * dtheta1
        rtstepScheme%dweightMatrixRHS = -3.0_DP * dtstep * dbeta * dtheta1
        rtstepScheme%dweightNewRHS    =  3.0_DP * dtstep * dalpha * dtheta1
        rtstepScheme%dweightOldRHS    =  3.0_DP * dtstep * dbeta * dtheta1
        rtstepScheme%dweightStationaryRHS = 3.0_DP * dtstep * dtheta1

      else
      
        ! 2nd substep

        rtstepScheme%dtstep           =  3.0_DP * dtstep * dthetp1
        
        rtstepScheme%dweightMatrixLHS =  3.0_DP * dtstep * dalpha * dtheta1
        rtstepScheme%dweightMatrixRHS = -3.0_DP * dtstep * dalpha * dthetp1
        rtstepScheme%dweightNewRHS    =  3.0_DP * dtstep * dbeta * dthetp1
        rtstepScheme%dweightOldRHS    =  3.0_DP * dtstep * dalpha * dthetp1
        rtstepScheme%dweightStationaryRHS = 3.0_DP * dtstep * dthetp1
      
      end if
      
    end if
    
  end subroutine

  !****************************************************************************
  
!<subroutine>

  subroutine timstp_nextSubstep (rtstepScheme)
  
!<description>
  ! Advances in time. In rtstepScheme, the current simulation time is increased
  ! to the next point in time. The weights for the terms in the differential
  ! equation are updated according to the next substep in the time stepping scheme.
!</description>
  
!<inputoutput>
  ! The time stepping structure. Is modified to represent the next time step.
  type(t_explicitTimeStepping), intent(inout) :: rtstepScheme
!</inputoutput>

!</subroutine>

    ! Remember the last step length.
    rtstepScheme%dtlaststep = rtstepScheme%dtstep

    ! Update simulation time and weights, right after each other.
    call timstp_nextSubstepTime (rtstepScheme)
    call timstp_nextSubstepWeights (rtstepScheme)

  end subroutine

  !****************************************************************************
  
!<subroutine>

  subroutine timstp_nextSubstepTime (rtstepScheme)
  
!<description>
  ! Advances in time. In rtstepScheme, the current simulation time is increased
  ! to the next point in time.
  ! Note that the weights and the number of the current (sub)step are not 
  ! changed! These have to be updated with an additional call to
  ! timstp_nextSubstepWeights!
!</description>
  
!<inputoutput>
  ! The time stepping structure. Is modified to represent the next time step.
  type(t_explicitTimeStepping), intent(inout) :: rtstepScheme
!</inputoutput>

!</subroutine>

    if (rtstepScheme%ctimestepType .lt. 0) then
      call output_line ('timstp_nextSubstep: Time stepping structure not initialised!')
    end if

    ! Increase the simulation time
    rtstepScheme%dcurrentTime = rtstepScheme%dcurrentTime + rtstepScheme%dtstep
    
  end subroutine

  !****************************************************************************
  
!<subroutine>

  subroutine timstp_nextSubstepWeights (rtstepScheme)
  
!<description>
  ! Advances the weights of the time step scheme in time. The number of the 
  ! current time (sub)step in increased and the weights for the terms in the 
  ! differential equation are updated according to the next substep in 
  ! the time stepping scheme.
!</description>
  
!<inputoutput>
  ! The time stepping structure. Is modified to represent the next time step.
  type(t_explicitTimeStepping), intent(inout) :: rtstepScheme
!</inputoutput>

!</subroutine>

    if (rtstepScheme%ctimestepType .lt. 0) then
      call output_line ('timstp_nextSubstep: Time stepping structure not initialised!')
    end if

    ! Increase number of current substep
    rtstepScheme%isubstep = rtstepScheme%isubstep + 1
    
    if (rtstepScheme%isubstep .gt. rtstepScheme%nsubsteps) then
      rtstepScheme%isubstep = 1
      ! Set the time with a slightly different approach to prevent rounding errors.
      rtstepScheme%dtimeMacrostep = rtstepScheme%dtimeMacrostep + &
          real(rtstepScheme%nsubsteps,DP) * rtstepScheme%dtstepFixed
      rtstepScheme%dcurrentTime = rtstepScheme%dtimeMacrostep
    end if
    
    ! Initialise all weights by calling the timstp_setBaseSteplength routine
    ! with time step length = current time step length, saved in the dtstepFixed
    ! variable of rtstepScheme.
    call timstp_setBaseSteplength (rtstepScheme, rtstepScheme%dtstepFixed)
    
  end subroutine

end module
