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
!# 1.) tdiscr_initOneStepTheta
!#     -> Initialise a time discretisation structure according to a one-step
!#        scheme.
!#
!# 2.) tdiscr_initFSTheta
!#     -> Initialise a time discretisation structure according to a 
!#        FS-theta scheme.
!#
!# 3.) tdiscr_initdG0
!#     -> Initialise a time discretisation structure according to the dG(0)
!#        scheme.
!#
!# 4.) tdiscr_igetNDofGlob
!#     -> Calculate the number of DOF's in time.
!#
!# 5.) tdiscr_done
!#     -> Releases a time discretisation structure.
!#
!# 6.) tdiscr_refineRegular
!#     -> Regular refinement of a time mesh.
!#
!# 7.) tdiscr_getTimestep
!#     -> Calculates the timestep for timestep interval iinterval.
!#
!# 8.) tdiscr_getTimestepWeights
!#     -> Calculate the weihts in a timestep.
!#
!# 9.) tdiscr_infoStatistics
!#     -> Print statistics about a time mesh.
!#
!# 10.) tdiscr_copy
!#      -> Copies a time discretisation structure.
!# </purpose>
!##############################################################################

module timediscretisation

  use fsystem
  use genoutput

  implicit none
  
  private
  
  public :: t_timeDiscretisation
  public :: tdiscr_initOneStepTheta
  public :: tdiscr_initFSTheta
  public :: tdiscr_initdG0
  public :: tdiscr_done
  public :: tdiscr_igetNDofGlob
  public :: tdiscr_refineRegular
  public :: tdiscr_getTimestep
  public :: tdiscr_getTimestepWeights
  public :: tdiscr_getOrder
  public :: tdiscr_infoStatistics
  public :: tdiscr_copy

!<constants>
!<constantblock description="Time discretisation type identifiers">

  ! Theta-scheme (implicit/explicit Euler, Crank Nicolson)
  integer, parameter, public :: TDISCR_ONESTEPTHETA    = 0
  
  ! Fractional step Theta, 3-step scheme
  integer, parameter, public :: TDISCR_FSTHETA         = 1

  ! The dG(0) scheme.
  integer, parameter, public :: TDISCR_DG0             = 2

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! This block realises a time discretisation structure which contains all
  ! parameters that are necessary for the time discretisation of an ODE.
  type t_timeDiscretisation
  
    ! A TDISCR_xxxx-Flag that identifies the time discretisation scheme.
    integer :: ctype = TDISCR_ONESTEPTHETA
    
    ! Number of time intervals
    integer :: nintervals          = 0

    ! Absolute start time of the simulation
    real(DP) :: dtimeInit          = 0.0_DP     
    
    ! Maximum time of the simulation
    real(DP) :: dtimeMax           = 0.0_DP
    
    ! Time step length of the time discretisation
    real(DP) :: dtstep             = 0.0_DP

    ! If ctype=TDISCR_ONESTEPTHETA: theta parameter that identifies the Theta scheme.
    !  =0.0: Forward Euler,
    !  =1.0: Backward Euler (Standard),
    !  =0.5: Crank Nicolson.
    real(DP) :: dtheta = 1.0_DP
  
  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine tdiscr_getOrder (rtimediscr,iorder)
  
!<description>
  ! Calculate the (theoretical) order of the time stepping scheme.
!</description>

!<input>
  ! The time discrtisation structure.
  type(t_timeDiscretisation), intent(in) :: rtimediscr
!</input>

!<output>
  ! Order.
  integer, intent(out) :: iorder
!</output>

    iorder = 1
    select case (rtimediscr%ctype)
    case (TDISCR_ONESTEPTHETA)
      if (rtimediscr%dtheta .eq. 0.5_DP) iorder = 2
    case (TDISCR_FSTHETA)
      iorder = 2
    end select

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine tdiscr_initOneStepTheta (rtimediscr, dstart, dend, nintervals, dtheta)
  
!<description>
  ! Initialises a time discretisation structure for the discretisation with
  ! a one-step Theta scheme.
!</description>

!<input>
  ! Start point of the time interval.
  real(DP), intent(in) :: dstart
  
  ! End point of the time interval
  real(DP), intent(in) :: dend
  
  ! Number of subintervals
  integer, intent(in) :: nintervals
  
  ! Theta scheme parameter
  !  =0.0: Forward Euler,
  !  =1.0: Backward Euler (Standard),
  !  =0.5: Crank Nicolson.
  real(DP), intent(in) :: dtheta
!</input>

!<output>
  ! The time discrtisation structure to be initialised.
  type(t_timeDiscretisation), intent(out) :: rtimediscr
!</output>

    ! Initialise the parameters of the structure.
    rtimediscr%ctype = TDISCR_ONESTEPTHETA
    rtimediscr%dtimeInit = dstart
    rtimediscr%dtimeMax = dend
    rtimediscr%nintervals = nintervals
    rtimediscr%dtstep = (dend-dstart)/real(nintervals,DP)
    rtimediscr%dtheta = dtheta

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine tdiscr_initFSTheta (rtimediscr, dstart, dend, nintervals)
  
!<description>
  ! Initialises a time discretisation structure for the discretisation with
  ! a 3-step FS-Theta scheme.
!</description>

!<input>
  ! Start point of the time interval.
  real(DP), intent(in) :: dstart
  
  ! End point of the time interval
  real(DP), intent(in) :: dend
  
  ! Number of subintervals; is rounded up to a multiple of 3.
  integer, intent(in) :: nintervals
!</input>

!<output>
  ! The time discrtisation structure to be initialised.
  type(t_timeDiscretisation), intent(out) :: rtimediscr
!</output>

    ! Initialise the parameters of the structure.
    rtimediscr%ctype = TDISCR_FSTHETA
    rtimediscr%dtimeInit = dstart
    rtimediscr%dtimeMax = dend
    rtimediscr%nintervals = 3*((nintervals+2)/3)
    rtimediscr%dtstep = 0.0_DP
    rtimediscr%dtheta = 0.0_DP

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine tdiscr_initdG0(rtimediscr, dstart, dend, nintervals)
  
!<description>
  ! Initialises a time discretisation structure for the discretisation with
  ! dG(0) scheme.
!</description>

!<input>
  ! Start point of the time interval.
  real(DP), intent(in) :: dstart
  
  ! End point of the time interval
  real(DP), intent(in) :: dend
  
  ! Number of subintervals
  integer, intent(in) :: nintervals
!</input>

!<output>
  ! The time discrtisation structure to be initialised.
  type(t_timeDiscretisation), intent(out) :: rtimediscr
!</output>

    ! Initialise the parameters of the structure.
    rtimediscr%ctype = TDISCR_DG0
    rtimediscr%dtimeInit = dstart
    rtimediscr%dtimeMax = dend
    rtimediscr%nintervals = nintervals
    rtimediscr%dtstep = (dend-dstart)/real(nintervals,DP)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine tdiscr_done(rtimediscr)
  
!<description>
  ! Releases a time discretisation structure
!</description>

!<inputoutput>
  ! The time discrtisation structure to be initialised.
  type(t_timeDiscretisation), intent(inout) :: rtimediscr
!</inputoutput>

    ! Initialise the parameters of the structure.
    rtimediscr%ctype = TDISCR_ONESTEPTHETA
    rtimediscr%dtimeInit = 0.0_DP
    rtimediscr%dtimeMax = 0.0_DP
    rtimediscr%nintervals = 0.0_DP
    rtimediscr%dtstep = 0.0_DP
    rtimediscr%dtheta = 1.0_DP

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine tdiscr_copy(rtimediscr,rtimeDiscrDest)
  
!<description>
  ! Copies a time discretisation.
!</description>

!<inputoutput>
  ! Source time discretisation
  type(t_timeDiscretisation), intent(in) :: rtimediscr
!</inputoutput>

!<inputoutput>
  ! The time discrtisation structure to be initialised.
  type(t_timeDiscretisation), intent(out) :: rtimediscrDest
!</inputoutput>
!</subroutine>

    ! As long as we do not do anything fancy (like allocating memory),
    ! we can do a flat copy.
    rtimeDiscrDest = rtimediscr

  end subroutine
  
  ! ***************************************************************************

!<function>  

  integer function tdiscr_igetNDofGlob(rtimediscr)

!<description>
  ! This function returns for a given time discretisation the number of 
  ! degrees of freedom in time.
!</description>

!<input>    
  ! The tiome discretisation structure that specifies the time discretisation.
  type(t_timeDiscretisation), intent(in) :: rtimediscr
!</input>

!<result>
  ! Global number of equations on current time grid
!</result>

!</function>

    tdiscr_igetNDofGlob = 0
    
    ! Cancel if the structure is not initialised.
    if (rtimediscr%nintervals .eq. 0) return

    select case(rtimediscr%ctype)
    case (TDISCR_ONESTEPTHETA,TDISCR_FSTHETA,TDISCR_DG0)
      ! One step scheme. We have as many DOF's as intervals + 1.
      tdiscr_igetNDofGlob = rtimediscr%nintervals + 1
      return

!    case (TDISCR_DG0)
!      ! dG(0). The DOF's are in the interval midpoints.
!      tdiscr_igetNDofGlob = rtimediscr%nintervals
!      return

    end select
    
    call output_line ('Unsupported time discretisation.', &
                      OU_CLASS_ERROR,OU_MODE_STD,'tdiscr_igetNDofGlob')
    call sys_halt()
    
  end function

  ! ***************************************************************************

!<subroutine>

  subroutine tdiscr_getTimestep(rtimediscr,iinterval,dtimeend,dtstep,dtimestart)
  
!<description>
  ! Calculates the timestep size for timestep interval iinterval.
!</description>

!<input>
  ! The time discrtisation structure to be initialised.
  type(t_timeDiscretisation), intent(in) :: rtimediscr
  
  ! Number of the subinterval to compute the length of.
  ! =1..nintervals, where iinterval=1 computes the lst
  ! interval from dtimeInit..dtimeInit+dtstep.
  ! A value < 1 returns the start time, a value > nintervals the end time
  ! of the time interval, both with time step length = 0.
  integer, intent(in) :: iinterval
!</input>

!<output>
  ! OPTIONAL: End time of the time interval
  real(dp), intent(out), optional :: dtimeend

  ! OPTIONAL: Timestep size for that step.
  real(dp), intent(out), optional :: dtstep

  ! OPTIONAL: Start time of the time interval
  real(dp), intent(out), optional :: dtimestart
!</output>

    integer :: isubstep
    real(DP) :: dtheta1,dthetp1
    real(DP) :: dlocalstart,dlocalend,dlocalstep

    if (iinterval .lt. 1) then
      dlocalstart = rtimediscr%dtimeinit
      dlocalend = rtimediscr%dtimeinit
      dlocalstep = 0.0_DP
    else if (iinterval .gt. rtimediscr%nintervals) then
      dlocalstart = rtimediscr%dtimeMax
      dlocalend = rtimediscr%dtimeMax
      dlocalstep = 0.0_DP
    else
      ! Type of the time stepping?
      select case (rtimediscr%ctype)
      case (TDISCR_ONESTEPTHETA,TDISCR_DG0)
        ! Constant timestep size
        dlocalstep = rtimediscr%dtstep
        dlocalstart = rtimediscr%dtimeinit + real(iinterval-1,dp) * dlocalstep
        dlocalend = dlocalstart + dlocalstep
        
      case (TDISCR_FSTHETA)
        ! That is slightly more involving. Step length depends on
        ! which of the three substep we have.
        isubstep = mod(iinterval-1-1,3)

        dtheta1 = 1.0_DP-sqrt(0.5_DP)
        dthetp1 = 1.0_DP-2.0_DP*dtheta1
        
        select case (isubstep)
        case (0)
          dlocalstep = 3.0_DP * rtimediscr%dtstep * dtheta1
          dlocalstart = rtimediscr%dtimeinit + &
              real(iinterval/3,dp) * 3.0_DP * rtimediscr%dtstep
        case (1)
          dlocalstep = 3.0_DP * rtimediscr%dtstep * dthetp1
          dlocalstart = rtimediscr%dtimeinit + &
              real(iinterval/3,dp) * 3.0_DP * rtimediscr%dtstep + dlocalstep
        case (2)
          dlocalstep = 3.0_DP * rtimediscr%dtstep * dtheta1
          dlocalstart = rtimediscr%dtimeinit + &
              real(1+(iinterval/3),dp) * 3.0_DP * rtimediscr%dtstep - dlocalstep
        end select
        
        dlocalend = dlocalstart + dlocalstep
      end select
    end if

    if (present(dtimestart)) &
        dtimestart = dlocalstart
    if (present(dtimeend)) &
        dtimeend = dlocalend
    if (present(dtstep)) &
        dtstep = dlocalstep

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine tdiscr_getTimestepWeights(rtimediscr,iinterval,dweightold,dweightnew)
  
!<description>
  ! Calculates the weights in the timestep scheme for the current interval.
  ! iinterval is allowed in the range 0..rtimediscr%nintervals. The weight
  ! computed here corresponds to the weights used to calculate u_(iinterval+1)
  ! from u_iinterval (so to say, how to come to the "end" of the interval). 
  ! The value iinterval=0 returns the weights used to calculate u_1.
!</description>

!<input>
  ! The time discrtisation structure to be initialised.
  type(t_timeDiscretisation), intent(in) :: rtimediscr
  
  ! Whether the timestepping is forward or backward in time.
  logical :: bforward
  
  ! Number of the subintervall to compute the weights for.
  integer, intent(in) :: iinterval
!</input>

!<output>
  ! Weight for the "old" solution.
  real(dp), intent(out) :: dweightold

  ! Weight for the "new" solution.
  real(dp), intent(out) :: dweightnew
!</output>

    integer :: isubstep
    real(DP) :: dtheta1,dthetp1,dalpha,dbeta
    real(DP) :: dlocalstart,dlocalend,dlocalstep

    ! We assume that the timestepping has the following form:
    !
    !  d/dt u = f(u)
    !
    ! For explicit timestepping methods, this formula is
    ! decomposed into
    !
    !  (u_n+1 - u_n) / delta_t = dweightold f(u_n)  +  dweightnew f(u_n+1)
    !
    ! We return the two weights here.

    ! Type of the time stepping?
    select case (rtimediscr%ctype)
    case (TDISCR_ONESTEPTHETA,TDISCR_DG0)
      ! Constant timestep size
      dweightnew = rtimeDiscr%dtheta
      dweightold = (1-rtimeDiscr%dtheta)
      
      if (iinterval .lt. 1) then
        dweightold = 0.0_DP
        dweightnew = 1.0_DP 
      end if
      if (iinterval .gt. rtimediscr%nintervals) then
        dweightold = 1.0_DP
        dweightnew = 0.0_DP
      end if
      
    case (TDISCR_FSTHETA)
      ! Fractional step.
      ! That is slightly more involving. Step length depends on
      ! which of the three substep we have.
      isubstep = mod(iinterval-1,3)

      dtheta1 = 1.0_DP-sqrt(0.5_DP)
      dthetp1 = 1.0_DP-2.0_DP*dtheta1
      dalpha  = dthetp1 / (1.0_DP-dtheta1)
      dbeta   = dtheta1 / (1.0_DP-dtheta1)
      
      select case (isubstep)
      case (0,2)
        dweightnew = 3.0_DP * dalpha * dtheta1
        dweightold  = 3.0_DP * dbeta * dtheta1
      case (1)
        dweightnew = 3.0_DP * dbeta * dthetp1
        dweightold  = 3.0_DP * dalpha * dthetp1
      end select
    
    end select
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine tdiscr_refineRegular (rtimediscr,nreflevels)

!<description>
  ! Refines a time mesh nreflevel times by subdividing every timestep.
  ! Note that for fractional step, every 3rd timestep is refined!
!</description>
 
!<input>
  ! Number of refinement levels
  integer, intent(in) :: nreflevels
!</input>

!<inputoutput>
  ! A new time discretisation structure level structure. Is replaced by a refined
  ! time discretisation.
  type(t_timeDiscretisation), intent(inout) :: rtimediscr
!</inputoutput>
  
!</subroutine>

    ! Cancel if there is no refinement.
    if (nreflevels .le. 0) return

    ! Type of timestep scheme?
    select case (rtimediscr%ctype)
    case (TDISCR_ONESTEPTHETA)
      ! Increase the number of iterations. Every intervall gets a new timestep
      ! in the midpoint.
      rtimediscr%nintervals = rtimediscr%nintervals * 2**nreflevels
      rtimediscr%dtstep = (rtimediscr%dtimeMax-rtimediscr%dtimeInit)/real(rtimediscr%nintervals,DP)

    case (TDISCR_FSTHETA)
      ! Fractional step. Get the number of macro-timesteps, refine and insert
      ! the substeps again.
      rtimediscr%nintervals = ((rtimediscr%nintervals/3) * 2**nreflevels) * 3
      
    end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine tdiscr_infoStatistics (rtimediscr,bheadline,ilevel)
  
!<description>
  ! Prints out statistical information of the given triangulation to the
  ! terminal. The output is formatted as a table with an optional headline
  ! and an optional level identifier.
!</description>

!<input>
  ! Time discretisation structure.
  type(t_timeDiscretisation), intent(in) :: rtimediscr
  
  ! OPTIONAL: Print out a headline above the statistical data.
  ! =FALSE: do not print = standard.
  ! =TRUE: print a headline.
  logical, intent(in), optional :: bheadline
  
  ! OPTIONAL: Level identifier.
  ! If specified, an additional column 'Level' is added to the front of
  ! the statistics table. ilevel is printed to this column.
  integer, intent(in), optional :: ilevel
!</input>

!</subroutine>

      if (present(bheadline)) then
        if (bheadline) then
          ! Print a headline
          if (present(ilevel)) then
            call output_line(&
              'Lv. type #intervals         start          stop        dtstep theta')
            call output_line(&
              '-------------------------------------------------------------------')
          else
            call output_line(&
              ' type #intervals         start          stop        dtstep theta')
            call output_line(&
              '----------------------------------------------------------------')
          end if
        end if
      end if

      ! Print out the statistics
      if (present(ilevel)) &
        call output_line (trim(sys_si(ilevel,3))//' ',bnolinebreak=.true.)
      
      call output_line (&
          trim(sys_si(rtimediscr%ctype,5)) &
        //trim(sys_si(rtimediscr%nintervals,11)) &
        //sys_adjustr(sys_sd(rtimediscr%dtimeInit,6),14) &
        //sys_adjustr(sys_sd(rtimediscr%dtimeMax,6),14) &
        //sys_adjustr(sys_sd(rtimediscr%dtstep,6),14) &
        //sys_adjustr(sys_sd(rtimediscr%dtheta,3),6) )

  end subroutine 

end module
