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

module timediscretisation

  use fsystem
  use genoutput

  implicit none

!<constants>

!<constantblock description="Time discretisation type identifiers">

  ! Theta-scheme (implicit/explicit Euler, Crank Nicolson)
  integer, parameter :: TDISCR_THETA           = 0
  
  ! The dG(0) scheme.
  integer, parameter :: TDISCR_DG0             = 1

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! This block realises a time discretisation structure which contains all
  ! parameters that are necessary for the time discretisation of an ODE.
  type t_timeDiscretisation
  
    ! A TDISCR_xxxx-Flag that identifies the time discretisation scheme.
    integer :: ctype = TDISCR_THETA
    
    ! Number of time intervals
    integer :: nintervals          = 0

    ! Absolute start time of the simulation
    real(DP) :: dtimeInit          = 0.0_DP     
    
    ! Maximum time of the simulation
    real(DP) :: dtimeMax           = 0.0_DP
    
    ! Time step length of the time discretisation
    real(DP) :: dtstep             = 0.0_DP

    ! If ctype=TDISCR_THETA: theta parameter that identifies the Theta scheme.
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

  subroutine tdiscr_initTheta (dstart, dend, nintervals, dtheta, rtimediscr)
  
!<description>
  ! Initialises a time discretisation structure for the discretisation with
  ! a one-step Theta scheme.
!</description>

!<input>
  ! Start point of the time interval.
  real(DP), intent(IN) :: dstart
  
  ! End point of the time interval
  real(DP), intent(IN) :: dend
  
  ! Number of subintervals
  integer, intent(IN) :: nintervals
  
  ! Theta scheme parameter
  !  =0.0: Forward Euler,
  !  =1.0: Backward Euler (Standard),
  !  =0.5: Crank Nicolson.
  real(DP), intent(IN) :: dtheta
!</input>

!<output>
  ! The time discrtisation structure to be initialised.
  type(t_timeDiscretisation), intent(OUT) :: rtimediscr
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
    rtimediscr%dtstep = (dend-dstart)/real(nintervals,DP)
    rtimediscr%dtheta = dtheta

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine tdiscr_initdG0(dstart, dend, nintervals, rtimediscr)
  
!<description>
  ! Initialises a time discretisation structure for the discretisation with
  ! dG(0) scheme.
!</description>

!<input>
  ! Start point of the time interval.
  real(DP), intent(IN) :: dstart
  
  ! End point of the time interval
  real(DP), intent(IN) :: dend
  
  ! Number of subintervals
  integer, intent(IN) :: nintervals
!</input>

!<output>
  ! The time discrtisation structure to be initialised.
  type(t_timeDiscretisation), intent(OUT) :: rtimediscr
!</output>

    ! Initialise the parameters of the structure.
    rtimediscr%ctype = TDISCR_DG0
    rtimediscr%dtimeInit = dstart
    rtimediscr%dtimeMax = dend
    rtimediscr%nintervals = nintervals
    rtimediscr%dtstep = (dend-dstart)/real(nintervals,DP)

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
  type(t_timeDiscretisation), intent(IN) :: rtimediscr
!</input>

!<result>
  ! Global number of equations on current time grid
!</result>

!</function>

    tdiscr_igetNDofGlob = 0
    
    ! Cancel if the structure is not initialised.
    if (rtimediscr%nintervals .eq. 0) return

    select case(rtimediscr%ctype)
    case (TDISCR_THETA)
      ! One step scheme. We have as many DOF's as intervals + 1.
      tdiscr_igetNDofGlob = rtimediscr%nintervals + 1
      return

    case (TDISCR_DG0)
      ! dG(0). The DOF's are in the interval midpoints.
      tdiscr_igetNDofGlob = rtimediscr%nintervals
      return

    end select
    
    call output_line ('Unsupported time discretisation.', &
                      OU_CLASS_ERROR,OU_MODE_STD,'tdiscr_igetNDofGlob')
    call sys_halt()
    
  end function

end module
