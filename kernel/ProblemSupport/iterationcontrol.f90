!##############################################################################
!# ****************************************************************************
!# <name> iterationcontrol </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a data structure as well as some routines which can
!# be used for iterative methods such as user-written non-linear solvers to
!# manage stopping criterions and maintain iteration statistics.
!#
!# The following routines can be found here:
!#
!#  1.) itc_initIteration
!#      -> Basic initialisation of an iteration control structure.
!#
!#  2.) itc_initResidual
!#      -> Initialises an iteration control structure.
!#
!#  3.) itc_pushResidual
!#      -> Pushes a new residual and updates the iteration control structure.
!#
!#  4.) itc_calcConvRateReal
!#      -> Calculates and returns the real convergence rate.
!#
!#  5.) itc_calcConvRateAsymptotic
!#      -> Calculates and returns the asymptotic convergence rate.
!#
!#  6.) itc_printStatistics
!#      -> Prints iteration statistics like e.g. number of iterations or
!#         convergence rates.
!#
!#  7.) itc_getParamsFromParlist
!#      -> Reads iteration parameters from a section in a parameter list.
!#
!#  8.) itc_copyStatistics
!#      -> Copies the statistics data to another iteration structure
!#
!# ****************************************************************************
!# Information regarding stopping criterions
!# -----------------------------------------
!# The main purpose of the iteration control is to manage stopping criterions
!# for any kind of user-defined iterative method. In addition to the obvious
!# option of stopping when a maximum number of allowed iterations has been
!# performed, the iteration control module offers three additional stopping
!# criterions: convergence-check, divergence-check, stagnation-check and
!# orbiting-check.
!#
!#
!# Convergence-Check
!# -----------------
!# The convergence check compares the residual against a user-defined relative
!# and/or absolute tolerance parameters, namely dtolAbs and dtolRel, to decide
!# whether the iteration can be interpreted as 'converged'.
!#
!# -> Absolute Tolerance: If at some point the residual drops below the dtolAbs
!#    constant, the iteration has fulfilled the 'absolute' convergence check.
!#
!# -> Relative Tolerance: If the residual drops below dtolRel multiplied by the
!#    inital residual (passed to the itc_initResidual routine), then the
!#    iteration has fulfilled the 'relative' convergence check.
!#
!# By specifying a corresponding ITC_STOP_MODE_* constant to the ctolMode
!# member of the control structure, the user can specify under which
!# circumstances the iteration is to be treated as 'converged'.
!#
!# There are five options available:
!#
!# -> ITC_STOP_MODE_NONE
!#    Do not check for convergence at all.
!#
!# -> ITC_STOP_MODE_REL
!#    Treat the iteration as converged if the relative convergence check is
!#    fulfilled.
!#
!# -> ITC_STOP_MODE_ABS
!#    Treat the iteration as converged if the absolute convergence check is
!#    fulfilled.
!#
!# -> ITC_STOP_MODE_BOTH
!#    Treat the iteration as converged if both the relative and absolute
!#    convergence checks are fulfilled.
!#
!# -> ITC_STOP_MODE_ONE_OF
!#    Treat the iteration as converged if at least one of the relative or
!#    absolute convergence checks are fulfilled.
!#
!#
!# Divergence-Check
!# ----------------
!# The divergence check compares the residual against user-defined relative
!# and/or absolute divergence bound parameters, namely ddivRel and ddivAbs,
!# to decide whether the iterations is to be interpreted as 'diverged'.
!# In analogy to the convergence checks, the user can specify a corresponding
!# mode to decide how divergence checks are to be performed.
!#
!#
!# Stagnation-Check
!# ----------------
!# The stagnation check tests whether the iteration has stagnated and therefore
!# the iteration can be cancelled since it is not making any progress anymore.
!# The stagnation check is configured by the 'dstagRate' and 'nstagIter'
!# parameters of the iteration control data structure. Let r_k denote the
!# residual in iteration k >= nstagIter, then the stagnation check is fulfilled
!# if
!#     for 0 <= i < nstagIter:    r_{k-i} >= dstagRate * r_{k-i-1},
!#
!# i.e. the iteration is treated as 'stagnated' if the residual did not
!# decrease by at least a user-defined factor for a specified number of
!# consecutive iterations.
!#
!#
!# Orbiting-Check
!# --------------
!# The orbiting check tests whether the iteration is caught in an orbit, i.e.
!# if the iterative method is alternating through a set of cluster points
!# instead of converging to a unique limit. The orbiting check is controlled
!# via three parameters, namely dorbitTol, nminPeriods and nmaxClusters.
!# nmaxClusters specifies the maximum number of cluster points that the
!# orbit may have - for nonlinear solvers, an orbit usually consists of two
!# cluster points. nminPeriods specifies how often a residual has to repeat
!# to identify the current residual as a cluster point. It is recommended to
!# set the minimal period to at least 2 to avoid that only two equal residuals
!# are accidently recognised as a period. The dorbitTol parameter specifies
!# the maximal relative difference between two residuals, which is used to
!# test whether two residuals are treated as equal.
!#
!#
!# Status-Codes
!# ------------
!# In each call of itc_pushResidual(), the iteration control structure is
!# update and a new status code is set as feedback to the caller indicating
!# what action has to be performed. The status code is stored in the 'cstatus'
!# member of the iteration control data structure and specifies whether any
!# of the stopping criterions has been fulfilled.
!#
!# There are five status codes:
!#
!# -> ITC_STATUS_CONTINUE
!#    Indicates that no stopping criterion has been fulfilled yet, therefore
!#    the iterative method can continue.
!#
!# -> ITC_STATUS_MAX_ITER
!#    Indicates that the maximum number of allowed iterations was performed
!#    without fulfilling any convergence, divergence or stagnation criterion.
!#
!# -> ITC_STATUS_CONVERGED
!#    Indicates that the convergence criterion has been fulfilled.
!#
!# -> ITC_STATUS_DIVERGED
!#    Indicates that the divergence criterion has been fulfilled.
!#
!# -> ITC_STATUS_STAGNATED
!#    Indicates that the stagnation criterion has been fulfilled.
!#
!# -> ITC_STATUS_ORBITING
!#    Indicates that the iteration is (probably) caught in an orbit.
!#
!# </purpose>
!##############################################################################

module iterationcontrol

use fsystem
use genoutput
use paramlist

implicit none

private

  public :: t_iterationControl
  public :: itc_initIteration
  public :: itc_initResidual
  public :: itc_pushResidual
  public :: itc_calcConvRateReal
  public :: itc_calcConvRateAsymptotic
  public :: itc_printStatistics
  public :: itc_getParamsFromParlist
  public :: itc_copyStatistics

!<constants>

!<constantblock description="Iteration control status codes">

  ! Undefined status, iteration not started/running.
  integer, parameter, public :: ITC_STATUS_UNDEFINED    = 0

  ! No stopping criterions fulfilled; continue iteration
  integer, parameter, public :: ITC_STATUS_CONTINUE     = 1

  ! Maximum number of allowed iterations performed
  integer, parameter, public :: ITC_STATUS_MAX_ITER     = 2

  ! Convergence criterion(s) fulfilled
  integer, parameter, public :: ITC_STATUS_CONVERGED    = 3

  ! Divergence criterion(s) fulfilled
  integer, parameter, public :: ITC_STATUS_DIVERGED     = 4

  ! Stagnation criterion fulfilled
  integer, parameter, public :: ITC_STATUS_STAGNATED    = 5

  ! Orbiting criterion fulfilled
  integer, parameter, public :: ITC_STATUS_ORBITING     = 6

!</constantblock>

!<constantblock description="Stopping modes">

  ! Disable stopping criterion
  integer, parameter, public :: ITC_STOP_MODE_NONE = 0

  ! Check only relative criterion
  integer, parameter, public :: ITC_STOP_MODE_REL = 1

  ! Check only absolute criterion
  integer, parameter, public :: ITC_STOP_MODE_ABS = 2

  ! Check for both absolute and relative criterions
  integer, parameter, public :: ITC_STOP_MODE_BOTH = 3

  ! Chheck for at least one of absolute or relative criterions
  integer, parameter, public :: ITC_STOP_MODE_ONE_OF = 4

!</constantblock>

!<constantblock description="Generic iteration control constants">

  ! Number of residuals to be stored in queue
  integer, parameter, public :: ITC_RES_QUEUE_SIZE = 32

!</constantblock>

!</constants>



!<types>

!<typeblock description="Iteration control data structure"

  type t_iterationControl

    ! <!-- INPUT PARAMETERS -->

    ! Input: Minimum number of iterations to be performed
    integer :: nminIterations = 0

    ! Input: Maximum number of iterations to be performed
    integer :: nmaxIterations = 10

    ! Input: Absolute stopping criterion
    real(DP) :: dtolAbs = 1E-5_DP

    ! Input: Relative stopping criterion
    real(DP) :: dtolRel = 1E-5_DP

    ! Input: Stopping mode for convergence
    integer :: ctolMode = ITC_STOP_MODE_REL

    ! Input: Absolute divergence criterion
    real(DP) :: ddivAbs = 1E+99_DP

    ! Input: Relative divergence criterion
    real(DP) :: ddivRel = 1E+10_DP

    ! Input: Stopping mode for divergence
    integer :: cdivMode = ITC_STOP_MODE_ONE_OF

    ! Input: Stagnation rate criterion
    real(DP) :: dstagRate = 1.0_DP

    ! Input: Minimum number of consecutive stagnated iterations required.
    ! Set to 0 to disable stagnation check.
    integer :: nstagIter = 0

    ! Input: Relative orbiting tolerance
    real(DP) :: dorbitTol = 1E-2_DP

    ! Input: Minimal orbit period
    ! Set to 0 to disable orbit check.
    integer :: nminPeriods = 0

    ! Input: Maximal cluster points in orbit
    integer :: nmaxClusters = 2
    
    ! Maximum number of residuals which are used fo the calculation
    ! of asymptotic convergence rates.
    integer :: nasympConvRateResiduals = 3

    ! <!-- OUTPUT PARAMETERS -->

    ! Output: current status code
    integer :: cstatus = ITC_STATUS_UNDEFINED

    ! Output: Total number of iterations performed
    integer :: niterations = 0

    ! Output: Initial residual
    real(DP) :: dresInitial = 0.0_DP

    ! Output: Final residual
    real(DP) :: dresFinal = 0.0_DP

    ! Residual queue.
    real(DP), dimension(ITC_RES_QUEUE_SIZE) :: Dresiduals

  end type

!</typeblock>

!</types>

contains

! *************************************************************************************************

!<subroutine>

  subroutine itc_initIteration(riter)

!<description>
  ! Basic initialisation of the structure. The number of iterations
  ! is reset to zero. All residuals are set to zero.
  ! The status is set to "continue".
!</description>

!<inputoutput>
  ! The iteration structure to be initialised.
  type(t_iterationControl), intent(inout) :: riter
!</inputoutput>

!<remark>
  ! Alternatively, the structure can be initialised by itc_initResidual
  ! which calls this routine automatically.
!</remark>

!</subroutine>

    ! store initial residual and reset anything necessary
    riter%niterations = 0
    riter%dresInitial = 0.0_DP
    riter%dresFinal = 0.0_DP
    riter%Dresiduals(1) = 0.0_DP
    riter%cstatus = ITC_STATUS_UNDEFINED

  end subroutine

! *************************************************************************************************

!<subroutine>

  subroutine itc_initResidual(riter, dres)

!<description>
  ! Initialises an iteration control data structure using the initial residual.
!</description>

!<inputoutput>
  ! The iteration structure to be initialised.
  type(t_iterationControl), intent(inout) :: riter
!</inputoutput>

!<remark>
  ! The routine internally calls itc_initIteration.
  ! Therefore, the structure can be initialised either with itc_initIteration
  ! or with itc_initResidual.
!</remark>

!<input>
  ! The initial residual
  real(DP), intent(in) :: dres
!</input>

!</subroutine>

    ! store initial residual and reset anything necessary
    call itc_initIteration(riter)
    riter%dresInitial = dres
    riter%dresFinal = dres
    riter%Dresiduals(1) = dres

    ! We're now going to check whether we have some absolute convergence or divergence criterions
    ! given. First of all, ensure that we do not have to perform a minimum number of iterations.
    if(riter%nminIterations .gt. 0) return

    ! Check for absolute convergence if desired
    if((riter%dtolAbs .gt. 0.0_DP) .and. (dres .le. riter%dtolAbs)) then

      ! Check the convergence mode
      select case(riter%ctolMode)
      case (ITC_STOP_MODE_ABS, ITC_STOP_MODE_ONE_OF)
        riter%cstatus = ITC_STATUS_CONVERGED
        return
      end select

    end if

    ! Check for absolute divergence if desired
    if((riter%ddivAbs .gt. 0.0_DP) .and. (dres .ge. riter%ddivAbs)) then

      ! Check the divergence mode
      select case(riter%cdivMode)
      case (ITC_STOP_MODE_ABS, ITC_STOP_MODE_ONE_OF)
        riter%cstatus = ITC_STATUS_DIVERGED
        return
      end select

    end if

  end subroutine

! *************************************************************************************************

!<subroutine>

  subroutine itc_pushResidual(riter, dres)

!<description>
  ! Pushes a new residual and updates the iteration control structure.
!</description>

!<inputoutput>
  ! The iteration structure to be updated.
  type(t_iterationControl), intent(inout) :: riter
!</inputoutput>

!<input>
  ! The new residual
  real(DP), intent(in) :: dres
!</input>

!</subroutine>

  real(DP) :: dx
  integer :: i, j, k, n, cmode
  logical :: btolRel, btolAbs, bdivRel, bdivAbs, bstag, bconv, bdiv, borbit

    ! set current defect and increase iteration count
    riter%niterations = riter%niterations+1
    riter%dresFinal = dres
    riter%Dresiduals(mod(riter%nIterations, ITC_RES_QUEUE_SIZE)+1) = dres

    ! Check against absolute criterions
    btolAbs = (riter%dtolAbs .gt. 0.0_DP) .and. (dres .le. riter%dtolAbs)
    bdivAbs = (riter%ddivAbs .gt. 0.0_DP) .and. (dres .ge. riter%ddivAbs)

    ! Divide current residual by initial residual
    dx = dres / max(riter%dresInitial, SYS_EPSREAL_DP)

    ! Check against relative criterions
    btolRel = (riter%dtolRel .gt. 0.0_DP) .and. (dx .le. riter%dtolRel)
    bdivRel = (riter%ddivRel .gt. 0.0_DP) .and. (dx .ge. riter%ddivRel)

    ! If the number of minimal iterations was not reached, continue
    if(riter%niterations .lt. riter%nminIterations) then

      riter%cstatus = ITC_STATUS_CONTINUE
      return

    end if

    ! Do we check for convergence?
    if(riter%ctolMode .ne. ITC_STOP_MODE_NONE) then

      if((riter%dtolRel .le. 0.0_DP) .and. (riter%dtolAbs .le. 0.0_DP)) then

        ! no stopping criterion
        cmode = ITC_STOP_MODE_NONE

      else if(riter%dtolRel .le. 0.0_DP) then

        ! No relative tolerance given
        select case(riter%ctolMode)
        case (ITC_STOP_MODE_ABS,ITC_STOP_MODE_BOTH,ITC_STOP_MODE_ONE_OF)
          cmode = ITC_STOP_MODE_ABS
        case default
          cmode = ITC_STOP_MODE_NONE
        end select

      else if(riter%dtolAbs .le. 0.0_DP) then

        ! No absolute tolerance given
        select case(riter%ctolMode)
        case (ITC_STOP_MODE_REL,ITC_STOP_MODE_BOTH,ITC_STOP_MODE_ONE_OF)
          cmode = ITC_STOP_MODE_REL
        case default
          cmode = ITC_STOP_MODE_NONE
        end select

      else

        ! Both tolerances given
        cmode = riter%ctolMode

      end if

      ! Check for tolerances
      select case(cmode)
      case (ITC_STOP_MODE_NONE)
        bconv = .false.
      case (ITC_STOP_MODE_REL)
        bconv = btolRel
      case (ITC_STOP_MODE_ABS)
        bconv = btolAbs
      case (ITC_STOP_MODE_BOTH)
        bconv = btolRel .and. btolAbs
      case (ITC_STOP_MODE_ONE_OF)
        bconv = btolRel .or. btolAbs
      end select

      ! Converged?
      if(bconv) then

        riter%cstatus = ITC_STATUS_CONVERGED
        return

      end if

    end if

    ! Do we check for divergence?
    if(riter%cdivMode .ne. ITC_STOP_MODE_NONE) then

      if((riter%ddivRel .le. 0.0_DP) .and. (riter%ddivAbs .le. 0.0_DP)) then

        ! no stopping criterion
        cmode = ITC_STOP_MODE_NONE

      else if(riter%ddivRel .le. 0.0_DP) then

        ! No relative tolerance given
        select case(riter%cdivMode)
        case (ITC_STOP_MODE_ABS,ITC_STOP_MODE_BOTH,ITC_STOP_MODE_ONE_OF)
          cmode = ITC_STOP_MODE_ABS
        case default
          cmode = ITC_STOP_MODE_NONE
        end select

      else if(riter%ddivAbs .le. 0.0_DP) then

        ! No absolute tolerance given
        select case(riter%cdivMode)
        case (ITC_STOP_MODE_REL,ITC_STOP_MODE_BOTH,ITC_STOP_MODE_ONE_OF)
          cmode = ITC_STOP_MODE_REL
        case default
          cmode = ITC_STOP_MODE_NONE
        end select

      else

        ! Both tolerances given
        cmode = riter%cdivMode

      end if

      ! Check for tolerances
      select case(cmode)
      case (ITC_STOP_MODE_NONE)
        bdiv = .false.
      case (ITC_STOP_MODE_REL)
        bdiv = bdivRel
      case (ITC_STOP_MODE_ABS)
        bdiv = bdivAbs
      case (ITC_STOP_MODE_BOTH)
        bdiv = bdivRel .and. bdivAbs
      case (ITC_STOP_MODE_ONE_OF)
        bdiv = bdivRel .or. bdivAbs
      end select

      ! Diverged?
      if(bdiv) then

        riter%cstatus = ITC_STATUS_CONVERGED
        return

      end if

    end if

    ! Check for stagnation?
    if((riter%nstagIter .gt. 0) .and. (riter%nstagIter .lt. riter%niterations)) then

      ! assume stagnation
      bstag = .true.

      ! loop over the last n iterations and check for stagnation
      n = min(riter%nstagIter, ITC_RES_QUEUE_SIZE-1)
      do i = 1, n

        ! calculate indices in residual queue
        j = mod(riter%niterations + ITC_RES_QUEUE_SIZE - i    , &
                ITC_RES_QUEUE_SIZE) + 1
        k = mod(riter%niterations + ITC_RES_QUEUE_SIZE - i - 1, &
                ITC_RES_QUEUE_SIZE) + 1

        ! check for stagnation
        bstag = bstag .and. (riter%Dresiduals(j) .ge. (riter%dstagRate * riter%Dresiduals(k)))

      end do

      ! Stagnated?
      if(bstag) then

        riter%cstatus = ITC_STATUS_STAGNATED
        return

      end if

    end if

    ! Check for orbiting?
    if((riter%nminPeriods .gt. 0) .and. (riter%nminPeriods .lt. riter%niterations)) then

      ! loop over all allowed cluster lengths
      do i = 2, riter%nmaxClusters

        ! assume no orbit
        borbit = .false.

        ! loop over the minimum period length
        do j = 1, riter%nminPeriods

          ! calculate offset from current iteration
          k = j*i - 1
          if(k .ge. min(riter%niterations,ITC_RES_QUEUE_SIZE)) then
            borbit = .false.
            exit
          end if
          k = mod(riter%niterations+ITC_RES_QUEUE_SIZE-k, ITC_RES_QUEUE_SIZE)

          ! check difference
          dx = (riter%Dresiduals(k) - riter%dresFinal)
          if(abs(riter%Dresiduals(k) - riter%dresFinal) .gt. (riter%dorbitTol*riter%dresFinal)) then
            borbit = .false.
            exit
          end if

          ! If we come out here, we have a satellite candidate...
          borbit = .true.

        end do

        ! Orbiting?
        if(borbit) then

          riter%cstatus = ITC_STATUS_ORBITING
          return

        end if

      end do

    end if

    ! Check against maximum iterations
    if(riter%niterations .ge. riter%nmaxIterations) then

      ! maximum number of iterations fulfilled
      riter%cstatus = ITC_STATUS_MAX_ITER
      return

    end if

    ! continue iteration
    riter%cstatus = ITC_STATUS_CONTINUE

  end subroutine

! *************************************************************************************************

!<function>

  real(DP) function itc_calcConvRateReal(riter) result(dcr)

!<description>
  ! Calculates the real convergence rate.
!</description>

!<result>
  ! The real convergence rate.
!</result>

!<input>
  ! The iteration control data structure.
  type(t_iterationControl), intent(in) :: riter
!</input>

!</function>

    ! Clear convergence rate
    dcr = 0.0_DP

    ! Test for anything that might blow up
    if(riter%dresInitial .le. SYS_EPSREAL_DP) return
    if(riter%dresFinal .le. 0.0_DP) return
    if(riter%niterations .le. 0) return

    ! Calculate convergence rate
    dcr = (riter%dresFinal / riter%dresInitial) ** (1.0_DP / real(riter%niterations,DP))

  end function

! *************************************************************************************************

!<function>

  real(DP) function itc_calcConvRateAsymptotic(riter, niter) result(dcr)

!<description>
  ! Calculates the asymptotic convergence rate.
!</description>

!<result>
  ! The asymptotic convergence rate.
!</result>

!<input>
  ! The iteration control data structure.
  type(t_iterationControl), intent(in) :: riter

  ! Optional: The number of last iterations from which the convergence rate is to be computed.
  ! If not specified, riter%nasympConvRateResiduals is used.
  integer, optional, intent(in) :: niter
!</input>

!</function>

  integer :: n, i, j

    ! Clear convergence rate
    dcr = 0.0_DP

    ! Test for anything that might blow up
    if(riter%dresInitial .le. SYS_EPSREAL_DP) return
    if(riter%dresFinal .le. 0.0_DP) return
    if(riter%niterations .le. 0) return

    ! calculate indices
    n = riter%nasympConvRateResiduals
    if(present(niter)) n = min(max(1,niter), ITC_RES_QUEUE_SIZE-1)
    n = min(n, riter%nIterations-1) + 1

    i = mod(riter%nIterations  , ITC_RES_QUEUE_SIZE) + 1
    j = mod(riter%nIterations-n, ITC_RES_QUEUE_SIZE) + 1

    ! Calculate convergence rate
    dcr = (riter%Dresiduals(i) / riter%Dresiduals(j)) ** (1.0_DP / real(n,DP))

  end function

! *************************************************************************************************

!<subroutine>

  subroutine itc_printStatistics(riter, coutputMode)

!<description>
  ! Prints statistics of the iteration.
!</description>

!<input>
  ! The iteration control data structure.
  type(t_iterationControl), intent(in) :: riter

  ! Optional: The output mode for the statistics, one of the OU_MODE_* constants defined in
  ! genoutput.f90. If not given, OU_MODE_STD is used.
  integer(I32), optional :: coutputMode
!</input>

!</subroutine>

  ! local variables
  integer(I32) :: coutput
  real(DP) :: drcr, dacr, dimp

    ! select output mode
    coutput = OU_MODE_STD
    if(present(coutputMode)) coutput = coutputMode

    ! calculate convergence rates
    drcr = itc_calcConvRateReal(riter)
    dacr = itc_calcConvRateAsymptotic(riter)

    ! calculate residual improvement
    if(riter%dresInitial .gt. SYS_EPSREAL_DP) then
      dimp = riter%dresFinal / riter%dresInitial
    else
      dimp = 0.0_DP
    end if

    ! Print result - but only if it is not ITC_STATUS_CONTINUE
    select case(riter%cstatus)
    case (ITC_STATUS_CONVERGED)
      call output_line('Status ................: converged', &
        OU_CLASS_MSG, coutput, 'itc_printStatistics')
    case (ITC_STATUS_DIVERGED)
      call output_line('Status ................: diverged', &
        OU_CLASS_MSG, coutput, 'itc_printStatistics')
    case (ITC_STATUS_MAX_ITER)
      call output_line('Status ................: maximum iterations reached', &
        OU_CLASS_MSG, coutput, 'itc_printStatistics')
    case (ITC_STATUS_STAGNATED)
      call output_line('Status ................: stagnated', &
        OU_CLASS_MSG, coutput, 'itc_printStatistics')
    case (ITC_STATUS_ORBITING)
      call output_line('Status ................: orbiting', &
        OU_CLASS_MSG, coutput, 'itc_printStatistics')
    end select

    ! Print other statistics
    call output_line('Iterations ............: ' // trim(sys_siL(riter%niterations,8)), &
        OU_CLASS_MSG, coutput, 'itc_printStatistics')
    call output_line('Final Residual ........: ' // trim(sys_sdEL(riter%dresFinal,12)), &
        OU_CLASS_MSG, coutput, 'itc_printStatistics')
    call output_line('Initial Residual ......: ' // trim(sys_sdEL(riter%dresInitial,12)), &
        OU_CLASS_MSG, coutput, 'itc_printStatistics')
    call output_line('Residual Improvement ..: ' // trim(sys_sdEL(dimp,12)), &
        OU_CLASS_MSG, coutput, 'itc_printStatistics')
    call output_line('Conv-Rate (Real) ......: ' // trim(sys_sdEL(drcr,12)), &
        OU_CLASS_MSG, coutput, 'itc_printStatistics')
    call output_line('Conv-Rate (Asympt) ....: ' // trim(sys_sdEL(dacr,12)), &
        OU_CLASS_MSG, coutput, 'itc_printStatistics')

  end subroutine

! *************************************************************************************************

!<subroutine>

  subroutine itc_getParamsFromParlist(riter,ssection,rparamList)
  
!<description>
  ! Reads iteration control parameters from a section in the parameter list rparamList.
!</description>
  
!<input>
  ! Parameter list with the parameters configuring the iteration.
  type(t_parlist), intent(in) :: rparamList

  ! Name of the section in the parameter list containing the parameters.
  character(LEN=*), intent(in) :: ssection
!</input>

!<output>
  ! The iteration control data structure which receives the parameters.
  type(t_iterationControl), intent(inout) :: riter
!</output>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    type(t_iterationControl) :: riterDefault

    ! Get the section, so we can access it faster.
    call parlst_querysection(rparamList, ssection, p_rsection)

    if (.not. associated(p_rsection)) then
      call output_line ("Cannot read parameters; no section ''"//&
          trim(ssection)//"''!", &
          OU_CLASS_ERROR,OU_MODE_STD,"itc_getParamsFromParlist")
      call sys_halt()
    end if

    ! Read the parameters
    call parlst_getvalue_int (p_rsection, "nminIterations", &
        riter%nminIterations,riterDefault%nminIterations)

    call parlst_getvalue_int (p_rsection, "nmaxIterations", &
        riter%nmaxIterations,riterDefault%nmaxIterations)

    call parlst_getvalue_double (p_rsection, "dtolAbs", &
        riter%dtolAbs,riterDefault%dtolAbs)

    call parlst_getvalue_double (p_rsection, "dtolRel", &
        riter%dtolRel,riterDefault%dtolRel)

    call parlst_getvalue_int (p_rsection, "ctolMode", &
        riter%ctolMode,riterDefault%ctolMode)

    call parlst_getvalue_double (p_rsection, "ddivAbs", &
        riter%ddivAbs,riterDefault%ddivAbs)

    call parlst_getvalue_double (p_rsection, "ddivRel", &
        riter%ddivRel,riterDefault%ddivRel)

    call parlst_getvalue_int (p_rsection, "cdivMode", &
        riter%cdivMode,riterDefault%cdivMode)

    call parlst_getvalue_double (p_rsection, "dstagRate", &
        riter%dstagRate,riterDefault%dstagRate)

    call parlst_getvalue_int (p_rsection, "nstagIter", &
        riter%nstagIter,riterDefault%nstagIter)

    call parlst_getvalue_int (p_rsection, "nasympConvRateResiduals", &
        riter%nasympConvRateResiduals,riterDefault%nasympConvRateResiduals)
        
    riter%nasympConvRateResiduals = &
        max(1,min(riter%nasympConvRateResiduals,ITC_RES_QUEUE_SIZE))

  end subroutine

! *************************************************************************************************

!<subroutine>

  subroutine itc_copyStatistics(riterSource,riterDest)
  
!<description>
  ! Copies the statistics data from riterSource to riterDest.
  ! Can be used, e.g., if the statistics from a subsolver should be
  ! passed to the solver.
!</description>
  
!<input>
  ! The iteration control data structure to be copied.
  type(t_iterationControl), intent(inout) :: riterSource
!</input>

!<output>
  ! Destination structure.
  type(t_iterationControl), intent(inout) :: riterDest
!</output>

!</subroutine>

    ! Copy statistics data
    riterDest%cstatus       = riterSource%cstatus      
    riterDest%niterations   = riterSource%niterations  
    riterDest%dresInitial   = riterSource%dresInitial  
    riterDest%dresFinal     = riterSource%dresFinal    
    riterDest%Dresiduals(:) = riterSource%Dresiduals(:)

  end subroutine

end module
