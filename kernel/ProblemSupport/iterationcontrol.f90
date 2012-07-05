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
!#  1.) itc_initResidual
!#      -> Initialises an iteration control structure.
!#
!#  2.) itc_pushResidual
!#      -> Pushes a new residual and updates the iteration control structure.
!#
!#  3.) itc_calcConvRateReal
!#      -> Calculates and returns the real convergence rate.
!#
!#  4.) itc_calcConvRateAsymptotic
!#      -> Calculates and returns the asymptotic convergence rate.
!#
!#  5.) itc_printStatistics
!#      -> Prints iteration statistics like e.g. number of iterations or
!#         convergence rates.
!#
!# ****************************************************************************
!# Information regarding stopping criterions
!# -----------------------------------------
!# The main purpose of the iteration control is to manage stopping criterions
!# for any kind of user-defined iterative method. In addition to the obvious
!# option of stopping when a maximum number of allowed iterations has been
!# performed, the iteration control module offers three additional stopping
!# criterions: convergence-check, divergence-check and stagnation-check.
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
!# The stagnation checks whether the iteration has stagnated and therefore
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
!# </purpose>
!##############################################################################

module iterationcontrol

use fsystem
use genoutput

implicit none

private

public :: t_iterationControl
public :: itc_initResidual
public :: itc_pushResidual
public :: itc_calcConvRateReal
public :: itc_calcConvRateAsymptotic
public :: itc_printStatistics

!<constants>

!<constantblock description="Iteration control status codes">

  ! No stopping criterions fulfilled; continue iteration
  integer, parameter, public :: ITC_STATUS_CONTINUE     = 0

  ! Maximum number of allowed iterations performed
  integer, parameter, public :: ITC_STATUS_MAX_ITER     = 1

  ! Convergence criterion(s) fulfilled
  integer, parameter, public :: ITC_STATUS_CONVERGED    = 2

  ! Divergence criterion(s) fulfilled
  integer, parameter, public :: ITC_STATUS_DIVERGED     = 3

  ! Stagnation criterion fulfilled
  integer, parameter, public :: ITC_STATUS_STAGNATED    = 4

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

    ! <!-- OUTPUT PARAMETERS -->

    ! Output: current status code
    integer :: cstatus = ITC_STATUS_CONTINUE

    ! Output: Total number of iterations performed
    integer :: nIterations = 0

    ! Output: Initial residual
    real(DP) :: dresInitial = 0.0_DP

    ! Output: Final residual
    real(DP) :: dresFinal = 0.0_DP

    ! Residual queue
    real(DP), dimension(ITC_RES_QUEUE_SIZE) :: Dresiduals

  end type

!</typeblock>

!</types>

contains

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

!<input>
  ! The initial residual
  real(DP), intent(in) :: dres
!</input>

!</subroutine>

    ! store initial residual and reset anything necessary
    riter%nIterations = 0
    riter%dresInitial = dres
    riter%dresFinal = dres
    riter%Dresiduals(1) = dres
    riter%cstatus = ITC_STATUS_CONTINUE

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
  logical :: btolRel, btolAbs, bdivRel, bdivAbs, bstag, bconv, bdiv

    ! set current defect and increase iteration count
    riter%nIterations = riter%nIterations+1
    riter%dresFinal = dres
    riter%Dresiduals(imod(riter%nIterations, ITC_RES_QUEUE_SIZE)+1) = dres

    ! Check against absolute criterions
    btolAbs = (riter%dtolAbs .gt. 0.0_DP) .and. (dres .le. riter%dtolAbs)
    bdivAbs = (riter%ddivAbs .gt. 0.0_DP) .and. (dres .ge. riter%ddivAbs)

    ! Divide current residual by initial residual
    dx = dres / max(riter%dresInitial, SYS_EPSREAL_DP)

    ! Check against relative criterions
    btolRel = (riter%dtolRel .gt. 0.0_DP) .and. (dx .le. riter%dtolRel)
    bdivRel = (riter%ddivRel .gt. 0.0_DP) .and. (dx .ge. riter%ddivRel)

    ! If the number of minimal iterations was not reached, continue
    if(riter%nIterations .lt. riter%nminIterations) then

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
    if((riter%nstagIter .gt. 0) .and. (riter%nstagIter .lt. riter%nIterations)) then

      ! assume stagnation
      bstag = .true.

      ! loop over the last n iterations and check for stagnation
      n = min(riter%nstagIter, ITC_RES_QUEUE_SIZE-1)
      do i = 1, n

        ! calculate indices in residual queue
        j = imod(riter%nIterations + ITC_RES_QUEUE_SIZE - i    , ITC_RES_QUEUE_SIZE) + 1
        k = imod(riter%nIterations + ITC_RES_QUEUE_SIZE - i - 1, ITC_RES_QUEUE_SIZE) + 1

        ! check for stagnation
        bstag = bstag .and. (riter%Dresiduals(j) .ge. (riter%dstagRate * riter%Dresiduals(k)))

      end do

      ! Stagnated?
      if(bstag) then

        riter%cstatus = ITC_STATUS_STAGNATED
        return

      end if

    end if

    ! Check against maximum iterations
    if(riter%nIterations .ge. riter%nmaxIterations) then

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
    if(riter%nIterations .le. 0) return

    ! Calculate convergence rate
    dcr = (riter%dresFinal / riter%dresInitial) ** (1.0_DP / real(riter%nIterations,DP))

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
  integer, optional, intent(in) :: niter
!</input>

!</function>

  integer :: n, i, j

    ! Clear convergence rate
    dcr = 0.0_DP

    ! Test for anything that might blow up
    if(riter%dresInitial .le. SYS_EPSREAL_DP) return
    if(riter%dresFinal .le. 0.0_DP) return
    if(riter%nIterations .le. 0) return

    ! calculate indices
    n = 3
    if(present(niter)) n = min(max(1,niter), ITC_RES_QUEUE_SIZE-1)
    n = min(n, riter%nIterations-1) + 1

    i = imod(riter%nIterations  , ITC_RES_QUEUE_SIZE) + 1
    j = imod(riter%nIterations-n, ITC_RES_QUEUE_SIZE) + 1

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
    end select

    ! Print other statistics
    call output_line('Iterations ............: ' // trim(sys_siL(riter%nIterations,8)), &
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

end module