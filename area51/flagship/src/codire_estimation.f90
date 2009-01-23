!##############################################################################
!# ****************************************************************************
!# <name> codire_estimation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all routines which are required to perform
!# error estimation for a conservation law for a scalar variable
!#
!# The following routines are available:
!#
!# 1.) codire_prepareErrorEstimator
!#     -> prepare the error estimator for a given solution
!#
!# 2.) codire_performErrorEstimation
!#     -> extract scalar variable used for error estimation
!#
!# </purpose>
!##############################################################################

module codire_estimation

  use fsystem
  use genoutput
  use geometryaux
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use pprocerror
  use pprocgradients
  use spatialdiscretisation
  use statistics
  use stdoperators
  use storage
  use triangulation

  use codire_basic
  use errorestimation
  use problem


  implicit none

  private

  public :: codire_prepareErrorEstimator
  public :: codire_performErrorEstimation

contains

  !*****************************************************************************

!<subroutine>

  subroutine codire_prepareErrorEstimator(rproblemLevel, rsolution, rerrorEstimator)

!<description>
    ! This subroutine prepares the error estimator for a given solution
    ! vector and initializes the internal structures
!</description>

!<input>
    ! multigrid structure
    type(t_problemLevel), intent(IN), target :: rproblemLevel

    ! solution vector
    type(t_vectorBlock), intent(IN) :: rsolution
!</input>

!<inputoutput>
    ! error estimator
    type(t_errorEstimator), intent(INOUT) :: rerrorEstimator
!</inputoutput>
!</subroutine>


    ! Start time measurement for error estimation
    call stat_startTimer(rtimer_errorestimation, STAT_TIMERSHORT)

    ! Prepare error estimator
    select case(rerrorEstimator%ierrorestimator)

    case (ERREST_CSPR_FACE:ERREST_LIMAVR)
      
      ! Set error variable to solution
      call lsyssc_duplicateVector(rsolution%RvectorBlock(1),&
          rerrorEstimator%p_rgradientRecovery%rerrorVariable,&
          LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_SHARE)
      
      ! Set pointers to triangulation and boundary structure
      rerrorEstimator%p_rgradientRecovery%p_rtriangulation => rproblemLevel%rtriangulation
      rerrorEstimator%p_rgradientRecovery%p_rboundary      => rproblemLevel%p_rproblem%rboundary


    case (ERREST_FIRSTDIFF:ERREST_SECONDDIFF)

      ! Set error variable to solution
      call lsyssc_duplicateVector(rsolution%RvectorBlock(1),&
          rerrorEstimator%p_rdifferenceIndicator%rerrorVariable,&
          LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_SHARE)

      ! Set pointers to triangulation and boundary structure
      rerrorEstimator%p_rdifferenceIndicator%p_rtriangulation => rproblemLevel%rtriangulation

    case DEFAULT
      call output_line('Invalid type of error estimator!',&
          OU_CLASS_ERROR,OU_MODE_STD,'codire_prepareErrorEstimator')
      call sys_halt()
    end select

    ! Stop time measurement for error estimation
    call stat_stopTimer(rtimer_errorestimation)

  end subroutine codire_prepareErrorEstimator

   !*****************************************************************************

!<subroutine>

  subroutine codire_performErrorEstimation(rerrorEstimator, rgridIndicator,&
                                        dprotectionLayerTolerance)

!<description>
    ! This subroutine performs error estimation based on the
    ! given error estimator and it computes the grid indicator
!</description>

!<input>
    ! tolerance for protection layer
    real(DP), intent(IN) :: dprotectionLayerTolerance
!</input>

!<inputoutput>
    ! error estimator
    type(t_errorEstimator), intent(INOUT) :: rerrorEstimator   
!</inputoutput>

!<output>
    ! grid indicator for adaptation
    type(t_vectorScalar), intent(OUT) :: rgridIndicator
!</output>
!</subroutine>

    ! local variables
    real(DP) :: derrorH1

    ! Start time measurement for error estimation
    call stat_startTimer(rtimer_errorestimation, STAT_TIMERSHORT)

    ! What type of error estimator are we?
    select case(rerrorEstimator%ierrorestimator)

    case (ERREST_CSPR_FACE:ERREST_LIMAVR)
      ! Recovery based error estimation of the H1-error
      call errest_calcH1Error(rerrorEstimator%p_rgradientRecovery,&
                              rgridIndicator, derrorH1)

      ! Construct the grid indicator based on the specified strategy
      call errest_calcGridIndicator(rerrorEstimator, rgridIndicator, derrorH1)
      

    case (ERREST_FIRSTDIFF)
      ! First-difference indicator
      call errest_calcFirstDiffIndicator(&
          rerrorEstimator%p_rdifferenceIndicator, rgridIndicator)

      ! Construct the grid indicator based on the specified strategy
      call errest_calcGridIndicator(rerrorEstimator, rgridIndicator)
      

    case (ERREST_SECONDDIFF)
      ! Second-difference indicator
      call errest_calcSecondDiffIndicator(&
          rerrorEstimator%p_rdifferenceIndicator, rgridIndicator)
      
      ! Construct the grid indicator based on the specified strategy
      call errest_calcGridIndicator(rerrorEstimator, rgridIndicator)

      
    case DEFAULT
      call output_line('Invalid type of error estimator!',&
          OU_CLASS_ERROR,OU_MODE_STD,'codire_performErrorEstimation')
      call sys_halt()
    end select


    ! Calculate protection layers (if required)
    if (rerrorEstimator%nprotectlayers .gt. 0)&
        call errest_calcProtectionLayers(rerrorEstimator, rgridIndicator,&
                                         dprotectionLayerTolerance)

    ! Stop time measurement for error estimation
    call stat_stopTimer(rtimer_errorestimation)

  end subroutine codire_performErrorEstimation
end module codire_estimation
