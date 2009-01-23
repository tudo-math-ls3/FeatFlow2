!##############################################################################
!# ****************************************************************************
!# <name> euler_estimation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all routines which are required to perform
!# error estimation for the compressible Euler/Navier-Stokes equations
!#
!# The following routines are available:
!#
!# 1.) euler_prepareErrorEstimator
!#     -> prepare the error estimator for a given solution
!#
!# 2.) euler_performErrorEstimation
!#     -> extract scalar variable used for error estimation
!#
!# </purpose>
!##############################################################################

module euler_estimation

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

  use errorestimation
  use euler_basic
  use problem

  implicit none

  private

  public :: euler_prepareErrorEstimator
  public :: euler_performErrorEstimation

contains

  !*****************************************************************************

!<subroutine>

  subroutine euler_prepareErrorEstimator(rproblemLevel, rsolution, rerrorEstimator)

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

    case (ERREST_CSPR_FACE: ERREST_LIMAVR)
      
      ! We need a scalar error variable
      call euler_calcScalarErrorVariable(rsolution, euler_getNVAR(rproblemLevel),&
          rerrorEstimator%p_rgradientRecovery%ierrorvariable,&
          rerrorEstimator%p_rgradientRecovery%rerrorVariable)
      
      ! Set pointers to triangulation and boundary structure
      rerrorEstimator%p_rgradientRecovery%p_rtriangulation => rproblemLevel%rtriangulation
      rerrorEstimator%p_rgradientRecovery%p_rboundary      => rproblemLevel%p_rproblem%rboundary


    case (ERREST_FIRSTDIFF:ERREST_SECONDDIFF)

      ! We need a scalar error variable
      call euler_calcScalarErrorVariable(rsolution, euler_getNVAR(rproblemLevel),&
          rerrorEstimator%p_rdifferenceIndicator%ierrorvariable,&
          rerrorEstimator%p_rdifferenceIndicator%rerrorVariable)

      ! Set pointers to triangulation and boundary structure
      rerrorEstimator%p_rdifferenceIndicator%p_rtriangulation => rproblemLevel%rtriangulation

    case DEFAULT
      call output_line('Invalid type of error estimator!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_prepareErrorEstimator')
      call sys_halt()
    end select

    ! Stop time measurement for error estimation
    call stat_stopTimer(rtimer_errorestimation)

  end subroutine euler_prepareErrorEstimator

  !*****************************************************************************

!<subroutine>

  subroutine euler_performErrorEstimation(rerrorEstimator, rgridIndicator,&
                                       dprotectionLayerTolerance)

!<description>
    ! This subroutine performs error estimation based on the given
    ! error estimator and it computes the grid indicator.
!</description>

!<input>
    ! tolerance for protection layers
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
                              rgridIndicator)

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
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_performErrorEstimation')
      call sys_halt()
    end select
    

    ! Calculate protection layers (if required)
    if (rerrorEstimator%nprotectlayers .gt. 0)&
        call errest_calcProtectionLayers(rerrorEstimator, rgridIndicator,&
                                         dprotectionLayerTolerance)

    ! Stop time measurement for error estimation
    call stat_stopTimer(rtimer_errorestimation)

  end subroutine euler_performErrorEstimation

  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************

  ! Here, the auxiliary routines follow

  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************

!<subroutine>

  subroutine euler_calcScalarErrorVariable(rsolution, NVAR, &
                                        ierrorvariable, rerrorVariable)

!<description>
    ! This subroutine calculates a scalar error variable from the global solution
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(IN) :: rsolution

    ! number of variables
    integer, intent(IN) :: NVAR

    ! variable for error estimator
    integer, intent(IN) :: ierrorvariable
!</input>

!<output>
    ! nodal error variable
    type(t_vectorScalar), intent(OUT) :: rerrorVariable
!</output>
!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rspatialdiscr
    real(DP), dimension(:), pointer :: p_Dsolution, p_DerrorVariable

    ! Set pointer
    p_rspatialdiscr => rsolution%RvectorBlock(1)%p_rspatialdiscr
    
    ! Set indicator variable
    select case(isystemFormat)
      
    case (SYSTEM_INTERLEAVEFORMAT)
      
      ! Create scalar vector for indicator variable
      call lsyssc_createVecByDiscr(p_rspatialdiscr, rerrorVariable,&
                                   .false., rsolution%cdataType)
      
      ! Set pointers
      call lsyssc_getbase_double(rerrorVariable, p_DerrorVariable)
      call lsysbl_getbase_double(rsolution,      p_Dsolution)
      
      ! Get indicator variable
      call euler_getVariableNodewise(rerrorVariable%NEQ, NVAR, p_Dsolution,&
                                  ierrorvariable, p_DerrorVariable)
      
    case (SYSTEM_BLOCKFORMAT)
      
      ! Create scalar vector for indicator variable
      call lsyssc_duplicateVector(rsolution%RvectorBlock(1), rerrorVariable,&
                                  LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_EMPTY)
      
      ! Set pointers
      call lsyssc_getbase_double(rerrorVariable, p_DerrorVariable)
      call lsysbl_getbase_double(rsolution,      p_Dsolution)
      
      ! Get indicator variables
      call euler_getVariableBlockwise(rerrorVariable%NEQ, NVAR, p_Dsolution,&
                                   ierrorvariable, p_DerrorVariable)

    case DEFAULT
      call output_line('Unsupproted system format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_calcScalarErrorVariable')
      call sys_halt()
    end select

  end subroutine euler_calcScalarErrorVariable
  
end module euler_estimation
