!##############################################################################
!# ****************************************************************************
!# <name> timestep </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic data structures and subroutines for
!# handling time-dependent problems. A calling program should create a
!# new object of time t_timestep and fill the required parameters
!# either by reading a user- supplied parameter file of via direct
!# adjustment.
!#
!# Furthermore, each object of type t_solver provides information
!# about convergence and norms after the solver has terminated.
!#
!# The following routines are available:
!#
!# 1.) tstep_createTimestep = tstep_createTimestepDirect /
!#                            tstep_createTimestepIndirect
!#     -> Creates a new time-stepping structure
!#
!# 2.) tstep_releaseTimestep
!#     -> Releases an existing time-stepping structure
!#
!# 3.) tstep_infoTimestep
!#     -> Outputs information about the time-stepping structure
!#
!# 4.) tstep_resetTimestep
!#     -> Resets the time-stepping structure to initial values
!#
!# 5.) tstep_removeTempFromTimestep
!#     -> Removes temporal storage from time-stepping structure
!#
!# 6.) tstep_performThetaStep = tstep_performThetaStepSc /
!#                              tstep_performThetaStepScCpl /
!#                              tstep_performThetaStepBl /
!#                              tstep_performThetaStepBlCpl
!#     -> Performs one step by the two-level theta scheme to compute the
!#        solution for the time interval dTime..dTime+dStep.
!#
!# 7.) tstep_performRKStep = tstep_performRKStepSc /
!#                           tstep_performRKStepScCpl /
!#                           tstep_performRKStepBl /
!#                           tstep_performRKStepBlCpl /
!#     -> Performs one step of an explicit Runge-Kutta scheme to compute the
!#        solution for the time interval dTime..dTime+dStep
!#
!# 8.) tstep_performPseudoStepping = tstep_performPseudoSteppingSc /
!#                                   tstep_performPseudoSteppingBl
!#     -> Performs pseudo time-stepping to compute the steady state solution
!#
!#
!# The following auxiliary routines are available:
!#
!# 1.) tstep_checkTimestep = tstep_checkTimestep /
!#                           tstep_checlTimestepCpl
!#     -> Checks the solution computed in one time step and adjust
!#        the size of the time step accordingly
!#
!# 2.) tstep_decodeOutputLevel
!#     -> Decodes the output level information into bitfields
!#
!# </purpose>
!##############################################################################

module timestep

#include "flagship.h"

!$ use omp_lib
  use collection
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use problem
  use timestepaux
  use solveraux
  use solverlinear
  use solvernonlinear

  implicit none

  private
  public :: tstep_createTimestep
  public :: tstep_releaseTimestep
  public :: tstep_infoTimestep
  public :: tstep_resetTimestep
  public :: tstep_removeTempFromTimestep
  public :: tstep_performThetaStep
  public :: tstep_performRKStep
  public :: tstep_performPseudoStepping

  ! *****************************************************************************

  interface tstep_createTimestep
    module procedure tstep_createTimestepDirect
    module procedure tstep_createTimestepIndirect
  end interface

  interface tstep_performThetaStep
    module procedure tstep_performThetaStepSc
    module procedure tstep_performThetaStepBl
    module procedure tstep_performThetaStepScCpl
    module procedure tstep_performThetaStepBlCpl
  end interface

  interface tstep_performRKStep
    module procedure tstep_performRKStepSc
    module procedure tstep_performRKStepBl
    module procedure tstep_performRKStepScCpl
    module procedure tstep_performRKStepBlCpl
  end interface

  interface tstep_performPseudoStepping
    module procedure tstep_performPseudoSteppingSc
    module procedure tstep_performPseudoSteppingBl
  end interface

  interface tstep_checkTimestep
    module procedure tstep_checkTimestep
    module procedure tstep_checkTimestepCpl
  end interface

  ! *****************************************************************************

!<constants>

!<constantblock description="Global Runge-Kutta weights">

  ! Runge-Kutta one-stage
  real(DP), dimension(1), parameter :: TSTEP_RK1 = (/ 1.0_DP /)

  ! Runge-Kutta two-stage
  real(DP), dimension(2), parameter :: TSTEP_RK2 = (/ 0.5_DP, 1.0_DP /)

  ! Runge-Kutta three-stage
  real(DP), dimension(3), parameter :: TSTEP_RK3 = (/ 0.6_DP, 0.6_DP, 1.0_DP /)

  ! Runge-Kutta four-stage
  real(DP), dimension(4), parameter :: TSTEP_RK4 = (/ 0.25_DP, 0.5_DP, 0.5_DP, 1.0_DP /)
!</constantblock>

!</constants>

contains

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_createTimestepDirect(rparlist, ssectionName,&
      rtimestep, nproblemsCoupled)

!<description>
    ! This subroutine creates a new time-stepping structure from a
    ! given parameter list
!</description>

!<input>
    ! Parameter list containing all data
    type(t_parlist), intent(in) :: rparlist

    ! Section name of the parameter list containing solver data
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: number of coupled problems
    ! If not given nproblemsCoupled = 1 is assumed.
    integer, intent(in), optional :: nproblemsCoupled
!</input>

!<output>
    ! Time-stepping object
    type(t_timestep), intent(out) :: rtimestep
!</output>
!</subroutine>

    ! local variables
    integer :: npc

    ! Set number of coupled problems
    npc = 1
    if (present(nproblemsCoupled)) npc = nproblemsCoupled


    ! The INTENT(out) already initialises rtimestep with the most
    ! important information.
    rtimestep%sName = trim(adjustl(ssectionName))

    ! Get mandatory configuration values from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
        "ctimestepType", rtimestep%ctimestepType)
    call parlst_getvalue_double(rparlist, ssectionName,&
        "dfinalTime", rtimestep%dfinalTime)
    call parlst_getvalue_double(rparlist, ssectionName,&
        "dinitialStep", rtimestep%dinitialStep)

    ! Get optional configuration values from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
        "ioutputlevel", rtimestep%ioutputlevel)
    call parlst_getvalue_int(rparlist, ssectionName,&
        "isolNorm", rtimestep%isolNorm)
    call parlst_getvalue_double(rparlist, ssectionName,&
        "dinitialTime", rtimestep%dinitialTime)
    call parlst_getvalue_double(rparlist, ssectionName,&
        "dminStep", rtimestep%dminStep, rtimestep%dinitialStep)
    call parlst_getvalue_double(rparlist, ssectionName,&
        "dmaxStep", rtimestep%dmaxStep, rtimestep%dinitialStep)
     call parlst_getvalue_double(rparlist, ssectionName,&
         "dstepReductionFactor", rtimestep%dstepReductionFactor)
    call parlst_getvalue_double(rparlist, ssectionName,&
        "depsSteady", rtimestep%depsSteady)
    call parlst_getvalue_int(rparlist, ssectionName,&
        "iadapttimestep", rtimestep%iadapttimestep)
    call parlst_getvalue_double(rparlist, ssectionName,&
        "dadaptTime", rtimestep%dadaptTime)

    ! Decode the output level
    call tstep_decodeOutputLevel(rtimestep)

    ! Get solver dependent configuration values from parameter list
    select case(rtimestep%ctimestepType)
    case (TSTEP_THETA_SCHEME)
      ! Two-level theta scheme
      call parlst_getvalue_double(rparlist, ssectionName,&
          "theta", rtimestep%theta)

      ! Allocate array of temporal vectors
      allocate(rtimestep%RtempVectors(npc))

    case (TSTEP_RK_SCHEME)
      ! Multilevel Runge-Kutta scheme
      call parlst_getvalue_int(rparlist, ssectionName,&
          "multisteps", rtimestep%multisteps)

      select case(rtimestep%multisteps)
      case (1)
        allocate(rtimestep%DmultistepWeights(npc))
        rtimestep%DmultistepWeights = TSTEP_RK1

      case (2)
        allocate(rtimestep%DmultistepWeights(2*npc))
        rtimestep%DmultistepWeights = TSTEP_RK2

      case (3)
        allocate(rtimestep%DmultistepWeights(3*npc))
        rtimestep%DmultistepWeights = TSTEP_RK3

      case (4)
        allocate(rtimestep%DmultistepWeights(4*npc))
        rtimestep%DmultistepWeights = TSTEP_RK4

      case default
        call output_line('Number of Runge-Kutta steps must be specified explicitly!',&
            OU_CLASS_ERROR,rtimestep%coutputModeError,&
            'tstep_createTimestepDirect')
        call sys_halt()
      end select

      ! Allocate array of temporal vectors
      allocate(rtimestep%RtempVectors(4))

    case default
      if (rtimestep%coutputModeError .gt. 0) then
        call output_line('Invalid type of time stepping algorithm!',&
            OU_CLASS_ERROR,rtimestep%coutputModeError,&
            'tstep_createTimestepDirect')
      end if
      call sys_halt()
    end select


    ! Perform adaptive time-stepping?
    select case(rtimestep%iadapttimestep)

    case (TSTEP_NOADAPT)
      ! No adaptive time-stepping

    case (TSTEP_SERADAPT)
      ! Adaptive time-stepping using switched evolution relaxation
      allocate(rtimestep%p_rserController)
      call parlst_getvalue_double(rparlist, ssectionName,&
          "dIncreaseFactor", rtimestep%p_rserController%dIncreaseFactor)
      call parlst_getvalue_double(rparlist, ssectionName,&
          "dDecreaseFactor", rtimestep%p_rserController%dDecreaseFactor)


    case (TSTEP_AUTOADAPT)
      ! Adaptive time-stepping using automatic time step control
      allocate(rtimestep%p_rautoController)

      ! Allocate array of temporal vectors
      allocate(rtimestep%p_rautoController%RtempVectors(2*npc))

      call parlst_getvalue_double(rparlist, ssectionName,&
          "dDecreaseFactor", rtimestep%p_rautoController%dDecreaseFactor)
      call parlst_getvalue_double(rparlist, ssectionName,&
          "depsRel", rtimestep%p_rautoController%depsRel)


    case (TSTEP_PIDADAPT)
      ! Adaptive time-stepping using the PID controller
      allocate(rtimestep%p_rpidController)

      ! Get configuration values from parameter list
      call parlst_getvalue_double(rparlist, ssectionName,&
          "dProportionalExponent", rtimestep%p_rpidController%dProportionalExponent)
      call parlst_getvalue_double(rparlist, ssectionName,&
          "dIntegralExponent", rtimestep%p_rpidController%dIntegralExponent)
      call parlst_getvalue_double(rparlist, ssectionName,&
          "dDerivativeExponent", rtimestep%p_rpidController%dDerivativeExponent)
      call parlst_getvalue_double(rparlist, ssectionName,&
          "dIncreaseFactor", rtimestep%p_rpidController%dIncreaseFactor)
      call parlst_getvalue_double(rparlist, ssectionName,&
          "dDecreaseFactor", rtimestep%p_rpidController%dDecreaseFactor)
      call parlst_getvalue_double(rparlist, ssectionName,&
          "depsRel", rtimestep%p_rpidController%depsRel)
      call parlst_getvalue_double(rparlist, ssectionName,&
          "dmaxRel", rtimestep%p_rpidController%dmaxRel)

    case default
      if (rtimestep%coutputModeError .gt. 0) then
        call output_line('Invalid type of adaptive time-stepping algorithm!',&
            OU_CLASS_ERROR,rtimestep%coutputModeError,&
            'tstep_createTimestepDirect')
      end if
      call sys_halt()
    end select


    ! Reset any other values
    call tstep_resetTimestep(rtimestep, .true.)

  end subroutine tstep_createTimestepDirect

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_createTimestepIndirect(rtimestep, rtimestepTemplate)

!<description>
    ! This subroutine creates a new time-stepping structure by cloning
    ! an existing time-stepping structure
!</description>

!<input>
    ! Template timestep structure
    type(t_timestep), intent(in) :: rtimestepTemplate
!</input>

!<output>
    ! Timestep structure
    type(t_timestep), intent(out) :: rtimestep
!</output>
!</subroutine>


    ! The INTENT(out) already initialises rtimestep with the most
    ! important information. The rest comes now
    rtimestep = rtimestepTemplate

    ! Create PID controller
    if (associated(rtimestepTemplate%p_rpidController)) then
      allocate(rtimestep%p_rpidController)
      rtimestep%p_rpidController = rtimestepTemplate%p_rpidController
    end if

    ! Create automatic controller
    if (associated(rtimestepTemplate%p_rautoController)) then
      allocate(rtimestep%p_rautoController)
      rtimestep%p_rautoController = rtimestepTemplate%p_rautoController
      if (associated(rtimestepTemplate%p_rautoController%RtempVectors)) then
        allocate(rtimestep%p_rautoController%RtempVectors(&
            size(rtimestepTemplate%p_rautoController%RtempVectors)))
      end if
    end if

    ! Create SER controller
    if (associated(rtimestepTemplate%p_rserController)) then
      allocate(rtimestep%p_rserController)
      rtimestep%p_rserController = rtimestepTemplate%p_rserController
    end if

    ! Create temporal vectors
    if (associated(rtimestepTemplate%RtempVectors)) then
      allocate(rtimestep%RtempVectors(&
          size(rtimestepTemplate%RtempVectors)))
    end if

    ! Create multistep weights
    if (associated(rtimestepTemplate%DmultistepWeights)) then
      allocate(rtimestep%DmultistepWeights(&
          size(rtimestepTemplate%DmultistepWeights)))
      rtimestep%DmultistepWeights = rtimestepTemplate%DmultistepWeights
    end if

  end subroutine tstep_createTimestepIndirect

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_releaseTimestep(rtimestep)

!<description>
    ! This subroutine releases an existing time-stepping structure
!</description>

!<inputoutput>
    ! Time-stepping object
    type(t_timestep), intent(inout) :: rtimestep
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i


    ! Release PID controller
    if (associated(rtimestep%p_rpidController)) then
      deallocate(rtimestep%p_rpidController)
      nullify(rtimestep%p_rpidController)
    end if

    ! Release automatic controller
    if (associated(rtimestep%p_rautoController)) then
      if (associated(rtimestep%p_rautoController%RtempVectors)) then
        do i = lbound(rtimestep%p_rautoController%RtempVectors, 1),&
               ubound(rtimestep%p_rautoController%RtempVectors, 1)
          call lsysbl_releaseVector(rtimestep%p_rautoController%RtempVectors(i))
        end do
        deallocate(rtimestep%p_rautoController%RtempVectors)
        nullify(rtimestep%p_rautoController%RtempVectors)
      end if
      deallocate(rtimestep%p_rautoController)
      nullify(rtimestep%p_rautoController)
    end if

    ! Release SER controller
    if (associated(rtimestep%p_rserController)) then
      deallocate(rtimestep%p_rserController)
      nullify(rtimestep%p_rserController)
    end if

    ! Release temporal vector
    if (associated(rtimestep%RtempVectors)) then
      do i = lbound(rtimestep%RtempVectors, 1),&
             ubound(rtimestep%RtempVectors, 1)
        call lsysbl_releaseVector(rtimestep%RtempVectors(i))
      end do
      deallocate(rtimestep%RtempVectors)
      nullify(rtimestep%RtempVectors)
    end if

    ! Release multistep weights
    if (associated(rtimestep%DmultistepWeights)) then
      deallocate(rtimestep%DmultistepWeights)
      nullify(rtimestep%DmultistepWeights)
    end if

  end subroutine tstep_releaseTimestep

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_infoTimestep(rtimestep, bprintInternal)

!<description>
    ! This subroutine prints information about the time-stepping structure
!</description>

!<input>
    ! OPTIONAL: Print internal data?
    logical, intent(in), optional :: bprintInternal
!</input>

!<inputoutput>
    ! Time-stepping object
    type(t_timestep), intent(inout) :: rtimestep
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i


    ! Output general information
    call output_line('Timestep:')
    call output_line('---------')
    call output_line('Name:                          '//trim(rtimestep%sName))
    call output_line('Number of time steps:          '//trim(sys_siL(rtimestep%nSteps,15)))
    call output_line('Number of rejected time steps: '//trim(sys_siL(rtimestep%nrejectedSteps,15)))

    ! Output detailed information
    if (present(bprintInternal)) then
      if (bprintInternal) then

        call output_line('ctimestepType:                 '//trim(sys_siL(rtimestep%ctimestepType,3)))
        call output_line('ioutputLevel:                  '//trim(sys_siL(rtimestep%ioutputLevel,3)))
        call output_line('isolNorm:                      '//trim(sys_siL(rtimestep%isolNorm,3)))
        call output_line('iadaptTimestep:                '//trim(sys_siL(rtimestep%iadaptTimestep,3)))
        call output_line('multiSteps:                    '//trim(sys_siL(rtimestep%multiSteps,3)))
        call output_line('theta:                         '//trim(sys_sdL(rtimestep%theta,5)))
        call output_line('dinitialTime:                  '//trim(sys_sdL(rtimestep%dinitialTime,5)))
        call output_line('dfinalTime:                    '//trim(sys_sdL(rtimestep%dfinalTime,5)))
        call output_line('dminStep:                      '//trim(sys_sdL(rtimestep%dminStep,5)))
        call output_line('dmaxStep:                      '//trim(sys_sdL(rtimestep%dmaxStep,5)))
        call output_line('dinitialStep:                  '//trim(sys_sdL(rtimestep%dinitialStep,5)))
        call output_line('dstepReductionFactor:          '//trim(sys_sdL(rtimestep%dstepReductionFactor,5)))
        call output_line('drelChange:                    '//trim(sys_sdL(rtimestep%drelChange,5)))
        call output_line('depsSteady:                    '//trim(sys_sdL(rtimestep%depsSteady,5)))
        call output_line('dTime:                         '//trim(sys_sdL(rtimestep%dTime,5)))
        call output_line('dStep:                         '//trim(sys_sdL(rtimestep%dstep,5)))
        call output_line('dStep1:                        '//trim(sys_sdL(rtimestep%dstep1,5)))

        ! Output information about weights of multistep method
        if (associated(rtimestep%DmultistepWeights)) then
          do i = lbound(rtimestep%DmultistepWeights, 1),&
                 ubound(rtimestep%DmultistepWeights, 1)
            call output_line('multistep weight['//trim(sys_siL(i,1))//']:           '//&
                             trim(sys_sdL(rtimestep%DmultistepWeights(i),5)))
          end do
        end if

        ! Output information about the evolutionary PID controller
        if (associated(rtimestep%p_rpidController)) then
          call output_lbrk()
          call output_line('Evolutionary PID controller:')
          call output_line('----------------------------')
          call output_line('dProportionalExponent:         '//&
              trim(sys_sdL(rtimestep%p_rpidController%dProportionalExponent,5)))
          call output_line('dIntegralExponent:             '//&
              trim(sys_sdL(rtimestep%p_rpidController%dIntegralExponent,5)))
          call output_line('dDerivativeExponent:           '//&
              trim(sys_sdL(rtimestep%p_rpidController%dDerivativeExponent,5)))
          call output_line('dIncreaseFactor:               '//&
              trim(sys_sdL(rtimestep%p_rpidController%dIncreaseFactor,5)))
          call output_line('dDecreaseFactor:               '//&
              trim(sys_sdL(rtimestep%p_rpidController%dDecreaseFactor,5)))
          call output_line('depsRel:                       '//&
              trim(sys_sdL(rtimestep%p_rpidController%depsRel,5)))
          call output_line('dmaxRel:                       '//&
              trim(sys_sdL(rtimestep%p_rpidController%dmaxRel,5)))
          call output_line('dcontrolValue1:                '//&
              trim(sys_sdL(rtimestep%p_rpidController%dcontrolValue1,5)))
          call output_line('dcontrolValue2:                '//&
              trim(sys_sdL(rtimestep%p_rpidController%dcontrolValue2,5)))
        end if

        ! Output information about the automatic time step controller
        if (associated(rtimestep%p_rautoController)) then
          call output_lbrk()
          call output_line('Automatic time step controller:')
          call output_line('-------------------------------')
          call output_line('dDecreaseFactor:               '//&
              trim(sys_sdL(rtimestep%p_rautoController%dDecreaseFactor,5)))
          call output_line('depsRel:                       '//&
              trim(sys_sdL(rtimestep%p_rautoController%depsRel,5)))

          if (associated(rtimestep%p_rautoController%RtempVectors)) then
            call output_lbrk()
            call output_line('Temporal vectors:')
            call output_line('-----------------')
            do i = lbound(rtimestep%p_rautoController%RtempVectors, 1),&
                   ubound(rtimestep%p_rautoController%RtempVectors, 1)
              call lsysbl_infoVector(rtimestep%p_rautoController%RtempVectors(i))
            end do
          end if
        end if

        ! Output information about the switched evolution relaxation controller
        if (associated(rtimestep%p_rserController)) then
          call output_lbrk()
          call output_line('Switched evolution relaxation (SER) controller:')
          call output_line('-----------------------------------------------')
          call output_line('dIncreaseFactor:               '//&
              trim(sys_sdL(rtimestep%p_rserController%dIncreaseFactor,5)))
          call output_line('dDecreaseFactor:               '//&
              trim(sys_sdL(rtimestep%p_rserController%dDecreaseFactor,5)))
          call output_line('dsteadyDefect:                 '//&
              trim(sys_sdL(rtimestep%p_rserController%dsteadyDefect,5)))
          call output_line('dsteadyDefect1:                '//&
              trim(sys_sdL(rtimestep%p_rserController%dsteadyDefect1,5)))
        end if

        ! Output information about auxiliary vectors
        if (associated(rtimestep%RtempVectors)) then
          call output_lbrk()
          call output_line('Temporal vectors:')
          call output_line('-----------------')
          do i = lbound(rtimestep%RtempVectors, 1),&
                 ubound(rtimestep%RtempVectors, 1)
            call lsysbl_infoVector(rtimestep%RtempVectors(i))
          end do
        end if
      end if
    end if

    call output_lbrk()

  end subroutine tstep_infoTimestep

  ! ***************************************************************************

!<subroutine>

  subroutine tstep_resetTimestep(rtimestep, bresetStatistics)

!<description>
    ! This subroutine resets the time-stepping structure
!</description>

!<input>
    ! If true, the statistical output parameters are reset
    logical, intent(in) :: bresetStatistics
!</input>

!<inputoutput>
    ! time-stepping object
    type(t_timestep), intent(inout) :: rtimestep
!</inputoutput>
!</subroutine>

    ! Reset time step
    rtimestep%dTime  = rtimeStep%dinitialTime
    rtimestep%dStep  = rtimeStep%dinitialStep
    rtimestep%dStep1 = rtimeStep%dinitialStep

    ! Reset statistical data (if required)
    if (bresetStatistics) then
      rtimestep%drelChange     = 0.0_DP
      rtimestep%nSteps         = 0
      rtimestep%nrejectedSteps = 0
    end if

    ! Reset PID controller
    if (associated(rtimestep%p_rpidController)) then
      rtimestep%p_rpidController%dcontrolValue1 = 1.0_DP
      rtimestep%p_rpidController%dcontrolValue2 = 1.0_DP
    end if

    ! Reset SER controller
    if (associated(rtimestep%p_rserController)) then
      rtimestep%p_rserController%dsteadyDefect  = 0.0_DP
      rtimestep%p_rserController%dsteadyDefect1 = 0.0_DP
    end if

  end subroutine tstep_resetTimestep

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_removeTempFromTimestep(rtimestep)

!<description>
    ! This subroutines removes all temporal memory from the time-stepping structure
!</description>

!<inputoutput>
    ! timestep structure
    type(t_timestep), intent(inout) :: rtimestep
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i

    ! Release temporal vectors in automatic time step control
    if (associated(rtimestep%p_rautoController)) then
      if (associated(rtimestep%p_rautoController%RtempVectors)) then
        do i = lbound(rtimestep%p_rautoController%RtempVectors, 1),&
               ubound(rtimestep%p_rautoController%RtempVectors, 1)
          call lsysbl_releaseVector(rtimestep%p_rautoController%RtempVectors(i))
        end do
      end if
    end if

    ! Release temporal vectors in time step structure
    if (associated(rtimestep%RtempVectors)) then
      do i = lbound(rtimestep%RtempVectors, 1),&
             ubound(rtimestep%RtempVectors, 1)
        call lsysbl_releaseVector(rtimestep%RtempVectors(i))
      end do
    end if

  end subroutine tstep_removeTempFromTimestep

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performThetaStepSc(rproblemLevel, rtimestep,&
      rsolver, rsolution, fcb_nlsolverCallback, rcollection, rsource)

!<description>
    ! This subroutine performs one step by means of the two-level
    ! theta-scheme
    !
    !   $$ Mu^{n+1}-\theta\Delta t N(u^{n+1}) =&
    !      Mu^n+(1-\theta)\Delta t N(u^n) + b $$
    !
    ! to advance the solution from time level $t^n$ to time level
    ! $t^{n+1}$.  In the above equation $M$ denotes the mass matrix,
    ! and $N(\cdot)$ is an abstract operator which depends in the
    ! solution $u$. The constant right-hand side vector is optional an
    ! is assumed to be zero if not present.
    !
    ! This routine can handle the forward Euler (theta=0), backward
    ! Euler (theta=1) and the Crank-Nicolson (theta=0.5) time stepping
    ! algorithm. However, the fully explicit scheme is known to be
    ! unstable. It is therefore highly recommended to use an explicit
    ! Runge-Kutta scheme instead.
    !
    ! Note that this subroutine serves as wrapper for scalar solution
    ! vectors only.
!</description>

!<input>
    ! Callback routines
    include 'intf_solvercallback.inc'

    ! OPTIONAL: constant right-hand side vector
    type(t_vectorScalar), intent(in), optional :: rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! solution vector
    type(t_vectorScalar), intent(inout) :: rsolution

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rsolutionBlock, rsourceBlock


    if (present(rsource)) then

      ! Convert scalar vectors to 1-block vectors
      call lsysbl_createVecFromScalar(rsolution, rsolutionBlock)
      call lsysbl_createVecFromScalar(rsource, rsourceBlock)
      
      ! Call block version of this subroutine
      call tstep_performThetaStepBl(rproblemLevel, rtimestep,&
          rsolver, rsolutionBlock, fcb_nlsolverCallback, rcollection,&
          rsourceBlock)

      ! Deallocate temporal 1-block vectors
      call lsysbl_releaseVector(rsolutionBlock)
      call lsysbl_releaseVector(rsourceBlock)

    else

      ! Convert scalar vector to 1-block vector
      call lsysbl_createVecFromScalar(rsolution, rsolutionBlock)

      ! Call block version of this subroutine
      call tstep_performThetaStepBl(rproblemLevel, rtimestep,&
          rsolver, rsolutionBlock, fcb_nlsolverCallback, rcollection)

      ! Deallocate temporal 1-block vector
      call lsysbl_releaseVector(rsolutionBlock)

    end if

  end subroutine tstep_performThetaStepSc

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performThetaStepBl(rproblemLevel, rtimestep,&
      rsolver, rsolution, fcb_nlsolverCallback, rcollection, rsource)

!<description>
    ! This subroutine performs one step by means of the two-level
    ! theta-scheme
    !
    !   $$ Mu^{n+1}-\theta\Delta t N(u^{n+1}) =
    !      Mu^n+(1-\theta)\Delta t N(u^n) + b $$
    !
    ! to advance the solution from time level $t^n$ to time level
    ! $t^{n+1}$.  In the above equation $M$ denotes the mass matrix,
    ! and $N(\cdot)$ is an abstract operator which depends in the
    ! solution $u$. The constant right-hand side vector is optional an
    ! is assumed to be zero if not present.
    !
    ! This routine can handle the forward Euler (theta=0), backward
    ! Euler (theta=1) and the Crank-Nicolson (theta=0.5) time stepping
    ! algorithm. However, the fully explicit scheme is known to be
    ! unstable. It is therefore highly recommended to use an explicit
    ! Runge-Kutta scheme instead.
!</description>

!<input>
    ! Callback routines
    include 'intf_solvercallback.inc'

    ! OPTIONAL: source vector
    type(t_vectorBlock), intent(in), optional :: rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver structure
    type(t_solver), intent(inout), target :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_solver), pointer :: p_rsolver
    type(t_vectorBlock), pointer :: p_rsolutionRef
    type(t_vectorBlock), pointer :: p_rsolutionAux
    type(t_vectorBlock), pointer :: p_rsolutionOld
    logical :: bcompatible, breject

    ! Set pointer to nonlinear solver
    p_rsolver => solver_getNextSolverByTypes(rsolver,&
        (/SV_NONLINEARMG, SV_NONLINEAR/))

    if (.not. associated(p_rsolver)) then
      if (rtimestep%coutputModeError .gt. 0) then
        call output_line('Unsupported/invalid solver type!',&
            OU_CLASS_ERROR,rtimestep%coutputModeError,&
            'tstep_performThetaStepBl')
      end if
      call sys_halt()
    end if

    ! Set pointers to temporal vectors
    p_rsolutionOld => rtimestep%RtempVectors(1)

    ! ... and check if vectors are compatible
    call lsysbl_isVectorCompatible(rsolution,&
        p_rsolutionOld, bcompatible)
    if (.not.bcompatible) then
      call lsysbl_resizeVectorBlock(p_rsolutionOld,&
          rsolution, .false.)
    end if

    ! Save the given solution vector to the temporal vector. If the
    ! computed time step is not accepted, then the backup of the
    ! given solution vector is used to recalculate the time step
    call lsysbl_copyVector(rsolution, p_rsolutionOld)


    ! Set pointers to temporal vectors for the computation of an
    ! auxiliary solution by means of local substepping
    if (rtimestep%iadaptTimestep .eq. TSTEP_AUTOADAPT) then
      p_rsolutionRef => rtimestep%p_rautoController%RtempVectors(1)
      p_rsolutionAux => rtimestep%p_rautoController%RtempVectors(2)

      ! ... and check if vectors are compatible
      call lsysbl_isVectorCompatible(rsolution, p_rsolutionRef, bcompatible)
      if (.not.bcompatible) then
        call lsysbl_resizeVectorBlock(p_rsolutionRef, rsolution, .false.)
        call lsysbl_resizeVectorBlock(p_rsolutionAux, rsolution, .false.)
      end if

      ! Make a backup copy of the solution vector
      call lsysbl_copyVector(rsolution, p_rsolutionRef)
    end if


    ! Adaptive time-stepping loop
    timeadapt: do

      ! Increase simulation time provisionally
      rtimestep%dStep = min(rtimestep%dStep, rtimestep%dfinalTime - rtimestep%dTime)
      rtimestep%dTime = rtimestep%dTime  + rtimestep%dStep
      rtimestep%nSteps= rtimestep%nSteps + 1

      ! Output information
      if (rtimestep%coutputModeInfo .gt. 0) then
        call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('Two-level theta-scheme, Time = '//&
                         trim(sys_sdEL(rtimestep%dTime,5))//&
                         ' Stepsize = '//trim(sys_sdEL(rtimestep%dStep,5)),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
      end if

      ! Solve the nonlinear algebraic system in the time interval (t^n, t^{n+1})
      call nlsol_solveMultigrid(rproblemLevel, rtimestep, p_rsolver,&
          rsolution, p_rsolutionOld, fcb_nlsolverCallback,&
          rcollection, rsource)

      ! Adjust status information of top-most solver
      call solver_copySolver(p_rsolver, rsolver, .false., .true.)


      ! Compute reference solution by performing two time steps of size Dt/2
      if (rtimestep%iadaptTimestep .eq. TSTEP_AUTOADAPT) then

        ! Set time step to smaller value
        rtimestep%dStep = rtimestep%dStep/2.0_DP

        ! Output information
        if (rtimestep%coutputModeInfo .gt. 0) then
          call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_line('First substep in automatic time step control',&
                           OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
        end if

        ! Solve the nonlinear algebraic system for time step t^n -> t^{n+1/2}
        call nlsol_solveMultigrid(rproblemLevel, rtimestep, p_rsolver,&
            p_rsolutionRef, p_rsolutionOld, fcb_nlsolverCallback,&
            rcollection, rsource)

        ! Save intermediate solution
        call lsysbl_copyVector(p_rsolutionRef, p_rsolutionAux)

        ! Output information
        if (rtimestep%coutputModeInfo .gt. 0) then
          call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_line('Second substep in automatic time step control',&
                           OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
        end if

        ! Solve the nonlinear algebraic system for time step t^{n+1/2} -> t^{n+1}
        call nlsol_solveMultigrid(rproblemLevel, rtimestep, p_rsolver,&
            p_rsolutionRef, p_rsolutionAux, fcb_nlsolverCallback,&
            rcollection, rsource)

        ! Check if solution from this time step can be accepted and
        ! adjust the time step size automatically if this is required
        breject = tstep_checkTimestep(rtimestep, p_rsolver,&
            p_rsolutionRef, p_rsolutionOld)

        ! Prepare time step size for next "large" time step
        rtimestep%dStep = rtimestep%dStep*2.0_DP

      else

        ! Check if solution from this time step can be accepted and
        ! adjust the time step size automatically if this is required
        breject = tstep_checkTimestep(rtimestep, p_rsolver,&
            rsolution, p_rsolutionOld)

      end if

      ! Output information
      if (rtimestep%coutputModeInfo .gt. 0) then
        call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_separator(OU_SEP_TILDE,OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('Time step was     '//merge('!!! rejected !!!','accepted        ',breject),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('New stepsize:     '//trim(sys_sdEL(rtimestep%dStep,5)),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('Last stepsize:    '//trim(sys_sdEL(rtimestep%dStep1,5)),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('Relative changes: '//trim(sys_sdEL(rtimestep%drelChange,5)),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_separator(OU_SEP_TILDE,OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
      end if

      ! Do we have to reject to current solution?
      if (breject) then
        ! Yes, so restore the old solution and
        ! repeat the adaptive time-stepping loop
        call lsysbl_copyVector(p_rsolutionOld, rsolution)
        
        if (rtimestep%coutputModeWarning .gt. 0) then
          call output_line('Time step was rejected!',&
              OU_CLASS_WARNING,rtimestep%coutputModeWarning,&
              'tstep_performThetaStepBl')
        end if
        
      else
        ! No, accept current solution and
        ! exit adaptive time-stepping loop
        exit timeadapt
      end if
    end do timeadapt

  end subroutine tstep_performThetaStepBl

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performThetaStepScCpl(rproblemLevel, rtimestep,&
      rsolver, Rsolution, fcb_nlsolverCallback, rcollection, Rsource)

!<description>
    ! This subroutine performs one step by means of the two-level
    ! theta-scheme
    !
    !   $$ Mu^{n+1}-\theta\Delta t N(u^{n+1}) =
    !      Mu^n+(1-\theta)\Delta t N(u^n)$$
    !
    ! to advance the solution from time level $t^n$ to time level
    ! $t^{n+1}$.  In the above equation $M$ denotes the mass matrix,
    ! and $N(\cdot)$ is an abstract operator which depends in the
    ! solution $u$. The constant right-hand side vector is optional an
    ! is assumed to be zero if not present.
    !
    ! This routine can handle the forward Euler (theta=0), backward
    ! Euler (theta=1) and the Crank-Nicolson (theta=0.5) time stepping
    ! algorithm. However, the fully explicit scheme is known to be
    ! unstable. It is therefore highly recommended to use an explicit
    ! Runge-Kutta scheme instead.
    !
    ! In contrast to routine tstep_performThetaStepSc this routine
    ! is applicable to a sequence of coupled problems which are
    ! solved repeatedly until convergence.
    !
    ! Note that this subroutine serves as wrapper for scalar solution
    ! vectors only.
!</description>

!<input>
    ! Callback routines
    include 'intf_solvercallback.inc'

    ! OPTIONAL: constant right-hand side vector
    type(t_vectorScalar), dimension(:), intent(in), optional :: Rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! solution vector
    type(t_vectorScalar), dimension(:), intent(inout) :: Rsolution



    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), dimension(:), allocatable :: RsolutionBlock
    type(t_vectorBlock), dimension(:), allocatable :: RsourceBlock
    integer :: icomponent


    if (present(Rsource)) then

      ! Allocate temporal memory
      allocate(RsolutionBlock(size(Rsolution)))
      allocate(RsourceBlock(size(Rsource)))

      ! Convert scalar vectors to 1-block vectors
      do icomponent = 1, size(Rsolution)
        call lsysbl_createVecFromScalar(Rsolution(icomponent),&
            RsolutionBlock(icomponent))
      end do

      ! Convert scalar vectors to 1-block vectors
      do icomponent = 1, size(Rsource)
            call lsysbl_createVecFromScalar(Rsource(icomponent),&
                RsourceBlock(icomponent))
      end do

      ! Call block version of this subroutine
      call tstep_performThetaStepBlCpl(rproblemLevel, rtimestep,&
          rsolver, RsolutionBlock, fcb_nlsolverCallback, rcollection,&
          RsourceBlock)

      ! Deallocate temporal 1-block vectors
      do icomponent = 1, size(Rsolution)
        call lsysbl_releaseVector(RsolutionBlock(icomponent))
      end do

      ! Deallocate temporal 1-block vectors
      do icomponent = 1, size(Rsource)
        call lsysbl_releaseVector(RsourceBlock(icomponent))
      end do

      ! Release temporal memory
      deallocate(RsolutionBlock)
      deallocate(RsourceBlock)

    else

      ! Allocate temporal memory
      allocate(RsolutionBlock(size(Rsolution)))

      ! Convert scalar vectors to 1-block vectors
      do icomponent = 1, size(Rsolution)
        call lsysbl_createVecFromScalar(Rsolution(icomponent),&
            RsolutionBlock(icomponent))
      end do

      ! Call block version of this subroutine
      call tstep_performThetaStepBlCpl(rproblemLevel, rtimestep,&
          rsolver, RsolutionBlock, fcb_nlsolverCallback, rcollection)

      ! Deallocate temporal 1-block vectors
      do icomponent = 1, size(Rsolution)
        call lsysbl_releaseVector(RsolutionBlock(icomponent))
      end do

      ! Release temporal memory
      deallocate(rsolutionBlock)

    end if

  end subroutine tstep_performThetaStepScCpl

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performThetaStepBlCpl(rproblemLevel, rtimestep,&
      rsolver, Rsolution, fcb_nlsolverCallback, rcollection, Rsource)

!<description>
    ! This subroutine performs one step by means of the two-level
    ! theta-scheme
    !
    !   $$ Mu^{n+1}-\theta\Delta t N(u^{n+1}) =
    !      Mu^n+(1-\theta)\Delta t N(u^n)$$
    !
    ! to advance the solution from time level $t^n$ to time level
    ! $t^{n+1}$.  In the above equation $M$ denotes the mass matrix,
    ! and $N(\cdot)$ is an abstract operator which depends in the
    ! solution $u$. The constant right-hand side vector is optional an
    ! is assumed to be zero if not present.
    !
    ! This routine can handle the forward Euler (theta=0), backward
    ! Euler (theta=1) and the Crank-Nicolson (theta=0.5) time stepping
    ! algorithm. However, the fully explicit scheme is known to be
    ! unstable. It is therefore highly recommended to use an explicit
    ! Runge-Kutta scheme instead.
    !
    ! In contrast to routine tstep_performThetaStepBl this routine
    ! is applicable to a sequence of coupled problems which are
    ! solved repeatedly until convergence.
!</description>

!<input>
    ! Callback routines
    include 'intf_solvercallback.inc'

    ! OPTIONAL: constant right-hand side vector
    type(t_vectorBlock), dimension(:), intent(in), optional :: Rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver structure
    type(t_solver), intent(inout), target :: rsolver

    ! solution vector
    type(t_vectorBlock), dimension(:), intent(inout) :: Rsolution

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_solver), pointer :: p_rsolver
    type(t_vectorBlock), dimension(:), pointer :: p_RsolutionRef
    type(t_vectorBlock), dimension(:), pointer :: p_RsolutionAux
    type(t_vectorBlock), dimension(:), pointer :: p_RsolutionOld
    logical :: bcompatible, breject
    integer :: icomponent,ncomponent

    ! Get number of components
    ncomponent = size(Rsolution)

    ! Check size of multi-component source vector (if any)
    if (present(Rsource)) then
      if (size(Rsource) .ne. ncomponent) then
        if (rtimestep%coutputModeError .gt. 0) then
          call output_line('Dimension of coupled problem mismatch!',&
              OU_CLASS_ERROR,rtimestep%coutputModeError,&
              'tstep_performThetaStepBlCpl')
        end if
        call sys_halt()
      end if
    end if

    ! Set pointer to coupled solver
    p_rsolver => solver_getNextSolverByType(rsolver, SV_COUPLED)

    if (.not. associated(p_rsolver)) then
      if (rtimestep%coutputModeError .gt. 0) then
        call output_line('Unsupported/invalid solver type!',&
            OU_CLASS_ERROR,rtimestep%coutputModeError,&
            'tstep_performThetaStepBlCpl')
      end if
      call sys_halt()
    end if

    ! Set pointer to temporal vectors
    p_RsolutionOld => rtimestep%RtempVectors(1:ncomponent)

    ! ... and check if vectors are compatible
    do icomponent = 1, ncomponent
      call lsysbl_isVectorCompatible(Rsolution(icomponent),&
          p_RsolutionOld(icomponent), bcompatible)
      if (.not.bcompatible) then
        call lsysbl_resizeVectorBlock(p_RsolutionOld(icomponent),&
            Rsolution(icomponent), .false.)
      end if

      ! Save the given solution vector to the temporal vector. If the
      ! computed time step is not accepted, then the backup of the
      ! given solution vector is used to recalculate the time step
      call lsysbl_copyVector(Rsolution(icomponent),&
          p_RsolutionOld(icomponent))
    end do

    ! Set pointers to temporal vectors for the computation of an
    ! auxiliary solution by means of local substepping
    if (rtimestep%iadaptTimestep .eq. TSTEP_AUTOADAPT) then
      p_RsolutionRef =>&
          rtimestep%p_rautoController%RtempVectors(1:ncomponent)
      p_RsolutionAux =>&
          rtimestep%p_rautoController%RtempVectors(ncomponent+1:2*ncomponent)

      ! ... and check if vectors are compatible
      do icomponent = 1, ncomponent
        call lsysbl_isVectorCompatible(Rsolution(icomponent),&
            p_RsolutionRef(icomponent), bcompatible)
        if (.not.bcompatible) then
          call lsysbl_resizeVectorBlock(p_RsolutionRef(icomponent),&
              Rsolution(icomponent), .false.)
          call lsysbl_resizeVectorBlock(p_RsolutionAux(icomponent),&
              Rsolution(icomponent), .false.)
        end if

        ! Make a backup copy of the solution vector
        call lsysbl_copyVector(Rsolution(icomponent),&
            p_RsolutionRef(icomponent))
      end do
    end if


    ! Adaptive time-stepping loop
    timeadapt: do

      ! Increase simulation time provisionally
      rtimestep%dStep = min(rtimestep%dStep, rtimestep%dfinalTime - rtimestep%dTime)
      rtimestep%dTime = rtimestep%dTime  + rtimestep%dStep
      rtimestep%nSteps= rtimestep%nSteps + 1


      ! Output information
      if (rtimestep%coutputModeInfo .gt. 0) then
        call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('Two-level theta-scheme, Time = '//&
                         trim(sys_sdEL(rtimestep%dTime,5))//&
                         ' Stepsize = '//trim(sys_sdEL(rtimestep%dStep,5)),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
      end if

      ! Solve the coupled system in the time interval (t^n, t^{n+1})
      call nlsol_solveCoupled(rproblemLevel, rtimestep, p_rsolver,&
          Rsolution, p_RsolutionOld, fcb_nlsolverCallback,&
          rcollection, Rsource)

      ! Adjust status information of top-most solver
      call solver_copySolver(p_rsolver, rsolver, .false., .true.)

      ! Compute reference solution by performing two time steps of size Dt/2
      if (rtimestep%iadaptTimestep .eq. TSTEP_AUTOADAPT) then

        ! Set time step to smaller value
        rtimestep%dStep = rtimestep%dStep/2.0_DP

      if (rtimestep%coutputModeInfo .gt. 0) then
          call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_line('First substep in automatic time step control',&
                           OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
        end if


        ! Solve the nonlinear algebraic system for time step t^n -> t^{n+1/2}
        call nlsol_solveCoupled(rproblemLevel, rtimestep, p_Rsolver,&
            p_RsolutionRef, p_RsolutionOld, fcb_nlsolverCallback,&
            rcollection, Rsource)

        ! Save intermediate solution
        do icomponent = 1, ncomponent
          call lsysbl_copyVector(p_RsolutionRef(icomponent),&
              p_RsolutionAux(icomponent))
        end do

        if (rtimestep%coutputModeInfo .gt. 0) then
          call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_line('Second substep in automatic time step control',&
                           OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
        end if

        ! Solve the nonlinear algebraic system for time step t^{n+1/2} -> t^{n+1}
        call nlsol_solveCoupled(rproblemLevel, rtimestep, p_rsolver,&
            p_RsolutionRef, p_RsolutionAux, fcb_nlsolverCallback,&
            rcollection, Rsource)

        ! Check if solution from this time step can be accepted and
        ! adjust the time step size automatically if this is required
        breject = tstep_checkTimestep(rtimestep, p_rsolver,&
            p_RsolutionRef, p_RsolutionOld)

        ! Prepare time step size for next "large" time step
        rtimestep%dStep = rtimestep%dStep*2.0_DP

      else

        ! Check if solution from this time step can be accepted and
        ! adjust the time step size automatically if this is required
        breject = tstep_checkTimestep(rtimestep, p_rsolver,&
            Rsolution, p_RsolutionOld)

      end if

      if (rtimestep%coutputModeInfo .gt. 0) then
        call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_separator(OU_SEP_TILDE,OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('Time step was     '//merge('!!! rejected !!!','accepted        ',breject),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('New stepsize:     '//trim(sys_sdEL(rtimestep%dStep,5)),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('Last stepsize:    '//trim(sys_sdEL(rtimestep%dStep1,5)),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('Relative changes: '//trim(sys_sdEL(rtimestep%drelChange,5)),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_separator(OU_SEP_TILDE,OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
      end if

      ! Do we have to reject to current solution?
      if (breject) then
        ! Yes, so restore the old solution and
        ! repeat the adaptive time-stepping loop
        do icomponent = 1, ncomponent
          call lsysbl_copyVector(p_RsolutionOld(icomponent),&
              Rsolution(icomponent))
        end do
      else
        ! No, accept current solution and
        ! exit adaptive time-stepping loop
        exit timeadapt
      end if
    end do timeadapt

  end subroutine tstep_performThetaStepBlCpl

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performRKStepSc(rproblemLevel, rtimestep,&
      rsolver, rsolution, fcb_nlsolverCallback, rcollection, rsource)

!<description>
    ! This subroutine performs one step by means of an explicit
    ! Runge-Kutta scheme to advance the solution from time level $t^n$
    ! to time level $t^{n+1}$. Note that this subroutine serves as
    ! wrapper for scalar solution vectors only.
!</description>

!<input>
    ! Callback routines
    include 'intf_solvercallback.inc'

    ! OPTIONAL: constant right-hand side vector
    type(t_vectorScalar), intent(in), optional :: rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! solution vector
    type(t_vectorScalar), intent(inout) :: rsolution

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rsourceBlock
    type(t_vectorBlock) :: rsolutionBlock


    if (present(rsource)) then

      ! Convert scalar vectors to 1-block vectos
      call lsysbl_createVecFromScalar(rsource, rsourceBlock)
      call lsysbl_createVecFromScalar(rsolution, rsolutionBlock)
      
      ! Call block version of this routine
      call tstep_performRKStepBl(rproblemLevel, rtimestep, rsolver,&
          rsolutionBlock, fcb_nlsolverCallback, rcollection, rsourceBlock)
      
      ! Free temporal 1-block vectors
      call lsysbl_releaseVector(rsourceBlock)
      call lsysbl_releaseVector(rsolutionBlock)
      
    else

      ! Convert scalar vector to 1-block vecto
      call lsysbl_createVecFromScalar(rsolution, rsolutionBlock)

      ! Call block version of this routine
      call tstep_performRKStepBl(rproblemLevel, rtimestep, rsolver,&
          rsolutionBlock, fcb_nlsolverCallback, rcollection)

      ! Free temporal 1-block vector
      call lsysbl_releaseVector(rsolutionBlock)

    end if

  end subroutine tstep_performRKStepSc

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performRKStepBl(rproblemLevel, rtimestep,&
      rsolver, rsolution,  fcb_nlsolverCallback, rcollection, rsource)

!<description>
    ! This subroutine performs one step by means of an explicit
    ! Runge-Kutta scheme to advance the solution from time level $t^n$
    ! to time level $t^{n+1}$.  In particular, the two step version
    ! suggested by Richtmyer and Morton is implemented so that no
    ! artificial diffusion needs to be constructed.
!</description>

!<input>
    ! Callback routines
    include 'intf_solvercallback.inc'

    ! OPTIONAL: constant right-hand side vector
    type(t_vectorBlock), intent(in), optional :: rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(inout), target :: rtimestep

    ! solver structure
    type(t_solver), intent(inout), target :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_solver), pointer :: p_rsolver
    type(t_vectorBlock), pointer :: p_rsolutionInitial, p_rsolutionRef, p_rsolutionAux
    type(t_vectorBlock), pointer :: p_rconstB, p_raux
    integer(I32) :: ioperationSpec
    integer :: istep
    logical :: bcompatible,breject


    ! Set pointer to linear solver
    p_rsolver => solver_getNextSolverByTypes(rsolver, (/SV_LINEARMG, SV_LINEAR/))

    if (.not. associated(p_rsolver)) then
      if (rtimestep%coutputModeError .gt. 0) then
        call output_line('Unsupported/invalid solver type!',&
            OU_CLASS_ERROR,rtimestep%coutputModeError,&
            'tstep_performRKStepBl')
      end if
      call sys_halt()
    end if

    ! Set pointers to temporal vector
    p_rsolutionInitial => rtimestep%RtempVectors(1)
    p_rconstB          => rtimestep%RtempVectors(2)
    p_raux             => rtimestep%RtempVectors(3)

    ! Check if vectors are compatible
    call lsysbl_isVectorCompatible(rsolution, p_rsolutionInitial, bcompatible)
    if (.not.bcompatible) then
      call lsysbl_resizeVectorBlock(p_rsolutionInitial, rsolution, .false.)
      call lsysbl_resizeVectorBlock(p_rconstB, rsolution, .false.)
      call lsysbl_resizeVectorBlock(p_raux, rsolution, .false.)
    end if

    ! Save the given solution vector
    call lsysbl_copyVector(rsolution, p_rsolutionInitial)


    ! Set pointer to second temporal vector for the solution computed by local substepping
    if (rtimestep%iadaptTimestep .eq. TSTEP_AUTOADAPT) then
      p_rsolutionRef => rtimestep%p_rautoController%RtempVectors(1)
      p_rsolutionAux => rtimestep%p_rautoController%RtempVectors(2)

      ! And check if vectors are compatible
      call lsysbl_isVectorCompatible(rsolution, p_rsolutionRef, bcompatible)
      if (.not.bcompatible) then
        call lsysbl_resizeVectorBlock(p_rsolutionRef, rsolution, .false.)
        call lsysbl_resizeVectorBlock(p_rsolutionAux, rsolution, .false.)
      end if
      call lsysbl_copyVector(rsolution, p_rsolutionRef)
    end if


    ! Adaptive time-stepping loop
    timeadapt: do

      ! Increase simulation time provisionally
      rtimestep%dStep = min(rtimestep%dStep, rtimestep%dfinalTime-rtimestep%dTime)
      rtimestep%dTime = rtimestep%dTime  + rtimestep%dStep
      rtimestep%nSteps= rtimestep%nSteps + 1

      if (rtimestep%coutputModeInfo .gt. 0) then
        call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('Explicit Runge-Kutta scheme, Time = '//trim(sys_sdEL(rtimestep%dTime,5))//&
                         ' Stepsize = '//trim(sys_sdEL(rtimestep%dStep,5)),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
      end if


      ! Perform multi-step Runge-Kutta method
      do istep = 1, rtimestep%multisteps

        if (rtimestep%coutputModeInfo .gt. 0) then
          call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_line('Explicit Runge-Kutta step '//trim(sys_siL(istep,5)),&
                           OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
        end if

        ! Initialise the constant right-hand side vector
        if (present(rsource)) then
          if (rsource%NEQ .gt. 0) then
            call lsysbl_copyVector(rsource, p_rconstB)
          else
            call lsysbl_clearVector(p_rconstB)
          end if
        else
          call lsysbl_clearVector(p_rconstB)
        end if

        ! Compute the explicit right-hand side vector and the residual
        ioperationSpec = NLSOL_OPSPEC_CALCRHS +&
                         NLSOL_OPSPEC_CALCRESIDUAL


        print *, "Not working"
        stop

!!$        call fcb_nlsolverCallback(rproblemLevel, rtimestep, p_rsolver,&
!!$                                  rsolution, p_rsolutionInitial, p_rconstB,&
!!$                                  0, ioperationSpec, rcollection, istatus)

!!$         ! Compute the new right-hand side
!!$         call fcb_calcB(rproblemLevel, rtimestep, p_rsolver, rsolution,&
!!$                          p_rsolutionInitial, p_rb, istep, rcollection)
!!$
!!$         ! Apply given right-hand side vector
!!$         if (present(rsource)) call lsysbl_vectorLinearComb(rsource, p_rb, 1._DP, 1._DP)
!!$
!!$         ! Impose boundary conditions
!!$         call fcb_setBoundary(rproblemLevel, rtimestep, rsolver,&
!!$                              rsolution, p_rb, p_rsolutionInitial, rcollection)

         ! Call linear multigrid to solve A*u=b
         call lsysbl_clearVector(p_raux)
         call linsol_solveMultigrid(rproblemLevel, p_rsolver, p_raux, p_rconstB)

         ! Add increment to solution vector
         call lsysbl_vectorLinearComb(p_raux, p_rsolutionInitial, 1._DP, 1._DP, rsolution)
      end do


      ! Compute reference solution by performing two time steps of size Dt/2
      if (rtimestep%iadaptTimestep .eq. TSTEP_AUTOADAPT) then

        ! Set time step to smaller value
        rtimestep%dStep = rtimestep%dStep/2.0_DP

        if (rtimestep%coutputModeInfo .gt. 0) then
          call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_line('First substep in automatic time step control',&
                           OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
        end if


        ! Perform multi-step Runge-Kutta method for first step
        do istep = 1, rtimestep%multisteps

          if (rtimestep%coutputModeInfo .gt. 0) then
            call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
            call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
            call output_line('Explicit Runge-Kutta step '//trim(sys_siL(istep,5)),&
                             OU_CLASS_MSG,rtimestep%coutputModeInfo)
            call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
            call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
          end if

!!$          ! Compute the new right-hand side
!!$          call fcb_calcB(rproblemLevel, rtimestep, p_rsolver, p_rsolutionRef,&
!!$                           p_rsolutionInitial, p_rb, istep, rcollection)
!!$
!!$          ! Apply given right-hand side vector
!!$          if (present(rsource)) call lsysbl_vectorLinearComb(rsource, p_rb, 1._DP, 1._DP)
!!$
!!$          ! Impose boundary conditions
!!$          call fcb_setBoundary(rproblemLevel, rtimestep, rsolver,&
!!$                               p_rsolutionRef, p_rb, p_rsolutionInitial, rcollection)

          ! Call linear multigrid to solve A*u=b
          call lsysbl_clearVector(p_raux)
          call linsol_solveMultigrid(rproblemLevel, p_rsolver, p_raux, p_rconstB)

          ! Add increment to solution vector
          call lsysbl_vectorLinearComb(p_raux, p_rsolutionInitial, 1._DP, 1._DP, p_rsolutionRef)
        end do


        ! Save intermediate solution
        call lsysbl_copyVector(p_rsolutionRef, p_rsolutionAux)

        if (rtimestep%coutputModeInfo .gt. 0) then
          call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_line('Second substep in automatic time step control',&
                           OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
          call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
        end if


        ! Perform multi-step Runge-Kutta method for first step
        do istep = 1, rtimestep%multisteps

          if (rtimestep%coutputModeInfo .gt. 0) then
            call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
            call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
            call output_line('Explicit Runge-Kutta step '//trim(sys_siL(istep,5)),&
                             OU_CLASS_MSG,rtimestep%coutputModeInfo)
            call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
            call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
          end if


!!$          ! Compute the new right-hand side
!!$          call fcb_calcB(rproblemLevel, rtimestep, p_rsolver, p_rsolutionRef,&
!!$                           p_rsolutionAux, p_rb, istep, rcollection)
!!$
!!$          ! Apply given right-hand side vector
!!$          if (present(rsource)) call lsysbl_vectorLinearComb(rsource, p_rb, 1._DP, 1._DP)
!!$
!!$          ! Impose boundary conditions
!!$          call fcb_setBoundary(rproblemLevel, rtimestep, rsolver,&
!!$                               p_rsolutionRef, p_rb, p_rsolutionAux, rcollection)

          ! Call linear multigrid to solve A*u=b
          call lsysbl_clearVector(p_raux)
          call linsol_solveMultigrid(rproblemLevel, p_rsolver, p_raux, p_rconstB)

          ! Add increment to solution vector
          call lsysbl_vectorLinearComb(p_raux, p_rsolutionAux, 1._DP, 1._DP, p_rsolutionRef)
        end do

        ! Check if solution from this time step can be accepted and
        ! adjust the time step size automatically if this is required
        breject = tstep_checkTimestep(rtimestep, p_rsolver,&
                                      p_rsolutionRef, p_rsolutionInitial)

        ! Prepare time step size for next "large" time step
        rtimestep%dStep = rtimestep%dStep*2.0_DP

      else

        ! Check if solution from this time step can be accepted and
        ! adjust the time step size automatically if this is required
        breject = tstep_checkTimestep(rtimestep, p_rsolver,&
                                      rsolution, p_rsolutionInitial)

      end if

      if (rtimestep%coutputModeInfo .gt. 0) then
        call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_separator(OU_SEP_TILDE,OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('Time step was     '//merge('!!! rejected !!!','accepted        ',breject),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('New stepsize:     '//trim(sys_sdEL(rtimestep%dStep,5)),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('Last stepsize:    '//trim(sys_sdEL(rtimestep%dStep1,5)),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('Relative changes: '//trim(sys_sdEL(rtimestep%drelChange,5)),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_separator(OU_SEP_TILDE,OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
      end if

      ! Do we have to reject to current solution?
      if (breject) then
        ! Yes, so restore the old solution and
        ! repeat the adaptive time-stepping loop
        call lsysbl_copyVector(p_rsolutionInitial, rsolution)
      else
        ! No, accept current solution and
        ! exit adaptive time-stepping loop
        exit timeadapt
      end if
    end do timeadapt
  end subroutine tstep_performRKStepBl

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performRKStepScCpl(rproblemLevel, rtimestep,&
      rsolver, Rsolution, fcb_nlsolverCallback, rcollection, Rsource)

!<description>
    ! This subroutine performs one step by means of an explicit
    ! Runge-Kutta scheme to advance the solution from time level $t^n$
    ! to time level $t^{n+1}$.  Note that this subroutine serves as
    ! wrapper for scalar solution vectors only.
    !
    ! In contrast to routine tstep_performRKStepSc this routine
    ! is applicable to a sequence of coupled problems which are
    ! solved repeatedly until convergence.
!</description>

!<input>
    ! Callback routines
    include 'intf_solvercallback.inc'

    ! OPTIONAL: constant right-hand side vector
    type(t_vectorScalar), dimension(:), intent(in), optional :: Rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! solution vector
    type(t_vectorScalar), dimension(:), intent(inout) :: Rsolution

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), dimension(:), allocatable :: RsolutionBlock
    type(t_vectorBlock), dimension(:), allocatable :: RsourceBlock
    integer :: icomponent

    if (present(Rsource)) then

      ! Allocate temporal memory
      allocate(RsolutionBlock(size(Rsolution)))
      allocate(RsourceBlock(size(Rsource)))

      ! Convert scalar vectors to 1-block vectors
      do icomponent = 1, size(Rsolution)
        call lsysbl_createVecFromScalar(Rsolution(icomponent),&
            RsolutionBlock(icomponent))
      end do

      ! Convert scalar vectors to 1-block vectors
      do icomponent = 1, size(Rsource)
        call lsysbl_createVecFromScalar(Rsource(icomponent),&
            RsourceBlock(icomponent))
      end do

      ! Call block version of this subroutine
      call tstep_performRKStepBlCpl(rproblemLevel, rtimestep,&
          rsolver, RsolutionBlock, fcb_nlsolverCallback, rcollection,&
          RsourceBlock)

      ! Deallocate temporal 1-block vectors
      do icomponent = 1, size(Rsolution)
        call lsysbl_releaseVector(RsolutionBlock(icomponent))
      end do

      ! Deallocate temporal 1-block vectors
      do icomponent = 1, size(Rsource)
        call lsysbl_releaseVector(RsourceBlock(icomponent))
      end do

      ! Release temporal memory
      deallocate(RsolutionBlock)
      deallocate(RsourceBlock)

    else

      ! Allocate temporal memory
      allocate(RsolutionBlock(size(Rsolution)))

      ! Convert scalar vectors to 1-block vectors
      do icomponent = 1, size(Rsolution)
        call lsysbl_createVecFromScalar(Rsolution(icomponent),&
            RsolutionBlock(icomponent))
      end do

      ! Call block version of this subroutine
      call tstep_performRKStepBlCpl(rproblemLevel, rtimestep,&
          rsolver, RsolutionBlock, fcb_nlsolverCallback, rcollection)

      ! Deallocate temporal 1-block vectors
      do icomponent = 1, size(Rsolution)
        call lsysbl_releaseVector(RsolutionBlock(icomponent))
      end do

      ! Release temporal memory
      deallocate(RsolutionBlock)

    end if
    
  end subroutine tstep_performRKStepScCpl

! *****************************************************************************

!<subroutine>

  subroutine tstep_performRKStepBlCpl(rproblemLevel, rtimestep,&
      rsolver, Rsolution,  fcb_nlsolverCallback, rcollection, Rsource)

!<description>
    ! This subroutine performs one step by means of an explicit
    ! Runge-Kutta scheme to advance the solution from time level $t^n$
    ! to time level $t^{n+1}$.  In particular, the two step version
    ! suggested by Richtmyer and Morton is implemented so that no
    ! artificial diffusion needs to be constructed.
    !
    ! In contrast to routine tstep_performRKStepBl this routine
    ! is applicable to a sequence of coupled problems which are
    ! solved repeatedly until convergence.
!</description>

!<input>
    ! Callback routines
    include 'intf_solvercallback.inc'

    ! OPTIONAL: constant right-hand side vector
    type(t_vectorBlock), dimension(:), intent(in), optional :: Rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(inout), target :: rtimestep

    ! solver structure
    type(t_solver), intent(inout), target :: rsolver

    ! solution vector
    type(t_vectorBlock), dimension(:), intent(inout) :: Rsolution

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    print *, "Subroutine tstep_performRKStepBlCpl is not implemented"
    stop

  end subroutine tstep_performRKStepBlCpl

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performPseudoSteppingBl(rproblemLevel,&
      rtimestep, rsolver, rsolution, fcb_nlsolverCallback,&
      rcollection, rsource)

!<description>
    ! This subroutine performs pseudo time-stepping to compute the
    ! steady state solution to a stationary problem. Each pseudo
    ! time step can be performed by the forward Euler (theta=0),
    ! backward Euler (theta=1), the Crank-Nicolson (theta=0.5)
    ! time stepping algorithm and an explicit Runge-Kutta scheme.
!</description>

!<input>
    ! Callback routines
    include 'intf_solvercallback.inc'

    ! OPTIONAL: source vector
    type(t_vectorBlock), intent(in), optional :: rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver structure
    type(t_solver), intent(inout), target :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>


    ! Infinite time loop
    timeloop: do

      ! What time-stepping scheme should be used?
      select case(rtimestep%ctimestepType)

      case (TSTEP_RK_SCHEME)

        ! Adopt explicit Runge-Kutta scheme
        call tstep_performRKStep(rproblemLevel, rtimestep, rsolver,&
            rsolution, fcb_nlsolverCallback, rcollection, rsource)

      case (TSTEP_THETA_SCHEME)

        ! Adopt two-level theta-scheme
        call tstep_performThetaStep(rproblemLevel, rtimestep, rsolver,&
            rsolution, fcb_nlsolverCallback, rcollection, rsource)

      case default
        if (rtimestep%coutputModeError .gt. 0) then
          call output_line('Unsupported time-stepping algorithm!',&
              OU_CLASS_ERROR,rtimestep%coutputModeError,&
              'tstep_performPseudoSteppingBl')
        end if
        call sys_halt()
      end select

      ! Reached final time, then exit infinite time loop?
      if (rtimestep%dTime .ge. rtimestep%dfinalTime) exit timeloop

      ! Reached steady state limit?
      if (rtimestep%depsSteady > 0.0_DP) then

        ! Check if steady-state residual exceeds tolerance
        if ((rsolver%dfinalDefect   < rsolver%dinitialDefect) .and.&
            (rsolver%dinitialDefect < rtimestep%dStep*rtimestep%depsSteady)) exit timeloop
      end if

    end do timeloop

  end subroutine tstep_performPseudoSteppingBl

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performPseudoSteppingSc(rproblemLevel,&
      rtimestep, rsolver, rsolution, fcb_nlsolverCallback,&
      rcollection, rsource)

!<description>
    ! This subroutine performs pseudo time-stepping to compute the
    ! steady state solution to a stationary problem. Each pseudo
    ! time step can be performed by the forward Euler (theta=0),
    ! backward Euler (theta=1), the Crank-Nicolson (theta=0.5)
    ! time stepping algorithm and an explicit Runge-Kutta scheme.
!</description>

!<input>
    ! Callback routines
    include 'intf_solvercallback.inc'

    ! OPTIONAL: source vector
    type(t_vectorScalar), intent(in), optional :: rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver structure
    type(t_solver), intent(inout), target :: rsolver

    ! solution vector
    type(t_vectorScalar), intent(inout) :: rsolution

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! Infinite time loop
    timeloop: do

      ! What time-stepping scheme should be used?
      select case(rtimestep%ctimestepType)

      case (TSTEP_RK_SCHEME)

        ! Adopt explicit Runge-Kutta scheme
        call tstep_performRKStep(rproblemLevel, rtimestep, rsolver,&
            rsolution, fcb_nlsolverCallback, rcollection, rsource)

      case (TSTEP_THETA_SCHEME)

        ! Adopt two-level theta-scheme
        call tstep_performThetaStep(rproblemLevel, rtimestep, rsolver,&
            rsolution, fcb_nlsolverCallback, rcollection, rsource)

      case default
        if (rtimestep%coutputModeError .gt. 0) then
          call output_line('Unsupported time-stepping algorithm!',&
              OU_CLASS_ERROR,rtimestep%coutputModeError,&
              'tstep_performPseudoSteppingSc')
        end if
        call sys_halt()
      end select

      ! Reached final time, then exit infinite time loop?
      if (rtimestep%dTime .ge. rtimestep%dfinalTime) exit timeloop

      ! Reached steady state limit?
      if (rtimestep%depsSteady > 0.0_DP) then

        ! Check if steady-state residual exceeds tolerance
        if ((rsolver%dfinalDefect   < rsolver%dinitialDefect) .and.&
            (rsolver%dinitialDefect < rtimestep%dStep*rtimestep%depsSteady)) exit timeloop
      end if

    end do timeloop

  end subroutine tstep_performPseudoSteppingSc


  ! *****************************************************************************
  ! AUXILIARY ROUTINES
  ! *****************************************************************************


!<function>

  function tstep_checkTimestep(rtimestep, rsolver, rsolution1,&
      rsolution2) result(breject)

!<description>
    ! This functions checks the result of the current time step and
    ! performs adaptive time-step controlling if this is required
    !
    ! First, it checks if the solution failed due to some critical error.
    ! If this is the case, the time step is reduced unless the smallest
    ! admissible time step size has been reached. In this case, the
    ! simulation is stopped unconditionally since a critical error occured.
    !
    ! If the time step size should be computed adaptively the corresponding
    ! time-stepping algorithm is used to adjust the time step.
!</description>

!<input>
    ! solver
    type(t_solver), intent(in):: rsolver

    ! first solution vector
    type(t_vectorBlock), intent(in) :: rsolution1

    ! second solution vector
    type(t_vectorBlock), intent(in) :: rsolution2
!</input>

!<inputoutput>
    ! time-stepping algorithm
    type(t_timestep), intent(inout) :: rtimestep
!</inputoutput>

!<result>
    ! reject solution
    logical :: breject
!</result>
!</function>

    ! local variables
    type(t_pidController), pointer  :: p_rpidController
    type(t_autoController), pointer :: p_rautoController
    type(t_serController), pointer  :: p_rserController
    real(DP), dimension(:), pointer :: p_Ddata1, p_Ddata2
    real(DP) :: dChange, dStepOpt, dPaux, dIaux, dDaux

    !---------------------------------------------------------------------------
    ! If the solver failed due to a critical error, we must reject the
    ! time step and recompute the solution adopting a smaller time step
    !---------------------------------------------------------------------------
    breject = (rsolver%istatus .eq. SV_INF_DEF  .or.&
               rsolver%istatus .eq. SV_INCR_DEF .or.&
               rsolver%istatus .eq. SV_NAN_RHS  .or.&
               rsolver%istatus .eq. SV_NAN_DEF)

    if (breject) then

      ! If the time step is already equal to the smallest
      ! admissible time step, then the simulation is terminated.
      if (rtimestep%dStep .le. rtimestep%dminStep + SYS_EPSREAL_DP) then
        if (rtimestep%coutputModeError .gt. 0) then
          call output_line('Time step reached smallest admissible value!',&
              OU_CLASS_ERROR,rtimestep%coutputModeError,&
              'tstep_checkTimestep')
        end if
        call sys_halt()
      end if

      ! Undo the last time step and adjust time step size
      rtimestep%nrejectedSteps = rtimestep%nrejectedSteps+1
      rtimestep%nSteps         = rtimestep%nSteps-1
      rtimestep%dTime          = rtimeStep%dTime-rtimeStep%dStep
      rtimestep%dStep          = max(rtimestep%dminStep,&
                                     rtimestep%dstepReductionFactor*rtimestep%dStep)

      ! That is all, recompute the last step
      return
    end if


    !---------------------------------------------------------------------------
    ! Perform adaptive time-stepping?
    !---------------------------------------------------------------------------

    select case(rtimestep%iadaptTimestep)

    case (TSTEP_AUTOADAPT)
      !-------------------------------------------------------------------------
      ! Automatic time step control
      !
      ! The vector rsolution1 has been computed by performing ONE time step
      ! of size Dt and vector rsolution2 has been computed by performing TWO
      ! time steps of size Dt/2.
      !
      ! The 'optimal' time step is computed from the following formula:
      !
      !    Dt_opt = SQRT( TOL * (Dt^2*(m^2-1)) / !! u1-u2 !! )
      !
      ! If the time step decreased too much, then the computed solution is
      ! rejected. Otherwise, it is accepted and the new time step is adopted.
      !-------------------------------------------------------------------------

      ! Set pointer
      p_rautoController => rtimestep%p_rautoController

      ! Compute norm of solution changes
      call lsysbl_getbase_double(rsolution1, p_Ddata1)
      call lsysbl_getbase_double(rsolution2, p_Ddata2)
      dChange = lalg_errorNorm(p_Ddata1, p_Ddata2, rtimestep%isolNorm)

      ! Calculate the 'optimal' time step size
      dStepOpt = sqrt(3.0_DP*p_rautoController%depsRel*rtimestep%dStep/dChange)

      ! Check if optimal time step decreased too much
      if (dStepOpt .le. p_rautoController%dDecreaseFactor*rtimestep%dStep) then

        ! Adjust time step accordingly
        rtimestep%nrejectedSteps = rtimestep%nrejectedSteps+1
        rtimestep%nSteps         = rtimestep%nSteps-1
        rtimestep%dTime          = rtimeStep%dTime-rtimeStep%dStep
        rtimestep%dStep          = max(rtimestep%dminStep, dStepOpt)

        ! That is it, recompute the last step
        breject = .true.
      else

        ! Accept the current solution
        breject = .false.

        ! Impose upper/lower bounds on absolute values
        rtimestep%dStep1 = rtimestep%dStep
        rtimestep%dStep  = max(rtimestep%dminStep, min(rtimestep%dmaxStep, dStepOpt))

        ! Calculate the relative changes for statistical information
        rtimestep%drelChange = dChange/max(SYS_EPSREAL_DP,&
                                           lalg_norm(p_Ddata1, rtimestep%isolNorm))

      end if


    case (TSTEP_PIDADAPT)
      !-------------------------------------------------------------------------
      ! Evolutionary PID controller
      !
      ! The vectors rvector1 and rvector2 represent the solution from the
      ! current and the previous time step. The relative changes
      !
      !    dChange = !!u1-u2!! / !!u1!!
      !
      ! are computed and the current solution is rejected if dChange > epsRel.
      ! In this case, the time step is recomputed with the new time step
      !
      !    Dt_opt = epsRel/dChange*Dt
      !
      ! If the computed solution is accepted, then the new time step is
      !
      !    Dt_opt = (e_{n-1}/e_n)^kp * (depsRel/e_n)^ki * (e_n^2/e_n/e_{n-2})^kd * Dt
      !
      ! where {e_n-i} denotes the relative changes of the i-th previous step.
      !
      !-------------------------------------------------------------------------

      ! Set pointer
      p_rpidController => rtimestep%p_rpidController

      ! Compute norm of solution changes
      call lsysbl_getbase_double(rsolution1, p_Ddata1)
      call lsysbl_getbase_double(rsolution2, p_Ddata2)
      dChange = lalg_errorNorm(p_Ddata1, p_Ddata2, rtimestep%isolNorm) /&
                     lalg_norm(p_Ddata1, rtimestep%isolNorm)
      rtimestep%drelChange = dChange

      ! Check if solver did not converge, that is, the convergence rate equals
      ! unity or the normalised changes exceed maximum tolerance
      breject = (dChange .gt. p_rpidController%dmaxRel) .or.&
                (rsolver%istatus .eq. SV_DIVERGED)

      ! If the solver failed due to a non-critical error, we may accept
      ! the time step in case it cannot be further reduces
      if (breject .and. (rtimestep%dStep .gt. rtimestep%dminStep + SYS_EPSREAL_DP)) then

        ! Adjust time step accordingly
        rtimestep%nrejectedSteps = rtimestep%nrejectedSteps+1
        rtimestep%nSteps         = rtimestep%nSteps-1
        rtimestep%dTime          = rtimeStep%dTime-rtimeStep%dStep

        ! Ensure the time step is reduced by a significant factor
        dStepOpt = min(rtimestep%dStepReductionFactor,&
                       p_rpidController%dmaxRel/dChange) * rtimestep%dStep

        ! Adjust time step and recompute last step
        rtimestep%dStep = max(rtimestep%dminStep, dStepOpt)
      else

        ! Accept the current solution
        breject = .false.

        ! Adopt previous time step if solution did not change
        if (dChange .le. SYS_EPSREAL_DP) return

        if (rtimestep%dTime .gt. rtimestep%dadaptTime) then

          ! Calculate the 'optimal' time step size
          dPaux = (p_rpidController%dcontrolValue1/dChange)**p_rpidController%dProportionalExponent
          dIaux = (p_rpidController%depsRel/dChange)**p_rpidController%dIntegralExponent
          dDaux = (p_rpidController%dcontrolValue1*p_rpidController%dcontrolValue1/&
                      dChange/p_rpidController%dcontrolValue2)**p_rpidController%dDerivativeExponent
          dStepOpt = dPaux*dIaux*dDaux*rtimestep%dStep

          ! Limit the growth and reduction of the time step
          dStepOpt = max(p_rpidController%dDecreaseFactor*rtimestep%dStep1,&
                         min(p_rpidController%dIncreaseFactor*rtimestep%dStep1,&
                             dStepOpt))
        else

          ! Adopt value from previous time step
          dstepOpt = rtimestep%dStep
        end if

        ! Impose upper/lower bounds on absolute values
        rtimestep%dStep1 = rtimestep%dStep
        rtimestep%dStep  = max(rtimestep%dminStep, min(rtimestep%dmaxStep, dStepOpt))

        ! Update history of relative changes
        p_rpidController%dcontrolValue2 = p_rpidController%dcontrolValue1
        p_rpidController%dcontrolValue1 = dChange

      end if


    case (TSTEP_SERADAPT)
      !-------------------------------------------------------------------------
      ! Switched evolution relaxation (SER)
      !
      ! From: W. Mulder and B. Van Leer. `Experiments with typical upwind
      ! methods for the Euler equations`, J. Comp. Phys., 59(1985), 232-246.
      ! The actual relaxation of the time step is performed by the
      ! relaxation subroutine which is called from within the nonlinear
      ! solver. Here, we only update data for the next time step.
      !-------------------------------------------------------------------------

      ! Set pointer
      p_rserController => rtimestep%p_rserController

      ! Accept the current solution
      breject = .false.

      ! Store stationary defect from previous loop
      p_rserController%dsteadyDefect1 = p_rserController%dsteadyDefect

      ! Store time step size from previous loop
      rtimestep%dStep1 = rtimestep%dStep

      ! Calculate the relative changes for statistical information
      call lsysbl_getbase_double(rsolution1, p_Ddata1)
      call lsysbl_getbase_double(rsolution2, p_Ddata2)
      dChange = lalg_errorNorm(p_Ddata1, p_Ddata2, rtimestep%isolNorm)
      rtimestep%drelChange = dChange/max(SYS_EPSREAL_DP,&
                                     lalg_norm(p_Ddata1, rtimestep%isolNorm))


    case default
      ! Accept time step unconditionally
      breject = .false.

      ! Calculate the relative changes for statistical information
      call lsysbl_getbase_double(rsolution1, p_Ddata1)
      call lsysbl_getbase_double(rsolution2, p_Ddata2)
      dChange = lalg_errorNorm(p_Ddata1, p_Ddata2, rtimestep%isolNorm)
      rtimestep%drelChange = dChange/max(SYS_EPSREAL_DP,&
                                         lalg_norm(p_Ddata1, rtimestep%isolNorm))
    end select

  end function tstep_checkTimestep

  ! *****************************************************************************

!<function>

  function tstep_checkTimestepCpl(rtimestep, rsolver, Rsolution1,&
      Rsolution2) result(breject)

!<description>
    ! This functions checks the result of the current time step and
    ! performs adaptive time-step controlling if this is required
    !
    ! First, it checks if the solution failed due to some critical error.
    ! If this is the case, the time step is reduced unless the smallest
    ! admissible time step size has been reached. In this case, the
    ! simulation is stopped unconditionally since a critical error occured.
    !
    ! If the time step size should be computed adaptively the corresponding
    ! time-stepping algorithm is used to adjust the time step.
!</description>

!<input>
    ! solver
    type(t_solver), intent(in):: rsolver

    ! first solution vector
    type(t_vectorBlock), dimension(:), intent(in) :: Rsolution1

    ! second soluion vector
    type(t_vectorBlock), dimension(:), intent(in) :: Rsolution2
!</input>

!<inputoutput>
    ! time-stepping algorithm
    type(t_timestep), intent(inout) :: rtimestep
!</inputoutput>

!<result>
    ! reject solution
    logical :: breject
!</result>
!</function>

    ! local variable
    type(t_timestep) :: rtimestep0, rtimestep1
    integer :: icomponent

    breject = .false.

    ! Make a duplicate of the given time step structure
    rtimestep0 = rtimestep

    ! Loop over all componenets
    do icomponent = 1, size(rsolver%p_rsolverSubnode)

      ! Adopt values from original time step structure
      rtimestep1 = rtimestep0

      ! Check if time step can be accepted
      breject = breject .or. tstep_checkTimestep(&
          rtimestep1, rsolver%p_rsolverSubnode(icomponent),&
          Rsolution1(icomponent), Rsolution2(icomponent))

      ! Update original time step
      rtimestep%nrejectedSteps = rtimestep1%nrejectedSteps

      if (icomponent .eq. 1) then
        rtimestep%dTime      = rtimestep1%dTime
        rtimestep%dStep      = rtimestep1%dStep
        rtimestep%dStep1     = rtimestep1%dStep1
        rtimestep%drelChange = rtimestep1%drelChange
      elseif (rtimestep1%dStep .lt. rtimestep%dStep) then
        rtimestep%dTime      = rtimestep1%dTime
        rtimestep%dStep      = rtimestep1%dStep
        rtimestep%dStep1     = rtimestep1%dStep1
        rtimestep%drelChange = rtimestep1%drelChange
      end if

    end do

  end function tstep_checkTimestepCpl

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_decodeOutputLevel(rtimestep)

!<description>
    ! This subroutine decodes the output level information into bitfields
!</description>

!<inputoutput>
    type(t_timestep), intent(inout) :: rtimestep
!</inputoutput>
!</subroutine>

    ! local variable
    character(len=3) :: coutputLevel
    integer :: i

    ! Initialisation
    rtimestep%coutputModeError   = 0_I32
    rtimestep%coutputModeWarning = 0_I32
    rtimestep%coutputModeInfo    = 0_I32
    rtimestep%coutputModeVerbose = 0_I32

    coutputLevel = sys_i03(rtimestep%ioutputLevel)
    
    ! Check output level for benchmark log file
    read(coutputLevel(1:1),*) i
    select case(i)
    case (TSTEP_IOLEVEL_VERBOSE)
      rtimestep%coutputModeVerbose = rtimestep%coutputModeVerbose + OU_MODE_BENCHLOG
    case (TSTEP_IOLEVEL_INFO)
      rtimestep%coutputModeInfo    = rtimestep%coutputModeInfo    + OU_MODE_BENCHLOG
    case (TSTEP_IOLEVEL_WARNING)
      rtimestep%coutputModeWarning = rtimestep%coutputModeWarning + OU_MODE_BENCHLOG
    case (TSTEP_IOLEVEL_ERROR)
      rtimestep%coutputModeError   = rtimestep%coutputModeError   + OU_MODE_BENCHLOG
    case (TSTEP_IOLEVEL_SILENT)
    case default
      call output_line('Invalid output level.',&
          OU_CLASS_MSG, OU_MODE_STD,'tstep_decodeOutputLevel')
    end select

    ! Check output level for terminal
    read(coutputLevel(2:2),*) i
    select case(i)
    case (TSTEP_IOLEVEL_VERBOSE)
      rtimestep%coutputModeVerbose = rtimestep%coutputModeVerbose + OU_MODE_TERM
    case (TSTEP_IOLEVEL_INFO)
      rtimestep%coutputModeInfo    = rtimestep%coutputModeInfo    + OU_MODE_TERM
    case (TSTEP_IOLEVEL_WARNING)
      rtimestep%coutputModeWarning = rtimestep%coutputModeWarning + OU_MODE_TERM
    case (TSTEP_IOLEVEL_ERROR)
      rtimestep%coutputModeError   = rtimestep%coutputModeError   + OU_MODE_TERM
    case (TSTEP_IOLEVEL_SILENT)
    case default
      call output_line('Invalid output level.',&
          OU_CLASS_MSG, OU_MODE_STD,'tstep_decodeOutputLevel')
    end select

    ! Check output level for log file
    read(coutputLevel(3:3),*) i
    select case(i)
    case (TSTEP_IOLEVEL_VERBOSE)
      rtimestep%coutputModeVerbose = rtimestep%coutputModeVerbose + OU_MODE_LOG
    case (TSTEP_IOLEVEL_INFO)
      rtimestep%coutputModeInfo    = rtimestep%coutputModeInfo    + OU_MODE_LOG
    case (TSTEP_IOLEVEL_WARNING)
      rtimestep%coutputModeWarning = rtimestep%coutputModeWarning + OU_MODE_LOG
    case (TSTEP_IOLEVEL_ERROR)
      rtimestep%coutputModeError   = rtimestep%coutputModeError   + OU_MODE_LOG
    case (TSTEP_IOLEVEL_SILENT)
    case default
      call output_line('Invalid output level.',&
          OU_CLASS_MSG, OU_MODE_STD,'tstep_decodeOutputLevel')
    end select

    ! Adjust lower log levels
    rtimestep%coutputModeInfo    = ior(rtimestep%coutputModeInfo,&
                                       rtimestep%coutputModeVerbose)
    rtimestep%coutputModeWarning = ior(rtimestep%coutputModeWarning,&
                                       rtimestep%coutputModeInfo)
    rtimestep%coutputModeError   = ior(rtimestep%coutputModeError,&
                                       rtimestep%coutputModeWarning)

  end subroutine tstep_decodeOutputLevel

end module timestep
