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
!# Furthermore, each object of type t_timestep provides information
!# about convergence and norms after the time-stepping algorithm has
!# terminated.
!#
!# The following routines are available:
!#
!# 1.) tstep_createTimestep = tstep_createTimestepFromParlist /
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
!# 5.) tstep_performTimestep = tstep_performTimestepSc /
!#                             tstep_performTimestepScCpl /
!#                             tstep_performTimestepBl /
!#                             tstep_performTimestepBlCpl
!#     -> Performs one time step with the given time-stepping structure
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
  use solverbase
  use solverlinear
  use solvernonlinear
  use timestepbase

  implicit none
 
  private
  public :: tstep_createTimestep
  public :: tstep_releaseTimestep
  public :: tstep_infoTimestep
  public :: tstep_resetTimestep
  public :: tstep_performTimestep

  ! *****************************************************************************

  interface tstep_createTimestep
    module procedure tstep_createTimestepFromParlist
    module procedure tstep_createTimestepIndirect
  end interface

  interface tstep_performTimestep
    module procedure tstep_performTimestepSc
    module procedure tstep_performTimestepBl
    module procedure tstep_performTimestepScCpl
    module procedure tstep_performTimestepBlCpl
  end interface
  
  interface tstep_checkTimestep
    module procedure tstep_checkTimestep
    module procedure tstep_checkTimestepCpl
  end interface

  ! *****************************************************************************

contains
  
  ! *****************************************************************************

!<subroutine>

  subroutine tstep_createTimestepFromParlist(rparlist, ssectionName,&
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
    character(len=SYS_STRLEN) :: stimestepName
    integer :: npc
    logical :: bcheck

    ! Set number of coupled problems
    npc = 1
    if (present(nproblemsCoupled)) npc = nproblemsCoupled
    
    ! The INTENT(out) already initialises rtimestep with the most
    ! important information.
    rtimestep%sName = trim(adjustl(ssectionName))
    
    ! Get time-stepping algorithm (family + scheme type)
    call parlst_getvalue_string(rparlist, ssectionName,&
        "ctimestep", stimestepName)
    rtimestep%ctimestep   = tstep_igetID(stimestepName)
    
    ! Get mandatory configuration values from parameter list
    call parlst_getvalue_double(rparlist, ssectionName,&
        "dfinalTime", rtimestep%dfinalTime)
    call parlst_getvalue_double(rparlist, ssectionName,&
        "dinitialStep", rtimestep%dinitialStep)
    call parlst_getvalue_double(rparlist, ssectionName,&
        "dinitialTime", rtimestep%dinitialTime)
    call parlst_getvalue_double(rparlist, ssectionName,&
        "dstepReductionFactor", rtimestep%dstepReductionFactor)
    call parlst_getvalue_double(rparlist, ssectionName,&
        "depsSteady", rtimestep%depsSteady)
    call parlst_getvalue_double(rparlist, ssectionName,&
        "dadaptTime", rtimestep%dadaptTime)
    call parlst_getvalue_int(rparlist, ssectionName,&
        "ioutputlevel", rtimestep%ioutputlevel)
    call parlst_getvalue_int(rparlist, ssectionName,&
        "isolNorm", rtimestep%isolNorm)
    call parlst_getvalue_int(rparlist, ssectionName,&
        "iadapttimestep", rtimestep%iadapttimestep)
    
    ! Get optional configuration values from parameter list
    call parlst_getvalue_double(rparlist, ssectionName,&
        "dminStep", rtimestep%dminStep, rtimestep%dinitialStep)
    call parlst_getvalue_double(rparlist, ssectionName,&
        "dmaxStep", rtimestep%dmaxStep, rtimestep%dinitialStep)

    ! Decode the output level
    call tstep_decodeOutputLevel(rtimestep)

    ! Get solver dependent configuration values from parameter list
    ! and allocate array of temporal vectors
    select case(tstep_getFamily(rtimestep%ctimestep))
    case (TSTEP_THETA)
      ! Two-level theta scheme family
      allocate(rtimestep%p_rthetaScheme)
      
      select case(tstep_getType(rtimestep%ctimestep))
      case (TSTEP_FORWARD_EULER)
        rtimestep%p_rthetaScheme%theta = 0.0_DP
        
      case (TSTEP_BACKWARD_EULER)
        rtimestep%ctimestep = ior(rtimestep%ctimestep,TSTEP_IMPLICIT)
        rtimestep%p_rthetaScheme%theta = 1.0_DP
        
      case (TSTEP_CRANK_NICHOLSON)
        rtimestep%ctimestep = ior(rtimestep%ctimestep,TSTEP_IMPLICIT)
        rtimestep%p_rthetaScheme%theta = 0.5_DP
        
      case (TSTEP_THETA_GEN)
        ! Get theta value from from parameter list
        call parlst_getvalue_double(rparlist, ssectionName,&
            "theta", rtimestep%p_rthetaScheme%theta)

        if (rtimestep%p_rthetaScheme%theta .ne. 0.0) then
          rtimestep%ctimestep = ior(rtimestep%ctimestep,TSTEP_IMPLICIT)
        end if
                
      case (TSTEP_FSTHETA_GEN)
        ! Get theta, and alpha values from from parameter list
        call parlst_getvalue_double(rparlist, ssectionName,&
            "theta", rtimestep%p_rthetaScheme%theta)
        rtimestep%p_rthetaScheme%gamma = 1.0_DP-2.0_DP*rtimestep%p_rthetaScheme%theta
        rtimestep%p_rthetaScheme%alpha = rtimestep%p_rthetaScheme%gamma/(1.0_DP-rtimestep%p_rthetaScheme%theta)
        call parlst_getvalue_double(rparlist, ssectionName,&
            "alpha", rtimestep%p_rthetaScheme%alpha, rtimestep%p_rthetaScheme%alpha)

        if (rtimestep%p_rthetaScheme%alpha .ne. 0.0) then
          rtimestep%ctimestep = ior(rtimestep%ctimestep,TSTEP_IMPLICIT)
        end if
        
      case (TSTEP_FSTHETA_2)
        rtimestep%ctimestep = ior(rtimestep%ctimestep,TSTEP_IMPLICIT)
        rtimestep%p_rthetaScheme%theta = 1.0_DP-sqrt(2.0_DP)/2.0_DP
        rtimestep%p_rthetaScheme%gamma = 1.0_DP-2.0_DP*rtimestep%p_rthetaScheme%theta
        rtimestep%p_rthetaScheme%alpha = rtimestep%p_rthetaScheme%gamma/(1.0_DP-rtimestep%p_rthetaScheme%theta)
        
      case default
        if (rtimestep%coutputModeError .gt. 0) then
          call output_line('Invalid type of time stepping algorithm!',&
              OU_CLASS_ERROR,rtimestep%coutputModeError,&
              'tstep_createTimestepFromParlist')
        end if
        call sys_halt()
      end select
      
      ! Allocate array of temporal vectors
      allocate(rtimestep%p_rthetaScheme%RtempVectors(npc))
      
    case (TSTEP_RUNGE_KUTTA)
      ! Multi-stage Runge-Kutta time-stepping scheme
      allocate(rtimestep%p_rRungeKuttaScheme)

      ! Try Butcher tableau
      call tstep_genButcherTableau(rtimestep%ctimestep,&
          rtimestep%p_rRungeKuttaScheme%Da,&
          rtimestep%p_rRungeKuttaScheme%Db,&
          rtimestep%p_rRungeKuttaScheme%Dc,&
          .true.,bcheck)

      if (.not.bcheck) then
        ! Try alpha-beta tableau
        call tstep_genAlphaBetaTableau(rtimestep%ctimestep,&
          rtimestep%p_rRungeKuttaScheme%Dalpha,&
          rtimestep%p_rRungeKuttaScheme%Dbeta,&
          rtimestep%p_rRungeKuttaScheme%Dc,&
          .true.,bcheck)

        if (.not.bcheck) then
          if (rtimestep%coutputModeError .gt. 0) then
            call output_line('Neither Butcher nor alpha-beta tableau available!',&
                OU_CLASS_ERROR,rtimestep%coutputModeError,&
                'tstep_createTimestepFromParlist')
          end if
          call sys_halt()
        end if
      end if

      ! Get number of stages
      rtimestep%p_rRungeKuttaScheme%nstages = size(rtimestep%p_rRungeKuttaScheme%Dc)
      
      ! Allocate array of temporal vectors
      allocate(rtimestep%p_rRungeKuttaScheme%RtempVectors(npc*rtimestep%p_rRungeKuttaScheme%nstages))
      
    case default
      if (rtimestep%coutputModeError .gt. 0) then
        call output_line('Invalid family of time stepping algorithm!',&
            OU_CLASS_ERROR,rtimestep%coutputModeError,&
            'tstep_createTimestepFromParlist')
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
            'tstep_createTimestepFromParlist')
      end if
      call sys_halt()
    end select

    ! Reset any other values
    call tstep_resetTimestep(rtimestep, .true.)
    
  end subroutine tstep_createTimestepFromParlist

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

    ! Copy data for time stepping scheme
    select case(tstep_getFamily(rtimestep%ctimestep))
    case (TSTEP_THETA)
      ! Create two-level theta scheme
      allocate(rtimestep%p_rthetaScheme)
      rtimestep%p_rthetaScheme = rtimestepTemplate%p_rthetaScheme

      ! Allocate memory for temporal vectors
      allocate(rtimestep%p_rthetaScheme%RtempVectors(&
          lbound(rtimestepTemplate%p_rthetaScheme%RtempVectors,1):&
          ubound(rtimestepTemplate%p_rthetaScheme%RtempVectors,1)))

    case (TSTEP_RUNGE_KUTTA)    
      ! Create multi-step Runge Kutta time-stepping scheme
      allocate(rtimestep%p_rRungeKuttaScheme)
      rtimestep%p_rRungeKuttaScheme = rtimestepTemplate%p_rRungeKuttaScheme

      ! Allocate memory for temporal memory
      allocate(rtimestep%p_rRungeKuttaScheme%RtempVectors(&
          lbound(rtimestepTemplate%p_rRungeKuttaScheme%RtempVectors,1):&
          ubound(rtimestepTemplate%p_rRungeKuttaScheme%RtempVectors,1)))

      if (iand(rtimestep%ctimestep,TSTEP_BUTCHER_TABLEAU) .eq. TSTEP_BUTCHER_TABLEAU) then
        ! Allocate memory for Butcher tableau
        allocate(rtimestep%p_rRungeKuttaScheme%Da(&
            lbound(rtimestepTemplate%p_rRungeKuttaScheme%Da,1):&
            ubound(rtimestepTemplate%p_rRungeKuttaScheme%Da,1),&
            lbound(rtimestepTemplate%p_rRungeKuttaScheme%Da,2):&
            ubound(rtimestepTemplate%p_rRungeKuttaScheme%Da,2)))
        allocate(rtimestep%p_rRungeKuttaScheme%Db(&
            lbound(rtimestepTemplate%p_rRungeKuttaScheme%Db,1):&
            ubound(rtimestepTemplate%p_rRungeKuttaScheme%Db,1)))
        allocate(rtimestep%p_rRungeKuttaScheme%Dc(&
            lbound(rtimestepTemplate%p_rRungeKuttaScheme%Dc,1):&
            ubound(rtimestepTemplate%p_rRungeKuttaScheme%Dc,1)))

        ! Copy coefficients
        rtimestep%p_rRungeKuttaScheme%Da = rtimestepTemplate%p_rRungeKuttaScheme%Da
        rtimestep%p_rRungeKuttaScheme%Db = rtimestepTemplate%p_rRungeKuttaScheme%Db
        rtimestep%p_rRungeKuttaScheme%Dc = rtimestepTemplate%p_rRungeKuttaScheme%Dc
      else
        ! Allocate memory for alpha-beta form
        allocate(rtimestep%p_rRungeKuttaScheme%Dalpha(&
            lbound(rtimestepTemplate%p_rRungeKuttaScheme%Dalpha,1):&
            ubound(rtimestepTemplate%p_rRungeKuttaScheme%Dalpha,1),&
            lbound(rtimestepTemplate%p_rRungeKuttaScheme%Dalpha,2):&
            ubound(rtimestepTemplate%p_rRungeKuttaScheme%Dalpha,2)))
        allocate(rtimestep%p_rRungeKuttaScheme%Dbeta(&
            lbound(rtimestepTemplate%p_rRungeKuttaScheme%Dbeta,1):&
            ubound(rtimestepTemplate%p_rRungeKuttaScheme%Dbeta,1),&
            lbound(rtimestepTemplate%p_rRungeKuttaScheme%Dbeta,2):&
            ubound(rtimestepTemplate%p_rRungeKuttaScheme%Dbeta,2)))
        allocate(rtimestep%p_rRungeKuttaScheme%Dc(&
            lbound(rtimestepTemplate%p_rRungeKuttaScheme%Dc,1):&
            ubound(rtimestepTemplate%p_rRungeKuttaScheme%Dc,1)))
        
        ! Copy coefficients
        rtimestep%p_rRungeKuttaScheme%Dalpha = rtimestepTemplate%p_rRungeKuttaScheme%Dalpha
        rtimestep%p_rRungeKuttaScheme%Dbeta  = rtimestepTemplate%p_rRungeKuttaScheme%Dbeta
        rtimestep%p_rRungeKuttaScheme%Dc     = rtimestepTemplate%p_rRungeKuttaScheme%Dc
      end if

    case default
      if (rtimestep%coutputModeError .gt. 0) then
        call output_line('Invalid type of time-stepping algorithm!',&
            OU_CLASS_ERROR,rtimestep%coutputModeError,&
            'tstep_createTimestepIndirect')
      end if
      call sys_halt()
    end select

    ! Copy data for adaptive time stepping scheme
    select case(rtimestep%iadapttimestep)
    case (TSTEP_NOADAPT)
      ! Create no adaptive time-stepping
      
    case (TSTEP_PIDADAPT)
      ! Create PID controller
      allocate(rtimestep%p_rpidController)
      rtimestep%p_rpidController = rtimestepTemplate%p_rpidController

    case (TSTEP_AUTOADAPT)
      ! Create automatic controller
      allocate(rtimestep%p_rautoController)
      rtimestep%p_rautoController = rtimestepTemplate%p_rautoController
      
      ! Allocate memory for temporal memory
      allocate(rtimestep%p_rautoController%RtempVectors(&
          size(rtimestepTemplate%p_rautoController%RtempVectors)))

    case (TSTEP_SERADAPT)
      ! Create SER controller
      allocate(rtimestep%p_rserController)
      rtimestep%p_rserController = rtimestepTemplate%p_rserController

    case default
      if (rtimestep%coutputModeError .gt. 0) then
        call output_line('Invalid type of adaptive time-stepping algorithm!',&
            OU_CLASS_ERROR,rtimestep%coutputModeError,&
            'tstep_createTimestepIndirect')
      end if
      call sys_halt()
    end select

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

    select case(tstep_getFamily(rtimestep%ctimestep))
    case (TSTEP_THETA)
      ! Release two-step theta scheme
      do i = lbound(rtimestep%p_rthetaScheme%RtempVectors, 1),&
             ubound(rtimestep%p_rthetaScheme%RtempVectors, 1)
        call lsysbl_releaseVector(rtimestep%p_rthetaScheme%RtempVectors(i))
      end do
      deallocate(rtimestep%p_rthetaScheme%RtempVectors)
      nullify(rtimestep%p_rthetaScheme%RtempVectors)
      deallocate(rtimestep%p_rthetaScheme)
      nullify(rtimestep%p_rthetaScheme)

    case (TSTEP_RUNGE_KUTTA)
      ! Release multi-step Runge-Kutta time-stepping scheme
      do i = lbound(rtimestep%p_rRungeKuttaScheme%RtempVectors, 1),&
             ubound(rtimestep%p_rRungeKuttaScheme%RtempVectors, 1)
        call lsysbl_releaseVector(rtimestep%p_rRungeKuttaScheme%RtempVectors(i))
      end do
      deallocate(rtimestep%p_rRungeKuttaScheme%RtempVectors)
      nullify(rtimestep%p_rRungeKuttaScheme%RtempVectors)

      if (iand(rtimestep%ctimestep,TSTEP_BUTCHER_TABLEAU) .eq. TSTEP_BUTCHER_TABLEAU) then
        ! Release memory for Butcher tableau
        deallocate(rtimestep%p_rRungeKuttaScheme%Da)
        deallocate(rtimestep%p_rRungeKuttaScheme%Db)
        deallocate(rtimestep%p_rRungeKuttaScheme%Dc)
        nullify(rtimestep%p_rRungeKuttaScheme%Da)
        nullify(rtimestep%p_rRungeKuttaScheme%Db)
        nullify(rtimestep%p_rRungeKuttaScheme%Dc)
      else
        deallocate(rtimestep%p_rRungeKuttaScheme%Dalpha)
        deallocate(rtimestep%p_rRungeKuttaScheme%Dbeta)
        nullify(rtimestep%p_rRungeKuttaScheme%Dalpha)
        nullify(rtimestep%p_rRungeKuttaScheme%Dbeta)
      end if
      
      deallocate(rtimestep%p_rRungeKuttaScheme)
      nullify(rtimestep%p_rRungeKuttaScheme)

    case default
      if (rtimestep%coutputModeError .gt. 0) then
        call output_line('Invalid type of time-stepping algorithm!',&
            OU_CLASS_ERROR,rtimestep%coutputModeError,&
            'tstep_releaseTimestep')
      end if
      call sys_halt()
    end select

    select case(rtimestep%iadapttimestep)
    case (TSTEP_NOADAPT)
      ! Nothing to release
      
    case (TSTEP_PIDADAPT)
      ! Release PID controller
      deallocate(rtimestep%p_rpidController)
      nullify(rtimestep%p_rpidController)

    case (TSTEP_AUTOADAPT)
      ! Release automatic controller
      do i = lbound(rtimestep%p_rautoController%RtempVectors, 1),&
          ubound(rtimestep%p_rautoController%RtempVectors, 1)
        call lsysbl_releaseVector(rtimestep%p_rautoController%RtempVectors(i))
      end do
      deallocate(rtimestep%p_rautoController%RtempVectors)
      nullify(rtimestep%p_rautoController%RtempVectors)
      deallocate(rtimestep%p_rautoController)
      nullify(rtimestep%p_rautoController)

    case (TSTEP_SERADAPT)
      ! Release SER controller
      deallocate(rtimestep%p_rserController)
      nullify(rtimestep%p_rserController)
      
    case default
      if (rtimestep%coutputModeError .gt. 0) then
        call output_line('Invalid type of adaptive time-stepping algorithm!',&
            OU_CLASS_ERROR,rtimestep%coutputModeError,&
            'tstep_releaseTimestep')
      end if
      call sys_halt()
    end select
    
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
    integer :: i,j


    ! Output general information
    call output_line('Timestep:')
    call output_line('---------')
    call output_line('Name:                          '//trim(rtimestep%sName))
    call output_line('Time stepping family:          '//trim(tstep_getFamilyName(rtimestep%ctimestep)))
    call output_line('Time stepping type:            '//trim(tstep_getTypeName(rtimestep%ctimestep)))
    call output_line('Number of time steps:          '//trim(sys_siL(rtimestep%nSteps,15)))
    call output_line('Number of rejected time steps: '//trim(sys_siL(rtimestep%nrejectedSteps,15)))

    ! Output detailed information
    if (present(bprintInternal)) then
      if (bprintInternal) then
        call output_line('ioutputLevel:                  '//trim(sys_siL(rtimestep%ioutputLevel,3)))
        call output_line('isolNorm:                      '//trim(sys_siL(rtimestep%isolNorm,3)))
        call output_line('iadaptTimestep:                '//trim(sys_siL(rtimestep%iadaptTimestep,3)))      
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
        call output_line('dStepPrevious:                 '//trim(sys_sdL(rtimestep%dstepPrevious,5)))

        ! Output information about the two-step theta scheme
        if (associated(rtimestep%p_rthetaScheme)) then
          call output_lbrk()
          call output_line('Two-step theta scheme:')
          call output_line('----------------------')
          call output_line('theta:                         '//&
              trim(sys_sdL(rtimestep%p_rthetaScheme%theta,5)))
          
          if (iand(rtimestep%ctimestep, TSTEP_FRACTIONAL_STEP) .eq. TSTEP_FRACTIONAL_STEP) then
            call output_line('alpha:                         '//&
                trim(sys_sdL(rtimestep%p_rthetaScheme%alpha,5)))
            call output_line('gamma:                         '//&
                trim(sys_sdL(rtimestep%p_rthetaScheme%gamma,5)))
          end if

          ! Output information about auxiliary vectors
          call output_lbrk()
          call output_line('Temporal vectors:')
          call output_line('-----------------')
          do i = lbound(rtimestep%p_rthetaScheme%RtempVectors, 1),&
                 ubound(rtimestep%p_rthetaScheme%RtempVectors, 1)
            call lsysbl_infoVector(rtimestep%p_rthetaScheme%RtempVectors(i))
          end do
        end if

        ! Output information about the multi-step runge-Kutta time-stepping scheme
        if (associated(rtimestep%p_rRungeKuttaScheme)) then
          call output_lbrk()
          call output_line('Multi-step Runge-Kutta scheme:')
          call output_line('------------------------------')
          call output_line('nstages:                       '//&
              trim(sys_siL(rtimestep%p_rRungeKuttaScheme%nstages,3)))         

          if (iand(rtimestep%ctimestep,TSTEP_BUTCHER_TABLEAU) .eq. TSTEP_BUTCHER_TABLEAU) then
            ! a-coefficient matrix
            call output_line('a-coefficient matrix:          '//&
                trim(sys_sdL(rtimestep%p_rRungeKuttaScheme%Da(1,1),5)), bnolinebreak=.true.)
            do j=2,rtimestep%p_rRungeKuttaScheme%nstages
              call output_line(trim(sys_sdP(rtimestep%p_rRungeKuttaScheme%Da(1,j),10,5)),&
                  bnolinebreak=.true.)
            end do
            do i=2,rtimestep%p_rRungeKuttaScheme%nstages
              call output_line('')
              call output_line('                               '//&
                trim(sys_sdL(rtimestep%p_rRungeKuttaScheme%Da(i,1),5)), bnolinebreak=.true.)
              do j=2,rtimestep%p_rRungeKuttaScheme%nstages
                call output_line(trim(sys_sdP(rtimestep%p_rRungeKuttaScheme%Da(i,j),10,5)),&
                    bnolinebreak=.true.)
              end do
            end do
            call output_line('')

            ! b-coefficient vector
            call output_line('b-coefficient vector:          '//&
                trim(sys_sdL(rtimestep%p_rRungeKuttaScheme%Db(1),5)), bnolinebreak=.true.)
            do i=2,rtimestep%p_rRungeKuttaScheme%nstages
              call output_line(trim(sys_sdP(rtimestep%p_rRungeKuttaScheme%Db(i),10,5)),&
                  bnolinebreak=.true.)
            end do
            call output_line('')

            ! c-coefficient vector
            call output_line('c-coefficient vector:          '//&
                trim(sys_sdL(rtimestep%p_rRungeKuttaScheme%Dc(1),5)), bnolinebreak=.true.)
            do i=2,rtimestep%p_rRungeKuttaScheme%nstages
              call output_line(trim(sys_sdP(rtimestep%p_rRungeKuttaScheme%Dc(i),10,5)),&
                  bnolinebreak=.true.)
            end do
            call output_line('')
            
          else
            ! alpha-coefficient matrix
            call output_line('alpha-coefficient matrix:      '//&
                trim(sys_sdL(rtimestep%p_rRungeKuttaScheme%Dalpha(1,1),5)), bnolinebreak=.true.)
            do j=2,rtimestep%p_rRungeKuttaScheme%nstages
              call output_line(trim(sys_sdP(rtimestep%p_rRungeKuttaScheme%Dalpha(1,j),10,5)),&
                  bnolinebreak=.true.)
            end do
            do i=2,rtimestep%p_rRungeKuttaScheme%nstages
              call output_line('')
              call output_line('                               '//&
                trim(sys_sdL(rtimestep%p_rRungeKuttaScheme%Dalpha(i,1),5)), bnolinebreak=.true.)
              do j=2,rtimestep%p_rRungeKuttaScheme%nstages
                call output_line(trim(sys_sdP(rtimestep%p_rRungeKuttaScheme%Dalpha(i,j),10,5)),&
                    bnolinebreak=.true.)
              end do
            end do
            call output_line('')
            
            ! beta-coefficient matrix
            call output_line('beta-coefficient matrix:       '//&
                trim(sys_sdL(rtimestep%p_rRungeKuttaScheme%Dbeta(1,1),5)), bnolinebreak=.true.)
            do j=2,rtimestep%p_rRungeKuttaScheme%nstages
              call output_line(trim(sys_sdP(rtimestep%p_rRungeKuttaScheme%Dbeta(1,j),10,5)),&
                  bnolinebreak=.true.)
            end do
            do i=2,rtimestep%p_rRungeKuttaScheme%nstages
              call output_line('')
              call output_line('                               '//&
                trim(sys_sdL(rtimestep%p_rRungeKuttaScheme%Dbeta(i,1),5)), bnolinebreak=.true.)
              do j=2,rtimestep%p_rRungeKuttaScheme%nstages
                call output_line(trim(sys_sdP(rtimestep%p_rRungeKuttaScheme%Dbeta(i,j),10,5)),&
                    bnolinebreak=.true.)
              end do
            end do
            call output_line('')

            ! c-coefficient vector
            call output_line('c-coefficient vector:          '//&
                trim(sys_sdL(rtimestep%p_rRungeKuttaScheme%Dc(1),5)), bnolinebreak=.true.)
            do i=2,rtimestep%p_rRungeKuttaScheme%nstages
              call output_line(trim(sys_sdP(rtimestep%p_rRungeKuttaScheme%Dc(i),10,5)),&
                  bnolinebreak=.true.)
            end do
            call output_line('')
            
          end if
          
          ! Output information about auxiliary vectors
          call output_lbrk()
          call output_line('Temporal vectors:')
          call output_line('-----------------')
          do i = lbound(rtimestep%p_rRungeKuttaScheme%RtempVectors, 1),&
                 ubound(rtimestep%p_rRungeKuttaScheme%RtempVectors, 1)
            call lsysbl_infoVector(rtimestep%p_rRungeKuttaScheme%RtempVectors(i))
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

          ! Output information about auxiliary vectors
          call output_lbrk()
          call output_line('Temporal vectors:')
          call output_line('-----------------')
          do i = lbound(rtimestep%p_rautoController%RtempVectors, 1),&
                 ubound(rtimestep%p_rautoController%RtempVectors, 1)
            call lsysbl_infoVector(rtimestep%p_rautoController%RtempVectors(i))
          end do
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
    rtimestep%dTime         = rtimeStep%dinitialTime
    rtimestep%dStep         = rtimeStep%dinitialStep
    rtimestep%dStepPrevious = rtimeStep%dinitialStep

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

  subroutine tstep_performTimestepSc(rproblemLevel, rtimestep,&
      rsolver, rsolution, fcb_nlsolverCallback, rcollection, nsteps, rsource)

!<description>   
    ! This subroutine performs nsteps time step with the
    ! time-stepping scheme given in rtimestep.
    ! This routine serves as wrapper for scalar vectors.
!</description>
    
!<input>
    ! Callback routines
    include 'intf_solvercallback.inc'

    ! Number of time steps to be performed
    integer, intent(in) :: nsteps
    
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
      call tstep_performTimestepBl(rproblemLevel, rtimestep,&
          rsolver, rsolutionBlock, fcb_nlsolverCallback, rcollection, nsteps,&
          rsourceBlock)
      
      ! Deallocate temporal 1-block vectors
      call lsysbl_releaseVector(rsolutionBlock)
      call lsysbl_releaseVector(rsourceBlock)
      
    else

      ! Convert scalar vector to 1-block vector
      call lsysbl_createVecFromScalar(rsolution, rsolutionBlock)
      
      ! Call block version of this subroutine
      call tstep_performTimestepBl(rproblemLevel, rtimestep,&
          rsolver, rsolutionBlock, fcb_nlsolverCallback, rcollection, nsteps)
      
      ! Deallocate temporal 1-block vector
      call lsysbl_releaseVector(rsolutionBlock)
      
    end if
    
  end subroutine tstep_performTimestepSc

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performTimestepBl(rproblemLevel, rtimestep,&
      rsolver, rsolution, fcb_nlsolverCallback, rcollection, nsteps, rsource)

!<description>   
    ! This subroutine performs nsteps time step with the
    ! time-stepping scheme given in rtimestep.
!</description>
    
!<input>
    ! Callback routines
    include 'intf_solvercallback.inc'

    ! Number of time steps to be performed
    integer, intent(in) :: nsteps
    
    ! OPTIONAL: constant right-hand side vector
    type(t_vectorBlock), intent(in), optional :: rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! Local variables
    type(t_solver), pointer :: p_rsolver
    type(t_vectorBlock), pointer :: p_rsolutionAux
    type(t_vectorBlock), pointer :: p_rsolutionRef
    type(t_vectorBlock), pointer :: p_rsolutionPrevious
    real(DP) :: dstep,dtime
    logical :: bcompatible,breject
    integer :: istep
    
    ! Set pointer to top-most nonlinear solver
    p_rsolver => solver_getNextSolverByTypes(rsolver, (/SV_NONLINEARMG, SV_NONLINEAR/))
    
    if (.not. associated(p_rsolver)) then
      if (rtimestep%coutputModeError .gt. 0) then
        call output_line('Unsupported/invalid solver type!',&
            OU_CLASS_ERROR,rtimestep%coutputModeError,&
            'tstep_performTimestepBl')
      end if
      call sys_halt()
    end if
    
    ! Set pointers to temporal vectors for the computation of an
    ! auxiliary solution by means of local sub-stepping
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
    
    ! What time-stepping family should be used?
    select case(tstep_getFamily(rtimestep%ctimestep))
      
    case (TSTEP_THETA)
      !-------------------------------------------------------------------------
      ! Two-level theta time-stepping scheme family
      !-------------------------------------------------------------------------
    
      ! Set pointers to temporal vectors
      p_rsolutionPrevious => rtimestep%p_rthetaScheme%RtempVectors(1)
      
      ! ... and check if vectors are compatible
      call lsysbl_isVectorCompatible(rsolution, p_rsolutionPrevious, bcompatible)
      if (.not.bcompatible) then
        call lsysbl_resizeVectorBlock(p_rsolutionPrevious, rsolution, .false.)
      end if

      ! Save the given solution vector to the temporal vector. If the
      ! computed time step is not accepted, then the backup of the
      ! given solution vector is used to recalculate the time step
      call lsysbl_copyVector(rsolution, p_rsolutionPrevious)

      ! Prepare time-stepping structure
      rtimestep%dscaleExplicit = (1.0_DP-rtimestep%p_rthetaScheme%theta)
      rtimestep%dscaleImplicit =         rtimestep%p_rthetaScheme%theta
      
      ! Perform nsteps steps with two-level theta scheme
      theta_loop: do istep = 1, nsteps

        ! Save current time
        dtime = rtimestep%dTime
        
        ! Adaptive time-stepping loop
        theta_loop_adapt: do
          
          ! Compute admissible time step size
          dstep = min(rtimestep%dStep, rtimestep%dfinalTime - dtime)
          
          ! Increase simulation time provisionally
          rtimestep%dStep  = dstep
          rtimestep%dTime  = dtime + rtimestep%dStep
          rtimestep%nSteps = rtimestep%nSteps + 1
          
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
              rsolution, p_rsolutionPrevious, fcb_nlsolverCallback,&
              rcollection, rsource)
          
          ! Adjust status information of top-most solver
          call solver_copySolver(p_rsolver, rsolver, .false., .true.)
          
          ! Compute reference solution by performing two time steps of size Dt/2
          if (rtimestep%iadaptTimestep .eq. TSTEP_AUTOADAPT) then

            ! Perform two halve time-steps
            rtimestep%dStep = dstep/2.0_DP

            ! Output information
            if (rtimestep%coutputModeInfo .gt. 0) then
              call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_line('Two-level theta-scheme [first substep], Time = '//&
                               trim(sys_sdEL(rtimestep%dTime,5))//&
                               ' Stepsize = '//trim(sys_sdEL(rtimestep%dStep,5)),&
                               OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
            end if

            ! Solve the nonlinear algebraic system for time step t^n -> t^{n+1/2}
            call nlsol_solveMultigrid(rproblemLevel, rtimestep, p_rsolver,&
                p_rsolutionRef, p_rsolutionPrevious, fcb_nlsolverCallback,&
                rcollection, rsource)
            
            ! Save intermediate solution
            call lsysbl_copyVector(p_rsolutionRef, p_rsolutionAux)

            ! Output information
            if (rtimestep%coutputModeInfo .gt. 0) then
              call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_line('Two-level theta-scheme [second substep], Time = '//&
                               trim(sys_sdEL(rtimestep%dTime,5))//&
                               ' Stepsize = '//trim(sys_sdEL(rtimestep%dStep,5)),&
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
                p_rsolutionRef, p_rsolutionPrevious)
            
            ! Reset the time step size
            rtimestep%dStep = dstep
          else
            
            ! Check if solution from this time step can be accepted and
            ! adjust the time step size automatically if this is required
            breject = tstep_checkTimestep(rtimestep, p_rsolver,&
                rsolution, p_rsolutionPrevious)
          end if
            
          ! Output information
          if (rtimestep%coutputModeInfo .gt. 0) then
            call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
            call output_separator(OU_SEP_TILDE,OU_CLASS_MSG,rtimestep%coutputModeInfo)
            call output_line('Time step was     '//merge('!!! rejected !!!','accepted        ',breject),&
                OU_CLASS_MSG,rtimestep%coutputModeInfo)
            call output_line('New stepsize:     '//trim(sys_sdEL(rtimestep%dStep,5)),&
                OU_CLASS_MSG,rtimestep%coutputModeInfo)
            call output_line('Last stepsize:    '//trim(sys_sdEL(rtimestep%dStepPrevious,5)),&
                OU_CLASS_MSG,rtimestep%coutputModeInfo)
            call output_line('Relative changes: '//trim(sys_sdEL(rtimestep%drelChange,5)),&
                OU_CLASS_MSG,rtimestep%coutputModeInfo)
            call output_separator(OU_SEP_TILDE,OU_CLASS_MSG,rtimestep%coutputModeInfo)
            call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
          end if

          ! Do we have to reject to current solution?
          if (breject) then
            ! Yes, so restore the old solution and repeat the time step
            call lsysbl_copyVector(p_rsolutionPrevious, rsolution)
            if (rtimestep%coutputModeWarning .gt. 0) then
              call output_line('Time step was rejected!',&
                  OU_CLASS_WARNING,rtimestep%coutputModeWarning,&
                  'tstep_performTimestepBl')
            end if
          else
            ! No, accept current solution and
            ! exit adaptive time-stepping loop
            exit theta_loop_adapt
          end if
        end do theta_loop_adapt
        
      end do theta_loop

    case (TSTEP_RUNGE_KUTTA)
      !-------------------------------------------------------------------------
      ! Multi-stage Runge-Kutta time-stepping scheme family
      !-------------------------------------------------------------------------
      
      ! Set pointers to temporal vectors
      p_rsolutionPrevious => rtimestep%p_rRungeKuttaScheme%RtempVectors(1)
      
      ! Save the given solution vector to the temporal vector. If the
      ! computed time step is not accepted, then the backup of the
      ! given solution vector is used to recalculate the time step
      call lsysbl_copyVector(rsolution, p_rsolutionPrevious)

      select case(rtimestep%ctimestep)

      case (TSTEP_SSPERK_1_1)
        
        ! Perform nsteps steps with explicit SSP-RK(1,1) scheme
        ssperk_1_1_loop: do istep = 1, nsteps

          ! Save current time
          dtime = rtimestep%dTime
          
          ! Adaptive time-stepping loop
          ssperk_1_1_loop_adapt: do

            ! Compute admissible time step size
            dstep = min(rtimestep%dStep, rtimestep%dfinalTime - dtime)

            ! Increase simulation time provisionally
            rtimestep%dStep  = rtimestep%p_rRungeKuttaScheme%Dc(1)*dstep
            rtimestep%dTime  = dtime + rtimestep%dStep
            rtimestep%nSteps = rtimestep%nSteps + 1
            
            ! Prepare time-stepping structure
            rtimestep%dscaleExplicit = rtimestep%p_rRungeKuttaScheme%Dbeta(1,1)
            rtimestep%dscaleImplicit = 0.0_DP

            ! Output information
            if (rtimestep%coutputModeInfo .gt. 0) then
              call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_line('Explicit SSP-RK(1,1) scheme, Time = '//&
                  trim(sys_sdEL(rtimestep%dTime,5))//&
                  ' Stepsize = '//trim(sys_sdEL(rtimestep%dStep,5)),&
                  OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
            end if

            ! Solve the nonlinear algebraic system in the time interval (t^n, t^{n+1})
            call nlsol_solveMultigrid(rproblemLevel, rtimestep, p_rsolver,&
                rsolution, p_rsolutionPrevious, fcb_nlsolverCallback,&
                rcollection, rsource)

            ! Adjust status information of top-most solver
            call solver_copySolver(p_rsolver, rsolver, .false., .true.)

            ! Reset the time step size
            rtimestep%dStep = dstep
            
            ! Check if solution from this time step can be accepted and
            ! adjust the time step size automatically if this is required
            breject = tstep_checkTimestep(rtimestep, p_rsolver,&
                rsolution, p_rsolutionPrevious)
            
            ! Output information
            if (rtimestep%coutputModeInfo .gt. 0) then
              call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_separator(OU_SEP_TILDE,OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_line('Time step was     '//merge('!!! rejected !!!','accepted        ',breject),&
                  OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_line('New stepsize:     '//trim(sys_sdEL(rtimestep%dStep,5)),&
                  OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_line('Last stepsize:    '//trim(sys_sdEL(rtimestep%dStepPrevious,5)),&
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
              call lsysbl_copyVector(p_rsolutionPrevious, rsolution)
              if (rtimestep%coutputModeWarning .gt. 0) then
                call output_line('Time step was rejected!',&
                    OU_CLASS_WARNING,rtimestep%coutputModeWarning,&
                    'tstep_performTimestepBl')
              end if
            else
              ! No, accept current solution and
              ! exit adaptive time-stepping loop
              exit ssperk_1_1_loop_adapt
            end if
          end do ssperk_1_1_loop_adapt
          
        end do ssperk_1_1_loop
        
      case (TSTEP_SSPERK_2_2)
        
        ! Perform nsteps steps with explicit SSP-RK(2,2) scheme
        ssperk_2_2_loop: do istep = 1, nsteps

          ! Save current time
          dtime = rtimestep%dTime
          
          ! Adaptive time-stepping loop
          ssperk_2_2_loop_adapt: do

            ! Compute admissible time step size
            dstep = min(rtimestep%dStep, rtimestep%dfinalTime - dtime)

            ! -- first step --
            
            ! Increase simulation time provisionally
            rtimestep%dStep  = rtimestep%p_rRungeKuttaScheme%Dc(1)*dstep
            rtimestep%dTime  = dtime + rtimestep%dStep
            rtimestep%nSteps = rtimestep%nSteps + 1
            
            ! Prepare time-stepping structure
            rtimestep%dscaleExplicit = rtimestep%p_rRungeKuttaScheme%Dbeta(1,1)
            rtimestep%dscaleImplicit = 0.0_DP

            ! Output information
            if (rtimestep%coutputModeInfo .gt. 0) then
              call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_line('Explicit SSP-RK(2,2) scheme [first step], Time = '//&
                  trim(sys_sdEL(rtimestep%dTime,5))//&
                  ' Stepsize = '//trim(sys_sdEL(rtimestep%dStep,5)),&
                  OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
            end if

            ! Solve the nonlinear algebraic system in the time interval (t^n, t^{n+1})
            call nlsol_solveMultigrid(rproblemLevel, rtimestep, p_rsolver,&
                rsolution, p_rsolutionPrevious, fcb_nlsolverCallback,&
                rcollection, rsource)

            ! -- second step --
            
            ! Increase simulation time provisionally
            rtimestep%dStep  = rtimestep%p_rRungeKuttaScheme%Dc(2)*dstep
            rtimestep%dTime  = dtime + rtimestep%dStep

            ! Prepare time-stepping structure
            rtimestep%dscaleExplicit = rtimestep%p_rRungeKuttaScheme%Dbeta(2,1)
            rtimestep%dscaleImplicit = 0.0_DP

            ! Output information
            if (rtimestep%coutputModeInfo .gt. 0) then
              call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_line('Explicit SSP-RK(2,2) scheme [second step], Time = '//&
                  trim(sys_sdEL(rtimestep%dTime,5))//&
                  ' Stepsize = '//trim(sys_sdEL(rtimestep%dStep,5)),&
                  OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_separator(OU_SEP_AT,OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
            end if

            ! Solve the nonlinear algebraic system in the time interval (t^n, t^{n+1})
            call nlsol_solveMultigrid(rproblemLevel, rtimestep, p_rsolver,&
                rsolution, p_rsolutionPrevious, fcb_nlsolverCallback,&
                rcollection, rsource)

            ! Solve the nonlinear algebraic system in the time interval (t^n, t^{n+1})
            call nlsol_solveMultigrid(rproblemLevel, rtimestep, p_rsolver,&
                rsolution, p_rsolutionPrevious, fcb_nlsolverCallback,&
                rcollection, rsource)
            
            ! Adjust status information of top-most solver
            call solver_copySolver(p_rsolver, rsolver, .false., .true.)

            ! Reset the time step size
            rtimestep%dStep = dstep
            
            ! Check if solution from this time step can be accepted and
            ! adjust the time step size automatically if this is required
            breject = tstep_checkTimestep(rtimestep, p_rsolver,&
                rsolution, p_rsolutionPrevious)
            
            ! Output information
            if (rtimestep%coutputModeInfo .gt. 0) then
              call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_separator(OU_SEP_TILDE,OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_line('Time step was     '//merge('!!! rejected !!!','accepted        ',breject),&
                  OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_line('New stepsize:     '//trim(sys_sdEL(rtimestep%dStep,5)),&
                  OU_CLASS_MSG,rtimestep%coutputModeInfo)
              call output_line('Last stepsize:    '//trim(sys_sdEL(rtimestep%dStepPrevious,5)),&
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
              call lsysbl_copyVector(p_rsolutionPrevious, rsolution)
              if (rtimestep%coutputModeWarning .gt. 0) then
                call output_line('Time step was rejected!',&
                    OU_CLASS_WARNING,rtimestep%coutputModeWarning,&
                    'tstep_performTimestepBl')
              end if
            else
              ! No, accept current solution and
              ! exit adaptive time-stepping loop
              exit ssperk_2_2_loop_adapt
            end if
          end do ssperk_2_2_loop_adapt
          
        end do ssperk_2_2_loop
        
      case (TSTEP_SSPERK_3_3)

      end select
      
    case default
      call output_line('Unsupported time-stepping scheme!',&
          OU_CLASS_ERROR,rtimestep%coutputModeWarning,&
          'tstep_performTiemstepBl')
      call sys_halt()
    end select
    
  contains

    ! Here, the real working routines follow
    
    !*************************************************************
    ! This function checks steady-state convergence
    
    logical function checkConvergence(rtimestep)
      type(t_timestep), intent(in) :: rtimestep
      
      ! Reached final time, then exit infinite time loop?
      checkConvergence = (rtimestep%dTime .ge. rtimestep%dfinalTime)
      
      ! Reached steady state limit?
      if (rtimestep%depsSteady > 0.0_DP) then        
        checkConvergence = checkConvergence .or.&
            ((rsolver%dfinalDefect   < rsolver%dinitialDefect) .and.&
             (rsolver%dinitialDefect < rtimestep%dStep*rtimestep%depsSteady))
      end if
    end function checkConvergence
    
  end subroutine tstep_performTimestepBl

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performTimestepScCpl(rproblemLevel, rtimestep,&
      rsolver, Rsolution, fcb_nlsolverCallback, rcollection, nsteps, Rsource)

!<description>
    ! This subroutine performs one single time step with the
    ! time-stepping scheme given in rtimestep.
    ! This routine serves as wrapper for scalar vectors.
!</description>

!<input>
    ! Callback routines
    include 'intf_solvercallback.inc'

    ! Number of time steps to be performed
    integer, intent(in) :: nsteps
    
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
        call lsysbl_createVecFromScalar(Rsource(icomponent),&
            RsourceBlock(icomponent))
      end do
      
      ! Call block version of this subroutine
      call tstep_performTimestepBlCpl(rproblemLevel, rtimestep,&
          rsolver, RsolutionBlock, fcb_nlsolverCallback, rcollection, nsteps,&
          RsourceBlock)
      
      ! Deallocate temporal 1-block vectors
      do icomponent = 1, size(Rsolution)
        call lsysbl_releaseVector(RsolutionBlock(icomponent))
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
      call tstep_performTimestepBlCpl(rproblemLevel, rtimestep,&
          rsolver, RsolutionBlock, fcb_nlsolverCallback, rcollection, nsteps)
      
      ! Deallocate temporal 1-block vectors
      do icomponent = 1, size(Rsolution)
        call lsysbl_releaseVector(RsolutionBlock(icomponent))
      end do
      
      ! Release temporal memory
      deallocate(rsolutionBlock)
      
    end if
    
  end subroutine tstep_performTimestepScCpl

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performTimestepBlCpl(rproblemLevel, rtimestep,&
      rsolver, Rsolution, fcb_nlsolverCallback, rcollection, nsteps, Rsource)

!<description>
    ! This subroutine performs one single time step with the
    ! time-stepping scheme given in rtimestep.
!</description>

!<input>
    ! Callback routines
    include 'intf_solvercallback.inc'

    ! Number of time steps to be performed
    integer, intent(in) :: nsteps
    
    ! OPTIONAL: constant right-hand side vector
    type(t_vectorBlock), dimension(:), intent(in), optional :: Rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! solution vector
    type(t_vectorBlock), dimension(:), intent(inout) :: Rsolution

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_solver), pointer :: p_rsolver
    type(t_vectorBlock), dimension(:), pointer :: p_RsolutionAux
    type(t_vectorBlock), dimension(:), pointer :: p_RsolutionRef
    type(t_vectorBlock), dimension(:), pointer :: p_RsolutionPrevious
    logical :: bcompatible
    integer :: icomponent,ncomponent
    
    ! Get number of components
    ncomponent = size(Rsolution)

    ! Check size of multi-component source vector (if any)
    if (present(Rsource)) then
      if (size(Rsource) .ne. ncomponent) then
        if (rtimestep%coutputModeError .gt. 0) then
          call output_line('Dimension of coupled problem mismatch!',&
              OU_CLASS_ERROR,rtimestep%coutputModeError,&
              'tstep_performTimestepBlCpl')
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
            'tstep_performTimestepBlCpl')
      end if
      call sys_halt()
    end if

    ! Set pointer to temporal vectors
    p_RsolutionPrevious => rtimestep%p_rthetaScheme%RtempVectors(1:ncomponent)
    
    ! ... and check if vectors are compatible
    do icomponent = 1, ncomponent
      call lsysbl_isVectorCompatible(Rsolution(icomponent),&
          p_RsolutionPrevious(icomponent), bcompatible)
      if (.not.bcompatible) then
        call lsysbl_resizeVectorBlock(p_RsolutionPrevious(icomponent),&
            Rsolution(icomponent), .false.)
      end if
      
      ! Save the given solution vector to the temporal vector. If the
      ! computed time step is not accepted, then the backup of the
      ! given solution vector is used to recalculate the time step
      call lsysbl_copyVector(Rsolution(icomponent),&
          p_RsolutionPrevious(icomponent))
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

    ! What time-stepping family should be used?
    select case(tstep_getFamily(rtimestep%ctimestep))
      
    case (TSTEP_THETA)
      !-------------------------------------------------------------------------
      ! Two-level theta time-stepping scheme
      !-------------------------------------------------------------------------
      
    case (TSTEP_RUNGE_KUTTA)
      !-------------------------------------------------------------------------
      ! Multi-stage Runge-Kutta time-stepping scheme family
      !-------------------------------------------------------------------------
      
    case default
      call output_line('Unsupported time-stepping scheme!',&
          OU_CLASS_ERROR,rtimestep%coutputModeWarning,&
          'tstep_performTimestepBlCpl')
      call sys_halt()
    end select
    
  end subroutine tstep_performTimestepBlCpl
  
  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performThetaStepCpl(rproblemLevel, rtimestep,&
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
    ! Runge-Kutta time-stepping scheme instead.
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
    type(t_vectorBlock), dimension(:), pointer :: p_RsolutionPrevious
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
              'tstep_performThetaStepCpl')
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
            'tstep_performThetaStepCpl')
      end if
      call sys_halt()
    end if

    ! Set pointer to temporal vectors
    p_RsolutionPrevious => rtimestep%p_rthetaScheme%RtempVectors(1:ncomponent)

    ! ... and check if vectors are compatible
    do icomponent = 1, ncomponent
      call lsysbl_isVectorCompatible(Rsolution(icomponent),&
          p_RsolutionPrevious(icomponent), bcompatible)
      if (.not.bcompatible) then
        call lsysbl_resizeVectorBlock(p_RsolutionPrevious(icomponent),&
            Rsolution(icomponent), .false.)
      end if

      ! Save the given solution vector to the temporal vector. If the
      ! computed time step is not accepted, then the backup of the
      ! given solution vector is used to recalculate the time step
      call lsysbl_copyVector(Rsolution(icomponent),&
          p_RsolutionPrevious(icomponent))
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
          Rsolution, p_RsolutionPrevious, fcb_nlsolverCallback,&
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
            p_RsolutionRef, p_RsolutionPrevious, fcb_nlsolverCallback,&
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
            p_RsolutionRef, p_RsolutionPrevious)

        ! Prepare time step size for next "large" time step
        rtimestep%dStep = rtimestep%dStep*2.0_DP

      else

        ! Check if solution from this time step can be accepted and
        ! adjust the time step size automatically if this is required
        breject = tstep_checkTimestep(rtimestep, p_rsolver,&
            Rsolution, p_RsolutionPrevious)

      end if

      if (rtimestep%coutputModeInfo .gt. 0) then
        call output_lbrk(OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_separator(OU_SEP_TILDE,OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('Time step was     '//merge('!!! rejected !!!','accepted        ',breject),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('New stepsize:     '//trim(sys_sdEL(rtimestep%dStep,5)),&
                         OU_CLASS_MSG,rtimestep%coutputModeInfo)
        call output_line('Last stepsize:    '//trim(sys_sdEL(rtimestep%dStepPrevious,5)),&
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
          call lsysbl_copyVector(p_RsolutionPrevious(icomponent),&
              Rsolution(icomponent))
        end do
      else
        ! No, accept current solution and
        ! exit adaptive time-stepping loop
        exit timeadapt
      end if
    end do timeadapt

  end subroutine tstep_performThetaStepCpl

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
        rtimestep%dStepPrevious = rtimestep%dStep
        rtimestep%dStep = max(rtimestep%dminStep, min(rtimestep%dmaxStep, dStepOpt))

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
          dStepOpt = max(p_rpidController%dDecreaseFactor*rtimestep%dStepPrevious,&
                         min(p_rpidController%dIncreaseFactor*rtimestep%dStepPrevious,&
                             dStepOpt))
        else

          ! Adopt value from previous time step
          dstepOpt = rtimestep%dStep
        end if

        ! Impose upper/lower bounds on absolute values
        rtimestep%dStepPrevious = rtimestep%dStep
        rtimestep%dStep = max(rtimestep%dminStep, min(rtimestep%dmaxStep, dStepOpt))

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
      rtimestep%dStepPrevious = rtimestep%dStep

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
        rtimestep%dTime         = rtimestep1%dTime
        rtimestep%dStep         = rtimestep1%dStep
        rtimestep%dStepPrevious = rtimestep1%dStepPrevious
        rtimestep%drelChange    = rtimestep1%drelChange
      elseif (rtimestep1%dStep .lt. rtimestep%dStep) then
        rtimestep%dTime         = rtimestep1%dTime
        rtimestep%dStep         = rtimestep1%dStep
        rtimestep%dStepPrevious = rtimestep1%dStepPrevious
        rtimestep%drelChange    = rtimestep1%drelChange
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
