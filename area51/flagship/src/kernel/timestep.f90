!##############################################################################
!# ****************************************************************************
!# <name> timestep </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic data structures and subroutines for handling
!# time-dependent problems. A calling program should create a new object of
!# time t_timestep and fill the required parameters either by reading a user-
!# supplied parameter file of via direct adjustment.
!#
!# The following routines are available:
!#
!# 1.) tstep_performThetaStep = tstep_performThetaStepScalar /
!#                              tstep_performThetaStepBlock
!#     -> Performs one step by the two-level theta scheme to compute the
!#        solution for the time interval dTime..dTime+dStep. 
!#
!# 2.) tstep_performRKStep = tstep_performRKStepScalar /
!#                           tstep_performRKStepBlock
!#     -> Performs one step of an explicit Runge-Kutta scheme to compute the
!#        solution for the time interval dTime..dTime+dStep
!#
!# 3.) tstep_performPseudoStepping = tstep_performPseudoStepScalar /
!#                                   tstep_performPseudoStepBlock
!#     -> Performs pseudo time-stepping to compute the steady state solution
!#
!#
!# The following auxiliary routines are available: 
!#
!# 1.) tstep_checkTimestep
!#     -> Checks the solution computed in one time step and adjust
!#        the size of the time step accordingly
!#
!# </purpose>
!##############################################################################

module timestep

  use collection
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use problem
  use solver
  use solverlinear
  use solvernonlinear

  implicit none

  private
  public :: tstep_performThetaStep
  public :: tstep_performRKStep
  public :: tstep_performPseudoStepping

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  interface tstep_performThetaStep
    module procedure tstep_performThetaStepScalar
    module procedure tstep_performThetaStepBlock
  end interface

  interface tstep_performRKStep
    module procedure tstep_performRKStepScalar
    module procedure tstep_performRKStepBlock
  end interface

  interface tstep_performPseudoStepping
    module procedure tstep_performPseudoSteppingScalar
    module procedure tstep_performPseudoSteppingBlock
  end interface

contains

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  !<subroutine>
  
  subroutine tstep_performThetaStepScalar(rproblemLevel, rtimestep, rsolver, rsolution,&
                                          fcb_calcResidual, fcb_calcJacobian,&
                                          fcb_applyJacobian, fcb_setBoundary,&
                                          rcollection, rrhs)

!<description>
    ! This subroutine performs one step by means of the two-level theta-scheme
    ! to advance the solution from time level $t^n$ to time level $t^{n+1}$.
    ! Note that this subroutine serves as wrapper for scalar solution vectors only.
!</description>

!<input>
    ! Callback routines
    include 'intf_nlsolcallback.inc'

    ! OPTIONAL: right-hand side vector
    type(t_vectorScalar), intent(IN), optional :: rrhs
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! solution vector
    type(t_vectorScalar), intent(INOUT) :: rsolution

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rsolutionBlock
    type(t_vectorBlock) :: rrhsBlock

    if (present(rrhs)) then

      call lsysbl_createVecFromScalar(rsolution, rsolutionBlock)
      call lsysbl_createVecFromScalar(rrhs, rrhsBlock)
      call tstep_performThetaStepBlock(rproblemLevel, rtimestep, rsolver,&
                                       rsolutionBlock,&
                                       fcb_calcResidual, fcb_calcJacobian,&
                                       fcb_applyJacobian, fcb_setBoundary,&
                                       rcollection, rrhsBlock)
      call lsysbl_releaseVector(rsolutionBlock)
      call lsysbl_releaseVector(rrhsBlock)

    else
      
      call lsysbl_createVecFromScalar(rsolution, rsolutionBlock)
      call tstep_performThetaStepBlock(rproblemLevel, rtimestep, rsolver,&
                                       rsolutionBlock,&
                                       fcb_calcResidual, fcb_calcJacobian,&
                                       fcb_applyJacobian, fcb_setBoundary,&
                                       rcollection)
      call lsysbl_releaseVector(rsolutionBlock)

    end if
  end subroutine tstep_performThetaStepScalar

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performThetaStepBlock(rproblemLevel, rtimestep, rsolver, rsolution,&
                                         fcb_calcResidual, fcb_calcJacobian,&
                                         fcb_applyJacobian, fcb_setBoundary,&
                                         rcollection, rrhs)

!<description>
    ! This subroutine performs one step by means of the two-level theta-scheme
    ! to advance the solution from time level $t^n$ to time level $t^{n+1}$.
    ! This routine can handle the forward Euler (theta=0), backward Euler (theta=1)
    ! and the Crank-Nicolson (theta=0.5) time stepping algorithm. However, the
    ! fully explicit scheme is known to be unstable. It is therefore highly 
    ! recommended to use an explicit Runge-Kutta scheme instead.
!</description>

!<input>
    ! Callback routines
    include 'intf_nlsolcallback.inc'

    ! OPTIONAL: right-hand side vector
    type(t_vectorBlock), intent(IN), optional :: rrhs
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep

    ! solver structure
    type(t_solver), intent(INOUT), target :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_solver), pointer :: p_rsolver
    type(t_vectorBlock), pointer :: p_rsolutionOld, p_rsolutionRef, p_rsolutionAux
    logical :: bcompatible,breject
    
    ! Set pointer to solver structure
    p_rsolver => rsolver
    
    ! Walk down the solver structure until applicable solver is reached
    do while (associated(p_rsolver))
      if ((p_rsolver%csolverType .eq. SV_NONLINEARMG) .or. &
          (p_rsolver%csolverType .eq. SV_NONLINEAR  )) exit
      p_rsolver => solver_getNextSolver(p_rsolver)
    end do

    if (.not. associated(p_rsolver)) then
      call output_line('Unsupported/invalid solver type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'tstep_performThetaStepBlock')
      call sys_halt()
    end if

    ! Set pointers to temporal vector
    p_rsolutionOld => rtimestep%RtempVectors(1)

    ! Check if vectors are compatible
    call lsysbl_isVectorCompatible(rsolution, p_rsolutionOld, bcompatible)
    if (.not.bcompatible)&
        call lsysbl_resizeVectorBlock(p_rsolutionOld, rsolution, .false.)

    ! Save the given solution vector
    call lsysbl_copyVector(rsolution, p_rsolutionOld)

    
    ! Set pointer to second temporal vector for the solution computed by local substepping
    if (rtimestep%iadaptTimestep .eq. SV_TIMESTEP_AUTOADAPT) then
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
      rtimestep%dTime = rtimestep%dTime+rtimestep%dStep
      rtimestep%nSteps= rtimestep%nSteps+1
      
      if (rtimestep%ioutputLevel .ge. SV_IOLEVEL_INFO) then
        call output_lbrk()
        call output_separator(OU_SEP_AT)
        call output_line('Two-level theta-scheme, Time = '//trim(sys_sdEL(rtimestep%dTime,5))//&
                         ' Stepsize = '//trim(sys_sdEL(rtimestep%dStep,5)))
        call output_separator(OU_SEP_AT)
        call output_lbrk()
      end if


      ! Solve the nonlinear algebraic system for time step t^n -> t^{n+1}
      call nlsol_solveMultigrid(rproblemLevel, rtimestep, p_rsolver,&
                                rsolution, p_rsolutionOld,&
                                fcb_calcResidual, fcb_calcJacobian,&
                                fcb_applyJacobian, fcb_setBoundary,&
                                rcollection, rrhs)

      ! Adjust status information of top-most solver
      call solver_copySolver(p_rsolver, rsolver, .false., .true.)


      ! Compute reference solution by performing two time steps of size Dt/2
      if (rtimestep%iadaptTimestep .eq. SV_TIMESTEP_AUTOADAPT) then

        ! Set time step to smaller value
        rtimestep%dStep = rtimestep%dStep/2.0_DP

        if (rtimestep%ioutputLevel .ge. SV_IOLEVEL_INFO) then
          call output_lbrk()
          call output_separator(OU_SEP_AT)
          call output_line('First substep in automatic time step control')
          call output_separator(OU_SEP_AT)
          call output_lbrk()
        end if


        ! Solve the nonlinear algebraic system for time step t^n -> t^{n+1/2}
        call nlsol_solveMultigrid(rproblemLevel, rtimestep, p_rsolver,&
                                  p_rsolutionRef, p_rsolutionOld,&
                                  fcb_calcResidual, fcb_calcJacobian,&
                                  fcb_applyJacobian, fcb_setBoundary,&
                                  rcollection, rrhs)

        ! Save intermediate solution
        call lsysbl_copyVector(p_rsolutionRef, p_rsolutionAux)

        if (rtimestep%ioutputLevel .ge. SV_IOLEVEL_INFO) then
          call output_lbrk()
          call output_separator(OU_SEP_AT)
          call output_line('Second substep in automatic time step control')
          call output_separator(OU_SEP_AT)
          call output_lbrk()
        end if

        ! Solve the nonlinear algebraic system for time step t^{n+1/2} -> t^{n+1}
        call nlsol_solveMultigrid(rproblemLevel, rtimestep, p_rsolver,&
                                  p_rsolutionRef, p_rsolutionAux,&
                                  fcb_calcResidual, fcb_calcJacobian,&
                                  fcb_applyJacobian, fcb_setBoundary,&
                                  rcollection, rrhs)

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
      
      if (rtimestep%ioutputlevel .ge. SV_IOLEVEL_VERBOSE) then
        call output_lbrk()
        call output_separator(OU_SEP_TILDE)
        call output_line('Time step was     '//merge('!!! rejected !!!','accepted        ',breject))
        call output_line('New stepsize:     '//trim(sys_sdEL(rtimestep%dStep,5)))
        call output_line('Last stepsize:    '//trim(sys_sdEL(rtimestep%dStep1,5)))
        call output_line('Relative changes: '//trim(sys_sdEL(rtimestep%drelChange,5)))
        call output_separator(OU_SEP_TILDE)
        call output_lbrk()
      end if

      ! Write time step to file?
      if (rtimestep%ioutputLevel .ge. SV_IOLEVEL_FILE) then
        write(rtimestep%iunitLogfile, FMT='("(01),",F20.5,3(",",E16.8E3))')&
            rtimestep%dTime, rtimestep%dStep, rtimestep%dStep/rtimestep%dStep1,&
            rtimestep%drelChange
      end if

      
      ! Do we have to reject to current solution?
      if (breject) then
        ! Yes, so restore the old solution and
        ! repeat the adaptive time-stepping loop
        call lsysbl_copyVector(p_rsolutionOld, rsolution)
      else
        ! No, accept current solution and 
        ! exit adaptive time-stepping loop
        exit timeadapt
      end if
    end do timeadapt
  end subroutine tstep_performThetaStepBlock

  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performRKStepScalar(rproblemLevel, rtimestep, rsolver, rsolution,&
                                       fcb_calcRHS, fcb_setBoundary, rcollection, rrhs)

!<description>
    ! This subroutine performs one step by means of an explicit Runge-Kutta scheme
    ! to advance the solution from time level $t^n$ to time level $t^{n+1}$.
    ! Note that this subroutine serves as wrapper for scalar solution vectors only.
!</description>

!<input>
    ! Callback routines
    include 'intf_nlsolcallback.inc'

    ! OPTIONAL: right-hand side vector
    type(t_vectorScalar), intent(IN), optional :: rrhs
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! solution vector
    type(t_vectorScalar), intent(INOUT) :: rsolution

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rsolutionBlock
    type(t_vectorBlock) :: rrhsBlock

    if (present(rrhs)) then

      call lsysbl_createVecFromScalar(rsolution, rsolutionBlock)
      call lsysbl_createVecFromScalar(rrhs, rrhsBlock)
      call tstep_performRKStepBlock(rproblemLevel, rtimestep, rsolver, rsolutionBlock,&
                                    fcb_calcRHS, fcb_setBoundary, rcollection, rrhsBlock)
      call lsysbl_releaseVector(rsolutionBlock)
      call lsysbl_releaseVector(rrhsBlock)

    else
      
      call lsysbl_createVecFromScalar(rsolution, rsolutionBlock)
      call tstep_performRKStepBlock(rproblemLevel, rtimestep, rsolver, rsolutionBlock,&
                                    fcb_calcRHS, fcb_setBoundary, rcollection)
      call lsysbl_releaseVector(rsolutionBlock)

    end if
  end subroutine tstep_performRKStepScalar
  
  ! *****************************************************************************

!<subroutine>

  subroutine tstep_performRKStepBlock(rproblemLevel, rtimestep, rsolver, rsolution, &
                                      fcb_calcRHS, fcb_setBoundary, rcollection, rrhs)

!<description>
    ! This subroutine performs one step by means of an explicit Runge-Kutta scheme
    ! to advance the solution from time level $t^n$ to time level $t^{n+1}$.
    ! In particular, the two step version suggested by Richtmyer and Morton is
    ! implemented so that no artificial diffusion needs to be constructed.
!</description>

!<input>
    ! Callback routines
    include 'intf_nlsolcallback.inc'

    ! OPTIONAL: right-hand side vector
    type(t_vectorBlock), intent(IN), optional :: rrhs
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(INOUT), target :: rtimestep

    ! solver structure
    type(t_solver), intent(INOUT), target :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>


    ! local variables
    type(t_solver), pointer :: p_rsolver
    type(t_vectorBlock), pointer :: p_rsolutionOld,p_rsolutionRef,p_rsolutionAux
    type(t_vectorBlock), pointer :: p_rrhs
    type(t_vectorBlock), pointer :: p_raux
    integer :: istep
    logical :: bcompatible,breject
    

    ! Set pointer
    p_rsolver => rsolver

    ! Walk down the solver structure until applicable solver is reached
    do while (associated(p_rsolver))
      if ((p_rsolver%csolverType .eq. SV_LINEARMG) .or. &
          (p_rsolver%csolverType .eq. SV_LINEAR  )) exit
      p_rsolver => solver_getNextSolver(p_rsolver)
    end do

    if (.not. associated(p_rsolver)) then
      call output_line('Unsupported/invalid solver type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'tstep_performRKStepBlock')
      call sys_halt()
    end if
    
    ! Set pointers to temporal vector
    p_rsolutionOld => rtimestep%RtempVectors(1)
    p_rrhs         => rtimestep%RtempVectors(2)
    p_raux         => rtimestep%RtempVectors(3)
    
    ! Check if vectors are compatible
    call lsysbl_isVectorCompatible(rsolution, p_rsolutionOld, bcompatible)
    if (.not.bcompatible) then
      call lsysbl_resizeVectorBlock(p_rsolutionOld, rsolution, .false.)
      call lsysbl_resizeVectorBlock(p_rrhs, rsolution, .false.)
      call lsysbl_resizeVectorBlock(p_raux, rsolution, .false.)
    end if

    ! Save the given solution vector
    call lsysbl_copyVector(rsolution, p_rsolutionOld)


    ! Set pointer to second temporal vector for the solution computed by local substepping
    if (rtimestep%iadaptTimestep .eq. SV_TIMESTEP_AUTOADAPT) then
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
      rtimestep%dTime = rtimestep%dTime+rtimestep%dStep
      rtimestep%nSteps= rtimestep%nSteps+1

      if (rtimestep%ioutputLevel .ge. SV_IOLEVEL_INFO) then
        call output_lbrk()
        call output_separator(OU_SEP_AT)
        call output_line('Explicit Runge-Kutta scheme, Time = '//trim(sys_sdEL(rtimestep%dTime,5))//&
                         ' Stepsize = '//trim(sys_sdEL(rtimestep%dStep,5)))
        call output_separator(OU_SEP_AT)
        call output_lbrk()
      end if

      
      ! Perform multi-step Runge-Kutta method
      do istep = 1, rtimestep%multisteps
         
        if (rtimestep%ioutputLevel .ge. SV_IOLEVEL_VERBOSE) then
          call output_lbrk()
          call output_separator(OU_SEP_AT)
          call output_line('Explicit Runge-Kutta step '//trim(sys_siL(istep,5)))
          call output_separator(OU_SEP_AT)
          call output_lbrk()
        end if
         
         ! Compute the new right-hand side
         call fcb_calcRHS(rproblemLevel, rtimestep, p_rsolver, rsolution,&
                          p_rsolutionOld, p_rrhs, istep, rcollection)
         
         ! Apply given right-hand side vector
         if (present(rrhs)) call lsysbl_vectorLinearComb(rrhs, p_rrhs, 1._DP, 1._DP)
         
         ! Impose boundary conditions
         call fcb_setBoundary(rproblemLevel, rtimestep, rsolver,&
                              rsolution, p_rrhs, p_rsolutionOld, rcollection)

         ! Call linear multigrid to solve A*u=rhs
         call lsysbl_clearVector(p_raux)
         call linsol_solveMultigrid(rproblemLevel, p_rsolver, p_raux, p_rrhs)

         ! Add increment to solution vector
         call lsysbl_vectorLinearComb(p_raux, p_rsolutionOld, 1._DP, 1._DP, rsolution)
      end do


      ! Compute reference solution by performing two time steps of size Dt/2
      if (rtimestep%iadaptTimestep .eq. SV_TIMESTEP_AUTOADAPT) then

        ! Set time step to smaller value
        rtimestep%dStep = rtimestep%dStep/2.0_DP

        if (rtimestep%ioutputLevel .ge. SV_IOLEVEL_INFO) then
          call output_lbrk()
          call output_separator(OU_SEP_AT)
          call output_line('First substep in automatic time step control')
          call output_separator(OU_SEP_AT)
          call output_lbrk()
        end if


        ! Perform multi-step Runge-Kutta method for first step
        do istep = 1, rtimestep%multisteps
          
          if (rtimestep%ioutputLevel .ge. SV_IOLEVEL_VERBOSE) then
            call output_lbrk()
            call output_separator(OU_SEP_AT)
            call output_line('Explicit Runge-Kutta step '//trim(sys_siL(istep,5)))
            call output_separator(OU_SEP_AT)
            call output_lbrk()
          end if
          
          ! Compute the new right-hand side
          call fcb_calcRHS(rproblemLevel, rtimestep, p_rsolver, p_rsolutionRef,&
                           p_rsolutionOld, p_rrhs, istep, rcollection)
          
          ! Apply given right-hand side vector
          if (present(rrhs)) call lsysbl_vectorLinearComb(rrhs, p_rrhs, 1._DP, 1._DP)
          
          ! Impose boundary conditions
          call fcb_setBoundary(rproblemLevel, rtimestep, rsolver,&
                               p_rsolutionRef, p_rrhs, p_rsolutionOld, rcollection)
          
          ! Call linear multigrid to solve A*u=rhs
          call lsysbl_clearVector(p_raux)
          call linsol_solveMultigrid(rproblemLevel, p_rsolver, p_raux, p_rrhs)
          
          ! Add increment to solution vector
          call lsysbl_vectorLinearComb(p_raux, p_rsolutionOld, 1._DP, 1._DP, p_rsolutionRef)
        end do

        
        ! Save intermediate solution
        call lsysbl_copyVector(p_rsolutionRef, p_rsolutionAux)
        
        if (rtimestep%ioutputLevel .ge. SV_IOLEVEL_INFO) then
          call output_lbrk()
          call output_separator(OU_SEP_AT)
          call output_line('Second substep in automatic time step control')
          call output_separator(OU_SEP_AT)
          call output_lbrk()
        end if

        
        ! Perform multi-step Runge-Kutta method for first step
        do istep = 1, rtimestep%multisteps
          
          if (rtimestep%ioutputLevel .ge. SV_IOLEVEL_VERBOSE) then
            call output_lbrk()
            call output_separator(OU_SEP_AT)
            call output_line('Explicit Runge-Kutta step '//trim(sys_siL(istep,5)))
            call output_separator(OU_SEP_AT)
            call output_lbrk()
          end if

          
          ! Compute the new right-hand side
          call fcb_calcRHS(rproblemLevel, rtimestep, p_rsolver, p_rsolutionRef,&
                           p_rsolutionAux, p_rrhs, istep, rcollection)
          
          ! Apply given right-hand side vector
          if (present(rrhs)) call lsysbl_vectorLinearComb(rrhs, p_rrhs, 1._DP, 1._DP)
          
          ! Impose boundary conditions
          call fcb_setBoundary(rproblemLevel, rtimestep, rsolver,&
                               p_rsolutionRef, p_rrhs, p_rsolutionAux, rcollection)
          
          ! Call linear multigrid to solve A*u=rhs
          call lsysbl_clearVector(p_raux)
          call linsol_solveMultigrid(rproblemLevel, p_rsolver, p_raux, p_rrhs)
          
          ! Add increment to solution vector
          call lsysbl_vectorLinearComb(p_raux, p_rsolutionAux, 1._DP, 1._DP, p_rsolutionRef)
        end do

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

      if (rtimestep%ioutputlevel .ge. SV_IOLEVEL_VERBOSE) then
        call output_lbrk()
        call output_separator(OU_SEP_TILDE)
        call output_line('Time step was     '//merge('!!! rejected !!!','accepted        ',breject))
        call output_line('New stepsize:     '//trim(sys_sdEL(rtimestep%dStep,5)))
        call output_line('Last stepsize:    '//trim(sys_sdEL(rtimestep%dStep1,5)))
        call output_line('Relative changes: '//trim(sys_sdEL(rtimestep%drelChange,5)))
        call output_separator(OU_SEP_TILDE)
        call output_lbrk()
      end if

      ! Write time step to file?
      if (rtimestep%ioutputLevel .ge. SV_IOLEVEL_FILE) then
        write(rtimestep%iunitLogfile, FMT='("(01),",F20.5,3(",",E16.8E3))')&
            rtimestep%dTime, rtimestep%dStep, rtimestep%dStep/rtimestep%dStep1,&
            rtimestep%drelChange
      end if


      ! Do we have to reject to current solution?
      if (breject) then
        ! Yes, so restore the old solution and
        ! repeat the adaptive time-stepping loop
        call lsysbl_copyVector(p_rsolutionOld, rsolution)
      else
        ! No, accept current solution and 
        ! exit adaptive time-stepping loop
        exit timeadapt
      end if
    end do timeadapt
  end subroutine tstep_performRKStepBlock

  ! *****************************************************************************
  
!<subroutine>

  subroutine tstep_performPseudoSteppingBlock(rproblemLevel, rtimestep, rsolver, rsolution,&
                                              fcb_calcRHS, fcb_calcResidual, fcb_calcJacobian,&
                                              fcb_applyJacobian, fcb_setBoundary, rcollection, rrhs)

!<description>
    ! This subroutine performs pseudo time-stepping to compute the
    ! steady state solution to a stationary problem. Each pseudo
    ! time step can be performed by the forward Euler (theta=0), 
    ! backward Euler (theta=1), the Crank-Nicolson (theta=0.5) 
    ! time stepping algorithm and an explicit Runge-Kutta scheme.
!</description>

!<input>
    ! Callback routines
    include 'intf_nlsolcallback.inc'

    ! OPTIONAL: right-hand side vector
    type(t_vectorBlock), intent(IN), optional :: rrhs
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep

    ! solver structure
    type(t_solver), intent(INOUT), target :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>


    ! Infinite time loop
    timeloop: do

      ! What time-stepping scheme should be used?
      select case(rtimestep%ctimestepType)
        
      case (SV_RK_SCHEME)
        
        ! Adopt explicit Runge-Kutta scheme
        call tstep_performRKStep(rproblemLevel, rtimestep,&
                                 rsolver, rsolution, fcb_calcRHS,&
                                 fcb_setBoundary, rcollection, rrhs)
        
      case (SV_THETA_SCHEME)
        
        ! Adopt two-level theta-scheme
        call tstep_performThetaStep(rproblemLevel, rtimestep,&
                                    rsolver, rsolution, fcb_calcResidual,&
                                    fcb_calcJacobian, fcb_applyJacobian,&
                                    fcb_setBoundary, rcollection, rrhs)
          
      case DEFAULT
        call output_line('Unsupported time-stepping algorithm!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'tstep_performPseudoSteppingBlock')
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

  end subroutine tstep_performPseudoSteppingBlock

  ! *****************************************************************************
  
!<subroutine>

  subroutine tstep_performPseudoSteppingScalar(rproblemLevel, rtimestep, rsolver, rsolution,&
                                               fcb_calcRHS, fcb_calcResidual, fcb_calcJacobian,&
                                               fcb_applyJacobian, fcb_setBoundary, rcollection, rrhs)

!<description>
    ! This subroutine performs pseudo time-stepping to compute the
    ! steady state solution to a stationary problem. Each pseudo
    ! time step can be performed by the forward Euler (theta=0), 
    ! backward Euler (theta=1), the Crank-Nicolson (theta=0.5) 
    ! time stepping algorithm and an explicit Runge-Kutta scheme.
!</description>

!<input>
    ! Callback routines
    include 'intf_nlsolcallback.inc'

    ! OPTIONAL: right-hand side vector
    type(t_vectorScalar), intent(IN), optional :: rrhs
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep

    ! solver structure
    type(t_solver), intent(INOUT), target :: rsolver

    ! solution vector
    type(t_vectorScalar), intent(INOUT) :: rsolution

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! Infinite time loop
    timeloop: do

      ! What time-stepping scheme should be used?
      select case(rtimestep%ctimestepType)
        
      case (SV_RK_SCHEME)
        
        ! Adopt explicit Runge-Kutta scheme
        call tstep_performRKStep(rproblemLevel, rtimestep,&
                                 rsolver, rsolution, fcb_calcRHS,&
                                 fcb_setBoundary, rcollection, rrhs)
        
      case (SV_THETA_SCHEME)
        
        ! Adopt two-level theta-scheme
        call tstep_performThetaStep(rproblemLevel, rtimestep,&
                                    rsolver, rsolution, fcb_calcResidual,&
                                    fcb_calcJacobian, fcb_applyJacobian,&
                                    fcb_setBoundary, rcollection, rrhs)
          
      case DEFAULT
        call output_line('Unsupported time-stepping algorithm!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'tstep_performPseudoSteppingScalar')
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

  end subroutine tstep_performPseudoSteppingScalar


  ! *****************************************************************************
  ! AUXILIARY ROUTINES
  ! *****************************************************************************


!<function>
  
  function tstep_checkTimestep(rtimestep, rsolver,&
                               rsolution1, rsolution2) result(breject)

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
    type(t_solver), intent(IN):: rsolver

    ! first solution vector
    type(t_vectorBlock), intent(IN) :: rsolution1

    ! second soluion vector
    type(t_vectorBlock), intent(IN) :: rsolution2
!</input>

!<inputoutput>
    ! time-stepping algorithm
    type(t_timestep), intent(INOUT) :: rtimestep
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
               rsolver%istatus .eq. SV_INCR_DEF)
    
    if (breject) then
      
      ! If the time step is already equal to the smallest 
      ! admissible time step, then the simulation is terminated.
      if (rtimestep%dStep .le. rtimestep%dminStep + SYS_EPSREAL) then
        call output_line('Time step reached smallest admissible value!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'tstep_checkTimestep')
        call sys_halt()
      end if
      
      ! Undo the last time step and adjust time step size
      rtimestep%nrejectedSteps = rtimestep%nrejectedSteps+1
      rtimestep%nSteps         = rtimestep%nSteps-1
      rtimestep%dTime          = rtimeStep%dTime-rtimeStep%dStep
      rtimestep%dStep          = max(rtimestep%dminStep,&
                                     rtimestep%dstepReductionFactor*rtimestep%dStep)

      ! That's all, recompute the last step
      return
    end if

    
    !---------------------------------------------------------------------------
    ! Perform adaptive time-stepping?
    !---------------------------------------------------------------------------   

    select case(rtimestep%iadaptTimestep)

    case (SV_TIMESTEP_AUTOADAPT)
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
      dChange = lalg_errorNormDble(p_Ddata1, p_Ddata2, rtimestep%isolNorm)

      ! Calculate the 'optimal' time step size
      dStepOpt = sqrt(3.0_DP*p_rautoController%depsRel*rtimestep%dStep/dChange)
      
      ! Check if optimal time step decreased too much
      if (dStepOpt .le. p_rautoController%dDecreaseFactor*rtimestep%dStep) then

        ! Adjust time step accordingly
        rtimestep%nrejectedSteps = rtimestep%nrejectedSteps+1
        rtimestep%nSteps         = rtimestep%nSteps-1
        rtimestep%dTime          = rtimeStep%dTime-rtimeStep%dStep
        rtimestep%dStep          = max(rtimestep%dminStep, dStepOpt)

        ! That's it, recompute the last step
        breject = .true.
      else
        
        ! Accept the current solution 
        breject = .false.
        
        ! Impose upper/lower bounds on absolute values
        rtimestep%dStep1 = rtimestep%dStep
        rtimestep%dStep  = max(rtimestep%dminStep, min(rtimestep%dmaxStep, dStepOpt))

        ! Calculate the relative changes for statistical information
        rtimestep%drelChange = dChange/max(SYS_EPSREAL,&
                                           lalg_normDble(p_Ddata1, rtimestep%isolNorm))

      end if


    case (SV_TIMESTEP_PIDADAPT)
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
      dChange = lalg_errorNormDble(p_Ddata1, p_Ddata2, rtimestep%isolNorm) /&
                     lalg_NormDble(p_Ddata1, rtimestep%isolNorm)
      rtimestep%drelChange = dChange

      ! Check if solver did not converge, that is, the convergence rate equals
      ! unity or the normalized changes exceed maximum tolerance
      breject = (dChange .gt. p_rpidController%dmaxRel) .or.&
                (rsolver%istatus .eq. SV_DIVERGED)
      
      ! If the solver failed due to a non-critical error, we may accept
      ! the time step in case it cannot be further reduces
      if (breject .and. (rtimestep%dStep .gt. rtimestep%dminStep + SYS_EPSREAL)) then
      
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
        if (dChange .le. SYS_EPSREAL) return

        if (rtimestep%nSteps .gt. rtimestep%npreadaptSteps) then
          
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


    case (SV_TIMESTEP_SERADAPT)
      !-------------------------------------------------------------------------
      ! Switched evolution relaxation (SER)
      !
      ! From: W. Mulder and B. Van Leer. "Experiments with typical upwind 
      ! methods for the Euler equations", J. Comp. Phys., 59(1985), 232-246.
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
      dChange = lalg_errorNormDble(p_Ddata1, p_Ddata2, rtimestep%isolNorm)
      rtimestep%drelChange = dChange/max(SYS_EPSREAL,&
                                     lalg_normDble(p_Ddata1, rtimestep%isolNorm))
            
      
    case DEFAULT
      ! Accept time step unconditionally
      breject = .false.  

      ! Calculate the relative changes for statistical information
      call lsysbl_getbase_double(rsolution1, p_Ddata1)
      call lsysbl_getbase_double(rsolution2, p_Ddata2)
      dChange = lalg_errorNormDble(p_Ddata1, p_Ddata2, rtimestep%isolNorm)
      rtimestep%drelChange = dChange/max(SYS_EPSREAL,&
                                         lalg_normDble(p_Ddata1, rtimestep%isolNorm))
      
    end select

  end function tstep_checkTimestep

end module timestep
