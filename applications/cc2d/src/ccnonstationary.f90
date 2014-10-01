!##############################################################################
!# ****************************************************************************
!# <name> ccnonstationary </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises a time dependent solver for the coupled Navier-Stokes
!# system.
!#
!# The following routines can be found here:
!#
!# 1.) cc_initParTimeDependence
!#     -> Initialise the parameters of the time dependent solver from DAT file
!#        parameters.
!#
!# 2.) cc_initTimeSteppingScheme
!#     -> Initialise the time stepping scheme from from DAT file
!#        parameters.
!#
!# 3.) cc_performTimestep
!#     -> Calculates one time step. This routine is similar to the MGSTP
!#        routine in old FEAT applications.
!#
!# 4.) cc_solveNonstationary
!#     -> Realises a time-loop to solve multiple time steps. Solves the
!#        nonstationary Navier-Stokes system.
!# </purpose>
!##############################################################################

module ccnonstationary

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use linearsolverautoinitialise
  use matrixrestriction
  use paramlist
  use timestepping

  use collection
  use convection

  use ccbasic
  use cccallback

  use ccnonlinearcore
  use ccnonlinearcoreinit
  use ccstationary
  use adaptivetimestep
  use cctimeanalysis
  use ccgeneraldiscretisation
  use ccpostprocessing
  use ccboundarycondition

  implicit none

!<types>

!<typeblock>

  ! Saves data about a time step. Is used to make a snapshot of the current
  ! flow/mesh situation. THat way, the solver can make backups that can be
  ! restored if necessary.

  type t_timestepSnapshot

    ! Current point in time
    real(DP)            :: dtime

    ! Current timestep
    integer             :: itimeStep

    ! Current timestep configuration
    type(t_explicitTimeStepping) :: rtimeStepping

    ! Right hand side at the current point in time
    type(t_vectorBlock) :: rrhs

    ! Solution vector at the current point in time
    type(t_vectorBlock) :: rsolution

  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initParTimeDependence (rproblem,ssection,rparams)

!<description>
  ! Initialises parameters in the problem structure according to whether the
  ! simulation is time dependent or not. Initialises the time stepping scheme,
  ! start proceduce, error bounds, etc.
  !
  ! The current simulation time is set to the initial time.
  !
  ! Note: This will not allocate an memory but only initialise the parameters
  ! in rproblem according to the parameters in rparams from the DAT file.
!</description>

!<input>
  ! A parameter list from a DAT file containing parameters that configure the
  ! time dependence.
  type(t_parlist), intent(in) :: rparams

  ! The name of the section in the parameter list containing the parameters
  ! for the time dependent simulation.
  character(LEN=*), intent(in) :: ssection
!</input>

!<inputoutput>
  ! The problem structure to be initialised.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    ! Fetch the parameters. Initialise with standard settings if they do not
    ! exist.
    call parlst_getvalue_int (rparams,ssection,"itimedependence",  &
        rproblem%itimedependence, 0)
    call parlst_getvalue_int (rparams,ssection,"niterations",     &
        rproblem%rtimedependence%niterations, 1000)
    call parlst_getvalue_double (rparams,ssection,"dtimeInit",   &
        rproblem%rtimedependence%dtimeInit, 0.0_DP)
    call parlst_getvalue_double (rparams,ssection,"dtimemax",     &
        rproblem%rtimedependence%dtimeMax, 20.0_DP)
    call parlst_getvalue_double (rparams,ssection,"dtimeStep",    &
        rproblem%rtimedependence%dtimeStep, 0.01_DP)
    call parlst_getvalue_double (rparams,ssection,"dminTimeDerivative",    &
        rproblem%rtimedependence%dminTimeDerivative, 0.00001_DP)
    rproblem%rtimedependence%itimeStep = 0
    rproblem%rtimedependence%dtime = rproblem%rtimedependence%dtimeInit

    ! Call the initialisation routine for the adaptive time stepping to
    ! initialise the rest of the parameters.
    call adtstp_init (rparams,ssection,&
        rproblem%rtimedependence%radaptiveTimeStepping)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_backupTimestep (rsnapshot,rproblem,&
      rtimeStepping,rsolution,rrhs)

!<description>
  ! Makes a backup of the current flow configuration into rsnapshot.
  ! If there is any previous data in rsnapshot, it is overwritten/updated.
!</description>

!<inputoutput>
  ! A snapshot structure that receives the snapshot.
  type(t_timestepSnapshot), intent(inout) :: rsnapshot

  ! The problem structure with all information about the current problem.
  type(t_problem), intent(inout) :: rproblem

  ! A time stepping structure that configures the current time step.
  type(t_explicitTimeStepping), intent(inout), optional :: rtimeStepping

  ! Right hand side at the current point in time
  type(t_vectorBlock), intent(inout), optional :: rrhs

  ! Solution vector at the current point in time
  type(t_vectorBlock), intent(inout), optional :: rsolution
!</inputoutput>

!</subroutine>

    ! Back up all data into rsnapshot.
    rsnapshot%dtime = rproblem%rtimedependence%dtime
    rsnapshot%itimestep = rproblem%rtimedependence%itimestep

    if (present(rtimeStepping)) &
      rsnapshot%rtimeStepping = rtimeStepping

    if (present(rrhs)) &
      call lsysbl_copyVector (rrhs,rsnapshot%rrhs)

    if (present(rsolution)) &
      call lsysbl_copyVector (rsolution,rsnapshot%rsolution)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_restoreTimestep (rsnapshot,rproblem,rtimeStepping,&
      rsolution,rrhs)

!<description>
  ! Restores a backup of a time step from rsnapshot.
!</description>

!<inputoutput>
  ! A snapshot structure that receives the snapshot.
  type(t_timestepSnapshot), intent(inout) :: rsnapshot

  ! The problem structure with all information about the current problem.
  type(t_problem), intent(inout) :: rproblem

  ! A time stepping structure that configures the current time step.
  type(t_explicitTimeStepping), intent(inout), optional :: rtimeStepping

  ! Right hand side at the current point in time
  type(t_vectorBlock), intent(inout), optional :: rrhs

  ! Solution vector at the current point in time
  type(t_vectorBlock), intent(inout), optional :: rsolution

!</inputoutput>

!</subroutine>

    ! Restore all data from rsnapshot.
    rproblem%rtimedependence%dtime = rsnapshot%dtime
    rproblem%rtimedependence%itimestep = rsnapshot%itimestep

    if (present(rtimeStepping)) &
      rtimeStepping = rsnapshot%rtimeStepping

    if (present(rrhs)) &
      call lsysbl_copyVector (rsnapshot%rrhs,rrhs)

    if (present(rsolution)) &
      call lsysbl_copyVector (rsnapshot%rsolution,rsolution)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_releaseSnapshot (rsnapshot)

!<description>
  ! Releases the rsnapshot structure. All allocated memory is released.
!</description>

!<inputoutput>
  ! A snapshot structure that receives the snapshot.
  type(t_timestepSnapshot), intent(inout) :: rsnapshot
!</inputoutput>

!</subroutine>

    ! Release allocated memory -- if memory is allocated.
    if (rsnapshot%rrhs%NEQ .ne. 0) &
      call lsysbl_releaseVector (rsnapshot%rrhs)
    if (rsnapshot%rsolution%NEQ .ne. 0) &
      call lsysbl_releaseVector (rsnapshot%rsolution)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initTimeSteppingScheme (rparams,rtimeStepping,ipressureFullyImplicit)

!<description>
  ! Initialises the time stepping scheme according to the parameters in the DAT file.
!</description>

!<input>
  ! A parameter list from a DAT file containing parameters that configure the
  ! time dependence.
  type(t_parlist), intent(in) :: rparams
!</input>

!<inputoutput>
  ! The time stepping scheme structure to be initialised.
  type(t_explicitTimeStepping), intent(inout) :: rtimeStepping

  ! Returns 1 if the pressure should be discretised fully implicitely, 0 otherwise.
  integer, intent(out) :: ipressureFullyImplicit
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: cscheme
    real(DP) :: dtheta,dtstep,dtimemin

    ! Get the parameters for the time stepping scheme from the parameter list
    call parlst_getvalue_int (rparams, &
        "TIME-DISCRETISATION", "ITIMESTEPSCHEME", cscheme, 0)
    call parlst_getvalue_double (rparams, &
        "TIME-DISCRETISATION", "DTIMESTEPTHETA", dtheta, 1.0_DP)
    call parlst_getvalue_double (rparams, &
        "TIME-DISCRETISATION", "DTIMESTEP", dtstep, 0.1_DP)
    call parlst_getvalue_double (rparams, &
        "TIME-DISCRETISATION", "DTIMEINIT", dtimemin, 0.0_DP)

    ! Initialise the time stepping in the problem structure
    if (cscheme .lt. 0) then
      cscheme = 0
    end if
    call timstp_init (rtimestepping, &
                      cscheme, dtimemin, dtstep, dtheta)

    ! Discretisation of the pressure.
    call parlst_getvalue_int (rparams, &
        "TIME-DISCRETISATION", "IPRESSUREFULLYIMPLICIT", ipressureFullyImplicit, 1)

    ! Force-off flag for DIRK schemes: treat pressure semi-implicitly, as the velocity
    if (rtimestepping%ctimestepType .eq. TSCHM_FS_DIRK  .or. &
        rtimestepping%ctimestepType .eq. TSCHM_DIRK23L  .or. &
        rtimestepping%ctimestepType .eq. TSCHM_DIRK34La .or. &
        rtimestepping%ctimestepType .eq. TSCHM_DIRK34Lb .or. &
        rtimestepping%ctimestepType .eq. TSCHM_DIRK44L  .or. &
        rtimestepping%ctimestepType .eq. TSCHM_DIRK54L  .or. &
        rtimestepping%ctimestepType .eq. TSCHM_SDIRK2   .or. &
        rtimestepping%ctimestepType .eq. TSCHM_SDIRK3PR) then
      ipressureFullyImplicit = 0
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_performTimestep (rproblem,rvector,rrhs,&
      rtimestepping,ipressureFullyImplicit,rnonlinearIteration,rnlSolver,&
      rtempVector,rtempVectorRhs)

!<description>
  ! Performs one time step: $t^n -> t^n+1$.
  ! Assembles system matrix and RHS vector.
  ! Solves the corresponding time-step equation and returns the solution vector
  ! at the end of the time step.
  ! Solves the given problem by applying a nonlinear solver.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem

  ! The current solution vector at time $t^n$. Is replaced by the
  ! solution vector at time $t^{n+1}.
  type(t_vectorBlock), intent(inout) :: rvector

  ! The RHS vector at time $t^n$. Is replaced by the RHS at time $t^{n+1}$.
  type(t_vectorBlock), intent(inout) :: rrhs

  ! Configuration block of the time stepping scheme.
  type(t_explicitTimeStepping)        :: rtimestepping

  ! Set to 1 if the pressure should be discretised fully implicitely, 0 otherwise.
  integer, intent(in) :: ipressureFullyImplicit

  ! Structure for the nonlinear iteration for solving the core equation.
  type(t_ccnonlinearIteration), intent(inout) :: rnonlinearIteration

  ! A configuration stucture for the nonlinear solver
  type(t_nlsolNode) :: rnlSolver

  ! Temporary vectors for the nonlinear solver
  type(t_vectorBlock), intent(inout) :: rtempVector,rtempVectorRhs
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    type(t_ccnonlinearIteration) :: rnonlinearIterationTmp
    type(t_nonlinearCCMatrix) :: rnonlinearCCMatrix
    type(t_timer) :: rtimerRHSgeneration

    ! auxiliary vector for Glowinski2003-FS-theta scheme
    type(t_vectorBlock), save :: rsolutionAux

    ! auxiliary vectors for diagonally implicit Runge-Kutta schemes:
    !
    ! <tex>$ f(t_{m} + c_{j} k) - N(U_{j}) U_{j} - B P_{j} $,
    ! see \cite[p. 517]{Rang200702}</tex> needed in right hand side assembly
    type(t_vectorBlock), dimension(5), save :: rvectorAuxFmjUjplusBPj
    ! U_0, solution from last macro time step
    type(t_vectorBlock), save :: rsolutionFromLastMacroTimeStep

    ! DEBUG!!!
    !REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2


    if (rtimestepping%ctimestepType .eq. TSCHM_FS_DIRK  .or. &
        rtimestepping%ctimestepType .eq. TSCHM_DIRK23L  .or. &
        rtimestepping%ctimestepType .eq. TSCHM_DIRK34La .or. &
        rtimestepping%ctimestepType .eq. TSCHM_DIRK34Lb .or. &
        rtimestepping%ctimestepType .eq. TSCHM_DIRK44L  .or. &
        rtimestepping%ctimestepType .eq. TSCHM_DIRK54L  .or. &
        rtimestepping%ctimestepType .eq. TSCHM_SDIRK2   .or. &
        rtimestepping%ctimestepType .eq. TSCHM_SDIRK3PR) then
      ! explicit first stage: use rtimestepping%isubstep - 1
      ! implicit first stage: use rtimestepping%isubstep
      select case (rtimestepping%isubstep - &
                   merge(1,0,rtimestepping%bexplicitFirstStage .eqv. .TRUE.))
      case (1)
        ! Store (velocity) solution from macro time step t_n for later use in subsequent
        ! stages of DIRK scheme where it is needed to assembly the right hand side
        ! vector.
        call lsysbl_createVector (rvector,rsolutionFromLastMacroTimeStep,bclear=.false.)
        call lsysbl_copyVector (rvector, rsolutionFromLastMacroTimeStep)

        call lsysbl_createVector (rvector, rvectorAuxFmjUjplusBPj(1), bclear=.false.)
      case (2)
        call lsysbl_createVector (rvector, rvectorAuxFmjUjplusBPj(2), bclear=.false.)
      case (3)
        call lsysbl_createVector (rvector, rvectorAuxFmjUjplusBPj(3), bclear=.false.)
      case (4)
        call lsysbl_createVector (rvector, rvectorAuxFmjUjplusBPj(4), bclear=.false.)
      case (5)
        call lsysbl_createVector (rvector, rvectorAuxFmjUjplusBPj(5), bclear=.false.)
      end select
    end if

    if (rtimestepping%ctimestepType .eq. TSCHM_FS_GLOWINSKI .and. &
        rtimestepping%isubstep .eq. 1) then
      ! Create auxiliary vector we need fir steps 2 and 3b
      call lsysbl_createVector (rvector, rsolutionAux, .false.)
      call lsysbl_copyVector (rvector, rsolutionAux)
    end if

    if (rtimestepping%ctimestepType .eq. TSCHM_FS_GLOWINSKI .and. &
        rtimestepping%isubstep .eq. 2) then

      ! Switch to the next point in time.
      rproblem%rtimedependence%dtime = rtimestepping%dcurrenttime + rtimestepping%dtstep

      ! u^{n+1-theta} = (1-theta)/theta u^{n+theta} + (2theta - 1)/theta u^{n}
      call lsyssc_vectorLinearComb (&
           rsolutionAux%RvectorBlock(1), &
           rvector%RvectorBlock(1), &
           (2.0_DP*rtimestepping%dtheta - 1.0_DP)/rtimestepping%dtheta, &
           (1.0_DP - rtimestepping%dtheta)/rtimestepping%dtheta)
      call lsyssc_vectorLinearComb (&
           rsolutionAux%RvectorBlock(2), &
           rvector%RvectorBlock(2), &
           (2.0_DP*rtimestepping%dtheta - 1.0_DP)/rtimestepping%dtheta, &
           (1.0_DP - rtimestepping%dtheta)/rtimestepping%dtheta)
      ! No solving needed, mark solver as successful
      rnlSolver%iresult     = 0
      rnlSolver%iiterations = 0
      rnlSolver%DinitialDefect      = 0.0_DP
      rnlSolver%DfinalDefect        = 0.0_DP
      rnlSolver%dinitialDefectTotal = 0.0_DP
      rnlSolver%dfinalDefectTotal   = 0.0_DP
      rnlSolver%dtimeTotal          = 0.0_DP
      rnlSolver%dtimeNLdefect       = 0.0_DP
      rnlSolver%dtimeNLpreconditioning = 0.0_DP

    else
      if (rtimestepping%ctimestepType .eq. TSCHM_FS_DIRK  .or. &
          rtimestepping%ctimestepType .eq. TSCHM_DIRK23L  .or. &
          rtimestepping%ctimestepType .eq. TSCHM_DIRK34La .or. &
          rtimestepping%ctimestepType .eq. TSCHM_DIRK34Lb .or. &
          rtimestepping%ctimestepType .eq. TSCHM_DIRK44L  .or. &
          rtimestepping%ctimestepType .eq. TSCHM_DIRK54L  .or. &
          rtimestepping%ctimestepType .eq. TSCHM_SDIRK2   .or. &
          rtimestepping%ctimestepType .eq. TSCHM_SDIRK3PR) then
! Note: inhomogeneous Neumann boundary conditions not supported (as of yet)!

        ! The new RHS will be set up in rtempVectorRhs. Assign the discretisation/
        ! boundary conditions of rrhs to that vector so that rtempVectorRhs
        ! acts as a RHS vector.
        call lsysbl_assignDiscrIndirect(rrhs, rtempVectorRhs)

        ! Algorithm:
        ! in substep 1:
        !     re-use pre-assembled f_n, assemble N(u_n) and store the result of
        !               f_n - N(u_n) u_n - B p_n
        !     for re-use in subsequent substeps. Then assemble f(t_n + c2 k) and calculate
        !        M u_n + k a21 (f(t_n) - N(u_n) u_n - B p_n)
        !              + k a22 f(t_n + c2 k)
        !     as right hand side vector
        ! in substep 2:
        !     re-use pre-assembled f(t_n + c2 k), assemble N(u(t_n + c2 k)) and store the
        !     result of
        !         f(t_n + c2 k) - N(u(t_n + c2 k)) u(t_n + c2 k) - B p(t_n + c2 k)
        !     for re-use in subsequent substeps. Then assemble f(t_n + c3 k) and calculate
        !        M u_n + k a31 (f(t_n) - N(u_n) u_n - B p_n)
        !              + k a32 (f(t_n + c2 k)
        !                       - N(u(t_n + c2 k)) u(t_n + c2 k) - B p(t_n + c2 k))
        !              + k a33 f(t_n + c3 k)
        ! in substep i:
        !     re-use pre-assembled f(t_n + c_{i-1} k), assemble N(u(t_n + c_{i-1} k)).
        !     Then assemble f(t_n + c_i k) = f(t_n+1) and calculate
        !        M u_n + k a_{i,1} (f(t_n + c1 k)
        !                           - N(u(t_n + c1 k)) u(t_n + c1 k) - B p(t_n + c1 k))
        !              + ...
        !              + k a_{i,i-1} (f(t_n + c_i k)
        !                           - N(u(t_n + c_i k)) u(t_n + c_i k) - B p(t_n + c_i k))
        !              + k a_{i,i} f(t_n+1)

        ! f = M u_{time before stage 1} +
        !     \tau \sum_{j=1}^{i-1} a_{i,j} (f(t_n + c_j k)
        !                                    - N(u(t_n + c_j k)) u(t_n + c_j k)
        !                                    - B p(t_n + c_j k))
        ! (i = current substep aka stage).
        !
        ! step 1:
        !    f = M u_{time before stage 1}
        call cc_initNonlinMatrix (rnonlinearCCMatrix,rproblem,&
             rproblem%RlevelInfo(rproblem%NLMAX)%rdiscretisation,&
             rproblem%RlevelInfo(rproblem%NLMAX)%rasmTempl,&
             rproblem%RlevelInfo(rproblem%NLMAX)%rdynamicInfo)
        call lsyssc_matVec (rnonlinearCCMatrix%p_rasmTempl%rmatrixMass, &
             rsolutionFromLastMacroTimeStep%RvectorBlock(1), &
             rtempVectorRhs%RvectorBlock(1), &
             1.0_DP, 0.0_DP)
        call lsyssc_matVec (rnonlinearCCMatrix%p_rasmTempl%rmatrixMass, &
             rsolutionFromLastMacroTimeStep%RvectorBlock(2), &
             rtempVectorRhs%RvectorBlock(2), &
             1.0_DP, 0.0_DP)
        if (rtimestepping%isubstep .gt. 1) then
          ! step 2:
          !   Determine (f(t_n + c_{i-1} k) - N(u(t_n + c_{i-1} k)) u(t_n + c_{i-1} k)
          !                                 - B p(t_n + c_{i-1} k))
          ! (i = current substep aka stage).
          !
          ! rvectorAuxFmjUjplusBPj(rtimestepping%isubstep-1) = rrhs
          call lsysbl_vectorLinearComb(&
               rrhs, rvectorAuxFmjUjplusBPj(rtimestepping%isubstep-1), &
               1.0_DP, 0.0_DP)
          ! Set up N(u(t_n + c_{i-1} k)) in order to calculate
          !      f_(t_n + c_{i-1} k) - N(u(t_n + c_{i-1} k)) u(t_n + c_{i-1} k)
          !                          - B p(t_n + c_{i-1} k),
          ! (i = current substep aka stage).
          ! Note that there is no need to scale the pressure here explicitly, as opposed
          ! to the case for the general theta and fractional-step scheme.
          ! Do not implement any boundary conditions when assembling this, they are
          ! implemented at the end when the full RHS is finished...
          rnonlinearCCMatrix%dmass = 0.0_DP
          rnonlinearCCMatrix%dstokes = -1.0_DP
          rnonlinearCCMatrix%dconvection = -real(1-rproblem%rphysics%iequation,DP)
          rnonlinearCCMatrix%dgradient = -1.0_DP
          rnonlinearCCMatrix%ddivergence = 0.0_DP
          ! Now store the final result of
          !      f_(t_n + c_{i-1} k) - N(u(t_n + c_{i-1} k)) u(t_n + c_{i-1} k)
          !                          - B p(t_n + c_{i-1} k)
          ! (i = current substep aka stage).
          call cc_nonlinearMatMul (&
               rnonlinearCCMatrix, rvector, &
               rvectorAuxFmjUjplusBPj(rtimestepping%isubstep-1), 1.0_DP, 1.0_DP, rproblem)

          ! step 3:
          ! f = M u_{time before stage 1} +
          !     \tau \sum_{j=1}^{i-1} a_ij (f(t_n + c_j k)
          !                                 - N(u(t_n + c_j k)) u(t_n + c_j k)
          !                                 - B p(t_n + c_j k))
          ! (i = current substep aka stage).
          do i = 1, rtimestepping%isubstep-1
            ! stage 1: do nothing
            ! stage 2: rvectorAuxFmjUjplusBPj(1) got filled in this stage, use it
            ! stage i, i>=3: rvectorAuxFmjUjplusBPj(i) got filled in this stage, use
            !          rvectorAuxFmjUjplusBPj(<=i)
            call lsysbl_vectorLinearComb(&
                 rvectorAuxFmjUjplusBPj(i), rtempVectorRhs, &
                 rtimestepping%dcoeffA(rtimestepping%isubstep,i) * rtimestepping%dtau, &
                 1.0_DP)
          end do
        end if

        ! -------------------------------------------
        ! Switch to the next point in time where the right hand side gets evaluated
        rproblem%rtimedependence%dtime = &
             rtimestepping%dtimeDIRKstage1 + &
             rtimestepping%dcoeffC(rtimestepping%isubstep) * rtimestepping%dtau

        ! Discretise the boundary conditions at the new point in time --
        ! if the boundary conditions are nonconstant in time!
        if (rproblem%iboundary .ne. 0) then
          call cc_updateDiscreteBC (rproblem)
        end if

        ! -------------------------------------------

        ! Generate (f_n+1, phi) into the rrhs overwriting the previous rhs.
        ! Do not implement any BC`s! We need the "raw" RHS for the next timestep.
        call stat_clearTimer(rtimerRHSgeneration)
        call stat_startTimer(rtimerRHSgeneration)
        call cc_generateBasicRHS (rproblem,&
             rproblem%RlevelInfo(rproblem%NLMAX)%rasmTempl,&
             rproblem%rrhsassembly,rrhs)

        ! Add (k a_ii * f_{n+1}, phi) with i = substep+1 to the current RHS.
        call lsysbl_vectorLinearComb(rrhs,rtempVectorRhs,&
             rtimestepping%dcoeffA(rtimestepping%isubstep,rtimestepping%isubstep) * &
             rtimestepping%dtau, 1.0_DP)

        ! Ensure f(3) is zero
        call lsyssc_clearVector(rtempVectorRhs%RvectorBlock(3))

        call stat_stopTimer(rtimerRHSgeneration)
        rproblem%rstatistics%dtimeRHSAssembly = &
             rproblem%rstatistics%dtimeRHSAssembly + rtimerRHSgeneration%delapsedReal

      else

        ! The new RHS will be set up in rtempVectorRhs. Assign the discretisation/
        ! boundary conditions of rrhs to that vector so that rtempVectorRhs
        ! acts as a RHS vector.
        call lsysbl_assignDiscrIndirect(rrhs,rtempVectorRhs)

        ! DEBUG!!!
        !CALL lsysbl_getbase_double (rvector,p_Ddata)
        !CALL lsysbl_getbase_double (rtempVectorRhs,p_Ddata2)

        ! We have an equation of the type
        !
        !   d/dt u(x,t)  +  N(u(x,t)) u(x,t) =  f(x,t)
        !
        ! Which is discretised in time with a Theta scheme, leading to
        !
        !   $$ u_{n+1} + w_1*N(u_n+1) u_{n+1}
        !      =   u_n + w_2*N(u_n) u_n  +  w_3*f_{n+1}  +  w_4*f_n $$
        !
        ! with k=time step size, u_{n+1} = u(.,t_{n+1}),etc., c.f. timestepping.f90.
        !
        ! The RHS of that equation therefore contains parts of the solution
        ! u_n, of the old RHS f_n and the new RHS f_{n+1}. At first, we make
        ! a weighted copy of the current RHS f_n to the "global" RHS vector
        ! according to the time stepping scheme.
        !
        !
        ! If we have inhomogeneous Neumann boundary conditions, the situation
        ! is slightly more complicated. The weak formulation of, e.g., the
        ! Stokes equations reads (note that there is a "k" included in the
        ! coefficients w_i, so w_i/k is constant!):
        !
        !    ( (u_n+1 - u_n)/k , phi )  +  nu w_1/k (grad u_n+1, grad phi)       -  nu w_2/k (grad u_n, grad phi)
        !                               -  nu w_1/k (du_n+1/dn , phi     )_Gamma +  nu w_2/k (du_n/dn , phi     )_Gamma
        !                               -           (p         , grad phi)
        !                               +           (p n       , phi     )_Gamma
        !  =   w_3/k ( f_n+1, phi)
        !    + w_4/k ( f_n  , phi)
        !
        ! The pressure is used fully implicitely, so the meaning of the pressure
        ! depends on the timestepping scheme used. For the CN scheme, e.g.,
        ! the above formula reads (keep in mind that FF2 stores the time step k
        ! in the pressure variable, i.e. 'p' is in fact 'k p')
        !
        !    ( (u_n+1 - u_n)/k , phi )  +  nu/2  (grad u_n+1, grad phi)       +  nu/2 (grad u_n, grad phi)
        !                               -  nu/2  (du_n+1/dn , phi     )_Gamma -  nu/2 (du_n/dn , phi     )_Gamma
        !                               -        (p_n+1/2   , grad phi)
        !                               +        (p_n+1/2 n , phi     )
        !  =   1/2 ( f_n+1  , phi)  +  1/2 ( f_n  , phi)
        !
        ! Some terms can be combined. For example, in the CN method, one could write
        !
        !    ( (u_n+1 - u_n)/k , phi )  +  nu/2   (grad u_n+1, grad phi)       +  nu/2 (grad u_n, grad phi)
        !                               -  nu     (du_n+1/2 / dn , phi )_Gamma
        !                               -         (p_n+1/2       , grad phi)
        !                               +         (p_n+1/2 n     , phi     )_Gamma
        !  =  ( f_n+1/2, phi)
        !
        ! which gives
        !
        !    ( (u_n+1 - u_n)/k , phi )  +  nu/2  (grad u_n+1, grad phi)  +  nu/2 (grad u_n, grad phi)
        !                               -        (p_n+1/2, grad phi)
        !  =  ( f_n+1/2, phi )          +        (nu du_n+1/2 / dn - p_n+1/2 n, phi)_Gamma
        !
        ! To implement inhomogeneous Neumann boundary conditions, one replaces
        ! the inhomogeneity on the RHS by the data, which results in
        !
        !  =  ( f_n+1/2  , phi           +         (g_n+1/2 , phi)_Gamma
        !
        ! so one has "g_n+1/2 = nu du_n+1/2 / dn - p_n+1/2 n", and as a consequence, the
        ! inhomogeneity has to be evaluated at the midpoint in time. Alternatively, both
        ! parts can be calculated with the trapezoidal rule (approximating the midpoint
        ! rule), so one ends up with
        !
        !  =  ( (f_n+1 + f_n)/2  , phi)  +  ( (g_n+1 + g_n)/2 , phi)_Gamma
        !
        ! Similar arguments can also be used in the general case. Here, one has to
        ! assemble
        !
        !  =  ( w_3 f_n+1 + w_4 f_n  , phi)  +  ( w_1 g_n+1 - w_2 g_n , phi)_Gamma
        !
        ! where "w_1 g_n+1 - w_2 g_n" approximates "( du/dn - p n, phi)" at the
        ! point in time corresponding to p.
        !
        !
        ! So what to do? We have "(f_n,phi)" from the last timestep and calculate
        ! "( w_3 f_n+1 + w_4 f_n  , phi )  +  ( w_1 g_n+1 - w_2 g_n , phi )_Gamma"
        ! in the following.

        ! Set up w_4*f_n.
        call lsysbl_vectorLinearComb(rrhs,rtempVectorRhs,&
             rtimestepping%dweightOldRHS,0.0_DP)

        ! Inhomogeneous Neumann part: "( -w_2 g_n, phi )_Gamma"
        call cc_assembleInhomNeumann (rproblem,&
             rproblem%rcollection,rtempVectorRhs,-rtimestepping%dweightMatrixRHS)

        ! For setting up M(u_n) + w_2*N(u_n), switch the sign of w_2 and call the method
        ! to calculate the Convection/Diffusion part of the nonlinear defect. This builds
        ! rtempVectorRhs = rtempVectorRhs - (-Mass)*u - (-w_2) (nu*Laplace*u + grad(u)u).
        ! Switch off the B-matrices as we do not need them for this defect.
        !
        ! Do not implement any boundary conditions when assembling this -- it is not
        ! a defect vector!
        ! The BC`s are implemented at the end when the full RHS is finished...

        call cc_initNonlinMatrix (rnonlinearCCMatrix,rproblem,&
             rproblem%RlevelInfo(rproblem%NLMAX)%rdiscretisation,&
             rproblem%RlevelInfo(rproblem%NLMAX)%rasmTempl,&
             rproblem%RlevelInfo(rproblem%NLMAX)%rdynamicInfo)

        if (rtimestepping%ctimestepType .eq. TSCHM_FS_GLOWINSKI) then
          rnonlinearCCMatrix%dmass = 1.0_DP
          rnonlinearCCMatrix%dstokes = 0.0_DP
          rnonlinearCCMatrix%dconvection = 0.0_DP
          rnonlinearCCMatrix%dgradient = 0.0_DP
          rnonlinearCCMatrix%ddivergence = 0.0_DP
        else
          rnonlinearCCMatrix%dmass = 1.0_DP
          rnonlinearCCMatrix%dstokes = rtimestepping%dweightMatrixRHS
          rnonlinearCCMatrix%dconvection = rtimestepping%dweightMatrixRHS * &
               real(1-rproblem%rphysics%iequation,DP)

          rnonlinearCCMatrix%dgradient = 0.0_DP
          rnonlinearCCMatrix%ddivergence = 0.0_DP
        end if

        ! Fully implicit pressure? There is only a difference if Crank-Nicolson
        ! is used.
        if (ipressureFullyImplicit .ne. 1) then
          rnonlinearCCMatrix%dgradient = rtimestepping%dweightMatrixRHS
        end if

        ! Calculate   rtempVectorRhs := rnonlinearCCMatrix rvector + rtempVectorRhs
        call cc_nonlinearMatMul (rnonlinearCCMatrix, rvector, rtempVectorRhs, &
             1.0_DP, 1.0_DP, rproblem)

        ! -------------------------------------------
        ! Switch to the next point in time.
        rproblem%rtimedependence%dtime = rtimestepping%dcurrenttime + rtimestepping%dtstep

        ! Discretise the boundary conditions at the new point in time --
        ! if the boundary conditions are nonconstant in time!
        if (rproblem%iboundary .ne. 0) then
          call cc_updateDiscreteBC (rproblem)
        end if

        ! -------------------------------------------

        ! Generate (f_n+1, phi) into the rrhs overwriting the previous rhs.
        ! Do not implement any BC`s! We need the "raw" RHS for the next timestep.
        call stat_clearTimer(rtimerRHSgeneration)
        call stat_startTimer(rtimerRHSgeneration)
        call cc_generateBasicRHS (rproblem,&
             rproblem%RlevelInfo(rproblem%NLMAX)%rasmTempl,&
             rproblem%rrhsassembly,rrhs)

        ! Add (w_3 * f_{n+1}, phi) to the current RHS.
        call lsysbl_vectorLinearComb(rrhs,rtempVectorRhs,&
             rtimestepping%dweightNewRHS,1.0_DP)

        ! Add the inhomogeneous Neumann BCs to the RHS: "(w_1 g_n , phi)_Gamma"
        call cc_assembleInhomNeumann (rproblem,&
             rproblem%rcollection,rtempVectorRhs,rtimestepping%dweightMatrixLHS)

        call stat_stopTimer(rtimerRHSgeneration)
        rproblem%rstatistics%dtimeRHSAssembly = &
             rproblem%rstatistics%dtimeRHSAssembly + rtimerRHSgeneration%delapsedReal
      end if  ! DIRK scheme or general theta/fractional-step theta scheme

      ! Implement boundary conditions into the RHS and solution vector, not
      ! into the matrices; the latter is done during the nonlinear iteration.
      call cc_implementBC (rproblem,rvector,rtempVectorRhs,.true.,.true.)

      ! That is it for the RHS and solution vector.
      !
      ! The LHS is "u_{n+1} + w_1*N(u_n+1)" which results in the system matrix
      ! "M + w_1 N(.)" for the next linear system to solve.
      ! Set up the corresponding core equation in a temporary core-equation
      ! structure.

      rnonlinearIterationTmp = rnonlinearIteration

      rnonlinearIterationTmp%dmass = 1.0_DP
      rnonlinearIterationTmp%dstokes = rtimestepping%dweightMatrixLHS
      rnonlinearIterationTmp%dconvection = rtimestepping%dweightMatrixLHS * &
          real(1-rproblem%rphysics%iequation,DP)
      rnonlinearIterationTmp%dgradient   = 1.0_DP
      rnonlinearIterationTmp%ddivergence   = 1.0_DP

      ! Scale the pressure by the length of the time step. The core equation routines
      ! handle the equation
      !   dmass*M*u + dstokes*nu*Laplace*u + dconvection*N(u)u + B*p = ...
      ! but we want to solve
      !   dmass*M*u + dstokes*nu*Laplace*u + dconvection*N(u)u + tstep*B*p = ...
      !
      ! So the trick is to scale p by tstep, solve the core equation
      !   dmass*M*u + dstokes*nu*Laplace*u + dconvection*N(u)u + B*(tstep*p) = ...
      ! and scale it back afterwards.
      !
      ! Note that there is an error in the book of [Turek] describing the factor
      ! in front of the pressure in the Crank Nicolson scheme! The pressure is
      ! handled fully implicitely. There is no part of the pressure on the RHS
      ! of the time step scheme and so the factor in front of the pressure
      ! is always the length of the current (sub)step!
      !
      ! For fully implicit pressure, just scale by the timestep size.
      ! For semi-implicit pressure, scale by the full weight of the LHS.
      ! There is only a difference if Crank-Nicolson or similar is used.
      if (ipressureFullyImplicit .ne. 1) then
        call lsyssc_scaleVector (rvector%RvectorBlock(NDIM2D+1),&
            rtimestepping%dweightMatrixLHS)
      else
        call lsyssc_scaleVector (rvector%RvectorBlock(NDIM2D+1),&
            rtimestepping%dtstep)
      end if

      ! Update the preconditioner for the case, something changed (e.g.
      ! the boundary conditions).
      ! Note: The bstructuralChange-parameter is set to FALSE here.
      ! In case the template matrices changed (e.g. during a mesh adaption),
      ! the routine must be called with bstructuralChange=true!
      call cc_updatePreconditioner (rproblem,rnonlinearIterationTmp,&
         rvector,rtempVectorRhs,.false.,.false.)

      ! Call the solver of the core equation to solve it using a nonlinear
      ! iteration.
      call cc_solveCoreEquation (rproblem,rnonlinearIterationTmp,rnlSolver,&
          rvector,rtempVectorRhs,rtempVector)

      ! scale the pressure back, then we have again the correct solution vector.
      if (ipressureFullyImplicit .ne. 1) then
        call lsyssc_scaleVector (rvector%RvectorBlock(NDIM2D+1),&
            1.0_DP/rtimestepping%dweightMatrixLHS)
      else
        call lsyssc_scaleVector (rvector%RvectorBlock(NDIM2D+1),&
            1.0_DP/rtimestepping%dtstep)
      end if

      ! rvector is the solution vector u^{n+1}.
    end if


    if (rtimestepping%ctimestepType .eq. TSCHM_FS_GLOWINSKI .and. &
        rtimestepping%isubstep .eq. 3) then
      ! Step 3b:
      !    p^{n+1} = (1-theta) p^{n+theta} + theta \tilde{p}^{n+1}
      ! with \tilde{p}^{n+1} being the pressure solution from 3rd substep of Glowinski`s
      ! fractional step scheme.
      call lsyssc_vectorLinearComb (&
           rsolutionAux%RvectorBlock(3),  rvector%RvectorBlock(3), &
           1.0_DP - rtimestepping%dtheta, rtimestepping%dtheta)

      call lsysbl_releaseVector (rsolutionAux)
    end if

    if (rtimestepping%ctimestepType .eq. TSCHM_FS_DIRK  .or. &
        rtimestepping%ctimestepType .eq. TSCHM_DIRK23L  .or. &
        rtimestepping%ctimestepType .eq. TSCHM_DIRK34La .or. &
        rtimestepping%ctimestepType .eq. TSCHM_DIRK34Lb .or. &
        rtimestepping%ctimestepType .eq. TSCHM_DIRK44L  .or. &
        rtimestepping%ctimestepType .eq. TSCHM_DIRK54L  .or. &
        rtimestepping%ctimestepType .eq. TSCHM_SDIRK2   .or. &
        rtimestepping%ctimestepType .eq. TSCHM_SDIRK3PR) then
      if (rtimestepping%isubstep .eq. rtimestepping%nsubsteps) then
        ! Free up (velocity) solution from macro time step t_n
        call lsysbl_releaseVector (rsolutionFromLastMacroTimeStep)
        do i = 5, 1, -1
          if (i .le. rtimestepping%nsubsteps - &
                     merge(1,0,rtimestepping%bexplicitFirstStage .eqv. .TRUE.)) then
            call lsysbl_releaseVector (rvectorAuxFmjUjplusBPj(i))
          end if
        end do
      end if
    end if

    ! Finally tell the time stepping scheme that we completed the time step.
    call timstp_nextSubstep (rtimestepping)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_interpolateTimesteps (rtimestepping,rvectorOld,dtimeOld,&
      rvectorNew,dtimeNew,rvectorInt,dtimeInt,ipressureFullyImplicit)

!<description>
  ! Interpolates the solution of the old and new timestep to a common
  ! point in time. May be necessary for example for a CN time discretisation
  ! where the pressure exists between the velocity timesteps.
!</description>

!<input>
  ! Configuration block of the time stepping scheme.
  type(t_explicitTimeStepping)        :: rtimestepping

  ! Solution vector of the previous timestep
  type(t_vectorBlock), intent(in) :: rvectorOld

  ! Point in time of the previous timestep
  real(dp), intent(in) :: dtimeOld

  ! Solution vector of the current timestep
  type(t_vectorBlock), intent(in) :: rvectorNew

  ! Point in time of the current timestep
  real(dp), intent(in) :: dtimeNew

  ! Set to 1 if the pressure should be discretised fully implicitely, 0 otherwise.
  integer, intent(in) :: ipressureFullyImplicit
!</input>

!<inputoutput>
  ! Interpolated solution vector
  type(t_vectorBlock), intent(inout) :: rvectorInt

  ! Point in time of the interpolated solution vector
  real(dp), intent(out) :: dtimeInt
!</inputoutput>

!</subroutine>

    real(DP) :: dfactor


    if (rtimestepping%ctimestepType .eq. TSCHM_ONESTEP) then
      ! With the implicit Euler scheme, the pressure solution lives in the endpoints in
      ! time (of every time intervall spanned by subsequent time steps). So, velocity and
      ! pressure live in the same points in time such that we can just take over the
      ! solution.
      if (rtimestepping%dtheta .eq. 1.0_DP) then
        call lsysbl_copyVector (rvectorNew, rvectorInt)
        dtimeInt = dtimeNew

      ! With the general theta scheme, the pressure solution (being treated fully
      ! implicitly, as opposed to the velocity variable) lives between the points in time
      ! where the velocity solution is calculated. See e.g. page 750 of paper
      !    @article{Rang2008747,
      !       author  = "J. Rang",
      !       title   = "Pressure corrected implicit $\theta$-schemes for %!" fix compiler
      !                  the incompressible Navier--Stokes equations",    %!" warnings
      !       journal = "Applied Mathematics and Computation",
      !       volume  = "201",
      !       number  = "1--2",
      !       pages   = "747--761",
      !       year    = "2008",
      !       issn    = "0096-3003",
      !       doi     = "http://dx.doi.org/10.1016/j.amc.2008.01.010",
      !       url  = "http://www.sciencedirect.com/science/article/pii/S0096300308000428",
      !       note    = "",
      !    }
      ! for a proof. Take an appropriate mean in order to have both variables live at a
      ! common point in time (mandatory for subsequent calculations like lift and drag
      ! calculation) and interpolate *not* the pressure variable, but the velocity
      ! variable. Because only for the latter do we always have a start solution and as
      ! such can use it in time step 1 for the first interpolation step.
      else
        dfactor = rtimestepping%dtheta
        dtimeInt = (1.0_DP-dfactor)*dtimeOld + dfactor*dtimeNew

        call lsysbl_copyVector (rvectorNew, rvectorInt)
        call lsyssc_vectorLinearComb (rvectorOld%RvectorBlock(1),&
             rvectorInt%RvectorBlock(1),&
             1.0_DP - dfactor, dfactor)
        call lsyssc_vectorLinearComb (rvectorOld%RvectorBlock(2),&
             rvectorInt%RvectorBlock(2),&
             1.0_DP - dfactor, dfactor)
      end if

    else if (rtimestepping%ctimestepType .eq. TSCHM_FRACTIONALSTEP .or. &
             rtimestepping%ctimestepType .eq. TSCHM_FS_DIRK) then
      if (ipressureFullyImplicit .eq. 1) then
        ! For the fractional-step theta scheme, the paper
        !    @article{Rang2008747,
        !       author  = "J. Rang",
        !       title   = "Pressure corrected implicit $\theta$-schemes   %!" fix compiler
        !                  for the incompressible Navier--Stokes
        !                  equations",                                    %!" warnings
        !       journal = "Applied Mathematics and Computation",
        !       volume  = "201",
        !       number  = "1--2",
        !       pages   = "747--761",
        !       year    = "2008",
        !       issn    = "0096-3003",
        !       doi     = "http://dx.doi.org/10.1016/j.amc.2008.01.010",
        !       url ="http://www.sciencedirect.com/science/article/pii/S0096300308000428",
        !       note    = "",
        !    }
        ! proves on page 751 that the pressure solution (being treated fully implicitly,
        ! as opposed to the velocity variable) lives in three distinct points in time:
        !   substep 1:  t^{n} + alpha theta k                    = t^{n} + (4 theta - 1) k
        !   substep 2:  t^{n} + theta k + beta theta' k          = t^{n} + (5 theta - 1) k
        !   substep 3:  t^{n} + (theta+theta') k + alpha theta k = t^{n} + 3 theta k
        ! So, again pressure and velocity do definitely not live in the same point in
        ! time.  Moreover, the same paper proves that the fractional-step theta scheme
        ! suffers from order reduction for stiff ODEs (Lemma 5.2) and due to that the
        ! pressure is only approximated with first (!) order accuracy along with second
        ! order approximation of the velocity variable. That can also be easily determined
        ! experimantally (turn on ierrorAnalysisTimeSpace in postprocessing.dat and
        ! compare respective L2 errors for a sequence of time steps for a right hand side
        ! chosen according to an analytic, known solution).
        ! Again, take an appropriate mean in order to have both variables live at a common
        ! point in time (mandatory for subsequent calculations like lift and drag
        ! calculation) and interpolate *not* the pressure variable, but the velocity
        ! variable. Because only for the latter do we always have a start solution and as
        ! such can use it in time step 1 for the first interpolation step.

        dfactor = 0.0_DP
        ! explicit first stage: rtimestepping%isubstep == 2,3 or 4?
        ! implicit first stage: rtimestepping%isubstep == 1,2 or 3?
        select case (mod(rtimestepping%isubstep + 1 - &
                     merge(1,0,rtimestepping%bexplicitFirstStage .eqv. .TRUE.), &
                     rtimestepping%nsubsteps)+1)
          ! (Note: substep gets incremented *before* postprocessing starts. Even though we
          !        are still in the 3rd substep, the counter points already to the first
          !        substep of the next macro time step.)
        case (1)  ! 1st substep
          dfactor = 4.0_DP * rtimestepping%dtheta - 1.0_DP

        case (2)  ! 2nd substep
          dfactor = 5.0_DP * rtimestepping%dtheta - 1.0_DP

        case (3)  ! 3rd substep
          dfactor = 3.0_DP * rtimestepping%dtheta

        end select

        dtimeInt = (1.0_DP-dfactor)*dtimeOld + dfactor*dtimeNew

        call lsysbl_copyVector (rvectorNew, rvectorInt)
        call lsyssc_vectorLinearComb (rvectorOld%RvectorBlock(1),&
             rvectorInt%RvectorBlock(1),&
             1.0_DP - dfactor, dfactor)
        call lsyssc_vectorLinearComb (rvectorOld%RvectorBlock(2),&
             rvectorInt%RvectorBlock(2),&
             1.0_DP - dfactor, dfactor)

      else
        ! Velocity and pressure are treated in the same way, namely semi-implicitly. They
        ! should both live in the same point in time.
        call lsysbl_copyVector (rvectorNew, rvectorInt)
        dtimeInt = dtimeNew
      end if ! treat pressure implicitly or semi-implicitly?

    else if (rtimestepping%ctimestepType .eq. TSCHM_FS_GLOWINSKI) then
      ! For the time being, it is unknown at which point in time the pressure variable
      ! lives for this time stepping scheme. Do not interpolate just yet.
      call lsysbl_copyVector (rvectorNew, rvectorInt)
      dtimeInt = dtimeNew

    else if (rtimestepping%ctimestepType .eq. TSCHM_DIRK23L  .or. &
             rtimestepping%ctimestepType .eq. TSCHM_DIRK34La .or. &
             rtimestepping%ctimestepType .eq. TSCHM_DIRK34Lb .or. &
             rtimestepping%ctimestepType .eq. TSCHM_DIRK44L  .or. &
             rtimestepping%ctimestepType .eq. TSCHM_DIRK54L  .or. &
             rtimestepping%ctimestepType .eq. TSCHM_SDIRK2   .or. &
             rtimestepping%ctimestepType .eq. TSCHM_SDIRK3PR) then
      ! Velocity and pressure are treated in the same way by these two time stepping
      ! schemes, namely semi-implicitly. They should both live in the same point in time.
      call lsysbl_copyVector (rvectorNew, rvectorInt)
      dtimeInt = dtimeNew

    else
      call output_line ("unknown time stepping scheme with index" // &
                        trim(sys_siL(rtimestepping%ctimestepType, 2)), &
                        OU_CLASS_ERROR, OU_MODE_STD, "cc_interpolateTimesteps")
      call sys_halt()

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_solveNonstationary (rproblem,rvector,rrhs,rpostprocessing)

!<description>
  ! Starts the time discretisation. Proceeds in time until the final time
  ! or the maximum number of time steps is reached.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem

  ! The initial solution vector. Is replaced by the final solution vector.
  ! Must be unsorted.
  type(t_vectorBlock), intent(inout) :: rvector

  ! The initial RHS vector. Is replaced by the final RHS vector.
  ! Boundary conditions must not be implemented into that vector, this is done
  ! internally depending on the time.
  type(t_vectorBlock), intent(inout) :: rrhs

  ! Postprocessing structure. Defines what to do with solution vectors.
  type(t_c2d2postprocessing), intent(inout) :: rpostprocessing

!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j
    type(t_vectorBlock) :: rtempBlock1,rtempBlock2
    type(t_ccNonlinearIteration) :: rnonlinearIteration

    ! Configuration block of the time stepping scheme
    type(t_explicitTimeStepping)        :: rtimestepping,rtimesteppingPredictor

    ! The nonlinear solver configuration
    type(t_nlsolNode) :: rnlSol

    real(DP) :: dtimederivative,dtmp,dtmperror,dtimeratio,dcpuTime

    ! Time-interpolated solution vector
    type(t_vectorBlock) :: rvectorInt
    real(DP) :: dtimeInt

    ! Time error analysis and adaptive time stepping variables
    type(t_timestepSnapshot) :: rsnapshotLastMacrostep
    type(t_vectorBlock) :: rpredictedSolution
    type(t_vectorBlock) :: rprevTimeStepSolution, rprevMacroTimeStepSolution
    real(DP) :: dprevTimeStepTime, dprevMacroTimeStepTime

    logical :: babortTimestep
    integer :: ipressureFullyImplicit

    integer :: irepetition
    integer(I32) :: isolverStatus,isolverStatusPredictor
    type(t_timeError) :: rtimeError
    type(t_timeDerivatives) :: rtimeDerivative

    ! Timer for the current timestep and the total time.
    type(t_timer) :: rtimerTimestep,rtimerAllTimesteps

    ! Backup postprocessing structure
    type(t_c2d2postprocessing) :: rpostprocessingBackup

    integer :: cscheme
    real(DP) :: daux1, daux2, daux3


    ! Some preparations for the nonlinear solver.
    !
    ! Initialise the nonlinear solver node rnlSol with parameters from
    ! the INI/DAT files.
    call cc_getNonlinearSolver (rnlSol, rproblem%rparamList, "CC2D-NONLINEAR")

    ! Initialise the nonlinear loop. This is to prepare everything for
    ! our callback routines that are called from the nonlinear solver.
    ! The preconditioner in that structure is initialised later.
    call cc_initNonlinearLoop (rproblem,rproblem%NLMIN,rproblem%NLMAX,&
        rvector,rrhs,rnonlinearIteration,"CC2D-NONLINEAR")

    ! Initialise the time stepping scheme according to the problem configuration
    call cc_initTimeSteppingScheme (rproblem%rparamList,rtimestepping,&
        ipressureFullyImplicit)

    ! Initialise the preconditioner for the nonlinear iteration
    call cc_initPreconditioner (rproblem,&
        rnonlinearIteration,rvector,rrhs)

    ! Create temporary vectors we need for the nonlinear iteration.
    call lsysbl_createVector (rrhs, rtempBlock1, .false.)
    call lsysbl_createVector (rrhs, rtempBlock2, .false.)

    ! Initial time step
    rproblem%rtimedependence%itimeStep = 0
    rproblem%rtimedependence%dtime = rproblem%rtimedependence%dtimeInit
    dtimederivative = rproblem%rtimedependence%dminTimeDerivative

    ! Discretise the boundary conditions at the initial time.
    call cc_updateDiscreteBC (rproblem)

    ! Implement the initial boundary conditions into the solution vector.
    ! Do not implement anything to matrices or RHS vector as these are
    ! maintained in the timeloop.
    call cc_implementBC (rproblem,rvector,rrhs,.true.,.false.)

    ! Postprocessing. Write out the initial solution.
    call output_line ("Starting postprocessing of initial solution...")
    call cc_postprocessingNonstat (rproblem,&
        rvector,rproblem%rtimedependence%dtimeInit,&
        rvector,rproblem%rtimedependence%dtimeInit,&
        rvector,rproblem%rtimedependence%dtimeInit,&
        0,rpostprocessing,bisIntermediateStage=.FALSE.)

    ! Reset counter of current macro step repetitions.
    irepetition = 0

    !----------------------------------------------------
    ! Timeloop
    !----------------------------------------------------

    ! Snapshot the current solution for the later calculation of the
    ! time derivative.
    call lsysbl_copyVector (rvector, rprevMacroTimeStepSolution)
    dprevMacroTimeStepTime = rproblem%rtimedependence%dtime

    ! Start with the 1st timestep
    rproblem%rtimedependence%itimeStep = 1

    ! Start timers that calculate the current time.
    call stat_clearTimer(rtimerAllTimesteps)
    call stat_startTimer(rtimerAllTimesteps)

    do while ((rproblem%rtimedependence%itimeStep .le. &
               rproblem%rtimedependence%niterations) .and. &
              (rproblem%rtimedependence%dtime .lt. &
               rproblem%rtimedependence%dtimemax-100.0_DP*SYS_EPSREAL_DP) .and. &
              (dtimederivative .ge. &
               rproblem%rtimedependence%dminTimeDerivative) .or. &
              ! When using a multi-stage time stepping scheme (like classic
              ! fractional-step theta scheme or diagonally implicit Runge-Kutta scheme)
              ! do not exit before all internal stages have completed.
              ! (Note that for DIRK schemes with an explicit first stage the first stage
              ! really performed is stage 2.)
              (rtimestepping%isubstep .gt. &
               1 + merge(1,0,rtimestepping%bexplicitFirstStage .eqv. .TRUE.)))

      if (rproblem%rtimedependence%itimeStep .eq. 3) then
        ! Re-check time stepping scheme to use. Possibly we only needed a few steps with
        ! scheme A to generate a suitable start solution for scheme B.
        ! (Why the choice of 3? Why not re-check after 1 time step with scheme A? Because
        !  the scheme DIRK23L consists of 2 time (sub)steps and changing from scheme A to
        !  DIRK23L after 2 time steps with A is consistent with how we treat other DIRK
        !  schemes consisting of more than 2 substeps.)
        call parlst_getvalue_int (rproblem%rparamList, &
             "TIME-DISCRETISATION", "ITIMESTEPSCHEME", cscheme, 0)
        if (cscheme .eq. -TSCHM_DIRK23L) then
          ! Re-initialise the time stepping in the problem structure
          daux1 = rtimestepping%dcurrentTime
          daux2 = rtimestepping%dtstep
          daux3 = rtimestepping%dtheta
          call timstp_init (rtimestepping, abs(cscheme), daux1, daux2, daux3)

          ! Force-off flag for DIRK schemes: treat pressure semi-implicitly, as the
          ! velocity
          ipressureFullyImplicit = 0
        end if
      end if

      if (rproblem%rtimedependence%itimeStep .eq. 4) then
        ! Re-check time stepping scheme to use. Possibly we only needed a few steps with
        ! scheme A to generate a suitable start solution for scheme B.
        ! (Why the choice of 4? Why not re-check after 1 time step with scheme A? Because
        !  all the schemes B in question consist of 4 stages with the first one being
        !  explicit, i.e. 3 time (sub)steps and changing from scheme A to B after 3 time
        !  steps makes it possible to keep the algorithm simple that determines the time
        !  plus space discretisation error in cc_errorAnalysis(): the algorithm does not
        !  need to consider whether initially a different time stepping scheme (A) is
        !  used.)
        call parlst_getvalue_int (rproblem%rparamList, &
             "TIME-DISCRETISATION", "ITIMESTEPSCHEME", cscheme, 0)
        select case (-cscheme)
        case (TSCHM_FS_DIRK, &
              TSCHM_DIRK34La, &
              TSCHM_DIRK34Lb, &
              TSCHM_DIRK44L, &
              TSCHM_DIRK54L)
          ! Re-initialise the time stepping in the problem structure
          daux1 = rtimestepping%dcurrentTime
          daux2 = rtimestepping%dtstep
          daux3 = rtimestepping%dtheta
          call timstp_init (rtimestepping, abs(cscheme), daux1, daux2, daux3)

          ! Force-off flag for DIRK schemes: treat pressure semi-implicitly, as the
          ! velocity
          if (rtimestepping%ctimestepType .eq. TSCHM_FS_DIRK  .or. &
              rtimestepping%ctimestepType .eq. TSCHM_DIRK34La .or. &
              rtimestepping%ctimestepType .eq. TSCHM_DIRK34Lb .or. &
              rtimestepping%ctimestepType .eq. TSCHM_DIRK44L  .or. &
              rtimestepping%ctimestepType .eq. TSCHM_DIRK54L) then
            ipressureFullyImplicit = 0
          end if
        end select
      end if

      if (rproblem%rtimedependence%itimeStep .eq. 5) then
        ! Re-check time stepping scheme to use. Possibly we only needed a few steps with
        ! scheme A to generate a suitable start solution for scheme B.
        ! (Why the choice of 5? Why not re-check after 1 time step with scheme A? Because
        !  the scheme SDIRK2 consists of 4 stages (the first one being implicit, too),
        !  i.e. 4 time (sub)steps and changing from scheme A to SDIRK3PR after 5 time
        !  steps with A is consistent with how we treat other DIRK schemes consisting of
        !  less than 5 substeps.)
        call parlst_getvalue_int (rproblem%rparamList, &
             "TIME-DISCRETISATION", "ITIMESTEPSCHEME", cscheme, 0)
        if (cscheme .eq. -TSCHM_SDIRK2) then
          ! Re-initialise the time stepping in the problem structure
          daux1 = rtimestepping%dcurrentTime
          daux2 = rtimestepping%dtstep
          daux3 = rtimestepping%dtheta
          call timstp_init (rtimestepping, abs(cscheme), daux1, daux2, daux3)

          ! Force-off flag for DIRK schemes: treat pressure semi-implicitly, as the
          ! velocity
          ipressureFullyImplicit = 0
        end if
      end if

      if (rproblem%rtimedependence%itimeStep .eq. 6) then
        ! Re-check time stepping scheme to use. Possibly we only needed a few steps with
        ! scheme A to generate a suitable start solution for scheme B.
        ! (Why the choice of 6? Why not re-check after 1 time step with scheme A? Because
        !  the scheme SDIRK3PR consists of 5 stages (the first one being implicit, too),
        !  i.e. 5 time (sub)steps and changing from scheme A to SDIRK3PR after 5 time
        !  steps with A is consistent with how we treat other DIRK schemes consisting of
        !  less than 5 substeps.)
        call parlst_getvalue_int (rproblem%rparamList, &
             "TIME-DISCRETISATION", "ITIMESTEPSCHEME", cscheme, 0)
        if (cscheme .eq. -TSCHM_SDIRK3PR) then
          ! Re-initialise the time stepping in the problem structure
          daux1 = rtimestepping%dcurrentTime
          daux2 = rtimestepping%dtstep
          daux3 = rtimestepping%dtheta
          call timstp_init (rtimestepping, abs(cscheme), daux1, daux2, daux3)

          ! Force-off flag for DIRK schemes: treat pressure semi-implicitly, as the
          ! velocity
          ipressureFullyImplicit = 0
        end if
      end if

      ! Time counter
      call stat_clearTimer(rtimerTimestep)
      call stat_startTimer(rtimerTimestep)

      ! The babortTimestep is normally FALSE. If set to TRUE, the computation
      ! of the next time step is aborted because of an error.

      babortTimestep = .false.

      !----------------------------------------------------
      ! Predictor step
      !----------------------------------------------------
      !
      ! When adaptive time stepping is activated, we proceed as follows:
      ! 1.) Calculate one large time (macro-)step with step size 3*dtstepFixed
      ! 2.) Calculate three small time substeps with step size dtstepFixed
      ! 3.) Compare the solutions, calculate a new step size and/or repeat
      !     the time step.
      ! The "3*dtstepFixed" step size just "coincidentally" coincides with the
      ! step size of three steps in the Fractional Step Theta scheme :-)
      ! So after each three substeps, the simulation time of the small
      ! substeps will always coincide with the simulation time after the
      ! large macrostep
      !
      ! Do we use adaptive time stepping?
      select case (rproblem%rtimedependence%radaptiveTimeStepping%ctype)
      case (TADTS_USERDEF)
        ! Nothing to be done

      case (TADTS_FIXED)
        ! Nothing to be done

      case (TADTS_PREDICTION,TADTS_PREDICTREPEAT,TADTS_PREDREPTIMECONTROL)
        ! Adaptive time stepping. Is this the first substep?
        if (mod(rproblem%rtimedependence%itimeStep,3) .eq. 1) then

          ! If this is not a repetition of a (macro-) timestep,
          ! create a backup of the current flow situation, so we can repeat the time
          ! step if necessary.
          if (irepetition .eq. 0) then
            call cc_backupTimestep (rsnapshotLastMacrostep,rproblem,rtimeStepping,&
                rvector,rrhs)
            call cc_copyPostprocessing (rpostprocessing,rpostprocessingBackup)
          end if

          ! At first, perform one predictor step with 3x time step size
          ! using implicit Euler.
          ! For this purpose, create a new time stepping structure
          ! based on the standard one (for small time steps).
          call timstp_init (rtimesteppingPredictor, TSCHM_ONESTEP, &
                rproblem%rtimedependence%dtime, &
                3.0_DP*rtimestepping%dtstepFixed, 1.0_DP)

          call output_lbrk (coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
          call output_separator(OU_SEP_EQUAL,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
          call output_line ("Predictor-Step at Time-Step "// &
              trim(sys_siL(rproblem%rtimedependence%itimeStep,6))// &
              ", Time = "// &
              trim(sys_sdL(rproblem%rtimedependence%dtime,5))// &
              ", Stepsize: DT1 = "// &
              trim(sys_sdL(rtimesteppingPredictor%dtstep,5)),&
              coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
          call output_separator(OU_SEP_EQUAL,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
          !CALL output_line (&
          !  "Macro step "//TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6))// &
          !  "      Time "//TRIM(sys_sdL(rproblem%rtimedependence%dtime,5)) // &
          !  "      Step size "//TRIM(sys_sdL(rtimesteppingPredictor%dtstep,5)))
          !CALL output_lbrk ()

          ! Proceed in time, calculate the predicted solution.
          call lsysbl_copyVector (rvector,rpredictedSolution)
          call cc_performTimestep (rproblem,rpredictedSolution,rrhs,&
              rtimesteppingPredictor,ipressureFullyImplicit,rnonlinearIteration,rnlSol,&
              rtempBlock1,rtempBlock2)

          ! Gather statistics
          rproblem%rstatistics%ntimesteps = rproblem%rstatistics%ntimesteps + 1

          ! Remember the status of the nonlinear solver for later use.
          isolverStatusPredictor = rnlSol%iresult

          ! Did the nonlinear solver break down?
          if (rnlSol%iresult .gt. 0) then
            ! Oops, not really good.
            babortTimestep = .true.

            ! Calculate a new time step size.
            ! Set bit 0 and not bit 2/3 in isolverStatus as we want to compute the
            ! new time step based on the solver status of the last computation and
            ! not on a combined analysis of predictor step and actual solution!
            isolverStatus = 0
            select case (rnlSol%iresult)
            case (:-1)
              ! Nonlinear solver worked, but could not reach convergence criterion
              isolverStatus = ior(isolverStatus,TADTS_SST_NLINCOMPLETE)
            case (0)
              ! Everything fine
            case (1)
              ! Nonlinear solver diverged
              isolverStatus =  ior(isolverStatus,TADTS_SST_NLFAIL)
            case (2)
              ! Nonlinear solver diverged because of error in the preconditioner
              isolverStatus =  ior(isolverStatus,&
                                   TADTS_SST_NLFAIL + TADTS_SST_NLPRECFAIL)
            case (3)
              ! General error
              isolverStatus = not(0)
            end select
            dtmp = adtstp_calcTimeStep (&
                rproblem%rtimedependence%radaptiveTimeStepping, &
                0.0_DP, &
                rproblem%rtimedependence%dtimeInit,&
                rproblem%rtimedependence%dtime, &
                rtimeStepping%dtstepFixed, &
                timstp_getOrder(rtimeStepping), &
                isolverStatus,irepetition)

            ! Tell the user that we have a new time step size.
            call output_separator(OU_SEP_AT,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
            call output_line ("Timestepping by "&
                //trim(sys_siL(irepetition,2)) &
                //" (" &
                //trim(sys_siL(rproblem%rtimedependence%itimeStep,6)) &
                //"), New Stepsize = " &
                //trim(sys_sdEP(dtmp,9,2)) &
                //", Old Stepsize = " &
                //trim(sys_sdEP(rtimeStepping%dtstepFixed,9,2)),&
                coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
            call output_separator(OU_SEP_AT,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)

            ! Accept the new step size
            call timstp_setBaseSteplength (rtimeStepping, dtmp)

            ! The old RHS is restored later...

          else

            ! Restore the flow situation to the beginning of the macrostep.
            ! We did not change the solution vector and the time stepping structure,
            ! so we do not have to restore them.
            call cc_restoreTimestep (rsnapshotLastMacrostep,rproblem,rrhs=rrhs)
            call cc_copyPostprocessing (rpostprocessingBackup,rpostprocessing)

            ! The solver worked, rpredictedSolution contains the predicted
            ! solution. Now we can continue with the usual time stepping.

          end if

        end if
      end select

      !----------------------------------------------------
      ! Proceed one step in time
      !----------------------------------------------------
      if (.not. babortTimestep) then

        !CALL output_separator(OU_SEP_MINUS)
        !CALL output_line (&
        !  "Time step "//TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6))// &
        !  "     Time "//TRIM(sys_sdL(rproblem%rtimedependence%dtime,5)) // &
        !  "     Step size "//TRIM(sys_sdL(rtimestepping%dtstep,5)))
        call output_lbrk (coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
        call output_separator(OU_SEP_AT,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
        if (rtimestepping%nsubsteps .gt. 1) then
          call output_line ("Time-Step "// &
               trim(sys_siL(rproblem%rtimedependence%itimeStep,6))// &
               ", Time = "// &
               trim(sys_sdL(rproblem%rtimedependence%dtime,5))// &
               ", Stepsize: DT3 = "// &
               trim(sys_sdL(rtimestepping%dtstep,5)) // &
               " (stage " // &
               trim(sys_siL(rtimestepping%isubstep,1)) // &
               ")", &
               coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
        else
          call output_line ("Time-Step "// &
               trim(sys_siL(rproblem%rtimedependence%itimeStep,6))// &
               ", Time = "// &
               trim(sys_sdL(rproblem%rtimedependence%dtime,5))// &
               ", Stepsize: DT3 = "// &
               trim(sys_sdL(rtimestepping%dtstep,5)),&
               coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
        end if
        call output_separator(OU_SEP_AT,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)

        ! Snapshot the current solution for the later calculation of an interpolated
        ! solution vector where the time of pressure and velocity matches.
        ! Given that the results in intermediate stages of most DIRK schemes have no
        ! physical meaning, avoid to use a snapshot of these when trying to determine
        ! (later on) whether a transient solution has reached a steady state.
        !
        ! explicit first stage: rtimestepping%isubstep == 2?
        ! implicit first stage: rtimestepping%isubstep == 1?
        if (rtimestepping%isubstep - &
            merge(1,0,rtimestepping%bexplicitFirstStage .eqv. .TRUE.) .eq. 1) then
          call lsysbl_copyVector (rvector,rprevTimeStepSolution)
          dprevTimeStepTime = rproblem%rtimedependence%dtime
        end if

        ! Proceed to next time step -- if we are allowed to.
        call cc_performTimestep (rproblem,rvector,rrhs,&
            rtimestepping,ipressureFullyImplicit,rnonlinearIteration,rnlSol,&
            rtempBlock1,rtempBlock2)

        ! Gather statistics
        rproblem%rstatistics%ntimesteps = rproblem%rstatistics%ntimesteps + 1

        ! Do we count in steps a 1 or in steps a 3?
        ! Respecting this, i is assigned the number of the substep in the
        ! macrostep.
        select case (rproblem%rtimedependence%radaptiveTimeStepping%ctype)
        case (TADTS_FIXED,TADTS_USERDEF)
          i = 1
          j = 1
        case (TADTS_PREDICTION,TADTS_PREDICTREPEAT,TADTS_PREDREPTIMECONTROL)
          i = mod(rproblem%rtimedependence%itimeStep-1,3)+1
          j = 3
        end select

        call output_separator(OU_SEP_AT,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
        call output_line ("Time-Step " &
            //trim(sys_siL(rproblem%rtimedependence%itimeStep,6)) &
            //" (Repetition = "//trim(sys_siL(irepetition,2)) &
            //", Substep = " &
            //trim(sys_siL(i,6))//" of "//trim(sys_siL(j,6)) &
            //") at time = " &
            //trim(sys_sdL(rproblem%rtimedependence%dtime,5)) &
            //" finished. ",coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )

        if (rproblem%rtimedependence%radaptiveTimeStepping%ctype .eq. TADTS_USERDEF) then

          ! User defined time stepping.
          !
          ! Calculate a new time step size.
          isolverStatus = 0
          select case (rnlSol%iresult)
          case (:-1)
            ! Nonlinear solver worked, but could not reach convergence criterion
            isolverStatus = ior(isolverStatus,TADTS_SST_NLINCOMPLETE)
          case (0)
            ! Everything fine
          case (1)
            ! Nonlinear solver diverged
            isolverStatus =  ior(isolverStatus,TADTS_SST_NLFAIL)
          case (2)
            ! Nonlinear solver diverged because of error in the preconditioner
            isolverStatus =  ior(isolverStatus,&
                                TADTS_SST_NLFAIL + TADTS_SST_NLPRECFAIL)
          case (3)
            ! General error
            isolverStatus = not(0)
          end select
          call cc_initCollectForAssembly (rproblem,rproblem%rcollection)
          dtmp = adtstp_calcTimeStep (&
              rproblem%rtimedependence%radaptiveTimeStepping, &
              0.0_DP, &
              rproblem%rtimedependence%dtimeInit,&
              rproblem%rtimedependence%dtime, &
              rtimeStepping%dtstepFixed, &
              timstp_getOrder(rtimeStepping), &
              isolverStatus,irepetition,calcAdaptiveTimestep,rproblem%rcollection)
          call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

          ! Tell the user that we have a new time step size.
          call output_line ("New Stepsize = " &
              //trim(sys_sdEP(dtmp,9,2)) &
              //", Old Stepsize = " &
              //trim(sys_sdEP(rtimeStepping%dtstepFixed,9,2)),&
              coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
          call output_separator(OU_SEP_AT,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)

          ! Accept the new step size
          call timstp_setBaseSteplength (rtimeStepping, dtmp)

        else

          ! Did the solver break down?
          if (rnlSol%iresult .lt. 0) then
            call output_line ("Accuracy notice: Nonlinear solver did not reach "// &
                              "the convergence criterion!",&
                              coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
          else if (rnlSol%iresult .gt. 0) then
            ! Oops, not really good.
            babortTimestep = .true.

            ! Do we have a time stepping algorithm that allows recomputation?
            select case (rproblem%rtimedependence%radaptiveTimeStepping%ctype)
            case (TADTS_FIXED,TADTS_PREDICTION)
              ! That is bad. Our solution is garbage!
              ! We do not do anything in this case. The repetition technique will
              ! later decide on whether to repeat the step or to stop the
              ! computation.
              call output_line("Nonlinear solver broke down. Solution probably garbage!",&
                  coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)

            case (TADTS_PREDICTREPEAT,TADTS_PREDREPTIMECONTROL)
              ! Yes, we have.
              call output_line ("Nonlinear solver broke down. "// &
                                "Calculating new time step size...",&
                                coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)

              ! Calculate a new time step size.
              isolverStatus = 0
              select case (rnlSol%iresult)
              case (:-1)
                ! Nonlinear solver worked, but could not reach convergence criterion
                isolverStatus = ior(isolverStatus,TADTS_SST_NLINCOMPLETE)
              case (0)
                ! Everything fine
              case (1)
                ! Nonlinear solver diverged
                isolverStatus =  ior(isolverStatus,TADTS_SST_NLFAIL)
              case (2)
                ! Nonlinear solver diverged because of error in the preconditioner
                isolverStatus =  ior(isolverStatus,&
                                    TADTS_SST_NLFAIL + TADTS_SST_NLPRECFAIL)
              case (3)
                ! General error
                isolverStatus = not(0)
              end select
              dtmp = adtstp_calcTimeStep (&
                  rproblem%rtimedependence%radaptiveTimeStepping, &
                  0.0_DP, &
                  rproblem%rtimedependence%dtimeInit,&
                  rproblem%rtimedependence%dtime, &
                  rtimeStepping%dtstepFixed, &
                  timstp_getOrder(rtimeStepping), &
                  isolverStatus,irepetition)

              ! Tell the user that we have a new time step size.
              !CALL output_line ("Timestepping by "&
              !    //TRIM(sys_siL(irepetition,2)) &
              !    //" (" &
              !    //TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6)) &
              !    //"), New Stepsize = " &
              !    //TRIM(sys_sdEP(dtmp,9,2)) &
              !    //", Old Stepsize = " &
              !    //TRIM(sys_sdEP(rtimeStepping%dtstepFixed,9,2)) )
              call output_line ("New Stepsize = " &
                  //trim(sys_sdEP(dtmp,9,2)) &
                  //", Old Stepsize = " &
                  //trim(sys_sdEP(rtimeStepping%dtstepFixed,9,2)),&
                  coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
              call output_separator(OU_SEP_AT,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)

              ! Accept the new step size
              call timstp_setBaseSteplength (rtimeStepping, dtmp)

            end select

          end if

        end if

      end if

      !----------------------------------------------------
      ! Time error control
      !----------------------------------------------------
      if (.not. babortTimestep) then

        ! Ok, everything worked fine, we have a valid solution!
        !
        ! Time step control. Do we have a time stepping algorithm that
        ! adapts the time step size and probably wants to repeat the
        ! calculation?
        select case (rproblem%rtimedependence%radaptiveTimeStepping%ctype)
        case (TADTS_USERDEF)
          ! No, continue as usual.

        case (TADTS_FIXED)
          ! No, continue as usual.

        case (TADTS_PREDICTION,TADTS_PREDICTREPEAT,TADTS_PREDREPTIMECONTROL)

          ! Yes. Time step adaption takes place after 3 steps, where we reached
          ! the same point in time like the big macrostep.

          if (mod(rproblem%rtimedependence%itimeStep,3) .eq. 0) then

            call output_separator (OU_SEP_MINUS,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
            call output_line ("Macrostep completed. Analysing time error and "// &
                              "computing new time step size...",&
                              coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)

            ! Calculate the new time step size.
            ! This is based on the solution, the predicted solution, the order
            ! of the time stepping algorithm and on the status of the solvers.
            !
            ! At first, calculate the time error.
            dtmperror =  cc_timeErrorByPredictor (&
                rproblem%rtimedependence%radaptiveTimeStepping%cadTimeStepErrorControl,&
                rvector,rpredictedSolution,rtempBlock1,rtimeError)

            ! Evaluate everything that went wrong in the solvers
            isolverStatus = 0
            select case (rnlSol%iresult)
            case (:-1)
              ! Nonlinear solver worked, but could not reach convergence criterion
              isolverStatus = ior(isolverStatus,TADTS_SST_NLINCOMPLETE)
            case (0)
              ! Everything fine
            case (1)
              ! Nonlinear solver diverged
              isolverStatus =  ior(isolverStatus,TADTS_SST_NLFAIL)
            case (2)
              ! Nonlinear solver diverged because of error in the preconditioner
              isolverStatus =  ior(isolverStatus,&
                                   TADTS_SST_NLFAIL + TADTS_SST_NLPRECFAIL)
            case (3)
              ! General error
              isolverStatus = not(0)
            end select

            select case (isolverStatusPredictor)
            case (:-1)
              ! Nonlinear solver worked, but could not reach convergence criterion
              isolverStatus = ior(isolverStatus,TADTS_SST_NLPREDINCOMPLETE)
            case (0)
              ! Everything fine
            case (1)
              ! Nonlinear solver diverged
              isolverStatus =  ior(isolverStatus,TADTS_SST_NLPREDFAIL)
            case (2)
              ! Nonlinear solver diverged because of error in the preconditioner
              isolverStatus =  ior(isolverStatus,&
                                   TADTS_SST_NLPREDFAIL + TADTS_SST_NLPREDPRECFAIL)
            case (3)
              ! General error
              isolverStatus = not(0)
            end select

            ! Calculate the new time step size
            dtmp = adtstp_calcTimeStep (&
                rproblem%rtimedependence%radaptiveTimeStepping, dtmperror, &
                rproblem%rtimedependence%dtimeInit,&
                rproblem%rtimedependence%dtime, &
                rtimeStepping%dtstepFixed, &
                timstp_getOrder(rtimeStepping), &
                isolverStatus,irepetition)

            ! Calculate the relation of the previous and new step size
            dtimeratio = dtmp / rtimeStepping%dtstepFixed

            ! Forces us the time stepping algorithm to recompute?
            if (rproblem%rtimedependence%radaptiveTimeStepping%ctype .eq. &
                TADTS_PREDREPTIMECONTROL) then

              ! When the new time step is much smaller than the old one,
              ! set babortTimestep to TRUE to indicate that
              ! the time-step has to be repeated.

              if (dtimeratio .lt. rproblem%rtimedependence%radaptiveTimeStepping% &
                                  depsAdaptiveRelTimeStep) then
                babortTimestep = .true.
              end if

            end if

            ! Assign i either 2 or 8, depending on whether we have a time stepping
            ! algorithm of order one or two in time!
            i = 2 + 6*(timstp_getOrder(rtimeStepping)-1)

            ! Print the result of the time error analysis
            call output_line ("Time error: " &
                //" U(L2)=" &
                //trim(sys_sdEP(rtimeError%drelUL2/real(i,DP),8,2)) &
                //"  U(MX)=" &
                //trim(sys_sdEP(rtimeError%drelUmax/real(i,DP),8,2)) &
                //"  P(L2)=" &
                //trim(sys_sdEP(rtimeError%drelPL2/real(i,DP),8,2)) &
                //"  P(MX)=" &
                //trim(sys_sdEP(rtimeError%drelPmax/real(i,DP),8,2)),&
                coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )

            ! Tell the user that we have a new time step size.
            call output_line ("New Stepsize = " &
                //trim(sys_sdEP(dtmp,9,2)) &
                //", Old Stepsize = " &
                //trim(sys_sdEP(rtimeStepping%dtstepFixed,9,2)),&
                coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
            !CALL output_line ("Timestepping by "&
            !    //TRIM(sys_siL(irepetition,2)) &
            !    //" (" &
            !    //TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6)) &
            !    //"), New Stepsize = " &
            !    //TRIM(sys_sdEP(dtmp,9,2)) &
            !    //", Old Stepsize = " &
            !    //TRIM(sys_sdEP(rtimeStepping%dtstepFixed,9,2)) )

            ! Accept the new step size
            call timstp_setBaseSteplength (rtimeStepping, dtmp)

          end if

        end select
      end if

      !----------------------------------------------------
      ! Postprocessing
      !----------------------------------------------------
      if (.not. babortTimestep) then

        ! Ok, everything worked fine, we have a valid solution of
        ! our current substep.

        call output_separator(OU_SEP_MINUS,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
        call output_line ("Starting postprocessing of the time step...",&
            coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)

        ! Calculate an interpolated solution vector where the time of
        ! pressure and velocity matches.
        call cc_interpolateTimesteps (rtimestepping, &
             rprevTimeStepSolution, dprevTimeStepTime,&
             rvector, rproblem%rtimedependence%dtime, rvectorInt, dtimeInt, &
             ipressureFullyImplicit)

        ! Postprocessing. Write out the solution if it was calculated successfully and
        ! measure time errors for test problems (with analytically given solution).
        rpostprocessing%ctimestepType = rtimestepping%ctimestepType
        call cc_postprocessingNonstat (rproblem, &
            rprevTimeStepSolution, dprevTimeStepTime, &
            rvector, rproblem%rtimedependence%dtime, &
            rvectorInt, dtimeInt, &
            rproblem%rtimedependence%itimeStep, rpostprocessing, &
            ! (Note: substep gets incremented *before* postprocessing starts. So, even
            !        though we are, e.g., still in the 3rd substep, the counter does
            !        already point to the first substep of the next macro time step.)
            ! explicit first stage: rtimestepping%isubstep > 2?
            ! implicit first stage: rtimestepping%isubstep > 1?
            bisIntermediateStage=&
               (rtimestepping%isubstep - &
                merge(1,0,rtimestepping%bexplicitFirstStage .eqv. .TRUE.) .gt. 1))

        call lsysbl_releaseVector (rvectorInt)

        ! Calculate time derivative, ignoring results from intermediate stages in a
        ! multi-stage time stepping scheme (FS, DIRK). This is actually non-optional as
        ! some DIRK schemes (DIRK34La, DIRK34Lb, DIRK44L) have multiple stages where the
        ! solution lives in the same point in time such that naive calculation of the time
        ! derivate fails due to coinciding interpolation points.
        ! (Note: substep gets incremented *before* postprocessing starts. So, even though
        !        we are, e.g., still in the 3rd substep, the counter does already point to
        !        the first substep of the next macro time step.)
        !
        ! explicit first stage: rtimestepping%isubstep == 2?
        ! implicit first stage: rtimestepping%isubstep == 1?
        if (rtimestepping%isubstep - &
            merge(1,0,rtimestepping%bexplicitFirstStage .eqv. .TRUE.) .eq. 1) then
          call output_separator(OU_SEP_MINUS,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
          call output_line ("Analysing time derivative...",&
               coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)

          ! Calculate the norm of the time derivative. This allows the DO-loop
          ! above to check if the solution got stationary.
          dtimeDerivative = cc_timeDerivative (&
               rproblem%rtimedependence%radaptiveTimeStepping%cadTimeStepErrorControl, &
               rvector, rprevTimeStepSolution, rtimestepping%dtstep, rtempBlock1, &
               rtimeDerivative)

          ! Print the results of the time analysis.

          call output_line ("Time derivative:  " &
               //" RELU(L2)=" &
               //trim(sys_sdEP(rtimeDerivative%drelUL2,9,2)) &
               //"  RELP(L2)=" &
               //trim(sys_sdEP(rtimeDerivative%drelPL2,9,2)) &
               //"  REL=" &
               //trim(sys_sdEP(dtimeDerivative,9,2)),&
               coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )

          if (dtimederivative .lt. rproblem%rtimedependence%dminTimeDerivative) then
            call output_line ("Solution reached stationary status. " // &
                 "Stopping simulation...", coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
          end if
        end if

        !CALL output_line ("#"&
        !    //TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6))&
        !    //" (" &
        !    //TRIM(sys_siL(i,6)) &
        !    //") TIME=" &
        !    //TRIM(sys_sdL(rproblem%rtimedependence%dtime,5)) &
        !    //" RELU(L2)=" &
        !    //TRIM(sys_sdEP(rtimeDerivative%drelUL2,9,2)) &
        !    //" RELP(L2)=" &
        !    //TRIM(sys_sdEP(rtimeDerivative%drelPL2,9,2)) &
        !    //" REL=" &
        !    //TRIM(sys_sdEP(dtimeDerivative,9,2)) )

      end if

      call output_separator(OU_SEP_MINUS,coutputMode=OU_MODE_STD)
      call stat_stopTimer(rtimerTimestep)
      call output_line ("Time for processing of this timestep: "// &
        trim(sys_sdL(rtimerTimestep%delapsedReal,10)))

      call stat_sampleTimer(rtimerAllTimesteps,dcpuTime)
      call output_line ("Time for all timesteps until now:     "// &
        trim(sys_sdL(dcpuTime,10)))

      !----------------------------------------------------
      ! Check if the time step must be repeated
      !----------------------------------------------------
      if (babortTimestep) then

        ! Uh, oh, something went wrong. Let us hope we can repeat the timestep!
        !
        ! We have to repeat the time step if
        !  a) we have repetitions left and
        !  b) The user has activated any kind of adaptive time stepping and
        !  c1) the parameters force us to repeat a step or
        !  c2) the solver for the nonlinear iteration (predictor
        !      or corrector step) broke down
        ! c2) is an exception to the rule that we do not repeat
        !     time steps in adaptive time stepping technique 1!

        ! Do we have a time stepping algorithm that allows recomputation?
        select case (rproblem%rtimedependence%radaptiveTimeStepping%ctype)
        case (TADTS_FIXED,TADTS_USERDEF)
          ! That is bad. Our solution is most probably garbage!
          ! We cancel the timeloop, it does not make any sense to continue.
          call output_line ("Solution garbage! Stopping simulation.",&
              coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
          call output_separator(OU_SEP_AT,&
              coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
          exit

        case (TADTS_PREDICTION,TADTS_PREDICTREPEAT,TADTS_PREDREPTIMECONTROL)
          ! So far, so good. Are we allowed to recompute the time step?
          if (irepetition .lt. &
              rproblem%rtimedependence%radaptiveTimeStepping%nrepetitions) then

            ! Lucky, we are allowed to recompute :-)
            ! The previous situation was restored before -- so increase the
            ! repetition counter and cycle the loop without increasing
            ! the number of the current time step.

            irepetition = irepetition + 1

            ! Restore the flow situation to the beginning of the macrostep.
            ! Use the newly calculated time step size for the next time step.
            dtmp = rtimeStepping%dtstepFixed
            call cc_restoreTimestep (rsnapshotLastMacrostep,rproblem,&
                rtimeStepping,rvector,rrhs)
            call cc_copyPostprocessing (rpostprocessingBackup,rpostprocessing)
            call timstp_setBaseSteplength (rtimeStepping, dtmp)

            call output_line ("Repeating macrostep. Returning to timestep " &
                //trim(sys_siL(rproblem%rtimedependence%itimeStep,6))//".",&
                coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
            call output_separator(OU_SEP_AT,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)

            ! Repeat the time step
            cycle

          else

            ! No repetitions left. Reset the repetition counter and
            ! continue with the next timestep, although the
            ! solution might be not really meaningful anymore.
            irepetition = 0

            call output_line ("No repetitions left. Cannot repeat macrostep. " &
                //"Continuing with the next one...",&
                coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)

          end if

        end select

      else  ! IF (babortTimestep) THEN

        ! Snapshot the current solution for the later calculation of the time derivative
        ! (for which intermediate stages in a multi-stage time stepping scheme (FS, DIRK)
        ! are to be ignored).
        ! (Note: substep gets incremented *before* postprocessing starts. So, even though
        !        we are, e.g., still in the 3rd substep, the counter does already point to
        !        the first substep of the next macro time step.)
        !
        ! explicit first stage: rtimestepping%isubstep == 2?
        ! implicit first stage: rtimestepping%isubstep == 1?
        if (rtimestepping%isubstep - &
            merge(1,0,rtimestepping%bexplicitFirstStage .eqv. .TRUE.) .eq. 1) then
          call lsysbl_copyVector (rvector, rprevMacroTimeStepSolution)
          dprevMacroTimeStepTime = rproblem%rtimedependence%dtime
        end if

        ! No, time step is ok. If this is the last time step of the macro step,
        ! reset the repetition counter such that it starts with 0 for the
        ! next macrostep.
        if (mod(rproblem%rtimedependence%itimeStep, rtimestepping%nsubsteps) .eq. 0) then
          irepetition = 0
        end if

      end if

      rproblem%rtimedependence%itimeStep = rproblem%rtimedependence%itimeStep + 1

      call output_separator(OU_SEP_AT,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)

    end do

    ! Clean up the stuff of/for the nonlinear solver.
    !
    ! Release the temporary vectors
    if (rpredictedSolution%NEQ .ne. 0) &
      call lsysbl_releaseVector (rpredictedSolution)
    if (rprevTimeStepSolution%NEQ .ne. 0) &
      call lsysbl_releaseVector (rprevTimeStepSolution)
    if (rprevMacroTimeStepSolution%NEQ .ne. 0) &
      call lsysbl_releaseVector (rprevMacroTimeStepSolution)
    call lsysbl_releaseVector (rtempBlock2)
    call lsysbl_releaseVector (rtempBlock1)

    ! Release vectors (possibly) allocated in cc_errorAnalysis
    if (rpostprocessing%icalcTimeSpaceDiscrErrors .ne. 0) then
      if (rpostprocessing%itimeSpaceDiscrErrorMethod .ge. 4) then
        call lsysbl_releaseVector (rpostprocessing%roldestSolution)
      end if
      if (rpostprocessing%itimeSpaceDiscrErrorMethod .ge. 2) then
        call lsysbl_releaseVector (rpostprocessing%rolderSolution)
        call lsysbl_releaseVector (rpostprocessing%roldSolution)
      end if
    end if

    ! Release existing snapshots
    call cc_releaseSnapshot (rsnapshotLastMacrostep)

    ! Release the preconditioner
    call cc_releasePreconditioner (rnonlinearIteration)

    ! Release parameters of the nonlinear loop, final clean up
    call cc_doneNonlinearLoop (rnonlinearIteration)

  end subroutine

end module
