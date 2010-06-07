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
    call parlst_getvalue_int (rparams,ssection,'itimedependence',  &
        rproblem%itimedependence, 0)
    call parlst_getvalue_int (rparams,ssection,'niterations',     &
        rproblem%rtimedependence%niterations, 1000)
    call parlst_getvalue_double (rparams,ssection,'dtimeInit',   &
        rproblem%rtimedependence%dtimeInit, 0.0_DP)
    call parlst_getvalue_double (rparams,ssection,'dtimemax',     &
        rproblem%rtimedependence%dtimeMax, 20.0_DP)
    call parlst_getvalue_double (rparams,ssection,'dtimeStep',    &
        rproblem%rtimedependence%dtimeStep, 0.01_DP)
    call parlst_getvalue_double (rparams,ssection,'dminTimeDerivative',    &
        rproblem%rtimedependence%dminTimeDerivative, 0.00001_DP)
    rproblem%rtimedependence%itimeStep = 0

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
  ! A snapshop structure that receives the snapshot.
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
  ! A snapshop structure that receives the snapshot.
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

  subroutine cc_releaseSnapshop (rsnapshot)

!<description>
  ! Releases the rsnapshot structure. All allocated memory is released.
!</description>

!<inputoutput>
  ! A snapshop structure that receives the snapshot.
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

  subroutine cc_initTimeSteppingScheme (rparams,rtimeStepping)
  
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
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: cscheme
    real(DP) :: dtheta,dtstep,dtimemin

    ! Get the parameters for the time stepping scheme from the parameter list
    call parlst_getvalue_int (rparams, &
        'TIME-DISCRETISATION', 'ITIMESTEPSCHEME', cscheme, 0)
    call parlst_getvalue_double (rparams, &
        'TIME-DISCRETISATION', 'DTIMESTEPTHETA', dtheta, 1.0_DP)
    call parlst_getvalue_double (rparams, &
        'TIME-DISCRETISATION', 'DTIMESTEP', dtstep, 0.1_DP)
    call parlst_getvalue_double (rparams, &
        'TIME-DISCRETISATION', 'DTIMEINIT', dtimemin, 0.0_DP)
    
    ! Initialise the time stepping in the problem structure
    call timstp_init (rtimestepping, &
                      cscheme, dtimemin, dtstep, dtheta)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine cc_performTimestep (rproblem,rvector,rrhs,&
      rtimestepping,rnonlinearIteration,rnlSolver,rtempVector,rtempVectorRhs)

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

  ! Structure for the nonlinear iteration for solving the core equation.
  type(t_ccnonlinearIteration), intent(inout) :: rnonlinearIteration
  
  ! A configuration stucture for the nonlinear solver
  type(t_nlsolNode) :: rnlSolver
  
  ! Temporary vectors for the nonlinear solver
  type(t_vectorBlock), intent(inout) :: rtempVector,rtempVectorRhs
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_ccnonlinearIteration) :: rnonlinearIterationTmp
    type(t_nonlinearCCMatrix) :: rnonlinearCCMatrix
    
    ! DEBUG!!!
    !REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2
  
    ! The new RHS will be set up in rtempVectorRhs. Assign the discretisation/
    ! boundary conditions of rrhs to that vector so that rtempVectorRhs
    ! acts as a RHS vector.
    call lsysbl_assignDiscrIndirect(rrhs,rtempVectorRhs)
    
    ! DEBUG!!!
    !CALL lsysbl_getbase_double (rvector,p_Ddata)
    !CALL lsysbl_getbase_double (rtempVectorRhs,p_Ddata2)
  
    ! We have an equation of the type
    !
    !   d/dt u(x,t)  +  N(u(x,t))  =  f(x,t)
    !
    ! Which is discretised in time with a Theta scheme, leading to
    !
    !   $$ u_{n+1} + w_1*N(u_n+1) 
    !      =   u_n + w_2*N(u_n)  +  w_3*f_{n+1}  +  w_4*f_n $$
    !
    ! with k=time step size, u_{n+1} = u(.,t_{n+1}),etc., c.f. timestepping.f90.
    !
    ! The RHS of that equation therefore contains parts of the solution
    ! u_n, of the old RHS f_n and the new RHS f_{n+1}. At first, we make
    ! a weighted copy of the current RHS f_n to the 'global' RHS vector
    ! according to the time stepping scheme.
    
    ! Set up w_4*f_n.
    call lsysbl_vectorLinearComb(rrhs,rtempVectorRhs,&
         rtimestepping%dweightOldRHS,0.0_DP)
    
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
        rproblem%RlevelInfo(rproblem%NLMAX)%rstaticInfo,&
        rproblem%RlevelInfo(rproblem%NLMAX)%rdynamicInfo)
    
    rnonlinearCCMatrix%dalpha = -1.0_DP
    rnonlinearCCMatrix%dtheta = -rtimestepping%dweightMatrixRHS
    rnonlinearCCMatrix%dgamma = -rtimestepping%dweightMatrixRHS * real(1-rproblem%iequation,DP)
    rnonlinearCCMatrix%deta = 0.0_DP
    rnonlinearCCMatrix%dtau = 0.0_DP
    
    call cc_nonlinearMatMul (rnonlinearCCMatrix,rvector,rtempVectorRhs,-1.0_DP,1.0_DP)

    ! -------------------------------------------    
    ! Switch to the next point in time.
    rproblem%rtimedependence%dtime = rtimestepping%dcurrenttime + rtimestepping%dtstep
          
    ! Discretise the boundary conditions at the new point in time -- 
    ! if the boundary conditions are nonconstant in time!
    if (rproblem%iboundary .ne. 0) then
      call cc_updateDiscreteBC (rproblem)
    end if

    ! -------------------------------------------    

    ! generate f_n+1 into the rrhs overwriting the previous rhs.
    ! Do not implement any BC`s! We need the "raw" RHS for the next timestep.
    call cc_generateBasicRHS (rproblem,rrhs)
    
    ! Add w_3*f_{n+1} to the current RHS.     
    call lsysbl_vectorLinearComb(rrhs,rtempVectorRhs,&
         rtimestepping%dweightNewRHS,1.0_DP)

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
    
    rnonlinearIterationTmp%dalpha = 1.0_DP
    rnonlinearIterationTmp%dtheta = rtimestepping%dweightMatrixLHS
    rnonlinearIterationTmp%dgamma = rtimestepping%dweightMatrixLHS * real(1-rproblem%iequation,DP)
    rnonlinearIterationTmp%deta   = 1.0_DP
    rnonlinearIterationTmp%dtau   = 1.0_DP

    ! Scale the pressure by the length of the time step. The core equation routines
    ! handle the equation
    !   alpha*M*u + theta*nu*Laplace*u + gamma*N(u)u + B*p = ...
    ! but we want to solve
    !   alpha*M*u + theta*nu*Laplace*u + gamma*N(u)u + tstep*B*p = ...
    !
    ! So the trick is to scale p by tstep, solve the core equation
    !   alpha*M*u + theta*nu*Laplace*u + gamma*N(u)u + B*(tstep*p) = ...
    ! and scale it back afterwards.
    !
    ! Note that there is an error in the book of [Turek] describing the factor
    ! in front of the pressure in the Crank Nicolson scheme! The pressure is
    ! handled fully implicitely. There is no part of the pressure on the RHS
    ! of the time step scheme and so the factor in front of the pressure
    ! is always the length of the current (sub)step!
    call lsyssc_scaleVector (rvector%RvectorBlock(NDIM2D+1),&
        rtimestepping%dtstep)
    
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
    call lsyssc_scaleVector (rvector%RvectorBlock(NDIM2D+1),&
        1.0_DP/rtimestepping%dtstep)

    ! rvector is the solution vector u^{n+1}.    
  
    ! Finally tell the time stepping scheme that we completed the time step.
    call timstp_nextSubstep (rtimestepping)

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
    
    ! Time error analysis and adaptive time stepping variables
    type(t_timestepSnapshot) :: rsnapshotLastMacrostep
    type(t_vectorBlock) :: rpredictedSolution,roldSolution
    real(dp) :: doldtime
    
    logical :: babortTimestep
    
    integer :: irepetition
    integer(I32) :: isolverStatus,isolverStatusPredictor
    type(t_timeError) :: rtimeerror
    type(t_timeDerivatives) :: rtimeDerivative
    
    ! Timer for the current timestep and the total time.
    type(t_timer) :: rtimerTimestep,rtimerAllTimesteps

    ! Backup postprocessing structure
    type(t_c2d2postprocessing) :: rpostprocessingBackup

    ! Some preparations for the nonlinear solver.
    !
    ! Initialise the nonlinear solver node rnlSol with parameters from
    ! the INI/DAT files.
    call cc_getNonlinearSolver (rnlSol, rproblem%rparamList, 'CC2D-NONLINEAR')

    ! Initialise the nonlinear loop. This is to prepare everything for
    ! our callback routines that are called from the nonlinear solver.
    ! The preconditioner in that structure is initialised later.
    call cc_initNonlinearLoop (rproblem,rproblem%NLMIN,rproblem%NLMAX,&
        rvector,rrhs,rnonlinearIteration,'CC2D-NONLINEAR')
    
    ! Initialise the time stepping scheme according to the problem configuration
    call cc_initTimeSteppingScheme (rproblem%rparamList,rtimestepping)
    
    ! Initialise the preconditioner for the nonlinear iteration
    call cc_initPreconditioner (rproblem,&
        rnonlinearIteration,rvector,rrhs)

    ! Create temporary vectors we need for the nonlinear iteration.
    call lsysbl_createVecBlockIndirect (rrhs, rtempBlock1, .false.)
    call lsysbl_createVecBlockIndirect (rrhs, rtempBlock2, .false.)
    ! First time step
    rproblem%rtimedependence%itimeStep = 1
    rproblem%rtimedependence%dtime = rproblem%rtimedependence%dtimeInit
    dtimederivative = rproblem%rtimedependence%dminTimeDerivative
    
    ! Discretise the boundary conditions at the initial time.
    call cc_updateDiscreteBC (rproblem)

    ! Implement the initial boundary conditions into the solution vector.
    ! Do not implement anything to matrices or RHS vector as these are
    ! maintained in the timeloop.

    call cc_implementBC (rproblem,rvector,rrhs,.true.,.false.)

    ! Postprocessing. Write out the initial solution.
    call output_line ('Starting postprocessing of initial solution...')
    call cc_postprocessingNonstat (rproblem,&
        rvector,rproblem%rtimedependence%dtimeInit,&
        rvector,rproblem%rtimedependence%dtimeInit,0,rpostprocessing,rtimestepping)

    
    ! Reset counter of current macro step repetitions.
    irepetition = 0
    
    !----------------------------------------------------
    ! Timeloop
    !----------------------------------------------------
    
    ! Start timers that calculate the current time.
    call stat_clearTimer(rtimerAllTimesteps)
    call stat_startTimer(rtimerAllTimesteps)
    
    do while ((rproblem%rtimedependence%itimeStep .le. &
               rproblem%rtimedependence%niterations) .and. &
              (rproblem%rtimedependence%dtime .lt. &
               rproblem%rtimedependence%dtimemax-100.0_DP*SYS_EPSREAL) .and. &
              (dtimederivative .ge. &
               rproblem%rtimedependence%dminTimeDerivative))

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
      ! The '3*dtstepFixed' step size just 'coincidentally' coincides with the 
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
          call output_line ('Predictor-Step at Time-Step '// &
              trim(sys_siL(rproblem%rtimedependence%itimeStep,6))// &
              ', Time = '// &
              trim(sys_sdL(rproblem%rtimedependence%dtime,5))// &
              ', Stepsize: DT1 = '// &
              trim(sys_sdL(rtimesteppingPredictor%dtstep,5)),&
              coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
          call output_separator(OU_SEP_EQUAL,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
          !CALL output_line (&
          !  'Macro step '//TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6))// &
          !  '      Time '//TRIM(sys_sdL(rproblem%rtimedependence%dtime,5)) // &
          !  '      Step size '//TRIM(sys_sdL(rtimesteppingPredictor%dtstep,5)))
          !CALL output_lbrk ()
          
          ! Proceed in time, calculate the predicted solution. 
          call lsysbl_copyVector (rvector,rpredictedSolution)              
          call cc_performTimestep (rproblem,rpredictedSolution,rrhs,&
              rtimesteppingPredictor,rnonlinearIteration,rnlSol,&
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
            call output_line ('Timestepping by '&
                //trim(sys_siL(irepetition,2)) &
                //' (' &
                //trim(sys_siL(rproblem%rtimedependence%itimeStep,6)) &
                //'), New Stepsize = ' &
                //trim(sys_sdEP(dtmp,9,2)) &
                //', Old Stepsize = ' &
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
        !  'Time step '//TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6))// &
        !  '     Time '//TRIM(sys_sdL(rproblem%rtimedependence%dtime,5)) // &
        !  '     Step size '//TRIM(sys_sdL(rtimestepping%dtstep,5)))
        call output_lbrk (coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
        call output_separator(OU_SEP_AT,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
        call output_line ('Time-Step '// &
            trim(sys_siL(rproblem%rtimedependence%itimeStep,6))// &
            ', Time = '// &
            trim(sys_sdL(rproblem%rtimedependence%dtime,5))// &
            ', Stepsize: DT3 = '// &
            trim(sys_sdL(rtimestepping%dtstep,5)),&
            coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
        call output_separator(OU_SEP_AT,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
        
        ! Snapshot the current solution for the later calculation of the
        ! time derivative.
        call lsysbl_copyVector (rvector,roldSolution)
        doldtime = rproblem%rtimedependence%dtime
        
        ! Proceed the next time step -- if we are allowed to.
        call cc_performTimestep (rproblem,rvector,rrhs,&
            rtimestepping,rnonlinearIteration,rnlSol,rtempBlock1,rtempBlock2)
            
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
        call output_line ('Time-Step ' &
            //trim(sys_siL(rproblem%rtimedependence%itimeStep,6)) &
            //' (Repetition = '//trim(sys_siL(irepetition,2)) &
            //', Substep = ' &
            //trim(sys_siL(i,6))//' of '//trim(sys_siL(j,6)) &
            //') at time = ' &
            //trim(sys_sdL(rproblem%rtimedependence%dtime,5)) &
            //' finished. ',coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )

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
          call output_line ('New Stepsize = ' &
              //trim(sys_sdEP(dtmp,9,2)) &
              //', Old Stepsize = ' &
              //trim(sys_sdEP(rtimeStepping%dtstepFixed,9,2)),&
              coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
          call output_separator(OU_SEP_AT,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)

          ! Accept the new step size
          call timstp_setBaseSteplength (rtimeStepping, dtmp)
        
        else

        ! Did the solver break down?
        if (rnlSol%iresult .lt. 0) then
          call output_line ('Accuracy notice: Nonlinear solver did not reach '// &
                            'the convergence criterion!',&
                            coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
        else if (rnlSol%iresult .gt. 0) then
          ! Oops, not really good. 
          babortTimestep = .true.

          ! Do we have a time stepping algorithm that allowes recomputation?          
          select case (rproblem%rtimedependence%radaptiveTimeStepping%ctype)
          case (TADTS_FIXED,TADTS_PREDICTION)
            ! That is bad. Our solution is garbage!
            ! We do not do anything in this case. The repetition technique will
            ! later decide on whether to repeat the step or to stop the 
            ! computation.
            call output_line ('Nonlinear solver broke down. Solution probably garbage!',&
                coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
            
          case (TADTS_PREDICTREPEAT,TADTS_PREDREPTIMECONTROL)
            ! Yes, we have. 
            call output_line ('Nonlinear solver broke down. '// &
                              'Calculating new time step size...',&
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
            !CALL output_line ('Timestepping by '&
            !    //TRIM(sys_siL(irepetition,2)) &
            !    //' (' &
            !    //TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6)) &
            !    //'), New Stepsize = ' &
            !    //TRIM(sys_sdEP(dtmp,9,2)) &
            !    //', Old Stepsize = ' &
            !    //TRIM(sys_sdEP(rtimeStepping%dtstepFixed,9,2)) )
            call output_line ('New Stepsize = ' &
                //trim(sys_sdEP(dtmp,9,2)) &
                //', Old Stepsize = ' &
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
            call output_line ('Macrostep completed. Analysing time error and '// &
                              'computing new time step size...',&
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
            call output_line ('Time error: ' &
                //' U(L2)=' &
                //trim(sys_sdEP(rtimeError%drelUL2/real(i,DP),8,2)) &
                //'  U(MX)=' &
                //trim(sys_sdEP(rtimeError%drelUmax/real(i,DP),8,2)) &
                //'  P(L2)=' &
                //trim(sys_sdEP(rtimeError%drelPL2/real(i,DP),8,2)) &
                //'  P(MX)=' &
                //trim(sys_sdEP(rtimeError%drelPmax/real(i,DP),8,2)),&
                coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
            
            ! Tell the user that we have a new time step size.
            call output_line ('New Stepsize = ' &
                //trim(sys_sdEP(dtmp,9,2)) &
                //', Old Stepsize = ' &
                //trim(sys_sdEP(rtimeStepping%dtstepFixed,9,2)),&
                coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
            !CALL output_line ('Timestepping by '&
            !    //TRIM(sys_siL(irepetition,2)) &
            !    //' (' &
            !    //TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6)) &
            !    //'), New Stepsize = ' &
            !    //TRIM(sys_sdEP(dtmp,9,2)) &
            !    //', Old Stepsize = ' &
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
        call output_line ('Starting postprocessing of the time step...',&
            coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
        
        ! Postprocessing. Write out the solution if it was calculated successfully.
        call cc_postprocessingNonstat (rproblem,&
            roldSolution,doldtime,&
            rvector,rproblem%rtimedependence%dtime,&
            rproblem%rtimedependence%itimeStep,rpostprocessing,rtimeStepping)
        
        call output_separator(OU_SEP_MINUS,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
        call output_line ('Analysing time derivative...',&
            coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
        
        ! Calculate the norm of the time derivative. This allowes the DO-loop
        ! above to check if the solution got stationary.
        dtimeDerivative = cc_timeDerivative (&
            rproblem%rtimedependence%radaptiveTimeStepping%cadTimeStepErrorControl,&
            rvector,roldSolution,rtimestepping%dtstep,rtempBlock1,rtimeDerivative)
            
        ! Print the results of the time analysis.
        
        call output_line ('Time derivative:  ' &
            //' RELU(L2)=' &
            //trim(sys_sdEP(rtimeDerivative%drelUL2,9,2)) &
            //'  RELP(L2)=' &
            //trim(sys_sdEP(rtimeDerivative%drelPL2,9,2)) &
            //'  REL=' &
            //trim(sys_sdEP(dtimeDerivative,9,2)),&
            coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
            
        if (dtimederivative .lt. rproblem%rtimedependence%dminTimeDerivative) then
          call output_line ('Solution reached stationary status. Stopping simulation...',&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
        end if
        !CALL output_line ('#'&
        !    //TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6))&
        !    //' (' &
        !    //TRIM(sys_siL(i,6)) &
        !    //') TIME=' &
        !    //TRIM(sys_sdL(rproblem%rtimedependence%dtime,5)) &
        !    //' RELU(L2)=' &
        !    //TRIM(sys_sdEP(rtimeDerivative%drelUL2,9,2)) &
        !    //' RELP(L2)=' &
        !    //TRIM(sys_sdEP(rtimeDerivative%drelPL2,9,2)) &
        !    //' REL=' &
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
        ! time steps in adaptive time stepping technique 1!
        
        ! Do we have a time stepping algorithm that allowes recomputation?          
        select case (rproblem%rtimedependence%radaptiveTimeStepping%ctype)
        case (TADTS_FIXED,TADTS_USERDEF) 
          ! That is bad. Our solution is most probably garbage!
          ! We cancel the timeloop, it does not make any sense to continue.
          call output_line ('Solution garbage! Stopping simulation.',&
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
            
            call output_line ('Repeating macrostep. Returning to timestep ' &
                //trim(sys_siL(rproblem%rtimedependence%itimeStep,6))//'.',&
                coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
            call output_separator(OU_SEP_AT,coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
            
            ! Repeat the time step            
            cycle
              
          else
          
            ! No repetitions left. Reset the repetition counter and
            ! continue with the next timestep, although the
            ! solution might be not really meaningful anymore.
            irepetition = 0
            
            call output_line ('No repetitions left. Cannot repeat macrostep. ' &
                //'Continuing with the next one...',&
                coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
            
          end if
          
        end select
      
      else  ! IF (babortTimestep) THEN
      
        ! No, time step is ok. If this is the last time step of the macrostep,
        ! reset the repetition counter such that it starts with 0 for the
        ! next macrostep.
        if (mod(rproblem%rtimedependence%itimeStep,3) .eq. 0) then
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
    if (roldSolution%NEQ .ne. 0) &
      call lsysbl_releaseVector (roldSolution)
    call lsysbl_releaseVector (rtempBlock2)
    call lsysbl_releaseVector (rtempBlock1)
    
    ! Release existing snapshots
    call cc_releaseSnapshop (rsnapshotLastMacrostep)
    
    ! Release the preconditioner
    call cc_releasePreconditioner (rnonlinearIteration)
    
    ! Release parameters of the nonlinear loop, final clean up
    call cc_doneNonlinearLoop (rnonlinearIteration)
             
  end subroutine
  
end module
