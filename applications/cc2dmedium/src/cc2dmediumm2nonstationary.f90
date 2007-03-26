!##############################################################################
!# ****************************************************************************
!# <name> cc2dmediumm2nonstationary </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises a time dependent solver for the coupled Navier-Stokes
!# system.
!#
!# The following routines can be found here:
!#
!# 1.) c2d2_initParTimeDependence
!#     -> Initialise the parameters of the time dependent solver from DAT file
!#        parameters.
!#
!# 2.) c2d2_initTimeSteppingScheme
!#     -> Initialise the time stepping scheme from from DAT file
!#        parameters.
!#
!# 3.) c2d2_performTimestep
!#     -> Calculates one time step.
!#
!# 4.) c2d2_solveNonstationary
!#     -> Realises a time-loop to solve multiple time steps. Solves the
!#        nonstationary Navier-Stokes system.
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2nonstationary

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE triangulation
  USE spatialdiscretisation
  USE coarsegridcorrection
  USE spdiscprojection
  USE nonlinearsolver
  USE paramlist
  USE linearsolverautoinitialise
  USE matrixrestriction
  USE paramlist
  USE timestepping
  
  USE collection
  USE convection
    
  USE cc2dmediumm2basic
  USE cc2dmedium_callback

  USE cc2dmediumm2nonlinearcore
  USE cc2dmediumm2nonlinearcoreinit
  USE cc2dmediumm2stationary
  USE adaptivetimestep
  USE cc2dmediumm2timeanalysis
  USE cc2dmediumm2boundary
  USE cc2dmediumm2discretisation
  USE cc2dmediumm2postprocessing
    
  IMPLICIT NONE

!<types>

!<typeblock>

  ! Saves data about a time step. Is used to make a snapshot of the current
  ! flow/mesh situation. THat way, the solver can make backups that can be
  ! restored if necessary.

  TYPE t_timestepSnapshot
  
    ! Current point in time
    REAL(DP)            :: dtime
    
    ! Current timestep
    INTEGER             :: itimeStep
    
    ! Current timestep configuration
    TYPE(t_explicitTimeStepping) :: rtimeStepping
  
    ! Right hand side at the current point in time
    TYPE(t_vectorBlock) :: rrhs
    
    ! Solution vector at the current point in time
    TYPE(t_vectorBlock) :: rsolution
  
  END TYPE

!</typeblock>

!</types>

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_initParTimeDependence (rproblem,ssection,rparams)
  
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
  TYPE(t_parlist), INTENT(IN) :: rparams
  
  ! The name of the section in the parameter list containing the parameters
  ! for the time dependent simulation.
  CHARACTER(LEN=*), INTENT(IN) :: ssection
!</input>

!<inputoutput>
  ! The problem structure to be initialised.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>
  
!</subroutine>

    ! Fetch the parameters. Initialise with standard settings if they don't 
    ! exist.
    CALL parlst_getvalue_int (rparams,ssection,'itimedependence',  &
        rproblem%itimedependence, 0)
    CALL parlst_getvalue_int (rparams,ssection,'niterations',     &
        rproblem%rtimedependence%niterations, 1000)
    CALL parlst_getvalue_double (rparams,ssection,'dtimestart',   &
        rproblem%rtimedependence%dtimeInit, 0.0_DP)
    CALL parlst_getvalue_double (rparams,ssection,'dtimemax',     &
        rproblem%rtimedependence%dtimeMax, 20.0_DP)
    CALL parlst_getvalue_double (rparams,ssection,'dtimeStep',    &
        rproblem%rtimedependence%dtimeStep, 0.01_DP)
    CALL parlst_getvalue_double (rparams,ssection,'dminTimeDerivative',    &
        rproblem%rtimedependence%dminTimeDerivative, 0.00001_DP)
    rproblem%rtimedependence%itimeStep = 0

    ! Call the initialisation routine for the adaptive time stepping to 
    ! initialise the rest of the parameters.
    CALL adtstp_init (rparams,ssection,&
        rproblem%rtimedependence%radaptiveTimeStepping)

  END SUBROUTINE  
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_backupTimestep (rsnapshot,rproblem,&
      rtimeStepping,rsolution,rrhs)

!<description>
  ! Makes a backup of the current flow configuration into rsnapshot.
  ! If there is any previous data in rsnapshot, it's overwritten/updated.
!</description>

!<inputoutput>
  ! A snapshop structure that receives the snapshot.
  TYPE(t_timestepSnapshot), INTENT(INOUT) :: rsnapshot
  
  ! The problem structure with all information about the current problem.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
  
  ! A time stepping structure that configures the current time step.
  TYPE(t_explicitTimeStepping), INTENT(INOUT), OPTIONAL :: rtimeStepping
  
  ! Right hand side at the current point in time
  TYPE(t_vectorBlock), INTENT(INOUT), OPTIONAL :: rrhs
  
  ! Solution vector at the current point in time
  TYPE(t_vectorBlock), INTENT(INOUT), OPTIONAL :: rsolution
!</inputoutput>

!</subroutine>

    ! Back up all data into rsnapshot.
    rsnapshot%dtime = rproblem%rtimedependence%dtime
    rsnapshot%itimestep = rproblem%rtimedependence%itimestep
    
    IF (PRESENT(rtimeStepping)) &
      rsnapshot%rtimeStepping = rtimeStepping
    
    IF (PRESENT(rrhs)) &
      CALL lsysbl_copyVector (rrhs,rsnapshot%rrhs)

    IF (PRESENT(rsolution)) &
      CALL lsysbl_copyVector (rsolution,rsnapshot%rsolution)
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_restoreTimestep (rsnapshot,rproblem,rtimeStepping,&
      rsolution,rrhs)

!<description>
  ! Restores a backup of a time step from rsnapshot.
!</description>

!<inputoutput>
  ! A snapshop structure that receives the snapshot.
  TYPE(t_timestepSnapshot), INTENT(INOUT) :: rsnapshot
  
  ! The problem structure with all information about the current problem.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
  
  ! A time stepping structure that configures the current time step.
  TYPE(t_explicitTimeStepping), INTENT(INOUT), OPTIONAL :: rtimeStepping
  
  ! Right hand side at the current point in time
  TYPE(t_vectorBlock), INTENT(INOUT), OPTIONAL :: rrhs
  
  ! Solution vector at the current point in time
  TYPE(t_vectorBlock), INTENT(INOUT), OPTIONAL :: rsolution
  
!</inputoutput>

!</subroutine>

    ! Restore all data from rsnapshot.
    rproblem%rtimedependence%dtime = rsnapshot%dtime 
    rproblem%rtimedependence%itimestep = rsnapshot%itimestep
    
    IF (PRESENT(rtimeStepping)) &
      rtimeStepping = rsnapshot%rtimeStepping 
    
    IF (PRESENT(rrhs)) &
      CALL lsysbl_copyVector (rsnapshot%rrhs,rrhs)
      
    IF (PRESENT(rsolution)) &
      CALL lsysbl_copyVector (rsnapshot%rsolution,rsolution)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_releaseSnapshop (rsnapshot)

!<description>
  ! Releases the rsnapshot structure. All allocated memory is released.
!</description>

!<inputoutput>
  ! A snapshop structure that receives the snapshot.
  TYPE(t_timestepSnapshot), INTENT(INOUT) :: rsnapshot
!</inputoutput>

!</subroutine>

    ! Release allocated memory -- if memory is allocated.
    IF (rsnapshot%rrhs%NEQ .NE. 0) &
      CALL lsysbl_releaseVector (rsnapshot%rrhs)
    IF (rsnapshot%rsolution%NEQ .NE. 0) &
      CALL lsysbl_releaseVector (rsnapshot%rsolution)
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_initTimeSteppingScheme (rparams,rtimeStepping)
  
!<description>
  ! Initialises the time stepping scheme according to the parameters in the DAT file.
!</description>

!<input>
  ! A parameter list from a DAT file containing parameters that configure the
  ! time dependence.
  TYPE(t_parlist), INTENT(IN) :: rparams
!</input>

!<inputoutput>
  ! The time stepping scheme structure to be initialised.
  TYPE(t_explicitTimeStepping), INTENT(INOUT) :: rtimeStepping
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: cscheme
    REAL(DP) :: dtheta,dtstep,dtimemin

    ! Get the parameters for the time stepping scheme from the parameter list
    CALL parlst_getvalue_int (rparams, &
        'TIME-DISCRETISATION', 'ITIMESTEPSCHEME', cscheme, 0)
    CALL parlst_getvalue_double (rparams, &
        'TIME-DISCRETISATION', 'DTIMESTEPTHETA', dtheta, 1.0_DP)
    CALL parlst_getvalue_double (rparams, &
        'TIME-DISCRETISATION', 'DTIMESTEP', dtstep, 0.1_DP)
    CALL parlst_getvalue_double (rparams, &
        'TIME-DISCRETISATION', 'DTIMEINIT', dtimemin, 0.0_DP)
    
    ! Initialise the time stepping in the problem structure
    CALL timstp_init (rtimestepping, &
                      cscheme, dtimemin, dtstep, dtheta)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_performTimestep (rproblem,rvector,rrhs,&
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
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! The current solution vector at time $t^n$. Is replaced by the
  ! solution vector at time $t^{n+1}.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector

  ! The RHS vector at time $t^n$. Is replaced by the RHS at time $t^{n+1}$.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhs
  
  ! Configuration block of the time stepping scheme.
  TYPE(t_explicitTimeStepping)        :: rtimestepping

  ! Structure for the nonlinear iteration for solving the core equation.
  TYPE(t_ccnonlinearIteration), INTENT(INOUT) :: rnonlinearIteration
  
  ! A configuration stucture for the nonlinear solver
  TYPE(t_nlsolNode) :: rnlSolver
  
  ! Temporary vectors for the nonlinear solver
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector,rtempVectorRhs
!</inputoutput>

!</subroutine>

    ! local variables
    TYPE(t_ccnonlinearIteration) :: rnonlinearIterationTmp
    
    ! DEBUG!!!
    !REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2
    
    ! Restore the standard matrix structure in case the matrices had been
    ! modified by c2d2_finaliseMatrices for the preconditioner -- i.e. 
    ! temporarily switch to the matrix structure that is compatible to 
    ! the discretisation.
    CALL c2d2_unfinaliseMatrices (rnonlinearIteration, &
        rnonlinearIteration%rfinalAssembly,.FALSE.)
  
    ! The new RHS will be set up in rtempVectorRhs. Assign the discretisation/
    ! boundary conditions of rrhs to that vector so that rtempVectorRhs
    ! acts as a RHS vector.
    CALL lsysbl_assignDiscretIndirect(rrhs,rtempVectorRhs)
    
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
    
    ! Set up w_2*N(u_n) + w_4*f_n.
    CALL lsysbl_vectorLinearComb(rrhs,rtempVectorRhs,&
         rtimestepping%dweightOldRHS,0.0_DP)
    
    ! For setting up M(u_n) + w_2*N(u_n), switch the sign of w_2 and call the method
    ! to calculate the Convection/Diffusion part of the nonlinear defect. This builds 
    ! rtempVectorRhs = rtempVectorRhs - (-Mass)*u - (-w_2) (nu*Laplace*u + grad(u)u).
    !
    ! Don't implement any boundary conditions when assembling this -- it's not
    ! a defect vector!
    ! The BC's are implemented at the end when the full RHS is finished...
    
    rnonlinearIterationTmp = rnonlinearIteration
    CALL c2d2_setupCoreEquation (rnonlinearIterationTmp, &
      -1.0_DP, &
      -rtimestepping%dweightMatrixRHS, &
      -rtimestepping%dweightMatrixRHS * &
       REAL(1-collct_getvalue_int (rproblem%rcollection,'ISTOKES'),DP))

    CALL c2d2_assembleConvDiffDefect (rnonlinearIterationTmp,rvector,rtempVectorRhs,&
        rproblem%rcollection)    

    ! -------------------------------------------    
    ! Switch to the next point in time.
    rproblem%rtimedependence%dtime = rtimestepping%dcurrenttime + rtimestepping%dtstep
          
    ! Discretise the boundary conditions at the new point in time -- 
    ! if the boundary conditions are nonconstant in time!
    IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
      CALL c2d2_updateDiscreteBC (rproblem, .FALSE.)
    END IF
    ! -------------------------------------------    

    ! generate f_n+1 into the rrhs overwriting the previous rhs.
    ! Don't implement any BC's! We need the "raw" RHS for the next timestep.
    CALL c2d2_generateBasicRHS (rproblem,rrhs)
    
    ! Add w_3*f_{n+1} to the current RHS.     
    CALL lsysbl_vectorLinearComb(rrhs,rtempVectorRhs,&
         rtimestepping%dweightNewRHS,1.0_DP)

    ! Implement boundary conditions into the RHS and solution vector, not
    ! into the matrices; the latter is done during the nonlinear iteration.
    CALL c2d2_implementBC (rproblem,rvector,rtempVectorRhs,.FALSE.,.TRUE.,.TRUE.)

    ! That's it for the RHS and solution vector.
    !    
    ! The LHS is "u_{n+1} + w_1*N(u_n+1)" which results in the system matrix
    ! "M + w_1 N(.)" for the next linear system to solve. 
    ! Set up the corresponding core equation in a temporary core-equation
    ! structure.

    rnonlinearIterationTmp = rnonlinearIteration
    CALL c2d2_setupCoreEquation (rnonlinearIterationTmp, &
      1.0_DP, &
      rtimestepping%dweightMatrixLHS, &
      rtimestepping%dweightMatrixLHS * &
       REAL(1-collct_getvalue_int (rproblem%rcollection,'ISTOKES'),DP))
    
    ! Using rfinalAssembly, make the matrices compatible 
    ! to our preconditioner if they are not. So we switch again to the matrix
    ! representation that is compatible to our preconditioner.
    CALL c2d2_finaliseMatrices (rnonlinearIteration)
    
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
    ! Note that there is an error in the book of [Turel] describing the factor
    ! in front of the pressure in the Crank Nicolson scheme! The pressure is
    ! handled fully implicitely. There is no part of the pressure on the RHS
    ! of the time step scheme and so the factor in front of the pressure
    ! is always the length of the current (sub)step!
    CALL lsyssc_scaleVector (rvector%RvectorBlock(NDIM2D+1),&
        rtimestepping%dtstep)
    
    ! Call the solver of the core equation to solve it using a nonlinear
    ! iteration.
    CALL c2d2_solveCoreEquation (rnlSolver,rnonlinearIterationTmp,&
        rvector,rtempVectorRhs,rproblem%rcollection,rtempVector)             

    ! scale the pressure back, then we have again the correct solution vector.
    CALL lsyssc_scaleVector (rvector%RvectorBlock(NDIM2D+1),&
        1.0_DP/rtimestepping%dtstep)

    ! rvector is the solution vector u^{n+1}.    
  
    ! Finally tell the time stepping scheme that we completed the time step.
    CALL timstp_nextSubstep (rtimestepping)

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_solveNonstationary (rproblem,rvector,rrhs)
  
!<description>
  ! Starts the time discretisation. Proceeds in time until the final time
  ! or the maximum number of time steps is reached.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

  ! The initial solution vector. Is replaced by the final solution vector.
  ! Must be unsorted.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
  
  ! The initial RHS vector. Is replaced by the final RHS vector.
  ! Boundary conditions must not be implemented into that vector, this is done
  ! internally depending on the time.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhs
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: i,j
    TYPE(t_vectorBlock) :: rtempBlock1,rtempBlock2
    TYPE(t_ccNonlinearIteration) :: rnonlinearIteration
    
    ! Configuration block of the time stepping scheme 
    TYPE(t_explicitTimeStepping)        :: rtimestepping,rtimesteppingPredictor

    ! The nonlinear solver configuration
    TYPE(t_nlsolNode) :: rnlSol

    REAL(DP) :: dtimederivative,dtmp,dtmperror,dtimeratio
    
    ! Time error analysis and adaptive time stepping variables
    TYPE(t_timestepSnapshot) :: rsnapshotLastMacrostep
    TYPE(t_vectorBlock) :: rpredictedSolution,roldSolution
    
    LOGICAL :: babortTimestep
    
    INTEGER :: irepetition
    INTEGER(I32) :: isolverStatus,isolverStatusPredictor
    TYPE(t_timeError) :: rtimeerror
    TYPE(t_timeDerivatives) :: rtimeDerivative

    ! Some preparations for the nonlinear solver.
    !
    ! Initialise the nonlinear solver node rnlSol with parameters from
    ! the INI/DAT files.
    CALL c2d2_getNonlinearSolver (rnlSol, rproblem%rparamList, 'CC2D-NONLINEAR')

    ! Initialise the nonlinear loop. This is to prepare everything for
    ! or callback routines that are called from the nonlinear solver.
    ! The preconditioner in that structure is initialised later.
    CALL c2d2_initNonlinearLoop (rproblem,rvector,rrhs,rnonlinearIteration,&
        'CC2D-NONLINEAR')
    
    ! Initialise the time stepping scheme according to the problem configuration
    CALL c2d2_initTimeSteppingScheme (rproblem%rparamList,rtimestepping)
    
    ! Check the matrices if they are compatible to our
    ! preconditioner. If not, we later have to modify the matrices a little
    ! bit to make it compatible. 
    ! The result of this matrix analysis is saved to the rfinalAssembly structure 
    ! in rnonlinearIteration and allows us later to switch between these two
    ! matrix representations: Compatibility to the discretisation routines
    ! and compatibity to the preconditioner.
    ! The c2d2_checkAssembly routine below uses this information to perform
    ! the actual modification in the matrices.
    CALL c2d2_checkAssembly (rproblem,rrhs,rnonlinearIteration%rfinalAssembly)
    
    ! Using rfinalAssembly as computed above, make the matrices compatible 
    ! to our preconditioner if they are not.
    CALL c2d2_finaliseMatrices (rnonlinearIteration)
    
    ! Initialise the preconditioner for the nonlinear iteration
    CALL c2d2_preparePreconditioner (rproblem,&
        rnonlinearIteration,rvector,rrhs)

    ! Create temporary vectors we need for the nonlinear iteration.
    CALL lsysbl_createVecBlockIndirect (rrhs, rtempBlock1, .FALSE.)
    CALL lsysbl_createVecBlockIndirect (rrhs, rtempBlock2, .FALSE.)

    ! Implement the initial boundary conditions into the solution vector.
    ! Don't implement anything to matrices or RHS vector as these are
    ! maintained in the timeloop.
    ! Afterwards, we can start the timeloop.
    CALL c2d2_implementBC (rproblem,rvector,rrhs,.FALSE.,.TRUE.,.FALSE.)

    ! First time step
    rproblem%rtimedependence%itimeStep = 1
    rproblem%rtimedependence%dtime = rproblem%rtimedependence%dtimeInit
    dtimederivative = rproblem%rtimedependence%dminTimeDerivative
    
    ! Reset counter of current macro step repetitions.
    irepetition = 0
    
    !----------------------------------------------------
    ! Timeloop
    !----------------------------------------------------
    
    DO WHILE ((rproblem%rtimedependence%itimeStep .LE. &
               rproblem%rtimedependence%niterations) .AND. &
              (rproblem%rtimedependence%dtime .LT. &
               rproblem%rtimedependence%dtimemax) .AND. &
              (dtimederivative .GE. &
               rproblem%rtimedependence%dminTimeDerivative))
              
      ! The babortTimestep is normally FALSE. If set to TRUE, the computation
      ! of the next time step is aborted because of an error. 
      
      babortTimestep = .FALSE.
      
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
      ! large mactostep
      !
      ! Do we use adaptive time stepping?
      SELECT CASE (rproblem%rtimedependence%radaptiveTimeStepping%ctype)
      CASE (TADTS_FIXED) 
        ! Nothing to be done
        
      CASE (TADTS_PREDICTION,TADTS_PREDICTREPEAT,TADTS_PREDREPTIMECONTROL)
        ! Adaptive time stepping. Is this the first substep?
        IF (MOD(rproblem%rtimedependence%itimeStep,3) .EQ. 1) THEN
        
          ! If this is not a repetition of a (macro-) timestep, 
          ! create a backup of the current flow situation, so we can repeat the time
          ! step if necessary.
          IF (irepetition .eq. 0) THEN
            CALL c2d2_backupTimestep (rsnapshotLastMacrostep,rproblem,rtimeStepping,&
                rvector,rrhs)
          END IF

          ! At first, perform one predictor step with 3x time step size
          ! using implicit Euler.
          ! For this purpose, create a new time stepping structure
          ! based on the standard one (for small time steps).
          CALL timstp_init (rtimesteppingPredictor, TSCHM_ONESTEP, &
                rproblem%rtimedependence%dtime, &
                3.0_DP*rtimestepping%dtstepFixed, 1.0_DP)

          CALL output_lbrk ()
          CALL output_separator(OU_SEP_EQUAL)
          CALL output_line ('Predictor-Step at Time-Step '// &
              TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6))// &
              ', Time = '// &
              TRIM(sys_sdL(rproblem%rtimedependence%dtime,5))// &
              ', Stepsize: DT1 = '// &
              TRIM(sys_sdL(rtimesteppingPredictor%dtstep,5)) )
          CALL output_separator(OU_SEP_EQUAL)
          !CALL output_line (&
          !  'Macro step '//TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6))// &
          !  '      Time '//TRIM(sys_sdL(rproblem%rtimedependence%dtime,5)) // &
          !  '      Step size '//TRIM(sys_sdL(rtimesteppingPredictor%dtstep,5)))
          !CALL output_lbrk ()
          
          ! Proceed in time, calculate the predicted solution. 
          CALL lsysbl_copyVector (rvector,rpredictedSolution)              
          CALL c2d2_performTimestep (rproblem,rpredictedSolution,rrhs,&
              rtimesteppingPredictor,rnonlinearIteration,rnlSol,&
              rtempBlock1,rtempBlock2)
              
          ! Remember the status of the nonlinear solver for later use.
          isolverStatusPredictor = rnlSol%iresult
          
          ! Did the nonlinear solver break down?
          IF (rnlSol%iresult .GT. 0) THEN
            ! Oops, not really good. 
            babortTimestep = .TRUE.
          
            ! Calculate a new time step size.
            ! Set bit 0 and not bit 2/3 in isolverStatus as we want to compute the 
            ! new time step based on the solver status of the last computation and
            ! not on a combined analysis of predictor step and actual solution!
            isolverStatus = 0
            SELECT CASE (rnlSol%iresult)
            CASE (:-1)
              ! Nonlinear solver worked, but could not reach convergence criterion
              isolverStatus = IOR(isolverStatus,TADTS_SST_NLINCOMPLETE)
            CASE (0)
              ! Everything fine
            CASE (1)
              ! Nonlinear solver diverged
              isolverStatus =  IOR(isolverStatus,TADTS_SST_NLFAIL)
            CASE (2)
              ! Nonlinear solver diverged because of error in the preconditioner
              isolverStatus =  IOR(isolverStatus,&
                                   TADTS_SST_NLFAIL + TADTS_SST_NLPRECFAIL)
            CASE (3)
              ! General error
              isolverStatus = NOT(0)
            END SELECT
            dtmp = adtstp_calcTimeStep (&
                rproblem%rtimedependence%radaptiveTimeStepping, &
                0.0_DP, &
                rproblem%rtimedependence%dtimeInit,&
                rproblem%rtimedependence%dtime, &
                rtimeStepping%dtstepFixed, &
                timstp_getOrder(rtimeStepping), &
                isolverStatus,irepetition) 
                
            ! Tell the user that we have a new time step size.
            CALL output_separator(OU_SEP_AT)
            CALL output_line ('Timestepping by '&
                //TRIM(sys_siL(irepetition,2)) &
                //' (' &
                //TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6)) &
                //'), New Stepsize = ' &
                //TRIM(sys_sdEP(dtmp,9,2)) &
                //', Old Stepsize = ' &
                //TRIM(sys_sdEP(rtimeStepping%dtstepFixed,9,2)) )
            CALL output_separator(OU_SEP_AT)

            ! Accept the new step size
            CALL timstp_setBaseSteplength (rtimeStepping, dtmp)
            
            ! The old RHS is restored later...
            
          ELSE

            ! Restore the flow situation to the beginning of the macrostep.
            ! We didn't change the solution vector and the time stepping structure,
            ! so we don't have to restore them.
            CALL c2d2_restoreTimestep (rsnapshotLastMacrostep,rproblem,rrhs=rrhs)

            ! The solver worked, rpredictedSolution contains the predicted
            ! solution. Now we can continue with the usual time stepping.
            
          END IF
            
        END IF
      END SELECT
      
      !----------------------------------------------------
      ! Proceed one step in time
      !----------------------------------------------------
      IF (.NOT. babortTimestep) THEN
      
        !CALL output_separator(OU_SEP_MINUS)
        !CALL output_line (&
        !  'Time step '//TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6))// &
        !  '     Time '//TRIM(sys_sdL(rproblem%rtimedependence%dtime,5)) // &
        !  '     Step size '//TRIM(sys_sdL(rtimestepping%dtstep,5)))
        CALL output_lbrk ()
        CALL output_separator(OU_SEP_AT)
        CALL output_line ('Time-Step '// &
            TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6))// &
            ', Time = '// &
            TRIM(sys_sdL(rproblem%rtimedependence%dtime,5))// &
            ', Stepsize: DT3 = '// &
            TRIM(sys_sdL(rtimestepping%dtstep,5)) )
        CALL output_separator(OU_SEP_AT)
        
        ! Snapshot the current solution for the later calculation of the
        ! time derivative.
        CALL lsysbl_copyVector (rvector,roldSolution)
        
        ! Proceed the next time step -- if we are allowed to.
        CALL c2d2_performTimestep (rproblem,rvector,rrhs,&
            rtimestepping,rnonlinearIteration,rnlSol,rtempBlock1,rtempBlock2)
            
        ! Do we count in steps a 1 or in steps a 3?
        ! Respecting this, i is assigned the number of the substep in the
        ! macrostep.
        SELECT CASE (rproblem%rtimedependence%radaptiveTimeStepping%ctype)
        CASE (TADTS_FIXED) 
          i = 1
          j = 1
        CASE (TADTS_PREDICTION,TADTS_PREDICTREPEAT,TADTS_PREDREPTIMECONTROL)
          i = MOD(rproblem%rtimedependence%itimeStep-1,3)+1
          j = 3
        END SELECT
        
        CALL output_separator(OU_SEP_AT)
        CALL output_line ('Time-Step ' &
            //TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6)) &
            //' (Repetition = '//TRIM(sys_siL(irepetition,2)) &
            //', Substep = ' &
            //TRIM(sys_siL(i,6))//' of '//TRIM(sys_siL(j,6)) &
            //') at time = ' &
            //TRIM(sys_sdL(rproblem%rtimedependence%dtime,5)) &
            //' finished. ' )

        ! Did the solver break down?
        IF (rnlSol%iresult .LT. 0) THEN
          CALL output_line ('Accuracy notice: Nonlinear solver did not reach '// &
                            'the convergence criterion!')
        ELSE IF (rnlSol%iresult .GT. 0) THEN
          ! Oops, not really good. 
          babortTimestep = .TRUE.

          ! Do we have a time stepping algorithm that allowes recomputation?          
          SELECT CASE (rproblem%rtimedependence%radaptiveTimeStepping%ctype)
          CASE (TADTS_FIXED,TADTS_PREDICTION)
            ! That's bad. Our solution is garbage!
            ! We don't do anything in this case. The repetition technique will
            ! later decide on whether to repeat the step or to stop the 
            ! computation.
            CALL output_line ('Nonlinear solver broke down. Solution probably garbage!')
            
          CASE (TADTS_PREDICTREPEAT,TADTS_PREDREPTIMECONTROL)
            ! Yes, we have. 
            CALL output_line ('Nonlinear solver broke down. '// &
                              'Calculating new time step size...')
            
            ! Calculate a new time step size.
            SELECT CASE (rnlSol%iresult)
            CASE (:-1)
              ! Nonlinear solver worked, but could not reach convergence criterion
              isolverStatus = IOR(isolverStatus,TADTS_SST_NLINCOMPLETE)
            CASE (0)
              ! Everything fine
            CASE (1)
              ! Nonlinear solver diverged
              isolverStatus =  IOR(isolverStatus,TADTS_SST_NLFAIL)
            CASE (2)
              ! Nonlinear solver diverged because of error in the preconditioner
              isolverStatus =  IOR(isolverStatus,&
                                   TADTS_SST_NLFAIL + TADTS_SST_NLPRECFAIL)
            CASE (3)
              ! General error
              isolverStatus = NOT(0)
            END SELECT
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
            CALL output_line ('New Stepsize = ' &
                //TRIM(sys_sdEP(dtmp,9,2)) &
                //', Old Stepsize = ' &
                //TRIM(sys_sdEP(rtimeStepping%dtstepFixed,9,2)) )
            CALL output_separator(OU_SEP_AT)

            ! Accept the new step size
            CALL timstp_setBaseSteplength (rtimeStepping, dtmp)

          END SELECT
          
        END IF  
            
      END IF
      
      !----------------------------------------------------
      ! Time error control
      !----------------------------------------------------
      IF (.NOT. babortTimestep) THEN
      
        ! Ok, everything worked fine, we have a valid solution!
        !        
        ! Time step control. Do we have a time stepping algorithm that
        ! adapts the time step size and probably wants to repeat the 
        ! calculation?
        SELECT CASE (rproblem%rtimedependence%radaptiveTimeStepping%ctype)
        CASE (TADTS_FIXED) 
          ! No, continue as usual.
         
        CASE (TADTS_PREDICTION,TADTS_PREDICTREPEAT,TADTS_PREDREPTIMECONTROL)
          
          ! Yes. Time step adaption takes place after 3 steps, where we reached
          ! the same point in time like the big macrostep.
          
          IF (MOD(rproblem%rtimedependence%itimeStep,3) .EQ. 0) THEN
          
            CALL output_separator (OU_SEP_MINUS)
            CALL output_line ('Macrostep completed. Analysing time error and '// &
                              'computing new time step size...')
          
            ! Calculate the new time step size.
            ! This is based on the solution, the predicted solution, the order
            ! of the time stepping algorithm and on the status of the solvers.
            !
            ! At first, calculate the time error.
            dtmperror =  c2d2_timeErrorByPredictor (&
                rproblem%rtimedependence%radaptiveTimeStepping%cadTimeStepErrorControl,&
                rvector,rpredictedSolution,rtempBlock1,rtimeError)            

            ! Evaluate everything that went wrong in the solvers
            isolverStatus = 0
            SELECT CASE (rnlSol%iresult)
            CASE (:-1)
              ! Nonlinear solver worked, but could not reach convergence criterion
              isolverStatus = IOR(isolverStatus,TADTS_SST_NLINCOMPLETE)
            CASE (0)
              ! Everything fine
            CASE (1)
              ! Nonlinear solver diverged
              isolverStatus =  IOR(isolverStatus,TADTS_SST_NLFAIL)
            CASE (2)
              ! Nonlinear solver diverged because of error in the preconditioner
              isolverStatus =  IOR(isolverStatus,&
                                   TADTS_SST_NLFAIL + TADTS_SST_NLPRECFAIL)
            CASE (3)
              ! General error
              isolverStatus = NOT(0)
            END SELECT

            SELECT CASE (isolverStatusPredictor)
            CASE (:-1)
              ! Nonlinear solver worked, but could not reach convergence criterion
              isolverStatus = IOR(isolverStatus,TADTS_SST_NLPREDINCOMPLETE)
            CASE (0)
              ! Everything fine
            CASE (1)
              ! Nonlinear solver diverged
              isolverStatus =  IOR(isolverStatus,TADTS_SST_NLPREDFAIL)
            CASE (2)
              ! Nonlinear solver diverged because of error in the preconditioner
              isolverStatus =  IOR(isolverStatus,&
                                   TADTS_SST_NLPREDFAIL + TADTS_SST_NLPREDPRECFAIL)
            CASE (3)
              ! General error
              isolverStatus = NOT(0)
            END SELECT

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
            IF (rproblem%rtimedependence%radaptiveTimeStepping%ctype .EQ. &
                TADTS_PREDREPTIMECONTROL) THEN

              ! When the new time step is much smaller than the old one, 
              ! set babortTimestep to TRUE to indicate that
              ! the time-step has to be repeated.
                
              IF (dtimeratio .LT. rproblem%rtimedependence%radaptiveTimeStepping% &
                                  depsAdaptiveRelTimeStep) THEN
                babortTimestep = .TRUE.
              END IF
              
            END IF

            ! Assign i either 2 or 8, depending on whether we have a time stepping
            ! algorithm of order one or two in time!
            i = 2 + 6*(timstp_getOrder(rtimeStepping)-1)

            ! Print the result of the time error analysis
            CALL output_line ('Time error: ' &
                //' U(L2)=' &
                //TRIM(sys_sdEP(rtimeError%drelUL2/REAL(i,DP),8,2)) &
                //'  U(MX)=' &
                //TRIM(sys_sdEP(rtimeError%drelUmax/REAL(i,DP),8,2)) &
                //'  P(L2)=' &
                //TRIM(sys_sdEP(rtimeError%drelPL2/REAL(i,DP),8,2)) &
                //'  P(MX)=' &
                //TRIM(sys_sdEP(rtimeError%drelPmax/REAL(i,DP),8,2)) )
            
            ! Tell the user that we have a new time step size.
            CALL output_line ('New Stepsize = ' &
                //TRIM(sys_sdEP(dtmp,9,2)) &
                //', Old Stepsize = ' &
                //TRIM(sys_sdEP(rtimeStepping%dtstepFixed,9,2)) )
            !CALL output_line ('Timestepping by '&
            !    //TRIM(sys_siL(irepetition,2)) &
            !    //' (' &
            !    //TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6)) &
            !    //'), New Stepsize = ' &
            !    //TRIM(sys_sdEP(dtmp,9,2)) &
            !    //', Old Stepsize = ' &
            !    //TRIM(sys_sdEP(rtimeStepping%dtstepFixed,9,2)) )

            ! Accept the new step size
            CALL timstp_setBaseSteplength (rtimeStepping, dtmp)

          END IF
        
        END SELECT
      END IF
        
      !----------------------------------------------------
      ! Postprocessing
      !----------------------------------------------------
      IF (.NOT. babortTimestep) THEN
      
        ! Ok, everything worked fine, we have a valid solution of 
        ! our current substep.
        
        CALL output_separator(OU_SEP_MINUS)
        CALL output_line ('Starting postprocessing of the time step...')
        
        ! Postprocessing. Write out the solution if it was calculated successfully.
        CALL c2d2_postprocessingNonstat (rproblem,rvector)
        
        CALL output_separator(OU_SEP_MINUS)
        CALL output_line ('Analysing time derivative...')
        
        ! Calculate the norm of the time derivative. This allowes the DO-loop
        ! above to check if the solution got stationary.
        dtimeDerivative = c2d2_timeDerivative (&
            rproblem%rtimedependence%radaptiveTimeStepping%cadTimeStepErrorControl,&
            rvector,roldSolution,rtimestepping%dtstep,rtempBlock1,rtimeDerivative)
            
        ! Print the results of the time analysis.
        
        CALL output_line ('Time derivative:  ' &
            //' RELU(L2)=' &
            //TRIM(sys_sdEP(rtimeDerivative%drelUL2,9,2)) &
            //'  RELP(L2)=' &
            //TRIM(sys_sdEP(rtimeDerivative%drelPL2,9,2)) &
            //'  REL=' &
            //TRIM(sys_sdEP(dtimeDerivative,9,2)) )
            
        IF (dtimederivative .LT. rproblem%rtimedependence%dminTimeDerivative) THEN
          CALL output_line ('Solution reached stationary status. Stopping simulation...')
        END IF
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
        
      END IF
      
      !----------------------------------------------------
      ! Check if the time step must be repeated
      !----------------------------------------------------
      IF (babortTimestep) THEN
      
        ! Uh, oh, somthing went wrong. Let's hope we can repeat the timestep!
        !
        ! We have to repeat the time step if
        !  a) we have repetitions left and
        !  b) The user has activated any kind of adaptive time stepping and
        !  c1) the parameters force us to repeat a step or
        !  c2) the solver for the nonlinear iteration (predictor
        !      or corrector step) broke down
        ! c2) is an exception to the rule that we don't repeat
        ! time steps in adaptive time stepping technique 1!
        
        ! Do we have a time stepping algorithm that allowes recomputation?          
        SELECT CASE (rproblem%rtimedependence%radaptiveTimeStepping%ctype)
        CASE (TADTS_FIXED) 
          ! That's bad. Our solution is most probably garbage!
          ! We cancel the timeloop, it doesn't make any sense to continue.
          CALL output_line ('Solution garbage! Stopping simulation.')
          CALL output_separator(OU_SEP_AT)
          EXIT
          
        CASE (TADTS_PREDICTION,TADTS_PREDICTREPEAT,TADTS_PREDREPTIMECONTROL)
          ! So far, so good. Are we allowed to recompute the time step?
          IF (irepetition .LT. &
              rproblem%rtimedependence%radaptiveTimeStepping%nrepetitions) THEN
            
            ! Lucky, we are allowed to recompute :-)
            ! The previous situation was restored before -- so increase the
            ! repetition counter and cycle the loop without increasing
            ! the number of the current time step.
            
            irepetition = irepetition + 1

            ! Restore the flow situation to the beginning of the macrostep.
            ! Use the newly calculated time step size for the next time step.
            dtmp = rtimeStepping%dtstepFixed
            CALL c2d2_restoreTimestep (rsnapshotLastMacrostep,rproblem,&
                rtimeStepping,rvector,rrhs)
            CALL timstp_setBaseSteplength (rtimeStepping, dtmp)
            
            CALL output_line ('Repeating macrostep. Returning to timestep ' &
                //TRIM(sys_siL(rproblem%rtimedependence%itimeStep,6))//'.')
            CALL output_separator(OU_SEP_AT)
            
            ! Repeat the time step            
            CYCLE
              
          ELSE
          
            ! No repetitions left. Reset the repetition counter and
            ! continue with the next timestep, although the
            ! solution might be not really meaningful anymore.
            irepetition = 0
            
            CALL output_line ('No repetitions left. Cannot repeat macrostep. ' &
                //'Continuing with the next one...')
            
          END IF
          
        END SELECT
      
      ELSE
      
        ! No, time step is ok. If this is the last time step of the macrostep,
        ! reset the repetition counter such that it starts with 0 for the
        ! next mactostep.
        IF (MOD(rproblem%rtimedependence%itimeStep,3) .EQ. 0) THEN
          irepetition = 0
        END IF
      
      END IF
           
      rproblem%rtimedependence%itimeStep = rproblem%rtimedependence%itimeStep + 1

      CALL output_separator(OU_SEP_AT)

    END DO

    ! Clean up the stuff of/for the nonlinear solver.
    !
    ! Release the temporary vectors
    IF (rpredictedSolution%NEQ .NE. 0) &
      CALL lsysbl_releaseVector (rpredictedSolution)
    IF (roldSolution%NEQ .NE. 0) &
      CALL lsysbl_releaseVector (roldSolution)
    CALL lsysbl_releaseVector (rtempBlock2)
    CALL lsysbl_releaseVector (rtempBlock1)
    
    ! Release existing snapshots
    CALL c2d2_releaseSnapshop (rsnapshotLastMacrostep)
    
    ! Release parameters of the nonlinear loop
    CALL c2d2_doneNonlinearLoop (rnonlinearIteration)
             
    ! Release the preconditioner
    CALL c2d2_releasePreconditioner (rnonlinearIteration)
    
  END SUBROUTINE
  
END MODULE
