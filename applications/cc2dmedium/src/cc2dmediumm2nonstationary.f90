!##############################################################################
!# ****************************************************************************
!# <name> cc2dmmediumm2timeerror </name>
!# ****************************************************************************
!#
!# <purpose>
!# 
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
  
  USE collection
  USE convection
    
  USE cc2dmediumm2basic
  USE cc2dmedium_callback

  USE cc2dmediumm2stationary
  USE adaptivetimestep
  USE cc2dmediumm2timeerror
  USE cc2dmediumm2boundary
  USE cc2dmediumm2discretisation
  USE cc2dmediumm2postprocessing
  
    
  IMPLICIT NONE

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

  SUBROUTINE c2d2_backupTimestep ()

!</subroutine>
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_restoreTimestep ()

!</subroutine>
  
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
    ! switch to the next point in time.
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
    
    ! Save the core equation structure to the collection to inform
    ! the callback routines how the core equation looks like.
    CALL c2d2_saveNonlinearLoop (rnonlinearIterationTmp,rproblem%rcollection)
    
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
    ! Call the nonlinear solver. For preconditioning and defect calculation,
    ! the solver calls our callback routines.
    CALL nlsol_performSolve(rnlSolver,rvector,rtempVectorRhs,rtempVector,&
                            c2d2_getDefect,c2d2_precondDefect,c2d2_resNormCheck,&
                            rcollection=rproblem%rcollection)

    ! scale the pressure back, then we have again the correct solution vector.
    CALL lsyssc_scaleVector (rvector%RvectorBlock(NDIM2D+1),&
        1.0_DP/rtimestepping%dtstep)

    ! rvector is the solution vector u^{n+1}.    
    !
    ! Clean up the temporary core equation. Only remove the elements from
    ! the collection, i.e. don't clean up the structure as it's a copy
    ! of the one of our parent -- we mustn't destroy that!
    CALL c2d2_doneNonlinearLoop (rnonlinearIterationTmp,rproblem%rcollection,.FALSE.)
  
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
    TYPE(t_vectorBlock) :: rtempBlock1,rtempBlock2
    TYPE(t_ccNonlinearIteration) :: rnonlinearIteration
    
    ! Configuration block of the time stepping scheme 
    TYPE(t_explicitTimeStepping)        :: rtimestepping,rtimesteppingPredictor

    ! The nonlinear solver configuration
    TYPE(t_nlsolNode) :: rnlSol

    INTEGER :: iiteration
    REAL(DP) :: dtimeerror

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
    CALL c2d2_checkAssembly (rproblem,rrhs,rnonlinearIteration)
    
    ! Using rfinalAssembly as computed above, make the matrices compatible 
    ! to our preconditioner if they are not.
    CALL c2d2_finaliseMatrices (rnonlinearIteration)
    
    ! Initialise the preconditioner for the nonlinear iteration
    CALL c2d2_preparePreconditioner (rproblem,&
        rnonlinearIteration%rpreconditioner,rvector,rrhs)

    ! Create temporary vectors we need for the nonlinear iteration.
    CALL lsysbl_createVecBlockIndirect (rrhs, rtempBlock1, .FALSE.)
    CALL lsysbl_createVecBlockIndirect (rrhs, rtempBlock2, .FALSE.)

    ! Implement the initial boundary conditions into the solution vector.
    ! Don't implement anything to matrices or RHS vector as these are
    ! maintained in the timeloop.
    ! Afterwards, we can start the timeloop.
    CALL c2d2_implementBC (rproblem,rvector,rrhs,.FALSE.,.TRUE.,.FALSE.)

    ! Let's start the timeloop
  
    iiteration = 1
    dtimeerror = rproblem%rtimedependence%dminTimeDerivative
    
    rproblem%rtimedependence%dtime = rproblem%rtimedependence%dtimeInit
    
    DO WHILE ((iiteration .LE. rproblem%rtimedependence%niterations) .AND. &
        (rproblem%rtimedependence%dtime .LT. rproblem%rtimedependence%dtimemax) .AND. &
        (dtimeerror .GE. rproblem%rtimedependence%dminTimeDerivative))
              
      rproblem%rtimedependence%itimeStep = iiteration
      
      CALL output_separator(OU_SEP_MINUS)
      CALL output_line ('Time step '//TRIM(sys_siL(iiteration,6))// &
                        '     Time '//TRIM(sys_sdL(rproblem%rtimedependence%dtime,5)))
      CALL output_lbrk ()
              
      ! Proceed to the next time step
      CALL c2d2_performTimestep (rproblem,rvector,rrhs,&
          rtimestepping,rnonlinearIteration,rnlSol,rtempBlock1,rtempBlock2)      
      
      ! Postprocessing. Write out the solution.
      CALL c2d2_postprocessingNonstat (rproblem,rvector)
           
      iiteration = iiteration+1
    END DO

    ! Clean up the stuff of/for the nonlinear solver.
    !
    ! Release the temporary vectors
    CALL lsysbl_releaseVector (rtempBlock2)
    CALL lsysbl_releaseVector (rtempBlock1)
    
    ! Release parameters of the nonlinear loop
    CALL c2d2_doneNonlinearLoop (rnonlinearIteration)
             
    ! Release the preconditioner
    CALL c2d2_releasePreconditioner (rproblem,rnonlinearIteration%rpreconditioner)
    
  END SUBROUTINE
  
END MODULE
