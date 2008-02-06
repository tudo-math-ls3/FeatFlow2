!##############################################################################
!# ****************************************************************************
!# <name> cc2dmediumm2nonstationary </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises a fully coupled time dependent solver for the 
!# coupled Navier-Stokes optimal control problem.
!#
!# The following routines can be found here:
!#
!# 1.) c2d2_initParTimeDependence
!#     -> Initialise the parameters of the time dependent solver from DAT file
!#        parameters.
!#
!# 2.) c2d2_solveNonstationaryDirect
!#     -> Solves the space-time coupled system.
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
  USE cc2dmediumm2timeanalysis
  USE cc2dmediumm2boundary
  USE cc2dmediumm2discretisation
  USE cc2dmediumm2postprocessing
  USE cc2dmediumm2oneshot
  USE cc2dmediumm2discretisation
    
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
    CALL parlst_getvalue_int (rparams, ssection,'ctimeStepScheme', &
          rproblem%rtimedependence%ctimeStepScheme, 0)
    CALL parlst_getvalue_double (rparams,ssection,'dtimeStepTheta',    &
        rproblem%rtimedependence%dtimeStepTheta, 1.0_DP)
    rproblem%rtimedependence%itimeStep = 0

  END SUBROUTINE  

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_solveNonstationaryDirect (rproblem)
  
!<description>
  ! Solve the nonstationary optimal control problem. This allocates
  ! memory for the solution vector and solves the space-time coupled
  ! system.
  ! This is the 'single-grid' variant of the coupled solver, i.e.
  ! it calls the solver once without any multigrid schemes in space/time.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    TYPE(t_ccoptSpaceTimeDiscretisation), DIMENSION(:), ALLOCATABLE :: RspaceTimeDiscr
    TYPE(t_spacetimeVector) :: rx,rd,rb
    INTEGER :: i,ispacelevelcoupledtotimelevel,cspaceTimeSolverType
    INTEGER :: ctypePreconditioner
    INTEGER(I32) :: TIMENLMIN,TIMENLMAX

    ! Get the minimum and maximum time level from the parameter list    
    CALL parlst_getvalue_int (rproblem%rparamList,'TIME-DISCRETISATION',&
                              'TIMENLMIN',TIMENLMIN,1)
    CALL parlst_getvalue_int (rproblem%rparamList,'TIME-DISCRETISATION',&
                              'TIMENLMAX',TIMENLMAX,1)
                              
    ! Values <= 0 result in a relativ minimum level.
    IF (TIMENLMIN .LE. 0) TIMENLMIN = TIMENLMAX+TIMENLMIN

    ! Allocate memory fo rthe space-time discretisation structures
    ! on all levels. Initialis the maximum level such that it represents
    ! the finest time level as well as the finest space level.
    ALLOCATE(RspaceTimeDiscr(TIMENLMIN:TIMENLMAX))
    
    CALL sptidis_initDiscretisation (rproblem,TIMENLMAX,&
        rproblem%NLMAX,RspaceTimeDiscr(TIMENLMAX))
        
    CALL sptivec_initVector (rx,RspaceTimeDiscr(TIMENLMAX)%rtimeDiscr,&
        RspaceTimeDiscr(TIMENLMAX)%p_rlevelInfo%p_rdiscretisation)
    CALL sptivec_initVector (rb,RspaceTimeDiscr(TIMENLMAX)%rtimeDiscr,&
        RspaceTimeDiscr(TIMENLMAX)%p_rlevelInfo%p_rdiscretisation)
    CALL sptivec_initVector (rd,RspaceTimeDiscr(TIMENLMAX)%rtimeDiscr,&
        RspaceTimeDiscr(TIMENLMAX)%p_rlevelInfo%p_rdiscretisation)

    ! Read the target flow -- stationary or nonstationary
    CALL c2d2_initTargetFlow (rproblem,&
        RspaceTimeDiscr(TIMENLMAX)%rtimeDiscr%dtimeInit,&
        RspaceTimeDiscr(TIMENLMAX)%rtimeDiscr%dtimeMax,&
        RspaceTimeDiscr(TIMENLMAX)%rtimeDiscr%nintervals)
        
    ! Figure out which type of solver we should use to solve
    ! the problem. 
    CALL parlst_getvalue_int (rproblem%rparamList,'TIME-SOLVER',&
                              'cspaceTimeSolverType',cspaceTimeSolverType,0)
                              
    ! Get the preconditioner type which is to be used in the nonlinear
    ! solver.
    CALL parlst_getvalue_int (rproblem%rparamList,'TIME-SOLVER',&
                              'ctypePreconditioner',ctypePreconditioner,1)

    ! Call the solver for the space/time coupled system.
    SELECT CASE (cspaceTimeSolverType)
    CASE (0)
      ! 1-level (in time) Gauss elimination solver.
      CALL c2d2_solveSupersystemDirect (rproblem, &
          RspaceTimeDiscr(TIMENLMAX), rx, rb, rd)

    CASE (1)
      ! 1-level (in time) defect correction solver with a preconditioner
      ! in space. 
      !
      ! Call the defect correction solver      
      CALL c2d2_solveSupersystemDefCorr (rproblem, &
          RspaceTimeDiscr(TIMENLMAX), rx, rb, rd, ctypePreconditioner)
      
    CASE (2)
      ! Space-time multigrid solver in space and/or time.
      !
      ! Should we couple space and time coarsening/refinement?
      CALL parlst_getvalue_int (rproblem%rparamList,'TIME-MULTIGRID',&
          'ispacelevelcoupledtotimelevel',ispacelevelcoupledtotimelevel,1)

      ! Initialise the supersystem on all levels below the topmost one
      SELECT CASE (ispacelevelcoupledtotimelevel)
      CASE (0)
        ! Everywhere the same space level
        DO i=TIMENLMIN,TIMENLMAX-1
          CALL sptidis_initDiscretisation (rproblem,i,&
              rproblem%NLMAX,RspaceTimeDiscr(i))
        END DO
        
      CASE (1) 
        ! Space level NLMAX-i = Time level TIMENLMAX-i
        DO i=TIMENLMIN,TIMENLMAX-1
          CALL sptidis_initDiscretisation (rproblem,i,&
              MAX(rproblem%NLMIN,rproblem%NLMAX-(TIMENLMAX-i)),&
              RspaceTimeDiscr(i))
        END DO
        
      CASE (2)
        ! Space level NLMAX-i = Time level TIMENLMAX-i,
        ! but no restriction in time
        DO i=TIMENLMIN,TIMENLMAX-1
          CALL sptidis_initDiscretisation (rproblem,TIMENLMAX,&
              MAX(rproblem%NLMIN,rproblem%NLMAX-(TIMENLMAX-i)),&
              RspaceTimeDiscr(i))
        END DO
      END SELECT
      
      ! Print some statistical information about the discretisation
      CALL output_lbrk ()
      CALL output_line ('Triangulation')
      CALL output_line ('-------------')
      DO i=rproblem%NLMIN,rproblem%NLMAX
        CALL tria_infoStatistics (rproblem%RlevelInfo(i)%rtriangulation,&
          i .EQ. rproblem%NLMIN,i)
      END DO
      
      CALL output_lbrk ()
      CALL output_line ('Time discretisation')
      CALL output_line ('-------------------')
      DO i=TIMENLMIN,TIMENLMAX
        CALL output_line ('Level: '//sys_siL(i,10))
        CALL sptidis_infoDiscretisation (RspaceTimeDiscr(i))
        CALL output_lbrk ()
      END DO
      
      ! Call the multigrid solver to solve on all these levels
      CALL c2d2_solveSupersystemMultigrid (rproblem, &
          RspaceTimeDiscr(TIMENLMIN:TIMENLMAX), rx, rb, rd, ctypePreconditioner)
    END SELECT
    
    ! POSTPROCESSING
    !
    CALL c2d2_spacetimepostproc (rproblem,RspaceTimeDiscr(TIMENLMAX),rx)
    
    CALL c2d2_doneTargetFlow (rproblem)
    
    DO i=TIMENLMIN,TIMENLMAX-1
      CALL sptidis_doneDiscretisation (RspaceTimeDiscr(i))
    END DO
    CALL sptivec_releaseVector (rb)
    CALL sptivec_releaseVector (rd)
    CALL sptivec_releaseVector (rx)

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_spacetimepostproc (rproblem,rspacetimediscr,rvector)
  
!<description>
  ! Performs postprocessing of the space time solution vector rvector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!<input>
  ! The space-time solution vector on which postprocessing should be applied.
  TYPE(t_spacetimevector), INTENT(IN) :: rvector
  
  ! The space-time discretisation structure belonging to that vector.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr
!</input>

!</subroutine>

    ! local variables
    INTEGER :: i
    CHARACTER(SYS_STRLEN) :: stemp,sfileName
    TYPE(t_vectorBlock) :: rvectorTmp
    
    ! Create a temp vector
    CALL lsysbl_createVecBlockByDiscr (&
        rspacetimediscr%p_rlevelInfo%p_rdiscretisation,rvectorTmp,.TRUE.)

    ! Attach the boundary conditions to that vector
    rvectorTmp%p_rdiscreteBC => rspacetimediscr%p_rlevelInfo%p_rdiscreteBC
    rvectorTmp%p_rdiscreteBCfict => rspacetimediscr%p_rlevelInfo%p_rdiscreteFBC

    ! Postprocessing of all solution vectors.
    DO i = 0,tdiscr_igetNDofGlob(rvector%p_rtimeDiscretisation)-1
    
      rproblem%rtimedependence%dtime = &
          rvector%p_rtimeDiscretisation%dtimeInit + &
          i * rvector%p_rtimeDiscretisation%dtstep
      rproblem%rtimedependence%itimeStep = i
    
      ! Evaluate the space time function in rvector in the point
      ! in time dtime. Independent of the discretisation in time,
      ! this will give us a vector in space.
      !CALL sptivec_getTimestepData (rx, i, rvectorTmp)
      CALL tmevl_evaluate(rvector,rproblem%rtimedependence%dtime,rvectorTmp)
    
      CALL c2d2_postprocessingNonstat (rproblem,rvectorTmp)  
      
    END DO
    
    ! If there's a ´filename for the solution file given,
    ! write out the current solution.
    CALL parlst_getvalue_string_direct (rproblem%rparamList, 'TIME-POSTPROCESSING', &
                                        'sfinalSolutionFileName', stemp, '')
    READ (stemp,*) sfileName
    IF (sfilename .NE. '') THEN
      stemp = '('''//TRIM(ADJUSTL(sfilename))//'.'',I5.5)'
      CALL sptivec_saveToFileSequence(rvector,stemp,.TRUE.)
    END IF
    
    CALL lsysbl_releaseVector(rvectorTmp);
    
  END SUBROUTINE
  
END MODULE
