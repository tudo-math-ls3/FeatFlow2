!##############################################################################
!# ****************************************************************************
!# <name> nonstationaryoptcsolver </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises a fully coupled time dependent solver for the
!# coupled Navier-Stokes optimal control problem.
!#
!# The following routines can be found here:
!#
!# 1.) cc_initParTimeDependence
!#     -> Initialise the parameters of the time dependent solver from DAT file
!#        parameters.
!#
!# 2.) cc_solveNonstationaryDirect
!#     -> Solves the space-time coupled system.
!# </purpose>
!##############################################################################

module nonstationaryoptcsolver

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
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
    
  use basicstructures
  use user_callback

  use spacepreconditioner
  use spacepreconditionerinit
  use stationaryoptcsolver
  use timeanalysis
  use spatialbc
  use spacediscretisation
  use postprocessing
  use nonlinearoneshotspacetimesolver
  use spacediscretisation
    
  implicit none

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
  type(t_parlist), intent(IN) :: rparams
  
  ! The name of the section in the parameter list containing the parameters
  ! for the time dependent simulation.
  character(LEN=*), intent(IN) :: ssection
!</input>

!<inputoutput>
  ! The problem structure to be initialised.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>
  
!</subroutine>

    ! Fetch the parameters. Initialise with standard settings if they don't
    ! exist.
    call parlst_getvalue_int (rparams,ssection,'itimedependence',  &
        rproblem%itimedependence, 0)
    call parlst_getvalue_int (rparams,ssection,'niterations',     &
        rproblem%rtimedependence%niterations, 1000)
    call parlst_getvalue_double (rparams,ssection,'dtimestart',   &
        rproblem%rtimedependence%dtimeInit, 0.0_DP)
    call parlst_getvalue_double (rparams,ssection,'dtimemax',     &
        rproblem%rtimedependence%dtimeMax, 20.0_DP)
    call parlst_getvalue_int (rparams, ssection,'ctimeStepScheme', &
          rproblem%rtimedependence%ctimeStepScheme, 0)
    call parlst_getvalue_double (rparams,ssection,'dtimeStepTheta',    &
        rproblem%rtimedependence%dtimeStepTheta, 1.0_DP)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine cc_solveNonstationaryDirect (rproblem,rvector,ipropagateStartVector)
  
!<description>
  ! Solve the nonstationary optimal control problem. This allocates
  ! memory for the solution vector and solves the space-time coupled
  ! system.
  ! This is the 'single-grid' variant of the coupled solver, i.e.
  ! it calls the solver once without any multigrid schemes in space/time.
!</description>

!<input>
  ! OPTIONAL: A block vector that serves as stationary initial condition.
  type(t_vectorBlock), intent(IN), optional :: rvector
  
  ! OPTIONAL: If a start vector is present, this specifies if the start vector
  ! is propagated to all timesteps.
  integer, intent(in), optional :: ipropagateStartVector
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    type(t_ccoptSpaceTimeDiscretisation), dimension(:), pointer :: RspaceTimeDiscr
    type(t_spacetimeVector) :: rx,rd,rb
    integer :: i,ispacelevelcoupledtotimelevel,cspaceTimeSolverType
    integer :: nmaxSimulRefLevel,icurrentspace,icurrenttime
    real(DP) :: dspacetimeRefFactor,dcurrentfactor
    integer :: ctypePreconditioner
    integer(I32) :: TIMENLMIN,TIMENLMAX,istep
    integer, dimension(:,:), allocatable :: Ispacetimelevel

    call output_lbrk()

    ! Get the minimum and maximum time level from the parameter list
    call parlst_getvalue_int (rproblem%rparamList,'TIME-DISCRETISATION',&
                              'TIMENLMIN',TIMENLMIN,1)
    call parlst_getvalue_int (rproblem%rparamList,'TIME-DISCRETISATION',&
                              'TIMENLMAX',TIMENLMAX,1)
                              
    ! Values <= 0 result in a relativ minimum level.
    if (TIMENLMIN .le. 0) TIMENLMIN = TIMENLMAX+TIMENLMIN

    ! Allocate memory for the space-time discretisation structures
    ! on all levels. Initialis the maximum level such that it represents
    ! the finest time level as well as the finest space level.
    allocate(RspaceTimeDiscr(1:TIMENLMAX))
    
    call output_line ('Initialising space-time discretisation of the maximum level...')
    
    call sptidis_initDiscretisation (rproblem,TIMENLMAX,&
        rproblem%NLMAX,RspaceTimeDiscr(TIMENLMAX))
        
    call sptivec_initVector (rx,RspaceTimeDiscr(TIMENLMAX)%rtimeDiscr,&
        RspaceTimeDiscr(TIMENLMAX)%p_rlevelInfo%rdiscretisation)
    call sptivec_initVector (rb,RspaceTimeDiscr(TIMENLMAX)%rtimeDiscr,&
        RspaceTimeDiscr(TIMENLMAX)%p_rlevelInfo%rdiscretisation)
    call sptivec_initVector (rd,RspaceTimeDiscr(TIMENLMAX)%rtimeDiscr,&
        RspaceTimeDiscr(TIMENLMAX)%p_rlevelInfo%rdiscretisation)
        
    ! If a start vector is given, propagate it to all timesteps.
    if (present(rvector)) then
      ! Initialise the solution at timestep 1.
      call sptivec_setTimestepData(rx,1,rvector)
      
      ! Eventually propagate the initial solution to all timesteps.
      if (present(ipropagateStartVector)) then
        if (ipropagateStartVector .gt. 0) then
          call output_line ('Propagating start vector...')
          do istep = 2,rx%NEQtime
            call sptivec_setTimestepData(rx,istep,rvector)
          end do
        end if
      end if
    end if

    call output_line ('Reading target flow...')

    ! Read the target flow -- stationary or nonstationary
    call cc_initTargetFlow (rproblem,&
        RspaceTimeDiscr(TIMENLMAX)%rtimeDiscr%dtimeInit,&
        RspaceTimeDiscr(TIMENLMAX)%rtimeDiscr%dtimeMax,&
        RspaceTimeDiscr(TIMENLMAX)%rtimeDiscr%nintervals)
        
    ! Figure out which type of solver we should use to solve
    ! the problem.
    call parlst_getvalue_int (rproblem%rparamList,'TIME-SOLVER',&
                              'cspaceTimeSolverType',cspaceTimeSolverType,0)
                              
    ! Get the preconditioner type which is to be used in the nonlinear
    ! solver.
    call parlst_getvalue_int (rproblem%rparamList,'TIME-SOLVER',&
                              'ctypePreconditioner',ctypePreconditioner,1)

    ! Call the solver for the space/time coupled system.
    select case (cspaceTimeSolverType)
    case (0)
      ! 1-level (in time) Gauss elimination solver.
      call output_line ('Invoking Gauss elimination solver...')
      call cc_solveSupersystemDirect (rproblem, &
          RspaceTimeDiscr(TIMENLMAX), rx, rb, rd)

    case (1)
      ! 1-level (in time) defect correction solver with a preconditioner
      ! in space.
      !
      ! Call the defect correction solver
      call output_line ('Invoking defect correction solver...')
      call cc_solveSupersystemDefCorr (rproblem, &
          RspaceTimeDiscr(TIMENLMAX), rx, rb, rd, ctypePreconditioner)
      
    case (2)
      ! Space-time multigrid solver in space and/or time.
      !
      ! Should we couple space and time coarsening/refinement?
      call parlst_getvalue_int (rproblem%rparamList,'TIME-MULTIGRID',&
          'ispacelevelcoupledtotimelevel',ispacelevelcoupledtotimelevel,1)

      ! Parameters of the refinement?
      call parlst_getvalue_int (rproblem%rparamList,'TIME-MULTIGRID',&
          'nmaxSimulRefLevel',nmaxSimulRefLevel,0)
      
      if (nmaxSimulRefLevel .le. 0) then
        nmaxSimulRefLevel = TIMENLMAX+nmaxSimulRefLevel
      end if
      nmaxSimulRefLevel = min(TIMENLMAX,max(TIMENLMIN,nmaxSimulRefLevel))
          
      call parlst_getvalue_double (rproblem%rparamList,'TIME-MULTIGRID',&
          'dspacetimeRefFactor',dspacetimeRefFactor,1.0_DP)

      call output_line ('Initialising space-time discretisation of the lower levels...')
      
      ! Initialise the supersystem on all levels below the topmost one
      select case (ispacelevelcoupledtotimelevel)
      case (0)
        ! Everywhere the same space level
        do i=TIMENLMIN,TIMENLMAX-1
          call sptidis_initDiscretisation (rproblem,i,&
              rproblem%NLMAX,RspaceTimeDiscr(i))
        end do
        
      case (1)
        ! Space level NLMAX-i = Time level TIMENLMAX-i.
        !
        ! But this is modified by nmaxSimulRefLevel and dspacetimeRefFactor!
        ! From level nmaxSimulRefLevel on, we have to use the maximum
        ! space level. Below of that, on every space refinement we have
        ! to do dspacetimeRefFactor time refinements (rounded).
        !
        ! To get a formula for that, one might get crazy. So we calculate
        ! the levels in advance. Allocate an array Ispacetimelevel
        ! and fill it: Ispacetimelevel(1,:) all the space levels,
        ! Ispacetimelevel(2,:) all the corresponding time levels.
        allocate(Ispacetimelevel(2,TIMENLMAX))
        
        icurrentspace = rproblem%NLMAX
        icurrenttime = TIMENLMAX
        dcurrentfactor = 0.0_DP
        do i=TIMENLMAX, TIMENLMIN,-1
          Ispacetimelevel(1,i) = icurrentspace
          Ispacetimelevel(2,i) = icurrenttime
          
          ! Reduce space and/or time level.
          if (dspacetimeRefFactor .ge. 1.0_DP) then
            ! Time coarsened as usual, space probably delayed
            icurrenttime = icurrenttime - 1
            if (i .le. nmaxSimulRefLevel) then
              ! Space coarsening allowed.
              dcurrentfactor = dcurrentfactor + 1.0_DP/dspacetimeRefFactor
              if (dcurrentfactor .ge. 1.0_DP) then
                ! Coarsening in space
                if (icurrentspace .gt. rproblem%NLMIN) &
                    icurrentspace = icurrentspace - 1
                dcurrentfactor = mod(dcurrentfactor,1.0_DP)
              end if
              ! Otherwise no coarsening in space.
            end if
          else
            ! Space coarsened as usual, time probably delayed
            if (icurrentspace .gt. rproblem%NLMIN) &
                icurrentspace = icurrentspace - 1
            if (i .le. nmaxSimulRefLevel) then
              ! Time coarsening allowed.
              dcurrentfactor = dcurrentfactor + dspacetimeRefFactor
              if (dcurrentfactor .ge. 1.0_DP) then
                ! Coarsening in space
                if (icurrenttime .gt. TIMENLMIN) icurrenttime = icurrenttime - 1
                dcurrentfactor = mod(dcurrentfactor,1.0_DP)
              end if
              ! Otherwise no coarsening in time.
            end if
          end if
        end do
        
        ! Now initialise the level hierarchy
        do i=TIMENLMIN,TIMENLMAX-1
          !CALL sptidis_initDiscretisation (rproblem,i,&
          !    MAX(rproblem%NLMIN,rproblem%NLMAX-(TIMENLMAX-i)),&
          !    RspaceTimeDiscr(i))
          call sptidis_initDiscretisation (rproblem,Ispacetimelevel(2,i),&
              Ispacetimelevel(1,i),RspaceTimeDiscr(i))
        end do
        
        deallocate (Ispacetimelevel)
        
      case (2)
        ! Space level NLMAX-i = Time level TIMENLMAX-i,
        ! but no restriction in time
        do i=TIMENLMIN,TIMENLMAX-1
          call sptidis_initDiscretisation (rproblem,TIMENLMAX,&
              max(rproblem%NLMIN,rproblem%NLMAX-(TIMENLMAX-i)),&
              RspaceTimeDiscr(i))
        end do
      end select
      
      ! Print some statistical information about the discretisation
      call output_lbrk ()
      call output_line ('Triangulation')
      call output_line ('-------------')
      do i=rproblem%NLMIN,rproblem%NLMAX
        call tria_infoStatistics (rproblem%RlevelInfo(i)%rtriangulation,&
          i .eq. rproblem%NLMIN,i)
      end do
      
      call output_lbrk ()
      call output_line ('Time discretisation')
      call output_line ('-------------------')
      do i=TIMENLMIN,TIMENLMAX
        call output_line ('Level: '//sys_siL(i,10))
        call sptidis_infoDiscretisation (RspaceTimeDiscr(i))
        call output_lbrk ()
      end do
      
      ! Call the multigrid solver to solve on all these levels
      call output_line ('Invoking nonlinear space-time solver...')
      call cc_solveSupersystemMultigrid (rproblem, &
          RspaceTimeDiscr(TIMENLMIN:TIMENLMAX), rx, rb, rd, ctypePreconditioner)
    end select
    
    call output_line ('Postprocessing...')
    
    ! POSTPROCESSING
    !
    call cc_spacetimepostproc (rproblem,RspaceTimeDiscr(TIMENLMAX),rx)
    
    call output_line ('Cleaning up...')
    
    call cc_doneTargetFlow (rproblem)
    
    do i=TIMENLMIN,TIMENLMAX-1
      call sptidis_doneDiscretisation (RspaceTimeDiscr(i))
    end do
    
    deallocate(RspaceTimeDiscr)
    
    call sptivec_releaseVector (rb)
    call sptivec_releaseVector (rd)
    call sptivec_releaseVector (rx)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_spacetimepostproc (rproblem,rspacetimediscr,rvector)
  
!<description>
  ! Performs postprocessing of the space time solution vector rvector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!<input>
  ! The space-time solution vector on which postprocessing should be applied.
  type(t_spacetimevector), intent(IN) :: rvector
  
  ! The space-time discretisation structure belonging to that vector.
  type(t_ccoptSpaceTimeDiscretisation), intent(IN) :: rspaceTimeDiscr
!</input>

!</subroutine>

    ! local variables
    integer :: i
    character(SYS_STRLEN) :: stemp,sfileName
    type(t_vectorBlock) :: rvectorTmp
    real(dp) :: dtime
    
    ! Create a temp vector
    call lsysbl_createVecBlockByDiscr (&
        rspacetimediscr%p_rlevelInfo%rdiscretisation,rvectorTmp,.true.)

    ! Postprocessing of all solution vectors.
    do i = 0,tdiscr_igetNDofGlob(rvector%p_rtimeDiscretisation)-1
    
      dtime = &
          rvector%p_rtimeDiscretisation%dtimeInit + i * rvector%p_rtimeDiscretisation%dtstep
    
      ! Evaluate the space time function in rvector in the point
      ! in time dtime. Independent of the discretisation in time,
      ! this will give us a vector in space.
      !CALL sptivec_getTimestepData (rx, i, rvectorTmp)
      call tmevl_evaluate(rvector,dtime,rvectorTmp)
    
      call cc_postprocessingNonstat (rproblem,i,dtime,rvectorTmp)
      
    end do
    
    ! If there's a ´filename for the solution file given,
    ! write out the current solution.
    call parlst_getvalue_string (rproblem%rparamList, 'TIME-POSTPROCESSING', &
                                 'sfinalSolutionFileName', stemp, '')
    read (stemp,*) sfileName
    if (sfilename .ne. '') then
      stemp = '('''//trim(adjustl(sfilename))//'.'',I5.5)'
      call sptivec_saveToFileSequence(rvector,stemp,.true.)
    end if
    
    call lsysbl_releaseVector(rvectorTmp);
    
  end subroutine
  
end module
