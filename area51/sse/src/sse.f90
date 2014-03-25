!##############################################################################
!# ****************************************************************************
!# <name> sse </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the elliptic
!# equation for sea surface elevation
!#
!#    - \Nabla \dot (A \Nabla N) - i\omega N = 0   in \Omega
!#
!#                                         N = N_D on \Gamma_D
!#
!#                     - (A \Nabla N)\cdot n = 0   on \Gamma_N
!#
!# for a scalar complex function N.
!#
!# </purpose>
!##############################################################################

program sse

  use collection
  use convergencetable
  use fsystem
  use genoutput
  use linearsystemblock
  use paramlist
  use statistics
  use storage

  use sse_main

  ! local variables
  type(t_timer) :: rtimerTria,rtimerDiscr,rtimerMatVec,&
                   rtimerBC,rtimerSolver,rtimerPostproc,rtimerFree
  type(t_parlist) :: rparlist
  type(t_problem) :: rproblem
  type(t_convergenceTable) :: rtable
  character(len=SYS_STRLEN) :: slogdir,slogfile
  integer :: NLMIN,NLMAX,i

  ! Initialise system-wide settings
  call sys_init()

  ! Initialise the output system.
  !
  ! Normally, we write all the output to the screen and to a file
  ! "./log/output.txt".
  ! In the case that environment variables "$logdir"/"$resultsfile" exists,
  ! we write all the output to that file. This can be used e.g. in
  ! regression tests to compare results to reference results.
  if (sys_getenv_string("LOGDIR",slogdir) .and. &
      sys_getenv_string("RESULTFILE",slogfile)) then
    call output_init(trim(slogdir)//"/"//trim(slogfile))
  else
    call output_init("./log/output.txt")
  end if
  
  ! Initialise the FEAT 2.0 storage management
  call storage_init(999, 100)

  ! Read parameter file
  call parlst_init(rparlist)
  call parlst_readfromfile(rparlist, 'data/master.dat', 'data', .true.)

  ! Get levels from parameter list
  call parlst_getvalue_int(rparlist, '', 'NLMIN', NLMIN)
  call parlst_getvalue_int(rparlist, '', 'NLMAX', NLMAX)
  call parlst_getvalue_int(rparlist, '', 'PROBLEMTYPE',&
      rproblem%cproblemtype, rproblem%cproblemtype)

  ! Initialise the collection
  call collct_init (rproblem%rcollection)
  do i=1,NLMAX
    call collct_addlevel (rproblem%rcollection)
  end do

  ! Initialise the convergence table
  call ctab_init(rtable)
  
  do i=NLMIN,NLMAX

    ! Initialisation
    call output_line('Initialising triangulation...')
    call stat_clearTimer(rtimerTria)
    call stat_startTimer(rtimerTria,STAT_TIMERSHORT)
    call sse_initParamTriang(rparlist,NLMIN,i,rproblem)
    call stat_stopTimer(rtimerTria)
    call output_line(&
        '............................................................'//&
        trim(sys_sdEL(rtimerTria%delapsedReal,3))//'sec')

    call output_line('Initialising discretisation...')
    call stat_clearTimer(rtimerDiscr)
    call stat_startTimer(rtimerDiscr,STAT_TIMERSHORT)
    call sse_initDiscretisation(rparlist,rproblem)
    call stat_stopTimer(rtimerDiscr)
    call output_line(&
        '............................................................'//&
        trim(sys_sdEL(rtimerDiscr%delapsedReal,3))//'sec')
    
    call output_line('Initialising matrices/vectors...')
    call stat_clearTimer(rtimerMatVec)
    call stat_startTimer(rtimerMatVec,STAT_TIMERSHORT)
    call sse_initMatVec(rparlist,rproblem)
    call stat_stopTimer(rtimerMatVec)
    call output_line(&
        '............................................................'//&
        trim(sys_sdEL(rtimerMatVec%delapsedReal,3))//'sec')
    
    call output_line('Initialising/implementing discrete boundary conditions...')
    call stat_clearTimer(rtimerBC)
    call stat_startTimer(rtimerBC,STAT_TIMERSHORT)
    call sse_initDiscreteBC(rproblem)
    
    ! Implementation of boundary conditions
    call sse_implementBC(rproblem)
    call stat_stopTimer(rtimerBC)
    call output_line(&
        '............................................................'//&
        trim(sys_sdEL(rtimerBC%delapsedReal,3))//'sec')
    
    ! Solve the problem
    call output_line('Solving problem...')
    call stat_clearTimer(rtimerSolver)
    call stat_startTimer(rtimerSolver,STAT_TIMERSHORT)
    call sse_solve(rproblem)
    call stat_stopTimer(rtimerSolver)
    call output_line(&
        '............................................................'//&
        trim(sys_sdEL(rtimerSolver%delapsedReal,3))//'sec')

    ! Post-processing
    call output_line('Postprocessing solution...')
    call stat_clearTimer(rtimerPostproc)
    call stat_startTimer(rtimerPostproc,STAT_TIMERSHORT)
    call sse_postprocessing(rproblem,rtable)
    call stat_stopTimer(rtimerPostproc)
    call output_line(&
        '............................................................'//&
        trim(sys_sdEL(rtimerPostproc%delapsedReal,3))//'sec')
    
    ! Cleanup
    call output_line('Freeing memory...')
    call stat_clearTimer(rtimerFree)
    call stat_startTimer(rtimerFree,STAT_TIMERSHORT)
    call sse_doneMatVec(rproblem)
    call sse_doneBC(rproblem)
    call sse_doneDiscretisation(rproblem)
    call sse_doneParamTriang(rproblem)
    call stat_stopTimer(rtimerFree)
    call output_line(&
        '............................................................'//&
        trim(sys_sdEL(rtimerFree%delapsedReal,3))//'sec')
  end do
  
  ! Export convergence table
  call sse_outputTable(rproblem,rtable)

  ! Clear convergence table
  call ctab_done(rtable)

  ! Clear parameter list
  call parlst_done(rparlist)
  
  ! Print out heap statistics
  call output_lbrk()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program sse
