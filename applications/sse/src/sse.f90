!##############################################################################
!# ****************************************************************************
!# <name> sse </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program to solve several types of
!# elliptic problems:
!#
!# 1) Poisson problem
!#
!#         -div(grad u) = f   in Omega
!#
!#                    u = u_D on Gamma_D
!#
!#                du/dn = g   on Gamma_N
!#
!# This equation can be solved either as is or cast into a first-order
!# system, whereby an inf-sup stable finite element pair is used to
!# approximate $u$ and its gradient $sigma=grad u$.
!#
!#          -div(sigma) = f   in Gamma
!#
!#       sigma - grad u = 0   in Gamma
!#
!#                    u = u_D on Gamma_D
!#
!#          d(sigma)/dn = g   on Gamma_N
!#
!#
!# 2) Sea surface elevation problem
!#
!#    div( A(0) * (grad N)) + i * omega * N = 0   in Omega
!#
!#                                        N = N_D on Gamma_D
!#
!#                             d(A(0)*N)/dn = 0   on Gamma_N
!#
!# for a scalar complex function $N$. This equation can be solved
!# either as is with high-order finite elements or cast into a
!# first-order system, whereby an inf-sup stable finite element pair
!# is used to approximate the sea surface elevation $N$ and its
!# gradient $\nabla N$. Since Feat2 does not support complex numbers
!# the unknown solution and the complex-valued system matrices are
!# split into their real and imaginary parts and solved as a coupled
!# but real-valued problem.
!#
!#
!# 3) Corin'e problem
!#
!#    TBA...
!#
!# </purpose>
!##############################################################################

program sse

  use collection
  use convergencetable
  use fsystem
  use genoutput
!!$  use linearsystemblock
  use paramlist
  use statistics
  use storage

  use sse_main
  use sse_base
  use sse_base_poisson
  use sse_base_sse
  use sse_base_corine  

  ! local variables
  type(t_timer) :: rtimerTria,rtimerDiscr,rtimerMatVec,&
                   rtimerBC,rtimerSolver,rtimerPostproc,rtimerFree
  type(t_parlist) :: rparlist
  type(t_problem) :: rproblem
  type(t_convergenceTable) :: rtable
  character(len=SYS_STRLEN) :: slogdir,slogfile
  integer :: ILMIN,NLMIN,NLMAX,i

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
  call parlst_getvalue_int(rparlist, '', 'ILMIN', ILMIN, NLMIN)
  call parlst_getvalue_int(rparlist, '', 'PROBLEMTYPE', rproblem%cproblemtype)
  call parlst_getvalue_int(rparlist, sse_getSection(rproblem%cproblemtype),&
                                         'PROBLEMSUBTYPE', rproblem%cproblemsubtype)
  
  ! Write configuration to screen
  call output_separator (OU_SEP_STAR)
  
  select case(rproblem%cproblemtype)
  case (POISSON_SCALAR, POISSON_SYSTEM)
    call sse_initParamPoisson(rproblem%cproblemtype,rparlist)
    call sse_infoPoisson(rproblem%cproblemtype,rproblem%cproblemsubtype)
    
  case (SSE_SCALAR, SSE_SYSTEM1, SSE_SYSTEM2)
    call sse_initParamSSE(rproblem%cproblemtype,rparlist)
    call sse_infoSSE(rproblem%cproblemtype,rproblem%cproblemsubtype)
    
  case (CORINE_1D, CORINE_2D)
    call sse_initParamCorine(rproblem%cproblemtype,rparlist)
    call sse_infoCorine(rproblem%cproblemtype,rproblem%cproblemsubtype)

  case default
    call output_line("Invalid problem type", &
        OU_CLASS_ERROR,OU_MODE_STD,"sse")
    call sys_halt()
  end select
  
  call output_separator (OU_SEP_STAR)

  ! Initialise the collection
  call collct_init (rproblem%rcollection)
  do i=1,NLMAX
    call collct_addlevel (rproblem%rcollection)
  end do
  
  ! Initialise the convergence table
  call ctab_init(rtable)
  
  do i=ILMIN,NLMAX

    ! Initialisation
    call output_line('Initialising triangulation')
    call stat_clearTimer(rtimerTria)
    call stat_startTimer(rtimerTria,STAT_TIMERSHORT)
    call sse_initParamTriang(rparlist,NLMIN,i,rproblem)
    call stat_stopTimer(rtimerTria)
    call output_line(&
        '............................................................'//&
        trim(sys_sdEL(rtimerTria%delapsedReal,3))//'sec')
    
    call output_line('Initialising discretisation')
    call stat_clearTimer(rtimerDiscr)
    call stat_startTimer(rtimerDiscr,STAT_TIMERSHORT)
    call sse_initDiscretisation(rparlist,rproblem)
    call stat_stopTimer(rtimerDiscr)
    call output_line(&
        '............................................................'//&
        trim(sys_sdEL(rtimerDiscr%delapsedReal,3))//'sec')

    call output_line('Initialising matrices/vectors')
    call stat_clearTimer(rtimerMatVec)
    call stat_startTimer(rtimerMatVec,STAT_TIMERSHORT)
    call sse_initMatVec(rparlist,rproblem)
    call stat_stopTimer(rtimerMatVec)
    call output_line(&
        '............................................................'//&
        trim(sys_sdEL(rtimerMatVec%delapsedReal,3))//'sec')

    call output_line('Initialising/implementing discrete boundary conditions')
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
    call output_line('Solving problem')
    call stat_clearTimer(rtimerSolver)
    call stat_startTimer(rtimerSolver,STAT_TIMERSHORT)
    call sse_solve(rparlist,rproblem)
    call stat_stopTimer(rtimerSolver)
    call output_line(&
        '............................................................'//&
        trim(sys_sdEL(rtimerSolver%delapsedReal,3))//'sec')

    ! Post-processing
    call output_line('Postprocessing solution')
    call stat_clearTimer(rtimerPostproc)
    call stat_startTimer(rtimerPostproc,STAT_TIMERSHORT)
    call sse_postprocessing(rparlist,rproblem,rtable)
    call stat_stopTimer(rtimerPostproc)
    call output_line(&
        '............................................................'//&
        trim(sys_sdEL(rtimerPostproc%delapsedReal,3))//'sec')

    ! Cleanup
    call output_line('Freeing memory')
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
