!##############################################################################
!# ****************************************************************************
!# <name> renum2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program uses a simple test the equation
!#
!#              - Laplace(u) = f
!#
!# for a scalar function u to test the MFLOP rate in matrix-vector
!# multiplication and the speed of different solvers on different levels.
!# </purpose>
!##############################################################################

program poisson
   
  use renum2d_test1
  use renum2d_test2
  
  implicit none
  
  ! local variables
  character(len=SYS_STRLEN) :: slogdir,slogfile,sdatafile
  type(t_parlist) :: rparlist
  logical :: bexists
  
  ! The very first thing in every application:
  ! Initialise system-wide settings:
  
  call system_init()
  
  ! Initialise the output system.
  !
  ! Normally, we write all the output to the screen and to a file
  ! './log/output.txt'.
  ! In the case that environment variables "$logdir"/"$resultsfile" exists,
  ! we write all the output to that file. This can be used e.g. in
  ! regression tests to compare results to reference results.
  if (sys_getenv_string('LOGDIR',slogdir) .and. &
      sys_getenv_string('RESULTFILE',slogfile)) then
    call output_init (trim(slogdir)//'/'//trim(slogfile))
  else
    call output_init ('./log/output.txt')
  end if

  ! The very second thing in every program:
  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)
  
  ! Print the configuration
  call sys_getcommandLineArg(1,sdatafile,sdefault='./data/renum.dat')
  inquire(file=sdatafile,exist=bexists)
  if (.not. bexists) then
    call output_line ('Data file not found: '//trim(sdatafile), OU_CLASS_ERROR)
    call sys_halt()
  end if
  
  call parlst_init(rparlist)
  call parlst_readfromfile (rparlist, sdatafile)
  call parlst_info(rparlist)
  call parlst_done(rparlist)
 
  !!!!!!! poisson2d_0_simple = MFLOP-Rate test !!!!!!!!

  ! Call the problem to solve. Poisson 2D method 1 - simple:
  call output_lbrk ()
  call output_line ('Calculating Performance-test MatVec')
  call output_line ('-----------------------------------')
  call renum2d_matvec_mfloptest
  
  !!!!!!! poisson2d_0_simple = Solver test vor different renumbering strategies  !!!!!!!!
  
  ! Call the problem to solve. Poisson 2D method 1 - multigrid:
  call output_lbrk ()
  call output_line ('Calculating Performance-test SolverTest')
  call output_line ('---------------------------------------')
  call renum2d_solvertest

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk ()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
