!##############################################################################
!# ****************************************************************************
!# <name> codire </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the equation
!#
!#     - div ( alpha*grad(u) )  +  beta*grad(u)  +  gamma*u  =  f
!#
!# on a 2D domain for a scalar function u.
!#
!# The main example (module codire_method3) separates the different stages
!# of the solution process into different subroutines like example 2.
!# The communication is done using a problem-related structure. For the
!# communication with callback routines during the assembly, a
!# collection structure is set up.
!#
!# The second example (module codire_method3) shows how to use a Multigrid
!# solver to solve the problem. Besides that, the module is similar to
!# codire_method3.
!#
!# The equation itself is set up using parameters of a .DAT file
!# named 'data/codire.dat'. Here all parameters (alpha, beta, gamma, level,...)
!# can be found and changed by the user without recompiling the program.
!# </purpose>
!##############################################################################

program codire

  use codire_method3
  use codire_method5
  
  implicit none
  
  ! local variables
  character(len=SYS_STRLEN) :: slogdir,slogfile

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
  ! Initialise the storage management:
  !
  ! Initialise FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Call the problem to solve.
  print *
  print *,'Calculating CoDiRe-Problem with method 3'
  print *,'----------------------------------------'
  call codire3

  ! Call the problem to solve - method 5 = multigrid
  print *
  print *,'Calculating CoDiRe-Problem with method 5'
  print *,'----------------------------------------'
  call codire5
  
  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call output_done()
  call storage_done()
  
end program
