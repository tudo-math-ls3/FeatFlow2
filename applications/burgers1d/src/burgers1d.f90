!##############################################################################
!# ****************************************************************************
!# <name> burgers1d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the equation
!#
!#    $$  u_t  +  u*u_x  -  \nu u_xx  =  0,   u(x,0) = sin(Pi*x)  $$
!#
!# on a 2D domain $\Omega=[0,1] \times [0,1]$ for a scalar function $u(x,t)$.
!# </purpose>
!##############################################################################

program burgers1d

  use burgers1d_method5
  use burgers1d_method6
  
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
  ! 2.) Initialise FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Call the problem to solve.
  print *
  print *,'Calculating Burgers1D-Problem with method 5'
  print *,'-------------------------------------------'
  call burgers1d5

  ! Call the problem to solve.
  print *
  print *,'Calculating Burgers1D-Problem with method 6'
  print *,'-------------------------------------------'
  call burgers1d6

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk();
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call output_done()
  call storage_done()
  
end program
