!##############################################################################
!# ****************************************************************************
!# <name> elasticity </name>
!# ****************************************************************************
!#
!# <purpose>
!#   Program for computing elasticity problems.
!#   This module is only responsible for choosing the desired application (see below).
!# </purpose>
!##############################################################################

program elasticity
   
  use elasticity_2d_disp_smallDeform_static
  
  implicit none
  
  ! local variables
  character(len=SYS_STRLEN) :: slogdir,slogfile
  
  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  call system_init()
  
  ! Initialise the output system.
  !
  ! Normally, we write all the output to the screen and to a file './log/output.txt'.
  ! In the case that environment variables "$logdir"/"$resultsfile" exists,
  ! we write all the output to that file. This can be used e.g. in
  ! regression tests to compare results to reference results.
  if (sys_getenv_string('LOGDIR',slogdir) .and. &
      sys_getenv_string('RESULTFILE',slogfile)) then
    call output_init(trim(slogdir)//'/'//trim(slogfile))
  else
    call output_init('./log/output.txt')
  end if

  ! The very second thing in every program: 
  ! Initialise the FEAT 2.0 storage management: 
  call storage_init(999, 100)

  ! Call the problem to solve.
  call output_lbrk()
  call output_line('-----------------------------------------------------')
  call output_line('Calculating elasticity problem:')
  call output_line(' * 2D')
  call output_line(' * pure displacements')
  call output_line(' * small deformation')
  call output_line(' * static')
  call output_line('-----------------------------------------------------')
  call elasticity_2d_disp_smallDeform_static_run

  ! further applications to follow...
  ! (transient, finite deformation, mixed formulation, ...)
  
  ! Print out heap statistics - just to check if everything is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish!
  call storage_done()
  
end program

