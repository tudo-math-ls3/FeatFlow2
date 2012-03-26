!##############################################################################
!# ****************************************************************************
!# <name> DG </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the equation
!#
!#          u_t + div f(u) = 0
!#
!# discontinuous Galerkin method.
!#
!# </purpose>
!##############################################################################

program dg

  use dg2d_method0_simple
  use dg2d_systems
  use dg2d_multigridscalar
  use dg2d_multigridsystem
  use dg2d_multigridpoisson
  use dg2d_multigridscalaradvectiondeal
  use dg2d_mgscpoissonde

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
  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)

!  ! Call the problem to solve. Scalar hyperbolic conservation law:
!  call output_lbrk ()
!  call output_line ('Calculating scalar hyperbolic conservation problem with DG')
!  call output_line ('----------------------------------------------------------')
!  call dg2d_0_simple
  
!  ! Call the problem to solve. Hyperbolic conservation law, system case:
!  call output_lbrk ()
!  call output_line ('Calculating hyperbolic system conservation problem with DG')
!  call output_line ('----------------------------------------------------------')
!  call dg2d_sys
  
!  ! Call the problem to solve. Linear scalar equation, multigrid solver:
!  call output_lbrk ()
!  call output_line ('Linear scalar equation, multigrid solver')
!  call output_line ('----------------------------------------')
!  call dg2d_mgsc
  
!  ! Call the problem to solve. Hyperbolic system conservation problem, multigrid solver:
!  call output_lbrk ()
!  call output_line ('Calculating hyperbolic system conservation problem')
!  call output_line ('with implicit DG discretisation and multigrid solver')
!  call output_line ('----------------------------------------------------')
!  call dg2d_mgsys

!  ! Call the problem to solve. Poisson equation, multigrid solver:
!  call output_lbrk ()
!  call output_line ('Poisson equation, multigrid solver')
!  call output_line ('----------------------------------')
!  call dg2d_mgpoisson
  
!  ! Call the problem to solve. Linear convection, as in DealII:
!  call output_lbrk ()
!  call output_line ('Linear convection equation, multigrid solver')
!  call output_line ('--------------------------------------------')
!  call dg2d_multigridscalaradvectiondeal_run

  ! Call the problem to solve. Poisson equation, as in DealII:
  call output_lbrk ()
  call output_line ('Poisson equation, multigrid solver')
  call output_line ('--------------------------------------------')
  call dg2d_dgmpd_run

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk ()
  call storage_info(.true.)

  ! Clean up the storage management, finish
  call storage_done()
  
  pause

end program
