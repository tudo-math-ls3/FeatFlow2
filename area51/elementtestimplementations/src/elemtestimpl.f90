!##############################################################################
!# ****************************************************************************
!# <name> stokes </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the Stokes equation
!#
!#              $$- \nu Laplace(u) + \nabla p = f $$
!#              $$ \nabla \cdot u = 0$$
!#
!# on a 2D domain for a 2D function $u=(u_1,u_2)$ or
!# on a 3D domain for a 3D function $u=(u_1,u_2,u_3)$.
!#
!# A linear block system with 3 (for 2D) or 4 (for 3D) solution components
!# is set up, discretised and solved with multigrid.
!#
!# This program has a set of 2D (stokes2d_XXXX) and 3D (stokes3d_XXXX) stokes
!# examples. The naming convention is based on the one from the poisson
!# examples: stokesXd_method0_simple is the most simple stokes example using
!# a single-grid solver. All other examples use multi-grid solvers.
!# The stokesXd_method1_XXXX modules extend the stokesXd_method0_simple example
!# to more general situations, but still in the same style of a one-in-all
!# solver. The stokesXd_method2_XXXX modules decompose the different tasks
!# of the problem into different subroutines, thus bringing more structure
!# into the code and highlighting more and more features of the kernel.
!# </purpose>
!##############################################################################

program stokes

  use stokes2d_EL_QPW4P1TVDF_2D
  
  implicit none
  
  ! local variables
  character(len=SYS_STRLEN) :: slogdir,slogfile
  
  ! The very first thing in every application:
  ! Initialise system-wide settings:
  
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
    call output_init (trim(slogdir)//"/"//trim(slogfile))
  else
    call output_init ("./log/output.txt")
  end if

  ! The very second thing in every program:
  ! Initialise the storage management:
  !
  ! 2.) Initialise FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Call the problem to solve.
  call output_lbrk()
  call output_line("Calculating 2D Stokes-Problem 3 - Divergence-free element")
  call output_line("---------------------------------------------------------")
  call test_stokes2d_EL_QPW4P1TVDF_2D

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display "Handles in use=0" and "Memory in use=0"!
  call output_lbrk()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
