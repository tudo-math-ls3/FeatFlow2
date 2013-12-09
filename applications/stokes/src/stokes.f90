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

  use stokes2d_method0_simple
  use stokes2d_method0_resort
  use stokes2d_method1_mg
  use stokes2d_method1_schur
  use stokes2d_method2_sv
  use stokes2d_method2_gv
  use stokes2d_method3_block
  
  use stokes3d_method0_simple
  use stokes3d_method1_mg
  
  use navst2d_method1_mg
  use navst3d_method1_mg
  
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

  ! Call the problem to solve. 2d stokes 0:
  call output_lbrk()
  call output_line("Calculating 2D Stokes-Problem 0 - simple")
  call output_line("----------------------------------------")
  call stokes2d_0_simple

  ! Call the problem to solve. 2d stokes 0:
  call output_lbrk()
  call output_line("Calculating 2D Stokes-Problem 0 - CM-sorting")
  call output_line("--------------------------------------------")
  call stokes2d_0_resort

  ! Call the problem to solve. 2d stokes 1:
  call output_lbrk()
  call output_line("Calculating 2D Stokes-Problem 1 - multigrid")
  call output_line("-------------------------------------------")
  call stokes2d_1_mg

  ! Call the problem to solve. 2d stokes 1:
  call output_lbrk()
  call output_line("Calculating 2D Stokes-Problem 1 - Schur-Complement")
  call output_line("--------------------------------------------------")
  call stokes2d_1_schur

  ! Call the problem to solve. 2d stokes 2:
  call output_lbrk()
  call output_line("Calculating 2D Stokes-Problem 2 - simple VANKA")
  call output_line("----------------------------------------------")
  call stokes2d_2_sv

  ! Call the problem to solve. 2d stokes 2:
  call output_lbrk()
  call output_line("Calculating 2D Stokes-Problem 2 - general VANKA")
  call output_line("-----------------------------------------------")
  call stokes2d_2_gv

  ! Call the problem to solve. 2d stokes 2:
  call output_lbrk()
  call output_line("Calculating 2D Stokes-Problem 3 - general VANKA, block-assembly")
  call output_line("---------------------------------------------------------------")
  call stokes2d_3_block

  ! Call the problem to solve. 3d stokes 0:
  call output_lbrk()
  call output_line("Calculating 3D Stokes-Problem 0 - simple")
  call output_line("----------------------------------------")
  call stokes3d_0_simple

  ! Call the problem to solve. 3d stokes 1:
  call output_lbrk()
  call output_line("Calculating 3D Stokes-Problem 1 - multigrid")
  call output_line("-------------------------------------------")
  call stokes3d_1_mg

  ! As the following two examples are Navier-Stokes examples rather than
  ! Stokes, they are commented out by default.

  ! Call the problem to solve. 2d navier-stokes 1:
  ! call output_lbrk()
  ! call output_line("Calculating 2D Navier-Stokes-Problem 1 - multigrid")
  ! call output_line("--------------------------------------------------")
  ! call navst2d_1_mg

  ! Call the problem to solve. 3d navier-stokes 1:
  ! call output_lbrk()
  ! call output_line("Calculating 3D Navier-Stokes-Problem 1 - multigrid")
  ! call output_line("--------------------------------------------------")
  ! call navst3d_1_mg

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display "Handles in use=0" and "Memory in use=0"!
  call output_lbrk()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
