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

PROGRAM stokes

  USE stokes2d_method0_simple
  USE stokes2d_method1_mg
  USE stokes2d_method2_sv
  USE stokes2d_method2_gv
  USE stokes3d_method0_simple
  USE stokes3d_method1_mg
  
  USE navst3d_method1_mg
  
  IMPLICIT NONE
  
  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  CALL system_init()
  CALL output_init()

  ! The very second thing in every program: 
  ! Initialise the storage management: 
  !
  ! 2.) Initialise FEAT 2.0 storage management:
  CALL storage_init(999, 100)

  ! Call the problem to solve. 2d stokes 0:
  CALL output_lbrk()
  CALL output_line('Calculating 2D Stokes-Problem 0 - simple')
  CALL output_line('----------------------------------------')
  CALL stokes2d_0_simple

  ! Call the problem to solve. 2d stokes 1:
  CALL output_lbrk()
  CALL output_line('Calculating 2D Stokes-Problem 1 - multigrid')
  CALL output_line('-------------------------------------------')
  CALL stokes2d_1_mg

  ! Call the problem to solve. 2d stokes 2:
  CALL output_lbrk()
  CALL output_line('Calculating 2D Stokes-Problem 2 - simple VANCA')
  CALL output_line('----------------------------------------------')
  CALL stokes2d_2_sv

  ! Call the problem to solve. 2d stokes 2:
  CALL output_lbrk()
  CALL output_line('Calculating 2D Stokes-Problem 2 - general VANCA')
  CALL output_line('-----------------------------------------------')
  CALL stokes2d_2_gv

  ! Call the problem to solve. 3d stokes 0:
  CALL output_lbrk()
  CALL output_line('Calculating 3D Stokes-Problem 0 - simple')
  CALL output_line('----------------------------------------')
  CALL stokes3d_0_simple

  ! Call the problem to solve. 3d stokes 1:
  CALL output_lbrk()
  CALL output_line('Calculating 3D Stokes-Problem 1 - multigrid')
  CALL output_line('-------------------------------------------')
  CALL stokes3d_1_mg

  ! Call the problem to solve. 3d navier-stokes 1:
  !CALL output_lbrk()
  !CALL output_line('Calculating 3D Navier-Stokes-Problem 1 - multigrid')
  !CALL output_line('--------------------------------------------------')
  !CALL navst3d_1_mg

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  CALL output_lbrk()
  CALL storage_info(.TRUE.)
  
  ! Clean up the storage management, finish
  CALL storage_done()
  
END PROGRAM
