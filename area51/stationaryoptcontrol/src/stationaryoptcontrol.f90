!##############################################################################
!# ****************************************************************************
!# <name> poisson </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the equation
!#
!#              - Laplace(u) = f
!#
!# for a scalar function u.
!#
!# There are a couple of examples provided how to solve this problem.
!# Each example has its own characteristics how to solve the problem.
!#
!# There are examples for 1D, 2D and 3D triangulations. The
!# poissonXd_method0_simple modules show the very basic steps how to solve the
!# Poisson problem. On top of that, the poissonXd_method1_XXXX modules
!# extend the poissonXd_method0_simple to more general situations (mmultigrid,
!# fictitious boundary,...), still in the style of poissonXd_method0_simple.
!# The poissonXd_method2_XXXX routines then decompose the different
!# tasks of the problem into different subroutines, thus bringing more
!# structure into the code and highlighting more and more features of the
!# kernel.
!#
!# </purpose>
!##############################################################################

program poisson
   
  use poisson2d_method0_simple
  use stokes2d_method0_simple
  use navstokes2d_method0_simple
  
  implicit none
  
  ! The very first thing in every application:
  ! Initialise system-wide settings:
  
  call sys_init()
  
  ! Initialise the output system. Write the program output to screen as
  ! well as to the file 'log/output.txt'.
  call output_init ('./log/output.txt')

  ! The very second thing in every program:
  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)
 
  ! Call the problem to solve. Poisson 2D method 1 - simple:
  call output_lbrk ()
  call output_line ('Calculating Poisson-2D-Problem with method 0 - simple')
  call output_line ('-----------------------------------------------------')
  !call poisson2d_0_simple

  ! Call the problem to solve. Stokes 2D method 1 - simple:
  call output_lbrk ()
  call output_line ('Calculating Stokes-2D-Problem with method 0 - simple')
  call output_line ('-----------------------------------------------------')
  call stokes2d_0_simple

  ! Call the problem to solve. Navier-Stokes 2D method 1 - simple:
  call output_lbrk ()
  call output_line ('Calculating Navier-Stokes-2D-Problem with method 0 - simple')
  call output_line ('-----------------------------------------------------------')
  !call navstokes2d_0_simple

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk ()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
