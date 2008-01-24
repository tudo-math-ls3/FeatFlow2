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
!# on a 2D domain for a scalar function u.
!#
!# There are a couple of examples provided how to solve this problem:
!#
!# The first example (module poisson_method1) discretises and solves this 
!# equation in a direct way, just listing all commands necessary for 
!# initialisation, discretisation, solving and cleanup.
!#
!# The second example (module poisson_method2) separates the different stages
!# of the solution process into different subroutines, which communicate
!# via a collection-structure.
!#
!# The third example (module poisson_method3) separates the different stages
!# of the solution process into different subroutines like example 2.
!# The communication is done using a problem-related structure. For the
!# communication with callback routines during the assembly, a
!# collection structure is set up.
!#
!# The fourth example (module poisson_method4) demonstrates the use of a
!# BiCGStab-solver with ILU(0) preconditioner.
!#
!# The fifth example (module poisson_method5) demonstrates the use of a
!# Multigrid solver with with ILU(0) smoother and UMFPACK-coarse grid solver.
!#
!# The 6th example (module poisson_method6) demonstrates like poisson_method1
!# a simple way of solving the Poisson-equation, but this time introduces
!# an additional fictitious boundary object inside of the domain.
!#
!# The 7th example is based again on example 1 and shows how to transform
!# a solution vector such that it gan be read by GMV. For this purpose,
!# an appropriate interpolation routine is called, that transforms an
!# arbitrary solution vector into a solution in the $Q_1$ space.
!# </purpose>
!##############################################################################

PROGRAM poisson
   
  USE poisson1d_method0_simple
  USE poisson1d_method1_mg
  USE poisson2d_method0_simple
  USE poisson2d_method1_mg
  USE poisson2d_method1_em30
  USE poisson2d_method1_fbc
  USE poisson2d_method1_hadapt
  USE poisson2d_method2
  USE poisson2d_method3
  USE poisson2d_method3_cmsort
  USE poisson2d_method3_mg
  USE poisson3d_method0_simple
  USE poisson3d_method1_mg
  USE poisson3d_method1_em30
  
  IMPLICIT NONE
  
  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  CALL system_init()
  
  ! Initialise the output system. Write the program output to screen as
  ! well as to the file 'log/output.txt'.
  CALL output_init ('./log/output.txt')

  ! The very second thing in every program: 
  ! Initialise the FEAT 2.0 storage management: 
  CALL storage_init(999, 100)
  
  ! Call the problem to solve. Poisson 1D method 1 - simple:
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-1D-Problem with method 0 - simple')
  CALL output_line ('-----------------------------------------------------')
  CALL poisson1d_0_simple

  ! Call the problem to solve. Poisson 1D method 1 - multigrid:
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-1D-Problem with method 1 - multigrid')
  CALL output_line ('--------------------------------------------------------')
  CALL poisson1d_1_mg

  ! Call the problem to solve. Poisson 2D method 1 - simple:
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-2D-Problem with method 0 - simple')
  CALL output_line ('-----------------------------------------------------')
  CALL poisson2d_0_simple
  
  ! Call the problem to solve. Poisson 2D method 1 - multigrid:
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-2D-Problem with method 1-mg')
  CALL output_line ('-----------------------------------------------')
  CALL poisson2d_1_mg

  ! Call the problem to solve. Poisson 1: Support for nonconforming elements
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-2D-Problem with method 1 - EM30')
  CALL output_line ('---------------------------------------------------')
  CALL poisson2d_1_em30

  ! Call the problem to solve. Poisson 1: Fictitious boundary support
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-2D-Problem with method 1 - FBC')
  CALL output_line ('--------------------------------------------------')
  CALL poisson2d_1_fbc

  ! Call the problem to solve. Poisson 1: h-adaptivity
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-2D-Problem with method 1 - hadapt')
  CALL output_line ('-----------------------------------------------------')
  CALL poisson2d_1_hadapt

  ! Call the problem to solve. Poisson 2D method 2:
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-2D-Problem with method 2')
  CALL output_line ('--------------------------------------------')
  CALL poisson2d_2
  
  ! Call the problem to solve. Poisson 3:
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-2D-Problem with method 3')
  CALL output_line ('--------------------------------------------')
  CALL poisson2d_3
  
  ! Call the problem to solve. Poisson 3: Sorting with Cuthill McKee
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-2D-Problem with method 3 - CM-sorting')
  CALL output_line ('---------------------------------------------------------')
  CALL poisson2d_3_cmsort
  
  ! Call the problem to solve. Poisson 5:
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-2D-Problem with method 3 - Multigrid')
  CALL output_line ('--------------------------------------------------------')
  CALL poisson2d_3_mg

  ! Call the problem to solve. Poisson3D-1:
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-3D-Problem with method 0 - simple')
  CALL output_line ('-----------------------------------------------------')
  CALL poisson3d_0_simple

  ! Call the problem to solve. Poisson3D-1:
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-3D-Problem with method 1 - multigrid')
  CALL output_line ('--------------------------------------------------------')
  CALL poisson3d_1_mg

  ! Call the problem to solve. Poisson3D-7:
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-3D-Problem with method 1 - EM30')
  CALL output_line ('---------------------------------------------------')
  CALL poisson3d_1_em30

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  CALL output_lbrk ()
  CALL storage_info(.TRUE.)
  
  ! Clean up the storage management, finish
  CALL storage_done()
  
END PROGRAM
