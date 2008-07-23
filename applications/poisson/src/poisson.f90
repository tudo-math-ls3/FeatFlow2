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

PROGRAM poisson
   
  USE poisson1d_method0_simple
  USE poisson1d_method1_mg
  USE poisson2d_method0_simple
  USE poisson2d_method1_mg
  USE poisson2d_method1_em30
  USE poisson2d_method1_fbc
  USE poisson2d_method1_hadapt
  USE poisson2d_method1_l2prj
  USE poisson2d_method2
  USE poisson2d_method2_collect
  USE poisson2d_method2_cmsort
  USE poisson2d_method2_mg
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
  CALL output_line ('Calculating Poisson-2D-Problem with method 1 - multigrid')
  CALL output_line ('--------------------------------------------------------')
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
  
  ! Call the problem to solve. Poisson 2D method 1 - L2-projection:
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-2D-Problem with method 1 - L2-projection')
  CALL output_line ('------------------------------------------------------------')
  CALL poisson2d_1_l2prj

  ! Call the problem to solve. Poisson 2D method 2:
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-2D-Problem with method 2')
  CALL output_line ('--------------------------------------------')
  CALL poisson2d_2
  
  ! Call the problem to solve. Poisson 3: Sorting with Cuthill McKee
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-2D-Problem with method 2 - CM-sorting')
  CALL output_line ('---------------------------------------------------------')
  CALL poisson2d_2_cmsort
  
  ! Call the problem to solve. Poisson 5:
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-2D-Problem with method 2 - multigrid')
  CALL output_line ('--------------------------------------------------------')
  CALL poisson2d_2_mg

  ! Call the problem to solve. Poisson 3: Collection support
  CALL output_lbrk ()
  CALL output_line ('Calculating Poisson-2D-Problem with method 2 - collection')
  CALL output_line ('---------------------------------------------------------')
  CALL poisson2d_2_collect
  
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
