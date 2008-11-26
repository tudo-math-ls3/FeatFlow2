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
   
  use poisson1d_method0_simple
  use poisson1d_method1_mg
  use poisson2d_method0_simple
  use poisson2d_method1_mg
  use poisson2d_method1_em30
  use poisson2d_method1_fbc
  use poisson2d_method1_hadapt
  use poisson2d_method1_l2prj
  use poisson2d_method1_prolmat
  use poisson2d_method2
  use poisson2d_method2_collect
  use poisson2d_method2_cmsort
  use poisson2d_method2_mg
  use poisson3d_method0_simple
  use poisson3d_method1_mg
  use poisson3d_method1_em30
  use poisson2d_method1_ncc
  
  implicit none
  
  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  call system_init()
  
  ! Initialise the output system. Write the program output to screen as
  ! well as to the file 'log/output.txt'.
  call output_init ('./log/output.txt')

  ! The very second thing in every program: 
  ! Initialise the FEAT 2.0 storage management: 
  call storage_init(999, 100)
 
  ! Call the problem to solve. Poisson 1D method 1 - simple:
  call output_lbrk ()
  call output_line ('Calculating Poisson-1D-Problem with method 0 - simple')
  call output_line ('-----------------------------------------------------')
  call poisson1d_0_simple

  ! Call the problem to solve. Poisson 1D method 1 - multigrid:
  call output_lbrk ()
  call output_line ('Calculating Poisson-1D-Problem with method 1 - multigrid')
  call output_line ('--------------------------------------------------------')
  call poisson1d_1_mg

  ! Call the problem to solve. Poisson 2D method 1 - simple:
  call output_lbrk ()
  call output_line ('Calculating Poisson-2D-Problem with method 0 - simple')
  call output_line ('-----------------------------------------------------')
  call poisson2d_0_simple

  ! Call the problem to solve. Poisson 2D method 1 - nonconstant coefficients:
  call output_lbrk ()
  call output_line ('Calculating Poisson-2D-Problem with method 1 - ncc')
  call output_line ('--------------------------------------------------')
  call poisson2d_1_ncc
  
  ! Call the problem to solve. Poisson 2D method 1 - multigrid:
  call output_lbrk ()
  call output_line ('Calculating Poisson-2D-Problem with method 1 - multigrid')
  call output_line ('--------------------------------------------------------')
  call poisson2d_1_mg

  ! Call the problem to solve. Poisson 1: Support for nonconforming elements
  call output_lbrk ()
  call output_line ('Calculating Poisson-2D-Problem with method 1 - EM30')
  call output_line ('---------------------------------------------------')
  call poisson2d_1_em30

  ! Call the problem to solve. Poisson 1: Fictitious boundary support
  call output_lbrk ()
  call output_line ('Calculating Poisson-2D-Problem with method 1 - FBC')
  call output_line ('--------------------------------------------------')
  call poisson2d_1_fbc

  ! Call the problem to solve. Poisson 1: h-adaptivity
  call output_lbrk ()
  call output_line ('Calculating Poisson-2D-Problem with method 1 - hadapt')
  call output_line ('-----------------------------------------------------')
  call poisson2d_1_hadapt
  
  ! Call the problem to solve. Poisson 2D method 1 - L2-projection:
  call output_lbrk ()
  call output_line ('Calculating Poisson-2D-Problem with method 1 - L2-projection')
  call output_line ('------------------------------------------------------------')
  call poisson2d_1_l2prj
  
  ! Call the problem to solve. Poisson 2D method 1 - Prolongation matrix:
  call output_lbrk ()
  call output_line ('Calculating Poisson-2D-Problem with method 1 - Prol.-Matrix')
  call output_line ('-----------------------------------------------------------')
  call poisson2d_1_prolmat

  ! Call the problem to solve. Poisson 2D method 2:
  call output_lbrk ()
  call output_line ('Calculating Poisson-2D-Problem with method 2')
  call output_line ('--------------------------------------------')
  call poisson2d_2
  
  ! Call the problem to solve. Poisson 3: Sorting with Cuthill McKee
  call output_lbrk ()
  call output_line ('Calculating Poisson-2D-Problem with method 2 - CM-sorting')
  call output_line ('---------------------------------------------------------')
  call poisson2d_2_cmsort
  
  ! Call the problem to solve. Poisson 5:
  call output_lbrk ()
  call output_line ('Calculating Poisson-2D-Problem with method 2 - multigrid')
  call output_line ('--------------------------------------------------------')
  call poisson2d_2_mg

  ! Call the problem to solve. Poisson 3: Collection support
  call output_lbrk ()
  call output_line ('Calculating Poisson-2D-Problem with method 2 - collection')
  call output_line ('---------------------------------------------------------')
  call poisson2d_2_collect
  
  ! Call the problem to solve. Poisson3D-1:
  call output_lbrk ()
  call output_line ('Calculating Poisson-3D-Problem with method 0 - simple')
  call output_line ('-----------------------------------------------------')
  call poisson3d_0_simple

  ! Call the problem to solve. Poisson3D-1:
  call output_lbrk ()
  call output_line ('Calculating Poisson-3D-Problem with method 1 - multigrid')
  call output_line ('--------------------------------------------------------')
  call poisson3d_1_mg

  ! Call the problem to solve. Poisson3D-7:
  call output_lbrk ()
  call output_line ('Calculating Poisson-3D-Problem with method 1 - EM30')
  call output_line ('---------------------------------------------------')
  call poisson3d_1_em30

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk ()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
