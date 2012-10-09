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
  use poisson2d_method0_neumann
  use poisson2d_method0_block
  use poisson2d_method0_smart
  use poisson2d_method0_cmsort
  use poisson2d_method1_mg
  use poisson2d_method0_agmg
  use poisson2d_method1_em30
  use poisson2d_method1_robin
  use poisson2d_method1_fbc
  use poisson2d_method1_hadapt
  use poisson2d_method1_l2prj
  use poisson2d_method1_prolmat
  use poisson2d_method2
  use poisson2d_method2_collect
  use poisson2d_method2_cmsort
  use poisson2d_method2_mg
  use poisson3d_method0_simple
  use poisson3d_method0_agmg
  use poisson3d_method1_mg
  use poisson3d_method1_em30
  use poisson2d_method1_ncc

  implicit none

  ! local variables
  character(len=SYS_STRLEN) :: slogdir,slogfile

  ! The very first thing in every application:
  ! Initialise system-wide settings:

  call system_init()

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
  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Call the problem to solve. Poisson 1D method 1 - simple:
  call output_lbrk ()
  call output_line ("Calculating Poisson-1D-Problem with method 0 - simple")
  call output_line ("-----------------------------------------------------")
  call poisson1d_0_simple

  ! Call the problem to solve. Poisson 1D method 1 - multigrid:
  call output_lbrk ()
  call output_line ("Calculating Poisson-1D-Problem with method 1 - multigrid")
  call output_line ("--------------------------------------------------------")
  call poisson1d_1_mg

  ! Call the problem to solve. Poisson 2D method 1 - simple:
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 0 - simple")
  call output_line ("-----------------------------------------------------")
  call poisson2d_0_simple

  ! Call the problem to solve. Poisson 2D method 1 - smart:
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 0 - smart")
  call output_line ("----------------------------------------------------")
  call poisson2d_0_smart

  ! Call the problem to solve. Poisson 2D method 1 - pure Neumann problem:
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 0 - Neumann")
  call output_line ("------------------------------------------------------")
  call poisson2d_0_neumann

  ! Call the problem to solve. Poisson 2D method 1 - CM-sorting:
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 0 - CM-sorting")
  call output_line ("---------------------------------------------------------")
  call poisson2d_0_cmsort

  ! Call the problem to solve. Poisson 2D method 0 - block variant:
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 0 - block")
  call output_line ("-----------------------------------------------------")
  call poisson2d_0_block

#ifdef ENABLE_AGMG
  ! Call the problem to solve. Poisson 2D method 1 - algebraic multigrid:
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 0 - AGMG")
  call output_line ("--------------------------------------------------------")
  call poisson2d_0_agmg
#endif

  ! Call the problem to solve. Poisson 2D method 1 - nonconstant coefficients:
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 1 - ncc")
  call output_line ("--------------------------------------------------")
  call poisson2d_1_ncc

  ! Call the problem to solve. Poisson 2D method 1 - multigrid:
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 1 - multigrid")
  call output_line ("--------------------------------------------------------")
  call poisson2d_1_mg

  ! Call the problem to solve. Poisson 1: Support for nonconforming elements
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 1 - EM30")
  call output_line ("---------------------------------------------------")
  call poisson2d_1_em30

  ! Call the problem to solve. Poisson 1: Support for nonconforming elements
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 1 - Robin BC")
  call output_line ("-------------------------------------------------------")
  call poisson2d_1_robin

  ! Call the problem to solve. Poisson 1: Fictitious boundary support
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 1 - FBC")
  call output_line ("--------------------------------------------------")
  call poisson2d_1_fbc

  ! Call the problem to solve. Poisson 1: h-adaptivity
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 1 - hadapt")
  call output_line ("-----------------------------------------------------")
  call poisson2d_1_hadapt

  ! Call the problem to solve. Poisson 2D method 1 - L2-projection:
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 1 - L2-projection")
  call output_line ("------------------------------------------------------------")
  call poisson2d_1_l2prj

  ! Call the problem to solve. Poisson 2D method 1 - Prolongation matrix:
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 1 - Prol.-Matrix")
  call output_line ("-----------------------------------------------------------")
  call poisson2d_1_prolmat

  ! Call the problem to solve. Poisson 2D method 2:
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 2")
  call output_line ("--------------------------------------------")
  call poisson2d_2

  ! Call the problem to solve. Poisson 3: Sorting with Cuthill McKee
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 2 - CM-sorting")
  call output_line ("---------------------------------------------------------")
  call poisson2d_2_cmsort

  ! Call the problem to solve. Poisson 5:
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 2 - multigrid")
  call output_line ("--------------------------------------------------------")
  call poisson2d_2_mg

  ! Call the problem to solve. Poisson 3: Collection support
  call output_lbrk ()
  call output_line ("Calculating Poisson-2D-Problem with method 2 - collection")
  call output_line ("---------------------------------------------------------")
  call poisson2d_2_collect

  ! Call the problem to solve. Poisson3D-1:
  call output_lbrk ()
  call output_line ("Calculating Poisson-3D-Problem with method 0 - simple")
  call output_line ("-----------------------------------------------------")
  call poisson3d_0_simple

#ifdef ENABLE_AGMG
  ! Call the problem to solve. Poisson3D-1:
  call output_lbrk ()
  call output_line ("Calculating Poisson-3D-Problem with method 0 - AGMG")
  call output_line ("-----------------------------------------------------")
  call poisson3d_0_agmg
#endif

  ! Call the problem to solve. Poisson3D-1:
  call output_lbrk ()
  call output_line ("Calculating Poisson-3D-Problem with method 1 - multigrid")
  call output_line ("--------------------------------------------------------")
  call poisson3d_1_mg

  ! Call the problem to solve. Poisson3D-7:
  call output_lbrk ()
  call output_line ("Calculating Poisson-3D-Problem with method 1 - EM30")
  call output_line ("---------------------------------------------------")
  call poisson3d_1_em30

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display "Handles in use=0" and "Memory in use=0"!
  call output_lbrk ()
  call storage_info(.true.)

  ! Clean up the storage management, finish
  call storage_done()

end program
