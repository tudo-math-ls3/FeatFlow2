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

  use fsystem
  use genoutput
  use paramlist

  use poisson1d_method0_simple
  use poisson1d_method1_mg
  use poisson2d_method0_agmg
  use poisson2d_method0_block
  use poisson2d_method0_cmsort
  use poisson2d_method0_neumann
  use poisson2d_method0_ravthomas
  use poisson2d_method0_simple
  use poisson2d_method0_smart
  use poisson2d_method1_em30
  use poisson2d_method1_fbc
  use poisson2d_method1_hadapt
  use poisson2d_method1_l2prj
  use poisson2d_method1_mg
  use poisson2d_method1_ncc
  use poisson2d_method1_prolmat
  use poisson2d_method1_robin
  use poisson2d_method2
  use poisson2d_method2_cmsort
  use poisson2d_method2_collect
  use poisson2d_method2_mg
  use poisson3d_method0_agmg
  use poisson3d_method0_simple
  use poisson3d_method1_em30
  use poisson3d_method1_mg

  implicit none

  ! local parameter list
  type(t_parlist) :: rparamlist

  ! local variables
  character(len=SYS_STRLEN) :: slogdir,slogfile,smaster
  integer :: i

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
  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Initialise the parameter list
  call parlst_init(rparamlist)
  
  ! Read the master file
  call sys_getcommandLineArg(1,smaster,sdefault='./data/master.dat')
  call parlst_readfromfile(rparamlist, smaster)

  ! Call the problem to solve. Poisson 1D method 1 - simple:
  call parlst_getvalue_int(rparamlist, "", "POISSON1D_METHOD0_SIMPLE", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-1D-Problem with method 0 - simple")
    call output_line ("-----------------------------------------------------")
    call poisson1d_0_simple
  end if

  ! Call the problem to solve. Poisson 1D method 1 - multigrid:
  call parlst_getvalue_int(rparamlist, "", "POISSON1D_METHOD1_MG", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-1D-Problem with method 1 - multigrid")
    call output_line ("--------------------------------------------------------")
    call poisson1d_1_mg
  end if

  ! Call the problem to solve. Poisson 2D method 1 - simple:
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD0_SIMPLE", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 0 - simple")
    call output_line ("-----------------------------------------------------")
    call poisson2d_0_simple
  end if

  ! Call the problem to solve. Poisson 2D method 1 - smart:
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD0_SMART", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 0 - smart")
    call output_line ("----------------------------------------------------")
    call poisson2d_0_smart
  end if

  ! Call the problem to solve. Poisson 2D method 1 - pure Neumann problem:
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD0_NEUMANN", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 0 - Neumann")
    call output_line ("------------------------------------------------------")
    call poisson2d_0_neumann
  end if

  ! Call the problem to solve. Poisson 2D method 1 - CM-sorting:
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD0_CMSORT", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 0 - CM-sorting")
    call output_line ("---------------------------------------------------------")
    call poisson2d_0_cmsort
  end if

  ! Call the problem to solve. Poisson 2D method 0 - block variant:
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD0_BLOCK", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 0 - block")
    call output_line ("-----------------------------------------------------")
    call poisson2d_0_block
  end if

  ! Call the problem to solve. Poisson 2D method 0 - Raviart-Thomas variant:
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD0_RAVTHOMAS", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 0 - RavThomas")
    call output_line ("--------------------------------------------------------")
    call poisson2d_0_ravthomas
  end if

  ! Call the problem to solve. Poisson 2D method 1 - algebraic multigrid:
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD0_AGMG", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 0 - AGMG")
    call output_line ("---------------------------------------------------")
    call poisson2d_0_agmg
  end if

  ! Call the problem to solve. Poisson 2D method 1 - nonconstant coefficients:
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD1_NCC", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 1 - ncc")
    call output_line ("--------------------------------------------------")
    call poisson2d_1_ncc
  end if

  ! Call the problem to solve. Poisson 2D method 1 - multigrid:
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD1_MG", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 1 - multigrid")
    call output_line ("--------------------------------------------------------")
    call poisson2d_1_mg
  end if

  ! Call the problem to solve. Poisson 1: Support for nonconforming elements
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD1_EM30", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 1 - EM30")
    call output_line ("---------------------------------------------------")
    call poisson2d_1_em30
  end if

  ! Call the problem to solve. Poisson 1: Support for nonconforming elements
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD1_ROBIN", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 1 - Robin BC")
    call output_line ("-------------------------------------------------------")
    call poisson2d_1_robin
  end if

  ! Call the problem to solve. Poisson 1: Fictitious boundary support
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD1_FBC", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 1 - FBC")
    call output_line ("--------------------------------------------------")
    call poisson2d_1_fbc
  end if

  ! Call the problem to solve. Poisson 1: h-adaptivity
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD1_HADAPT", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 1 - hadapt")
    call output_line ("-----------------------------------------------------")
    call poisson2d_1_hadapt
  end if

  ! Call the problem to solve. Poisson 2D method 1 - L2-projection:
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD1_L2PRJ", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 1 - L2-projection")
    call output_line ("------------------------------------------------------------")
    call poisson2d_1_l2prj
  end if

  ! Call the problem to solve. Poisson 2D method 1 - Prolongation matrix:
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD1_PROLMAT", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 1 - Prol.-Matrix")
    call output_line ("-----------------------------------------------------------")
    call poisson2d_1_prolmat
  end if

  ! Call the problem to solve. Poisson 2D method 2:
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD2", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 2")
    call output_line ("--------------------------------------------")
    call poisson2d_2
  end if

  ! Call the problem to solve. Poisson 3: Sorting with Cuthill McKee
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD2_CMSORT", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 2 - CM-sorting")
    call output_line ("---------------------------------------------------------")
    call poisson2d_2_cmsort
  end if

  ! Call the problem to solve. Poisson 5:
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD2_MG", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 2 - multigrid")
    call output_line ("--------------------------------------------------------")
    call poisson2d_2_mg
  end if

  ! Call the problem to solve. Poisson 3: Collection support
  call parlst_getvalue_int(rparamlist, "", "POISSON2D_METHOD2_COLLECT", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-2D-Problem with method 2 - collection")
    call output_line ("---------------------------------------------------------")
    call poisson2d_2_collect
  end if

  ! Call the problem to solve. Poisson3D-1:
  call parlst_getvalue_int(rparamlist, "", "POISSON3D_METHOD0_SIMPLE", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-3D-Problem with method 0 - simple")
    call output_line ("-----------------------------------------------------")
    call poisson3d_0_simple
  end if

  ! Call the problem to solve. Poisson3D-1:
  call parlst_getvalue_int(rparamlist, "", "POISSON3D_METHOD0_AGMG", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-3D-Problem with method 0 - AGMG")
    call output_line ("---------------------------------------------------")
    call poisson3d_0_agmg
  end if

  ! Call the problem to solve. Poisson3D-1:
  call parlst_getvalue_int(rparamlist, "", "POISSON3D_METHOD1_MG", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-3D-Problem with method 1 - multigrid")
    call output_line ("--------------------------------------------------------")
    call poisson3d_1_mg
  end if

  ! Call the problem to solve. Poisson3D-7:
  call parlst_getvalue_int(rparamlist, "", "POISSON3D_METHOD1_EM30", i, 0)
  if (i.ne. 0) then
    call output_lbrk ()
    call output_line ("Calculating Poisson-3D-Problem with method 1 - EM30")
    call output_line ("---------------------------------------------------")
    call poisson3d_1_em30
  end if

  ! Release parameter list
  call parlst_done(rparamlist)

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display "Handles in use=0" and "Memory in use=0"!
  call output_lbrk ()
  call storage_info(.true.)

  ! Clean up the storage management, finish
  call storage_done()

end program
