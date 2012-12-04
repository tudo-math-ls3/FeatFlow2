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

!  use poisson1d_method0_simple
!  use poisson1d_method1_mg
!  use poisson2d_method0_simple
  use poisson2d_method1_mg
  use poisson2d_method1_deflgmres
!  use poisson2d_method1_em30
!  use poisson2d_method1_robin
!  use poisson2d_method1_fbc
!  use poisson2d_method1_hadapt
!  use poisson2d_method1_l2prj
!  use poisson2d_method1_prolmat
!  use poisson2d_method2
!  use poisson2d_method2_collect
!  use poisson2d_method2_cmsort
!  use poisson2d_method2_mg
!  use poisson3d_method0_simple
!  use poisson3d_method1_mg
!  use poisson3d_method1_em30
!  use poisson2d_method1_ncc

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

  ! Comment in/out here what to test!
  !call poisson2d_1_deflgmres_quad
  !call poisson2d_1_mg_quad
  call poisson2d_1_deflgmres_bench1
  !call poisson2d_1_mg_bench1

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk ()
  call storage_info(.true.)

  ! Clean up the storage management, finish
  call storage_done()

end program
