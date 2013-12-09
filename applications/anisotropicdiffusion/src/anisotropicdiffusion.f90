!##############################################################################
!# ****************************************************************************
!# <name> anisotropicdiffusion </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the equation
!#
!#              - div (A grad(u) = f
!#
!# on a 2D domain for a scalar function u and a matrix
!#
!#   A = ( a11 a12 )
!#       ( a21 a22 )
!#
!# a11,...,a22 are configured with a .dat file. The matrix may be
!# rotated according to
!#
!#   A~ =  ( cos(t) -sin(t) )  ( a11 a12 )  (  cos(t) sin(t) )
!#         ( sin(t)  cos(t) )  ( a21 a22 )  ( -sin(t) cos(t) )
!#
!# The example (module anisotropicdiffusion_method1) discretises and
!# solves this equation in a direct way, just listing all commands necessary
!# for initialisation, discretisation, solving and cleanup.
!##############################################################################

program anisotropicdiffusion

  use anisotropicdiffusion_method1
  use anisotropicdiffusion_method2
  
  implicit none
  
  ! local variables
  character(len=SYS_STRLEN) :: slogdir,slogfile

  ! The very first thing in every application:
  ! Initialise system-wide settings:
  
  call sys_init()
  
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

  ! Call the problem to solve.
  call output_lbrk ()
  call output_line ('Calculating anisotropic diffusion problem with method 1')
  call output_line ('-------------------------------------------------------')
  call anisotropicdiffusion1

  ! Call the problem to solve.
  call output_lbrk ()
  call output_line ('Calculating anisotropic diffusion problem with method 2')
  call output_line ('-------------------------------------------------------')
  call anisotropicdiffusion2
  
  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk ()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call output_done()
  call storage_done()
  
end program
