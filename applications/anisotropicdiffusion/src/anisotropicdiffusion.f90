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
  
  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  call system_init(); sys_haltmode = SYS_HALT_THROWFPE
  
  ! Initialise the output system. Write the program output to screen as
  ! well as to the file 'log/output.txt'.
  call output_init ('./log/output.txt')

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
  call storage_done()
  
end program
