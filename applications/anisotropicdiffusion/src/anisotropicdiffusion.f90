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

PROGRAM anisotropicdiffusion

  USE anisotropicdiffusion_method1
  USE anisotropicdiffusion_method2
  
  IMPLICIT NONE
  
  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  CALL system_init(); sys_haltmode = SYS_HALT_THROWFPE
  
  ! Initialise the output system. Write the program output to screen as
  ! well as to the file 'log/output.txt'.
  CALL output_init ('./log/output.txt')

  ! The very second thing in every program: 
  ! Initialise the FEAT 2.0 storage management: 
  CALL storage_init(999, 100)

  ! Call the problem to solve. 
  CALL output_lbrk ()
  CALL output_line ('Calculating anisotropic diffusion problem with method 1')
  CALL output_line ('-------------------------------------------------------')
  CALL anisotropicdiffusion1

  ! Call the problem to solve. 
  CALL output_lbrk ()
  CALL output_line ('Calculating anisotropic diffusion problem with method 2')
  CALL output_line ('-------------------------------------------------------')
  CALL anisotropicdiffusion2
  
  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  CALL output_lbrk ()
  CALL storage_info(.TRUE.)
  
  ! Clean up the storage management, finish
  CALL storage_done()
  
END PROGRAM
