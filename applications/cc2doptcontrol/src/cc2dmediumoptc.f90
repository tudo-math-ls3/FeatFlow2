!##############################################################################
!# ****************************************************************************
!# <name> cc2dmediumoptc </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising a 
!# stationary and nonstationary Navier-Stokes optimal control problem 
!#
!#  $$ min J(y,u) = 1/2||y-z||_{L^2} + \gamma/2||y(T)-z(T)||_{L^2} + \alphga/2||u||^2 $$
!#
!#  $$- \nu Delta(y) + y*\Nabla(y) + \Nabla p = f $$
!#  $$ \Nabla \cdot y = 0$$
!#  $$- \nu Delta(\lambda) - y*\Nabla(\lambda) + \lambda\Nabla y + \Nabla \xi = y-z $$
!#  $$ \Nabla \cdot \lambda = 0$$
!#              
!#
!# on a 2D domain for a 2D function $y=(y_1,y_2)$, a pressure $p$,
!# a dual velocity $\lambda$ and a dual pressure $\xi$. $u$ is the control
!# and $z$ a desired flow field.
!#
!# on a 2D domain for a 2D velocity function $y=(y_1,y_2)$, a pressure $p$,
!# a 2D dual velocity $\lambda=(\lambda_1,\lambda_2)$ and a dual pressure
!# $\xi$.
!#
!# A linear block system with 6 solution components is set up,
!# discretised and solved.
!#
!# </purpose>
!##############################################################################

PROGRAM cc2dmediumoptc

  USE cc2dmedium_method2
  
  IMPLICIT NONE
  
  INCLUDE 'cmem.inc'
  INCLUDE 'cout.inc'
  INCLUDE 'cerr.inc'
  INCLUDE 'cfileout.inc'

  ! As we still use some FEAT 1.x routines, we have to initialise some
  ! output variables.
  M = 0
  MT = 2
  MTERM = 6
  MFILE = 0
  ICHECK = 0

  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  CALL system_init()
  
  ! Initialise log file for output.
  CALL output_init ('log/output.txt')
  
  ! The very second thing in every program: 
  ! Initialise the storage management: 
  !
  ! 2.) Initialise FEAT 2.0 storage management:
  CALL storage_init(999, 100)
  
  ! 3.) Initialise old FEAT 1.x storage management for compatibility.
  CALL ZINIT(NNWORK,'feat.msg','log/feat1.err','log/feat1.prt',&
             'log/feat1.sys','log/feat1.trc') 
  
  ! Call the problem to solve. 
  CALL output_lbrk ()
  CALL output_line ('Calculating cc2dmedium-Problem')
  CALL output_separator (OU_SEP_MINUS)
  
  CALL cc2dmedium2optc ()

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  CALL output_lbrk ()
  CALL storage_info(.TRUE.)
  
  ! Clean up the storage management, finish
  CALL storage_done()
  
END PROGRAM
