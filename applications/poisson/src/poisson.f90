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
!# There are two examples provided how to solve this problem:
!#
!# The first example (module poissonmeth1) discretises and solves this 
!# equation in a direct way, just listing all commands necessary for 
!# initialisation, discretisation, solving and cleanup.
!#
!# The second example (module poissonmeth2) separates the different stages
!# of the solution process into different subroutines, which communicate
!# via a collection-structure.
!# </purpose>
!##############################################################################

PROGRAM poisson

  USE poissonmeth1
  USE poissonmeth2
  
  IMPLICIT NONE
  
  INCLUDE 'cmem.inc'
  INCLUDE 'cout.inc'
  INCLUDE 'cerr.inc'

  ! As we still use some FEAT 1.x routines, we have to initialise some
  ! output variables.

  M = 0
  ICHECK = 0

  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  CALL system_init()

  ! The very second thing in every program: 
  ! Initialise the storage management: 
  !
  ! 2.) Initialise FEAT 2.0 storage management:
  CALL storage_init(999, 100)

  ! 3.) Initialise old FEAT 1.x storage management for compatibility.
  CALL ZINIT(NNWORK,'feat.msg','data/cc2d.err','data/cc2d.prt',&
             'data/cc2d.sys','data/cc2d.trc') 
  
  ! Call the problem to solve. Poisson 1:
  PRINT *,'Calculating Laplace-Problem with method 1'
  PRINT *,'-----------------------------------------'
  CALL poisson_method1
  
  PRINT *

  ! Call the problem to solve. Poisson 2:
  PRINT *,'Calculating Laplace-Problem with method 2'
  PRINT *,'-----------------------------------------'
  CALL poisson_method2
  
  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  PRINT *
  CALL storage_info()
  
END PROGRAM
