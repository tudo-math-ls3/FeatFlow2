!##############################################################################
!# ****************************************************************************
!# <name> stokes </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the equation
!#
!#              - Laplace(u) = f
!#
!# on a 2D domain for a scalar function u.
!#
!# There are five examples provided how to solve this problem:
!#
!# The first example (module stokes_method1) discretises and solves this 
!# equation in a direct way, just listing all commands necessary for 
!# initialisation, discretisation, solving and cleanup.
!#
!# The second example (module stokes_method2) separates the different stages
!# of the solution process into different subroutines, which communicate
!# via a collection-structure.
!#
!# The third example (module stokes_method3) separates the different stages
!# of the solution process into different subroutines like example 2.
!# The communication is done using a problem-related structure. For the
!# communication with callback routines during the assembly, a
!# collection structure is set up.
!#
!# The fourth example (module stokes_method4) demonstrates the use of a
!# BiCGStab-solver with ILU(0) preconditioner.
!#
!# The fith example (module stokes_method5) demonstrates the use of a
!# Multigrid solver with with ILU(0) smoother and BiCGStab-coarse grid solver.
!# </purpose>
!##############################################################################

PROGRAM stokes

  USE stokes_method1
  
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
  CALL ZINIT(NNWORK,'feat.msg','log/feat1.err','log/feat1.prt',&
             'log/feat1.sys','log/feat1.trc') 
  
  ! Call the problem to solve. stokes 1:
  PRINT *
  PRINT *,'Calculating Stokes-Problem with method 1'
  PRINT *,'-----------------------------------------'
  CALL stokes1

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  PRINT *
  CALL storage_info(.TRUE.)
  
  ! Clean up the storage management, finish
  CALL storage_done()
  
END PROGRAM
