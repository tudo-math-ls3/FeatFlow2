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
!# There are a couple of examples provided how to solve this problem:
!#
!# The first example (module poisson_method1) discretises and solves this 
!# equation in a direct way, just listing all commands necessary for 
!# initialisation, discretisation, solving and cleanup.
!#
!# The second example (module poisson_method2) separates the different stages
!# of the solution process into different subroutines, which communicate
!# via a collection-structure.
!#
!# The third example (module poisson_method3) separates the different stages
!# of the solution process into different subroutines like example 2.
!# The communication is done using a problem-related structure. For the
!# communication with callback routines during the assembly, a
!# collection structure is set up.
!#
!# The fourth example (module poisson_method4) demonstrates the use of a
!# BiCGStab-solver with ILU(0) preconditioner.
!#
!# The fifth example (module poisson_method5) demonstrates the use of a
!# Multigrid solver with with ILU(0) smoother and UMFPACK-coarse grid solver.
!#
!# The 6th example (module poisson_method6) demonstrates like poisson_method1
!# a simple way of solving the Laplace-equation, but this time introduces
!# an additional fictitious boundary object inside of the domain.
!# </purpose>
!##############################################################################

PROGRAM poisson

  USE poisson_method1
  USE poisson_method2
  USE poisson_method3
  USE poisson_method4
  USE poisson_method5
  USE poisson_method6
  USE poisson_method7
  
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

  ! The very second thing in every program: 
  ! Initialise the storage management: 
  !
  ! 2.) Initialise FEAT 2.0 storage management:
  CALL storage_init(999, 100)

  ! 3.) Initialise old FEAT 1.x storage management for compatibility.
  CALL ZINIT(NNWORK,'feat.msg','log/feat1.err','log/feat1.prt',&
             'log/feat1.sys','log/feat1.trc') 
  
  ! Call the problem to solve. Poisson 1:
  CALL output_lbrk ()
  CALL output_line ('Calculating Laplace-Problem with method 1')
  CALL output_line ('-----------------------------------------')
  CALL poisson1
  
  ! Call the problem to solve. Poisson 2:
  CALL output_lbrk ()
  CALL output_line ('Calculating Laplace-Problem with method 2')
  CALL output_line ('-----------------------------------------')
  CALL poisson2
  
  ! Call the problem to solve. Poisson 3:
  CALL output_lbrk ()
  CALL output_line ('Calculating Laplace-Problem with method 3')
  CALL output_line ('-----------------------------------------')
  CALL poisson3
  
  ! Call the problem to solve. Poisson 4:
  CALL output_lbrk ()
  CALL output_line ('Calculating Laplace-Problem with method 4')
  CALL output_line ('-----------------------------------------')
  CALL poisson4
  
  ! Call the problem to solve. Poisson 5:
  CALL output_lbrk ()
  CALL output_line ('Calculating Laplace-Problem with method 5')
  CALL output_line ('-----------------------------------------')
  CALL poisson5

  ! Call the problem to solve. Poisson 6:
  CALL output_lbrk ()
  CALL output_line ('Calculating Laplace-Problem with method 6')
  CALL output_line ('-----------------------------------------')
  CALL poisson6

  ! Call the problem to solve. Poisson 7:
  CALL output_lbrk ()
  CALL output_line ('Calculating Laplace-Problem with method 7')
  CALL output_line ('-----------------------------------------')
  CALL poisson7

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  CALL output_lbrk ()
  CALL storage_info(.TRUE.)
  
  ! Clean up the storage management, finish
  CALL storage_done()
  
END PROGRAM
