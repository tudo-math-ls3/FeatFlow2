!##############################################################################
!# ****************************************************************************
!# <name> burgers1d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the equation
!#
!#    $$  u_t  +  u*u_x  -  \nu u_xx  =  0,   u(x,0) = sin(Pi*x)  $$
!#
!# on a 2D domain $\Omega=[0,1] \times [0,1]$ for a scalar function $u(x,t)$.
!# </purpose>
!##############################################################################

PROGRAM burgers1d

  USE burgers1d_method5
  USE burgers1d_method6
  
  IMPLICIT NONE
  
  INCLUDE 'cmem.inc'
  INCLUDE 'cout.inc'
  INCLUDE 'cerr.inc'
  INCLUDE 'cfileout.inc'

  ! As we still use some FEAT 1.x routines, we have to initialise some
  ! output variables.

  M = 0
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
  
  ! Call the problem to solve. 
  PRINT *
  PRINT *,'Calculating Burgers1D-Problem with method 5'
  PRINT *,'-------------------------------------------'
  CALL burgers1d5

  ! Call the problem to solve. 
  PRINT *
  PRINT *,'Calculating Burgers1D-Problem with method 6'
  PRINT *,'-------------------------------------------'
  CALL burgers1d6

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  PRINT *
  CALL storage_info(.TRUE.)
  
  ! Clean up the storage management, finish
  CALL storage_done()
  
END PROGRAM
