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

program burgers1d

  use burgers1d_method5
  use burgers1d_method6
  
  implicit none
  
  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  call system_init()

  ! The very second thing in every program: 
  ! Initialise the storage management: 
  !
  ! 2.) Initialise FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Call the problem to solve. 
  print *
  print *,'Calculating Burgers1D-Problem with method 5'
  print *,'-------------------------------------------'
  call burgers1d5

  ! Call the problem to solve. 
  print *
  print *,'Calculating Burgers1D-Problem with method 6'
  print *,'-------------------------------------------'
  call burgers1d6

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  print *
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
