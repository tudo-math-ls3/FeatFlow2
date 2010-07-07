!##############################################################################
!# ****************************************************************************
!# <name> chemotaxis </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the equation
!#
!#   d/dt u  - div ( alpha*grad(u) )  +  beta*grad(u)  +  gamma*u  =  f
!#
!# on a 2D domain for a scalar function u=u(x,t).
!#
!# The first example (module heatcond_method5) is a very simple one-routine
!# solver for the head condition equation without any specials. The equation
!# is discretised with explicit Euler for a fixed number of time steps.
!# Boundary conditions and RHS are constant in time.
!#
!# The second example (module heatcond_method5) separates the different steps
!# of the solution process into different subroutines.
!# The communication is done using a problem-related structure. For the
!# communication with callback routines during the assembly, a
!# collection structure is set up. As solver, a simple Multigrid solver is
!# used.
!#
!# The equation itself is set up using parameters of a .DAT file
!# named 'data/heatcond.dat'. Here all parameters (alpha, beta, gamma, level,...)
!# can be found and changed by the user without recompiling the program.
!# </purpose>
!##############################################################################

program chemotaxis

      use chemotaxis_blowup

  implicit none

  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  call system_init()

  ! The very second thing in every program: 
  ! Initialise the storage management: 
  !
  ! 2.) Initialise FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Call the problem to solve - explicit method  = very very basic
!   call output_lbrk()
!   call output_line('Calculating heatcond-Problem with explicit euler method')
!   call output_line('------------------------------------------')
!    call  heatcond1
!   call output_line('....Skip this model....')

  
!   Call the problem to solve - chemotaxis ( very very basic )
  call output_lbrk()
  call output_line('Calculating chemotaxis_blowup')
  call output_line('------------------------------------------')
  print*, "..."
  call chemotaxisBlowup


  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
!   call output_lbrk()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
