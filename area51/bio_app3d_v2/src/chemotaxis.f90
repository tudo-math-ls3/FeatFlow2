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

!   use chemotaxis_coupldRHSnontrivialtest
!    use chemotaxis_coupldRHSadapttest
!   use chemotaxis_coupldRHSconst
!   use chemotaxis_coupldRHS_ss_l2
!   use chemotaxis_coupldRHS
!   use chemotaxis_coupldRHStimeconsttest
!   use chemotaxis_coupld1_5
!   use heatcond_method1
!   use chemotaxis_implicit
!   use chemotaxis_coupld
!   use heatcond_method5
!   use heatcond_1
!   use chemotaxis_coupld2
!   use heatcond_methodExpl
!    use chemotaxis_cherkur
!     use chemotaxis_cherkur_corner
!     use chemotaxis_tystlev
!     use chemotaxis_cherkur_3
!     use chemotaxis_cherkur_opt
!     use chemotaxis_cherkur_bilf
!     use chemotaxis_cherkur_bilf_stat
!     use chemotaxis_cherkur_bilf_nonlin
!     use chemotaxis_cherkur_bilf_nonlin_cpct
!     use chemotaxis_cherkur_bilf_nonlin_mdl3RHS
!     use chemotaxis_cherkur_bilf_nonlin_mdl3
!      use chemotaxis_cherkur_TVD_RHS
!    use chemotaxis_cherkur_bilf_nonlin_cpct_RHS
!     use chemotaxis_cherkur_TVD_l2_test_RHS
!     use chemotaxis_cherkur_template_lumped
!     use chemotaxis_cherkur_template_star
      use chemotaxis_pattern_FCT

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
  call output_line('Calculating a simple implicit chemotaxis-Problem')
  call output_line('------------------------------------------')
    print*, "..."
!  call chemotaxischerkurbilfnonlincpctrhs
     call chemotaxispatternFCT
 
  ! Call the problem to solve - method 5 = multigrid
 ! call output_lbrk()
 ! call output_line('Calculating heatcond-Problem with method 5')
 ! call output_line('------------------------------------------')
 ! call heatcond5

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
!   call output_lbrk()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
