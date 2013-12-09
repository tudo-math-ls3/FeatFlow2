!##############################################################################
!# ****************************************************************************
!# <name> cc3d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising a simple
!# stationary Navier-Stokes equation
!#
!#              $$- \nu Laplace(u) + u*grad(u) + \Nabla p = f $$
!#              $$ \Nable \cdot p = 0$$
!#
!# on a 3D domain for a 3D function $u=(u_1,u_2,u_3)$ and a pressure $p$.
!#
!# A linear block system with 4 solution components is set up,
!# discretised and solved.
!#
!# </purpose>
!##############################################################################

program cc3d

  use ccmainproblem
  
  implicit none
  
  character(LEN=SYS_STRLEN) :: slogfile,serrorfile
  
  ! The very first thing in every application:
  ! Initialise system-wide settings:
  call sys_init()
  
  ! General output init - temporary until we read in the output settings
  call output_init ()
  
  ! Read the name of the message and error log file.
  call cc_getLogFiles (slogfile,serrorfile)
  
  ! Release output stuff
  call output_done()
  
  ! Initialise log file for output.
  call output_init (slogfile,serrorfile)
  
  ! Now we can really start!
  !
  ! Initialise the storage management:
  call storage_init(999, 100)
  
  ! Initialise the parser
  call fparser_init ()

  ! Call the problem to solve.
  call output_lbrk ()
  call output_line ('Calculating cc3d-Problem')
  call output_separator (OU_SEP_MINUS)
  
  call cc3dmain ()

  ! Release the parser
  call fparser_done ()

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk ()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
