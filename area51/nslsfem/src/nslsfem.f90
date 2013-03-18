!##############################################################################
!# ****************************************************************************
!# <name> nslsfem </name>
!# ****************************************************************************
!# <purpose>
!# This program solves the 2D Navier-Stokes (NS) eq. using LSFEM.
!#
!# The second-order elliptic NS equations are reformulated into 
!# first-order equations using:
!# a- the definition of the vorticity:
!#   vorticity:  w = curl(u).
!#    http://www.student.math.uwaterloo.ca/
!#                 ~amat361/Fluid%20Mechanics/topics/vorticity.htm
!# The resulting system called Vorticity-Velocity-Pressure (VVP). 
!# 
!# b- the definition of the stresses:
!#   stress: -pI + \nu \nabla(\bu)
!# The resulting system in this case is Stress-Velocity-Pressure
!#   (SVP).
!# The problem is solved in a coupled manner for the solution of:
!#   a1- velocity component    u1, u2 (in 2D)
!#   a2- pressure              p
!#   a3- vorticity function    w
!#          OR
!#   b1- velocity component    u1, u2 (in 2D)
!#   b2- pressure              p
!#   b3-stress component       s1, s2, s3 (in 2D)
!# variables.
!#
!# </purpose>
!##############################################################################

program nslsfem

  use LS_NS_VVP_MG2D
  use LS_NS_VVP_Time_MG2D
  
  use LS_NS_SVP_MG2D
  use LS_NS_SVP_RT2D
  use LS_NS_SVP
  
  use LS_LS
  
  implicit none
  
  ! local variables
  character(len=SYS_STRLEN) :: slogdir,slogfile
  
  ! The very first thing in every application:
  ! Initialise system-wide settings:
  
  call system_init()

  ! Initialise the output system.
  !
  ! Normally, we write all the output to the screen and to a file
  ! './log/output.txt'.
  ! In the case that environment variables "$logdir"/"$resultsfile" exists,
  ! we write all the output to that file. This can be used e.g. in
  ! regression tests to compare results to reference results.
  if (sys_getenv_string('LOGDIR',slogdir) .and. &
      sys_getenv_string('RESULTFILE',slogfile)) then
    call output_init (trim(slogdir)//'/'//trim(slogfile))
  else
    call output_init ('./log/output.txt','','./log/error.txt')
  end if

  ! The very second thing in every program:
  ! Initialise the storage management:
  !
  ! 2.) Initialise FEAT 2.0 storage management:
  call storage_init(999, 100)

!  ! Call the problem to solve 2D Navier-stokes:
!  call output_lbrk()
!  call output_line('Calculating 2D Navier-Stokes-LSFEM')
!  call output_lbrk()
!  call output_line('Vorticity-Velocity-Pressure-Multigrid')  
!  call output_line('-------------------------------------')
!  call ls_vvp_mg2d

!  ! Call the problem to solve 2D transient Navier-stokes:
!  call output_lbrk()
!  call output_line('Calculating 2D transient Navier-Stokes-LSFEM')
!  call output_lbrk()
!  call output_line('Vorticity-Velocity-Pressure-Multigrid')
!  call ls_vvp_time_mg2d


  ! Call the problem to solve 2D Navier-stokes:
  call output_lbrk()
  call output_line('Calculating 2D Navier-Stokes-LSFEM')
  call output_lbrk()
  call output_line('Stress-Velocity-Pressure-Multigrid')  
  call output_line('----------------------------------')
  call ls_svp_mg2d


!  ! Call the problem to solve 2D Navier-stokes:
!  call output_lbrk()
!  call output_line('Calculating 2D Navier-Stokes-LSFEM')
!  call output_lbrk()
!  call output_line('Stress-Velocity-Pressure-Multigrid')  
!  call output_line('       Defect-based Version       ') 
!  call output_line('----------------------------------')
!  call ls_svp_2d

!  ! Call the problem to solve 2D Navier-stokes:
!  call output_lbrk()
!  call output_line('Calculating 2D Navier-Stokes-LSFEM')
!  call output_lbrk()
!  call output_line('Stress-Velocity-Pressure-Raviart-Thomas')  
!  call output_line('---------------------------------------')
!  call ls_svp_rt2d

!  ! Call the problem to solve 2D Navier-stokes:
!  call output_lbrk()
!  call output_line('Level Set 2D-LSFEM')
!  call output_line('------------------')
!  call ls_levelset

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program
