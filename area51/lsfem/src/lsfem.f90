!##############################################################################
!# ****************************************************************************
!# <name> lsfem </name>
!# ****************************************************************************
!# <purpose>
!#
!# This package solves 2D problems using the least-squares FEM.
!# The following PDEs are considered in this code:
!#  - Advection-diffusion-reaction equations (Poisson, diffusion-reaction)
!#  - Navier-Stokes equations (vorticity-based, stress-based)
!#  - non-Newtonian fluid flows (stress-based)
!#
!# </purpose> 
!##############################################################################

program lsfem

  use LS_Poisson_MG2D

  use LS_NS_VVP_MG2D
  use LS_NS_VVP_Time_MG2D
  
  use LS_NS_SVP_MG2D
  use LS_NS_SVP_MG_MAT2D
  use LS_NS_SVP_RT2D
  use LS_NS_SVP
  use LS_SVP_NNew_MG2D
  
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


  ! Call the problem to solve 2D Poisson equation:
  call output_lbrk()
  call output_line('Calculating 2D Poisson')
  call output_lbrk()
  call output_line('Div-Grad-(Curl)-Multigrid')  
  call output_line('-------------------------')
  call LS_Poisson_mg

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


!  ! Call the problem to solve 2D Navier-stokes:
!  call output_lbrk()
!  call output_line('----------------------------------')  
!  call output_line('Calculating 2D Navier-Stokes-LSFEM')
!  call output_line('Stress-Velocity-Pressure-Multigrid')  
!  call output_line('----------------------------------')
!  call ls_svp_mg2d


!  ! Call the problem to solve 2D Navier-stokes:
!  call output_lbrk()
!  call output_line('----------------------------------')  
!  call output_line('Calculating 2D Navier-Stokes-LSFEM')
!  call output_line('Stress-Velocity-Pressure-Multigrid') 
!  call output_line('     Matrix-Based Prol./Rest.     ') 
!  call output_line('----------------------------------')
!  call ls_svp_mg_mat2d


!  ! Call the problem to solve Non-Newtonian fluid flow:
!  call output_lbrk()
!  call output_line('----------------------------------')  
!  call output_line('Calculating 2D Navier-Stokes-LSFEM')
!  call output_line('Stress-Velocity-Pressure-Multigrid') 
!  call output_line('       Non-Newtonian Fluid        ') 
!  call output_line('----------------------------------')
!  call ls_svp_nn_mg2d

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


!  ! Call the problem to solve level-set equation:
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
