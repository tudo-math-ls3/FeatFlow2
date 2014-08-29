!##############################################################################
!# ****************************************************************************
!# <name> meshadapt </name>
!# ****************************************************************************
!#
!# <purpose>
!# This application performs local mesh adaptation based on a given
!# indicator function and prescribed refinement/recoarsening tolerances.
!#
!# Calling example:
!#
!# meshadapt -read2d mymesh -indicator myindicator.dat -refmax 3 \
!#           -reftol 0.8 -crstol 0.2
!#
!# The initial 2D mesh is read from the TRI/PRM file mymesh.tri/prm
!# and refined based on the error element-wise indicator given in file
!# myerror.dat. Cells with error larger than 0.8 are refined and cells
!# with error smalled than 0.2 are re-coarsened. The maximum
!# refinement level is 3, that is, an element of the initial grid can
!# be refined 3 times at most.
!#
!# </purpose>
!##############################################################################

program meshadapt

  use fsystem
  use genoutput
  use signals
  use storage
  
  use meshadaptbase

  implicit none
  
  ! Initialise system-wide settings
  call sys_init()
  
  ! Initialise the output system
  call output_init()

  ! Initialise the FEAT 2.0 storage management
  call storage_init(100, 100)
  
  if (madapt_signalhandler(SIGUSR1) .eq. 0) then
    
    ! We are in daemon mode, hence, register signal handler
    call fsignal(SIGUSR2, madapt_signalhandler) ! export mesh
    call fsignal(SIGINT,  madapt_signalhandler) ! single step
    call fsignal(SIGQUIT, madapt_signalhandler) ! finalise
    call fsignal(SIGHUP,  madapt_signalhandler) ! finalise
    call fsignal(SIGTERM, madapt_signalhandler) ! finalise

    daemon: do
      ! Perform mesh adaptation step?
      call fsignal(SIGINT, madapt_signalhandler)

      ! Export mesh to file?
      call fsignal(SIGUSR2, madapt_signalhandler)

      ! Finalise and exit?
      call fsignal(SIGQUIT, madapt_signalhandler)
      call fsignal(SIGTERM, madapt_signalhandler)
      call fsignal(SIGHUP,  madapt_signalhandler)
    end do daemon

  else
    
    ! We are not in daemon mode, hence, perform a single mesh
    ! adaptation step, export the mesh to file and exit
    if (madapt_signalhandler(SIGINT) .ne. 0) then
      call output_line("An error occured during mesh adaptation!",&
          OU_CLASS_ERROR,OU_MODE_STD,"meshadapt")
      call sys_halt()
    end if

    if (madapt_signalhandler(SIGUSR2) .ne. 0) then
      call output_line("An error occured while exporting mesh to file!",&
          OU_CLASS_ERROR,OU_MODE_STD,"meshadapt")
      call sys_halt()
    end if

    if (madapt_signalhandler(SIGQUIT) .ne. 0) then
      call output_line("An error occured during finalisation!",&
          OU_CLASS_ERROR,OU_MODE_STD,"meshadapt")
      call sys_halt()
    end if
  end if
  
  ! Clean up the storage management, finish
  call storage_done()

end program meshadapt
