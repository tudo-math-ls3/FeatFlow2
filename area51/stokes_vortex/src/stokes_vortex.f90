!##############################################################################
!# ****************************************************************************
!# <name> stokes_vortex </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a benchmark for solving the Stokes equation with a steady
!# standing vortex solution.
!# </purpose>
!##############################################################################

program stokes_vortex

use fsystem
use genoutput
use storage
use stokes2d_vortex_slip
use stokes2d_vortex_noslip
use stokes2d_no_flow

implicit none

! local variables
character(len=SYS_STRLEN) :: slogdir,slogfile

  ! Initialise system-wide settings:
  call system_init()

  ! Initialise the output system.
  if (sys_getenv_string('LOGDIR',slogdir) .and. &
      sys_getenv_string('RESULTFILE',slogfile)) then
    call output_init (trim(slogdir)//'/'//trim(slogfile))
  else
    call output_init ('./log/output.txt')
  end if

  ! Initialise FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Run Slip-BC driver
  call output_separator(OU_SEP_STAR)
  call output_line('2D Stokes: Vortex Slip BCs')
  call output_line('--------------------------')
  call stokes2d_vtx_slip
  call output_lbrk()

  ! Run No-Slip-BC driver
  call output_separator(OU_SEP_STAR)
  call output_line('2D Stokes: Vortex No-Slip BCs')
  call output_line('-----------------------------')
  call stokes2d_vtx_noslip
  call output_lbrk()

  ! Run No-Flow driver
  call output_separator(OU_SEP_STAR)
  call output_line('2D Stokes: No-Flow')
  call output_line('------------------')
  call stokes2d_noflow
  call output_lbrk()

  ! Print out heap statistics
  call output_lbrk()
  call storage_info(.true.)

  ! Clean up the storage management, finish
  call storage_done()

end program
