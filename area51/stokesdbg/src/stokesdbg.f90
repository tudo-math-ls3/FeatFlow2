!##############################################################################
!# ****************************************************************************
!# <name> stokesdbg </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# </purpose>
!##############################################################################

program stokesdbg

use fsystem
use storage
use genoutput
use paramlist

use stokesdbg2d

implicit none

type(t_parlist) :: rparam
character(LEN=256) :: sLogFile, sDatFile
integer :: idriver
  
  ! The very first thing in every application:
  ! Initialise system-wide settings:
  call sys_init()
  
  ! The very second thing in every program:
  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Fetch parameter file name from command line
  sDatFile = './data/stokesdbg.dat'
  if(sys_ncommandLineArgs() .ge. 1) &
    call getarg(1, sDatFile)

  ! Read in parameter list
  call parlst_init(rparam)
  call parlst_readfromfile(rparam, sDatFile)

  ! Get log file name
  call parlst_getvalue_string(rparam, '', 'SLOGFILE', sLogFile, '')

  ! Initialise the output system.
  call output_init (sLogFile)
  
  ! Get dimension
  call parlst_getvalue_int(rparam, '', 'DRIVER', idriver, -1)

  ! Call the corresponding debugger
  select case(idriver)
  case(2000:2999)
    ! 2D Stokes Debugger
    call stokesdbg_run2d(rparam, idriver)
  
  case default
    ! Error
    call output_line('Invalid DRIVER parameter', OU_CLASS_ERROR, OU_MODE_STD, 'stokesdbg')
    
  end select
  
  ! Release parameter list
  call parlst_done (rparam)

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk ()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()

end program
