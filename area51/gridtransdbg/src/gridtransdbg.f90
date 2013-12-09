!##############################################################################
!# ****************************************************************************
!# <name> gridtransdbg </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# </purpose>
!##############################################################################

program gridtransdbg

  use fsystem
  use genoutput
  use storage
  use paramlist
  
  use gridtransdbg_test1

  implicit none
  
  type(t_parlist) :: rparam
  character(LEN=64) :: sConfigSection
  character(LEN=256) :: sLogFile
  integer :: itest
  
  ! The very first thing in every application:
  ! Initialise system-wide settings:
  call sys_init()
  
  ! The very second thing in every program:
  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Read in parameter list
  call parlst_init(rparam)
  call parlst_readfromfile(rparam, './data/gridtransdbg.dat')

  ! Get log file name
  call parlst_getvalue_string(rparam, '', 'SLOGFILE', sLogFile, '')

  ! Initialise the output system.
  call output_init (sLogFile)
  
  ! Get config section name
  call parlst_getvalue_string(rparam, '', 'SCONFIGSECTION', sConfigSection, '')
  
  ! Get dimension
  call parlst_getvalue_int(rparam, sConfigSection, 'ITEST', itest, -1)

  ! Call the corresponding debugger
  select case(itest)
  case(1)
    ! Grid-Transfer Debugger #1
    call gridtransdbg_1(rparam, sConfigSection, itest)
 
  case default
    ! Error
    call output_line('Invalid ITEST parameter or config section name',&
                     OU_CLASS_ERROR, OU_MODE_STD, 'gridtransdbg')
    
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