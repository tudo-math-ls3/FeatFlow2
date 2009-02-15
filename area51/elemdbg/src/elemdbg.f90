!##############################################################################
!# ****************************************************************************
!# <name> elemdbg </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# </purpose>
!##############################################################################

program elemdbg
   
  use elemdbg1d_test1
  use elemdbg2d_test1
  use elemdbg3d_test1
  use paramlist
  
  implicit none
  
  type(t_parlist) :: rparam
  character(LEN=64) :: sConfigSection
  integer :: itest
  
  ! The very first thing in every application: 
  ! Initialise system-wide settings:
  
  call system_init()
  
  ! Initialise the output system. Write the program output to screen as
  ! well as to the file 'log/output.txt'.
  call output_init ('./log/output.txt')

  ! The very second thing in every program: 
  ! Initialise the FEAT 2.0 storage management: 
  call storage_init(999, 100)
  
  ! Read in parameter list
  call parlst_init(rparam)
  call parlst_readfromfile(rparam, './data/elemdbg.dat')
  
  ! Get config section name
  call parlst_getvalue_string(rparam, '', 'SCONFIGSECTION', sConfigSection, '')
  
  ! Get dimension
  call parlst_getvalue_int(rparam, sConfigSection, 'ITEST', itest, -1)

  ! Call the corresponding debugger
  select case(itest)
  case(101,102)
    ! 1D Element-Debugger, test 1
    call elemdbg1d_1(rparam,sConfigSection,itest)

  case(201,202)
    ! 2D Element-Debugger, test 1
    call elemdbg2d_1(rparam,sConfigSection,itest)

  case(301,302)
    ! 3D Element-Debugger, test 1
    call elemdbg3d_1(rparam,sConfigSection,itest)
  
  case default
    ! Error
    call output_line('Invalid ITEST parameter or config section name',&
                     OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg')
    
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
