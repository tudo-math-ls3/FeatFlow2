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
  character(LEN=256) :: sLogDir,sLogFile,sDatFile
  integer :: itest
  
  ! The very first thing in every application:
  ! Initialise system-wide settings:
  call sys_init()
  
  ! The very second thing in every program:
  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Read in parameter list
  call parlst_init(rparam)
  if (sys_getenv_string("DATFILE",sDatFile)) then
    call parlst_readfromfile(rparam, trim(sDatFile))
  else
    call parlst_readfromfile(rparam, './data/elemdbg.dat')
  end if

  ! Initialise the output system.
  !
  ! Normally, we write all the output to the screen and to a file.
  ! In the case that environment variables "$logdir"/"$resultsfile" exists,
  ! we write all the output to that file. This can be used e.g. in
  ! regression tests to compare results to reference results.
  if (sys_getenv_string("LOGDIR",slogdir) .and. &
      sys_getenv_string("RESULTFILE",slogfile)) then
    call output_init (trim(slogdir)//"/"//trim(slogfile))
  else
    ! Get log file name
    call parlst_getvalue_string(rparam, '', 'SLOGFILE', sLogFile, '')
    
    ! Initialise the output system.
    call output_init (sLogFile)
  end if
  
  ! Get config section name
  call parlst_getvalue_string(rparam, '', 'SCONFIGSECTION', sConfigSection, '')
  
  ! Get dimension
  call parlst_getvalue_int(rparam, sConfigSection, 'ITEST', itest, -1)

  ! Call the corresponding debugger
  select case(itest)
  case(101,102,103)
    ! 1D Element-Debugger #1
    call elemdbg1d_1(rparam,sConfigSection,itest)

  case(201,202,203)
    ! 2D Element-Debugger #1
    call elemdbg2d_1(rparam,sConfigSection,itest)

  case(301,302,303)
    ! 3D Element-Debugger #1
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
