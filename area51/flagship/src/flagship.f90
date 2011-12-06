!##############################################################################
!# ****************************************************************************
!# <name> flagship </name>
!# ****************************************************************************
!#
!# <purpose>
!# This is the main program which calls the individual application modules.
!#
!# Frequently asked questions
!# --------------------------
!#
!# 1.) How to add a new application/model?
!#
!#     To add a new application, you should create a new subdirectory in
!#     src/models giving it the name of your application. Then you should
!#     implement your application so that it can be evoked by a single
!#     call to xxx_app(rparlist,'appname'), where xxx stands for the
!#     application name and user-defined parameters are supplied via
!#     the parameter list rparlist.
!#
!# </purpose>
!##############################################################################

program flagship

  use fparser
  use fsystem
  use genoutput
  use paramlist
  use signals
  use storage

  use flagship_basic
  use hydro_application
  use mhd_application
  use transport_application
  use zpinch_application

  implicit none

  ! global parameter list
  type(t_parlist) :: rparlist

  ! local variables
  character(LEN=SYS_STRLEN) :: cbuffer, hostname, hosttype, username
  character(LEN=SYS_STRLEN) :: slogdir, slogfile
  character(LEN=SYS_STRLEN) :: application
  character(LEN=SYS_STRLEN) :: sparameterfile, sperfconfigfile
  character(LEN=10) :: stime
  character(LEN=8)  :: sdate
  integer, external :: signal_SIGINT, signal_SIGQUIT


#ifdef ENABLE_COPROCESSOR_SUPPORT
  ! Initialise CUDA subsystem
  call coproc_init()
#endif

  ! Initialise Feat2 subsystem
  call system_init()

  ! Set system halt mode
#ifdef ENABLE_ERROR_TRACEBACK
  sys_haltmode = SYS_HALT_THROWFPE
#else
  sys_haltmode = SYS_HALT_STOP
#endif

  ! Initialise the output system
  ! Use $LOGDIR/$LOGFILE if set, otherwise hardcoded setting based on current time
  call date_and_time(sdate, stime)
  call getenv('LOGDIR',cbuffer); slogdir = adjustl(cbuffer)
  call getenv('LOGFILE',cbuffer); slogfile = adjustl(cbuffer)
  if (trim(slogdir) .eq. '') then
    call output_init('./log/flagship_'//sdate//'_'//stime(1:10)//'.log')
  else
    if (trim(slogfile) .eq. '') then
      slogfile = 'flagship_' // sdate // '_' // stime(1:10) // '.log'
    end if
    call output_init(trim(slogdir) // '/' // trim(slogfile))
  end if

  ! Initialise storage subsystem
  call storage_init(500, 100)

  ! Initialise function parser
  call fparser_init()

  ! Initialise signal handler for SIGINT and SIGQUIT
  call fsignal(SIGINT, signal_SIGINT)
  call fsignal(SIGQUIT, signal_SIGQUIT)

  ! Print welcome screen
  call output_lbrk()
  call output_separator(OU_SEP_STAR)
  call output_line('  FlAGSHiP: Flux-corrected Aerodynamics by Galerkin')
  call output_line('            Schemes with High Performance (2004-2011)')
  call output_lbrk()
  call output_line('  Authors:  Dmitri Kuzmin, Matthias Moeller')
  call output_line('            Institute of Applied Mathematics')
  call output_line('            Dortmund University of Technology')
  call output_line('            Vogelpothsweg 87, 44227 Dortmund, Germany')
  call output_separator(OU_SEP_STAR)
  call getenv('HOST',cbuffer); hostname = adjustl(cbuffer)
  call output_line('  Hostname:        '//trim(hostname))
  call getenv('HOSTTYPE',cbuffer); hosttype = adjustl(cbuffer)
  call output_line('  Hosttype:        '//trim(hosttype))
  call getenv('USER',cbuffer); username = adjustl(cbuffer)
  call output_line('  Username:        '//trim(username))
  call output_line('  Date:            '//sdate(7:8)//'.'//sdate(5:6)//'.'//sdate(1:4))
  call output_line('  Time:            '//stime(1:2)//':'//stime(3:4)//':'//stime(5:6))

  ! Get command line arguments
  if (command_argument_count() .eq. 0) then
    call output_lbrk()
    call output_line('  PARAMETERFILE missing!!!')
    call output_lbrk()
    call sys_halt()
  end if

  ! Initialise parameter list from file
  call get_command_argument(command_argument_count(), cbuffer)
  sparameterfile = adjustl(cbuffer)
  call parlst_init(rparlist)
  call parlst_readfromfile(rparlist, trim(sparameterfile))
  call parlst_getvalue_string(rparlist, '', "application", application)
  call sys_tolower(application)
  call output_line('  Application:     '//trim(application))
  call output_line('  Parameterfile:   '//trim(sparameterfile))
  call output_separator(OU_SEP_STAR)
  call output_lbrk()
  call output_line('The following settings are used for simulation')
  call output_separator(OU_SEP_MINUS)
  call parlst_info(rparlist)
  call output_separator(OU_SEP_MINUS)

  ! Initialise global performace configurations
  call parlst_getvalue_string(rparlist, '', 'sperfconfigfile', sperfconfigfile)
  call flagship_initPerfConfig(sperfconfigfile)

  ! Switch to application module
  if (trim(application) .eq. 'transport') then
    call transp_app(rparlist, 'transport')

  elseif (trim(application) .eq. 'hydro') then
    call hydro_app(rparlist, 'hydro')

  elseif (trim(application) .eq. 'zpinch') then
    call zpinch_app(rparlist, 'zpinch')

  elseif (trim(application) .eq. 'mhd') then
    call mhd_app(rparlist, 'mhd')
  else
    call output_line('Invalid application name!',&
        OU_CLASS_WARNING,OU_MODE_STD,'flagship')
    call sys_halt()
    
  end if

  ! Release parameter list
  call parlst_done(rparlist)

  ! Release function parser
  call fparser_done()

  ! Release storage
  call storage_info(.true.)
  call storage_done()
  call output_lbrk()

  ! Close logfile
  call output_done()

end program flagship

!*****************************************************************************

!<function>

function signal_SIGINT(signum) result(sigcount)

  use fsystem
  use genoutput
  use signals

!<description>
  ! This subroutine performs signal handling for SIGINT. In essence,
  ! it counts the number if SIGINT`s received and terminates if user
  ! sent SIGINT more than three times.
!</description>

!<input>
  integer, intent(in) :: signum
!</input>

!<result>
  ! signal
  integer :: sigcount
!</result>
!</function>

  ! local variables
  integer, save :: icount = 0

  sigcount = icount

  if (signum .eq. -1) then

    ! Reset counter
    icount = 0

  elseif (signum .eq. SIGINT) then

    ! Increase counter
    icount = icount+1
    if (icount .ge. 3) then
      call output_line('Simulation terminated due to user interruption (SIGINT)')
      call sys_halt()
    end if

  end if
end function signal_SIGINT

!*****************************************************************************

!<function>

function signal_SIGQUIT(signum) result(sigcount)

  use fsystem
  use genoutput
  use signals

!<description>
  ! This subroutine performs signal handling for SIGQUIT.
!</description>

!<input>
  integer, intent(in) :: signum
!</input>

!<result>
  ! signal
  integer :: sigcount
!</result>
!</function>

  call output_line('Simulation terminated due to user interruption (SIGQUIT)')
  stop
end function signal_SIGQUIT
