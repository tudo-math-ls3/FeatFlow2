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
!#     call to xxx_app(rparlist), where xxx stands for the application
!#     name and user-defined parameters are supplied via parameter list.
!#
!# TODO
!#
!# 1.) Remove the old splib (ILU-k) from the solver module and replace it
!#     by the routines from the iluk.f90 module.
!#
!# 2.) Prelimiting does not work for FEM-FCT algorithms except for the
!#     linearised version. We still need to think about how to assemble
!#     the fluxes for prelimiting.
!#
!# 3.) FEM-FCT algorithms for Euler model are not working except for the
!#     linearised version.
!#
!# 4.) Jacobian matrix for the semi-implicit FEM-FCT algorithm has to be
!#     fixed (it was based on the fluxes but now it has to be based on the
!#     edgewise correciton factors).
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

  use euler_application
  use transport_application
  use zpinch_application
  use eulerlagrange_application

  implicit none

  ! global parameter list
  type(t_parlist) :: rparlist

  ! local constants
  character(LEN=*), parameter  :: sversion = 'VERSION'
  character(LEN=*), parameter  :: sbuild = 'BUILD'

  ! local variables
  character(LEN=SYS_STRLEN) :: cbuffer, hostname, hosttype, username
  character(LEN=SYS_STRLEN) :: application, sparameterfileName
  character(LEN=10) :: stime
  character(LEN=8)  :: sdate
  integer, external :: signal_SIGINT, signal_SIGQUIT


  ! Initialize Feat2 subsystem
  call system_init()

  ! Set system halt mode
  sys_haltmode = SYS_HALT_THROWFPE

  ! Initialize the output system
  call date_and_time(sdate, stime)
  call output_init('./log/flagship_'//sdate//'_'//stime(1:10)//'.log')

  ! Initialize storage subsystem
  call storage_init(500, 100)

  ! Initialize function parser
  call fparser_init()

  ! Initialize signal handler for SIGINT and SIGQUIT
  call fsignal(SIGINT, signal_SIGINT)
  call fsignal(SIGQUIT, signal_SIGQUIT)

  ! Print welcome screen
  call output_lbrk()
  call output_separator(OU_SEP_STAR)
  call output_line('  FLAGSHIP: Version ' // sversion // ',   Build ' // sbuild)
  call output_line('            Date '//sdate(7:8)//'.'//sdate(5:6)//'.'//sdate(1:4)//&
                   ', Time '//stime(1:2)//':'//stime(3:4)//':'//stime(5:6))
  call output_separator(OU_SEP_STAR)
  call output_line('  FlAGSHiP: Flux-corrected Aerodynamics by Galerkin')
  call output_line('            Schemes with High Performance (2004-2009)')
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

  ! Get command line arguments
  if (command_argument_count() .eq. 0) then
    call output_lbrk()
    call output_line('  PARAMETERFILE missing!!!')
    call output_lbrk()
    call sys_halt()
  end if

  ! Initialize parameter list from file
  call get_command_argument(command_argument_count(), cbuffer)
  sparameterfileName = adjustl(cbuffer)
  call parlst_init(rparlist)
  call parlst_readfromfile(rparlist, trim(sparameterfileName))
  call parlst_getvalue_string(rparlist, '', "application", application)
  call sys_tolower(application)
  call output_line('  Application:     '//trim(application))
  call output_line('  Parameterfile:   '//trim(sparameterfileName))
  call output_separator(OU_SEP_STAR)
  call output_lbrk()
  call output_line('The following settings are used for simulation')
  call output_separator(OU_SEP_MINUS)
  call parlst_info(rparlist)
  call output_separator(OU_SEP_MINUS)


  ! Switch to application module
  if (trim(application) .eq. 'transport') then
    call transp_app(rparlist, 'transport')

  elseif (trim(application) .eq. 'euler') then
    call euler_app(rparlist, 'euler')

  elseif (trim(application) .eq. 'zpinch') then
    call zpinch_app(rparlist, 'zpinch')

  elseif (trim(application) .eq. 'eulerlagrange') then
    call eulerlagrange_app(rparlist, 'eulerlagrange')

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
