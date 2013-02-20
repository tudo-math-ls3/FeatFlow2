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

#include "flagship.h"

!$use omp_lib
  use fparser
  use fsystem
  use genoutput
  use io
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
  character(LEN=SYS_STRLEN) :: application,cbuffer
  character(LEN=SYS_STRLEN) :: slogfile,serrorfile,sbenchlogfile
  character(LEN=SYS_STRLEN) :: sparameterfile,sperfconfigfile
  character(LEN=10) :: stime
  character(LEN=8)  :: sdate
  integer, external :: signal_SIGINT, signal_SIGQUIT
  integer :: iunit
  logical :: bexit


#ifdef ENABLE_COPROCESSOR_SUPPORT
  ! Initialise CUDA subsystem
  call coproc_init(0)
#endif

  ! Initialise Feat2 subsystem
  call system_init()

  ! Set system halt mode
#ifdef ENABLE_ERROR_TRACEBACK
  sys_haltmode = SYS_HALT_THROWFPE
#else
  sys_haltmode = SYS_HALT_STOP
#endif

  ! Initialise the output system - temporary until we read in the output settings
  call output_init()

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
  call output_line(' /FFFFFFFF /LL  /AAAAAA   /GGGGGG   /SSSSSS  /HH   /HH /II /PPPPPPP ')
  call output_line('| FF_____/| LL /AA__  AA /GG__  GG /SS__  SS| HH  | HH|__/| PP__  PP')
  call output_line('| FF      | LL| AA  \ AA| GG  \__/| SS  \__/| HH  | HH /II| PP  \ PP')
  call output_line('| FFFFF   | LL| AAAAAAAA| GG /GGGG|  SSSSSS | HHHHHHHH| II| PPPPPPP/')
  call output_line('| FF__/   | LL| AA__  AA| GG|_  GG \____  SS| HH__  HH| II| PP____/ ')
  call output_line('| FF      | LL| AA  | AA| GG  \ GG /SS  \ SS| HH  | HH| II| PP      ')
  call output_line('| FF      | LL| AA  | AA|  GGGGGG/|  SSSSSS/| HH  | HH| II| PP      ')
  call output_line('|__/      |__/|__/  |__/ \______/  \______/ |__/  |__/|__/|__/      ')
  call output_lbrk()
  call output_line(' Flux-corrected Aerodynamics by Galerkin Schemes with High Performance')
    call output_separator(OU_SEP_STAR)
  call output_line(' Kernel  : Featflow')
  call output_line(' Version : 2.0')
  call output_line(' Web     : http://www.featflow.de')
  call output_lbrk()
  call output_line(' Author  : Matthias Moeller (2004-2013)')
  call output_line('           Institute of Applied Mathematics')
  call output_line('           Technical University of Dortmund')
  call output_line('           Vogelpothsweg 87')
  call output_line('           44227 Dortmund')
  call output_line('           Germany')
  call output_separator(OU_SEP_STAR)
  if (sys_getenv_string('HOST',cbuffer))&
      call output_line('  Hostname:        '//trim(adjustl(cbuffer)))
  if (sys_getenv_string('HOSTTYPE',cbuffer))&
      call output_line('  Hosttype:        '//trim(adjustl(cbuffer)))
  if (sys_getenv_string('USER',cbuffer))&
      call output_line('  Username:        '//trim(adjustl(cbuffer)))
  call date_and_time(sdate, stime)
  call output_line('  Date:            '//sdate(7:8)//'.'//sdate(5:6)//'.'//sdate(1:4))
  call output_line('  Time:            '//stime(1:2)//':'//stime(3:4)//':'//stime(5:6))

  ! Check presence of command line argument(s)
  bexit = (sys_ncommandLineArgs() .eq. 0)

  ! Check existence of parameter file
  if (.not.bexit) then
    call sys_getcommandLineArg(sys_ncommandLineArgs(), cbuffer)
    sparameterfile = trim(adjustl(cbuffer))
    iunit = sys_getFreeUnit()
    bexit = .not.(sys_fileExists(iunit,trim(sparameterfile)))
  end if
  
  ! Exit if application is not called correctly
  if (bexit) then
    call output_separator(OU_SEP_STAR)
    call output_lbrk()
    call output_line('Application must be called as:')
    call output_line('  appname <args> parameterfile')
    call output_lbrk()
    call output_line('Admissible arguments are as follows:')
    call output_line(' -D<section>.parameter:<entry>=value')
    call output_line('   Overwrites the entry of the parameter in section by the given value.')
    call output_line('   If <section> is missing then the default section is used.')
    call output_line('   If <entry> is missing then the first entry is used.')
    call output_lbrk()
    call exit(-1)
  end if

  ! Initialise parameter list from file
  call parlst_init(rparlist)
  call parlst_readfromfile(rparlist, trim(sparameterfile), bexpandVars=.false.)

  ! Check for command line arguments
  call flagship_parseCmdlArguments(rparlist)

  ! Expand all subvariables and environment variables to the actual values
  call parlst_expandEnvVariables(rparlist)
  call parlst_expandSubvars(rparlist)

  ! Get log files for output
  call parlst_getvalue_string (rparlist,'',&
      'slogfile', slogfile, '')
  call parlst_getvalue_string (rparlist,'',&
      'serrorfile', serrorfile, '')
  call parlst_getvalue_string (rparlist,'',&
      'sbenchlogfile', sbenchlogfile, '')
  call parlst_getvalue_string(rparlist, '',&
      'sperfconfigfile', sperfconfigfile)
  
  ! Release output stuff
  call output_done()

  ! Initialise log file for output
  call output_init (slogfile, serrorfile, sbenchlogfile)

  ! Initialise global performace configurations
  call flagship_initPerfConfig(sperfconfigfile)

  ! Switch to application module
  call parlst_getvalue_string(rparlist, '', "application", application, '')
  call sys_toupper(application)
  call output_line('  Application:     '//trim(application))
  call output_line('  Parameterfile:   '//trim(sparameterfile))
  call output_separator(OU_SEP_STAR)
  call output_lbrk()
  call output_line('The following settings are used for simulation')
  call output_separator(OU_SEP_MINUS)
  call parlst_info(rparlist)
  call output_separator(OU_SEP_MINUS)

  ! Switch to application module
  if (trim(application) .eq. 'TRANSPORT') then
    call transp_app(rparlist, 'TRANSPORT')

  elseif (trim(application) .eq. 'HYDRO') then
    call hydro_app(rparlist, 'HYDRO')

  elseif (trim(application) .eq. 'ZPINCH') then
    call zpinch_app(rparlist, 'ZPINCH')

  elseif (trim(application) .eq. 'MHD') then
    call mhd_app(rparlist, 'MHD')
    
  else
    call output_line('Invalid application name!',&
        OU_CLASS_ERROR,OU_MODE_STD,'flagship')
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

contains

!************************************************************************

!<subroutine>

  subroutine flagship_parseCmdlArguments(rparlist)

!<description>
    ! This subroutine parses the commandline arguments and modifies the
    ! parameter values in the global parameter list.
!</description>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(inout) :: rparlist
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    character(LEN=PARLST_MLSECTION) :: ssection
    character(LEN=PARLST_MLNAME) :: sparameter,sentry
    character(LEN=PARLST_MLDATA) :: svalue
    character(LEN=SYS_STRLEN) :: soption
    character(LEN=PARLST_MLDATA), dimension(:), pointer :: Svalues
    integer :: iarg,narg,iformat,itoken1,itoken2,isubstring,nsubstrings,i

    iarg = 1; narg = sys_ncommandLineArgs()-1

    cmdarg: do while(iarg .le. narg)
      ! Retrieve next command line argument
      call sys_getcommandLineArg(iarg,soption,svalue,iformat)

      select case(iformat)
      case (0)
        ! Options without parameter values

      case (1,2)
        ! Options with given parameter values

        ! What option are we?
        if (soption(1:1) .eq. 'D') then

          ! Tokenize string: soption=D<section>.variable:<entry>
          itoken1 = scan(soption,'.')
          if (itoken1 .ne. 0) then
            ssection = trim(adjustl(soption(2:itoken1-1)))
          else
            ssection = ''
          end if
          
          itoken2 = scan(soption,':')
          if (itoken2 .ne. 0) then
            sparameter = trim(adjustl(soption(max(2,itoken1+1):itoken2-1)))
            sentry     = trim(adjustl(soption(max(2,itoken2+1):)))
            read(sentry, fmt='(I10)') isubstring
          else
            sparameter = trim(adjustl(soption(max(2,itoken1+1):)))
            sentry     = ''
            isubstring = 0
          end if

          ! Query/add section in parameter list
          call parlst_querysection(rparlist, trim(adjustl(ssection)), p_rsection)
          if (.not.associated(p_rsection)) then
            call parlst_addsection (rparlist, trim(adjustl(ssection)))
            call parlst_querysection(rparlist, trim(adjustl(ssection)), p_rsection)
          end if

          ! We need to consider several cases
          if (parlst_queryvalue(p_rsection, trim(adjustl(sparameter))) .ne. 0) then

            ! Parameter exists already in section
            if (isubstring .eq. 0) then
              ! Overwrite existing value
              call parlst_addvalue(p_rsection, trim(adjustl(sparameter)),&
                  trim(adjustl(svalue)))
            else
              ! Get number of existing substrings
              nsubstrings = parlst_querysubstrings(p_rsection, trim(adjustl(sparameter)))
              if (isubstring .lt. nsubstrings) then
                ! Overwrite existing substring
                call parlst_setvalue(p_rsection, trim(adjustl(sparameter)),&
                    trim(adjustl(svalue)), isubstring=isubstring)
              else
                ! Make a backup of existing substrings
                allocate(Svalues(0:nsubstrings))
                do i=0,nsubstrings
                  call parlst_getvalue_string(p_rsection, trim(adjustl(sparameter)),&
                      Svalues(i), isubstring=i)
                end do
                ! Increase the number of substrings
                call parlst_addvalue(p_rsection, trim(adjustl(sparameter)),&
                    trim(adjustl(svalue)), isubstring)
                ! Copy existing substrings
                do i=0,nsubstrings
                  call parlst_setvalue(p_rsection, trim(adjustl(sparameter)),&
                      Svalues(i), isubstring=i)
                end do
                ! Add new substring
                call parlst_setvalue(p_rsection, trim(adjustl(sparameter)),&
                    trim(adjustl(svalue)), isubstring=isubstring)
                deallocate(Svalues)
              end if
            end if
          else
            ! Add new value to parameter list
            if (isubstring .eq. 0) then
              call parlst_addvalue(p_rsection, trim(adjustl(sparameter)),&
                  trim(adjustl(svalue)))
            else
              call parlst_addvalue(p_rsection, trim(adjustl(sparameter)),&
                  trim(adjustl(svalue)), isubstring)
              call parlst_setvalue(p_rsection, trim(adjustl(sparameter)),&
                  trim(adjustl(svalue)), isubstring=isubstring)
            end if
          end if
        else
          call output_line('Invalid option: '//trim(adjustl(soption))//'!',&
              OU_CLASS_WARNING,OU_MODE_STD,'flagship')
        end if
      end select
      
      ! Proceed with next command line argument
      iarg=iarg+1
    end do cmdarg
        
  end subroutine flagship_parseCmdlArguments

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
