module ExtFEcomparer_parameters

! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  use paramlist

  use element

  use ExtFEcomparer_typedefs

  implicit none

  private

  public :: ExtFEcomparer_get_parameters
  public :: ExtFEcomparer_init_parameters
  public :: ExtFEcomparer_getLogFiles
  public :: ExtFEcomparer_parseCmdlArguments
  public :: ExtFEcomparer_getCmdlMasterDat

contains


subroutine ExtFEcomparer_get_parameters(rparamList)

!<description>
  ! Reads in all DAT files into the parameter list rparamlist
!</description>

!<inputoutput>
  ! The parameter list where the values of the DAT files should be stored.
  ! The structure must have been initialised, the parameters are just added
  ! to the list.
  type(t_parlist), intent(inout) :: rparamList
!</inputoutput>

!</subroutine>

    logical :: bexists
    character(LEN=ExtFE_STRLEN) :: smaster

    ! Check if a command line parameter specifies the master.dat file.
    ! It always is the last argument
    smaster = ''
    call ExtFEcomparer_getCmdlMasterDat(smaster)

    ! Let's see
    if(smaster .eq. '') then
        smaster = './data/master.dat'
    end if

    ! Read the file "master.dat".
    ! If that does not exist, try to manually read files with parameters from a
    ! couple of files.
    inquire(file=smaster, exist=bexists)

    if (bexists) then
      ! Read the master file. That either one contains all parameters or
      ! contains references to subfiles with data.
      call parlst_readfromfile (rparamList, smaster)
    else
      ! Each "readfromfile" command adds the parameter of the specified file
      ! to the parameter list.
      call parlst_readfromfile (rparamList, "./data/generalsettings.dat")
      call parlst_readfromfile (rparamList, "./data/function_specific_settings.dat")
      call parlst_readfromfile (rparamList, "./data/postprocessing.dat")
      call parlst_readfromfile (rparamList, "./data/logfile_settings.dat")

    end if

    ! Now we parse the command-line parameters in the parameter list
    call ExtFEcomparer_parseCmdlArguments(rparamList)

    ! Thats it

end subroutine


subroutine ExtFEcomparer_init_parameters(rproblem, section)
!<description>
  ! Initialises the structure rproblem with data from the initialisation
  ! files.
!</description>

!<input>
  ! A string to define the section
  ! Must be "ExtFE-FIRST" or "ExtFE-SECOND" to specify
  ! which function we init.
  character(LEN=*), intent(in) :: section
!<input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

    integer :: i,ilvmax, tmpvalue, ielemPairId
    integer(I32) :: ielemtype

    ! Variable for a filename:
    character(LEN=ExtFE_STRLEN) :: sString
    character(LEN=ExtFE_STRLEN) :: sPRMFile, sTRIFile
    character(LEN=ExtFE_STRLEN) :: sVectorFile
    character(LEN=ExtFE_STRLEN) :: smessage

    ! ilvmax receives the level where our solution lives on

    call parlst_getvalue_int (rproblem%rparamlist,section,&
                              "NLMAX",ilvmax)
    rproblem%NLMAX=ilvmax

    ! find out the dimension of the problem
    call parlst_getvalue_int(rproblem%rparamlist, "ExtFE-DOMAININFO", &
                            "dim", tmpvalue)
    rproblem%iDimension = tmpvalue

    ! Get the path of the .prm and the .tri file from the parameter list.
    ! .prm is not there in the 1D-Case
    if (rproblem%iDimension .ne. ExtFE_NDIM1) then
        call parlst_getvalue_string (rproblem%rparamList,section,&
                                    "sParametrisation",sPRMFile,bdequote=.true.)
        rproblem%sPRMFile = sPRMFile
    end if

    ! Get the path of the mesh
    call parlst_getvalue_string (rproblem%rparamList,section,&
                                 "sMesh",sTRIFile,bdequote=.true.)
    rproblem%sTRIFile = sTRIFile

    ! Get the path of the vector file
    call parlst_getvalue_string (rproblem%rparamList,section,&
                                 "sVector",sVectorFile,bdequote=.true.)
    rproblem%sVectorFile = sVectorFile

    ! Get the format identifier telling us what we have
    ! in the vector file
    call parlst_getvalue_int(rproblem%rparamlist,section, &
                            "sFileFormat",tmpvalue)
    rproblem%vectorFileFormat = tmpvalue
    ! Written out by write_blockVectorHR or something else?
    call parlst_getvalue_int(rproblem%rparamlist,section, &
                            "sVectorType", tmpvalue)
    rproblem%vectorType = tmpvalue
    ! How many variables are in the vector?
    call parlst_getvalue_int(rproblem%rparamlist,section, &
                            "iNVAR", tmpvalue)
    rproblem%NVAR = tmpvalue

    ! get the element type. Standard is -1 so that it will crash when
    ! we try to recreate the discretisation and nothing is specified.
    ! however, there is one more option to specify the element - a string!
    ! Get the element setting: One element or an ID of a pair?

    call parlst_getvalue_int(rproblem%rparamlist,section, &
                            "iElementSetting",tmpvalue)
    rproblem%elementSetting = tmpvalue


    select case(rproblem%elementSetting)
        case(ExtFE_ElementPair)
            call output_line("Searching for a numerical identifier of the element pair")
            call parlst_getvalue_int(rproblem%rparamlist,section,"ielementType",&
                        ielemPairId)
            call output_line("Found a numerical identifier of the element pair")
            rproblem%iElemPairID = ielemPairId

        case(ExtFE_OneElement)
            call output_line("Searching for the name of an element")
            call parlst_getvalue_string(rproblem%rparamlist,section, &
                    "ElementName",sString,bdequote=.true.)
            ielemtype = elem_igetID(sString)
            call output_line("Found a name of an element")
            rproblem%ielemType = ielemtype

        case(ExtFE_ElementList)
            call output_line("Search for an element list")

            ! How long is the list we find?
            tmpvalue = parlst_querysubstrings(rproblem%rparamlist,&
                                    section,"ElementList")
            if(tmpvalue .ne. rproblem%NVAR) then
                write(smessage,*) 'Input error: List length &
                &does not match the number of variables for function ', section
                call output_line(smessage , &
                OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_parameters")
                call sys_halt()
            end if

            ! Allocate the array for the element list with
            ! one element per variable
            allocate(rproblem%iElemList(rproblem%NVAR))
            do i=1,rproblem%NVAR
                call parlst_getvalue_string(rproblem%rparamlist, &
                    section, "ElementList", &
                            sString,sdefault="",isubstring=i,bdequote=.TRUE.)
                rproblem%iElemList(i) = elem_igetID(sString)
            end do
            call output_line("Found an element list")

        case default
            write(smessage,*) "Input error: choice not allowed for &
            &the selection of iElementSetting for function ", section
            call output_line(smessage , &
            OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_parameters")
            call sys_halt()
    end select

end subroutine

!<subroutine>
subroutine ExtFEcomparer_getLogFiles(slogfile,serrorfile,sbenchlogfile)

    !<description>
    ! Temporarily reads the output DAT file to get the names of the output
    ! files.
    !</description>

    !<output>
    ! Name of the message log file.
    character(LEN=*), intent(out) :: slogfile

    ! Name of the error log file.
    character(LEN=*), intent(out) :: serrorfile

    ! Name of the benchmark log file.
    character(LEN=*), intent(out) :: sbenchlogfile
    !</output>

    type(t_parlist) :: rparlist
    character(LEN=ExtFE_STRLEN) :: smaster
    logical :: bexists

    ! Init parameter list that accepts parameters for output files
    call parlst_init (rparlist)

    smaster = ''
    ! Check if a command line parameter specifies the master.dat file.
    call ExtFEcomparer_getCmdlMasterDat(smaster)
    if(smaster .eq. '') then
        smaster = './data/master.dat'
    else
        call output_line("Using the following master.dat: "//smaster)
    end if

    ! Read parameters that configure the output
    inquire(file=smaster, exist=bexists)

    if (bexists) then
        ! Read the master file. That either one contains all parameters or
        ! contains references to subfiles with data.
        call parlst_readfromfile (rparlist, smaster)
    else
        call parlst_readfromfile (rparlist, "./data/logfile_settings.dat")
    end if

    ! Now the real initialisation of the output including log file stuff!
    call parlst_getvalue_string (rparlist,"ExtFE-LOGFILESETTINGS",&
        "smsgLog",slogfile,"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,"ExtFE-LOGFILESETTINGS",&
        "serrorLog",serrorfile,"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,"ExtFE-LOGFILESETTINGS",&
        "sbenchLog",sbenchlogfile,"",bdequote=.true.)

    ! That temporary parameter list is not needed anymore.
    call parlst_done (rparlist)

end subroutine


!************************************************************************

!<subroutine>

  subroutine ExtFEcomparer_parseCmdlArguments(rparlist)

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

    iarg = 1; narg = sys_ncommandLineArgs()

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
        elseif (soption(1:7) .eq. 'smaster') then
            ! This is a master.dat We search for it in another routine
            ! don't do anything
        else
          call output_line('Invalid option: '//trim(adjustl(soption))//'!',&
              OU_CLASS_WARNING,OU_MODE_STD,'ExtFEcomparer_parseCmdLineArg')
        end if
      end select

      ! Proceed with next command line argument
      iarg=iarg+1
    end do cmdarg

  end subroutine


  subroutine ExtFEcomparer_getCmdlMasterDat(masterDat)

!<description>
    ! This subroutine parses the commandline arguments and modifies the
    ! parameter values in the global parameter list.
!</description>

!<inputoutput>
    ! path of the master.dat
    character(LEN=ExtFE_STRLEN),intent(inout) :: masterDat
!</inputoutput>
!</subroutine>

    ! local variables
    character(LEN=SYS_STRLEN) :: soption
    character(LEN=PARLST_MLDATA) :: svalue
    integer :: iarg,narg,iformat
    integer :: ipathlen
    iarg = 1; narg = sys_ncommandLineArgs()

    cmdarg: do while(iarg .le. narg)
      ! Retrieve next command line argument
      call sys_getcommandLineArg(iarg,soption,svalue,iformat)

      select case(iformat)
      case (0)
        ! Options without parameter values

      case (1,2)
        ! Options with given parameter values

        ! What option are we?
        if (soption(1:7) .eq. 'smaster') then
            ! no need to tokenize the string. soption=smaster:somepath
            ! We just get everything after the smaster:
            call output_line("Found a master.dat on the command line")
            !ipathlen = len(adjustl(trim(soption(9:len(soption)))))
            masterDat = svalue
            masterDat = trim(masterDat)
            return
        end if
      end select

      ! Proceed with next command line argument
      iarg=iarg+1

    end do cmdarg

  end subroutine


end module
