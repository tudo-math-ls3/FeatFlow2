module ExtFEcomparer_parameters

! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  use paramlist

  use element

  use ExtFEcomparer_typedefs

  implicit none


contains


subroutine ExtFEcomparer_get_parameters(rparamList)

!<description>
  ! Reads in all DAT files into the parameter list rparamlist
  ! This is based on the cc2d-code
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
    call sys_getcommandLineArg(1,smaster,sdefault="./data/master.dat")

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



end subroutine


subroutine ExtFEcomparer_init_parameters(rproblem, section)

!<description>
  ! Initialises the structure rproblem with data from the initialisation
  ! files.
  !
  ! The parameters in rproblem\%rparameters are evaluated.
  ! Important parameters are written to the problem structure
  ! rproblem. The parametrisation is read in, generated and
  ! stored in rproblem.
  ! The code is based on the cc2d-code
!</description>


!<input>
  ! A string to define the section
  character(LEN=*), intent(in) :: section
!<input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    integer :: i,ilvmin,ilvmax, ielemtype, tmpvalue

    ! Variable for a filename:
    character(LEN=ExtFE_STRLEN) :: sString
    character(LEN=ExtFE_STRLEN) :: sPRMFile, sTRIFile
    character(LEN=ExtFE_STRLEN) :: sVectorFile
    character(LEN=ExtFE_STRLEN) :: smessage

    ! ilvmax receives the level where our solution lives on

    call parlst_getvalue_int (rproblem%rparamlist,section,&
                              "NLMAX",ilvmax)
    rproblem%NLMAX=ilvmax


    call parlst_getvalue_int(rproblem%rparamlist,section, &
                            "iElementSetting",tmpvalue)
    rproblem%elementSetting = tmpvalue

    ! get the element type. Standard is -1 so that it will crash when
    ! we try to recreate the discretisation and nothing is specified.
    ! however, there is one more option to specify the element - a string!


    if (rproblem%elementSetting .eq. ExtFE_ElementPair) then
        call output_line("Searching for a numerical identifier of the element pair")
        call parlst_getvalue_int(rproblem%rparamlist,section,"ielementType",&
                    ielemtype)
        call output_line("Found a numerical identifier of the element pair")
    else if (rproblem%elementSetting .eq. ExtFE_OneElement) then
        call output_line("Searching for the name of an element")
        call parlst_getvalue_string(rproblem%rparamlist,section, &
                "ElementName",sString,bdequote=.true.)
        ielemtype = elem_igetID(sString)
        call output_line("Found a name of an element")
    else
        write(smessage,*) 'Input error: choice not allowed for &
        &the selection of iElementSetting for function ', section
        call output_line(smessage , &
        OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_parameters")
        call sys_halt()
    end if

    rproblem%ielemtype=ielemtype


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

    call parlst_getvalue_string (rproblem%rparamList,section,&
                                 "sMesh",sTRIFile,bdequote=.true.)
    rproblem%sTRIFile = sTRIFile

    ! Get the path of the vector file
    call parlst_getvalue_string (rproblem%rparamList,section,&
                                 "sVector",sVectorFile,bdequote=.true.)
    rproblem%sVectorFile = sVectorFile

    call parlst_getvalue_int(rproblem%rparamlist,section, &
                            "sFileFormat",tmpvalue)
    rproblem%vectorFileFormat = tmpvalue

    call parlst_getvalue_int(rproblem%rparamlist,section, &
                            "sVectorType", tmpvalue)
    rproblem%vectorType = tmpvalue

    call parlst_getvalue_int(rproblem%rparamlist,section, &
                            "iNVAR", tmpvalue)
    rproblem%NVAR = tmpvalue


end subroutine


!<subroutine>

subroutine ExtFEcomparer_getLogFiles (slogfile,serrorfile,sbenchlogfile)

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

    !</subroutine>

    type(t_parlist) :: rparlist
    character(LEN=ExtFE_STRLEN) :: smaster
    logical :: bexists

    ! Init parameter list that accepts parameters for output files
    call parlst_init (rparlist)

    ! Check if a command line parameter specifies the master.dat file.
    call sys_getcommandLineArg(1,smaster,sdefault="./data/master.dat")

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

end module
