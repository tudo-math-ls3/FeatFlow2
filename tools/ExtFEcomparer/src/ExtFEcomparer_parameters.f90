module ExtFEcomparer_parameters

! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  use fparser
  use paramlist
  use collection

  use triangulation
  use meshgeneration

  use element
  use cubature
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use bilinearformevaluation

  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop

  use vectorio
  use ExtFEcomparer_typedefs
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
    character(LEN=SYS_STRLEN) :: smaster

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

    integer :: i,ilvmin,ilvmax, ielemtype

    ! Variable for a filename:
    character(LEN=SYS_STRLEN) :: sString
    character(LEN=SYS_STRLEN) :: sPRMFile, sTRIFile
    character(LEN=SYS_STRLEN) :: sVectorFile

    ! Local variable
    character(LEN=SYS_STRLEN) :: sVectorName


    ! Get min/max level from the parameter file.
    !
    ! ilvmin receives the minimal level - the level of the
    ! mesh we created.
    ! ilvmax receives the level where our solution lives on

    call parlst_getvalue_int (rproblem%rparamlist,section,&
                              "NLMIN",ilvmin,2)
    call parlst_getvalue_int (rproblem%rparamlist,section,&
                              "NLMAX",ilvmax,4)

    call parlst_getvalue_int(rproblem%rparamlist,section,"ielementType",ielemtype,3)

    rproblem%NLMAX=ilvmax
    rproblem%NLMIN=ilvmin
    rproblem%ielemtype=ielemtype


    ! Get the .prm and the .tri file from the parameter list.
    call parlst_getvalue_string (rproblem%rparamList,section,&
                                 "sParametrisation",sPRMFile,bdequote=.true.)

    call parlst_getvalue_string (rproblem%rparamList,section,&
                                 "sMesh",sTRIFile,bdequote=.true.)

    call parlst_getvalue_string (rproblem%rparamList,section,&
                                 "sVector",sVectorFile,bdequote=.true.)


    !Save these information in the problem structure
    rproblem%sVectorFile = sVectorFile
    rproblem%sTRIFile = sTRIFile
    rproblem%sPRMFile = sPRMFile


    ! read out which cubature rule to use
    call parlst_getvalue_string (rproblem%rparamList,section,&
                                 "sCubFormula",sString,"AUTO_G3")
    rproblem%I_Cubature_Formula = cub_igetID(sString)

     ! Load the vector:
    call vecio_readBlockVectorHR(rproblem%coeffVector,sVectorName,.FALSE.,0,sVectorFile,.TRUE.)

end subroutine

subroutine ExtFEcomparer_doneParameters (rproblem)

!<description>
  ! Cleans up parameters read from the DAT files. Removes all references to
  ! parameters from the collection rproblem\%rcollection that were
  ! set up in cc_initParameters.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    ! Deallocate memory
    deallocate(rproblem%RlevelInfo)

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
    character(LEN=SYS_STRLEN) :: smaster
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
    call parlst_getvalue_string (rparlist,"LOGFILESETTINGS",&
        "smsgLog",slogfile,"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,"LOGFILESETTINGS",&
        "serrorLog",serrorfile,"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,"LOGFILESETTINGS",&
        "sbenchLog",sbenchlogfile,"",bdequote=.true.)

    ! That temporary parameter list is not needed anymore.
    call parlst_done (rparlist)

end subroutine

end module
