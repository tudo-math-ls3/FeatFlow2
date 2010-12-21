! Contains the command handlers that encapsule ann FEAT commands.

module featcommandhandler

  use fsystem
  use genoutput
  use storage
  use collection
  use io
  use basicgeometry
  use boundary
  use triangulation
  
  use meshhierarchy
  use element
  use cubature

  use fespacehierarchybase
  use fespacehierarchy
  use spatialdiscretisation
  use multilevelprojection  
  use linearsystemscalar
  use linearsystemblock
  
  use multilevelprojection
  
  use derivatives
  use scalarpde
  use domainintegration
  use bilinearformevaluation
  use stdoperators
  use feevaluation
  use analyticprojection
  use spdiscprojection
  
  use vectorio
  use ucd
  
  use typedsymbol
  use commandparserbase
  use stdinoutparser

  implicit none
  
  private
  
  ! Execution mode: Standard execution
  integer, parameter, public :: FCMD_EXECSTD = 0

  ! Execution mode: Print short help message
  integer, parameter, public :: FCMD_EXECSHORTHELP = 1

  ! Execution mode: Print long help message
  integer, parameter, public :: FCMD_EXECLONGHELP = 2
  
  public :: fcmd_evalCommand
  
  ! Specification for a variable for the automatic variable check.
  type t_varspec
  
    ! Tag for the variable. May be ="", then the variable is an untagged, mandatory
    ! variable.
    character(len=SYS_NAMELEN) :: svartag
    
    ! Mandatory type for the argument. An STYPE_xxxx constant
    integer :: ctype = STYPE_INVALID
    
    ! Mandatory type for the variable in the collection. If not set to COLLCT_UNDEFINED,
    ! the variable must be part of a collection and have this type.
    integer :: ccollectionType = COLLCT_UNDEFINED
   
  end type
  
  
contains
  
  ! ***************************************************************************

  !<subroutine>
  
  subroutine fcmd_evalCommand (sstring,rcmdStatus,inestlevel,rvalue,bunknown,Rvalues)
  
  !<description>
    ! Evaluates a FEAT command.
  !</description>

  !<input>
    ! Name of the function
    character(len=*), intent(in) :: sstring
  
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus

    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! OPTIONAL: list of values or NULL if the list is empty.
    type(t_symbolValue), dimension(:), optional :: Rvalues
  !</input>
  
  !<output>
    ! Returns TRUE if the command is unknown.
    logical, intent(out) :: bunknown
  
    ! Return value of the function.
    type(t_symbolValue), intent(out) :: rvalue
  !</output>
  
  !</subroutine>
  
    bunknown = .false.
    rvalue%ctype = STYPE_INVALID
    
    if (sstring .eq. "help") then

      ! Call the FEAT command handler
      call fcmd_help (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "halt") then

      ! Stop the program.
      call fcmd_halt (rcmdStatus,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "meminfo") then

      ! Memory information
      call fcmd_meminfo (rcmdStatus,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "delete") then

      ! Delte a variabnle
      call fcmd_delete (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "destroy") then

      ! Destroy a variable
      call fcmd_destroy (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "show") then

      ! Information about a variable.
      call fcmd_show (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "info") then

      ! Information about a variable.
      call fcmd_info (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "printf") then

      ! Print text.
      call fcmd_printf (rcmdStatus,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "sprintf") then

      ! Put text to string.
      call fcmd_sprintf (rcmdStatus,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "read2dprm") then
      
      ! Call the FEAT command handler
      call fcmd_read2dprm (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "read2dtri") then
      
      ! Call the FEAT command handler
      call fcmd_read2dtri (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "meshrefine") then
      
      ! Call the FEAT command handler
      call fcmd_meshrefine (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "meshhierarchy") then
      
      ! Call the FEAT command handler
      call fcmd_meshhierarchy (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "fespace") then
      
      ! Call the FEAT command handler
      call fcmd_fespace (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "fehierarchy") then
      
      ! Call the FEAT command handler
      call fcmd_fehierarchy (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "mlevelprjhierarchy") then
      
      ! Call the FEAT command handler
      call fcmd_mlevelprjhierarchy (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "createblockvector") then
      
      ! Call the FEAT command handler
      call fcmd_createblockvector (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "readblockvector") then
      
      ! Call the FEAT command handler
      call fcmd_readblockvector (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "writeblockvector") then
      
      ! Call the FEAT command handler
      call fcmd_writeblockvector (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "copyvector") then
      
      ! Call the FEAT command handler
      call fcmd_copyvector (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "interpolatevector") then
      
      ! Call the FEAT command handler
      call fcmd_interpolatevector (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "extractfespacefromhier") then
      
      ! Call the FEAT command handler
      call fcmd_extractfespacefromhier (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "extractmeshfromhier") then
      
      ! Call the FEAT command handler
      call fcmd_extractmeshfromhier (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "extractsubvector") then
      
      ! Call the FEAT command handler
      call fcmd_extractsubvector (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "l2projection") then
      
      ! Call the FEAT command handler
      call fcmd_l2projection (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "writeucd") then
      
      ! Call the FEAT command handler
      call fcmd_writeucd (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "daxpyvector") then
      
      ! Call the FEAT command handler
      call fcmd_daxpyvector (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "clearvector") then
      
      ! Call the FEAT command handler
      call fcmd_clearvector (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else if (sstring .eq. "scalevector") then
      
      ! Call the FEAT command handler
      call fcmd_scalevector (rcmdStatus,inestlevel,rvalue,FCMD_EXECSTD,Rvalues)

    else
    
      ! Unknown command
      bunknown = .true.

    end if

  end subroutine

  ! ***************************************************************************

  !<subroutine>
  
  recursive subroutine fcmd_help (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: HELP.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>
  
  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! local variables
    character(len=COLLCT_MLNAME) :: scommand
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  help               - This help page.")
      
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("help - Shows help to a command.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("  help (""[command]"")")
      call output_lbrk ()
      call output_line ("Shows help for the command ""[command]"".")
    
      return
    end select

    if (.not. present(Rvalues)) then
      ! Short help
      call output_line ("FEAT2 command parser.")
      call output_line ("---------------------")
      call output_line ("The following commands are available. For a specific help")
      call output_line ("type 'help (""[command]"")'.")
      call output_lbrk ()
      
      ! Show the short help for all commands.
      call fcmd_help (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_internal_run (rcmdStatus,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_halt (rcmdStatus,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_meminfo (rcmdStatus,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_delete (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_destroy (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_show (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_info (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_printf (rcmdStatus,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_sprintf (rcmdStatus,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_read2dprm (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_read2dtri (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_meshrefine (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_meshhierarchy (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_fespace (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_fehierarchy (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_mlevelprjhierarchy (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_createblockvector (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_readblockvector (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_writeblockvector (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_copyvector (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_extractfespacefromhier (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_extractmeshfromhier (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_extractsubvector (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_l2projection (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_writeucd (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_daxpyvector (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_clearvector (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)
      call fcmd_scalevector (rcmdStatus,inestlevel,rreturn,FCMD_EXECSHORTHELP)

    else
      ! Get the command.
      call cmdprs_dequoteStd(Rvalues(1)%svalue,scommand)
    
      if (scommand .eq. "help") then
        call fcmd_help (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)
        
      else if (scommand .eq. "run") then
        call fcmd_internal_run (rcmdStatus,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "exit") then
        call fcmd_halt (rcmdStatus,rreturn,FCMD_EXECLONGHELP)
        
      else if (scommand .eq. "meminfo") then
        call fcmd_meminfo (rcmdStatus,rreturn,FCMD_EXECLONGHELP)
        
      else if (scommand .eq. "delete") then
        call fcmd_delete (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)
        
      else if (scommand .eq. "destroy") then
        call fcmd_destroy (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)
        
      else if (scommand .eq. "show") then
        call fcmd_show (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)
        
      else if (scommand .eq. "info") then
        call fcmd_info (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)
        
      else if (scommand .eq. "halt") then
        call fcmd_halt (rcmdStatus,rreturn,FCMD_EXECLONGHELP)
        
      else if (scommand .eq. "printf") then
        call fcmd_printf (rcmdStatus,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "sprintf") then
        call fcmd_sprintf (rcmdStatus,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "read2dprm") then
        call fcmd_read2dprm (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)
        
      else if (scommand .eq. "read2dtri") then
        call fcmd_read2dtri (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)
        
      else if (scommand .eq. "meshrefine") then
        call fcmd_meshrefine (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)
        
      else if (scommand .eq. "meshhierarchy") then
        call fcmd_meshhierarchy (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "fespace") then
        call fcmd_fespace (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "fehierarchy") then
        call fcmd_fehierarchy (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "mlevelprjhierarchy") then
        call fcmd_mlevelprjhierarchy (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "createblockvector") then
        call fcmd_createblockvector (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "readblockvector") then
        call fcmd_readblockvector (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "writeblockvector") then
        call fcmd_writeblockvector (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "copyvector") then
        call fcmd_copyvector (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "interpolatevector") then
        call fcmd_interpolatevector (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "extractfespacefromhier") then
        call fcmd_extractfespacefromhier (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)
      
      else if (scommand .eq. "extractmeshfromhier") then
        call fcmd_extractmeshfromhier (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)
      
      else if (scommand .eq. "extractsubvector") then
        call fcmd_extractsubvector (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "l2projection") then
        call fcmd_l2projection (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "writeucd") then
        call fcmd_writeucd (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "daxpyvector") then
        call fcmd_daxpyvector (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "clearvector") then
        call fcmd_clearvector (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)

      else if (scommand .eq. "scalevector") then
        call fcmd_scalevector (rcmdStatus,inestlevel,rreturn,FCMD_EXECLONGHELP)
      else
        call output_line ("help: Unknown command.")
      end if
    end if

    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************
  ! Auzxiliary routines
  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_getparameters (Rvarspec,Iindex,bsuccess,rcollection,inestlevel,Rvalues)
  
  !<description>
    ! Checks parameters and finds indices of optional parameters.
  !</description>
  
  !<input>
    ! Variable specification. Contains name and mandatory types vor the
    ! parameters in Rvalues. Parameters without a tag are mandatory
    ! and their number has to match the non-anynomous variables in Rvalues.
    ! Tagged variables are optional and searched for in Rvalues.
    ! For every variable in Snames, Itypes must specify a mandatory type,
    ! Iindex will receive the index of the named (optional) variables
    ! in Rvalues.
    type(t_varspec), dimension(:), intent(in) :: Rvarspec
    
    ! Collection object that contains all variables of the program.
    type(t_collection), intent(in) :: rcollection
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>
  
  !<output>
    ! For each named variable in Snames, this array receives an index
    ! where this variable can be found in Rvalues. The index is =0
    ! if the variable does not appear in Rvalues.
    ! The unnamed variables receive indices 1,2,3,...
    integer, dimension(:), intent(out) :: Iindex

    ! Returns TRUE if the parameters match in number and type the definitions
    ! in Snames/Itypes.
    logical, intent(out) :: bsuccess
  !</output>
  
  !</subroutine>
  
    ! local variables
    integer :: i,j,nuntagged
    integer :: ctype
    character(len=SYS_NAMELEN) :: ssection
    logical :: bexists
    
    bsuccess = .false.
    Iindex(:) = 0
    
    ! Figure out the number of untagged variables
    do nuntagged = 1,size(Rvarspec)
      if (Rvarspec(nuntagged)%svartag .ne. "") exit
    end do
    nuntagged = nuntagged-1
    
    ! This must match the number of arguments.
    if (nuntagged .gt. 0) then
      if (.not. present (Rvalues)) return
      if (ubound(Rvalues,1) .lt. nuntagged) then
        call output_line ("Not enough arguments.")
        return
      end if
      do i=1,nuntagged
        if (Rvalues(i)%svartag .ne. "") then
          call output_line ("Not enough arguments.")
          return
        end if
      end do
    else
      ! If there are no arguments and all are optional, that's ok.
      if (.not. present (Rvalues)) then
        bsuccess = .true.
        return
      end if
    end if
    
    ! Fill the indices of the unnamed variables. Stop if a type does not match.
    do i=1,nuntagged
      if (Rvarspec(i)%ctype .ne. STYPE_UNDEFINED) then
        if (Rvalues(i)%ctype .ne. Rvarspec(i)%ctype) then
          call output_line ("Invalid arguments.")
          return
        end if
      end if
      Iindex(i) = i
    end do
    
    ! Now find the optional parameters.
    do i=nuntagged+1,size(Rvarspec)
      do j=nuntagged+1,size(Rvalues)
        if (Rvarspec(i)%svartag .eq. Rvalues(j)%svartag) then
          if (Rvarspec(i)%ctype .ne. STYPE_UNDEFINED) then
            ! Cancel if the type is wrong.
            if (Rvarspec(i)%ctype .ne. Rvalues(j)%ctype) then
              call output_line ("Invalid arguments.")
              return
            end if
          end if
          
          ! Remember the index, parameter done.
          Iindex(i) = j
          exit
        end if
      end do
    end do
    
    ! Check the parameters which have to fulfil the type in the collection
    ! to be correct.
    do i=1,size(Rvarspec)
      
      if (Rvarspec(i)%ccollectiontype .ne. COLLCT_UNDEFINED) then
        
        ! Check if we have this parameter
        if (Iindex(i) .ne. 0) then
        
          ! Is this a parameter from the collection?
          if (Rvalues(Iindex(i))%svarname .eq. "") then
            call output_line ("Unknown variable")
            return
          end if
          
          call cmdprs_getSymbolSection (rcollection,Rvalues(Iindex(i))%svarname,&
              inestlevel,ssection,bexists)

          if (.not. bexists) then
            call output_line ("Unknown variable: "//trim(Rvalues(Iindex(i))%svarname))
            return
          end if
        
          ! Is this parameter present in the collection?
          ctype = collct_gettype (rcollection, Rvalues(Iindex(i))%svarname, ssectionName=ssection)
          if (ctype .ne. Rvarspec(i)%ccollectionType) then
            call output_line ("Invalid type: "//trim(Rvalues(Iindex(i))%svarname))
            return
          end if
        
        end if
        
      end if
      
    end do
    
    bsuccess = .true.
    
  end subroutine

  ! ***************************************************************************
  ! internal subroutines, only help messages
  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_internal_run (rcmdStatus,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: RUN.
    ! Internal command, this routine only provides help text.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! local variables
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  run                - Execute script.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("run() - Execute a script on hard disc.")
      call output_lbrk ()
      call output_line ("Usage on the command line:")
      call output_line ("  run ""[filename/path]""")
      call output_lbrk ()
      call output_line ("Usage in a script:")
      call output_line ("  run (""[filename/path]"")")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_internal_print (rcmdStatus,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: RUN.
    ! Internal command, this routine only provides help text.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! local variables
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  print              - Print some text.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("print - prints some text to the terminal.")
      call output_lbrk ()
      call output_line ("Usage on command line:")
      call output_line ("  PRINT ""[some text]""")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************
  ! General FEAT commands and expressions
  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_halt (rcmdStatus,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: HALT.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! local variables
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  exit               - Exits the command line/a running loop.")
      call output_line ("  halt()             - Exits the command line/Stop the program.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("halt - Stops a running program or close interactive command line.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("  halt()")
      call output_lbrk ()
      call output_lbrk ()
      call output_line ("exit - Exit a loop in a program.")
      call output_lbrk ()
      call output_line ("Usage in a script:")
      call output_line ("  exit")
      call output_lbrk ()
      call output_lbrk ()
      call output_line ("exit - Close interactive command line, return to terminal.")
      call output_lbrk ()
      call output_line ("Usage on command line:")
      call output_line ("  exit")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    rcmdStatus%bterminate = .true.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_meminfo (rcmdStatus,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: HALT.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! local variables
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  meminfo()          - Information about memory management.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("meminfo - Information about memory management.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("  meminfo()")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    call storage_info (.true.)
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_show (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: INFO.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(1), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_UNDEFINED, COLLCT_UNDEFINED) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    character(len=SYS_NAMELEN) :: svalname,ssection
    character(len=SYS_STRLEN) :: svalue
    integer :: ctype
    logical :: bexists
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  show()             - Show environment variable (if possible).")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("show - Show content about an environment vairbale.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("  show ()")
      call output_line ("  show ([name])")
      call output_lbrk ()
      call output_line ("If possible, this routine prints information about an")
      call output_line ("environment variable to the terminal. For standard types")
      call output_line ("(int, double, string,...), the value is returned.")
      call output_line ("For extended types (e.g. meshes), status information")
      call output_line ("is returned.")
      call output_line ("If no variable is specified, the environment information")
      call output_line ("for all variables is shown.")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    if (.not. present(Rvalues)) then
    
      ! Environment statistics.
      call output_lbrk ()
      call output_line ("Environment statistics:")
      call output_lbrk ()
      ! Print the current status of the collection.
      call collct_printStatistics (rcmdStatus%rcollection)
    
    else
  
      ! Check parameters
      call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
          rcmdStatus%rcollection,inestlevel,Rvalues)
      if (.not. bsuccess) return

      ! Get the identifiers
      svalname = Rvalues(1)%svarname
      
      if (svalname .eq. "") then
        call output_line ("Unknown environment variable!")
        return
      end if

      ! That's an awful task. Figure out the type at first.
      ! Then get the value and/or print information.
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,svalname,inestlevel,ssection)
      ctype = collct_gettype (rcmdStatus%rcollection, svalname, bexists=bexists, &
          ssectionName=ssection)
      
      if (.not. bexists) then
        call output_line ("Unknown environment variable!")
      else
        select case (ctype)
        case (COLLCT_INTEGER)
          call output_line (trim(sys_siL(&
            collct_getvalue_int (rcmdStatus%rcollection, svalname) ,10)))
            
        case (COLLCT_REAL)
          call output_line (trim(sys_sdEL(&
            collct_getvalue_real (rcmdStatus%rcollection, svalname) ,10)))
            
        case (COLLCT_STRING)
          call collct_getvalue_string (rcmdStatus%rcollection, svalname, svalue)
          call output_line (svalue)
            
        case default
          call output_line ("No information available.")
        end select
        
      end if
     
    end if
    
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_info (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: INFO.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(1), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_UNDEFINED, COLLCT_UNDEFINED) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    character(len=COLLCT_MLNAME) :: ssection
    character(len=SYS_NAMELEN) :: svalname
    character(len=SYS_STRLEN) :: svalue
    integer :: ctype
    logical :: bexists
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_meshhierarchy), pointer :: p_rmeshhierarchy
    type(t_fespaceLevel), pointer :: p_rfespace
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_vectorBlock), pointer :: p_rvectorBlock
    type(t_vectorScalar), pointer :: p_rvectorScalar
    
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  info()             - Detailed information about an environment variable.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("info - Detailed information about environment variables.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("  info [name]")
      call output_lbrk ()
      call output_line ("Shows detailed information about an environment variable.")
      call output_line ("This includes type, value,...")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return

    ! Get the identifiers
    svalname = Rvalues(1)%svarname
    
    if (svalname .eq. "") then
      call output_line ("Unknown environment variable!")
      return
    end if

    ! That's an awful task. Figure out the type at first.
    ! Then get the value and/or print information.
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,svalname,inestlevel,ssection)
    ctype = collct_gettype (rcmdStatus%rcollection, svalname, bexists=bexists, &
        ssectionName=ssection)
    
    if (.not. bexists) then
      call output_line ("Unknown environment variable!")
    else
      call output_line ("Variable: "//svalname)
      
      select case (ctype)
      case (COLLCT_INTEGER)
        call output_line ("Type    : Integer")
        call output_line ("Content : ",bnolinebreak=.true., bnotrim=.true.)
        call output_line (trim(sys_siL(&
          collct_getvalue_int (rcmdStatus%rcollection, svalname, ssectionName=ssection) ,10)))
          
      case (COLLCT_REAL)
        call output_line ("Type    : Double precision")
        call output_line ("Content : ",bnolinebreak=.true., bnotrim=.true.)
        call output_line (trim(sys_sdEL(&
          collct_getvalue_real (rcmdStatus%rcollection, svalname, ssectionName=ssection) ,10)))
          
      case (COLLCT_STRING)
        call output_line ("Type    : string")
        call output_line ("Content : ",bnolinebreak=.true., bnotrim=.true.)
        call collct_getvalue_string (rcmdStatus%rcollection, svalname, svalue, ssectionName=ssection)
        call output_line (svalue)
          
      case (COLLCT_BOUNDARY)
        call output_line ("Type    : Boundary object (2D)")

      case (COLLCT_TRIA)
        call output_line ("Type    : Triangulation object")
        call output_lbrk ()
        p_rtriangulation => collct_getvalue_tria (rcmdStatus%rcollection, svalname, ssectionName=ssection)
        call tria_infoStatistics (p_rtriangulation,.true.)

      case (COLLCT_MSHHIERARCHY)
        call output_line ("Type    : Mesh hierarchy")
        call output_lbrk ()
        p_rmeshhierarchy => collct_getvalue_mshh (rcmdStatus%rcollection, svalname, ssectionName=ssection)
        call mshh_printHierStatistics (p_rmeshhierarchy)

      case (COLLCT_FESPACE)
        call output_line ("Type    : FE space")
        call output_lbrk ()
        p_rfeSpace => collct_getvalue_fesp (rcmdStatus%rcollection, svalname, ssectionName=ssection)
        call fesph_infoStatistics (p_rfeSpace,.true.)
        call spdiscr_infoBlockDiscr (p_rfeSpace%p_rdiscretisation)

      case (COLLCT_FEHIERARCHY)
        call output_line ("Type    : Mesh hierarchy")
        call output_lbrk ()
        p_rfeHierarchy => collct_getvalue_feh (rcmdStatus%rcollection, svalname, ssectionName=ssection)
        call fesph_printHierStatistics (p_rfeHierarchy)

      case (COLLCT_MLPRJHIERARCHY)
        call output_line ("Type    : Multilevel projection hierarchy")
        call output_lbrk ()

      case (COLLCT_BLKVECTOR)
        call output_line ("Type    : Block vector")
        call output_lbrk ()
        p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, svalname, ssectionName=ssection)
        call lsysbl_infoVector (p_rvectorBlock)

      case (COLLCT_SCAVECTOR)
        call output_line ("Type    : Scalar vector")
        call output_lbrk ()
        p_rvectorScalar => collct_getvalue_vecsca (rcmdStatus%rcollection, svalname, ssectionName=ssection)
        call lsyssc_infoVector (p_rvectorScalar)

      case default
        call output_line ("Type    : "//sys_siL(ctype,10))
      end select
      
    end if
    
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_destroy (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: destroy.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(1), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_UNDEFINED, COLLCT_UNDEFINED) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    character(len=COLLCT_MLNAME) :: ssection
    character(len=SYS_NAMELEN) :: svalname
    logical :: bexists
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  info()             - Detailed information about an environment variable.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("info - Detailed information about environment variables.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("  info [name]")
      call output_lbrk ()
      call output_line ("Shows detailed information about an environment variable.")
      call output_line ("This includes type, value,...")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return

    ! Get the identifiers
    svalname = Rvalues(1)%svarname
    
    ! Destroy the variable
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,svalname,inestlevel,ssection,bexists)
    if (bexists) then
      call do_destroy (rcmdStatus,svalname,ssection,.true.)
    end if
    
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_delete (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: delete.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(1), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_UNDEFINED, COLLCT_UNDEFINED) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    character(len=COLLCT_MLNAME) :: ssection
    character(len=SYS_NAMELEN) :: svalname
    logical :: bexists
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  delete()           - Delete variable from environment.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("delete - Delete a variable from the environment.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("  delete ([name])")
      call output_lbrk ()
      call output_line ("WARNING: This does not release probably associated memory!")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return

    ! Get the identifiers
    svalname = Rvalues(1)%svarname
    
    ! Destroy the variable
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,svalname,inestlevel,ssection,bexists)
    if (bexists) then
      call collct_deletevalue (rcmdStatus%rcollection, svalname,ssectionName=ssection)
    end if
    
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_printf (rcmdStatus,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: PRINTF.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! local variables
    character (len=SYS_STRLEN) :: stemp,stemp2
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  printf()           - Print text to the terminal.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("printf - Print text to the termininal.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("  printf(""[format string]"",...)")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    rreturn%ctype = STYPE_INTEGER

    if (present(Rvalues)) then
      if (Rvalues(1)%ctype .eq. STYPE_STRING) then
      
        ! Evaluate the arguments.
        stemp = Rvalues(1)%svalue
        call sioprs_qualifyString (stemp,Rvalues(2:))
        
        ! Print to terminal
        call cmdprs_dequoteStd(stemp,stemp2)
        call output_line (stemp2)
        
        ! Worked.
        return
      end if
      
    end if

    ! If we come to here, there is something wrong.
    call output_line ("Invalid arguments!")

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_sprintf (rcmdStatus,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: SPRINTF.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(inout), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! local variables
    character (len=SYS_STRLEN) :: stemp
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  sprintf()          - Print text to a string.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("sprintf - Print text to a string.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("  sprintf([destination string],""[format string]"",...)")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    rreturn%ctype = STYPE_INTEGER

    if (present(Rvalues)) then
      if (size(Rvalues) .ge. 2) then
        if ((Rvalues(1)%ctype .eq. STYPE_STRING) .and.&
            (Rvalues(2)%ctype .eq. STYPE_STRING)) then
        
          ! Evaluate the arguments.
          stemp = Rvalues(2)%svalue
          if (size(Rvalues) .ge. 3) then
            call sioprs_qualifyString (stemp,Rvalues(3:))
          else
            call sioprs_qualifyString (stemp)
          end if
          
          ! Print to string
          Rvalues(1)%svalue = trim(stemp)
          Rvalues(1)%ilength = len_trim(stemp)
          call tpsym_saveSymbol (Rvalues(1))
          
          ! Worked.
          return
        end if
      end if
      
    end if

    ! If we come to here, there is something wrong.
    call output_line ("Invalid arguments!")

  end subroutine

  ! ***************************************************************************

  subroutine do_destroy (rcmdStatus,svalname,svalsection,bverbose)
  
  !<description>
    ! Command handler: DESTROY.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Name of the variable to destroy
    character(len=*), intent(in) :: svalname

    ! Name of the section containing the variable
    character(len=*), intent(in) :: svalsection
    
    ! Print vebose messages
    logical, intent(in) :: bverbose
  !</input>
  
    ! local variables
    integer :: ctype,itag
    logical :: bexists,bdelete
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_meshhierarchy), pointer :: p_rmeshhierarchy
    type(t_feSpaceLevel), pointer :: p_rfeSpace
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_interlevelProjectionHier), pointer :: p_rinterlevelProjectionHier
    type(t_vectorBlock), pointer :: p_rvectorBlock
    type(t_vectorScalar), pointer :: p_rvectorScalar
  
    ! That's an awful task. Figure out the type at first.
    ! Then get the value and/or print information.
    ctype = collct_gettype (rcmdStatus%rcollection, svalname, bexists=bexists, ssectionName=svalsection)
    
    if (.not. bexists) then
      if (bverbose) then
        call output_line ("Unknown environment variable!")
      end if
    else
      ! Get the user defined tag. If the tag is set to 1, do not destroy the variable
      ! but only delete it -- it is associated to another structure!
      itag = collct_gettag(rcmdStatus%rcollection, svalname, ssectionName=svalsection)
      
      bdelete = .true.

      if (itag .eq. 0) then
        select case (ctype)
        case (COLLCT_INTEGER,COLLCT_REAL,COLLCT_STRING)
          ! Primitive type
            
        case (COLLCT_BOUNDARY)
          ! Boundary object
          p_rboundary => collct_getvalue_bdry (rcmdStatus%rcollection, svalname,ssectionName=svalsection)
          call boundary_release(p_rboundary)
          deallocate(p_rboundary)

        case (COLLCT_TRIA)
          ! Boundary object
          p_rtriangulation => collct_getvalue_tria (rcmdStatus%rcollection, svalname,ssectionName=svalsection)
          call tria_done(p_rtriangulation)
          deallocate(p_rtriangulation)

        case (COLLCT_MSHHIERARCHY)
          ! Boundary object
          p_rmeshhierarchy => collct_getvalue_mshh (rcmdStatus%rcollection, svalname,ssectionName=svalsection)
          call mshh_releaseHierarchy(p_rmeshhierarchy)
          deallocate(p_rmeshhierarchy)

        case (COLLCT_FESPACE)
          ! Boundary object
          p_rfeSpace => collct_getvalue_fesp (rcmdStatus%rcollection, svalname,ssectionName=svalsection)
          call fesph_releaseFEspace(p_rfeSpace)
          deallocate(p_rfeSpace)

        case (COLLCT_FEHIERARCHY)
          ! Boundary object
          p_rfeHierarchy => collct_getvalue_feh (rcmdStatus%rcollection, svalname,ssectionName=svalsection)
          call fesph_releaseHierarchy(p_rfeHierarchy)
          deallocate(p_rfeHierarchy)

        case (COLLCT_MLPRJHIERARCHY)
          ! Multilevel projection hierarchy
          p_rinterlevelProjectionHier => collct_getvalue_mlprjh (rcmdStatus%rcollection, svalname,&
              ssectionName=svalsection)
          call mlprj_releasePrjHierarchy(p_rinterlevelProjectionHier)
          deallocate(p_rinterlevelProjectionHier)

        case (COLLCT_BLKVECTOR)
          ! Multilevel projection hierarchy
          p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, svalname,&
              ssectionName=svalsection)
          call lsysbl_releaseVector(p_rvectorBlock)
          deallocate(p_rvectorBlock)

        case (COLLCT_SCAVECTOR)
          ! Multilevel projection hierarchy
          p_rvectorScalar => collct_getvalue_vecsca (rcmdStatus%rcollection, svalname,&
              ssectionName=svalsection)
          call lsyssc_releaseVector(p_rvectorScalar)
          deallocate(p_rvectorScalar)

        case default
          call output_line ("Unknown type, variable cannot be destroyed!")         
          bdelete = .false.
        end select
      end if

      if (bdelete) then
        call collct_deletevalue (rcmdStatus%rcollection, svalname, ssectionName=svalsection)
      end if
      
    end if
    
  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_read2dprm (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: READ2DPRM.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(2), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &    ! Variable name
         t_varspec("", STYPE_STRING, COLLCT_UNDEFINED)  &    ! File name
       /)
    integer, dimension(size(Rvarspec)) :: Iindex

    ! local variables
    logical :: bsuccess
    logical :: bexists
    character(len=SYS_STRLEN) :: sfilename
    character(len=COLLCT_MLNAME) :: sname, ssection
    type(t_boundary), pointer :: p_rboundary
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  read2dprm()        - Read 2D .prm file.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("read2dprm - Read 2D .PRM file.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("  read2dprm ([varname],[filename])")
      call output_lbrk ()
      call output_line ("Read in a .prm file and store it using the name [variable].")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return

    ! Get the identifier and file name 
    sname = Rvalues(1)%svarname
    call cmdprs_dequoteStd(Rvalues(2)%svalue,sfilename)
    
    if ((sname .eq. "") .or. (sfilename .eq. "")) then
      call output_line ("Invalid arguments.")
      return
    end if

    ! Open the file and read.
    inquire(file=trim(sfilename), exist=bexists)
    if (.not. bexists) then
      call output_line ("File not found!")
      return
    else
      call output_line ("Reading file: "//trim(sfilename))
    
      ! Read
      allocate (p_rboundary)
      call boundary_read_prm(p_rboundary, sfilename)
      
      ! Remove old value from collection if present
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,sname,inestlevel,ssection)
      call do_destroy(rcmdStatus,sname,ssection,.false.)
      
      ! Add to the collection
      call collct_setvalue_bdry (rcmdStatus%rcollection, sname, p_rboundary, .true., &
          ssectionName=ssection)
    end if

    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_read2dtri (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: READ2DTRI.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(3), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_STRING, COLLCT_UNDEFINED), &
         t_varspec("boundary", STYPE_VAR, COLLCT_BOUNDARY)  &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex

    ! local variables
    logical :: bsuccess
    logical :: bexists
    character(len=SYS_STRLEN) :: sfilename, stoken
    character(len=COLLCT_MLNAME) :: sname,ssection
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  read2dtri()        - Read 2D .tri file.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("read2dtri - Read 2D .TRI file.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("  read2dtri ([varname],[filename] [,boundary=[varbd]]")
      call output_lbrk ()
      call output_line ("Read in a .tri file and store it using the name [variable].")
      call output_line ("If BOUNDARY is specified, the triangulation is connected")
      call output_line ("to a boundary object varbd.")
      call output_lbrk ()
      call output_line ("Example:")
      call output_line ("  read2dprm (myprm,""myprmfile.prm"";");
      call output_line ("  read2dtri (mytri,""mytrifile.tri"",boundary=myprm);")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check the parameters
    if (.not. present(Rvalues)) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    if (size(Rvalues) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if

    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return
    
    ! Get the identifier and file name 
    sname = Rvalues(1)%svarname
    call cmdprs_dequoteStd(Rvalues(2)%svalue,sfilename)
    nullify(p_rboundary)
    
    if (Iindex(3) .ne. 0) then
      ! Get the bondary object.
      stoken = Rvalues(Iindex(3))%svarname
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
      p_rboundary => collct_getvalue_bdry(rcmdStatus%rcollection,stoken,bexists=bexists,&
          ssectionName=ssection)
    end if
    
    ! Open the file and read.
    inquire(file=trim(sfilename), exist=bexists)
    if (.not. bexists) then
      call output_line ("File not found!")
    else
    
      call output_line ("Reading file: "//trim(sfilename))
    
      ! Read
      allocate (p_rtriangulation)
      
      if (associated(p_rboundary)) then
        call tria_readTriFile2D(p_rtriangulation, sfilename, p_rboundary)
        
        ! Create a standard mesh.
        call tria_initStandardMeshFromRaw (p_rtriangulation, p_rboundary)
      else
        call output_line("Warning. No boundary specified!")
      
        call tria_readTriFile2D(p_rtriangulation, sfilename)
        
        ! Create a standard mesh.
        call tria_initStandardMeshFromRaw (p_rtriangulation)
      end if
      
      ! Remove old value from collection if present
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,sname,inestlevel,ssection)
      call do_destroy(rcmdStatus,sname,ssection,.false.)
      
      ! Add to the collection
      call collct_setvalue_tria (rcmdStatus%rcollection, sname, p_rtriangulation, .true.,&
          ssectionName=ssection)
    end if

    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_meshrefine (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: MESHREFINE.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(4), parameter :: Rvarspec = &
      (/ t_varspec("",         STYPE_VAR,     COLLCT_TRIA), &
         t_varspec("levels",   STYPE_INTEGER, COLLCT_UNDEFINED), &
         t_varspec("boundary", STYPE_VAR,     COLLCT_BOUNDARY), &
         t_varspec("method",   STYPE_STRING,  COLLCT_UNDEFINED) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    character(len=SYS_STRLEN) :: stoken
    integer :: cmethod, ilevels
    integer :: i
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    character(len=COLLCT_MLNAME) :: ssection
    logical :: bexists
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  meshrefine()       - Refine a mesh.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("meshrefine - Refine a mesh.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   meshrefine ([varmesh] [, ...options...]);")
      call output_lbrk ()
      call output_line ("Refine a mesh with a given method one or multiple times.")
      call output_line ("[varmesh] identifies the mesh to refine.")
      call output_line ("The following options are possible in [...options...]:")
      call output_lbrk ()
      call output_line ("  ... levels=[varname] ...")
      call output_line ("      Refine the mesh [varname] times. Default is ""levels=1"".")
      call output_lbrk ()
      call output_line ("  ... boundary=[varname] ...")
      call output_line ("      Use the specified boundary object [varname] fdor boundary refinement.")
      call output_line ("      If not specified, no boundary object is used.")
      call output_lbrk ()
      call output_line ("  ... method=[varname] ...")
      call output_line ("      Use a specific method. Default is ""method ""2levelordered"""". Possible choices:")
      call output_line ("      ""2levelordered""   - Use 2-level-ordering refinement.")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return

    ! Get the identifier and file name 
    nullify(p_rboundary)
    nullify(p_rtriangulation)
    ilevels = 1
    cmethod = 0 ! 2-level ordering
    
    if (Iindex(2) .ne. 0) then
      ! Refinement levels
      ilevels = Rvalues(Iindex(2))%ivalue
    end if
    
    if (Iindex(3) .ne. 0) then
      ! Get the bondary object.
      stoken = Rvalues(Iindex(3))%svarname
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
      p_rboundary => collct_getvalue_bdry(rcmdStatus%rcollection,stoken,bexists=bexists,&
          ssectionName=ssection)
    end if
    
    if (Iindex(4) .ne. 0) then
      ! Method
      if (Rvalues(Iindex(4))%svalue .eq. "2levelordered") then
        cmethod = 0
      else
        call output_line ("Warning: Unknown refinement method! Using default.")
        cmethod = 0
      end if
    end if

    ! Get the triangulation object.
    stoken = Rvalues(Iindex(1))%svarname
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
    p_rtriangulation => collct_getvalue_tria(rcmdStatus%rcollection,stoken,ssectionName=ssection)

    call output_line ("Refining mesh... [",bnolinebreak=.true.)
  
    ! Refine the mesh
    select case (cmethod)
    case (0)
      ! 2-level ordered
      
      do i = 1,ilevels
        call output_line (" "//trim(sys_siL(i,10)),bnolinebreak=.true.)

        ! Boundary present?
        if (associated(p_rboundary)) then
          call tria_quickRefine2LevelOrdering(1,p_rtriangulation)
        else
          call tria_quickRefine2LevelOrdering(1,p_rtriangulation,p_rboundary)
        end if
      end do

      ! Create a standard mesh.
      if (associated(p_rboundary)) then
        call tria_initStandardMeshFromRaw (p_rtriangulation)
      else
        call tria_initStandardMeshFromRaw (p_rtriangulation,p_rboundary)
      end if
      
      call output_line ("]")
    end select

    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_meshhierarchy (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: MESHIHIERARCHY.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(5), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("mesh", STYPE_VAR, COLLCT_TRIA), &
         t_varspec("boundary", STYPE_VAR, COLLCT_BOUNDARY), &
         t_varspec("levels", STYPE_INTEGER, COLLCT_UNDEFINED), &
         t_varspec("method", STYPE_INTEGER, COLLCT_UNDEFINED)  &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    character(len=SYS_STRLEN) :: stoken
    integer :: cmethod, ilevels
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    character(len=COLLCT_MLNAME) :: sname,ssection
    type(t_meshHierarchy), pointer :: p_rmeshHierarchy
    logical :: bexists
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  meshhierarchy()    - Create a mesh hierarchy.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("meshhierarchy - Create a mesh hierarchy.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   meshhierarchy ([varname] [, ...options...])")
      call output_lbrk ()
      call output_line ("Refine a mesh with a given method one or multiple times.")
      call output_line ("[varname] identifies the variable to create.")
      call output_line ("The following options are possible in [...options...]:")
      call output_lbrk ()
      call output_line ("  ... mesh=[varname] ...")
      call output_line ("      Use mesh [varname] as coarse mesh of teh hierarchy.")
      call output_line ("      Mandatory argument")
      call output_lbrk ()
      call output_line ("  ... boundary=[varname] ...")
      call output_line ("      Use the specified boundary object [varname] for boundary refinement.")
      call output_line ("      If not specified, no boundary object is used.")
      call output_lbrk ()
      call output_line ("  ... levels=[varname] ...")
      call output_line ("      Creates a hierarchy of [varname] levels. Default is ""levels=1"".")
      call output_lbrk ()
      call output_line ("  ... method=""[varname]"" ...")
      call output_line ("      Use a specific method. Default is ""method ""2levelordered"""". Possible choices:")
      call output_line ("      2levelordered   - Use 2-level-ordering refinement.")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return

    ! Get the identifier and file name 
    nullify(p_rboundary)
    nullify(p_rtriangulation)
    sname = Rvalues(1)%svarname
    ilevels = 1
    cmethod = 0 ! 2-level ordering
    
    if (Iindex(2) .ne. 0) then
      ! Get the triangulation object.
      stoken = Rvalues(Iindex(2))%svarname
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
      p_rtriangulation => collct_getvalue_tria(rcmdStatus%rcollection,stoken,ssectionName=ssection)
    end if

    if (Iindex(3) .ne. 0) then
      ! Get the bondary object.
      stoken = Rvalues(Iindex(3))%svarname
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
      p_rboundary => collct_getvalue_bdry(rcmdStatus%rcollection,stoken,bexists=bexists,&
          ssectionName=ssection)
    end if
    
    if (Iindex(4) .ne. 0) then
      ! Refinement levels
      ilevels = Rvalues(Iindex(4))%ivalue
    end if
    
    if (Iindex(5) .ne. 0) then
      ! Method
      if (Rvalues(Iindex(5))%svalue .eq. "2levelordered") then
        cmethod = 0
      else
        call output_line ("Warning: Unknown refinement method! Using default.")
        cmethod = 0
      end if
    end if

    if (.not. associated (p_rtriangulation)) then
      call output_line ("Invalid triangulation!")
      return
    end if

    ! Create the hierarchy.
    allocate(p_rmeshHierarchy)
    
    ! Create the hierarchy
    select case (cmethod)
    case (0)
      ! 2-level ordered
  
      ! Boundary present?
      if (associated(p_rboundary)) then
        call output_line ("Creating mesh hierarchy... [1",bnolinebreak=.true.)
        call mshh_initHierarchy (p_rmeshHierarchy,p_rtriangulation,0,ilevels,&
            rboundary=p_rboundary)
            
        call mshh_refineHierarchy2lv (p_rmeshHierarchy,ilevels,&
            rboundary=p_rboundary,bprint=.true.)
      else
        call output_line ("Warning: No boundary present!")
        
        call output_line ("Creating mesh hierarchy... [1",bnolinebreak=.true.)
        call mshh_initHierarchy (p_rmeshHierarchy,p_rtriangulation,0,ilevels)

        call mshh_refineHierarchy2lv (p_rmeshHierarchy,ilevels,bprint=.true.)
      end if

      call output_line ("]")
    end select
    
    ! Remove old value from collection if present
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,sname,inestlevel,ssection)
    call do_destroy(rcmdStatus,sname,ssection,.false.)
    
    ! Add to the collection
    call collct_setvalue_mshh (rcmdStatus%rcollection, sname, p_rmeshHierarchy, .true.,&
        ssectionName=ssection)

    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

    subroutine fgetDiscr(rtriangulation,rdiscr,rboundary,rcollection)
    
    !<description>
      ! Callback routine to set up a block discretisation.
    !</description>
    
    !<input>
      ! Triangulation structure
      type(t_triangulation), intent(in) :: rtriangulation
  
      ! Definition of the boundary
      type(t_boundary), intent(in), optional :: rboundary
    
      ! Collection structure with information about the discretisation
      type(t_collection), intent(inout), optional :: rcollection
    !</input>
  
    !<output>
      ! Block discretisation structure
      type(t_blockDiscretisation), intent(out) :: rdiscr
    !</output>

      ! local variables
      integer :: i,nspaces
      logical :: bexists
      integer :: celtri,celquad,ccubtri,ccubquad
      integer, dimension(:), allocatable :: p_IelementIdsTri
      integer, dimension(:), allocatable :: p_IelementIdsQuad
      integer, dimension(:), allocatable :: p_IcubIdsTri
      integer, dimension(:), allocatable :: p_IcubIdsQuad

      nspaces = rcollection%IquickAccess(1)
      allocate(p_IelementIdsTri(nspaces))
      allocate(p_IelementIdsQuad(nspaces))
      allocate(p_IcubIdsTri(nspaces))
      allocate(p_IcubIdsQuad(nspaces))
      
      call collct_getvalue_intarr (rcollection, "IelementIdsTri", p_IelementIdsTri, bexists=bexists)
      call collct_getvalue_intarr (rcollection, "IelementIdsQuad", p_IelementIdsQuad, bexists=bexists)
      call collct_getvalue_intarr (rcollection, "IcubIdsTri", p_IcubIdsTri, bexists=bexists)
      call collct_getvalue_intarr (rcollection, "IcubIdsQuad", p_IcubIdsQuad, bexists=bexists)

      select case (rtriangulation%ndim)
      case (NDIM2D)
        ! Create a block discretisation of the specific size.
        call spdiscr_initBlockDiscr (rdiscr,nspaces,rtriangulation,rboundary)
        
        ! Create the sub-discretisation structures.
        celtri  = 0
        celquad = 0
        ccubtri  = 0
        ccubquad = 0
        
        do i=1,nspaces
          celtri   = p_IelementIdsTri(i)
          celquad  = p_IelementIdsQuad(i)
          ccubtri  = p_IcubIdsTri(i)
          ccubquad = p_IcubIdsQuad(i)
        
          if (celquad .eq. 0) then
            ! Pure tri space
            call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(i), &
                INT(celtri,I32),INT(ccubtri,I32),rtriangulation, &
                rboundary)
                
          else if (celtri .eq. 0) then
            ! Pure quad space
            call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(i), &
                INT(celquad,I32),INT(ccubquad,I32),rtriangulation, &
                rboundary)
          else
            ! Combined space
            call spdiscr_initDiscr_triquad (rdiscr%RspatialDiscr(i), &
                int(celtri,I32), INT(celquad,I32),&
                int(ccubtri,I32), INT(ccubquad,I32),rtriangulation, &
                rboundary)
          end if
        end do
        
      case default
        call output_line("Error. Dimension not supported!")
        call sys_halt()
      end select

    end subroutine

  ! ***************************************************************************

    subroutine prepare_fgetDiscr(rcollection,rcmdStatus,nspaces,&
        p_IelementIdsTri,p_IelementIdsQuad,p_IcubIdsTri,p_IcubIdsQuad)
    
    !<description>
      ! Prepares setting up a discretisation.
      ! Must be called in advance to fgetDiscr.
    !</description>
    
    !<input>
      ! Number of components
      integer, intent(in) :: nspaces
    
      ! List of element id's for tri/tetra elements
      integer, dimension(:), pointer :: p_IelementIdsTri

      ! List of element id's for quad/hexa elements
      integer, dimension(:), pointer :: p_IelementIdsQuad

      ! List of cubature id's for tri/tetra elements
      integer, dimension(:), pointer :: p_IcubIdsTri

      ! List of cubature id's for quad/hexa elements
      integer, dimension(:), pointer :: p_IcubIdsQuad
    !</input>
    
    !<inputoutput>
      ! Current status block.
      type(t_commandstatus), intent(inout), target :: rcmdStatus

      ! Collection structure with information about the discretisation
      type(t_collection), intent(inout), optional :: rcollection
    !</inputoutput>

      rcollection%p_rnextCollection => rcmdStatus%rcollection
      
      rcollection%IquickAccess(1) = nspaces
      if (associated(p_IelementIdsTri)) then
        call collct_setvalue_intarr (rcollection, "IelementIdsTri", p_IelementIdsTri, .true.)
      end if
      if (associated(p_IelementIdsQuad)) then
        call collct_setvalue_intarr (rcollection, "IelementIdsQuad", p_IelementIdsQuad, .true.)
      end if
      if (associated(p_IcubIdsTri)) then
        call collct_setvalue_intarr (rcollection, "IcubIdsTri", p_IcubIdsTri, .true.)
      end if
      if (associated(p_IcubIdsQuad)) then
        call collct_setvalue_intarr (rcollection, "IcubIdsQuad", p_IcubIdsQuad, .true.)
      end if
      
    end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_fespace (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: FESPACE.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(9), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("mesh", STYPE_VAR, COLLCT_TRIA), &
         t_varspec("boundary", STYPE_VAR, COLLCT_BOUNDARY), &
         t_varspec("components", STYPE_INTEGER, COLLCT_UNDEFINED), &
         t_varspec("trielements", STYPE_STRING, COLLCT_UNDEFINED), &
         t_varspec("quadelements", STYPE_STRING, COLLCT_UNDEFINED), &
         t_varspec("tricub", STYPE_STRING, COLLCT_UNDEFINED), &
         t_varspec("quadcub", STYPE_STRING, COLLCT_UNDEFINED), &
         t_varspec("concat", STYPE_STRING, COLLCT_UNDEFINED) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    logical :: bexists,berror
    character(len=SYS_STRLEN) :: stoken
    integer :: ncomponents
    integer :: i,j,istart,iend,ilength
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    character(len=COLLCT_MLNAME) :: sname,ssection
    character(len=SYS_STRLEN) :: strielements,squadelements,stricub,squadcub,sconcat
    type(t_feSpaceLevel), pointer :: p_rfeSpace
    type(t_feSpaceLevel), pointer :: p_rfeSpace1,p_rfeSpace2
    type (t_collection) :: rcollection
    integer, dimension(:), pointer :: p_IelementIdsTri,p_IelementIdsQuad
    integer, dimension(:), pointer :: p_IcubIdsTri,p_IcubIdsQuad
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  fespace()          - Create a FE space.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("fespace - Create an FE space.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   fespace ([varname] [, ...options...])")
      call output_lbrk ()
      call output_line ("Creates an FE space that can be used to set up matrices/vectors.")
      call output_line ("[varname] identifies the variable to create.")
      call output_line ("The following options are possible in [...options...]:")
      call output_lbrk ()
      call output_line ("  ... mesh=[varname] ...")
      call output_line ("      Use mesh [varname] as coarse mesh of the hierarchy.")
      call output_line ("      Mandatory argument")
      call output_lbrk ()
      call output_line ("  ... boundary=[varname] ...")
      call output_line ("      Use boundary definition [varname] for refinements.")
      call output_lbrk ()
      call output_line ("  ... components=[num] ...")
      call output_line ("      Creates a FE space with [num] components. Mandatory argument.")
      call output_line ("      Must be defined in advance to all element/cubature specifications.")
      call output_lbrk ()
      call output_line ("  ... trielements=""EL_xxx EL_xxx EL_xxx"" ...")
      call output_line ("      Specifies a list of TRI/TETRA elemnent types to be used for")
      call output_line ("      triangular/tetrahedral elements. There must be one ID per component.")
      call output_lbrk ()
      call output_line ("  ... quadelements=""EL_xxx EL_xxx EL_xxx"" ...")
      call output_line ("      Specifies a list of QUAD/HEXA elemnent types to be used for")
      call output_line ("      quadrilateral/hexahedral elements. There must be one ID per component.")
      call output_lbrk ()
      call output_line ("  ... tricub=""CUB_xxx CUB_xxx CUB_xxx"" ...")
      call output_line ("      Specifies a list of TRI/TETRA cubature formulas to be used for")
      call output_line ("      triangular/tetrahedral elements. There must be one ID per component.")
      call output_line ("      If not defined, the default cubature formula is used.")
      call output_lbrk ()
      call output_line ("  ... quadcub=""CUB_xxx CUB_xxx CUB_xxx"" ...")
      call output_line ("      Specifies a list of QUAD/HEXA cubature formulas to be used for")
      call output_line ("      quadrilateral/hexahedral elements. There must be one ID per component.")
      call output_line ("      If not defined, the default cubature formula is used.")
      call output_lbrk ()
      call output_line ("Alternative usage:")
      call output_line ("   fespace ([varname] , concat=""[var1] [var2] ..."")")
      call output_lbrk ()
      call output_line ("Forms a FE-space by concatenation of other FE-spaces..")
      call output_line ("""[var1] [var2] ..."" must be a list of FE spaces to be concatenated.")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return

    ! Get the identifiers
    sname = Rvalues(1)%svarname
    sconcat = ""
    ncomponents = 0
    nullify(p_rboundary)
    nullify(p_rtriangulation)
    nullify(p_rfeSpace)
    nullify(p_rfeSpace1)
    nullify(p_rfeSpace2)
    
    if (Iindex(2) .ne. 0) then
      ! Get the triangulation object.
      stoken = Rvalues(Iindex(2))%svarname
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
      p_rtriangulation => collct_getvalue_tria(rcmdStatus%rcollection,stoken,ssectionName=ssection)
    end if

    if (Iindex(3) .ne. 0) then
      ! Get the bondary object.
      stoken = Rvalues(Iindex(3))%svarname
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
      p_rboundary => collct_getvalue_bdry(rcmdStatus%rcollection,stoken,bexists=bexists,&
          ssectionName=ssection)
    end if

    if (Iindex(4) .ne. 0) then
      ! Number of components
      ncomponents = Rvalues(Iindex(4))%ivalue
    end if
    
    if (Iindex(5) .ne. 0) then
      ! tri-elements
      strielements = Rvalues(Iindex(5))%svalue
    end if

    if (Iindex(6) .ne. 0) then
      ! Quad-elements
      squadelements = Rvalues(Iindex(6))%svalue
    end if

    if (Iindex(7) .ne. 0) then
      ! Tri-cubature
      stricub = Rvalues(Iindex(7))%svalue
    end if

    if (Iindex(8) .ne. 0) then
      ! Quad-Cubature
      squadcub = Rvalues(Iindex(8))%svalue
    end if

    if (Iindex(9) .ne. 0) then
      ! Concatenation
      sconcat = Rvalues(Iindex(9))%svalue
    end if
    
    
    if (sconcat .ne. "") then
      ! Concatenation of FE-spaces.
      if (cmdprs_counttokens(sconcat," ") .lt. 2) then
        call output_line("Error. Not enough source spaces specified!")
        berror = .true.
      else
      
        ! Concatenate the first two.
        istart = 0
        iend = 0

        call cmdprs_nexttoken (sconcat,istart,iend,ilength,cseparator=" ")
        sname = sconcat(istart:iend)
        p_rfeSpace1 => collct_getvalue_fesp(rcmdStatus%rcollection,sname,bexists=bexists)
        if (.not. bexists) then
          call output_line("Error. 1st source FE space does not exist.")
          berror = .true.
        end if

        if (.not. berror) then
          call cmdprs_nexttoken (sconcat,istart,iend,ilength,cseparator=" ")
          sname = sconcat(istart:iend)
          p_rfeSpace2 => collct_getvalue_fesp(rcmdStatus%rcollection,sname,bexists=bexists)
          if (.not. bexists) then
            call output_line("Error. 1st source FE space does not exist.")
            berror = .true.
          end if
        end if
        
        if (.not. berror) then
        
          ! Concat.
          call fesph_concatFeSpaces (p_rfeSpace1,p_rfeSpace2,p_rfeSpace,p_rtriangulation)
          
          ! Now the case that there are even more spaces.
          do i=3,cmdprs_counttokens(sconcat," ")
            
            ! Shift the spaces.
            p_rfeSpace1 => p_rfeSpace

            ! Get the next one.            
            call cmdprs_nexttoken (sconcat,istart,iend,ilength,cseparator=" ")
            sname = sconcat(istart:iend)
            p_rfeSpace2 => collct_getvalue_fesp(rcmdStatus%rcollection,sname,bexists=bexists)
            if (.not. bexists) then
              call output_line("Warning. Invalid source space """//trim(sname)//""". Ignored.")
              ! Ignore the error.
            else
              ! Concatenate and remove the old one.
              call fesph_concatFeSpaces (p_rfeSpace1,p_rfeSpace2,p_rfeSpace,p_rtriangulation)
              call fesph_releaseFEspace (p_rfeSpace1)
            end if
            
          end do
        
        end if
        
      end if
      
    else
    
      if (.not. associated (p_rtriangulation)) then
        call output_line ("Invalid triangulation!")
        return
      end if

      if (ncomponents .le. 0) then
        call output_line("Error. Number of components undefined!")
        return
      end if

      berror = .false.
      
      if ((strielements .ne. "") .and. (cmdprs_counttokens(strielements," ") .eq. ncomponents)) then
        allocate (p_IelementIdsTri(ncomponents))
        istart = 0
        iend = 0
        ilength = len_trim(strielements)
        do i=1,ncomponents
          ! Parse the Id.
          call cmdprs_nexttoken (strielements,istart,iend,ilength,cseparator=" ")
          p_IelementIdsTri(j) = elem_igetID(strielements(istart:iend),.true.)
          
          if (p_IelementIdsTri(j) .eq. 0) then
            call output_line("Error. Invalid element ID: "//trim(stoken))
            berror = .true.
            exit
          end if

          if (elem_igetDimension(p_IelementIdsTri(j)) .ne. NDIM2D) then
            call output_line("Error. Not a 2D element: "//trim(stoken))
            berror = .true.
            exit
          end if

          if (elem_igetNVE(p_IelementIdsTri(j)) .ne. 3) then
            call output_line("Error. Not a tri element: "//trim(stoken))
            berror = .true.
            exit
          end if
        end do
      end if

      if ((.not. berror) .and. (squadelements .ne. "")) then 
        if (cmdprs_counttokens(squadelements," ") .eq. ncomponents) then
          allocate (p_IelementIdsQuad(ncomponents))
          istart = 0
          iend = 0
          ilength = len_trim(squadelements)
          do i=1,ncomponents
            ! Parse the Id.
            call cmdprs_nexttoken (squadelements,istart,iend,ilength,cseparator=" ")
            p_IelementIdsTri(j) = elem_igetID(squadelements(istart:iend),.true.)
            
            if (p_IelementIdsquad(j) .eq. 0) then
              call output_line("Error. Invalid element ID: "//trim(stoken))
              berror = .true.
              exit
            end if

            if (elem_igetDimension(p_IelementIdsquad(j)) .ne. NDIM2D) then
              call output_line("Error. Not a 2D element: "//trim(stoken))
              berror = .true.
              exit
            end if

            if (elem_igetNVE(p_IelementIdsquad(j)) .ne. 3) then
              call output_line("Error. Not a quad element: "//trim(stoken))
              berror = .true.
              exit
            end if
          end do
        else
          call output_line("Error. Invalid parameters!")
          berror = .true.
        end if
      end if

      if ((.not. berror) .and. (stricub .ne. "")) then 
        if (cmdprs_counttokens(stricub," ") .eq. ncomponents) then
          allocate (p_IcubIdsTri(ncomponents))
          istart = 0
          iend = 0
          ilength = len_trim(stricub)
          do i=1,ncomponents
            ! Parse the Id.
            call cmdprs_nexttoken (stricub,istart,iend,ilength,cseparator=" ")
            p_IcubIdsTri(j) = cub_igetID(stricub(istart:iend))
            
            if (p_IcubIdsTri(j) .eq. 0) then
              call output_line("Error. Invalid cubature ID: "//trim(stoken))
              berror = .true.
              exit
            end if
          end do
        else
          call output_line("Error. Invalid parameters!")
          berror = .true.
        end if
      end if

      if ((.not. berror) .and. (squadcub .ne. "")) then 
        if (cmdprs_counttokens(squadcub," ") .eq. ncomponents) then
          allocate (p_IcubIdsquad(ncomponents))
          istart = 0
          iend = 0
          ilength = len_trim(squadcub)
          do i=1,ncomponents
            ! Parse the Id.
            call cmdprs_nexttoken (squadcub,istart,iend,ilength,cseparator=" ")
            p_IcubIdsquad(j) = cub_igetID(squadcub(istart:iend))
            
            if (p_IcubIdsquad(j) .eq. 0) then
              call output_line("Error. Invalid cubature ID: "//trim(stoken))
              berror = .true.
              exit
            end if
          end do
        else
          call output_line("Error. Invalid parameters!")
          berror = .true.
        end if
      end if

      ! Create the FE space usinf fgetDiscr.
      call collct_init (rcollection)
      call prepare_fgetDiscr(rcollection,rcmdStatus,ncomponents,&
          p_IelementIdsTri,p_IelementIdsQuad,p_IcubIdsTri,p_IcubIdsQuad)
      if (associated(p_rboundary)) then
        call fesph_createFEspace (p_rfeSpace,1,&
            p_rtriangulation,1,fgetDiscr,rcollection,rboundary=p_rboundary)
      else
        call output_line ("Warning: No boundary present!")
        
        call fesph_createFEspace (p_rfeSpace,1,&
            p_rtriangulation,1,fgetDiscr,rcollection)
      end if
      call collct_done (rcollection)

    end if

    if (.not. berror) then
      ! Remove old value from collection if present
      call cmdprs_getSymbolSection (rcollection,sname,inestlevel,ssection)
      call do_destroy(rcmdStatus,sname,ssection,.false.)
      
      ! Add to the collection
      call collct_setvalue_fesp (rcmdStatus%rcollection, sname, p_rfeSpace, .true.,&
          ssectionName=ssection)
    end if

    ! Release allocated memory
   if (associated(p_IelementIdsTri)) deallocate(p_IelementIdsTri)
   if (associated(p_IelementIdsQuad)) deallocate(p_IelementIdsQuad)
   if (associated(p_IcubIdsTri)) deallocate(p_IcubIdsTri)
   if (associated(p_IcubIdsQuad)) deallocate(p_IcubIdsQuad)

    if (berror) then
      return
    end if

    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_fehierarchy (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: FEHIERARCHY.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(9), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("meshhierarchy", STYPE_VAR, COLLCT_MSHHIERARCHY), &
         t_varspec("boundary", STYPE_VAR, COLLCT_BOUNDARY), &
         t_varspec("components", STYPE_INTEGER, COLLCT_UNDEFINED), &
         t_varspec("trielements", STYPE_STRING, COLLCT_UNDEFINED), &
         t_varspec("quadelements", STYPE_STRING, COLLCT_UNDEFINED), &
         t_varspec("tricub", STYPE_STRING, COLLCT_UNDEFINED), &
         t_varspec("quadcub", STYPE_STRING, COLLCT_UNDEFINED), &
         t_varspec("concat", STYPE_STRING, COLLCT_UNDEFINED) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    logical :: bexists,berror
    character(len=SYS_STRLEN) :: stoken
    integer :: ncomponents
    integer :: i,istart,iend,ilength
    type(t_boundary), pointer :: p_rboundary
    type(t_meshhierarchy), pointer :: p_rmeshhierarchy
    character(len=COLLCT_MLNAME) :: sname,ssection
    character(len=SYS_STRLEN) :: strielements,squadelements,stricub,squadcub,sconcat
    type(t_feHierarchy), pointer :: p_rfeSpHier
    type(t_feHierarchy), pointer :: p_rfeSpHier1,p_rfeSpHier2
    type (t_collection) :: rcollection
    integer, dimension(:), pointer :: p_IelementIdsTri,p_IelementIdsQuad
    integer, dimension(:), pointer :: p_IcubIdsTri,p_IcubIdsQuad
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  fehierarchy()      - Create a FE hierarchy.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("fehierarchy - Create a FE hierarchy.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   fehierarchy ([varname] [, ...options...])")
      call output_lbrk ()
      call output_line ("Creates an FE hierarchy that can be used to set up matrices/vectors.")
      call output_line ("[varname] identifies the variable to create.")
      call output_line ("The following options are possible in [...options...]:")
      call output_lbrk ()
      call output_line ("  ... meshhierarchy=[varname] ...")
      call output_line ("      Use mesh hierarchy [varname] as basis mesh of the hierarchy.")
      call output_line ("      Mandatory argument.")
      call output_lbrk ()
      call output_line ("  ... boundary=[varname] ...")
      call output_line ("      Use boundary definition [varname] for refinements.")
      call output_lbrk ()
      call output_line ("  ... components=[num] ...")
      call output_line ("      Creates a FE space with [num] components. Mandatory argument.")
      call output_line ("      Must be defined in advance to all element/cubature specifications.")
      call output_lbrk ()
      call output_line ("  ... trielements=""EL_xxx EL_xxx EL_xxx"" ...")
      call output_line ("      Specifies a list of TRI/TETRA elemnent types to be used for")
      call output_line ("      triangular/tetrahedral elements. There must be one ID per component.")
      call output_lbrk ()
      call output_line ("  ... quadelements=""EL_xxx EL_xxx EL_xxx"" ...")
      call output_line ("      Specifies a list of QUAD/HEXA elemnent types to be used for")
      call output_line ("      quadrilateral/hexahedral elements. There must be one ID per component.")
      call output_lbrk ()
      call output_line ("  ... tricub=""CUB_xxx CUB_xxx CUB_xxx"" ...")
      call output_line ("      Specifies a list of TRI/TETRA cubature formulas to be used for")
      call output_line ("      triangular/tetrahedral elements. There must be one ID per component.")
      call output_line ("      If not defined, the default cubature formula is used.")
      call output_lbrk ()
      call output_line ("  ... quadcub=""CUB_xxx CUB_xxx CUB_xxx"" ...")
      call output_line ("      Specifies a list of QUAD/HEXA cubature formulas to be used for")
      call output_line ("      quadrilateral/hexahedral elements. There must be one ID per component.")
      call output_line ("      If not defined, the default cubature formula is used.")
      call output_lbrk ()
      call output_line ("Alternative usage:")
      call output_line ("   fespace ([varname] , concat=""[var1] [var2] ..."")")
      call output_lbrk ()
      call output_line ("Forms a FE-space hierarchy by concatenation of other FE-space hierarchies.")
      call output_line ("""[var1] [var2] ..."" must be a list of FE spaces to be concatenated.")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Get the identifiers
    sname = Rvalues(1)%svarname
    sconcat = ""
    stricub = ""
    squadcub = ""
    strielements = ""
    squadelements = ""
    ncomponents = 0
    nullify(p_rboundary)
    nullify(p_rmeshHierarchy)
    
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return
    
    if (Iindex(2) .ne. 0) then
      ! Get the triangulation object.
      stoken = Rvalues(Iindex(2))%svarname
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
      p_rmeshhierarchy => collct_getvalue_mshh(rcmdStatus%rcollection,stoken,ssectionName=ssection)
    end if

    if (Iindex(3) .ne. 0) then
      ! Get the bondary object.
      stoken = Rvalues(Iindex(3))%svarname
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
      p_rboundary => collct_getvalue_bdry(rcmdStatus%rcollection,stoken,bexists=bexists,&
          ssectionName=ssection)
    end if

    if (Iindex(4) .ne. 0) then
      ! Number of components
      ncomponents = Rvalues(Iindex(4))%ivalue
    end if
    
    if (Iindex(5) .ne. 0) then
      ! tri-elements
      call cmdprs_dequoteStd(Rvalues(Iindex(5))%svalue,strielements)
    end if

    if (Iindex(6) .ne. 0) then
      ! Quad-elements
      call cmdprs_dequoteStd(Rvalues(Iindex(6))%svalue,squadelements)
    end if

    if (Iindex(7) .ne. 0) then
      ! Tri-cubature
      call cmdprs_dequoteStd(Rvalues(Iindex(7))%svalue,stricub)
    end if

    if (Iindex(8) .ne. 0) then
      ! Quad-Cubature
      call cmdprs_dequoteStd(Rvalues(Iindex(8))%svalue,squadcub)
    end if

    if (Iindex(9) .ne. 0) then
      ! Concatenation
      call cmdprs_dequoteStd(Rvalues(Iindex(9))%svalue,sconcat)
    end if
    
    if (sconcat .ne. "") then
      ! Concatenation of FE-spaces.
      if (cmdprs_counttokens(sconcat," ") .lt. 2) then
        call output_line("Error. Not enough source spaces specified!")
        berror = .true.
      else
      
        ! Concatenate the first two.
        istart = 0
        iend = 0

        call cmdprs_nexttoken (sconcat,istart,iend,ilength,cseparator=" ")
        sname = sconcat(istart:iend)
        p_rfeSpHier1 => collct_getvalue_feh(rcmdStatus%rcollection,sname,bexists=bexists)
        if (.not. bexists) then
          call output_line("Error. 1st source FE space does not exist.")
          berror = .true.
        end if

        if (.not. berror) then
          call cmdprs_nexttoken (sconcat,istart,iend,ilength,cseparator=" ")
          sname = sconcat(istart:iend)
          p_rfeSpHier2 => collct_getvalue_feh(rcmdStatus%rcollection,sname,bexists=bexists)
          if (.not. bexists) then
            call output_line("Error. 1st source FE space does not exist.")
            berror = .true.
          end if
        end if
        
        if (.not. berror) then
        
          ! Concat.
          call fesph_concatFeHierarchies (p_rfeSpHier1,p_rfeSpHier2,p_rfeSpHier)
          
          ! Now the case that there are even more spaces.
          do i=3,cmdprs_counttokens(sconcat," ")
            
            ! Shift the spaces.
            p_rfeSpHier1 => p_rfeSpHier

            ! Get the next one.            
            call cmdprs_nexttoken (sconcat,istart,iend,ilength,cseparator=" ")
            sname = sconcat(istart:iend)
            p_rfeSpHier2 => collct_getvalue_feh(rcmdStatus%rcollection,sname,bexists=bexists)
            if (.not. bexists) then
              call output_line("Warning. Invalid source space """//trim(sname)//""". Ignored.")
              ! Ignore the error.
            else
              ! Concatenate and remove the old one.
              call fesph_concatFeHierarchies (p_rfeSpHier1,p_rfeSpHier2,p_rfeSpHier)
              call fesph_releaseHierarchy (p_rfeSpHier1)
            end if
            
          end do
        
        end if
        
      end if
      
    else
    
      if (ncomponents .le. 0) then
        call output_line("Error. Number of components undefined!")
        return
      end if

      berror = .false.
      
      if ((strielements .ne. "") .and. (cmdprs_counttokens(strielements," ") .eq. ncomponents)) then
        allocate (p_IelementIdsTri(ncomponents))
        istart = 0
        iend = 0
        ilength = len_trim(strielements)
        do i=1,ncomponents
          ! Parse the Id.
          call cmdprs_nexttoken (strielements,istart,iend,ilength,cseparator=" ")
          p_IelementIdsTri(i) = elem_igetID(strielements(istart:iend),.true.)
          
          if (p_IelementIdsTri(i) .eq. 0) then
            call output_line("Error. Invalid element ID: "//trim(strielements(istart:iend)))
            berror = .true.
            exit
          end if

          if (elem_igetDimension(p_IelementIdsTri(i)) .ne. NDIM2D) then
            call output_line("Error. Not a 2D element: "//trim(strielements(istart:iend)))
            berror = .true.
            exit
          end if

          if (elem_igetNVE(p_IelementIdsTri(i)) .ne. 3) then
            call output_line("Error. Not a tri element: "//trim(strielements(istart:iend)))
            berror = .true.
            exit
          end if
        end do
      end if

      if ((.not. berror) .and. (squadelements .ne. "")) then 
        if (cmdprs_counttokens(squadelements," ") .eq. ncomponents) then
          allocate (p_IelementIdsQuad(ncomponents))
          istart = 0
          iend = 0
          ilength = len_trim(squadelements)
          do i=1,ncomponents
            ! Parse the Id.
            call cmdprs_nexttoken (squadelements,istart,iend,ilength,cseparator=" ")
            p_IelementIdsQuad(i) = elem_igetID(squadelements(istart:iend),.true.)
            
            if (p_IelementIdsQuad(i) .eq. 0) then
              call output_line("Error. Invalid element ID: "//trim(squadelements(istart:iend)))
              berror = .true.
              exit
            end if

            if (elem_igetDimension(p_IelementIdsQuad(i)) .ne. NDIM2D) then
              call output_line("Error. Not a 2D element: "//trim(squadelements(istart:iend)))
              berror = .true.
              exit
            end if

            if (elem_igetNVE(p_IelementIdsQuad(i)) .ne. 4) then
              call output_line("Error. Not a quad element: "//trim(squadelements(istart:iend)))
              berror = .true.
              exit
            end if
          end do
        else
          call output_line("Error. Invalid parameters!")
          berror = .true.
        end if
      end if

      if ((.not. berror) .and. (stricub .ne. "")) then 
        if (cmdprs_counttokens(stricub," ") .eq. ncomponents) then
          allocate (p_IcubIdsTri(ncomponents))
          istart = 0
          iend = 0
          ilength = len_trim(stricub)
          do i=1,ncomponents
            ! Parse the Id.
            call cmdprs_nexttoken (stricub,istart,iend,ilength,cseparator=" ")
            p_IcubIdsTri(i) = cub_igetID(stricub(istart:iend))
            
            if (p_IcubIdsTri(i) .eq. 0) then
              call output_line("Error. Invalid cubature ID: "//trim(stricub(istart:iend)))
              berror = .true.
              exit
            end if
          end do
        else
          call output_line("Error. Invalid parameters!")
          berror = .true.
        end if
      end if

      if ((.not. berror) .and. (squadcub .ne. "")) then 
        if (cmdprs_counttokens(squadcub," ") .eq. ncomponents) then
          allocate (p_IcubIdsquad(ncomponents))
          istart = 0
          iend = 0
          ilength = len_trim(squadcub)
          do i=1,ncomponents
            ! Parse the Id.
            call cmdprs_nexttoken (squadcub,istart,iend,ilength,cseparator=" ")
            p_IcubIdsquad(i) = cub_igetID(squadcub(istart:iend))
            
            if (p_IcubIdsquad(i) .eq. 0) then
              call output_line("Error. Invalid cubature ID: "//trim(squadcub(istart:iend)))
              berror = .true.
              exit
            end if
          end do
        else
          call output_line("Error. Invalid parameters!")
          berror = .true.
        end if
      end if

      if (.not. berror) then
        ! Create the FE space using fgetDiscr.
        call collct_init (rcollection)
        call prepare_fgetDiscr(rcollection,rcmdStatus,ncomponents,&
            p_IelementIdsTri,p_IelementIdsQuad,p_IcubIdsTri,p_IcubIdsQuad)
        allocate(p_rfeSpHier)
        if (associated(p_rboundary)) then
          call fesph_createHierarchy (p_rfeSpHier,p_rmeshHierarchy%nlevels,&
              p_rmeshHierarchy,fgetDiscr,rcollection,rboundary=p_rboundary)
        else
          call output_line ("Warning: No boundary present!")
          
          call fesph_createHierarchy (p_rfeSpHier,p_rmeshHierarchy%nlevels,&
              p_rmeshHierarchy,fgetDiscr,rcollection)
        end if
        call collct_done (rcollection)
      end if

    end if

    if (.not. berror) then
      ! Remove old value from collection if present
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,sname,inestlevel,ssection)
      call do_destroy(rcmdStatus,sname,ssection,.false.)
      
      ! Add to the collection
      call collct_setvalue_feh (rcmdStatus%rcollection, sname, p_rfeSpHier, .true.,&
          ssectionName=ssection)
    end if

    ! Release allocated memory
   if (associated(p_IelementIdsTri)) deallocate(p_IelementIdsTri)
   if (associated(p_IelementIdsQuad)) deallocate(p_IelementIdsQuad)
   if (associated(p_IcubIdsTri)) deallocate(p_IcubIdsTri)
   if (associated(p_IcubIdsQuad)) deallocate(p_IcubIdsQuad)

    if (berror) then
      return
    end if

    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_extractfespacefromhier (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: EXTRACTFESPACEFROMHIER.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(3), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_INTEGER, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_VAR, COLLCT_FEHIERARCHY) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess
    integer :: ilevel

    ! local variables
    character(len=COLLCT_MLNAME) :: sname,ssection
    type(t_fehierarchy), pointer :: p_rfeHierarchy
    type(t_feSpaceLevel), pointer :: p_rfeSpace
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  extractfespacefromhier() - Extract an FE space from a hierarchy.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("extractfespacefromhier - Extract an FE space from a hierarchy.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   extractfespacefromhier ([varname], [varhier], [level])")
      call output_lbrk ()
      call output_line ("Extracts the FE-space of level [level] from the hierarchy [varhier]")
      call output_line ("and stores it to the variable [varname].")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return

    ! Get the identifiers
    sname = Rvalues(1)%svarname
    ilevel = Rvalues(2)%ivalue
    nullify(p_rfeHierarchy)
    p_rfeHierarchy => collct_getvalue_feh (rcmdStatus%rcollection, Rvalues(3)%svarname)
    
    ! Get the FE space
    if ((ilevel .lt. 1) .or. (ilevel .gt. p_rfeHierarchy%nlevels)) then
      call output_line ("Invalid level.")
      return
    end if
    
    p_rfeSpace => p_rfeHierarchy%p_rfeSpaces(ilevel)
    
    ! Remove old value from collection if present
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,sname,inestlevel,ssection)
    call do_destroy(rcmdStatus,sname,ssection,.false.)
    
    ! Add to the collection
    call collct_setvalue_fesp (rcmdStatus%rcollection, sname, p_rfeSpace, .true.,&
        ssectionName=ssection)
        
    ! Set the user defined tag to 1 to mark this as "shared copy" so that
    ! no memory is released upon a destroy.
    call collct_settag(rcmdStatus%rcollection, sname, 1, ssectionName=ssection)
  
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine
    
  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_extractmeshfromhier (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: EXTRACTMESHFROMHIER.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(3), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_INTEGER, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_VAR, COLLCT_MSHHIERARCHY) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess
    integer :: ilevel

    ! local variables
    character(len=COLLCT_MLNAME) :: sname,ssection
    type(t_meshhierarchy), pointer :: p_rmeshHierarchy
    type(t_triangulation), pointer :: p_rtriangulation
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  extractmeshfromhier() - Extract an FE space from a hierarchy.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("extractmeshfromhier - Extract a mesh from a hierarchy.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   extractmeshfromhier ([varname], [varhier], [level])")
      call output_lbrk ()
      call output_line ("Extracts the mesh of level [level] from the mesh hierarchy [varhier]")
      call output_line ("and stores it to the variable [varname].")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return

    ! Get the identifiers
    sname = Rvalues(1)%svarname
    ilevel = Rvalues(2)%ivalue
    nullify(p_rmeshHierarchy)
    p_rmeshHierarchy => collct_getvalue_mshh (rcmdStatus%rcollection, Rvalues(3)%svarname)
    
    ! Get the FE space
    if ((ilevel .lt. 1) .or. (ilevel .gt. p_rmeshHierarchy%nlevels)) then
      call output_line ("Invalid level.")
      return
    end if
    
    p_rtriangulation => p_rmeshHierarchy%p_Rtriangulations(ilevel)
    
    ! Remove old value from collection if present
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,sname,inestlevel,ssection)
    call do_destroy(rcmdStatus,sname,ssection,.false.)
    
    ! Add to the collection
    call collct_setvalue_tria (rcmdStatus%rcollection, sname, p_rtriangulation, .true.,&
        ssectionName=ssection)
  
    ! Set the user defined tag to 1 to mark this as "shared copy" so that
    ! no memory is released upon a destroy.
    call collct_settag(rcmdStatus%rcollection, sname, 1, ssectionName=ssection)
  
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine
    
  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_extractsubvector (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: extractsubvector.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(3), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_INTEGER, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_VAR, COLLCT_BLKVECTOR) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess
    integer :: icomponent

    ! local variables
    character(len=COLLCT_MLNAME) :: sname,ssection
    type(t_vectorBlock), pointer :: p_rvectorBlock
    type(t_vectorScalar), pointer :: p_rvectorScalar
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  extractsubvector() - Extract a subvector from a block vector.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("extractsubvector - Extract a subvector from a block vector.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   extractsubvector ([varname], [varblock], [vecidx])")
      call output_lbrk ()
      call output_line ("Extracts subvector [vecidx] from the block vector [varblock]")
      call output_line ("and stores it to the variable [varname].")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return

    ! Get the identifiers
    sname = Rvalues(1)%svarname
    icomponent = Rvalues(2)%ivalue
    nullify(p_rvectorBlock)
    p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, Rvalues(3)%svarname)
    
    ! Get the FE space
    if ((icomponent .lt. 1) .or. (icomponent .gt. p_rvectorBlock%nblocks)) then
      call output_line ("Invalid component.")
      return
    end if
    
    p_rvectorScalar => p_rvectorBlock%RvectorBlock(icomponent)
    
    ! Remove old value from collection if present
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,sname,inestlevel,ssection)
    call do_destroy(rcmdStatus,sname,ssection,.false.)
    
    ! Add to the collection
    call collct_setvalue_vecsca (rcmdStatus%rcollection, sname, p_rvectorScalar, .true.,&
        ssectionName=ssection)
  
    ! Set the user defined tag to 1 to mark this as "shared copy" so that
    ! no memory is released upon a destroy.
    call collct_settag(rcmdStatus%rcollection, sname, 1, ssectionName=ssection)
  
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_mlevelprjhierarchy (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: MLEVELPRJHIERARCHY.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>

    ! Command arguments
    type(t_varspec), dimension(2), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("fehierarchy", STYPE_VAR, COLLCT_FEHIERARCHY)  &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess
  
    ! local variables
    logical :: bexists
    integer :: i
    character(len=COLLCT_MLNAME) :: sname,ssection
    character(len=SYS_STRLEN) :: stoken
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_interlevelProjectionHier), pointer :: p_rprjHier
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  mlevelprjhierarchy() - Create a multilevel projection hierarchy.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("mlevelprjhierarchy - Create an multilevel projection hierarchy.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   mlevelprjhierarchy ([varname] [, ...options...])")
      call output_lbrk ()
      call output_line ("Creates a multilevel projection hierarchy that can be used to transfer")
      call output_line ("soltion/rhs vectors from one level to another.")
      call output_line ("[varname] identifies the variable to create.")
      call output_line ("The following options are possible in [...options...]:")
      call output_lbrk ()
      call output_line ("  ... fehierarchy=[varname] ...")
      call output_line ("      Create a multilevel projection hierarchy based on the")
      call output_line ("      FE space hierarchy [varname]. Mandatory argument")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return

    ! Get the identifiers
    sname = Rvalues(1)%svarname

    if (Iindex(2) .ne. 0) then
      ! Get the triangulation object.
      stoken = Rvalues(Iindex(2))%svarname
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
      p_rfeHierarchy => collct_getvalue_feh(rcmdStatus%rcollection,stoken,bexists=bexists,ssectionName=ssection)
    end if
    
    if (.not. associated (p_rfeHierarchy)) then
      call output_line ("Invalid FE hierarchy!")
      return
    end if
  
    ! Create the fe space.
    allocate(p_rprjHier)
    
    call output_line ("Creating multilevel projection hierarchy.")
  
    ! Create the FE space usinf fgetDiscr.
    call mlprj_initPrjHierarchy(p_rprjHier,1,p_rfeHierarchy%nlevels)
    do i = 1,p_rfeHierarchy%nlevels
      call mlprj_initPrjHierarchyLevel(p_rprjHier,i,&
          p_rfeHierarchy%p_rfeSpaces(i)%p_rdiscretisation)
    end do
    call mlprj_commitPrjHierarchy (p_rprjHier)

    ! Remove old value from collection if present
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,sname,inestlevel,ssection)
    call do_destroy(rcmdStatus,sname,ssection,.false.)
    
    ! Add to the collection
    call collct_setvalue_mlprjh (rcmdStatus%rcollection, sname, p_rprjHier, .true.,&
        ssectionName=ssection)
    
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_createblockvector (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: CREATEBLOCKVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(2), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_VAR, COLLCT_FESPACE) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    logical :: bexists
    character(len=COLLCT_MLNAME) :: sname,ssection
    character(len=SYS_STRLEN) :: stoken
    type(t_feSpaceLevel), pointer :: p_rfeSpace
    type(t_vectorBlock), pointer :: p_rvectorBlock
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  createblockvector() - Create an empty block vector.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("createblockvector - Create an empty block vector.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   createblockvector ([varname],[varfespace])")
      call output_lbrk ()
      call output_line ("Creates an empty block vector..")
      call output_line ("[varname] identifies the variable to create.")
      call output_line ("[varfespace] defines the FE space, the vector is based on.")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return

    ! Get the identifiers
    sname = Rvalues(1)%svarname
    
    stoken = Rvalues(2)%svarname
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
    p_rfeSpace => collct_getvalue_fesp(rcmdStatus%rcollection,stoken,bexists=bexists,ssectionName=ssection)
    
    ! Create the vector.
    allocate(p_rvectorBlock)
    call lsysbl_createVectorBlock (p_rfeSpace%p_rdiscretisation,p_rvectorBlock,.true.)

    ! Remove old value from collection if present
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,sname,inestlevel,ssection)
    call do_destroy(rcmdStatus,sname,ssection,.false.)
    
    ! Add to the collection
    call collct_setvalue_vec (rcmdStatus%rcollection, sname, p_rvectorBlock, .true., &
        ssectionName=ssection)
    
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_readblockvector (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: READBLOCKVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(3), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_BLKVECTOR), &
         t_varspec("", STYPE_STRING, COLLCT_UNDEFINED), &
         t_varspec("format", STYPE_STRING, COLLCT_UNDEFINED) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    character(len=SYS_STRLEN) :: sfilename
    character(len=COLLCT_MLNAME) :: sname,ssection
    character(len=SYS_STRLEN) :: stoken
    type(t_vectorBlock), pointer :: p_rvectorBlock
    logical :: bformatted
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  readblockvector()  - Read a block vector from a file.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
        call output_line ("readblockvector - Read block vector from file.")
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   readblockvector ([varname],[filename] [,...options...]")
        call output_lbrk ()
        call output_line ("Read a block vector from a file.")
        call output_line ("[varname] identifies the variable where to read data to.")
        call output_line ("[filename] identifies the filename.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... format=""unformatted"" ...")
        call output_line ("      Read a binary file. Default is formatted, human readable file.")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return

    ! Get the identifiers
    sname = Rvalues(1)%svarname
    call cmdprs_dequoteStd(Rvalues(2)%svalue,sfilename)
    bformatted = .true.
    
    if (Iindex(3) .ne. 0) then
      call cmdprs_dequoteStd(Rvalues(Iindex(3))%svalue,stoken)
      bformatted = sys_upcase(stoken) .ne. "UNFORMATTED"
    end if
    
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,sname,inestlevel,ssection)
    p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, sname, ssectionName=ssection)
    
    call output_line ("Reading vector: "//trim(sfilename))

    call vecio_readBlockVectorHR (p_rvectorBlock, stoken, .true.,&
        0, sfilename, bformatted)
    
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_writeblockvector (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: WRITEBLOCKVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(3), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_STRING, COLLCT_UNDEFINED), &
         t_varspec("format", STYPE_STRING, COLLCT_UNDEFINED) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    character(len=SYS_STRLEN) :: sfilename,sformat
    character(len=COLLCT_MLNAME) :: sname,ssection
    character(len=SYS_STRLEN) :: stoken
    type(t_vectorBlock), pointer :: p_rvectorBlock
    logical :: bformatted
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  writeblockvector() - Write a block vector to a file.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("writeblockvector - Write block vector from file.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   writeblockvector ([varname],[filename] [,...options...]")
      call output_lbrk ()
      call output_line ("Writes a block vector to a file.")
      call output_line ("[varname] identifies the variable where to read data to.")
      call output_line ("[filename] identifies the filename.")
      call output_line ("The following options are possible in [...options...]:")
      call output_lbrk ()
      call output_line ("  ... format=""unformatted"" ...")
      call output_line ("      Write a binary file. Default is formatted, human readable file.")
      call output_lbrk ()
      call output_line ("  ... format=""(E20.10)"" ...")
      call output_line ("      Write a human readable file with the specified number format.")
  
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return

    ! Get the identifiers
    sname = Rvalues(1)%svarname
    call cmdprs_dequoteStd(Rvalues(2)%svalue,sfilename)
    bformatted = .true.
    sformat = "(E20.10)"
    
    if (Iindex(3) .ne. 0) then
      bformatted = Rvalues(Iindex(3))%svalue .ne. "unformatted"
      if (.not. bformatted) sformat = Rvalues(Iindex(3))%svalue
    end if
    
    stoken = Rvalues(2)%svarname
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,sname,inestlevel,ssection)
    p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, sname, ssectionName=ssection)
    
    call output_line ("Writing vector: "//trim(sfilename))

    if (bformatted) then
      call vecio_writeBlockVectorHR (p_rvectorBlock, "vector", .true.,&
          0, sfilename, sformat)
    else
      call vecio_writeBlockVectorHR (p_rvectorBlock, "vector", .true.,&
          0, sfilename)
    end if
    
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_copyvector (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: copyvector.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(2), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_VAR, COLLCT_UNDEFINED) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    character(len=COLLCT_MLNAME) :: ssection
    character(len=SYS_STRLEN) :: stoken
    integer :: ctype1,ctype2
    type(t_vectorBlock), pointer :: p_rvectorBlock1,p_rvectorBlock2
    type(t_vectorScalar), pointer :: p_rvectorScalar1,p_rvectorScalar2
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  copyvector()       - Copy a vector to another.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("copyvector - Copy a scalar or block vector.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   copyvector ([varsource],[vardest])")
      call output_lbrk ()
      call output_line ("Copies vector [varsource] to vector [vardest].")
      call output_line ("The vector may be a full block vector or a single.")
      call output_line ("subvector.")
      call output_lbrk ()
      call output_line ("Example:")
      call output_line ("    copyvector (source,dest)")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return

    stoken = Rvalues(1)%svarname
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
    ctype1 = collct_gettype (rcmdStatus%rcollection, Rvalues(1)%svarname, ssectionName=ssection)

    if (ctype1 .eq. COLLCT_BLKVECTOR) then
    
      p_rvectorBlock1 => collct_getvalue_vec (rcmdStatus%rcollection, stoken,&
            ssectionName=ssection)
    
      stoken = Rvalues(2)%svarname
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
      ctype2 = collct_gettype (rcmdStatus%rcollection, Rvalues(2)%svarname, ssectionName=ssection)
      if (ctype2 .ne. COLLCT_BLKVECTOR) then
        call output_line ("Source and destination not compatible.")
        return
      end if
      
      p_rvectorBlock2 => collct_getvalue_vec (rcmdStatus%rcollection, stoken,&
            ssectionName=ssection)
            
      ! Copy!
      call lsysbl_copyVector (p_rvectorBlock1,p_rvectorBlock2)
      
    else if (ctype1 .eq. COLLCT_SCAVECTOR) then

      p_rvectorScalar1 => collct_getvalue_vecsca (rcmdStatus%rcollection, stoken,&
            ssectionName=ssection)
    
      stoken = Rvalues(2)%svarname
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
      ctype2 = collct_gettype (rcmdStatus%rcollection, Rvalues(2)%svarname, ssectionName=ssection)
      if (ctype2 .ne. COLLCT_BLKVECTOR) then
        call output_line ("Source and destination not compatible.")
        return
      end if

      p_rvectorScalar2 => collct_getvalue_vecsca (rcmdStatus%rcollection, stoken,&
            ssectionName=ssection)
            
      ! Copy!
      call lsyssc_copyVector (p_rvectorScalar1,p_rvectorScalar2)
      
    else
      call output_line ("Invalid source vector.")
      return
    end if
    
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_interpolatevector (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: interpolate.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(6), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_INTEGER, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_INTEGER, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_VAR, COLLCT_FEHIERARCHY), &
         t_varspec("", STYPE_VAR, COLLCT_MLPRJHIERARCHY) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    character(len=COLLCT_MLNAME) :: ssection
    character(len=SYS_STRLEN) :: stoken
    type(t_vectorBlock), pointer :: p_rvectorBlock1,p_rvectorBlock2
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_interlevelProjectionHier), pointer :: p_rprjHier
    type(t_vectorBlock), pointer :: p_rtemp1,p_rtemp2
    integer :: i,isource,idest
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  interpolatevector() - Interpolate a vector to another level.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("interpolatevector - Interpolate a vector to another level.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   interpolatevector ([varsource],[lvlsource],[vardest],[lvldest],[varhier],[mlprj])")
      call output_lbrk ()
      call output_line ("Interpolates vector [varsource] from level [lvlsource] to")
      call output_line ("level [lvldest] and writes the result to [vardest].")
      call output_line ("[varhier] specifies the FE hierarchy, the level refer to.")
      call output_line ("[mlprj] specifies the projection hierarchy to be used.")
      call output_line ("The method uses prolongation/interpolation for the level change.")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return
    
    ! Get the variables.
    stoken = Rvalues(1)%svarname
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
    p_rvectorBlock1 => collct_getvalue_vec (rcmdStatus%rcollection, stoken,&
          ssectionName=ssection)
          
    isource = Rvalues(2)%ivalue

    stoken = Rvalues(3)%svarname
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
    p_rvectorBlock2 => collct_getvalue_vec (rcmdStatus%rcollection, stoken,&
          ssectionName=ssection)
          
    idest = Rvalues(4)%ivalue
    
    stoken = Rvalues(5)%svarname
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
    p_rfeHierarchy => collct_getvalue_feh (rcmdStatus%rcollection, stoken,&
          ssectionName=ssection)

    stoken = Rvalues(6)%svarname
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
    p_rprjHier => collct_getvalue_mlprjh (rcmdStatus%rcollection, stoken,&
          ssectionName=ssection)
    
    ! Go either up or down - or copy.
    if (isource .eq. idest) then
    
      call lsysbl_copyVector (p_rvectorBlock1,p_rvectorBlock2)
      
    else if (isource .lt. idest) then
    
      if (isource .eq. idest-1) then
      
        ! One prolongation
        p_rtemp1 => p_rvectorBlock1
      
      else
      
        if (.not. associated (p_rfeHierarchy)) then
          call output_line ("Invalid FE hierarchy!")
          return
        end if
      
        ! Multiple prolongations
        allocate (p_rtemp1)
        call lsysbl_createVectorBlock (p_rvectorBlock1,p_rtemp1,.false.)
        call lsysbl_copyVector (p_rvectorBlock1,p_rtemp1)
        
        ! Create the temp vectors using the FE hierarchy.
        do i=isource+1,idest-1
          allocate (p_rtemp2)
          call lsysbl_createVectorBlock (p_rfeHierarchy%p_rfeSpaces(i)%p_rdiscretisation,&
              p_rtemp2,.false.)
          call mlprj_performProlongationHier (p_rprjHier,&
              i,p_rtemp1,p_rtemp2)
          call lsysbl_releaseVector (p_rtemp1)
          deallocate (p_rtemp1)
          p_rtemp1 => p_rtemp2
        end do
      end if

      ! Final prolongation
      call mlprj_performProlongationHier (p_rprjHier,&
          idest,p_rtemp1,p_rvectorBlock2)
          
      if (isource .eq. idest-1) then
        ! Cleanup
        call lsysbl_releaseVector (p_rtemp1)
        deallocate(p_rtemp1)
      end if
    
    else
      ! Interpolation. NOT RESTRICTION!!!
    
      if (isource-1 .eq. idest) then
      
        ! One prolongation
        p_rtemp1 => p_rvectorBlock1
      
      else
      
        if (.not. associated (p_rfeHierarchy)) then
          call output_line ("Invalid FE hierarchy!")
          return
        end if
        
        ! Multiple interpolations
        allocate (p_rtemp1)
        call lsysbl_createVectorBlock (p_rvectorBlock1,p_rtemp1,.false.)
        call lsysbl_copyVector (p_rvectorBlock1,p_rtemp1)
        
        ! Create the temp vectors using the FE hierarchy.
        do i=idest-1,isource+1,-1
          allocate (p_rtemp2)
          call lsysbl_createVectorBlock (p_rfeHierarchy%p_rfeSpaces(i)%p_rdiscretisation,&
              p_rtemp2,.false.)
          call mlprj_performInterpolationHier (p_rprjHier,&
              i+1,p_rtemp2,p_rtemp1)
          call lsysbl_releaseVector (p_rtemp1)
          deallocate (p_rtemp1)
          p_rtemp1 => p_rtemp2
        end do
      end if

      ! Final interpolation
      call mlprj_performInterpolationHier (p_rprjHier,&
          idest+1,p_rvectorBlock2,p_rtemp1)
          
      if (isource-1 .eq. idest) then
        ! Cleanup
        call lsysbl_releaseVector (p_rtemp1)
        deallocate(p_rtemp1)
      end if
      
    end if
    
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine
  
  ! ***************************************************************************

    subroutine fcoeff_analytPrj (rdiscretisation, rform, &
                  nelements, npointsPerElement, Dpoints, &
                  IdofsTest, rdomainIntSubset, &
                  Dcoefficients, rcollection)
    
    ! Returns values of an FE function in cubature points.
    
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    type(t_linearForm), intent(in) :: rform
    integer, intent(in) :: nelements
    integer, intent(in) :: npointsPerElement
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
    integer, dimension(:,:), intent(in) :: IdofsTest
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset    
    type(t_collection), intent(inout), optional :: rcollection
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients

      ! local variables
      integer :: icomponent
      type(t_vectorBlock), pointer :: p_rvectorBlock
      
      ! Get the component and the FE function
      p_rvectorBlock => rcollection%p_rvectorQuickAccess1
      icomponent = rcollection%IquickAccess(1)

      ! Evaluate the FE function      
      call fevl_evaluate_sim (p_rvectorBlock%RvectorBlock(icomponent), &
          rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
  
    end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_l2projection (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: L2PROJECTION.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(4), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("verbose", STYPE_INTEGER, COLLCT_UNDEFINED), &
         t_varspec("relerror", STYPE_DOUBLE, COLLCT_UNDEFINED) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    character(len=COLLCT_MLNAME) :: ssection
    character(len=SYS_STRLEN) :: stoken
    type(t_matrixScalar) :: rmatrixMass
    integer :: i,ctype
    logical :: bverbose
    type(t_collection) :: rcollection
    type(t_configL2ProjectionByMass) :: rL2ProjectionConfig
    type(t_vectorBlock), pointer :: p_rvectorBlock1,p_rvectorBlock2
    type(t_vectorScalar), pointer :: p_rvectorScalar1,p_rvectorScalar2
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  l2projection()     - Appies an L2 projection to a vector.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("l2projection       - Appies an L2 projection to a vector.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   l2projection ([varsource],[vardest] [,...options...])")
      call output_lbrk ()
      call output_line ("Interpolates vector [varsource] to [vardest] using an")
      call output_line ("L2-projection.")
      call output_lbrk ()
      call output_line ("Example:")
      call output_line ("    l2projection (source,dest)")
      call output_lbrk ()
      call output_line ("The following options are possible in [...options...]:")
      call output_lbrk ()
      call output_line ("  ... verbose==1 ...")
      call output_line ("      Activate/Deactivate verbose output.")
      call output_lbrk ()
      call output_line ("  ... relerror=[error] ...")
      call output_line ("      Defines the relative accuracy of the projection.")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return
    
    ! Get the parameters
    nullify(p_rvectorBlock1)
    nullify(p_rvectorBlock2)
    nullify(p_rvectorScalar1)
    nullify(p_rvectorScalar2)
    bverbose = .true.
    
    if (Iindex(3) .ne. 0) then
      bverbose = Rvalues(Iindex(3))%ivalue .ne. 0
    end if

    if (Iindex(4) .ne. 0) then
      rL2ProjectionConfig%depsrel = Rvalues(Iindex(4))%dvalue
    end if
        
    ! Vectors may be block or scalar
    stoken = Rvalues(1)%svarname
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
    ctype = collct_gettype (rcollection, stoken, ssectionName=ssection)
    
    if (ctype .eq. COLLCT_SCAVECTOR) then
      p_rvectorScalar1 => collct_getvalue_vecsca (rcmdStatus%rcollection, stoken,&
          ssectionName=ssection)
    else if (ctype .eq. COLLCT_BLKVECTOR) then
      p_rvectorBlock1 => collct_getvalue_vec (rcmdStatus%rcollection, stoken,&
          ssectionName=ssection)
    else
      call output_line ("Invalid source vector")
      return
    end if

    stoken = Rvalues(2)%svarname
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
    ctype = collct_gettype (rcollection, stoken, ssectionName=ssection)
    
    if (ctype .eq. COLLCT_SCAVECTOR) then
      p_rvectorScalar2 => collct_getvalue_vecsca (rcmdStatus%rcollection, stoken,&
          ssectionName=ssection)
    else if (ctype .eq. COLLCT_BLKVECTOR) then
      p_rvectorBlock2 => collct_getvalue_vec (rcmdStatus%rcollection, stoken,&
          ssectionName=ssection)
    else
      call output_line ("Invalid destination vector")
      return
    end if

    if ((associated(p_rvectorScalar1) .and. .not. associated(p_rvectorScalar2)) .or. &
        (associated(p_rvectorBlock1) .and. .not. associated(p_rvectorBlock2)) .or. &
        (associated(p_rvectorScalar2) .and. .not. associated(p_rvectorScalar1)) .or. &
        (associated(p_rvectorBlock2) .and. .not. associated(p_rvectorBlock1))) then
      call output_line ("Source and destination vector not compatible.")
      return
    end if

    if (associated(p_rvectorScalar1)) then
      ! Scalar projection
      call lsyssc_clearVector (p_rvectorScalar2)
      
      if (bverbose) then
        call output_line ("Creating mass matrix - structure...")
      end if

      ! Create a mass matrix in that space
      call bilf_createMatrixStructure (p_rvectorScalar2%p_rspatialDiscr,&
          LSYSSC_MATRIX9,rmatrixMass)

      if (bverbose) then
        call output_line ("Creating mass matrix - content...")
      end if

      call stdop_assembleSimpleMatrix (rmatrixMass,DER_FUNC,DER_FUNC,1.0_DP,.true.)
      
      if (bverbose) then
        call output_line ("Projecting...")
      end if

      ! Do the L2 projection. Put the vector in a temporary block vector.
      allocate (rcollection%p_rvectorQuickAccess1)
      call lsysbl_createVecFromScalar (p_rvectorScalar1,rcollection%p_rvectorQuickAccess1)
      rcollection%IquickAccess(1) = 1
      
      call anprj_analytL2projectionByMass (p_rvectorScalar2, rmatrixMass,&
          fcoeff_analytPrj, rcollection, rL2ProjectionConfig)
          
      call lsysbl_releaseVector (rcollection%p_rvectorQuickAccess1)
      deallocate (rcollection%p_rvectorQuickAccess1)

      if (bverbose) then
        call output_line ("Rel. error: "//trim(sys_sdEL(rL2ProjectionConfig%drelError,10)))
        call output_line ("Abs. error: "//trim(sys_sdEL(rL2ProjectionConfig%dabsError,10)))
        call output_line ("Iteraions : "//trim(sys_siL(rL2ProjectionConfig%iiterations,10)))
      end if
          
      ! Release the mass matrix
      call lsyssc_releaseMatrix (rmatrixMass)      
      
    else
      ! Block projection. All blocks separately.
        
      ! Clear the destination
      call lsysbl_clearVector (p_rvectorBlock2)

      ! Loop through the components
      do i=1,min(p_rvectorBlock1%nblocks,p_rvectorBlock2%nblocks)
        if (bverbose) then
          call output_lbrk ()
          call output_line ("Component : "//trim(sys_siL(i,10)))
          call output_line ("Creating mass matrix - structure...")
        end if

        ! Create a mass matrix in that space
        call bilf_createMatrixStructure (p_rvectorBlock2%p_rblockDiscr%RspatialDiscr(i),&
            LSYSSC_MATRIX9,rmatrixMass)

        if (bverbose) then
          call output_line ("Creating mass matrix - content...")
        end if

        call stdop_assembleSimpleMatrix (rmatrixMass,DER_FUNC,DER_FUNC,1.0_DP,.true.)
        
        if (bverbose) then
          call output_line ("Projecting...")
        end if

        ! Do the L2 projection
        rcollection%p_rvectorQuickAccess1 => p_rvectorBlock1
        rcollection%IquickAccess(1) = i
        call anprj_analytL2projectionByMass (p_rvectorBlock2%RvectorBlock(i), rmatrixMass,&
            fcoeff_analytPrj, rcollection, rL2ProjectionConfig)

        if (bverbose) then
          call output_line ("Rel. error: "//trim(sys_sdEL(rL2ProjectionConfig%drelError,10)))
          call output_line ("Abs. error: "//trim(sys_sdEL(rL2ProjectionConfig%dabsError,10)))
          call output_line ("Iteraions : "//trim(sys_siL(rL2ProjectionConfig%iiterations,10)))
        end if
            
        ! Release the mass matrix
        call lsyssc_releaseMatrix (rmatrixMass)
      end do
      
    end if

    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_writeucd (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: writeucd.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(3), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_STRING, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_VAR, COLLCT_TRIA), &
         t_varspec("", STYPE_STRING, COLLCT_UNDEFINED) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    character(len=COLLCT_MLNAME) :: ssection
    character(len=SYS_STRLEN) :: sfilename,stype,svecname,sname,stoken
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_vectorScalar), pointer :: p_rvectorScalar
    real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2,p_Ddata3
    type(t_ucdExport) :: rexport
    integer :: istart,iend,ilength,iparam,ncomponents
    logical :: bexists
    integer :: ctype
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  writeucd()         - Writes a postprocessing file.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("writeucd - Write a postprocessing file.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   writeucd ([filename],[mesh],[type] [,...options...])")
      call output_lbrk ()
      call output_line ("Writes a postproceessing file [filename] based on the")
      call output_line ("mesh [mesh]. [type] specifies the type of the output.")
      call output_lbrk ()
      call output_line ("Example:")
      call output_line ("    writeucd (""mypostproc.vtk"",mymesh,""vtk"")")
      call output_lbrk ()
      call output_line ("The following postprocessing output is possible:")
      call output_lbrk ()
      call output_line ("  [type] = ""vtk"": VTK-output")
      call output_line ("           ""gmv"": GMV output")
      call output_lbrk ()
      call output_line ("The following parameters are available in [...options...]:")
      call output_lbrk ()
      call output_line ("  pointdatascalar=""[name] [vector]""")
      call output_line ("          Writes out subvector of vector [vector]")
      call output_line ("          as scalar array with the name [name]. May be specified")
      call output_line ("          more than once. The data is interpreted in the corner")
      call output_line ("          points of the elements. Example:")
      call output_line ("              pointdatascalar=""vel_y myvec""")
      call output_lbrk ()
      call output_line ("  pointdatavec=""[name] [vector1] [vector2] ...""")
      call output_line ("          Writes out a set of subvectors [vector1] [vector2]...")
      call output_line ("          as a vector field with name [name]. May be specified more than once.")
      call output_line ("          Example:")
      call output_line ("              pointdatavec=""velocity myvec""")
      call output_lbrk ()
      call output_line ("  celldatascalar=""[name] [vector]""")
      call output_line ("          Writes out vector [vector] as scalar array with the name ")
      call output_line ("          [name]. May be specified more than once.")
      call output_line ("          The data is interpreted in the elements.")
      call output_line ("          Example:")
      call output_line ("              celldatascalar=""pressure myvec""")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return
    
    call cmdprs_dequoteStd(Rvalues(1)%svalue,sfilename)
    
    stoken = Rvalues(2)%svarname
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
    p_rtriangulation => collct_getvalue_tria (rcmdStatus%rcollection, stoken,&
          ssectionName=ssection)
    
    ! Start the output file
    call cmdprs_dequoteStd(sys_upcase(Rvalues(3)%svalue),stype)
    if (stype .eq. "VTK") then
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfilename)
    else if (stype .eq. "GMV") then
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfilename)
    else
      call output_line ("Unknown output format!")
      return
    end if
    
    ! Loop through the other arguments
    do iparam = 4,size(Rvalues)
    
      if (Rvalues(iparam)%svartag .eq. "pointdatascalar") then
      
        ! Get the data. Ignore if wrong.
        if (Rvalues(iparam)%ctype .ne. STYPE_STRING) then
          call output_line ("Ignoring invalid parameter.")
          cycle
        end if
        call cmdprs_dequoteStd(Rvalues(iparam)%svalue,stoken)
        
        ! Must have two entries
        if (cmdprs_counttokens (stoken," ") .ne. 2) then
          call output_line ("Ignoring invalid parameter.")
          cycle
        end if

        istart = 0
        iend = 0
        ilength = len_trim(stoken)
        call cmdprs_nexttoken (stoken,istart,iend,ilength,cseparator=" ")
        sname = stoken(istart:iend)
        call cmdprs_nexttoken (stoken,istart,iend,ilength,cseparator=" ")
        svecname = stoken(istart:iend)
        
        ! Check the type
        call cmdprs_getSymbolSection (rcmdStatus%rcollection,svecname,inestlevel,ssection,bexists)
        if (.not. bexists) then
          call output_line ("Ignoring unknown vector: "//trim(svecname))
          cycle
        end if

        ctype = collct_gettype (rcmdStatus%rcollection, svecname, ssectionName=ssection)
        if (ctype .ne. COLLCT_SCAVECTOR) then
          call output_line ("Ignoring unknown type: "//trim(svecname))
          cycle
        end if
        
        p_rvectorScalar => collct_getvalue_vecsca (rcmdStatus%rcollection, svecname, &
            ssectionName=ssection)

        ! Write the block
        call spdp_projectToVertices (p_rvectorScalar,p_Ddata1)
        call ucd_addVariableVertexBased (rexport,trim(sname),UCD_VAR_STANDARD,p_Ddata1)
        deallocate(p_Ddata1)
      
      else if (Rvalues(iparam)%svartag .eq. "pointdatavec") then

        ! Get the data. Ignore if wrong.
        if (Rvalues(iparam)%ctype .ne. STYPE_STRING) then
          call output_line ("Ignoring invalid parameter.")
          cycle
        end if
        call cmdprs_dequoteStd(Rvalues(iparam)%svalue,stoken)
        
        ! Must have at most 4 entries
        if (cmdprs_counttokens (stoken," ") .gt. 4) then
          call output_line ("Ignoring invalid parameter.")
          cycle
        end if

        ! Must have at most 4 entries
        if (cmdprs_counttokens (stoken," ") .lt. 2) then
          call output_line ("Ignoring invalid parameter.")
          cycle
        end if

        ! Up to three components
        nullify(p_Ddata1)
        nullify(p_Ddata2)
        nullify(p_Ddata3)
        ncomponents = 0

        ! Get the name and the components
        istart = 0
        iend = 0
        ilength = len_trim(stoken)
        call cmdprs_nexttoken (stoken,istart,iend,ilength,cseparator=" ")
        sname = stoken(istart:iend)
        
        ! 1st component
        call cmdprs_nexttoken (stoken,istart,iend,ilength,cseparator=" ")
        svecname = stoken(istart:iend)
        
        ! Check the type
        call cmdprs_getSymbolSection (rcmdStatus%rcollection,svecname,inestlevel,ssection,bexists)
        if (.not. bexists) then
          call output_line ("Ignoring unknown vector: "//trim(svecname))
          cycle
        end if

        ctype = collct_gettype (rcmdStatus%rcollection, svecname, ssectionName=ssection)
        if (ctype .ne. COLLCT_SCAVECTOR) then
          call output_line ("Ignoring unknown type: "//trim(svecname))
          cycle
        end if
        
        p_rvectorScalar => collct_getvalue_vecsca (rcmdStatus%rcollection, svecname, &
            ssectionName=ssection)

        ! Write the block
        call spdp_projectToVertices (p_rvectorScalar,p_Ddata1)
        
        ncomponents = ncomponents + 1
        
        call cmdprs_nexttoken (stoken,istart,iend,ilength,cseparator=" ")

        if (istart .ne. 0) then
          ! 2nd component

          svecname = stoken(istart:iend)
          
          ! Check the type
          call cmdprs_getSymbolSection (rcmdStatus%rcollection,svecname,inestlevel,ssection,bexists)
          if (.not. bexists) then
            call output_line ("Ignoring unknown vector: "//trim(svecname))
            cycle
          end if

          ctype = collct_gettype (rcmdStatus%rcollection, svecname, ssectionName=ssection)
          if (ctype .ne. COLLCT_SCAVECTOR) then
            call output_line ("Ignoring unknown type: "//trim(svecname))
            cycle
          end if
          
          p_rvectorScalar => collct_getvalue_vecsca (rcmdStatus%rcollection, svecname, &
              ssectionName=ssection)

          ! Write the block
          call spdp_projectToVertices (p_rvectorScalar,p_Ddata2)
          
          ncomponents = ncomponents + 1
          
          call cmdprs_nexttoken (stoken,istart,iend,ilength,cseparator=" ")
        end if

        if (istart .ne. 0) then
          ! 3rd component

          svecname = stoken(istart:iend)
          
          ! Check the type
          call cmdprs_getSymbolSection (rcmdStatus%rcollection,svecname,inestlevel,ssection,bexists)
          if (.not. bexists) then
            call output_line ("Ignoring unknown vector: "//trim(svecname))
            cycle
          end if

          ctype = collct_gettype (rcmdStatus%rcollection, svecname, ssectionName=ssection)
          if (ctype .ne. COLLCT_SCAVECTOR) then
            call output_line ("Ignoring unknown type: "//trim(svecname))
            cycle
          end if
          
          p_rvectorScalar => collct_getvalue_vecsca (rcmdStatus%rcollection, svecname, &
              ssectionName=ssection)

          ! Write the block
          call spdp_projectToVertices (p_rvectorScalar,p_Ddata3)
          
          ncomponents = ncomponents + 1
          
          call cmdprs_nexttoken (stoken,istart,iend,ilength,cseparator=" ")
        end if
      
        ! Write the block
        select case (ncomponents)
        case (1)
          call ucd_addVarVertBasedVec (rexport,trim(sname),p_Ddata1)
          deallocate(p_Ddata1)
        case (2)
          call ucd_addVarVertBasedVec (rexport,trim(sname),p_Ddata1,p_Ddata2)
          deallocate(p_Ddata2)
          deallocate(p_Ddata1)
        case (3:)
          call ucd_addVarVertBasedVec (rexport,trim(sname),p_Ddata1,p_Ddata2,p_Ddata3)
          deallocate(p_Ddata3)
          deallocate(p_Ddata2)
          deallocate(p_Ddata1)
        end select
      
      else if (Rvalues(iparam)%svartag .eq. "celldatascalar") then
    
        ! Get the data. Ignore if wrong.
        if (Rvalues(iparam)%ctype .ne. STYPE_STRING) then
          call output_line ("Ignoring invalid parameter.")
          cycle
        end if
        call cmdprs_dequoteStd(Rvalues(iparam)%svalue,stoken)
        
        ! Must have two entries
        if (cmdprs_counttokens (stoken," ") .ne. 2) then
          call output_line ("Ignoring invalid parameter.")
          cycle
        end if

        istart = 0
        iend = 0
        ilength = len_trim(stoken)
        call cmdprs_nexttoken (stoken,istart,iend,ilength,cseparator=" ")
        sname = stoken(istart:iend)
        call cmdprs_nexttoken (stoken,istart,iend,ilength,cseparator=" ")
        svecname = stoken(istart:iend)
        
        ! Check the type
        call cmdprs_getSymbolSection (rcmdStatus%rcollection,svecname,inestlevel,ssection,bexists)
        if (.not. bexists) then
          call output_line ("Ignoring unknown vector: "//trim(svecname))
          cycle
        end if

        ctype = collct_gettype (rcmdStatus%rcollection, svecname, ssectionName=ssection)
        if (ctype .ne. COLLCT_SCAVECTOR) then
          call output_line ("Ignoring unknown type: "//trim(svecname))
          cycle
        end if
        
        p_rvectorScalar => collct_getvalue_vecsca (rcmdStatus%rcollection, svecname, &
            ssectionName=ssection)

        ! Write the block
        call spdp_projectToCells (p_rvectorScalar,p_Ddata1)
        call ucd_addVariableElementBased (rexport,trim(sname),UCD_VAR_STANDARD,p_Ddata1)
        deallocate(p_Ddata1)

      end if
      
    end do
    
    ! Write the file, done.
    call ucd_write(rexport)
    call ucd_release(rexport)
    
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_daxpyvector (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: daxpyvector.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(4), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_DOUBLE, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_DOUBLE, COLLCT_UNDEFINED) &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    character(len=COLLCT_MLNAME) :: ssection
    character(len=SYS_STRLEN) :: stoken
    integer :: ctype1,ctype2
    type(t_vectorBlock), pointer :: p_rvectorBlock1,p_rvectorBlock2
    type(t_vectorScalar), pointer :: p_rvectorScalar1,p_rvectorScalar2
    real(DP) :: da,dp
    
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  daxpyvector()      - Linear combination of two vectors")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("daxpyvector - Do a linear combination of two vectors.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   daxpyvector ([varsource],[vardest],[da],[dp])")
      call output_lbrk ()
      call output_line ("Does a DAXPY operation in the following form:")
      call output_line ("  [vardest] = [da]*[varsource] + [dp]*[vardest]")
      call output_line ("The vectors may be scalar or block.")
      call output_lbrk ()
      call output_line ("Example:")
      call output_line ("    daxpyvector (source,dest,-1.0,1.0)")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return
    
    
    ! Get the parameters
    
    da = Rvalues(3)%dvalue
    dp = Rvalues(4)%dvalue
    
    stoken = Rvalues(1)%svarname
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
    ctype1 = collct_gettype (rcmdStatus%rcollection, Rvalues(1)%svarname, ssectionName=ssection)

    if (ctype1 .eq. COLLCT_BLKVECTOR) then
    
      p_rvectorBlock1 => collct_getvalue_vec (rcmdStatus%rcollection, stoken,&
            ssectionName=ssection)
    
      stoken = Rvalues(2)%svarname
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
      ctype2 = collct_gettype (rcmdStatus%rcollection, Rvalues(2)%svarname, ssectionName=ssection)
      if (ctype2 .ne. COLLCT_BLKVECTOR) then
        call output_line ("Source and destination not compatible.")
        return
      end if
      
      p_rvectorBlock2 => collct_getvalue_vec (rcmdStatus%rcollection, stoken,&
            ssectionName=ssection)
            
      ! Copy!
      call lsysbl_vectorLinearComb (p_rvectorBlock1,p_rvectorBlock2,da,dp)
      
    else if (ctype1 .eq. COLLCT_SCAVECTOR) then

      p_rvectorScalar1 => collct_getvalue_vecsca (rcmdStatus%rcollection, stoken,&
            ssectionName=ssection)
    
      stoken = Rvalues(2)%svarname
      call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
      ctype2 = collct_gettype (rcmdStatus%rcollection, Rvalues(2)%svarname, ssectionName=ssection)
      if (ctype2 .ne. COLLCT_BLKVECTOR) then
        call output_line ("Source and destination not compatible.")
        return
      end if

      p_rvectorScalar2 => collct_getvalue_vecsca (rcmdStatus%rcollection, stoken,&
            ssectionName=ssection)
            
      ! Copy!
      call lsyssc_vectorLinearComb (p_rvectorScalar1,p_rvectorScalar2,da,dp)
      
    else
      call output_line ("Invalid source vector.")
      return
    end if
    
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_clearvector (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: clearvector.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(2), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("value", STYPE_DOUBLE, COLLCT_UNDEFINED)  &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    character(len=COLLCT_MLNAME) :: ssection
    character(len=SYS_STRLEN) :: stoken
    integer :: ctype1
    type(t_vectorBlock), pointer :: p_rvectorBlock1
    type(t_vectorScalar), pointer :: p_rvectorScalar1
    real(DP) :: da
        
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  clearvector()      - Clears a vector or overwrites with a value.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("clearvector - Clears a vector or overwrites with a numbert.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   clearvector ([vector])")
      call output_line ("   clearvector ([vector],value=[val])")
      call output_lbrk ()
      call output_line ("Overwrites [vector] with zero. If [val] is specified,")
      call output_line ("overwrites the vector with [val].")
      call output_line ("The vectors may be scalar or block.")
      call output_lbrk ()
      call output_line ("Example:")
      call output_line ("    clearvector (source)")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return
    
    ! Get the parameters
    
    da = 0.0_DP
    if (Iindex(2) .ne. 0) da = Rvalues(Iindex(2))%dvalue
    
    stoken = Rvalues(1)%svarname
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
    ctype1 = collct_gettype (rcmdStatus%rcollection, Rvalues(1)%svarname, ssectionName=ssection)

    if (ctype1 .eq. COLLCT_BLKVECTOR) then
    
      p_rvectorBlock1 => collct_getvalue_vec (rcmdStatus%rcollection, stoken,&
            ssectionName=ssection)

      if (da .eq. 0.0_DP) then
        call lsysbl_clearVector (p_rvectorBlock1)
      else
        call lsysbl_clearVector (p_rvectorBlock1,da)
      end if
      
    else if (ctype1 .eq. COLLCT_SCAVECTOR) then

      p_rvectorScalar1 => collct_getvalue_vecsca (rcmdStatus%rcollection, stoken,&
            ssectionName=ssection)
    
      if (da .eq. 0.0_DP) then
        call lsyssc_clearVector (p_rvectorScalar1)
      else
        call lsyssc_clearVector (p_rvectorScalar1,da)
      end if
      
    else
      call output_line ("Invalid vector.")
      return
    end if
    
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine fcmd_scalevector (rcmdStatus,inestlevel,rreturn,cexecmode,Rvalues)
  
  !<description>
    ! Command: clearvector.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
    
    ! Level of nesting
    integer, intent(in) :: inestlevel

    ! Type of execution mode. One of the FCMD_EXECxxxx constants.
    integer, intent(in) :: cexecmode
  !</inputoutput>

  !<input>
    ! OPTIONAL: Command line arguments.
    type(t_symbolValue), dimension(:), intent(in), optional :: Rvalues
  !</input>

  !<output>
    ! Return value
    type(t_symbolValue), intent(inout) :: rreturn
  !</output>
  
  !</subroutine>
  
    ! Command arguments
    type(t_varspec), dimension(2), parameter :: Rvarspec = &
      (/ t_varspec("", STYPE_VAR, COLLCT_UNDEFINED), &
         t_varspec("", STYPE_DOUBLE, COLLCT_UNDEFINED)  &
       /)
    integer, dimension(size(Rvarspec)) :: Iindex
    logical :: bsuccess

    ! local variables
    character(len=COLLCT_MLNAME) :: ssection
    character(len=SYS_STRLEN) :: stoken
    integer :: ctype1
    type(t_vectorBlock), pointer :: p_rvectorBlock1
    type(t_vectorScalar), pointer :: p_rvectorScalar1
    real(DP) :: da
        
    select case (cexecmode)
    case (FCMD_EXECSHORTHELP)
      ! Print a short help message and return
      call output_line ("  scalevector()      - Scales a vector by a value.")
      
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    case (FCMD_EXECLONGHELP)
      ! Print a long help message and return
      call output_line ("scalevector - Scales a vector by a value.")
      call output_lbrk ()
      call output_line ("Usage:")
      call output_line ("   scalevector ([vector],[val])")
      call output_lbrk ()
      call output_line ("Scales the vector [vector] by the value [val].")
      call output_line ("The vectors may be scalar or block.")
      call output_lbrk ()
      call output_line ("Example:")
      call output_line ("    scalevector (source,-1.0)")
    
      ! Ok.
      rreturn%ctype = STYPE_INTEGER
      return
    end select
  
    ! Check parameters
    call fcmd_getparameters (Rvarspec,Iindex,bsuccess,&
        rcmdStatus%rcollection,inestlevel,Rvalues)
    if (.not. bsuccess) return
    
    ! Get the parameters
    
    da = Rvalues(2)%dvalue
    
    stoken = Rvalues(1)%svarname
    call cmdprs_getSymbolSection (rcmdStatus%rcollection,stoken,inestlevel,ssection)
    ctype1 = collct_gettype (rcmdStatus%rcollection, Rvalues(1)%svarname, ssectionName=ssection)

    if (ctype1 .eq. COLLCT_BLKVECTOR) then
    
      p_rvectorBlock1 => collct_getvalue_vec (rcmdStatus%rcollection, stoken,&
            ssectionName=ssection)

      call lsysbl_scaleVector (p_rvectorBlock1,da)
      
    else if (ctype1 .eq. COLLCT_SCAVECTOR) then

      p_rvectorScalar1 => collct_getvalue_vecsca (rcmdStatus%rcollection, stoken,&
            ssectionName=ssection)
    
      call lsyssc_scaleVector (p_rvectorScalar1,da)
      
    else
      call output_line ("Invalid vector.")
      return
    end if
    
    ! Ok.
    rreturn%ctype = STYPE_INTEGER

  end subroutine

end module
