module commandparser

  use ISO_VARYING_STRING
  
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
  
  use vectorio

  implicit none
  
  private
  
  ! Type that encapsules the current status.
  type t_commandstatus
    
    ! Set to TRUE to terminate
    logical :: bterminate = .false.
    
    ! Echo of the command line.
    logical :: becho = .false.
    
    ! Global collection object with all variables
    type(t_collection) :: rcollection
    
  end type
  
  public :: t_commandstatus
  public :: cmdprs_init
  public :: cmdprs_done
  public :: cmdprs_parsestream
  public :: cmdprs_parseterminal
  
contains

  ! ***************************************************************************

  subroutine cmdprs_init (rcmdStatus)

  !<description>
    ! Initialises a command line.
  !</description>
  
  !<output>
    ! Command line to initialise.
    type(t_commandstatus), intent(out) :: rcmdStatus
  !</output>
  
    call collct_init (rcmdStatus%rcollection)

  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_done (rcmdStatus)

  !<description>
    ! Clean up a command line.
  !</description>

  !<inputoutput>
    ! Command line to clean up.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>
  
    call collct_done (rcmdStatus%rcollection)
  
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_releaseargs (p_Sargs)
  
  !<description>
    ! Releases memory of a command line array.
  !</description>
  
  !<inputoutput>
    ! Command line to release.
    type(VARYING_STRING), dimension(:), pointer :: p_Sargs
  !</inputoutput>

    integer :: i
    
    do i=1,size(p_Sargs)
      p_Sargs(i) = ""
    end do
    
    deallocate (p_Sargs)

  end subroutine

  !************************************************************************

  elemental function cmdprs_dequote (ssource)
  !<description>
    ! Remove all quotation marks from a string.
  !</description>
  
  !<input>
    ! A source string.
    type(VARYING_STRING), intent(in) :: ssource
  !</input>
  
  !<output>
    ! Destination string. Non-escaped quotation marks are removed.
    type(VARYING_STRING) :: cmdprs_dequote
  !</output>
  
    type(VARYING_STRING) :: stemp,sdest
    integer :: i, idest
    character :: squotechar, scurrentchar
    logical :: bescape
    
    ! Copy the string, remove leading/trailing spaces.
    stemp = trim(ssource)
    
    ! Initialise destination
    sdest = stemp
    
    ! Loop through the string. Copy characters.
    ! Replace quotation marks.
    i = 1
    idest = 0
    bescape = .false.
    squotechar = ' '
    do 
      if (i .gt. len(stemp)) exit
      
      ! Get the current character.
      scurrentchar = getchar (stemp,i);
      
      if (bescape) then
      
        ! Switch off escape mode.
        bescape = .false.
        
      else
        
        if (scurrentchar .eq. "\") then
        
          ! take the next character as it is.
          bescape = .true.
          
          ! Ignore the escape char.
          i = i+1
          cycle
        
        else if ((scurrentchar .eq. "'") .or. (scurrentchar .eq. """")) then
          
          ! Start quoting or stop it. Only the current quote char
          ! triggers off the quoting.
          if (squotechar .eq. scurrentchar) then
            squotechar = " "
          else
            squotechar = scurrentchar
          end if
          
          ! Ignore the quotation char.
          i = i+1
          cycle
          
        end if
      end if

      ! Copy the character and go to the next char.
      idest = idest+1
      call setchar (sdest,idest,scurrentchar)

      i = i+1
    end do
    
    ! From the destination, take the actual substring.
    cmdprs_dequote = extract (sdest,1,idest)
    
    ! Release memory.
    sdest = ""
    stemp = ""
  
  end function 

  ! ***************************************************************************

  subroutine cmdprs_getparam (Sargs,iparam,sparam,bupcase,bdequote)
  
  !<description>
    ! Returns parameter iparam from a list of parameters.
  !</description>
  
  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
    
    ! Number of the argument.
    integer, intent(in) :: iparam
    
    ! Whether to upper-case the parameter or nor
    logical, intent(in) :: bupcase

    ! Whether to dequote the parameter or not.
    logical, intent(in) :: bdequote
  !</input>
  
  !<output>
    ! Returns the iparam'sth parameter or "" if that one does
    ! not exist.
    character(len=*), intent(out) :: sparam
  !</output>
  
    ! local variables
  
    if (iparam .gt. size(Sargs)) then
      ! Does not exist.
      sparam = ""
      return
    end if
    
    ! Get the parameter
    sparam = adjustl(char(Sargs(iparam)))

    if (bupcase) then
      call sys_toupper(sparam)
    end if
    
    if (bdequote) then
      call sys_dequote(sparam)
    end if

  end subroutine
  
  ! ***************************************************************************

  subroutine cmdprs_splitline (ssource,p_Sargs)
  
  !<description>
    ! Subduvides a line into an array of substrings.
    ! Quotation marks around substrings are removed.
  !</description>
  
  !<input>
    ! A source string.
    type(VARYING_STRING), intent(in) :: ssource
  !</input>
  
  !<output>
    ! A pointer pointing to an array of strings, each string referring to
    ! one parameter in ssource.
    ! The value NULL is returned if there are no arguments
    ! in the string.
    type(VARYING_STRING), dimension(:), pointer :: p_Sargs
  !</output>
  
    type(VARYING_STRING) :: stemp
    integer :: i, icount, istart
    character :: squotechar, scurrentchar
    logical :: bescape
    
    
    ! Copy the string, remove leading/trailing spaces.
    stemp = trim(ssource)
    
    ! Loop through the string. Replace spaces by #0 except
    ! for those being quoted.
    i = 1
    bescape = .false.
    squotechar = ' '
    do 
      if (i .gt. len(stemp)) exit
      
      if (bescape) then
      
        ! Switch off escape mode.
        bescape = .false.

      else
        ! Get the current character.
        scurrentchar = getchar (stemp,i);
        
        if (scurrentchar .eq. "\") then
        
          ! take the next character as it is.
          bescape = .true.
        
        else if ((scurrentchar .eq. "'") .or. (scurrentchar .eq. """")) then
          
          ! Start quoting or stop it. Only the current quote char
          ! triggers off the quoting.
          if (squotechar .eq. scurrentchar) then
            squotechar = " "
          else
            squotechar = scurrentchar
          end if
          
        else if ((squotechar .eq. " ") .and. (scurrentchar .eq. " ")) then
        
          ! Replace by #0. That's the separator.
          call setchar (stemp,i,char(0))
          
        end if
      end if

      ! next char
      i = i+1
    end do
    
    ! Count number of words.#
    icount = 1
    do i=1,len(stemp)
      if (getchar(stemp,i) .eq. char(0)) then
        icount = icount + 1
      end if
    end do
    
    ! Allocate memory
    allocate (p_Sargs(icount))

    ! Get the parameters
    istart = 1
    icount = 0
    do i=1,len(stemp)
      if (getchar(stemp,i) .eq. char(0)) then
        icount = icount + 1
        p_Sargs(icount) = extract(stemp,istart,i-1)
        
        ! Continue behind current position
        istart = i+1
      end if
    end do
    
    ! Last argument
    icount = icount + 1
    p_Sargs(icount) = extract(stemp,istart,len(stemp))
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_parsestream (rcmdStatus,istream)
  
  !<description>
    ! Parses an input stream.
  !</description>
  
  !<inputoutput>
    ! Command line status object
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! The input stream channel.
    integer, intent(in) :: istream
  !</input>

    type(VARYING_STRING) :: sinput,sin,sin2
    type(VARYING_STRING), dimension(:), pointer :: p_Sargs
    integer :: ierr

    do while (.not. rcmdStatus%bterminate)
    
      sin = ""
      sin2 = ""
      sinput = ""
      
      call get(istream,sin,IOSTAT=ierr) ! read next line of file
      
      ! Trim. Warning, this allocates memory!!!
      sin2 = adjustl(sin)
      sinput = trim(sin2)
      
      
      if (rcmdStatus%becho) then
        call output_line (CHAR(sinput))
      end if
      
      if (ierr == -1 .or. ierr > 0) then
        
        ! Stop, there's no more data.
        exit
        
      else
        if (len(sinput) .ne. 0) then
          ! Ignore comments
          if (getchar (sinput,1) .ne. "#") then
      
            ! Parse the line.
            call cmdprs_splitline (sinput,p_Sargs)
            
            ! A command involved?
            if (ubound(p_Sargs,1) .ne. 0) then
            
              ! Execute the command
              call cmdprs_docommand (rcmdStatus,p_Sargs)
            
            end if
            
            ! Release memory
            call cmdprs_releaseargs(p_Sargs)
          end if
        end if
      end if
    end do

    ! Release memory    
    sin = ""
    sin2 = ""
    sinput = ""
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_parseterminal (rcmdStatus)
  
  !<description>
    ! Parses an the terminal.
  !</description>

  !<inputoutput>
    ! Command line status object
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

    type(VARYING_STRING) :: sinput
    type(VARYING_STRING), dimension(:), pointer :: p_Sargs
    character(len=SYS_STRLEN) :: sstr

    do while (.not. rcmdStatus%bterminate)
      
      ! Read tghe next line from the termninal
      call output_line("> ",bnolinebreak=.true.,bnotrim=.true.)
      read (IN_TERMINAL,'(A)') sstr
      
      sinput = trim(adjustl(sstr))
      
      if (rcmdStatus%becho) then
        call output_line (CHAR(sinput))
      end if
      
      if (len(sinput) .ne. 0) then
      
        ! Ignore comments
        if (getchar (sinput,1) .ne. "#") then

          ! Parse the line.
          call cmdprs_splitline (sinput,p_Sargs)
          
          ! A command involved?
          if (associated(p_Sargs)) then
          
            ! Execute the command
            call cmdprs_docommand (rcmdStatus,p_Sargs)
          
          end if
          
          ! Release memory
          call cmdprs_releaseargs(p_Sargs)
          
        end if
      
      end if
      
      call output_lbrk()
     
    end do
    
    ! Release memory.
    sinput = ""
  
  end subroutine

  ! ***************************************************************************

  recursive subroutine cmdprs_docommand (rcmdStatus,Sargs)
  
  !<description>
    ! Executes a command.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus

    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</inputoutput>

    type(VARYING_STRING) :: scmd, sarg
    
    ! Get the command
    scmd = Sargs(1)
    call TOUPPER(scmd)
    
    if (scmd .eq. "EXIT") then
      rcmdStatus%bterminate = .true.
      return
    end if
    
    if (scmd .eq. "PRINT") then
      if (size (Sargs) .gt. 1) then
        sarg = cmdprs_dequote(Sargs(2))
        call output_line (char(sarg))
        sarg = ""
      else
        call output_lbrk ()
      end if
      return
    end if

    if (scmd .eq. "MEMINFO") then
      call output_lbrk ()
      call storage_info(.true.)
      return
    end if

    if (scmd .eq. "SET") then
      call cmdprs_do_set (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "ECHO") then
      if (size (Sargs) .gt. 1) then
        if (sys_upcase(char(Sargs(2))) .eq. "ON") then
          rcmdStatus%becho = .true.
        else if (sys_upcase(char(Sargs(2))) .eq. "OFF") then
          rcmdStatus%becho = .false.
        end if
      else
        if (rcmdStatus%becho) then
          call output_line ("Echo is on.")
        else
          call output_line ("Echo is off.")
        end if
      end if
      return
    end if

    if (scmd .eq. "EXEC") then
      call cmdprs_do_exec (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "SHOW") then
      call cmdprs_do_show (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "INFO") then
      call cmdprs_do_info (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "DELETE") then
      call cmdprs_do_delete (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "DESTROY") then
      call cmdprs_do_destroy (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "HELP") then
      call cmdprs_do_help (Sargs)
      return
    end if

    if (scmd .eq. "READ2DPRM") then
      call cmdprs_do_read2dprm (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "READ2DTRI") then
      call cmdprs_do_read2dtri (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "MESHREFINE") then
      call cmdprs_do_meshrefine (rcmdStatus,Sargs)
      return
    end if
    
    if (scmd .eq. "MESHHIERARCHY") then
      call cmdprs_do_meshhierarchy (rcmdStatus,Sargs)
      return
    end if
    
    if (scmd .eq. "FESPACE") then
      call cmdprs_do_fespace (rcmdStatus,Sargs)
      return
    end if
    
    if (scmd .eq. "FEHIERARCHY") then
      call cmdprs_do_fehierarchy (rcmdStatus,Sargs)
      return
    end if
    
    if (scmd .eq. "MLEVELPRJHIERARCHY") then
      call cmdprs_do_mlevelprjhierarchy (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "READBLOCKVECTOR") then
      call cmdprs_do_readblockvector (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "WRITEBLOCKVECTOR") then
      call cmdprs_do_writeblockvector (rcmdStatus,Sargs)
      return
    end if
    
    if (scmd .eq. "CREATEBLOCKVECTOR") then
      call cmdprs_do_createblockvector (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "COPYSUBVECTOR") then
      call cmdprs_do_copysubvector (rcmdStatus,Sargs)
      return
    end if
    
    if (scmd .eq. "COPYVECTOR") then
      call cmdprs_do_copyvector (rcmdStatus,Sargs)
      return
    end if

    if (scmd .eq. "INTERPOLATEVECTOR") then
      call cmdprs_do_interpolatevector (rcmdStatus,Sargs)
      return
    end if
    
    if (scmd .eq. "L2PROJECTION") then
      call cmdprs_do_l2projection (rcmdStatus,Sargs)
      return
    end if
    
    if (scmd .ne. "") then
      call output_line ("Unknown command!")
    end if
  
  end subroutine

  ! ***************************************************************************
  ! Extended command handlers
  ! ***************************************************************************

  ! ***************************************************************************

  subroutine cmdprs_do_help (Sargs)
  
  !<description>
    ! Command: HELP.
  !</description>
  
  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    character(len=20) :: sargformatted
  
    call output_lbrk ()

    if (size(Sargs) .eq. 1) then
      call output_line ("FEAT2 command parser.")
      call output_line ("---------------------")
      call output_line ("The following commands are available. For a specific help")
      call output_line ("type ""HELP [command]"".")
      call output_lbrk ()
      call output_line ("  help               - This help page.");
      call output_line ("  exit               - Exits the command line.");
      call output_line ("  meminfo            - Information about memory management.");
      call output_line ("  exec               - Execute script.");
      call output_line ("  print              - Print some text.");
      call output_line ("  set                - Modify/set/print environment.");
      call output_line ("  show               - Show environment variable (if possible).");
      call output_line ("  delete             - Delete variable from environment.");
      call output_line ("  destroy            - Destroys an environment structure (extended delete).");
      call output_line ("  read2dprm          - Read 2D .prm file.");
      call output_line ("  read2dtri          - Read 2D .tri file.");
      call output_line ("  meshrefine         - Refine a mesh.");
      call output_line ("  meshhierarchy      - Create a mesh hierarchy.");
      call output_line ("  fespace            - Create a FE space.");
      call output_line ("  fehierarchy        - Create a FE hierarchy.");
      call output_line ("  mlevelprjhierarchy - Create a multilevel projection hierarchy.");
      call output_line ("  createblockvector  - Create an empty block vector.");
      call output_line ("  readblockvector    - Read a block vector from a file.");
      call output_line ("  writeblockvector   - Write a block vector to a file.");
      call output_line ("  copyvector         - Copy a vector to another.");
      call output_line ("  copysubvector      - Copy a subvector to another.");
      call output_line ("  interpolatevector  - Interpolate a vector to another level.");
      call output_line ("  l2projection       - Appies an L2 projection to a vector.");
    
    else
    
      call cmdprs_getparam (Sargs,2,sargformatted,.true.,.false.)
      
      if (sargformatted .eq. "EXIT") then
        call output_line ("EXIT - Close interactive command line, return to terminal.");
        call output_lbrk ()
        call output_line ("Usage:");
        call output_line ("  EXIT");
      
        return
      end if

      if (sargformatted .eq. "EXEC") then
        call output_line ("EXEC - Execute a script on hard disc.");
        call output_lbrk ()
        call output_line ("Usage:");
        call output_line ("  EXEC ""[filename/path]""");
      
        return
      end if

      if (sargformatted .eq. "PRINT") then
        call output_line ("PRINT - prints some text to the terminal.");
        call output_lbrk ()
        call output_line ("Usage:");
        call output_line ("  PRINT ""[some text]""");
      
        return
      end if

      if (sargformatted .eq. "SET") then
        call output_line ("SET - Modify or show environment.");
        call output_lbrk ()
        call output_line ("Usage:");
        call output_line ("* To show statistics about the environment:");
        call output_line ("      SET");
        call output_lbrk ()
        call output_line ("*To set/modify an integer variable:");
        call output_line ("      SET INT [name] [value]");
        call output_lbrk ()
        call output_line ("*To set/modify an double precision variable:");
        call output_line ("      SET DOUBLE [name] [value]");
        call output_lbrk ()
        call output_line ("*To set/modify a string variable:");
        call output_line ("      SET STRING [name] [value]");
      
        return
      end if
    
      if (sargformatted .eq. "DELETE") then
        call output_line ("DELETE - Delete a variable from the environment.");
        call output_lbrk ()
        call output_line ("Usage:");
        call output_line ("  DELETE [name]");
        call output_lbrk ()
        call output_line ("WARNING: This does not release probably associated memory!");
      
        return
      end if
    
      if (sargformatted .eq. "DESTROY") then
        call output_line ("DESTROY - Releases memory associated to an environment variable.");
        call output_lbrk ()
        call output_line ("Usage:");
        call output_line ("  DESTROY [name]");
        call output_lbrk ()
        call output_line ("Automatically detects the type of an environment variable and");
        call output_line ("releases the associated structure from memory (if there is one).");
        call output_line ("The variable is removed from the environment");
      
        return
      end if
    
      if (sargformatted .eq. "SHOW") then
        call output_line ("SHOW - Show content about an environment vairbale.");
        call output_lbrk ()
        call output_line ("Usage:");
        call output_line ("  SHOW [name]");
        call output_lbrk ()
        call output_line ("If possible, this routine prints information about an");
        call output_line ("environment variable to the terminal. For standard types");
        call output_line ("(int, double, string,...), the value is returned.");
        call output_line ("For extended types (e.g. meshes), status information");
        call output_line ("is returned.");
      
        return
      end if
    
      if (sargformatted .eq. "INFO") then
        call output_line ("INFO - Detailed information about environment vairbales.");
        call output_lbrk ()
        call output_line ("Usage:");
        call output_line ("  SHOW [name]");
        call output_lbrk ()
        call output_line ("Shows detailed information about an environment variable.");
        call output_line ("This includes type, value,...");
      
        return
      end if
    
      if (sargformatted .eq. "READ2DPRM") then
        call output_line ("READ2DPRM - Read 2D .PRM file.");
        call output_lbrk ()
        call output_line ("Usage:");
        call output_line ("  READ2DPRM [varname] [filename]");
        call output_lbrk ()
        call output_line ("Read in a .prm file and store it using the name [variable].")
      
        return
      end if

      if (sargformatted .eq. "READ2DTRI") then
        call output_line ("READ2DTRI - Read 2D .TRI file.");
        call output_lbrk ()
        call output_line ("Usage:");
        call output_line ("  READ2DTRI [varname] [filename] [ BOUNDARY varbd ]");
        call output_lbrk ()
        call output_line ("Read in a .tri file and store it using the name [variable].")
        call output_line ("If BOUNDARY is specified, the triangulation is connected")
        call output_line ("to a boundary object varbd.")
        call output_lbrk ()
        call output_line ("Example:")
        call output_line ("  read2dprm myprm ""myprmfile.prm""");
        call output_line ("  read2dtri mytri ""mytrifile.tri"" boundary myprm");
      
        return
      end if

      if (sargformatted .eq. "MESHREFINE") then
        call output_line ("MESHREFINE - Refine a mesh.");
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   meshrefine [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Refine a mesh with a given method one or multiple times.")
        call output_line ("[varname] identifies the variable to create.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --TIMES [varname] ...")
        call output_line ("      Refine the mesh [varname] times. Default is ""--TIMES 1"".")
        call output_lbrk ()
        call output_line ("  ... --BOUNDARY [varname] ...")
        call output_line ("      Use the specified boundary object [varname] fdor boundary refinement.")
        call output_line ("      If not specified, no boundary object is used.")
        call output_lbrk ()
        call output_line ("  ... --METHOD [varname] ...")
        call output_line ("      Use a specific method. Default is ""--METHOD 2LEVELORDERED"". Possible choices:")
        call output_line ("      2LEVELORDERED   - Use 2-level-ordering refinement.")
      
        return
      end if

      if (sargformatted .eq. "MESHHIERARCHY") then
        call output_line ("MESHHIERARCHY - Create a mesh hierarchy.");
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   meshhierarchy [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Refine a mesh with a given method one or multiple times.")
        call output_line ("[varname] identifies the variable to create.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --MESH [varname] ...")
        call output_line ("      Use mesh [varname] as coarse mesh of teh hierarchy.")
        call output_line ("      Mandatory argument")
        call output_lbrk ()
        call output_line ("  ... --BOUNDARY [varname] ...")
        call output_line ("      Use boundary definition [varname] for refinements.")
        call output_lbrk ()
        call output_line ("  ... --LEVELS [varname] ...")
        call output_line ("      Creates a hierarchy of [varname] levels. Default is ""--LEVELS 1"".")
        call output_lbrk ()
        call output_line ("  ... --BOUNDARY [varname] ...")
        call output_line ("      Use the specified boundary object [varname] fdor boundary refinement.")
        call output_line ("      If not specified, no boundary object is used.")
        call output_lbrk ()
        call output_line ("  ... --METHOD [varname] ...")
        call output_line ("      Use a specific method. Default is ""--METHOD 2LEVELORDERED"". Possible choices:")
        call output_line ("      2LEVELORDERED   - Use 2-level-ordering refinement.")
      
        return
      end if

      if (sargformatted .eq. "FESPACE") then
        call output_line ("FESPACE - Create an FE space.");
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   fespace [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Creates an FE space that can be used to set up matrices/vectors.")
        call output_line ("[varname] identifies the variable to create.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --MESH [varname] ...")
        call output_line ("      Use mesh [varname] as coarse mesh of the hierarchy.")
        call output_line ("      Mandatory argument")
        call output_lbrk ()
        call output_line ("  ... --BOUNDARY [varname] ...")
        call output_line ("      Use boundary definition [varname] for refinements.")
        call output_lbrk ()
        call output_line ("  ... --COMPONENTS [varname] ...")
        call output_line ("      Creates a FE space with [varname] components. Mandatory argument.")
        call output_line ("      Must be defined in advance to all element/cubature specifications.")
        call output_lbrk ()
        call output_line ("  ... --TRIELEMENTS EL_xxx EL_xxx EL_xxx ...")
        call output_line ("      Specifies a list of TRI/TETRA elemnent types to be used for");
        call output_line ("      triangular/tetrahedral elements. There must be one ID per component.")
        call output_lbrk ()
        call output_line ("  ... --QUADELEMENTS EL_xxx EL_xxx EL_xxx ...")
        call output_line ("      Specifies a list of QUAD/HEXA elemnent types to be used for");
        call output_line ("      quadrilateral/hexahedral elements. There must be one ID per component.")
        call output_lbrk ()
        call output_line ("  ... --TRICUB CUB_xxx CUB_xxx CUB_xxx ...")
        call output_line ("      Specifies a list of TRI/TETRA cubature formulas to be used for");
        call output_line ("      triangular/tetrahedral elements. There must be one ID per component.")
        call output_line ("      If not defined, the default cubature formula is used.")
        call output_lbrk ()
        call output_line ("  ... --QUADCUB CUB_xxx CUB_xxx CUB_xxx ...")
        call output_line ("      Specifies a list of QUAD/HEXA cubature formulas to be used for");
        call output_line ("      quadrilateral/hexahedral elements. There must be one ID per component.")
        call output_line ("      If not defined, the default cubature formula is used.")
      
        return
      end if

      if (sargformatted .eq. "FEHIERARCHY") then
        call output_line ("FEHIERARCHY - Create an FE hierarchy.");
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   fehierarchy [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Creates an FE space that can be used to set up matrices/vectors.")
        call output_line ("[varname] identifies the variable to create.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --MESHHIERARCHY [varname] ...")
        call output_line ("      Use mesh hierarchy [varname] as coarse mesh of the hierarchy.")
        call output_line ("      Mandatory argument")
        call output_lbrk ()
        call output_line ("  ... --BOUNDARY [varname] ...")
        call output_line ("      Use boundary definition [varname] for refinements.")
        call output_lbrk ()
        call output_line ("  ... --COMPONENTS [varname] ...")
        call output_line ("      Creates a FE space with [varname] components. Mandatory argument.")
        call output_line ("      Must be defined in advance to all element/cubature specifications.")
        call output_lbrk ()
        call output_line ("  ... --TRIELEMENT EL_xxx EL_xxx EL_xxx ...")
        call output_line ("      Specifies a list of TRI/TETRA elemnent types to be used for");
        call output_line ("      triangular/tetrahedral elements. There must be one ID per component.")
        call output_lbrk ()
        call output_line ("  ... --QUADELEMENT EL_xxx EL_xxx EL_xxx ...")
        call output_line ("      Specifies a list of QUAD/HEXA elemnent types to be used for");
        call output_line ("      quadrilateral/hexahedral elements. There must be one ID per component.")
        call output_lbrk ()
        call output_line ("  ... --TRICUB CUB_xxx CUB_xxx CUB_xxx ...")
        call output_line ("      Specifies a list of TRI/TETRA cubature formulas to be used for");
        call output_line ("      triangular/tetrahedral elements. There must be one ID per component.")
        call output_lbrk ()
        call output_line ("  ... --QUADCUB CUB_xxx CUB_xxx CUB_xxx ...")
        call output_line ("      Specifies a list of QUAD/HEXA cubature formulas to be used for");
        call output_line ("      quadrilateral/hexahedral elements. There must be one ID per component.")
      
        return
      end if

      if (sargformatted .eq. "MLEVELPRJHIERARCHY") then
        call output_line ("MLEVELPRJHIERARCHY - Create an multilevel projection hierarchy.");
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   mlevelprjhierarchy [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Creates a multilevel projection hierarchy that can be used to transfer")
        call output_line ("soltion/rhs vectors from one level to another.")
        call output_line ("[varname] identifies the variable to create.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --FEHIERARCHY [varname] ...")
        call output_line ("      Create a multilevel projection hierarchy based on the");
        call output_line ("      FE space hierarchy [varname]. Mandatory argument");
      
        return
      end if

      if (sargformatted .eq. "CREATEBLOCKVECTOR") then
        call output_line ("CREATEBLOCKVECTOR - Create an empty block vector.");
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   createblockvector [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Creates an empty block vector..")
        call output_line ("[varname] identifies the variable to create.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --FEHIERARCHY [varname] ...")
        call output_line ("      Defines the FE hierarchy, the vector is based on.");
        call output_line ("      Mandatory argument");
        call output_lbrk ()
        call output_line ("  ... --LEVEL [varname] ...")
        call output_line ("      Level in the hierarchy specifying the discretisation of the vector.");
        call output_line ("      Default is max. level.");
      
        return
      end if

      if (sargformatted .eq. "READBLOCKVECTOR") then
        call output_line ("READBLOCKVECTOR - Read block vector from file.");
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   readblockvector [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Read a block vector from a file.")
        call output_line ("[varname] identifies the variable where to read data to.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --FILENAME [varname] ...")
        call output_line ("      Defines the filename of the file to read.");
        call output_line ("      Mandatory argument");
        call output_lbrk ()
        call output_line ("  ... --UNFORMATTED ...")
        call output_line ("      Read a binary file. Default is formatted, human readable file.");
      
        return
      end if

      if (sargformatted .eq. "WRITEBLOCKVECTOR") then
        call output_line ("WRITEBLOCKVECTOR - Write block vector to file.");
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   writeblockvector [varname] [...options...]")
        call output_lbrk ()
        call output_line ("Write a block vector to a file.")
        call output_line ("[varname] identifies the vector.")
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --FILENAME [varname] ...")
        call output_line ("      Defines the filename of the file to read.");
        call output_line ("      Mandatory argument");
        call output_lbrk ()
        call output_line ("  ... --FORMAT [varname] ...")
        call output_line ("      Defines the format of the numbers in the file.");
        call output_line ("      Applies only if --UNFORMATTED is not specified.");
        call output_line ("      Default: E20.10");
        call output_lbrk ()
        call output_line ("  ... --UNFORMATTED ...")
        call output_line ("      Read a binary file. Default is formatted, human readable file.");
      
        return
      end if

      if (sargformatted .eq. "COPYVECTOR") then
        call output_line ("COPYVECTOR - Copy a vector.");
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   copysubvector [varsource] [vardest]")
        call output_lbrk ()
        call output_line ("Copies vector [varsource] to vector [vardest].")
        call output_lbrk ()
        call output_line ("Example:")
        call output_line ("    copyvector source dest")
      
        return
      end if

      if (sargformatted .eq. "INTERPOLATEVECTOR") then
        call output_line ("INTERPOLATEVECTOR - Interpolate a vector.");
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   interpolatevector [varsource] [lvlsource] [vardest] [lvldest] [...options...]")
        call output_lbrk ()
        call output_line ("Interpolates vector [varsource] from level [lvlsource] to")
        call output_line ("level [lvldest] and writes the result to [vardest].")
        call output_line ("[mlprj] specifies the projection hierarchy to be used.")
        call output_line ("The method uses prolongation/interpolation for the level change.")
        call output_lbrk ()
        call output_line ("Example:")
        call output_line ("    interpolatevector source 3 dest 5 --MLPRJHIERARCHY mhier --FEHIERARCHY feh")
        call output_lbrk ()
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --MLPRJHIERARCHY [varname] ...")
        call output_line ("      Defines the underlying projection hierarchy.");
        call output_line ("      Mandatory argument.");
        call output_lbrk ()
        call output_line ("  ... --FEHIERARCHY [varname] ...")
        call output_line ("      Defines the underlying FE hierarchy.");
        call output_line ("      Mandatory argument if source and destination level differ");
        call output_line ("      by more than one level.");
      
        return
      end if

      if (sargformatted .eq. "L2PROJECTION") then
        call output_line ("L2PROJECTION - Applies an L2 projection to a vector.");
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   l2projection [varsource] [vardest] [...options...]")
        call output_lbrk ()
        call output_line ("Interpolates vector [varsource] to [vardest] using an")
        call output_line ("L2-projection.")
        call output_lbrk ()
        call output_line ("Example:")
        call output_line ("    interpolatevector source dest")
        call output_lbrk ()
        call output_line ("The following options are possible in [...options...]:")
        call output_lbrk ()
        call output_line ("  ... --VERBOSE ...")
        call output_line ("      Activate verbose output.")
        call output_lbrk ()
        call output_line ("  ... --RELERROR [error} ...")
        call output_line ("      Defines the relative accuracy of the projection.")
      
        return
      end if

      if (sargformatted .eq. "COPYSUBVECTOR") then
        call output_line ("COPYSUBVECTOR - Copy a subvector.");
        call output_lbrk ()
        call output_line ("Usage:")
        call output_line ("   copysubvector [varsource] [sourcecomp] [vardest] [destcomp]")
        call output_lbrk ()
        call output_line ("Copies subvector [sourcecomp] of vector [varsource] to")
        call output_line ("subvector [destcomp] of vector [vardest].")
        call output_lbrk ()
        call output_line ("Example:")
        call output_line ("    copysubvector source 1 dest 4")
      
        return
      end if

      call output_line ("No specific help available.")
    end if
  
  end subroutine

  ! ***************************************************************************

  recursive subroutine cmdprs_do_exec (rcmdStatus,Sargs)
  
  !<description>
    ! Command: EXEC.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    integer :: iunit
    logical :: bexists
    character(len=SYS_STRLEN) :: svalue
  
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the script name
    call cmdprs_getparam (Sargs,2,svalue,.false.,.true.)

    ! Open the file and read.
    inquire(file=trim(svalue), exist=bexists)
    if (.not. bexists) then
      call output_line ("File not found!")
    else
      
      ! Open
      call io_openFileForReading(svalue, iunit, .true.)

      ! Interpret the stream
      call cmdprs_parsestream (rcmdStatus,iunit)
      
      close (iunit)
    end if

  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_set (rcmdStatus,Sargs)
  
  !<description>
    ! Command: SET.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    integer :: iarg
    real(DP) :: darg
    character(len=COLLCT_MLNAME) :: svalname
    character(len=20) :: sformat
    character(len=SYS_STRLEN) :: svalue
  
    if (size(Sargs) .eq. 1) then
      call output_lbrk ()
      call output_line ("Environment statistics:")
      call output_lbrk ()
      ! Print the current status of the collection.
      call collct_printStatistics (rcmdStatus%rcollection)
      
    else
    
      if (size(Sargs) .lt. 4) then
        call output_line ("Not enough arguments.")
        return
      end if
      
      ! Get the next token(s)
      call cmdprs_getparam (Sargs,2,sformat,.true.,.false.)
      call cmdprs_getparam (Sargs,3,svalname,.false.,.true.)
      
      if (sformat .eq. "INT") then
        call cmdprs_getparam (Sargs,4,svalue,.false.,.true.)
        read (svalue,*) iarg
        call collct_setvalue_int (rcmdStatus%rcollection, &
            svalname, iarg, .true.)
        return
      end if
    
      if (sformat .eq. "DOUBLE") then
        call cmdprs_getparam (Sargs,4,svalue,.false.,.true.)
        read (svalue,*) darg
        call collct_setvalue_real (rcmdStatus%rcollection, &
            svalname, darg, .true.)
        return
      end if

      if (sformat .eq. "STRING") then
        call cmdprs_getparam (Sargs,4,svalue,.true.,.false.)
        call collct_setvalue_string (rcmdStatus%rcollection, &
            svalname, trim(svalue), .true.)
        return
      end if
    
    end if
  
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_delete (rcmdStatus,Sargs)
  
  !<description>
    ! Command: SET.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    character(len=COLLCT_MLNAME) :: svalname
  
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the next token(s)
    call cmdprs_getparam (Sargs,2,svalname,.true.,.false.)
    
    call collct_deletevalue (rcmdStatus%rcollection, svalname)
    
  end subroutine

  ! ***************************************************************************

  subroutine do_destroy (rcmdStatus,svalname,bverbose)
  
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
    
    ! Print vebose messages
    logical, intent(in) :: bverbose
  !</input>
  
    ! local variables
    integer :: ctype
    logical :: bexists
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_meshhierarchy), pointer :: p_rmeshhierarchy
    type(t_feSpaceLevel), pointer :: p_rfeSpace
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_interlevelProjectionHier), pointer :: p_rinterlevelProjectionHier
    type(t_vectorBlock), pointer :: p_rvectorBlock
  
    ! That's an awful task. Figure out the type at first.
    ! Then get the value and/or print information.
    ctype = collct_gettype (rcmdStatus%rcollection, svalname, bexists=bexists)
    
    if (.not. bexists) then
      if (bverbose) then
        call output_line ("Unknown environment variable!")
      end if
    else
      select case (ctype)
      case (COLLCT_INTEGER,COLLCT_REAL,COLLCT_STRING)
        ! Primitive type
        call collct_deletevalue (rcmdStatus%rcollection, svalname)
          
      case (COLLCT_BOUNDARY)
        ! Boundary object
        p_rboundary => collct_getvalue_bdry (rcmdStatus%rcollection, svalname)
        call boundary_release(p_rboundary)
        deallocate(p_rboundary)

        call collct_deletevalue (rcmdStatus%rcollection, svalname)
          
      case (COLLCT_TRIA)
        ! Boundary object
        p_rtriangulation => collct_getvalue_tria (rcmdStatus%rcollection, svalname)
        call tria_done(p_rtriangulation)
        deallocate(p_rtriangulation)

        call collct_deletevalue (rcmdStatus%rcollection, svalname)

      case (COLLCT_MSHHIERARCHY)
        ! Boundary object
        p_rmeshhierarchy => collct_getvalue_mshh (rcmdStatus%rcollection, svalname)
        call mshh_releaseHierarchy(p_rmeshhierarchy)
        deallocate(p_rmeshhierarchy)

        call collct_deletevalue (rcmdStatus%rcollection, svalname)

      case (COLLCT_FESPACE)
        ! Boundary object
        p_rfeSpace => collct_getvalue_fesp (rcmdStatus%rcollection, svalname)
        call fesph_releaseFEspace(p_rfeSpace)
        deallocate(p_rfeSpace)

        call collct_deletevalue (rcmdStatus%rcollection, svalname)

      case (COLLCT_FEHIERARCHY)
        ! Boundary object
        p_rfeHierarchy => collct_getvalue_feh (rcmdStatus%rcollection, svalname)
        call fesph_releaseHierarchy(p_rfeHierarchy)
        deallocate(p_rfeHierarchy)

        call collct_deletevalue (rcmdStatus%rcollection, svalname)

      case (COLLCT_MLPRJHIERARCHY)
        ! Multilevel projection hierarchy
        p_rinterlevelProjectionHier => collct_getvalue_mlprjh (rcmdStatus%rcollection, svalname)
        call mlprj_releasePrjHierarchy(p_rinterlevelProjectionHier)
        deallocate(p_rinterlevelProjectionHier)

        call collct_deletevalue (rcmdStatus%rcollection, svalname)

      case (COLLCT_BLKVECTOR)
        ! Multilevel projection hierarchy
        p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, svalname)
        call lsysbl_releaseVector(p_rvectorBlock)
        deallocate(p_rvectorBlock)

        call collct_deletevalue (rcmdStatus%rcollection, svalname)

      case default
        call output_line ("Unknown type, variable cannot be destroyed!")         
      end select
      
    end if
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_destroy (rcmdStatus,Sargs)
  
  !<description>
    ! Command: DESTROY.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    character(len=COLLCT_MLNAME) :: svalname
  
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the next token(s)
    call cmdprs_getparam (Sargs,2,svalname,.true.,.false.)
    
    call do_destroy (rcmdStatus,svalname,.true.)
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_show (rcmdStatus,Sargs)
  
  !<description>
    ! Command: SHOW.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    character(len=COLLCT_MLNAME) :: svalname
    character(len=SYS_STRLEN) :: svalue
    integer :: ctype
    logical :: bexists
  
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the next token(s)
    call cmdprs_getparam (Sargs,2,svalname,.true.,.false.)
    
    ! That's an awful task. Figure out the type at first.
    ! Then get the value and/or print information.
    ctype = collct_gettype (rcmdStatus%rcollection, svalname, bexists=bexists)
    
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
        call output_line ("No detailed information available.")
      end select
      
    end if
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_info (rcmdStatus,Sargs)
  
  !<description>
    ! Command: INFO.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    character(len=COLLCT_MLNAME) :: svalname
    character(len=SYS_STRLEN) :: svalue
    integer :: ctype
    logical :: bexists
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_meshhierarchy), pointer :: p_rmeshhierarchy
    type(t_feSpaceLevel), pointer :: p_rfeSpace
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_vectorBlock), pointer :: p_rvectorBlock
  
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the next token(s)
    call cmdprs_getparam (Sargs,2,svalname,.true.,.false.)
    
    ! That's an awful task. Figure out the type at first.
    ! Then get the value and/or print information.
    ctype = collct_gettype (rcmdStatus%rcollection, svalname, bexists=bexists)
    
    if (.not. bexists) then
      call output_line ("Unknown environment variable!")
    else
      call output_line ("Variable: "//svalname)
      select case (ctype)
      case (COLLCT_INTEGER)
        call output_line ("Type    : Integer")
        call output_line ("Content : ",bnolinebreak=.true., bnotrim=.true.)
        call output_line (trim(sys_siL(&
          collct_getvalue_int (rcmdStatus%rcollection, svalname) ,10)))
          
      case (COLLCT_REAL)
        call output_line ("Type    : Double precision")
        call output_line ("Content : ",bnolinebreak=.true., bnotrim=.true.)
        call output_line (trim(sys_sdEL(&
          collct_getvalue_real (rcmdStatus%rcollection, svalname) ,10)))
          
      case (COLLCT_STRING)
        call output_line ("Type    : string")
        call output_line ("Content : ",bnolinebreak=.true., bnotrim=.true.)
        call collct_getvalue_string (rcmdStatus%rcollection, svalname, svalue)
        call output_line (svalue)
          
      case (COLLCT_BOUNDARY)
        call output_line ("Type    : Boundary object (2D)")

      case (COLLCT_TRIA)
        call output_line ("Type    : Triangulation object")
        call output_lbrk ()
        p_rtriangulation => collct_getvalue_tria (rcmdStatus%rcollection, svalname)
        call tria_infoStatistics (p_rtriangulation,.true.)

      case (COLLCT_MSHHIERARCHY)
        call output_line ("Type    : Mesh hierarchy")
        call output_lbrk ()
        p_rmeshhierarchy => collct_getvalue_mshh (rcmdStatus%rcollection, svalname)
        call mshh_printHierStatistics (p_rmeshhierarchy)

      case (COLLCT_FESPACE)
        call output_line ("Type    : FE space")
        call output_lbrk ()
        p_rfeSpace => collct_getvalue_fesp (rcmdStatus%rcollection, svalname)
        call fesph_infoStatistics (p_rfeSpace,.true.)
        call spdiscr_infoBlockDiscr (p_rfeSpace%p_rdiscretisation)

      case (COLLCT_FEHIERARCHY)
        call output_line ("Type    : Mesh hierarchy")
        call output_lbrk ()
        p_rfeHierarchy => collct_getvalue_feh (rcmdStatus%rcollection, svalname)
        call fesph_printHierStatistics (p_rfeHierarchy)

      case (COLLCT_MLPRJHIERARCHY)
        call output_line ("Type    : Multilevel projection hierarchy")
        call output_lbrk ()

      case (COLLCT_BLKVECTOR)
        call output_line ("Type    : Block vector")
        call output_lbrk ()
        p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, svalname)
        call lsysbl_infoVector (p_rvectorBlock)

      case default
        call output_line ("Type    : "//sys_siL(ctype,10))
      end select
      
    end if
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_read2dprm (rcmdStatus,Sargs)
  
  !<description>
    ! Command: READ2DPRM.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sfilename
    character(len=COLLCT_MLNAME) :: sname
    type(t_boundary), pointer :: p_rboundary
  
    if (size(Sargs) .lt. 3) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the identifier and file name 
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    call cmdprs_getparam (Sargs,3,sfilename,.false.,.true.)

    ! Open the file and read.
    inquire(file=trim(sfilename), exist=bexists)
    if (.not. bexists) then
      call output_line ("File not found!")
    else
      call output_line ("Reading file: "//trim(sfilename))
    
      ! Read
      allocate (p_rboundary)
      call boundary_read_prm(p_rboundary, sfilename)
      
      ! Remove old value from collection if present
      call do_destroy(rcmdStatus,sname,.false.)
      
      ! Add to the collection
      call collct_setvalue_bdry (rcmdStatus%rcollection, sname, p_rboundary, .true.)
    end if

  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_read2dtri (rcmdStatus,Sargs)
  
  !<description>
    ! Command: READ2DTRI.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sfilename, stoken
    character(len=COLLCT_MLNAME) :: sname
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
  
    if (size(Sargs) .lt. 3) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the identifier and file name 
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    call cmdprs_getparam (Sargs,3,sfilename,.false.,.true.)

    ! Open the file and read.
    inquire(file=trim(sfilename), exist=bexists)
    if (.not. bexists) then
      call output_line ("File not found!")
    else
      nullify(p_rboundary)
    
      ! Is there a boundary object attached?
      if (size(Sargs) .ge. 5) then
        call cmdprs_getparam (Sargs,4,stoken,.true.,.false.)
        if (stoken .eq. "--BOUNDARY") then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,5,stoken,.false.,.true.)
          
          p_rboundary => collct_getvalue_bdry(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Boundary object does not exist!")
          end if
        end if
      end if
    
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
      call do_destroy(rcmdStatus,sname,.false.)
      
      ! Add to the collection
      call collct_setvalue_tria (rcmdStatus%rcollection, sname, p_rtriangulation, .true.)
    end if

  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_meshrefine (rcmdStatus,Sargs)
  
  !<description>
    ! Command: MESHREFINE.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: stoken
    integer :: cmethod, icount
    integer :: i
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
  
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Get the identifier and file name 
    call cmdprs_getparam (Sargs,2,stoken,.true.,.false.)
    p_rtriangulation => collct_getvalue_tria(rcmdStatus%rcollection,stoken,bexists=bexists)
    if (.not. bexists) then
      call output_line ("Unknown variable!")
      return
    end if
    
    ! Read the other parameters
    cmethod = 0 ! 2-level ordering
    icount = 1
    nullify(p_rboundary)
  
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--BOUNDARY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rboundary => collct_getvalue_bdry(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Boundary object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      
      else if (stoken .eq. "--TIMES") then
        if (i .le. size(Sargs)) then
          ! Get the method type
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          ! Number of refinements
          read (stoken,*) icount
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--METHOD") then
        if (i .le. size(Sargs)) then
          ! Get the method type
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          if (stoken .eq. "--2LEVELORDERED") then
            ! 2-level ordering
            cmethod = 0
          else
            call output_line ("Warning: Unknown refinement method! Using default.")
            cmethod = 0
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      end if
    end do
  
    call output_line ("Refining mesh... [",bnolinebreak=.true.)
  
    ! Refine the mesh
    select case (cmethod)
    case (0)
      ! 2-level ordered
      
      do i = 1,icount
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
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_meshhierarchy (rcmdStatus,Sargs)
  
  !<description>
    ! Command: MESHHIERARCHY.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sname, stoken
    integer :: cmethod, ilevels
    integer :: i
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_meshHierarchy), pointer :: p_rmeshHierarchy
  
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Name of the hierarchy
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    
    ! Parse the options.
    cmethod = 0 ! 2-level ordering
    ilevels = 1
    nullify(p_rboundary)
    nullify(p_rtriangulation)
  
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--BOUNDARY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rboundary => collct_getvalue_bdry(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Boundary object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      
      else if (stoken .eq. "--MESH") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rtriangulation => collct_getvalue_tria(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Mesh object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      
      else if (stoken .eq. "--LEVELS") then
        if (i .le. size(Sargs)) then
          ! Get the method type
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          ! Number of refinements
          read (stoken,*) ilevels
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--METHOD") then
        if (i .le. size(Sargs)) then
          ! Get the method type
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          if (stoken .eq. "2LEVELORDERED") then
            ! 2-level ordering
            cmethod = 0
          else
            call output_line ("Warning: Unknown refinement method! Using default.")
            cmethod = 0
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      end if
    end do
    
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
    call do_destroy(rcmdStatus,sname,.false.)
    
    ! Add to the collection
    call collct_setvalue_mshh (rcmdStatus%rcollection, sname, p_rmeshHierarchy, .true.)
    
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
      integer, dimension(:), allocatable :: p_IelementIdsTri
      integer, dimension(:), allocatable :: p_IelementIdsQuad
      integer, dimension(:), allocatable :: p_IcubIdsTri
      integer, dimension(:), allocatable :: p_IcubIdsQuad

      nspaces = rcollection%IquickAccess(1)
      allocate(p_IelementIdsTri(nspaces))
      allocate(p_IelementIdsQuad(nspaces))
      allocate(p_IcubIdsTri(nspaces))
      allocate(p_IcubIdsQuad(nspaces))
      
      call collct_getvalue_intarr (rcollection, "IelementIdsTri", p_IelementIdsTri)
      call collct_getvalue_intarr (rcollection, "IelementIdsQuad", p_IelementIdsQuad)
      call collct_getvalue_intarr (rcollection, "IcubIdsTri", p_IcubIdsTri)
      call collct_getvalue_intarr (rcollection, "IcubIdsQuad", p_IcubIdsQuad)

      select case (rtriangulation%ndim)
      case (NDIM2D)
        ! Create a block discretisation of the specific size.
        call spdiscr_initBlockDiscr (rdiscr,nspaces,rtriangulation,rboundary)
        
        ! Create the sub-discretisation structures.
        do i=1,nspaces
          if (p_IelementIdsQuad(i) .eq. 0) then
            ! Pure tri space
            call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(i), &
                INT(p_IelementIdsTri(i),I32),INT(p_IcubIdsTri(i),I32),rtriangulation, &
                rboundary)
                
          else if (p_IelementIdsTri(i) .eq. 0) then
            ! Pure quad space
            call spdiscr_initDiscr_simple (rdiscr%RspatialDiscr(i), &
                INT(p_IelementIdsQuad(i),I32),INT(p_IcubIdsQuad(i),I32),rtriangulation, &
                rboundary)
          else
            ! Combined space
            call spdiscr_initDiscr_triquad (rdiscr%RspatialDiscr(i), &
                int(p_IelementIdsTri(i),I32), INT(p_IelementIdsQuad(i),I32),&
                int(p_IcubIdsTri(i),I32), INT(p_IcubIdsQuad(i),I32),rtriangulation, &
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
        IelementIdsTri,IelementIdsQuad,IcubIdsTri,IcubIdsQuad)
    
    !<description>
      ! Prepares setting up a discretisation.
      ! Must be called in advance to fgetDiscr.
    !</description>
    
    !<input>
      ! Number of components
      integer, intent(in) :: nspaces
    
      ! List of element id's for tri/tetra elements
      integer, dimension(:), intent(in), target :: IelementIdsTri

      ! List of element id's for quad/hexa elements
      integer, dimension(:), intent(in), target :: IelementIdsQuad

      ! List of cubature id's for tri/tetra elements
      integer, dimension(:), intent(in), target :: IcubIdsTri

      ! List of cubature id's for quad/hexa elements
      integer, dimension(:), intent(in), target :: IcubIdsQuad
    !</input>
    
    !<inputoutput>
      ! Current status block.
      type(t_commandstatus), intent(inout), target :: rcmdStatus

      ! Collection structure with information about the discretisation
      type(t_collection), intent(inout), optional :: rcollection
    !</inputoutput>

      rcollection%p_rnextCollection => rcmdStatus%rcollection
      
      rcollection%IquickAccess(1) = nspaces
      call collct_setvalue_intarr (rcollection, "IelementIdsTri", IelementIdsTri, .true.)
      call collct_setvalue_intarr (rcollection, "IelementIdsQuad", IelementIdsQuad, .true.)
      call collct_setvalue_intarr (rcollection, "IcubIdsTri", IcubIdsTri, .true.)
      call collct_setvalue_intarr (rcollection, "IcubIdsQuad", IcubIdsQuad, .true.)
      
    end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_fespace (rcmdStatus,Sargs)
  
  !<description>
    ! Command: FESPACE.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sname, stoken
    integer :: cmethod, ilevels, ncomponents
    integer :: i,j
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_feSpaceLevel), pointer :: p_rfeSpace
    type(t_feSpaceLevel), pointer :: p_rfeSpace1,p_rfeSpace2
    type (t_collection) :: rcollection
    integer, dimension(:), allocatable :: p_IelementIdsTri,p_IelementIdsQuad
    integer, dimension(:), allocatable :: p_IcubIdsTri,p_IcubIdsQuad
    
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Name of the hierarchy
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    
    ! Parse the options.
    cmethod = 0 ! 2-level ordering
    ncomponents = 0
    ilevels = 1
    nullify(p_rboundary)
    nullify(p_rtriangulation)
    
    ! For concatenation
    nullify(p_rfeSpace1)
    nullify(p_rfeSpace2)
  
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--BOUNDARY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rboundary => collct_getvalue_bdry(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Boundary object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      
      else if (stoken .eq. "--MESH") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rtriangulation => collct_getvalue_tria(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Mesh object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--COMPONENTS") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          read (stoken,*) ncomponents
          
          if (ncomponents .le. 0) then
            call output_line("Error. Number of components invalid!")
            return
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--TRIELEMENTS") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IelementIdsTri(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IelementIdsTri(j) = elem_igetID(stoken,.true.)
            if (p_IelementIdsTri(j) .eq. 0) then
              call output_line("Error. Invalid element ID: "//trim(stoken))
              return
            end if

            if (elem_igetDimension(p_IelementIdsTri(j)) .ne. NDIM2D) then
              call output_line("Error. Not a 2D element: "//trim(stoken))
              return
            end if

            if (elem_igetNVE(p_IelementIdsTri(j)) .ne. 3) then
              call output_line("Error. Not a tri element: "//trim(stoken))
              return
            end if
          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--QUADELEMENTS") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IelementIdsQuad(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IelementIdsQuad(j) = elem_igetID(stoken,.true.)
            if (p_IelementIdsQuad(j) .eq. 0) then
              call output_line("Error. Invalid element ID: "//trim(stoken))
              return
            end if

            if (elem_igetDimension(p_IelementIdsQuad(j)) .ne. NDIM2D) then
              call output_line("Error. Not a 2D element: "//trim(stoken))
              return
            end if

            if (elem_igetNVE(p_IelementIdsQuad(j)) .ne. 4) then
              call output_line("Error. Not a quad element: "//trim(stoken))
              return
            end if
          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--TRICUB") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IcubIdsTri(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IcubIdsTri(j) = cub_igetID(stoken,.true.)
            if (p_IelementIdsTri(j) .eq. 0) then
              call output_line("Error. Invalid cubature ID: "//trim(stoken))
              return
            end if

          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--QUADCUB") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IcubIdsQuad(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IcubIdsQuad(j) = cub_igetID(stoken,.true.)
            if (p_IelementIdsQuad(j) .eq. 0) then
              call output_line("Error. Invalid cubature ID: "//trim(stoken))
              return
            end if

          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--CONCAT") then
        if (i .lt. size(Sargs)-1) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          if (stoken .eq. sname) then
            call output_line("Error. Destination variable must be different from source!")
            return
          end if

          p_rfeSpace1 => collct_getvalue_fesp(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            call output_line("Error. 1st source FE space does not exist.")
            return
          end if

          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          if (stoken .eq. sname) then
            call output_line("Error. Destination variable must be different from source!")
            return
          end if

          p_rfeSpace2 => collct_getvalue_fesp(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            call output_line("Error. 2nd source FE space does not exist.")
            return
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      end if      

    end do
    
    if (.not. associated (p_rtriangulation)) then
      call output_line ("Invalid triangulation!")
      return
    end if
  
    ! Create the fe space.
    allocate(p_rfeSpace)
    
    if (associated (p_rfeSpace1)) then
      ! Concat two predefined FE spaces.
      call fesph_concatFeSpaces (p_rfeSpace1,p_rfeSpace2,p_rfeSpace,p_rtriangulation)
    else
      if ((.not. allocated(p_IelementIdsTri)) .and. &
          (.not. allocated(p_IelementIdsQuad))) then
        call output_line("Error. No element ID's specified!")
        return
      end if
    
      ! Create missing id arrays.
      if (.not. allocated(p_IelementIdsTri)) then
        allocate(p_IelementIdsTri(ncomponents))
        p_IelementIdsTri(:) = 0
      end if
      
      if (.not. allocated(p_IelementIdsQuad)) then
        allocate(p_IelementIdsQuad(ncomponents))
        p_IelementIdsQuad(:) = 0
      end if
    
      if (.not. allocated(p_IcubIdsTri)) then
        allocate(p_IcubIdsTri(ncomponents))
        p_IcubIdsTri(:) = 0
      end if

      ! Probably take the standard cubature formula.      
      where ((p_IcubIdsTri .eq. 0) .and. (p_IelementIdsTri .ne. 0)) 
        p_IcubIdsTri = spdiscr_getStdCubature(int(p_IelementIdsTri,I32))
      end where

      ! Probably take the standard cubature formula.      
      where ((p_IcubIdsQuad .eq. 0) .and. (p_IelementIdsQuad .ne. 0)) 
        p_IcubIdsQuad = spdiscr_getStdCubature(int(p_IelementIdsQuad,I32))
      end where
    
      if (.not. allocated(p_IcubIdsQuad)) then
        allocate(p_IcubIdsQuad(ncomponents))
        p_IcubIdsQuad(:) = 0
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
    
    ! Remove old value from collection if present
    call do_destroy(rcmdStatus,sname,.false.)
    
    ! Add to the collection
    call collct_setvalue_fesp (rcmdStatus%rcollection, sname, p_rfeSpace, .true.)
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_fehierarchy (rcmdStatus,Sargs)
  
  !<description>
    ! Command: FEHIERARCHY.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sname, stoken
    integer :: cmethod, ilevels, ncomponents
    integer :: i,j
    type(t_boundary), pointer :: p_rboundary
    type(t_meshhierarchy), pointer :: p_rmeshhierarchy
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type (t_collection) :: rcollection
    integer, dimension(:), allocatable :: p_IelementIdsTri,p_IelementIdsQuad
    integer, dimension(:), allocatable :: p_IcubIdsTri,p_IcubIdsQuad
    
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Name of the hierarchy
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    
    ! Parse the options.
    cmethod = 0 ! 2-level ordering
    ncomponents = 0
    ilevels = 1
    nullify(p_rboundary)
    nullify(p_rmeshhierarchy)
    
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--BOUNDARY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rboundary => collct_getvalue_bdry(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Boundary object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      
      else if (stoken .eq. "--MESHHIERARCHY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rmeshhierarchy => collct_getvalue_mshh(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rboundary)
            call output_line("Warning. Mesh hierarchy object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--COMPONENTS") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          read (stoken,*) ncomponents
          
          if (ncomponents .le. 0) then
            call output_line("Error. Number of components invalid!")
            return
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--TRIELEMENTS") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IelementIdsTri(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IelementIdsTri(j) = elem_igetID(stoken,.true.)
            if (p_IelementIdsTri(j) .eq. 0) then
              call output_line("Error. Invalid element ID: "//trim(stoken))
              return
            end if

            if (elem_igetDimension(p_IelementIdsTri(j)) .ne. NDIM2D) then
              call output_line("Error. Not a 2D element: "//trim(stoken))
              return
            end if

            if (elem_igetNVE(p_IelementIdsTri(j)) .ne. 3) then
              call output_line("Error. Not a tri element: "//trim(stoken))
              return
            end if
          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--QUADELEMENTS") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IelementIdsQuad(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IelementIdsQuad(j) = elem_igetID(stoken,.true.)
            if (p_IelementIdsQuad(j) .eq. 0) then
              call output_line("Error. Invalid element ID: "//trim(stoken))
              return
            end if

            if (elem_igetDimension(p_IelementIdsQuad(j)) .ne. NDIM2D) then
              call output_line("Error. Not a 2D element: "//trim(stoken))
              return
            end if

            if (elem_igetNVE(p_IelementIdsQuad(j)) .ne. 4) then
              call output_line("Error. Not a quad element: "//trim(stoken))
              return
            end if
          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--TRICUB") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IcubIdsTri(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IcubIdsTri(j) = cub_igetID(stoken,.true.)
            if (p_IelementIdsTri(j) .eq. 0) then
              call output_line("Error. Invalid cubature ID: "//trim(stoken))
              return
            end if

          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--QUADCUB") then
        if (ncomponents .le. 0) then
          call output_line("Error. Number of components undefined!")
          return
        end if
        
        if (i .le. size(Sargs)-ncomponents+1) then
          ! Allocate memory for element types and read.
          allocate (p_IcubIdsQuad(ncomponents))
          do j=1,ncomponents
            call cmdprs_getparam (Sargs,i,stoken,.false.,.false.)
            i = i+1
            
            ! Parse the Id.
            p_IcubIdsQuad(j) = cub_igetID(stoken,.true.)
            if (p_IelementIdsQuad(j) .eq. 0) then
              call output_line("Error. Invalid cubature ID: "//trim(stoken))
              return
            end if

          end do
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      end if      

    end do
    
    if (.not. associated (p_rmeshhierarchy)) then
      call output_line ("Invalid mesh hierarchy!")
      return
    end if
  
    ! Create the fe space.
    allocate(p_rfeHierarchy)
    
    if ((.not. allocated(p_IelementIdsTri)) .and. &
        (.not. allocated(p_IelementIdsQuad))) then
      call output_line("Error. No element ID's specified!")
      return
    end if
  
    ! Create missing id arrays.
    if (.not. allocated(p_IelementIdsTri)) then
      allocate(p_IelementIdsTri(ncomponents))
      p_IelementIdsTri(:) = 0
    end if
    
    if (.not. allocated(p_IelementIdsQuad)) then
      allocate(p_IelementIdsQuad(ncomponents))
      p_IelementIdsQuad(:) = 0
    end if
  
    if (.not. allocated(p_IcubIdsTri)) then
      allocate(p_IcubIdsTri(ncomponents))
      p_IcubIdsTri(:) = 0
    end if

    ! Probably take the standard cubature formula.      
    where ((p_IcubIdsTri .eq. 0) .and. (p_IelementIdsTri .ne. 0)) 
      p_IcubIdsTri = spdiscr_getStdCubature(int(p_IelementIdsTri,I32))
    end where

    ! Probably take the standard cubature formula.      
    where ((p_IcubIdsQuad .eq. 0) .and. (p_IelementIdsQuad .ne. 0)) 
      p_IcubIdsQuad = spdiscr_getStdCubature(int(p_IelementIdsQuad,I32))
    end where
  
    if (.not. allocated(p_IcubIdsQuad)) then
      allocate(p_IcubIdsQuad(ncomponents))
      p_IcubIdsQuad(:) = 0
    end if
  
    call output_line ("Creating FE hierarchy.")
  
    ! Create the FE space usinf fgetDiscr.
    call collct_init (rcollection)
    call prepare_fgetDiscr(rcollection,rcmdStatus,ncomponents,&
        p_IelementIdsTri,p_IelementIdsQuad,p_IcubIdsTri,p_IcubIdsQuad)
    if (associated(p_rboundary)) then
      call fesph_createHierarchy (p_rfeHierarchy,p_rmeshHierarchy%nlevels,&
          p_rmeshHierarchy,fgetDiscr,rcollection,rboundary=p_rboundary)
    else
      call output_line ("Warning: No boundary present!")
      
      call fesph_createHierarchy (p_rfeHierarchy,p_rmeshHierarchy%nlevels,&
          p_rmeshHierarchy,fgetDiscr,rcollection)
    end if
    call collct_done (rcollection)
    
    ! Remove old value from collection if present
    call do_destroy(rcmdStatus,sname,.false.)
    
    ! Add to the collection
    call collct_setvalue_feh (rcmdStatus%rcollection, sname, p_rfeHierarchy, .true.)
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_mlevelprjhierarchy (rcmdStatus,Sargs)
  
  !<description>
    ! Command: MLEVELPRJHIERARCHY.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sname, stoken
    integer :: i
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_interlevelProjectionHier), pointer :: p_rprjHier
    
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Name of the hierarchy
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    
    ! Parse the options.
    nullify(p_rfeHierarchy)
    
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--FEHIERARCHY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rfeHierarchy => collct_getvalue_feh(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rfeHierarchy)
            call output_line("Warning. FE hierarchy object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      end if      

    end do
    
    if (.not. associated (p_rfeHierarchy)) then
      call output_line ("Invalid mesh hierarchy!")
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
    call do_destroy(rcmdStatus,sname,.false.)
    
    ! Add to the collection
    call collct_setvalue_mlprjh (rcmdStatus%rcollection, sname, p_rprjHier, .true.)
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_createblockvector (rcmdStatus,Sargs)
  
  !<description>
    ! Command: READBLOCKVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sname, stoken, sfilename
    integer :: i,ilevel
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_vectorBlock), pointer :: p_rvectorBlock
    logical :: bformatted
    
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Name of the hierarchy
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    
    ! Parse the options.
    ilevel = 0
    nullify(p_rfeHierarchy)
    bformatted = .true.
    
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--FEHIERARCHY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          p_rfeHierarchy => collct_getvalue_feh(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rfeHierarchy)
            call output_line("Warning. FE hierarchy object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--LEVEL") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          read (stoken,*) ilevel
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      end if      

    end do
    
    if (.not. associated (p_rfeHierarchy)) then
      call output_line ("Invalid mesh hierarchy!")
      return
    end if
    
    if (ilevel .eq. 0) ilevel = p_rfeHierarchy%nlevels
  
    inquire(file=trim(sfilename), exist=bexists)
    
    call output_line ("Creating block vector.")

    ! Create the vector.
    allocate(p_rvectorBlock)
    call lsysbl_createVectorBlock (p_rfeHierarchy%p_rfeSpaces(ilevel)%p_rdiscretisation,&
        p_rvectorBlock,.true.)

    ! Remove old value from collection if present
    call do_destroy(rcmdStatus,sname,.false.)
    
    ! Add to the collection
    call collct_setvalue_vec (rcmdStatus%rcollection, sname, p_rvectorBlock, .true.)
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_readblockvector (rcmdStatus,Sargs)
  
  !<description>
    ! Command: READBLOCKVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    character(len=SYS_STRLEN) :: sname, stoken, sfilename
    integer :: i,ilevel
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_vectorBlock), pointer :: p_rvectorBlock
    logical :: bformatted
    
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Name of the hierarchy
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    
    ! Parse the options.
    ilevel = 1
    nullify(p_rfeHierarchy)
    bformatted = .true.
    
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--FILENAME") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          sfilename = trim(stoken)
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--UNFORMATTED") then

        bformatted = .false.

      end if      

    end do
    
    if (sfilename .eq. "") then
      call output_line ("Invalid filename!")
      return
    end if

    inquire(file=trim(sfilename), exist=bexists)
    
    if (.not. bexists) then
      call output_line ("File not found!")
    end if
    
    p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, sname)
    if (.not. associated(p_rvectorBlock)) then
      call output_line ("Invalid vector!")
      return
    end if
    
    call output_line ("Reading vector: "//trim(sfilename))

    ! Create the vector.
    call vecio_readBlockVectorHR (p_rvectorBlock, stoken, .true.,&
        0, sfilename, bformatted)

  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_writeblockvector (rcmdStatus,Sargs)
  
  !<description>
    ! Command: WRITEBLOCKVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    character(len=SYS_STRLEN) :: sname, stoken, sfilename, sformat
    integer :: i,ilevel
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    type(t_vectorBlock), pointer :: p_rvectorBlock
    logical :: bformatted
    
    if (size(Sargs) .lt. 2) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Name of the hierarchy
    call cmdprs_getparam (Sargs,2,sname,.true.,.false.)
    
    ! Parse the options.
    ilevel = 1
    nullify(p_rfeHierarchy)
    bformatted = .true.
    sformat = "(E20.10)"
    
    i = 3
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--FILENAME") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          sfilename = trim(stoken)
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--UNFORMATTED") then

        bformatted = .false.

      else if (stoken .eq. "--FORMAT") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          sformat = "("//trim(adjustl(stoken))//")"
          
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      end if      

    end do
    
    if (sfilename .eq. "") then
      call output_line ("Invalid filename!")
      return
    end if
    
    p_rvectorBlock => collct_getvalue_vec (rcmdStatus%rcollection, sname)
    if (.not. associated(p_rvectorBlock)) then
      call output_line ("Invalid vector!")
      return
    end if
    
    call output_line ("Writing vector: "//trim(sfilename))
    
    if (bformatted) then
      call vecio_writeBlockVectorHR (p_rvectorBlock, "vector", .true.,&
          0, sfilename, sformat)
    else
      call vecio_writeBlockVectorHR (p_rvectorBlock, "vector", .true.,&
          0, sfilename)
    end if

  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_copysubvector (rcmdStatus,Sargs)
  
  !<description>
    ! Command: WRITEBLOCKVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    type(t_vectorBlock), pointer :: p_rvectorBlock1,p_rvectorBlock2
    character(len=COLLCT_MLNAME) :: ssource, sdest, stoken
    integer :: icomp1, icomp2
    
    if (size(Sargs) .lt. 5) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Source vector and component
    call cmdprs_getparam (Sargs,2,ssource,.true.,.false.)
    call cmdprs_getparam (Sargs,3,stoken,.false.,.false.)
    read (stoken,*) icomp1
    
    ! Destination
    call cmdprs_getparam (Sargs,4,sdest,.true.,.false.)
    call cmdprs_getparam (Sargs,5,stoken,.false.,.false.)
    read (stoken,*) icomp2

    p_rvectorBlock1 => collct_getvalue_vec (rcmdStatus%rcollection, ssource)
    if (.not. associated(p_rvectorBlock1)) then
      call output_line ("Invalid source vector!")
      return
    end if
    
    if ((icomp1 .lt. 1) .or. (icomp1 .gt. p_rvectorBlock1%nblocks)) then
      call output_line ("Invalid source component!")
      return
    end if

    p_rvectorBlock2 => collct_getvalue_vec (rcmdStatus%rcollection, sdest)
    if (.not. associated(p_rvectorBlock2)) then
      call output_line ("Invalid destination vector!")
      return
    end if

    if ((icomp2 .lt. 1) .or. (icomp2 .gt. p_rvectorBlock1%nblocks)) then
      call output_line ("Invalid destination component!")
      return
    end if
    
    call lsyssc_copyVector (p_rvectorBlock1%RvectorBlock(icomp1),&
        p_rvectorBlock2%RvectorBlock(icomp2))
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_copyvector (rcmdStatus,Sargs)
  
  !<description>
    ! Command: WRITEBLOCKVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    type(t_vectorBlock), pointer :: p_rvectorBlock1,p_rvectorBlock2
    character(len=COLLCT_MLNAME) :: ssource, sdest
    
    if (size(Sargs) .lt. 3) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Source vector and component
    call cmdprs_getparam (Sargs,2,ssource,.true.,.false.)
    call cmdprs_getparam (Sargs,3,sdest,.true.,.false.)

    p_rvectorBlock1 => collct_getvalue_vec (rcmdStatus%rcollection, ssource)
    if (.not. associated(p_rvectorBlock1)) then
      call output_line ("Invalid source vector!")
      return
    end if
    
    p_rvectorBlock2 => collct_getvalue_vec (rcmdStatus%rcollection, sdest)
    if (.not. associated(p_rvectorBlock2)) then
      call output_line ("Invalid destination vector!")
      return
    end if

    call lsysbl_copyVector (p_rvectorBlock1,p_rvectorBlock2)
    
  end subroutine

  ! ***************************************************************************

  subroutine cmdprs_do_interpolatevector (rcmdStatus,Sargs)
  
  !<description>
    ! Command: INTERPOLATEVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    logical :: bexists
    type(t_vectorBlock), pointer :: p_rvectorBlock1,p_rvectorBlock2
    character(len=COLLCT_MLNAME) :: ssource, sdest, stoken
    type(t_interlevelProjectionHier), pointer :: p_rprjHier
    type(t_feHierarchy), pointer :: p_rfeHierarchy
    integer :: i,isource,idest
    type(t_vectorBlock), pointer :: p_rtemp1,p_rtemp2
    
    if (size(Sargs) .lt. 6) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Source vector and component
    call cmdprs_getparam (Sargs,2,ssource,.true.,.false.)
    call cmdprs_getparam (Sargs,3,stoken,.true.,.false.)
    read (stoken,*) isource
    
    ! Destination
    call cmdprs_getparam (Sargs,4,sdest,.true.,.false.)
    call cmdprs_getparam (Sargs,5,stoken,.true.,.false.)
    read (stoken,*) idest

    p_rvectorBlock1 => collct_getvalue_vec (rcmdStatus%rcollection, ssource)
    if (.not. associated(p_rvectorBlock1)) then
      call output_line ("Invalid source vector!")
      return
    end if

    p_rvectorBlock2 => collct_getvalue_vec (rcmdStatus%rcollection, sdest)
    if (.not. associated(p_rvectorBlock2)) then
      call output_line ("Invalid destination vector!")
      return
    end if

    nullify(p_rfeHierarchy)
    nullify(p_rprjHier)

    i = 5
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--MLPRJHIERARCHY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1

          p_rprjHier => collct_getvalue_mlprjh(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rprjHier)
            call output_line("Warning. Multilevel projection hierarchy does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      else if (stoken .eq. "--FEHIERARCHY") then
        if (i .le. size(Sargs)) then
          ! Get the name of the attached boundary object -- and the object
          call cmdprs_getparam (Sargs,i,stoken,.false.,.true.)
          i = i+1
          
          p_rfeHierarchy => collct_getvalue_feh(rcmdStatus%rcollection,stoken,bexists=bexists)
          if (.not. bexists) then
            nullify(p_rfeHierarchy)
            call output_line("Warning. FE hierarchy object does not exist!")
          end if
        else
          call output_line("Error. Invalid parameters!")
          return
        end if

      end if      

    end do
    
    if (.not. associated (p_rprjHier)) then
      call output_line ("Invalid multilevel projection hierarchy!")
      return
    end if
    
    if ((isource .lt. p_rprjHier%nlmin) .or. (isource .gt. p_rprjHier%nlmax)) then
      call output_line ("Invalid source level!")
      return
    end if
    
    if ((idest .lt. p_rprjHier%nlmin) .or. (idest .gt. p_rprjHier%nlmax)) then
      call output_line ("Invalid destination level!")
      return
    end if
    
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

  subroutine cmdprs_do_l2projection (rcmdStatus,Sargs)
  
  !<description>
    ! Command: INTERPOLATEVECTOR.
  !</description>
  
  !<inputoutput>
    ! Current status block.
    type(t_commandstatus), intent(inout) :: rcmdStatus
  !</inputoutput>

  !<input>
    ! Command line arguments.
    type(VARYING_STRING), dimension(:), intent(in) :: Sargs
  !</input>
  
    ! local variables
    type(t_vectorBlock), pointer :: p_rvectorBlock1,p_rvectorBlock2
    character(len=COLLCT_MLNAME) :: ssource, sdest, stoken
    type(t_matrixScalar) :: rmatrixMass
    integer :: i
    logical :: bverbose
    type(t_collection) :: rcollection
    type(t_configL2ProjectionByMass) :: rL2ProjectionConfig
    
    if (size(Sargs) .lt. 3) then
      call output_line ("Not enough arguments.")
      return
    end if
    
    ! Source vector and destination
    call cmdprs_getparam (Sargs,2,ssource,.true.,.false.)
    call cmdprs_getparam (Sargs,3,sdest,.true.,.false.)

    p_rvectorBlock1 => collct_getvalue_vec (rcmdStatus%rcollection, ssource)
    if (.not. associated(p_rvectorBlock1)) then
      call output_line ("Invalid source vector!")
      return
    end if

    p_rvectorBlock2 => collct_getvalue_vec (rcmdStatus%rcollection, sdest)
    if (.not. associated(p_rvectorBlock2)) then
      call output_line ("Invalid destination vector!")
      return
    end if

    i = 4
    do 
      if (i .gt. size(Sargs)) exit
      
      ! Get the next token.
      call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
      i = i+1
      
      if (stoken .eq. "--VERBOSE") then
        bverbose = .true.
      end if      

      if (stoken .eq. "--RELERROR") then
        if (i .le. size(Sargs)) then
          call cmdprs_getparam (Sargs,i,stoken,.true.,.false.)
          i = i+1
          read (stoken,*) rL2ProjectionConfig%depsrel
        else
          call output_line("Error. Invalid parameters!")
          return
        end if
      end if      

    end do

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

  end subroutine

end module
