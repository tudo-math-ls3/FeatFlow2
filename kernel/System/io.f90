!##############################################################################
!# ****************************************************************************
!# <name> io </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains several routines for input/output purposes.
!#
!# The following routines can be found here:
!#
!# 1.) io_openFileForReading
!#     -> Opens a file for reading; the file handle is automatically determined
!#
!# 2.) io_openFileForWriting
!#     -> Opens a file for writing; the file handle is automatically determined
!#
!# 3.) io_readlinefromfile
!#     -> Reads one line from an opend file 
!#
!# 4.) io_deleteFile
!#     -> Deletes a file.
!#
!# 5.) io_pathExtract
!#     -> Extracts path information from a filename
!#
!# 6.) io_pathConcat
!#     -> Concatenates a filename to a path.
!#
!# </purpose>
!##############################################################################

module io

  use fsystem
  use error

  implicit none

!<constants>

  !<constantblock description="Input/output block type identifiers">
  
  ! defines the default value for files
  integer, parameter :: IO_UNKNOWN = 0

  ! defines that a file must already exist
  integer, parameter :: IO_OLD = 1
  
  ! defines that a file must not exist
  integer, parameter :: IO_NEW = 2

  ! defines that an existing file should be replaced
  integer, parameter :: IO_REPLACE = 3

  ! defines that a temporary file should be deleted when closed
  integer, parameter :: IO_SCRATCH = 4
    
  !</constantblock>

!</constants>


contains

!************************************************************************************

!<subroutine>
  subroutine io_openFileForReading(sfilename, iunit, bformatted)

!<description>
    !This routine tries to open a file for reading. If succesful, on can read from it
    !via unit "iunit". Otherwise, iunit is -1.
!</description>

!<input>

    !filename
    character(*), intent(in) :: sfilename

    ! OPTIONAL: 
    ! TRUE : Open the file formatted, i.e. in human readable form
    ! FALSE: Open the file in unformatted, machine dependent form
    ! If not specified, the default system dependent setting is used.
    logical, intent(IN), optional :: bformatted

!</input>

!<output>

    !number of unit
    integer, intent(out) :: iunit
!</output>
!</subroutine>

    logical :: bexists !true, if a file with name sfilename exists
    integer :: istatus !status variable for opening. 0, if opening succesful


    if (trim(sfilename) .eq. "") then
      call error_print(ERR_IO_EMPTYFILENAME, "io_openFileForReading", ERR_NOT_CRITICAL, &
                       sarg1 = "sfilename")
      return
    endif

    iunit = sys_getFreeUnit()
    if (iunit .eq. -1) then
      call error_print(ERR_IO_NOFREEUNIT, "io_openFileForReading", ERR_NOT_CRITICAL)

      !give it up
      return
    endif

    inquire(file=trim(sfilename), exist=bexists)

    if (bexists) then
    
      if (.not. present(bformatted)) then
        open(unit=iunit, file=trim(sfilename), iostat=istatus, action="read")
      else if (bformatted) then
        open(unit=iunit, file=trim(sfilename), iostat=istatus, action="read",&
             form="formatted")
      else
        open(unit=iunit, file=trim(sfilename), iostat=istatus, action="read",&
             form="unformatted")
      end if
      if (istatus .ne. 0) then
        write(unit=*,fmt=*) "*** Error while opening file '",trim(sfilename),"'. ***"
        iunit = -1
      end if
      
    else
    
      call error_print(ERR_IO_NOSUCHFILE, "io_openFileForReading", ERR_CRITICAL, &
                      sarg1 = sfilename)
                      
    endif

  end subroutine io_openFileForReading

!************************************************************************


!<subroutine>
  subroutine io_openFileForWriting(sfilename, iunit, cflag, bfileExists, bformatted)

!<description>
    !This routine tries to open a file for writing. If succesful, one can write to it
    !via unit "iunit". Otherwise, iunit is -1. cflag specifies if an already existing file
    !should be replaced or if the output should be appended. bfileExists is an optional
    !parameter, which will be set to true, if the file already existed, otherwise false.
!</description>

!<input>

    !filename
    character(*), intent(in) :: sfilename

    !mode: SYS_APPEND or SYS_REPLACE
    integer, intent(in) :: cflag

    ! OPTIONAL: 
    ! TRUE : Open the file formatted, i.e. in human readable form
    ! FALSE: Open the file in unformatted, machine dependent form
    ! If not specified, the default system dependent setting is used.
    logical, intent(IN), optional :: bformatted
!</input>

!<output>

    !unit of the opened file
    integer, intent(out) :: iunit

    !optional parameter (see description)
    logical, intent(out),optional :: bfileExists
    
!</output>
!</subroutine>

    logical :: bexists !true, if the file to be written in exists
    integer :: istatus !status variable for opening procedure

    if (trim(sfilename) .eq. "") then 
      call error_print(ERR_IO_EMPTYFILENAME, "io_openFileForWriting", ERR_NOT_CRITICAL, &
                       sarg1 = "sfilename")
      return
    endif

    iunit = sys_getFreeUnit()
    if (iunit .eq. -1) then
      call error_print(ERR_IO_NOFREEUNIT, "io_openFileForWriting", ERR_NOT_CRITICAL)
      return
    endif

    inquire(file=trim(sfilename), exist=bexists)
    if (.not. present(bformatted)) then
      if (bexists .and. cflag .eq. SYS_REPLACE) then
        open(unit=iunit, file=trim(sfilename), iostat=istatus, status="replace", &
            action="write")
      else
        open(unit=iunit, file=trim(sfilename), iostat=istatus, action="write", &
            position="append")
      endif
    else
      if (bexists .and. cflag .eq. SYS_REPLACE) then
        if (bformatted) then
          open(unit=iunit, file=trim(sfilename), iostat=istatus, status="replace", &
              action="write", form="formatted")
        else
          open(unit=iunit, file=trim(sfilename), iostat=istatus, status="replace", &
              action="write", form="unformatted")
        end if
      else
        if (bformatted) then
          open(unit=iunit, file=trim(sfilename), iostat=istatus, action="write", &
              position="append", form="formatted")
        else
          open(unit=iunit, file=trim(sfilename), iostat=istatus, action="write", &
              position="append", form="unformatted")
        end if
      endif    
    end if
    if (present(bfileExists)) then
      bfileExists = bexists
    endif
    if (istatus .ne. 0) then
      write(unit=*,fmt=*) "*** Error while opening file '", trim(sfilename), "'. ***"
      iunit = -1
    endif

  end subroutine io_openFileForWriting

! ***************************************************************************

!<subroutine>
  subroutine io_deleteFile(sfilename)

!<description>
  ! this routione deletes a file sfilename.
!</description>

!<input>
  ! filename
  character(*), intent(IN) :: sfilename
!</input>
!</subroutine>

    integer :: iunit
    
    ! Open the file for writing, overwrite the old one.
    call io_openFileForWriting(sfilename, iunit, SYS_REPLACE)
    
    ! Close the file and delete it.
    close (iunit, STATUS='DELETE')

  end subroutine io_deleteFile

! ***************************************************************************

!<subroutine>
  
  subroutine io_readlinefromfile (iunit, sdata, ilinelen, ios)

!<description>
  !This routine reads a line from a text file
!</description>

!<input>  
    ! The unit where to read from; must be connected to a file.
    integer, intent(IN) :: iunit
!</input>  

!<output>
    ! The string where to write data to
    character(LEN=*), intent(OUT) :: sdata
    
    ! Length of the output
    integer, intent(OUT) :: ilinelen
    
    ! Status of the reading process. Set to a value <> 0 if the end
    ! of the file is reached.
    integer, intent(OUT) :: ios
!</output>
!</subroutine>
    
    ! local variables
    integer :: eol
    character :: c
    
    sdata = ''
    ilinelen = 0
    
    ! Read the data - as long as the line/file does not end.
    eol = NO
    ios = 0
    do while ((ios .eq. 0) .and. (eol .eq. NO))
      
      ! Read a character.
      ! Unfortunately, Fortran forces me to use this dirty GOTO
      ! to decide processor-independently whether the line or
      ! the record ends.
      read (unit=iunit,fmt='(A1)',iostat=ios,advance='NO', end=10, eor=20) c
      goto 30
      
10    continue
      ! End of file. 
      ios = -1
      goto 30
      
20    continue
      ! End of record = END OF LINE.
      eol = YES
      
      ! Set error flag back to 0.
      ios = 0
      
30    continue    
      ! Don't do anything in case of an error
      if (ios .eq. 0) then
        
        ilinelen = ilinelen + 1
        sdata (ilinelen:ilinelen) = c
        
      end if
      
    end do
    
  end subroutine io_readlinefromfile

  ! ***************************************************************************

  !<subroutine>

  subroutine io_pathExtract (sfile, sfilepath, sfilename)
  
  !<description>
    ! Extracts the path of a file from a path+filename string.
  !</description>
  
  !<input>
    ! Filename + path of a specific file (or directory).
    character(len=*), intent(in) :: sfile
  !</input>
  
  !<output>
    ! Receives the directory that contains the specific file, 
    ! or "" if no directory was specified in sfile.
    character(len=*), intent(out) :: sfilepath

    ! OPTIONAL: Receives the name of the file without a probably preceding
    ! directory string.
    character(len=*), intent(out), optional :: sfilename
  !</output>
  
  !</subroutine>
  
    integer :: i
  
    ! Find the last "/" or "\" in sfile.
    ! Note that we specified "\\" and not "\" because the PGI compiler
    ! (stupid thing) would otherwise use the backslash to escape the quote
    ! character. So PGI sees "/\" and other compiler see "/\\", but this
    ! doesn't matter since the string must only contain a couple of
    ! delimiters which may occur more than once in the string.
    i = scan(sfile,"/\\",.true.)
    if (i .ne. 0) then
      ! Directory ends at position i.
      sfilepath = sfile(1:i-1)
      if (present(sfilename)) sfilename = sfile(i+1:)
    else
      ! No directory specified.
      sfilepath = ""
      if (present(sfilename)) sfilename = sfile
    end if
  
  end subroutine

  ! ***************************************************************************

  !<function>

  function io_pathConcat (spath,sfilename) result (sfile)
  
  !<description>
    ! Concatenates a filename to a path specifier.
  !</description>
  
  !<input>
    ! Path to the file.
    character(len=*), intent(in) :: spath

    ! Name of the file (or directory)
    character(len=*), intent(in) :: sfilename
  !</input>
  
  !<result>
    ! Path + filename to a specific file (or directory).
    character(len=LEN_TRIM(spath)+LEN_TRIM(sfilename)+1) :: sfile
  !</result>
  
  !</function>
  
    sfile = trim(spath)//"/"//trim(sfilename)
  
  end function

end module io
