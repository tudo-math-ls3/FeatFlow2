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
!# 7.) io_isDirectory
!#     -> Checks whether a given string is a directory
!#
!# 8.) io_makeDirectory
!#     -> Creates a directory (including parent directories if needed)
!#
!# </purpose>
!##############################################################################

module io

!$ use omp_lib
  use fsystem
  use genoutput

  implicit none

  private

!<constants>

  !<constantblock description="Input/output block type identifiers">

  ! defines the default value for files
  integer, parameter, public :: IO_UNKNOWN = 0

  ! defines that a file must already exist
  integer, parameter, public :: IO_OLD = 1

  ! defines that a file must not exist
  integer, parameter, public :: IO_NEW = 2

  ! defines that an existing file should be replaced
  integer, parameter, public :: IO_REPLACE = 3

  ! defines that a temporary file should be deleted when closed
  integer, parameter, public :: IO_SCRATCH = 4

  !</constantblock>

!</constants>

  public :: io_openFileForReading
  public :: io_openFileForWriting
  public :: io_readlinefromfile
  public :: io_deleteFile
  public :: io_pathExtract
  public :: io_pathConcat
  public :: io_isDirectory

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
    logical, intent(in), optional :: bformatted

!</input>

!<output>

    !number of unit
    integer, intent(out) :: iunit
!</output>
!</subroutine>

    logical :: bexists !true, if a file with name sfilename exists
    integer :: istatus !status variable for opening. 0, if opening succesful


    if (trim(sfilename) .eq. "") then
      call output_line('File name <'//trim(adjustl(sfilename))//'> empty.', &
                       OU_CLASS_ERROR,OU_MODE_STD,'io_openFileForReading')
      return
    endif

    iunit = sys_getFreeUnit()
    if (iunit .eq. -1) then
      call output_line('No free unit found. Not able to open the file.', &
                       OU_CLASS_ERROR,OU_MODE_STD,'io_openFileForReading')
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
        call output_line('Error while opening file <'//trim(adjustl(sfilename))//'> for reading.', &
                         OU_CLASS_ERROR,OU_MODE_STD,'io_openFileForReading')
        iunit = -1
      end if

    else

      call output_line('File <'//trim(adjustl(sfilename))//'> does not exist.', &
                       OU_CLASS_ERROR,OU_MODE_STD,'io_openFileForReading')
      call sys_halt()
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
    logical, intent(in), optional :: bformatted
!</input>

!<output>

    ! unit of the opened file
    integer, intent(out) :: iunit

    ! optional parameter (see description)
    logical, intent(out), optional :: bfileExists

!</output>
!</subroutine>

    logical :: bexists !true, if the file to be written in exists
    integer :: istatus !status variable for opening procedure

    ! the result "dirname(sfilename)" would yield
    character(len=len(sfilename)) :: sfilepath

    if (len_trim(sfilename) .eq. 0) then
      call output_line('File name <'//trim(adjustl(sfilename))//'> empty.', &
                       OU_CLASS_ERROR,OU_MODE_STD,'io_openFileForWriting')
      return
    end if

    iunit = sys_getFreeUnit()
    if (iunit .eq. -1) then
      call output_line('No free unit found. Not able to open the file.', &
                       OU_CLASS_ERROR,OU_MODE_STD,'io_openFileForWriting')
      return
    end if

    inquire(file=trim(sfilename), exist=bexists)
    if (.not. present(bformatted)) then
      if (bexists .and. cflag .eq. SYS_REPLACE) then
        open(unit=iunit, file=trim(sfilename), iostat=istatus, status="replace", &
            action="write")
      else
        ! Ensure that the path up to the given file name does exist
        call io_pathExtract(sfilename, sfilepath)
        if (.not. io_isDirectory(sfilepath)) then
          call io_makeDirectory(sfilepath)
        end if
        open(unit=iunit, file=trim(sfilename), iostat=istatus, action="write", &
            position="append")
      end if
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
        ! Ensure that the path up to the given file name does exist. If it does not,
        ! create it
        call io_pathExtract(sfilename, sfilepath)
        if (.not. io_isDirectory(sfilepath)) then
          call io_makeDirectory(sfilepath)
        end if
        if (bformatted) then
          open(unit=iunit, file=trim(sfilename), iostat=istatus, action="write", &
              position="append", form="formatted")
        else
          open(unit=iunit, file=trim(sfilename), iostat=istatus, action="write", &
              position="append", form="unformatted")
        end if
      end if
    end if
    if (present(bfileExists)) then
      bfileExists = bexists
    end if
    if (istatus .ne. 0) then
      call output_line('Error while opening file <'//trim(adjustl(sfilename))//'> for writing.', &
                         OU_CLASS_ERROR,OU_MODE_STD,'io_openFileForWriting')
      iunit = -1
    end if

  end subroutine io_openFileForWriting

! ***************************************************************************

!<subroutine>
  subroutine io_deleteFile(sfilename)

!<description>
  ! this routione deletes a file sfilename.
!</description>

!<input>
  ! filename
  character(*), intent(in) :: sfilename
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
    integer, intent(in) :: iunit
!</input>

!<output>
    ! The string where to write data to
    character(LEN=*), intent(out) :: sdata

    ! Length of the output
    integer, intent(out) :: ilinelen

    ! Status of the reading process. Set to a value <> 0 if the end
    ! of the file is reached.
    integer, intent(out) :: ios
!</output>
!</subroutine>

    ! local variables
    character :: c

    sdata = ''
    ilinelen = 0

    ! Read the data - as long as the line/file does not end.
    do

      ! Read a character.
      ! Unfortunately, Fortran forces me to use this dirty GOTO
      ! to decide processor-independently whether the line or
      ! the record ends.
      read (unit=iunit, fmt='(A1)', iostat=ios, advance='NO',&
            end=10, eor=20) c

      ! Do not do anything in case of an error
      if (ios .eq. 0) then

        ilinelen = ilinelen + 1
        sdata (ilinelen:ilinelen) = c

      end if

      ! Proceed to next character
      cycle

      ! End of file.
10    ios = -1
      exit

      ! End of record = END OF LINE.
20    ios = 0
      exit

    end do

  end subroutine io_readlinefromfile

  ! ***************************************************************************

  !<subroutine>

  subroutine io_pathExtract (sfile, sfilepath, sfilename, babsolute)

  !<description>
    ! Extracts the path of a file from a path+filename string.
  !</description>

  !<input>
    ! Filename + path of a specific file (or directory).
    character(len=*), intent(in) :: sfile
  !</input>

  !<output>
    ! OPTIONAL: Receives the directory that contains the specific file,
    ! or "" if no directory was specified in sfile.
    character(len=*), intent(out), optional :: sfilepath

    ! OPTIONAL: Receives the name of the file without a probably preceding
    ! directory string.
    character(len=*), intent(out), optional :: sfilename

    ! OPTIONAL: Returns TRUE if the path specification in sfile points to an
    ! absolute path. Returns FALSE if the path in sfile is relative.
    logical, intent(out), optional :: babsolute
  !</output>

  !</subroutine>

    integer :: i
    character(len=10) :: ssubpath

    ! Find the last "/" or "\" in sfile.                                (!" cpp fix)
    ! Note that we specified "\\" and not "\" because the PGI compiler  (!" cpp fix)
    ! (stupid thing) would otherwise use the backslash to escape the quote
    ! character. So PGI sees "/\" and other compiler see "/\\", but this (!" cpp fix)
    ! does not matter since the string must only contain a couple of
    ! delimiters which may occur more than once in the string.
    i = scan(sfile,"/\\",.true.)
    if (i .ne. 0) then
      ! Directory ends at position i.
      if (present(sfilepath)) sfilepath = sfile(1:i-1)
      if (present(sfilename)) sfilename = sfile(i+1:)
    else
      ! No directory specified.
      if (present(sfilepath)) sfilepath = ""
      if (present(sfilename)) sfilename = sfile
    end if

    if (present(babsolute)) then
      ! Take a look if this is an absolute or relative path.
      i = scan(trim(adjustl(sfile)),"/\\",.false.)
      babsolute = i .eq. 1

      ! In Windows environments, the path is also absolute if
      ! a volume descriptor like "C:" precedes the (back-)slash.
      if (.not. babsolute) then
        if (i .eq. 3) then
          ! Extract the first 10 characters and check
          ssubpath = trim(adjustl(sfile))
          if (ssubpath(2:2) .eq. ":") then
            babsolute = .true.
          end if
        end if
      end if
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
    character(len=len_trim(spath)+len_trim(sfilename)+1) :: sfile
  !</result>

  !</function>

    sfile = trim(spath)//"/"//trim(sfilename)

  end function


  ! ***************************************************************************


!<function>
  function io_isDirectory (spath) result (bisDir)

  !<description>
    ! Checks whether a given string is a directory
  !</description>

  !<input>
    ! Path to the file.
    character(len=*), intent(in) :: spath
  !</input>

  !<result>
    ! Is .TRUE. if the given string is an existing directory
    logical :: bisDir
  !</result>

!</function>

    integer(I32) :: ierr
    integer(I32) :: iisdir

    ! use wrapper for C system call stat() in kernel/System/isdirectory.c
    external isdirectory
    ierr = 0
    bisDir = .false.
    call isdirectory(trim(spath) // achar(0), iisdir, ierr)
    if (ierr .ge. 0) then
      bisDir = (iisdir .gt. 0)
    endif

  end function io_isDirectory


  ! ***************************************************************************


!<subroutine>
  subroutine io_makeDirectory (spath)

  !<description>
    ! Portably creates directories (even recursively)
  !</description>

  !<input>
    ! Path to the file.
    character(len=*), intent(in) :: spath
  !</input>

!</subroutine>

    external mkdir_recursive

    ! status variable for (recursive) directory creation
    ! > 0: number of subdirectories created
    ! < 0: -errno
    integer :: istatus


    if (len_trim(spath) .eq. 0) then
      call output_line('Path name <'//trim(adjustl(spath))//'> empty.', &
                       OU_CLASS_ERROR,OU_MODE_STD,'io_makeDirectory')
      return
    end if

    ! Ensure to explicitly NULL-terminate the string when calling C function!
    call mkdir_recursive(trim(spath) // achar(0), istatus)
    if (istatus .lt. 0) then
      call output_line('Could not (auto-)create path <'//trim(adjustl(spath))//'>.', &
                       OU_CLASS_ERROR,OU_MODE_STD,'io_makeDirectory')
      call sys_halt()
    end if

  end subroutine io_makeDirectory

end module io
