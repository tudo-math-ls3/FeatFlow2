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
!# </purpose>
!##############################################################################

MODULE io

  USE fsystem
  USE error

  IMPLICIT NONE

!<constants>

  !<constantblock description="Input/output block type identifiers">
  
  ! defines the default value for files
  INTEGER, PARAMETER :: IO_UNKNOWN = 0

  ! defines that a file must already exist
  INTEGER, PARAMETER :: IO_OLD = 1
  
  ! defines that a file must not exist
  INTEGER, PARAMETER :: IO_NEW = 2

  ! defines that an existing file should be replaced
  INTEGER, PARAMETER :: IO_REPLACE = 3

  ! defines that a temporary file should be deleted when closed
  INTEGER, PARAMETER :: IO_SCRATCH = 4
    
  !</constantblock>

!</constants>


CONTAINS

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
    LOGICAL, INTENT(IN), OPTIONAL :: bformatted

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
    
      IF (.NOT. PRESENT(bformatted)) THEN
        open(unit=iunit, file=trim(sfilename), iostat=istatus, action="read")
      ELSE IF (bformatted) THEN
        open(unit=iunit, file=trim(sfilename), iostat=istatus, action="read",&
             form="formatted")
      ELSE
        open(unit=iunit, file=trim(sfilename), iostat=istatus, action="read",&
             form="unformatted")
      END IF
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
    LOGICAL, INTENT(IN), OPTIONAL :: bformatted
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
    IF (.NOT. PRESENT(bformatted)) THEN
      if (bexists .and. cflag .eq. SYS_REPLACE) then
        open(unit=iunit, file=trim(sfilename), iostat=istatus, status="replace", &
            action="write")
      else
        open(unit=iunit, file=trim(sfilename), iostat=istatus, action="write", &
            position="append")
      endif
    ELSE
      if (bexists .and. cflag .eq. SYS_REPLACE) then
        IF (bformatted) THEN
          open(unit=iunit, file=trim(sfilename), iostat=istatus, status="replace", &
              action="write", form="formatted")
        ELSE
          open(unit=iunit, file=trim(sfilename), iostat=istatus, status="replace", &
              action="write", form="unformatted")
        END IF
      else
        IF (bformatted) THEN
          open(unit=iunit, file=trim(sfilename), iostat=istatus, action="write", &
              position="append", form="formatted")
        ELSE
          open(unit=iunit, file=trim(sfilename), iostat=istatus, action="write", &
              position="append", form="unformatted")
        END IF
      endif    
    END IF
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
  SUBROUTINE io_deleteFile(sfilename)

!<description>
  ! this routione deletes a file sfilename.
!</description>

!<input>
  ! filename
  CHARACTER(*), INTENT(IN) :: sfilename
!</input>
!</subroutine>

    INTEGER :: iunit
    
    ! Open the file for writing, overwrite the old one.
    CALL io_openFileForWriting(sfilename, iunit, SYS_REPLACE)
    
    ! Close the file and delete it.
    CLOSE (iunit, STATUS='DELETE')

  END SUBROUTINE 

! ***************************************************************************

!<subroutine>
  
  SUBROUTINE io_readlinefromfile (iunit, sdata, ilinelen, ios)

!<description>
  !This routine reads a line from a text file
!</description>

!<input>  
    ! The unit where to read from; must be connected to a file.
    INTEGER, INTENT(IN) :: iunit
!</input>  

!<output>
    ! The string where to write data to
    CHARACTER(LEN=*), INTENT(OUT) :: sdata
    
    ! Length of the output
    INTEGER, INTENT(OUT) :: ilinelen
    
    ! Status of the reading process. Set to a value <> 0 if the end
    ! of the file is reached.
    INTEGER, INTENT(OUT) :: ios
!</output>
!</subroutine>
    
    ! local variables
    INTEGER :: eol
    CHARACTER :: c
    
    sdata = ''
    ilinelen = 0
    
    ! Read the data - as long as the line/file does not end.
    eol = NO
    ios = 0
    DO WHILE ((ios .EQ. 0) .AND. (eol .EQ. NO))
      
      ! Read a character.
      ! Unfortunately, Fortran forces me to use this dirty GOTO
      ! to decide processor-independently whether the line or
      ! the record ends.
      READ (unit=iunit,fmt='(A1)',iostat=ios,advance='NO', end=10, eor=20) c
      GOTO 30
      
10    CONTINUE
      ! End of file. 
      ios = -1
      GOTO 30
      
20    CONTINUE
      ! End of record = END OF LINE.
      eol = YES
      
      ! Set error flag back to 0.
      ios = 0
      
30    CONTINUE    
      ! Don't do anything in case of an error
      IF (ios .EQ. 0) THEN
        
        ilinelen = ilinelen + 1
        sdata (ilinelen:ilinelen) = c
        
      END IF
      
    END DO
    
  END SUBROUTINE

END MODULE
