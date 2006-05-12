!##############################################################################
!# ****************************************************************************
!# <name> io </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains several routines for input/output purposes.
!# </purpose>
!##############################################################################

MODULE IO

  USE fsystem
  USE error

  IMPLICIT NONE


CONTAINS

!************************************************************************************

!<subroutine>
  subroutine io_openFileForReading(sfilename, iunit)

    !<description>
    !This routine tries to open a file for reading. If succesful, on can read from it
    !via unit "iunit". Otherwise, iunit is -1.
    !</description>

    !<input>

    !filename
    character(*), intent(in) :: sfilename

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
      open(unit=iunit, file=trim(sfilename), iostat=istatus, action="read")
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
  subroutine io_openFileForWriting(sfilename, iunit, cflag, bfileExists)

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
    if (bexists .and. cflag .eq. SYS_REPLACE) then
      open(unit=iunit, file=trim(sfilename), iostat=istatus, status="replace", &
           action="write")
    else
      open(unit=iunit, file=trim(sfilename), iostat=istatus, action="write", &
           position="append")
    endif
    if (present(bfileExists)) then
      bfileExists = bexists
    endif
    if (istatus .ne. 0) then
      write(unit=*,fmt=*) "*** Error while opening file '", trim(sfilename), "'. ***"
      iunit = -1
    endif

  end subroutine io_openFileForWriting

END MODULE