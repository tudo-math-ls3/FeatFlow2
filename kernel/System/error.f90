!##############################################################################
!# ****************************************************************************
!# <name> error </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides a mapping for error IDs to hardcoded error messages.
!# Currently, this is basically only used by the io module. One could argue to
!# resolve this artificial splitting in the future.
!# </purpose>
!##############################################################################
module error

!$use omp_lib
  use genoutput
  use fsystem

  implicit none

  private

  !constants to use for subroutine error_print(...) indicating if the
  !error is critical (causing program exit).
  logical, parameter, public :: ERR_CRITICAL     = .true.
  logical, parameter, public :: ERR_NOT_CRITICAL = .false.

!********************************* 00 global ************************************
  ! No error
  integer, parameter, public :: ERR_NO_ERROR = 0000

  !Yet not implemented
  integer, parameter, public :: ERR_YNI = 0001

  !Conversion error string to int
  integer, parameter, public :: ERR_STRING_TO_INT = 0002

  !Conversion error string to real
  integer, parameter, public :: ERR_STRING_TO_REAL = 0003

  integer, parameter, public :: ERR_PARLST_ENV_VAR_NOT_FOUND = 0004

!*********************************** 14 io *****************************************
  integer, parameter, public :: ERR_IO_FILEIO             = 1401
  integer, parameter, public :: ERR_IO_WRONGSTRUCT        = 1402

  !empty file name
  integer, parameter, public :: ERR_IO_EMPTYFILENAME      = 1403

  !'mkdir -p' failed
  integer, parameter, public :: ERR_IO_MKDIR_P_FAILED     = 1404

  !no free unit found in subroutine sys_getFreeunit
  integer, parameter, public :: ERR_IO_NOFREEUNIT         = 1405

  !file not found error
  integer, parameter, public :: ERR_IO_NOSUCHFILE         = 1406

  !trying to export a non-assembled matrix to disk
  integer, parameter, public :: ERR_IO_MATRIX_UNASSEMBLED = 1407

!*********************************** 29 solver *************************************

  !matrix factorisation with UMFPACK failed
  integer, parameter, public :: ERR_SLV_UMFPACKFACTORISESYMB    = 2901
  integer, parameter, public :: ERR_SLV_UMFPACKFACTORISENUM     = 2902

  !solving the linear system with UMFPACK failed
  integer, parameter, public :: ERR_SLV_UMFPACKSOLVEFAILED      = 2903
  integer, parameter, public :: ERR_SLV_READGLOBSOLSIZEMISMATCH = 2904
  integer, parameter, public :: ERR_SLV_READGLOBSOLOUTOFRANGE   = 2905

  !invalid keyword in ScaRC solver definition
  integer, parameter, public :: ERR_SLV_INVALID_KEYWORD         = 2907

  !invalid solver type in ScaRC solver definition
  integer, parameter, public :: ERR_SLV_INVALID_SOLVER          = 2908

  !scarc solver file not found
  integer, parameter, public :: ERR_SLV_SCARCFILE_NOT_FOUND     = 2909
  !wrong chlayer
  integer, parameter, public :: ERR_SLV_WRONG_CHLAYER           = 2910

  !wrong chlayer
  integer, parameter, public :: ERR_SLV_GLOBALGRIDDEFECT        = 2911


!*********************************** 31 storage ************************************

  !memory allocation error
  integer, parameter, public :: ERR_ST_ALLOCF     = 3101
  !no descriptor block free
  integer, parameter, public :: ERR_ST_NODES      = 3102
  !wrong descriptor block
  integer, parameter, public :: ERR_ST_DESF       = 3103
  !wrong descriptor type
  integer, parameter, public :: ERR_ST_DESTF      = 3104
  !not enough memory
  integer, parameter, public :: ERR_ST_NOMEM      = 3105
  !request for negative or zero amount of memory
  integer, parameter, public :: ERR_ST_INVREQUEST = 3106
  !There cannot be more FREE storage blocks than the maximum number of blocks
  integer, parameter, public :: ERR_ST_FREE_BLOCK_ERROR = 3107

  ! last error occurred
  integer, private :: ierror = ERR_NO_ERROR

  public :: error_print
  public :: error_print_aux
  public :: error_clearError

contains

!************************************************************************

!<subroutine>
  subroutine error_print(&
       icode, sroutine, bcritical, &
       iarg1, iarg2, &
       darg1, darg2, &
       sarg1, sarg2)

  !<description>
    ! This routine acts as wrapper routine for error_print_aux. The preprocessor f90cpp
    ! will replace all occurences of this routine with calls to error_print_aux,
    ! inserting corrected line numbers etc.
  !</description>
  !<input>

    !error code
    integer, intent(in) :: icode

    !name of the calling routine
    character (len = *), intent(in) :: sroutine

    !flag if the error is critical, then terminate the program
    logical, intent(in) :: bcritical

    !integer argument 1 (optional)
    integer, optional :: iarg1

    !integer argument 2 (optional)
    integer, optional :: iarg2

    !double argument 1 (optional)
    real(DP), optional :: darg1

    !double argument 2 (optional)
    real(DP), optional :: darg2

    !string argument 1 (optional)
    character(len = *), optional :: sarg1

    !string argument 2 (optional)
    character(len = *), optional :: sarg2
  !</input>
!</subroutine>

    call error_print_aux(0, "-", icode, sroutine, bcritical, &
                         iarg1, iarg2, &
                         darg1, darg2, &
                         sarg1, sarg2)

  end subroutine error_print

!************************************************************************

!<subroutine>

  subroutine error_print_aux(&
       iline, sfile, icode, sroutine, bcritical, &
       iarg1, iarg2, &
       darg1, darg2, &
       sarg1, sarg2)

  !<description>
    ! This routine prints the error message for the given error code. If the occured
    ! error is critical so that the further program execution is impossible the program
    ! terminates.
    ! Even though this routine is public, users should never call it directly, because
    ! it is used by the preprocessor f90cpp.
  !</description>

  !<input>

    ! line of calling error_print
    integer, intent(in) :: iline

    ! file of calling error_print
    character(len=*), intent(in) :: sfile

    !error code
    integer, intent(in) :: icode

    !name of the calling routine
    character (len = *), intent(in) :: sroutine

    !flag if the error is critical, then terminate the program
    logical, intent(in) :: bcritical

    !integer argument 1 (optional)
    integer, optional :: iarg1

    !integer argument 2 (optional)
    integer, optional :: iarg2

    !double argument 1 (optional)
    real(DP), optional :: darg1

    !double argument 2 (optional)
    real(DP), optional :: darg2

    !string argument 1 (optional)
    character(len = *), optional :: sarg1

    !string argument 2 (optional)
    character(len = *), optional :: sarg2
  !</input>

!</subroutine>

    character(len = SYS_STRLEN) :: sstring, sstring1, sstring2, sprefix

    ! strings to construct the actual error message
    ! hard-coded indices, please obey:
    ! Smessage(1)  = delimiter
    ! Smessage(2)  = error code and routine
    ! Smessage(3)  = line number and file
    ! Smessage(4)  = empty line for improved layout
    ! Smessage(5)  = message body
    ! Smessage(6)  = message body
    ! Smessage(7)  = message body
    ! Smessage(8)  = message body
    ! Smessage(9)  = message body
    ! Smessage(10) = message body
    ! Smessage(11) = message body
    ! Smessage(12) = message body
    ! Smessage(13) = message body
    ! Smessage(14) = message body
    ! Smessage(15) = message body
    ! Smessage(16) = delimiter
    integer(I32), parameter :: NMESSAGELINES = 16
    character(len=SYS_STRLEN), dimension(NMESSAGELINES) :: Smessage
    character(len=SYS_STRLEN), dimension(NMESSAGELINES) :: Smessage2
    integer :: i,inextline

    ! initialise all message strings
    Smessage(:) = ""
    Smessage2(:) = ""


    ! Set module-private error code
    ierror = icode


    ! Error or warning only?
    if (bcritical) then
      sstring1 = "ERROR"
      sprefix = "[ERR]"
    else
      sstring1 = "WARNING"
      sprefix = "[WARN]"
    endif

    ! message header
    Smessage(1) = "*************************************" // &
                  "*************************************"

    ! Line information available?
    if (iline .ne. 0 ) then
      Smessage(2) = trim(sstring1) // " " // trim(sys_siL(icode, 5)) // &
                  " in <" // sroutine // ">,"
      Smessage(3) = "raised in line " // trim(sys_siL(iline, 5)) // &
                  " in file <" // sfile // ">:"
      Smessage(4) = ""
    else
      Smessage(2) = trim(sstring1) // " " // trim(sys_siL(icode, 5)) // &
                  " in <" // sroutine // ">:"
      Smessage(3) = ""
      Smessage(4) = ""
    endif



    select case (icode)
    case (ERR_YNI)
      if (present(sarg1)) then
        Smessage(5) = trim(sarg1)
        if (present(sarg2)) Smessage(6) = trim(sarg2)
      else
        Smessage(5) = "Yet not implemented."
      endif

    case (ERR_STRING_TO_INT)
      Smessage(5) = "Error while converting string <" // trim(sarg1) // "> to int."

    case (ERR_STRING_TO_REAL)
      Smessage(5) = "Error while converting string <" // trim(sarg1) // "> to real."

    case (ERR_PARLST_ENV_VAR_NOT_FOUND)
      Smessage(5) = "Environment variable <" // trim(sarg1) // "> not set."

!*********************************** 14 io *****************************************
    case(ERR_IO_NOFREEUNIT)
      Smessage(5) = "No free unit found. Not able to open the file."

    case (ERR_IO_FILEIO)
      Smessage(5) = "File input/output error."
      if (present(sarg1)) then
        Smessage(6) = sarg1
      endif

    case (ERR_IO_EMPTYFILENAME)
      if (present(sarg1)) then
        Smessage(5) = "File name <" // trim(sarg1) // "> empty."
      else
        Smessage(5) = "File name empty."
      endif

    case (ERR_IO_MKDIR_P_FAILED)
      if (present(sarg1)) then
        Smessage(5) = "Could not (auto-)create path <" // trim(sarg1) // ">."
      else
        Smessage(5) = "Could not (auto-)create path to output file."
      endif

    case (ERR_IO_WRONGSTRUCT)
      Smessage(5) = "error during read: wrong structure"

    case(ERR_IO_NOSUCHFILE)
      Smessage(5) = "File " // trim(sarg1) // " does not exist."
      Smessage(6) = "Read from file failed."

    case(ERR_IO_MATRIX_UNASSEMBLED)
      if (present(iarg1)) then
        sstring = sys_siL(iarg1,1)
      else
        sstring = "[unknown]"
      endif

      if (present(iarg2)) then
        sstring1 = sys_siL(iarg2,1)
      else
        sstring1 = "[unknown]"
      endif

      if (present(sarg1)) then
        sstring2 = sarg1
      else
        sstring2 = "[unknown]"
      endif
      Smessage(5) = "Block (" // trim(sstring) // "," // &
                    trim(sstring1) // ") of matrix <" // &
                    trim(sstring2) // ">"
      Smessage(6) = "seems not properly assembled. Export skipped."


!********************************* 31 storage **********************************

    case (ERR_ST_NODES)
      Smessage(5) = "FEAT ran out of storage descriptors. Please increase the value of"
      Smessage(6) = "STORAGE_DES in the configuration file."

    case (ERR_ST_DESF)
      Smessage(5) = "Wrong descriptor block."

    case (ERR_ST_DESTF)
      Smessage(5) = "Wrong descriptor type."

    case (ERR_ST_NOMEM)
      Smessage(5) = "Not enough memory on heap " // trim(sys_siL(iarg1, 6)) // &
                    ". Needs additionally " // trim(sys_siL(iarg2, 12)) // &
                    " entries! "

    case (ERR_ST_ALLOCF)
      Smessage(5) = "Memory allocation error."

    case (ERR_ST_INVREQUEST)
      Smessage(5) = "Request for invalid amount of memory:"
      Smessage(6) = "For " // trim(sarg1) // " the size " // &
                    trim(sys_siL(iarg1, 12)) // " was requested."

    case (ERR_ST_FREE_BLOCK_ERROR)
      Smessage(5) = "There cannot be more FREE storage blocks " // &
                    "than the maximum number of blocks."

!********************************* default **********************************

    case default
      Smessage(5) = "Unknown error code raised! " // trim(sys_siL(icode, 12))

    end select

    ! message footer
    Smessage(NMESSAGELINES) = "******************************************" // &
                              "********************************"


!    call output_line(OU_CLASS_ERROR, "", "******************************************" // &
!                                   "********************************")

    !Force warnings/errors to appear on screen/in log file.
    call sys_flush(OU_LOG)
    call sys_flush(6)

    ! pretty-print error message (warning: semi-hardcoded based on index ranges)
    do i = 1, 4
      Smessage2(i) = trim(sprefix) // " " // trim(Smessage(i))
    enddo
    inextLine = 5
    do i = 5, NMESSAGELINES-1
      if (Smessage(i) .ne. "") then
        Smessage2(inextLine) = trim(sprefix) // " " // trim(Smessage(i))
        inextLine = inextLine + 1
      endif
    enddo
    Smessage2(inextLine) = trim(sprefix) // " " // trim(Smessage(NMESSAGELINES))

    do i = 1, inextLine
!      write(SYS_LOG, '(a)') trim(Smessage2(i))
      call output_line(OU_CLASS_ERROR, "", trim(Smessage2(i)))
    enddo

    !In case the thrown error has been critical, end program now
    if (bcritical) then
      write(OU_LOG, '(a)') "Fatal exit!"

      !Write error code to screen
      sstring = "Error code " // trim(sys_siL(icode, 5))
      write (6, '(a)') trim(sstring)

      !Write error code to log file
      write (OU_LOG, '(a)') trim(sstring)

    endif

  end subroutine error_print_aux

!************************************************************************

!<function>
  integer function error_askError()

!<description>
!This function returns the error constant of the last occurred error.
!If no error has occurred ERR_NO_ERROR will be given back.
!Naturally, this works only for noncritical errors!
!</description>


!<result>
! error code
!</result>
!</function>
    error_askError = ierror
  end function error_askError

!<subroutine>
  subroutine error_clearError()

!<description>
!This routine resets the internal error memory to ERR_NO_ERROR.
!</description>

!</subroutine>

    ierror = ERR_NO_ERROR
  end subroutine error_clearError

end module error
