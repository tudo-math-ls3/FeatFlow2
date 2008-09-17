!########################################################################
!# FINITE ELEMENT ANALYSIS & SOLUTION TOOLS  F E A S T  (Release 1.0)   #
!#                                                                      #
!# Authors: Ch.Becker, S.Kilian, S.Turek                                #
!#          Institute of Applied Mathematics & Simulation               #
!#          University of Dortmund                                      #
!#          D-44227 DORTMUND                                            #
!#                                                                      #
!########################################################################
!#                                                                      #
!# <name> error </name>                                                 #
!#                                                                      #
!#                                                                      #
!# <purpose>                                                            #
!#                                                                      #
!# This module contains the definition of the error codes and the       #
!# routine which prints the output text if an error has ocuured.        #
!# </purpose>                                                           #
!#                                                                      #
!########################################################################

!# Current version: $Id: error.f90,v 1.98 2006/03/21 16:29:28 goeddeke4 Exp $

!<!--
!#include <feastdefs.h>
! -->

module error

  use genoutput
  use fsystem

  !constants to use for subroutine error_print(...) indicating if the
  !error is critical (causing program exit).
  logical, parameter :: ERR_CRITICAL     = .true.
  logical, parameter :: ERR_NOT_CRITICAL = .false.

!********************************* 00 global ************************************
  ! No error
  integer, parameter :: ERR_NO_ERROR = 0000
  
  !Yet not implemented
  integer, parameter :: ERR_YNI = 0001

  !Conversion error string to int
  integer, parameter :: ERR_STRING_TO_INT = 0002

  !Conversion error string to real
  integer, parameter :: ERR_STRING_TO_REAL = 0003


!*********************************** 14 io *****************************************
  integer, parameter :: ERR_IO_FILEIO             = 1401
  integer, parameter :: ERR_IO_WRONGSTRUCT        = 1402

  !empty file name
  integer, parameter :: ERR_IO_EMPTYFILENAME      = 1403

  !no free unit found in subroutine sys_getFreeunit
  integer, parameter :: ERR_IO_NOFREEUNIT         = 1404

  !file not found error
  integer, parameter :: ERR_IO_NOSUCHFILE         = 1405

  !trying to export a non-assembled matrix to disk
  integer, parameter :: ERR_IO_MATRIX_UNASSEMBLED = 1406

!*********************************** 29 solver *************************************

  !matrix factorisation with UMFPACK failed
  integer, parameter :: ERR_SLV_UMFPACKFACTORISESYMB    = 2901
  integer, parameter :: ERR_SLV_UMFPACKFACTORISENUM     = 2902

  !solving the linear system with UMFPACK failed
  integer, parameter :: ERR_SLV_UMFPACKSOLVEFAILED      = 2903
  integer, parameter :: ERR_SLV_READGLOBSOLSIZEMISMATCH = 2904
  integer, parameter :: ERR_SLV_READGLOBSOLOUTOFRANGE   = 2905

  !invalid keyword in ScaRC solver definition
  integer, parameter :: ERR_SLV_INVALID_KEYWORD         = 2907

  !invalid solver type in ScaRC solver definition
  integer, parameter :: ERR_SLV_INVALID_SOLVER          = 2908

  !scarc solver file not found
  integer, parameter :: ERR_SLV_SCARCFILE_NOT_FOUND     = 2909
  !wrong chlayer
  integer, parameter :: ERR_SLV_WRONG_CHLAYER           = 2910

  !wrong chlayer
  integer, parameter :: ERR_SLV_GLOBALGRIDDEFECT        = 2911


!*********************************** 31 storage ************************************

  !memory allocation error
  integer, parameter :: ERR_ST_ALLOCF     = 3101
  !no descriptor block free
  integer, parameter :: ERR_ST_NODES      = 3102
  !wrong descriptor block
  integer, parameter :: ERR_ST_DESF       = 3103
  !wrong descriptor type
  integer, parameter :: ERR_ST_DESTF      = 3104
  !not enough memory
  integer, parameter :: ERR_ST_NOMEM      = 3105
  !request for negative or zero amount of memory
  integer, parameter :: ERR_ST_INVREQUEST = 3106
  !There cannot be more FREE storage blocks than the maximum number of blocks
  integer, parameter :: ERR_ST_FREE_BLOCK_ERROR = 3107

  ! last error occurred
  integer, private :: ierror = ERR_NO_ERROR

contains

!************************************************************************
!<subroutine>
  subroutine error_print(icode, sroutine, bcritical, iarg1, iarg2, darg1, darg2, &
                         sarg1, sarg2)
!<description>
!This routine prints the error message for the given error code. If the occured error is
!critical so that the further program execution is impossible the program terminates.
!</description>

!<input>

    !name of the calling routine
    character (len = *), intent(in) :: sroutine

    !flag if the error is critical, then terminate the program
    logical, intent(in) :: bcritical

    !error code
    integer, intent(in) :: icode

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
! </input>

!</subroutine>

    character(len = SYS_STRLEN) :: sstring, sstring1, sstring2
    
    ierror = icode
    
    call output_line(OU_CLASS_ERROR, "", "******************************************" // &
         "********************************")
    if (bcritical) then
      call output_line(OU_CLASS_ERROR, "", &
                       "ERROR " // trim(sys_siL(icode, 5)) // " in '" // sroutine // "': ")
    else
      call output_line(OU_CLASS_ERROR, "", &
                       "WARNING " // trim(sys_siL(icode, 5)) // " in '" // sroutine // "': ")
    endif

    select case (icode)
    case (ERR_YNI)
      if (present(sarg1)) then
        call output_line(OU_CLASS_ERROR, "", trim(sarg1))
      else
        write(OU_LOG, '(A)') "Yet not implemented."
      end if

    case (ERR_STRING_TO_INT)
      call output_line(OU_CLASS_ERROR, "", "Error while converting string '" // trim(sarg1) // &
                                     "' to int.")

    case (ERR_STRING_TO_REAL)
      call output_line(OU_CLASS_ERROR, "", "Error while converting string '" // trim(sarg1) // &
                                     "' to real.")

! Case never used in FEAST
!    case (ERR_ST_NOMB)
!      write(OU_LOG, '(A)') "No memory block free."

! Case never used in FEAST
!    case (ERR_ST_HEAPF)
!      write(OU_LOG, '(A)') "Wrong heap descriptor."

!*********************************** 14 io *****************************************
    case(ERR_IO_NOFREEUNIT)
      call output_line(OU_CLASS_ERROR, "", "No free unit found, not able to open the file.")

    case (ERR_IO_FILEIO)
      call output_line(OU_CLASS_ERROR, "", "File input/output error.")
      if (present(sarg1)) then
        call output_line(OU_CLASS_ERROR, "", sarg1)
      endif

    case (ERR_IO_EMPTYFILENAME)
      if (present(sarg1)) then
        call output_line(OU_CLASS_ERROR, "", "File name '" // trim(sarg1) // "' empty.")
      else
        call output_line(OU_CLASS_ERROR, "", "File name empty.")
      endif

    case (ERR_IO_WRONGSTRUCT)
      write(OU_LOG, '(A)') "error during read: wrong structure"

    case(ERR_IO_NOSUCHFILE)
      call output_line(OU_CLASS_ERROR, "", "File " // trim(sarg1) // " does not exist.")
      call output_line(OU_CLASS_ERROR, "", "Read from file failed.")

    case(ERR_IO_MATRIX_UNASSEMBLED)
      if (present(iarg1)) then
        sstring = sys_siL(iarg1,1)
      else
        sstring = "[unknown]"
      end if

      if (present(iarg1)) then
        sstring1 = sys_siL(iarg2,1)
      else
        sstring1 = "[unknown]"
      end if

      if (present(sarg1)) then
        sstring2 = sarg1
      else
        sstring2 = "[unknown]"
      end if
      call output_line(OU_CLASS_ERROR, "", "Block (" // trim(sstring) // "," // &
                                     trim(sstring1) // ") of matrix '" // &
                                     trim(sstring2) // "' seems not properly assembled.")
      call output_line(OU_CLASS_ERROR, "", "Export skipped.")


!********************************* 31 storage **********************************

    case (ERR_ST_NODES)
      call output_line(OU_CLASS_ERROR, "", "No descriptor block free.")

    case (ERR_ST_DESF)
      call output_line(OU_CLASS_ERROR, "", "Wrong descriptor block.")

    case (ERR_ST_DESTF)
      call output_line(OU_CLASS_ERROR, "", "Wrong descriptor type.")

    case (ERR_ST_NOMEM)
      call output_line(OU_CLASS_ERROR, "", "Not enough memory on heap " // &
                                     trim(adjustl(sys_i6(iarg1))) // &
                                     ". Needs additionally " // &
                                     trim(adjustl(sys_i12(iarg2))) // &
                                     " entries! Aborting program.")

    case (ERR_ST_ALLOCF)
      call output_line(OU_CLASS_ERROR, "", "Memory allocation error.")

    case (ERR_ST_INVREQUEST)
      call output_line(OU_CLASS_ERROR, "", "Request for invalid amount of memory:")
      call output_line(OU_CLASS_ERROR, "", "For " // sarg1 // " the size " // &
                                     trim(adjustl(sys_i12(iarg1))) // &
                                     " was requested.")
    case (ERR_ST_FREE_BLOCK_ERROR)
      call output_line(OU_CLASS_ERROR, "", "There cannot be more FREE storage blocks " // &
                                     "than the maximum number of blocks.")

!********************************* default **********************************

    case default
      call output_line(OU_CLASS_ERROR, "", "Unknown error code raised! " // sys_i12(icode))

    end select



    call output_line(OU_CLASS_ERROR, "", "******************************************" // &
                                   "********************************")

    !Force warnings/errors to appear on screen/in log file.
    call sys_flush(OU_LOG)
    call sys_flush(6)

    !In case the thrown error has been critical, end program now
    if (bcritical) then
      write(OU_LOG, '(a)') "Fatal exit!"

      !Write error code to screen
      sstring = "Error code " // trim(sys_siL(icode, 5))
      write (6, '(a)') trim(sstring)

      !Write error code to log file
      write (OU_LOG, '(a)') trim(sstring)

    endif

  end subroutine error_print
  
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
