!##############################################################################
!# Tutorial 003a: File IO
!##############################################################################

module tutorial003a

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage
  use io

  implicit none
  private
  
  public :: start_tutorial003a

contains

  ! ***************************************************************************

  subroutine start_tutorial003a

    ! Declare some variables
    integer :: iunit,ios,ilength
    character(LEN=SYS_STRLEN) :: sdata

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 003a")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Open a file, write, close.
    ! =================================
    call io_openFileForWriting("post/tutorial003a.txt", &
        iunit, SYS_REPLACE, bformatted=.true.)
    
    write (iunit,"(A)") "This is a"
    write (iunit,"(A)") "text file for"
    write (iunit,"(A)") "tutorial 003a"
    
    close (iunit)
    
    ! =================================
    ! Open the file, print, close.
    ! =================================
    call io_openFileForReading("post/tutorial003a.txt", iunit, .true.)
    
    ! With io_readlinefromfile, we can read line by line
    call io_readlinefromfile (iunit, sdata, ilength, ios)
    
    do while (ios .eq. 0)
      call output_line (sdata(1:ilength))
      call io_readlinefromfile (iunit, sdata, ilength, ios)
    end do
    
    close (iunit)
    
  end subroutine

end module
