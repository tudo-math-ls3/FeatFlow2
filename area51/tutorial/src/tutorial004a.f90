!##############################################################################
!# Tutorial 004a: poststructures - sorting of arrays
!##############################################################################

module tutorial004a

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use sort

  implicit none
  private
  
  public :: start_tutorial004a

contains

  ! ***************************************************************************

  subroutine start_tutorial004a

    ! Declare some variables
    integer :: i
    integer, dimension(9) :: Ipost
    real(DP), dimension(9) :: Dpost

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 004a")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Sort an integer array
    ! =================================
    
    Ipost = (/  5, 2, 7, 3, 8, 4, 1, 9, 6  /)
    call sort_int (Ipost,SORT_QUICK)
    
    ! =================================
    ! Sort an double precision array
    ! =================================

    Dpost = (/  5.0_DP, 2.0_DP, 7.0_DP, 3.0_DP, 8.0_DP, 4.0_DP, 1.0_DP, 9.0_DP, 6.0_DP  /)
    call sort_dp (Dpost,SORT_QUICK)
    
    ! =================================
    ! Output
    ! =================================

    do i=1,9
      call output_line (" " // trim(sys_siL(Ipost(i),10)), bnolinebreak=(i .ne. 9) )
    end do

    do i=1,9
      call output_line (" " // trim(sys_sdL(Dpost(i),2)), bnolinebreak=(i .ne. 9) )
    end do

  end subroutine

end module
