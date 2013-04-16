!##############################################################################
!# Tutorial 002b: Memory management, 2D arrays
!##############################################################################

module tutorial002b

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage

  implicit none
  private
  
  public :: start_tutorial002b

contains

  ! ***************************************************************************

  subroutine start_tutorial002b

    ! Declare some variables
    integer :: i,j
    integer :: ihandleInt2d, ihandleDouble2d
    integer, dimension(2) :: Isize
    
    integer, dimension(:,:), pointer :: p_Ipost
    real(DP), dimension(:,:), pointer :: p_Dpost

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 002b")
    call output_separator (OU_SEP_MINUS)

    ! ===============    
    ! Allocate memory
    ! ===============
    Isize = (/ 4,4 /)
    
    ! 2D integer; use "ST_NEWBLOCK_NOINIT" to not initialise with zero.
    call storage_new ("start_tutorial002b", &
        "arrayI", Isize, ST_INT, ihandleInt2D, ST_NEWBLOCK_ZERO)

    ! 2D double; use "ST_NEWBLOCK_NOINIT" to not initialise with zero.
    call storage_new ("start_tutorial002b", &
        "arrayI", Isize, ST_DOUBLE, ihandleDouble2D, ST_NEWBLOCK_ZERO)

    ! =====================
    ! Get pointers
    ! =====================
    call storage_getbase_int2D (ihandleInt2D,p_Ipost)
    call storage_getbase_double2D (ihandleDouble2D,p_Dpost)

    ! =====================
    ! Fill with post
    ! =====================
    do j=1,4
      do i=1,4
        p_Ipost(i,j) = i**2 * j
        p_Dpost(i,j) = real(i**2 * j,DP)
      end do
    end do

    ! =====================
    ! Print the post
    ! =====================
    do i=1,4
      do j=1,4
        call output_line (sys_si(p_Ipost(i,j),5), bnolinebreak=(j .ne. 4) )
      end do
    end do

    do i=1,4
      do j=1,4
        call output_line (sys_adjustr(sys_sdE(p_Dpost(i,j),2),10), bnolinebreak=(j .ne. 4) )
      end do
    end do
    
    ! =====================
    ! Release the post.
    ! Pointers get invalid.
    ! =====================
    call storage_free (ihandleInt2D)
    call storage_free (ihandleDouble2D)

    ! =====================
    ! Print information
    ! about the memory.
    ! =====================
    call output_lbrk()
    call storage_info()

  end subroutine

end module
