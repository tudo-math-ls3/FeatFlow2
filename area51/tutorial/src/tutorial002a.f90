!##############################################################################
!# Tutorial 002a: Memory management, 1D arrays
!##############################################################################

module tutorial002a

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage

  implicit none
  private
  
  public :: start_tutorial002a

contains

  ! ***************************************************************************

  subroutine start_tutorial002a

    ! Declare some variables
    integer :: i
    integer :: ihandleInt, ihandleDouble
    
    integer, dimension(:), pointer :: p_Ipost
    real(DP), dimension(:), pointer :: p_Dpost

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 002a")
    call output_separator (OU_SEP_MINUS)
    
    ! ===============
    ! Allocate memory
    ! ===============
    
    ! 1D integer
    call storage_new ("start_tutorial002a", &
        "arrayI", 10, ST_INT, ihandleInt, ST_NEWBLOCK_ZERO)

    ! 1D double
    call storage_new ("start_tutorial002a", &
        "arrayI", 10, ST_DOUBLE, ihandleDouble, ST_NEWBLOCK_ZERO)
    
    ! =====================
    ! Get pointers
    ! =====================
    call storage_getbase_int (ihandleInt,p_Ipost)
    call storage_getbase_double (ihandleDouble,p_Dpost)

    ! =====================
    ! Fill with post
    ! =====================
    do i=1,10
      p_Ipost(i) = i**2
      p_Dpost(i) = real(i,DP) ** 2
    end do

    ! =====================
    ! Print the post
    ! =====================
    do i=1,10
      call output_line (trim(sys_siL(p_Ipost(i),10)))
    end do

    do i=1,10
      call output_line (trim(sys_sdEL(p_Dpost(i),10)))
    end do
    
    ! =====================
    ! Release the post.
    ! Pointers get invalid.
    ! =====================
    call storage_free (ihandleInt)
    call storage_free (ihandleDouble)

    ! =====================
    ! Print information
    ! about the memory.
    ! =====================
    call output_lbrk()
    call storage_info()
    
  end subroutine

end module
