!##############################################################################
!# Tutorial 018e: Collection handling. Using sections.
!##############################################################################

module tutorial018e

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput

  use collection

  implicit none
  private
  
  public :: start_tutorial018e

contains

  ! ***************************************************************************

  subroutine printarray (rcollection,ssection)
  type(t_collection), intent(inout) :: rcollection
  character(len=*), intent(in) :: ssection
  
    integer, dimension(:), pointer :: p_Iarray
    real(DP), dimension(:), pointer :: p_Darray
    integer :: i
    
    call output_line ("Section: " // ssection)
    
    ! Get a pointer to the data from the collection
    call collct_getvalue_iarrp (rcollection, "MYIARRAY", p_Iarray, ssectionName=ssection)
    call collct_getvalue_rarrp (rcollection, "MYDARRAY", p_Darray, ssectionName=ssection)
    
    ! Print the content
    do i=1,ubound(p_Iarray,1)
      call output_line (sys_si(p_Iarray(i),5))
    end do

    do i=1,ubound(p_Darray,1)
      call output_line (sys_sd(p_Darray(i),5))
    end do
  
  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial018e

    ! Declare some variables.
    type(t_collection) :: rcollection
    integer, dimension(5), target :: Iarray1, Iarray2
    real(DP), dimension(5), target :: Darray1, Darray2
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 018e")
    call output_separator (OU_SEP_MINUS)

    ! For storing named items in the collection, the collection must be
    ! initialised.
    !
    ! Initialise the collection
    call collct_init (rcollection)
    
    ! Add two sections.
    call collct_addsection (rcollection,"SECTION1")
    call collct_addsection (rcollection,"SECTION2")
    
    ! Initialise arrays
    Iarray1(:) = (/ 1,2,3,4,5 /)
    Darray1(:) = (/ 5.0_DP,4.0_DP,3.0_DP,2.0_DP,1.0_DP /)

    Iarray2(:) = (/ 6,7,8,9,10 /)
    Darray2(:) = (/ 10.0_DP,9.0_DP,8.0_DP,7.0_DP,6.0_DP /)
    
    ! Save pointers to the arrays in the collection.
    !
    ! Section 1
    call collct_setvalue_iarrp (rcollection, "MYIARRAY", Iarray1, .true.,&
        ssectionName = "SECTION1")
    
    call collct_setvalue_rarrp (rcollection, "MYDARRAY", Darray1, .true.,&
        ssectionName = "SECTION1")
    
    ! Section 2
    call collct_setvalue_iarrp (rcollection, "MYIARRAY", Iarray2, .true.,&
        ssectionName = "SECTION2")
    
    call collct_setvalue_rarrp (rcollection, "MYDARRAY", Darray2, .true.,&
        ssectionName = "SECTION2")

    ! Pass the collection to the subroutine
    call printarray (rcollection, "SECTION1")
    call printarray (rcollection, "SECTION2")
    
    ! Remove the data
    call collct_deletevalue (rcollection,"MYIARRAY",ssectionName = "SECTION1")
    call collct_deletevalue (rcollection,"MYDARRAY",ssectionName = "SECTION1")
    
    call collct_deletevalue (rcollection,"MYIARRAY",ssectionName = "SECTION2")
    call collct_deletevalue (rcollection,"MYDARRAY",ssectionName = "SECTION2")
    
    ! Remove the sections
    call collct_deletesection (rcollection, "SECTION1")
    call collct_deletesection (rcollection, "SECTION2")
    
    ! Print statistics of the collection
    call collct_printStatistics (rcollection)

    ! Release the collection
    call collct_done (rcollection)
    
  end subroutine

end module
