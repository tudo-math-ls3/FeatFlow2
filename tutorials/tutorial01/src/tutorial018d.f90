!##############################################################################
!# Tutorial 018d: Collection handling. Using levels.
!##############################################################################

module tutorial018d

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput

  use collection

  implicit none
  private
  
  public :: start_tutorial018d

contains

  ! ***************************************************************************

  subroutine printarray (rcollection,ilevel)
  type(t_collection), intent(inout) :: rcollection
  integer, intent(in) :: ilevel
  
    integer, dimension(:), pointer :: p_Iarray
    real(DP), dimension(:), pointer :: p_Darray
    integer :: i
    
    call output_line ("Level: " // sys_si(ilevel,3))
    
    ! Get a pointer to the data from the collection
    call collct_getvalue_iarrp (rcollection, "MYIARRAY", p_Iarray, ilevel)
    call collct_getvalue_rarrp (rcollection, "MYDARRAY", p_Darray, ilevel)
    
    ! Print the content
    do i=1,ubound(p_Iarray,1)
      call output_line (sys_si(p_Iarray(i),5))
    end do

    do i=1,ubound(p_Darray,1)
      call output_line (sys_sd(p_Darray(i),5))
    end do
  
  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial018d

    ! Declare some variables.
    type(t_collection) :: rcollection
    integer, dimension(:), pointer :: p_Iarray
    real(DP), dimension(:), pointer :: p_Darray
    integer :: ilevel
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 018d")
    call output_separator (OU_SEP_MINUS)

    ! For storing named items in the collection, the collection must be
    ! initialised.
    !
    ! Initialise the collection
    call collct_init (rcollection)
    
    ! Create a hierarchy of 3 levels
    do ilevel = 1,3
      call collct_addlevel (rcollection)
    end do
    
    ! Put some data to all levels
    do ilevel = 1,3
      allocate (p_Iarray(2**ilevel))
      p_Iarray(:) = ilevel*2
      
      allocate (p_Darray(2**ilevel))
      p_Darray(:) = 1.0_DP/real(ilevel,DP)
      
      ! Put the data arrays to the collection at level ilevel
      call collct_setvalue_iarrp (rcollection, "MYIARRAY", p_Iarray, .true.,&
          ilevel)
    
      call collct_setvalue_rarrp (rcollection, "MYDARRAY", p_Darray, .true.,&
          ilevel)

    end do
    
    ! Pass the collection to the subroutine
    do ilevel = 1,3
      call printarray (rcollection, ilevel)
    end do
    
    ! Remove the data
    do ilevel = 1,3
    
      ! Get the pointers, 
      call collct_getvalue_iarrp (rcollection, "MYIARRAY", p_Iarray, ilevel)
      call collct_getvalue_rarrp (rcollection, "MYDARRAY", p_Darray, ilevel)
      
      ! Remove the data from the collection
      call collct_deletevalue (rcollection,"MYIARRAY",ilevel)
      call collct_deletevalue (rcollection,"MYDARRAY",ilevel)

      ! Deallocate
      deallocate(p_Iarray)
      deallocate(p_Darray)

    end do

    ! Print statistics of the collection
    call collct_printStatistics (rcollection)

    ! Release the collection
    call collct_done (rcollection)
    
  end subroutine

end module
