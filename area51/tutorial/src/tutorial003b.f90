!##############################################################################
!# Tutorial 003b: Reading of post files
!##############################################################################

module tutorial003b

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage
  use paramlist

  implicit none
  private
  
  public :: start_tutorial003b

contains

  ! ***************************************************************************

  subroutine start_tutorial003b

    ! Declare some variables
    type(t_parlist) :: rparlist
    integer :: ipost, i,n
    real(DP) :: dpost
    character(LEN=SYS_STRLEN) :: spost

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 003b")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Read the post file.
    ! =================================
    call parlst_init(rparlist)
    call parlst_readfromfile (rparlist, "data/tutorial003b_data.dat")

    ! =================================
    ! Print the unnamed section
    ! =================================
    call parlst_getvalue_int (rparlist,"","value1",ipost)
    call parlst_getvalue_double (rparlist,"","value2",dpost)
    
    call output_line (sys_siL(ipost,5))
    call output_line (sys_sdEL(dpost,10))

    ! =================================
    ! Print the section 1
    ! =================================
    call parlst_getvalue_int (rparlist,"PARAMETERSET1","value1",ipost)
    call parlst_getvalue_double (rparlist,"PARAMETERSET1","value2",dpost)
    call parlst_getvalue_string (rparlist,"PARAMETERSET1","value3",spost,bdequote=.true.)
    ! bdequote dequotes the string automatically...
    
    call output_line (sys_siL(ipost,5))
    call output_line (sys_sdEL(dpost,10))
    call output_line (trim(spost))

    ! =================================
    ! Print the section 2
    ! =================================
    n = parlst_querysubstrings (rparlist,"PARAMETERSET2","valuelist1")
    
    do i=1,n
      call parlst_getvalue_double (&
          rparlist,"PARAMETERSET2","valuelist1",dpost, iarrayindex = i)
      call output_line (sys_sdEL(dpost,10))
    end do

    ! =================================
    ! Cleanup
    ! =================================
    call parlst_done(rparlist)
    
  end subroutine

end module
