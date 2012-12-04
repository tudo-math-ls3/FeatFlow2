!##############################################################################
!# ****************************************************************************
!# <name> preprocstd </name>
!# ****************************************************************************
!#
!# <purpose>
!# Standard preprocessor commands. Support for ASSERT and ASSUME.
!#
!# </purpose>
!##############################################################################

module preprocstd

  use fsystem
  use genoutput
  use io
  
  private
  
  public assert_internal
  
contains

  subroutine assert_internal (bcondition,iline,sfile)
  
    logical, intent(in) :: bcondition
    integer, intent(in) :: iline
    character(len=*), intent(in) :: sfile
    character(len=80) :: sname
  
    if (.not. bcondition) then
      call io_pathExtract (sfile, sfilename=sname)
      call output_line ("Assumption failed; line "//trim(sys_siL(iline,10)), &
          OU_CLASS_ERROR,OU_MODE_STD,trim(sname))
      call sys_halt()
    end if

  end subroutine

end module
