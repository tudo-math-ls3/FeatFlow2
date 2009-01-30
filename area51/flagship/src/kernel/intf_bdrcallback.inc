  ! This file describes the interfaces of the callback routines of the
  ! boundary treatment routines.  
  
  ! ***************************************************************************

  interface

!<subroutine>

  subroutine fcb_calcBoundaryvalues(DbdrNormal, DpointNormal,&
                                    DbdrValue, ibdrCondType,&
                                    Du, Du0, istatus)

  use fsystem

!<description>
    ! This subroutine computes the boundary values
!</description>

!<input>
  real(DP), dimension(:), intent(IN) :: DbdrNormal
  real(DP), dimension(:), intent(IN) :: DpointNormal
  real(DP), dimension(:), intent(IN) :: DbdrValue
  real(DP), dimension(:), intent(IN) :: Du0
  integer, intent(IN) :: ibdrCondType
!</input>

!<inputoutput>
  real(DP), dimension(:), intent(INOUT) :: Du
  integer, intent(INOUT), optional :: istatus
!</inputoutput>

!</subroutine>

  end subroutine fcb_calcBoundaryvalues

  end interface
