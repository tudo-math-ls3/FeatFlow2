! -*- mode: f90; -*-

  ! This file describes the interfaces of the callback routines of the
  ! boundary treatment routines.

  ! ***************************************************************************

  interface

!<subroutine>

  subroutine fcb_calcBoundaryvalues(DbdrNormal, DpointNormal,&
      DbdrValue, ibdrCondType, Du, Du0, istatus)

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

  ! ***************************************************************************

  interface

!<subroutine>

  subroutine fcb_getBdrCondExprNumber(ibdrCondType, ndimension, nexpr)

!<description>
    ! This subroutine calculates the number of mathematical expressions
    ! required for boundary conditions of type ibdrCondType in the
    ! given spatial dimension ndimension.
!</description>

!<input>
    ! type of boundary condition
    integer, intent(in) :: ibdrCondType

    ! number of spatial dimensions
    integer, intent(in) :: ndimension
!</input>

!<output>
    integer, intent(out) :: nexpr
!</outpu>
!</subroutine>

  end subroutine fcb_getBdrCondExprNumber

  end interface
