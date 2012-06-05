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

  subroutine fcb_parseBoundaryCondition(cbdrCondType, ndimension,&
      ibdrCondType, nexpressions)

!<description>
    ! This subroutine parses the boundarry condition defined in the
    ! parameter file. That is, the string cbdrCondType is tansformed
    ! into an integer value ibdrCondType and the number of
    ! mathematical expressions corresponding to the given boundary
    ! type and the spatial dimension ndimension are returned.
!</description>

!<input>
    ! character string: type of boundary conditions
    character(len=*), intent(in) :: cbdrCondType

    ! number of spatial dimensions
    integer, intent(in) :: ndimension
!</input>

!<output>
    ! type of boundary condition
    integer, intent(out) :: ibdrCondType

    ! number of mathematical expressions
    integer, intent(out) :: nexpressions
!</outpu>
!</subroutine>

  end subroutine fcb_parseBoundaryCondition

  end interface
