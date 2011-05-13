  ! -*- mode: f90; -*-

  interface

!<subroutine>

  subroutine fcb_extractVariableScalar(rvectorSrc, cvariable, rvectorDest)

  use linearsystemscalar

!<description>
    ! This subroutine extracts a single variable identified by the
    ! string cvariable from the scalar vector rvectorSrc and stores
    ! the result in the scalar destination vector rvectorDest.
!</description>

!<input>
    ! block vector
    type(t_vectorScalar), intent(in) :: rvectorSrc

    ! identifier for the variable
    character(LEN=*), intent(in) :: cvariable
!</input>

!<inputoutput>
    ! scalar vector
    type(t_vectorScalar), intent(inout) :: rvectorDest
!</inputoutput>
!</subroutine>

  end subroutine fcb_extractVariableScalar

  end interface