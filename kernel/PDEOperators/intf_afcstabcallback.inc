! -*- mode: f90; -*-
  
  ! This file describes the interfaces of the callback routines
  ! for the discrete stabilisation by means of algebraic flux correction

  ! ***************************************************************************

  interface

  ! ***************************************************************************

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

  ! ***************************************************************************

!<subroutine>

  subroutine fcb_extractVariableBlock(rvectorSrc, cvariable, rvectorDest)

  use linearsystemscalar
  use linearsystemblock

!<description>
    ! This subroutine extracts a single variable identified by the
    ! string cvariable from the block vector rvectorSrc and stores
    ! the result in the scalar destination vector rvectorDest.
!</description>

!<input>
    ! block vector
    type(t_vectorBlock), intent(in) :: rvectorSrc

    ! identifier for the variable
    character(LEN=*), intent(in) :: cvariable
!</input>

!<inputoutput>
    ! scalar vector
    type(t_vectorScalar), intent(inout) :: rvectorDest
!</inputoutput>
!</subroutine>

  end subroutine fcb_extractVariableBlock

  ! ***************************************************************************

!<subroutine>

  subroutine fcb_extractVariableArray(RvectorSrc, cvariable, rvectorDest)

  use linearsystemscalar
  use linearsystemblock

!<description>
    ! This subroutine extracts a single variable identified by the
    ! string cvariable from the array of block vectors rvectorSrc and
    ! stores the result in the scalar destination vector rvectorDest.
!</description>

!<input>
    ! array of block vector
    type(t_vectorBlock), dimension(:), intent(in) :: RvectorSrc

    ! identifier for the variable
    character(LEN=*), intent(in) :: cvariable
!</input>

!<inputoutput>
    ! scalar vector
    type(t_vectorScalar), intent(inout) :: rvectorDest
!</inputoutput>
!</subroutine>

  end subroutine fcb_extractVariableArray

  end interface
