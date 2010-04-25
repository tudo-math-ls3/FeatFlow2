! -*- mode: f90; -*-
  
  ! This file describes the interfaces of the callback routines
  ! for the discrete stabilisation by means of algebraic flux correction

  ! ***************************************************************************

  interface

  ! ***************************************************************************

!<subroutine>

  subroutine fcb_extractVariable(rvectorBlock, cvariable, rvectorScalar)

  use linearsystemscalar
  use linearsystemblock

!<description>
    ! This subroutine extracts a single variable identified by the
    ! string cvariable from the block vector rvector.
!</description>

!<input>
    ! block vector
    type(t_vectorBlock), intent(in) :: rvectorBlock

    ! identifier for the variable
    character(LEN=*), intent(in) :: cvariable
!</input>

!<inputoutput>
    ! scalar vector
    type(t_vectorScalar), intent(inout) :: rvectorScalar
!</inputoutput>
!</subroutine>

  end subroutine fcb_extractVariable

  end interface
