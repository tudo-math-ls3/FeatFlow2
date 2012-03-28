  interface

!<subroutine>

  subroutine fcb_extractVariable(rvectorSrc, cvariable, rvectorDest, Imask)

  use linearsystemscalar
  use linearsystemblock

!<description>
    ! This subroutine extracts a single variable identified by the
    ! string cvariable from the block vector rvectorSrc and stores
    ! the result in the scalar destination vector rvectorDest.
!</description>

!<input>
    ! source block vector
    type(t_vectorBlock), intent(in) :: rvectorSrc

    ! identifier for the variable
    character(LEN=*), intent(in) :: cvariable

    ! OPTIONAL: integer mask array
    ! If present only those entries of the destination vector are
    ! computed which are given by the integer mask.
    integer, dimension(:), intent(in), optional :: Imask
!</input>

!<inputoutput>
    ! scalar destination vector
    type(t_vectorScalar), intent(inout) :: rvectorDest
!</inputoutput>
!</subroutine>

  end subroutine fcb_extractVariable

  end interface
