   ! -*- mode: f90; -*-

  interface

!<subroutine>

  subroutine fcb_calcDivVector(rafcstab, rx, ry, dscale, bclear, rcollection)

    use afcstabilisation
    use collection
    use fsystem
    use linearsystemblock

!<description>
    ! This subroutine assembles the divergence flux vector
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! divergence vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

  end subroutine fcb_calcDivVector

  end interface