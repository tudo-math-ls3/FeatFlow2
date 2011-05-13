  ! -*- mode: f90; -*-

  interface

!<subroutine>

  subroutine fcb_limitNodal(NEQ, NVAR, dscale,&
        ML, Dpp, Dpm, Dqp, Dqm, Drp, Drm, rcollection)

    use fsystem
    use collection

!<description>
    ! This subroutine computes the nodal correction factors
!</description>

!<input>
    ! Number of nodes
    integer, intent(in) :: NEQ

    ! Number of variables
    integer, intent(in) :: NVAR

    ! Sums of antidiffusive increments
    real(DP), dimension(NVAR,NEQ), intent(in) :: Dpp,Dpm

    ! Local bounds
    real(DP), dimension(NVAR,NEQ), intent(in) :: Dqp,Dqm

    ! Lumped mass matrix
    real(DP), dimension(:), intent(in) :: ML

    ! Scaling parameter
    real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! Nodal correction factors
    real(DP), dimension(NVAR,NEQ), intent(inout) :: Drp,Drm

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

  end subroutine fcb_limitNodal

  end interface