  interface

!<subroutine>

  subroutine fcb_calcFailsafe(IedgeListIdx, IedgeList,&
      NEDGE, NEQ, NVAR, ndim1, ndim2, ioperationSpec, dscale, dtol,&
      ML, Dx, Dalpha, Dflux, bisAccepted, Dbeta, rcollection)

!<description>
    ! This subroutine calculates the failsafe correction factors
!</description>

    use fsystem
    use collection

!<input>
    ! Number of edges
    integer, intent(in) :: NEDGE

    ! Number of nodes
    integer, intent(in) :: NEQ

    ! Number of variables
    integer, intent(in) :: NVAR

    ! Dimensions of the solution vector
    integer, intent(in) :: ndim1, ndim2

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Tolerance parameter
    real(DP), intent(in) :: dtol

    ! Operation specification tag. This is a bitfield coming from an
    ! OR combination of different AFCSTAB_FAILSAFE_xxxx constants and
    ! specifies which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! Lumped mass matrix
    real(DP), dimension(:), intent(in) :: ML

    ! Index pointer for edge data structure
    integer, dimension(:), intent(in) :: IedgeListIdx

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IedgeList
!</input>

!<inputoutput>
    ! Solution vector
    real(DP), dimension(ndim1,ndim2), intent(inout) :: Dx

    ! Edgewise correction factors
    real(DP), dimension(:), intent(inout) :: Dalpha

    ! Raw antidiffusive fluxes
    real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux

    ! OPTIONAL: Edgewise failsafe correction factors
    real(DP), dimension(:), intent(inout), optional :: Dbeta

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Flag which is TRUE if the computed solution vector
    ! is accepted by the failsafe correction procedure
    logical, intent(out) :: bisAccepted
!</output>
!</subroutine>

  end subroutine fcb_calcFailsafe

  end interface
