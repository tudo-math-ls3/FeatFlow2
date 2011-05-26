  ! -*- mode: f90; -*-

  interface

!<subroutine>

  subroutine fcb_calcCorrection(IverticesAtEdgeIdx, IverticesAtEdge, NEDGE,&
      NEQ, NVAR, dscale, Dalpha, Dflux, ndim1, ndim2, Dy, ML, rcollection)

    use fsystem
    use collection

!<description>
    ! This subroutine applied the corrected antidiffusive fluxes
!</description>

!<input>
    ! Number of edges
    integer, intent(in) :: NEDGE

    ! Number of nodes
    integer, intent(in) :: NEQ

    ! Number of variables
    integer, intent(in) :: NVAR

    ! Dimensions of the solution vector
    integer, intent(in) :: ndim1, ndim2

    ! Edgewise correction factors
    real(DP), dimension(:), intent(in) :: Dalpha

    ! Raw antidiffusive fluxes
    real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Index pointer for edge data structure
    integer, dimension(:), intent(in) :: IverticesAtEdgeIdx

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! OPTIONAL :: Lumped mass matrix
    real(DP), dimension(:), intent(in), optional :: ML
!</input>

!<inputoutput>
    ! Vector to which the correction is applied
    real(DP), dimension(ndim1,ndim2), intent(inout) :: Dy

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

  end subroutine fcb_calcCorrection

  end interface