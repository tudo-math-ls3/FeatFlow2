  interface

!<subroutine>

  subroutine fcb_calcADFluxes(IedgeListIdx, IedgeList, NEDGE, NEQ,&
      NVAR, NVARtransformed, ndim1, ndim2, dscale, Dx, Dflux, rcollection)

    use fsystem
    use collection

!<description>
    ! This subroutine assembles the antidiffusive fluxes
!</description>

!<input>
    ! Number of edges
    integer, intent(in) :: NEDGE

    ! Number of nodes
    integer, intent(in) :: NEQ

    ! Number of solution variables
    integer, intent(IN) :: NVAR

    ! Number of transformed variables
    integer, intent(IN) :: NVARtransformed

    ! Dimensions of the solution vector
    integer, intent(in) :: ndim1, ndim2

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Solution used for flux calculation
    real(DP), dimension(ndim1,ndim2), intent(in) :: Dx

    ! Index pointer for edge data structure
    integer, dimension(:), intent(in) :: IedgeListIdx

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IedgeList
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Raw antidiffusive flux
    real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux
!</output>
!</subroutine

  end subroutine fcb_calcADFluxes

  end interface
