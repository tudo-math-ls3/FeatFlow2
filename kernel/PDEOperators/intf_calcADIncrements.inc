  interface

!<subroutine>

  subroutine fcb_calcADIncrements(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVAR, NVARtransformed, ndim1, ndim2, Dx, Dflux, Dalpha,&
        Dpp, Dpm, fcb_calcFluxTransformation_sim, DfluxPrel, rcollection)

    use fsystem
    use collection

!<description>
    ! This subroutine assembles the sums of antidiffusive increments
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

    ! Solution used for flux transformation
    real(DP), dimension(ndim1,ndim2), intent(in) :: Dx

    ! Raw antidiffusive flux
    real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux

    ! Edgewise correction factor
    real(DP), dimension(:), intent(in) :: Dalpha

    ! Index pointer for edge data structure
    integer, dimension(:), intent(in) :: IedgeListIdx

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IedgeList

    ! OPTIONAL: callback function to compute variable transformation
    include 'intf_calcFluxTransform_sim.inc'
    optional :: fcb_calcFluxTransformation_sim

    ! OPTIONAL: Antidiffusive flux for prelimiting
    real(DP), dimension(NVAR,NEDGE), intent(in), optional :: DfluxPrel
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Sums of antidiffusive increments
    real(DP), dimension(NVARtransformed,NEQ), intent(out) :: Dpp,Dpm
!</output>
!</subroutine

  end subroutine fcb_calcADIncrements

  end interface
