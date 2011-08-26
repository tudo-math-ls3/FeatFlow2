  ! -*- mode: f90; -*-

  interface

!<subroutine>

  subroutine fcb_limitEdgewise(IverticesAtEdgeIdx, IverticesAtEdge,&
      NEDGE, NEQ, NVAR, NVARtransformed, ndim1, ndim2, Dx, Dflux, Dalpha,&
      Drp, Drm, fcb_calcFluxTransformation_sim, DfluxConstr, rcollection)

    use fsystem
    use collection

!<description>
    ! This subroutine computes the edgewise correction factors
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

    ! Nodal correction factors
    real(DP), dimension(NVARtransformed,NEQ), intent(in) :: Drp,Drm

    ! Index pointer for edge data structure
    integer, dimension(:), intent(in) :: IverticesAtEdgeIdx

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! OPTIONAL: callback function to compute variable transformation
    include 'intf_calcFluxTransformation_sim.inc'
    optional :: fcb_calcFluxTransformation_sim

    ! OPTIONAL: Antidiffusive flux for constraining
    real(DP), dimension(NVAR,NEDGE), intent(in), optional :: DfluxConstr
!</intput>

!<inputoutput>
    ! Edgewise correction factors
    real(DP), dimension(:), intent(inout) :: Dalpha

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

  end subroutine fcb_limitEdgewise

  end interface
