  ! -*- mode: f90; -*-
  
  ! This file describes the interfaces of the callback routines
  ! of the discrete operators for the group finite element method

  ! ***************************************************************************

  interface

  ! ***************************************************************************

!<subroutine>

  subroutine fcb_calcADIncrements(IverticesAtEdgeIdx, IverticesAtEdge,&
        NEDGE, NEQ, NVAR, NVARtransformed, ndim1, ndim2, Dx, Dflux, Dalpha,&
        Dpp, Dpm, fcb_calcFluxTransformation_sim, Dflux0, rcollection)

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
    integer, dimension(:), intent(in) :: IverticesAtEdgeIdx

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! OPTIONAL: callback function to compute variable transformation
    include 'intf_gfsyscallback.inc'
    optional :: fcb_calcFluxTransformation_sim

    ! OPTIONAL: Antidiffusive flux for prelimiting
    real(DP), dimension(NVAR,NEDGE), intent(in), optional :: Dflux0
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

  ! ***************************************************************************

!<subroutine>

  subroutine fcb_calcBounds(IverticesAtEdgeIdx, IverticesAtEdge,&
      NEDGE, NEQ, NVAR, NVARtransformed, ndim1, ndim2, Dx,&
      Dqp, Dqm, fcb_calcDiffTransformation_sim, rcollection)

    use fsystem
    use collection

!<description>
  ! This subroutine assembles the local bounds
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

    ! Solution values
    real(DP), dimension(ndim1,ndim2), intent(in) :: Dx

    ! Index pointer for edge data structure
    integer, dimension(:), intent(in) :: IverticesAtEdgeIdx

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! OPTIONAL: callback function to compute variable transformation
    include 'intf_gfsyscallback.inc'
    optional :: fcb_calcDiffTransformation_sim
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Local bounds
    real(DP), dimension(NVARtransformed,NEQ), intent(out) :: Dqp,Dqm
!</output>
!</subroutine>

  end subroutine fcb_calcBounds

  ! ***************************************************************************

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

  ! ***************************************************************************

!<subroutine>

  subroutine fcb_limitEdgewise(IverticesAtEdge, NEDGE, NEQ, NVAR,&
      NVARtransformed, ndim1, ndim2, Dx, Dflux, Drp, Drm,&
      Dalpha, fcb_calcFluxTransformation_sim, Dflux0, rcollection)

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

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! OPTIONAL: callback function to compute variable transformation
    include 'intf_gfsyscallback.inc'
    optional :: fcb_calcFluxTransformation_sim

    ! OPTIONAL: Antidiffusive flux for constraining
    real(DP), dimension(NVAR,NEDGE), intent(in), optional :: Dflux0
!</intput>

!<inputoutput>
    ! Edgewise correction factors
    real(DP), dimension(:), intent(inout) :: Dalpha

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

  end subroutine fcb_limitEdgewise

  ! ***************************************************************************

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
