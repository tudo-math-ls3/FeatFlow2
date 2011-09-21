  ! -*- mode: f90; -*-

  interface

!<subroutine>

  subroutine fcb_calcBounds(IedgeListIdx, IedgeList,&
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
    integer, dimension(:), intent(in) :: IedgeListIdx

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IedgeList

    ! OPTIONAL: callback function to compute variable transformation
    include 'intf_calcDiffTransformation_sim.inc'
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

  end interface
