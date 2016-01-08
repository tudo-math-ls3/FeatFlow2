  interface

!<subroutine>

  subroutine fcb_calcNodalTransformation_sim(DdataAtNode,&
      nnodes, DtransformedDataAtNode, rcollection)

    use fsystem
    use collection

!<description>
    ! This subroutine computes the transformed solution.
!</description>

!<input>
  ! Nodal solution values for all nodes under consideration
  !   DIMENSION(nvar,nnodes)
  ! with nvar the number of variables at each node
  real(DP), dimension(:,:), intent(in) :: DdataAtNode

  ! Number of nodes
  integer, intent(in) :: nnodes
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Transformed solution values for all nodes under consideration
  !   DIMENSION(nvar,nnodes)
  ! with nvar the number of variables at each node
  real(DP), dimension(:,:), intent(out) :: DtransformedDataAtNode
!</output>

!</subroutine>

  end subroutine fcb_calcNodalTransformation_sim

  end interface