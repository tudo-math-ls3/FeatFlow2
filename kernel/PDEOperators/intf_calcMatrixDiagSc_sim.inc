  interface

!<subroutine>

  subroutine fcb_calcMatrixDiagSc_sim(DdataAtNode, DcoeffsAtNode,&
      InodeList, dscale, nnodes, DmatrixAtNode, rcollection)

    use fsystem
    use collection

!<description>
    ! Given the solution data DdataAtNode and auxiliary coefficients
    ! DcoeffsAtNode this subroutine computes the local matrix entries
    ! DmatrixAtNode for the node $i$.
!</description>

!<input>
  ! Nodal solution values for all nodes under consideration
  !   DIMENSION(nnodes)
  real(DP), dimension(:), intent(in) :: DdataAtNode

  ! Entries of the coefficient matrices for all nodes under consideration
  !   DIMENSION(ndim,nnodes)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode

  ! Numbers of nodes and matrix entries for all nodes under consideration
  !   DIMENSION(2,nnodes)
  integer, dimension(:,:), intent(in) :: InodeList

  ! Scaling parameter
  real(DP), intent(in) :: dscale

  ! Number of nodes
  integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all nodes under consideration
  !   DIMENSION(ncoeffs,nnodes)
  ! with ncoeffs the number of matrix coefficients at the node
  real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>

!</subroutine>

  end subroutine fcb_calcMatrixDiagSc_sim

  end interface
