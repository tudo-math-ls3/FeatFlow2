   ! -*- mode: f90; -*-

  interface

!<subroutine>
    
  subroutine fcb_calcMatrixDiagSys_sim(DdataAtNode, DmatrixCoeffsAtNode,&
      IdofsAtNode, dscale, nnodes, DcoefficientsAtNode, rcollection)

    use fsystem
    use collection

!<description>
    ! This subroutine computes the local matrix entries for the node
    ! $i$. Given the conservative variables at states $i$ and the
    ! multi-component weighting coefficients $C_{ii}$, this subroutine
    ! computes the matrix coefficients $k_{ii} which are assembled
    ! into the global system matrix.
!</description>

!<input>
  ! Nodal solution values for all nodes under consideration
  !   DIMENSION(nvar,nnodes)
  ! with nvar the number of variables at the node
  real(DP), dimension(:,:), intent(in) :: DdataAtNode

  ! Entries of the coefficient matrices for all nodes under consideration
  !   DIMENSION(nmatrices,nnodes)
  ! with nmatrices the number of coefficient matrices
  real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode

  ! Numbers of vertices and matrix entries for all nodes under consideration
  !   DIMENSION(2,nnodes)
  integer, dimension(:,:), intent(in) :: IdofsAtNode

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
  !   DIMENSION(mvar,ncoeffs,nnodes)
  ! with ncoeffs the number of matrix coefficients at the node and
  ! mvar the size of the local matrix entry at the node.
  ! For local block-diagonal matrices mvar=nvar and for locally
  ! dense matrices mvar=nvar*nvar.
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtNode
!</output>

!</subroutine>

  end subroutine fcb_calcMatrixDiagSys_sim

  end interface
