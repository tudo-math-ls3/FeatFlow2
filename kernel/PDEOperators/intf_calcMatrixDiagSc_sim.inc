  ! -*- mode: f90; -*-

  interface

!<subroutine>
    
  subroutine fcb_calcMatrixDiagSc_sim(DdataAtNode, DmatrixCoeffsAtNode,&
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
  !   DIMENSION(nnodes)
  real(DP), dimension(:), intent(in) :: DdataAtNode

  ! Entries of the coefficient matrices for all nodes under consideration
  !   DIMENSION(ndim,nnodes)
  ! with ndim the number of spatial dimensions
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
  !   DIMENSION(ncoeffs,nnodes)
  ! with ncoeffs the number of matrix coefficients at the node
  real(DP), dimension(:,:), intent(out) :: DcoefficientsAtNode
!</output>

!</subroutine>

  end subroutine fcb_calcMatrixDiagSc_sim

  end interface
