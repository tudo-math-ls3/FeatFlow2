  ! -*- mode: f90; -*-
  
  ! This file describes the interfaces of the callback routines
  ! of the discrete operators for scalar conservation laws
  
  ! ***************************************************************************
  
  interface

!<subroutine>
    
  subroutine fcb_calcMatrixDiagonal_sim(DdataAtNode, DmatrixCoeffsAtNode,&
      IverticesAtNode, dscale, nnodes, DcoefficientsAtNode, rcollection)

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
  integer, dimension(:,:), intent(in) :: IverticesAtNode

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

  end subroutine fcb_calcMatrixDiagonal_sim

  end interface

  ! ***************************************************************************

  interface

!<subroutine>
    
  subroutine fcb_calcMatrix_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, nedges, DcoefficientsAtEdge, rcollection)

    use fsystem
    use collection

!<description>
    ! This subroutine computes the local matrix entries for the edge
    ! $ij$. Given the conservative variables at the two states $i$ and
    ! $j$ and the multi-component weighting coefficients $C_{ij}$ and
    ! $C_{ji}$, this subroutine computes the matrix coefficients
    ! $k_{ij} and $k_{ji} which are assembled into the global system
    ! matrix. Moreover, the diffusion coefficient $d_{ij}$ is computed.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(2,nedges)
  real(DP), dimension(:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IverticesAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale

  ! Number of edges
  integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  !   DIMENSION(ncoeffs,nedges)
  ! with ncoeffs the number of matrix coefficients at the edge
  real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
!</output>

!</subroutine>

  end subroutine fcb_calcMatrix_sim

  end interface
