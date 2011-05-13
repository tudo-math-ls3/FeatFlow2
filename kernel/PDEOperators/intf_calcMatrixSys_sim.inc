   ! -*- mode: f90; -*-

  interface

!<subroutine>
    
  subroutine fcb_calcMatrixSys_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
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
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(nmatrices,2,nedges)
  ! with nmatrices the number of coefficient matrices
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
  !   DIMENSION(mvar,ncoeffs,nedges)
  ! with ncoeffs the number of matrix coefficients at the edge and
  ! mvar the size of the local matrix entry at the edge.
  ! For local block-diagonal matrices mvar=nvar and for locally
  ! dense matrices mvar=nvar*nvar.
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>

!</subroutine>

  end subroutine fcb_calcMatrixSys_sim

  end interface
