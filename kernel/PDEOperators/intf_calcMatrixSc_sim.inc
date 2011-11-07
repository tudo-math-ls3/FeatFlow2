  ! -*- mode: f90; -*-

  interface

!<subroutine>
    
  subroutine fcb_calcMatrixSc_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DmatrixAtEdge, rcollection)

    use fsystem
    use collection

!<description>
    ! Given the solution data DdataAtEdge and auxiliary coefficients
    ! DcoeffsAtEdge this subroutine computes the local matrix entries
    ! DmatrixAtEdge for the edge $ij$.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(2,nedges)
  real(DP), dimension(:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of nodes and matrix entries for all edges under consideration
  !   DIMENSION(6,nedges)
  integer, dimension(:,:), intent(in) :: IedgeList

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
  real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>

!</subroutine>

  end subroutine fcb_calcMatrixSc_sim

  end interface
