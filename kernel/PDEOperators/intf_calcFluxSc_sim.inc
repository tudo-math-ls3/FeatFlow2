   ! -*- mode: f90; -*-

  interface

!<subroutine>
    
  subroutine fcb_calcFluxSc_sim(DdataAtEdge, DcoeffsAtEdge,&
      IedgeList, dscale, nedges, DfluxesAtEdge, Dcoefficients,&
      rcollection)

    use fsystem
    use collection

!<description>
    ! Given the solution data DdataAtEdge and auxiliary coefficients
    ! DcoeffsAtEdge this subroutine computes the internodal fluxes
    ! from node $j$ into node $i$: $F_{ij}$ and the internodal flux
    ! from node $i$ into node $j$: $F_{ji}$. In some cases it is
    ! useful to have skew-symmetry, i.e. $F_{ji}=-F_{ij}$ but this is
    ! not mandatory, and hence, this routine computes both.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(nmatrices,2,nedges)
  ! with nmatrices the number of coefficient matrices
  real(DP), dimension(:,:,:), intent(in) ::  DcoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
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
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:), intent(out) :: DfluxesAtEdge

  ! OPTIONAL: Coefficients of the matrix for all edges under consideration
  !   DIMENSION(ncoeffs,nedges)
  ! with ncoeffs the number of matrix coefficients at the edge
  real(DP), dimension(:,:), intent(out), optional :: Dcoefficients
!</output>

!</subroutine>

  end subroutine fcb_calcFluxSc_sim

  end interface
