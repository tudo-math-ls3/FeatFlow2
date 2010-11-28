  ! -*- mode: f90; -*-
  
  ! This file describes the interfaces of the callback routines
  ! of the discrete operators for hyperbolic systems

  ! ***************************************************************************
  
  interface

!<subroutine>
    
  subroutine fcb_calcMatrixDiagonal_sim(DdataAtNode, DmatrixCoeffsAtNode,&
      IverticesAtNode, dscale, DcoefficientsAtNode, rcollection)

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
  integer, dimension(:,:), intent(in) :: IverticesAtNode

  ! Scaling parameter
  real(DP), intent(in) :: dscale
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

  end subroutine fcb_calcMatrixDiagonal_sim

  end interface

  ! ***************************************************************************

  interface

!<subroutine>
    
  subroutine fcb_calcMatrix_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DcoefficientsAtEdge, rcollection)

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

  end subroutine fcb_calcMatrix_sim

  end interface

  ! ***************************************************************************

  interface

!<subroutine>
    
  subroutine fcb_calcFlux_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

    use fsystem
    use collection

!<description>
    ! This subroutine computes the internodal fluxes. Given the
    ! conservative variables at the two states $i$ and $j$ and the
    ! multi-component weighting coefficients $C_{ij}$ and $C_{ji}$,
    ! this subroutine computes the internodal flux from node $j$ into
    ! node $i$: $F_{ij}$ and the internodal flux from node $i$ into
    ! node $j$: $F_{ji}$. In some cases it is useful to have
    ! skew-symmetry, i.e. $F_{ji}=-F_{ij}$ but this is not mandatory,
    ! and hence, this routine computes both.
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
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

  end subroutine fcb_calcFlux_sim

  end interface

  ! ***************************************************************************

  interface

!<subroutine>
    
  subroutine fcb_calcFlux2_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IverticesAtEdge, dscale, DfluxesAtEdge, rcollection)

    use fsystem
    use collection

!<description>
    ! This subroutine computes the internodal fluxes. Given the
    ! conservative variables at the two states $i$ and $j$ and the
    ! multi-component weighting coefficients $C_{ij}$ and $C_{ji}$,
    ! this subroutine computes the internodal fluxes from node
    ! $j$ into node $i$ $F_{ij}$ and the internodal flux from node
    ! $i$ into node $j$ $F_{ji}=-F_{ij}$.
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
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

  end subroutine fcb_calcFlux2_sim

  end interface

  ! ***************************************************************************

  interface

!<subroutine>

  subroutine fcb_calcFluxTransformation_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

    use fsystem
    use collection

!<description>
    ! This subroutine computes the transformed fluxes.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Transformed internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>

!</subroutine>

  end subroutine fcb_calcFluxTransformation_sim

  end interface

  ! ***************************************************************************

  interface

!<subroutine>

  subroutine fcb_calcDiffTransformation_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

    use fsystem
    use collection

!<description>
    ! This subroutine computes the transformed solution difference.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Difference of transformed solution values for all edges under consideration
  !   DIMENSION(nvar,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>

!</subroutine>

  end subroutine fcb_calcDiffTransformation_sim

  end interface

  ! ***************************************************************************

  interface

!<subroutine>

  subroutine fcb_calcCharacteristics_sim(Dweight, DdataAtEdge,&
      DcharVariablesAtEdge, DeigenvaluesAtEdge,&
      DrightEigenvectorsAtEdge, DleftEigenvectorsAtEdge, rcollection)

    use fsystem
    use collection

!<description>
    ! This subroutine computes the characteristic transformation
    ! variables.  Given the conservative variables at the two states
    ! $i$ and $j$ and the multi-component weighting coefficient, this
    ! routine computes:
    !
    ! 1.) The characteristic solution difference
    !     $W_{ij}=R_{ij}(U_j-U_i)$
    !
    ! 2.) The diagonal matrix of real eigenvalues $\Lambda_{ij}$
    !
    ! 3.) The matrix of right eigenvectors $R_{ij}$
    !
    ! 4.) The matrix of left eigenvectors $L_{ij}$
    !
    ! For hyperbolic systems, the following characteristic
    ! decomposition holds
    !
    ! $$A_{ij} = L_{ij} \Lambda_{ij} R_{ij}$$
!</description>

!<input>
  ! Weighting coefficient for wave-decomposition
  real(DP), dimension(:), intent(in)  :: Dweight

  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! OPTIONAL: Characteristic variables for all edges under consideration
  !   DIMENSION(nvar,nedge)
  ! with nvar the number of variables at each edge
  real(DP), dimension(:,:), intent(out), optional :: DcharVariablesAtEdge

  ! OPTIONAL: Eigenvalues for all edges under consideration
  !   DIMENSION(nvar,nedge)
  ! with nvar the number of variables at each edge
  real(DP), dimension(:,:), intent(out), optional :: DeigenvaluesAtEdge

  ! OPTIONAL: Matrices of left eigenvectors for all edges under consideration
  !   DIMENSION(nvar*nvar,nedge)
  ! with nvar the number of variables at each edge
  real(DP), dimension(:,:), intent(out), optional :: DleftEigenvectorsAtEdge

  ! OPTIONAL: Matrices of right eigenvectors for all edges under consideration
  !   DIMENSION(nvar*nvar,nedge)
  ! with nvar the number of variables at each edge
  real(DP), dimension(:,:), intent(out), optional :: DrightEigenvectorsAtEdge
!</output>

!</subroutine>

  end subroutine fcb_calcCharacteristics_sim

  end interface
