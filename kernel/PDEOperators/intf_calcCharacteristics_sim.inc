   ! -*- mode: f90; -*-

  interface

!<subroutine>

  subroutine fcb_calcCharacteristics_sim(Dweight, DdataAtEdge,&
      nedges, DcharVariablesAtEdge, DeigenvaluesAtEdge,&
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

  ! Number of edges
  integer, intent(in) :: nedges
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! OPTIONAL: Characteristic variables for all edges under consideration
  !   DIMENSION(nvar,nedges)
  ! with nvar the number of variables at each edge
  real(DP), dimension(:,:), intent(out), optional :: DcharVariablesAtEdge

  ! OPTIONAL: Eigenvalues for all edges under consideration
  !   DIMENSION(nvar,nedges)
  ! with nvar the number of variables at each edge
  real(DP), dimension(:,:), intent(out), optional :: DeigenvaluesAtEdge

  ! OPTIONAL: Matrices of left eigenvectors for all edges under consideration
  !   DIMENSION(nvar*nvar,nedges)
  ! with nvar the number of variables at each edge
  real(DP), dimension(:,:), intent(out), optional :: DleftEigenvectorsAtEdge

  ! OPTIONAL: Matrices of right eigenvectors for all edges under consideration
  !   DIMENSION(nvar*nvar,nedges)
  ! with nvar the number of variables at each edge
  real(DP), dimension(:,:), intent(out), optional :: DrightEigenvectorsAtEdge
!</output>

!</subroutine>

  end subroutine fcb_calcCharacteristics_sim

  end interface