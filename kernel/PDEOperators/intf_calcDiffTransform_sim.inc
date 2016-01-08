  interface

!<subroutine>

  subroutine fcb_calcDiffTransformation_sim(DdataAtEdge,&
      nedges, DtransformedDataAtEdge, rcollection)

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

  ! Number of edges
  integer, intent(in) :: nedges
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