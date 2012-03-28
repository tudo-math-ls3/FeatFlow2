  interface

!<subroutine>

  subroutine fcb_calcVectorNodeSc(rgroupFEMSet, rx, ry, dscale, bclear,&
      fcb_calcVectorSc_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine assembles the discrete vector for scalar
    ! problems using the group finite element formulation by Fletcher.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the vector entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: callback function to compute local vectors
    include 'intf_calcVectorSc_sim.inc'
    optional :: fcb_calcVectorSc_sim
!</input>

!<inputoutput>
    ! Destination vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

  end subroutine fcb_calcVectorNodeSc

  end interface
