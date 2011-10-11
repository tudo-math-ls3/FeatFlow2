   ! -*- mode: f90; -*-

  interface

!<subroutine>

  subroutine fcb_calcOperatorNodeSc(rgroupFEMSet, rx, rmatrix, dscale, bclear,&
      fcb_calcMatrixDiagSc_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine assembles the discrete operator for scalar
    ! problems using the group finite element formulation by Fletcher.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the matrix entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: callback function to compute local matrices
    include 'intf_calcMatrixDiagSc_sim.inc'
    optional :: fcb_calcMatrixDiagSc_sim
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

  end subroutine fcb_calcOperatorNodeSc

  end interface
