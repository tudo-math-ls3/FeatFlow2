   ! -*- mode: f90; -*-

  interface

!<subroutine>

  subroutine fcb_calcOperatorEdgeSys(rgroupFEMSet, rx, rmatrix, dscale, bclear,&
      fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim, rcollection, rafcstab)

    use afcstabbase
    use collection
    use fsystem
    use groupfembase
    use linearsystemblock

!<description>
    ! This subroutine assembles the discrete operator for systems of
    ! equations using the group finite element formulation by Fletcher.
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
    include 'intf_calcMatrixDiagSys_sim.inc'
    include 'intf_calcMatrixSys_sim.inc'
    optional :: fcb_calcMatrixDiagSys_sim
    optional :: fcb_calcMatrixSys_sim
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

  end subroutine fcb_calcOperatorEdgeSys

  end interface
