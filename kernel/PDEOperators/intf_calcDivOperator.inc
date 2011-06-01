   ! -*- mode: f90; -*-

  interface

!<subroutine>

  subroutine fcb_calcDivOperator(rafcstab, rx, rmatrix, dscale, bclear,&
      fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim, rcollection)

    use afcstabilisation
    use collection
    use fsystem
    use linearsystemblock

!<description>
    ! This subroutine assembles the divergence flux operator
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: callback function to compute local fluxes
    include 'intf_calcMatrixDiagSys_sim.inc'
    include 'intf_calcMatrixSys_sim.inc'
    optional :: fcb_calcMatrixDiagSys_sim
    optional :: fcb_calcMatrixSys_sim
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! divergence operator
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

  end subroutine fcb_calcDivOperator

  end interface
