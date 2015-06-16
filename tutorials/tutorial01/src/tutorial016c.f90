!##############################################################################
!# Tutorial 016c: Solve simple linear system with CG solver, SSOR prec.
!##############################################################################

module tutorial016c

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput

  use linearalgebra
  use triangulation
  use meshgeneration

  use element
  use cubature
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use bilinearformevaluation

  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop

  use discretebc
  use bcassembly
  use meshregion
  use vectorfilters
  use matrixfilters

  use linearsolver
  use collection

  use ucd

  implicit none
  private

  public :: start_tutorial016c

contains

  ! ***************************************************************************
!<subroutine>

  subroutine fgetBoundaryValuesMR (Icomponents,rdiscretisation,rmeshRegion,&
      cinfoNeeded,Iwhere,Dwhere,Dcoords,Dvalues,rcollection)

!<description>
  ! Auxiliary function that returns values on the boundary.
!</description>

!<input>
  ! Component specifier.
  integer, dimension(:), intent(in) :: Icomponents

  ! Discretisation structure of the underlying equation.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation

  ! Mesh region that is currently being processed.
  type(t_meshRegion), intent(in) :: rmeshRegion

  ! The type of information, the routine should calculate.
  integer, intent(in) :: cinfoNeeded

  ! Information about the DOF to be processed
  integer, dimension(4), intent(in) :: Iwhere

  ! The coordinates of the point which is currently processed.
  real(DP), dimension(:), intent(in) :: Dwhere

  ! The coordinates of the point for which the boundary values are to be calculated.
  real(DP), dimension(:), intent(in) :: Dcoords
!</input>

!<inputoutput>
  ! Optional: A collection structure to provide additional information.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Dvalues(1) receives the calculated value on the boundary.
  real(DP), dimension(:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! Zero boundary values.
    Dvalues(1) = 0.0_DP

  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial016c

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    type(t_matrixBlock) :: rmatrix
    type(t_vectorBlock) :: rrhs, rsolution, rtemp
    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_discreteBC) :: rdiscreteBC
    type(t_meshRegion) :: rmeshRegion
    type(t_ucdExport) :: rexport
    character(LEN=SYS_STRLEN) :: spostdir

    type(t_linsolNode), pointer :: p_rsolverNode, p_rpreconditioner
    integer :: ierror

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 016c")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Create a brick mesh
    ! =================================

    ! The mesh must always be in "standard" format.
    ! First create a 51x51-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (rtriangulation, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 50, 50)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Discretise with Q1.
    !
    ! Create a structure rspatialDiscr
    ! which describes the discretisation.
    ! We generate a 1x1 system with Q1.
    ! =================================

    ! Create a spatial discretisation with Q1
    call spdiscr_initDiscr_simple (rspatialDiscr,EL_Q1_2D,rtriangulation)

    ! Create a block discretisation with 1 block Q1.
    call spdiscr_initBlockDiscr (rblockDiscr,rtriangulation)
    call spdiscr_appendBlockComponent (rblockDiscr,rspatialDiscr)
    call spdiscr_commitBlockDiscr (rblockDiscr)

    ! =================================
    ! Assemble a matrix, a RHS and create
    ! an empty solution vector.
    ! =================================

    ! Use a Gauss 3x3 formula for the discretisation.
    call spdiscr_createDefCubStructure (rspatialDiscr,rcubatureInfo,CUB_GEN_AUTO_G3)

    ! -----------------------------------------------------
    ! Initialise a block matrix and block vectors
    ! corresponding to the block discretisation.
    call lsysbl_createMatrix (rblockDiscr,rmatrix)
    call lsysbl_createVector (rblockDiscr,rrhs)
    call lsysbl_createVector (rblockDiscr,rsolution)
    call lsysbl_createVector (rblockDiscr,rtemp)

    ! Create an empty matrix in block (1,1) in CSR format
    call bilf_createMatrixStructure (rmatrix,1,1,LSYSSC_MATRIX9)
    call lsysbl_allocEmptyMatrix (rmatrix)

    ! -----------------------------------------------------
    ! Discretise a simple Laplace matrix into (1,1)
    call lsysbl_clearMatrix (rmatrix)
    call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD,&
        bma_fcalc_laplace,rcubatureInfo=rcubatureInfo)

    ! -----------------------------------------------------
    ! Create a RHS vector into block (1)
    call lsysbl_clearVector (rrhs)
    call bma_buildVector (rrhs,BMA_CALC_STANDARD,&
        bma_fcalc_rhsBubble,rcubatureInfo=rcubatureInfo)

    ! =================================
    ! Discretise boundary conditions
    ! =================================

    ! Initialise a boundary condition structure
    call bcasm_initDiscreteBC(rdiscreteBC)

    ! Get a mesh region for the complete boundary
    call mshreg_createFromNodalProp(rmeshRegion, rtriangulation, MSHREG_IDX_ALL)

    ! Discretise Dirichlet boundary conditions
    call bcasm_newDirichletBConMR (rblockDiscr, 1, rdiscreteBC, &
        rmeshRegion, fgetBoundaryValuesMR)

    call mshreg_done(rmeshregion)

    ! =================================
    ! Impose boundary conditions
    ! =================================

    ! Impose the BC into the matrix and the vectors.
    call matfil_discreteBC (rmatrix,rdiscreteBC)
    call vecfil_discreteBCrhs (rrhs,rdiscreteBC)
    call vecfil_discreteBCsol (rsolution,rdiscreteBC)

    ! =================================
    ! Solve the system with Gauss elimination
    ! =================================

    call output_line ("Solving linear system...")

    ! ----------------
    ! Solver preparation
    ! ----------------

    ! Create a CG solver with SSOR preconditioner
    call linsol_initSSOR (p_rpreconditioner)
    call linsol_initCG (p_rsolverNode,p_rpreconditioner)

    ! Attach the system matrix
    call linsol_setMatrix (p_rsolverNode, rmatrix)

    ! Symbolic factorisation
    call linsol_initStructure (p_rsolverNode, ierror)

    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line ("Error during symbolic factorisation.")
      call sys_halt()
    end if

    ! Numeric factorisation
    call linsol_initData (p_rsolverNode, ierror)

    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line ("Error during numeric factorisation.")
      call sys_halt()
    end if

    ! ----------------
    ! Set stopping criteria, etc.
    ! ----------------

    p_rsolverNode%depsRel        = 1E-8_DP        ! ||final res|| <= 1E-8 * ||initial res||
    p_rsolverNode%depsAbs        = 1E-5_DP        ! ||final res|| <= 1E-5
    p_rsolverNode%istoppingCriterion = LINSOL_STOP_STANDARD  ! Check for abs. + rel. res.
    p_rsolverNode%nmaxIterations = 1000           ! <= 1000 iterations
    p_rsolverNode%iresNorm       = LINALG_NORML2  ! <= check l2-norm of the vector
    p_rsolverNode%ioutputlevel   = 2              ! Print some output about the iteration
    p_rsolverNode%niteResOutput  = 10             ! Print every 10th residual

    ! ----------------
    ! Solve the system
    ! ----------------

    ! Clear the solution
    call lsysbl_clearVector (rsolution)

    ! Solve
    call linsol_solveAdaptively (p_rsolverNode,rsolution,rrhs,rtemp)

    ! ----------------
    ! Cleanup
    ! ----------------

    ! Numeric data
    call linsol_doneData (p_rsolverNode)

    ! Symbolic data
    call linsol_doneStructure (p_rsolverNode)

    ! Remaining solver data
    call linsol_releaseSolver (p_rsolverNode)

    ! =================================
    ! Create output files
    ! =================================

    call output_line ("Writing postprocessing files...")

    ! Open / write / close; write the solution to a VTK file.
    if (sys_getenv_string("POSTDIR",spostdir)) then
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       trim(spostdir)//"/tutorial016c.vtk")
    else
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       "post/tutorial016c.vtk")
    end if

    call ucd_addVectorByVertex (rexport, "solution", &
        UCD_VAR_STANDARD, rsolution%RvectorBlock(1))
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Cleanup
    ! =================================

    ! Release the BC
    call bcasm_releaseDiscreteBC(rdiscreteBC)

    ! Release the matrix/vectors
    call lsysbl_releaseMatrix (rmatrix)
    call lsysbl_releaseVector (rtemp)
    call lsysbl_releaseVector (rrhs)
    call lsysbl_releaseVector (rsolution)

    ! Release the discretisation
    call spdiscr_releaseBlockDiscr (rblockDiscr)
    call spdiscr_releaseDiscr (rspatialDiscr)

    ! Cubature done.
    call spdiscr_releaseCubStructure (rcubatureInfo)

    ! Release the triangulation
    call tria_done (rtriangulation)

  end subroutine

end module
