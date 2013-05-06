!##############################################################################
!# Tutorial 015c: Impose Dirichlet boundary conditions on bench1
!##############################################################################

module tutorial015c

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  
  use boundary
  use triangulation
  use meshgeneration
  use mprimitives
  
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
  use vectorfilters
  use matrixfilters
  
  use collection
  
  use matrixio
  use vectorio
  use ucd

  implicit none
  private
  
  public :: start_tutorial015c

contains

  ! ***************************************************************************
  !<subroutine>

    subroutine fgetBoundaryValues (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                   cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)

    use fsystem
    use boundary
    use collection
    use spatialdiscretisation
    use discretebc

  !<description>
    ! Auxiliary function that returns values on the boundary.
  !</description>

  !<input>
    ! Component specifier.
    !   Icomponents(1) defines the number of the solution component.
    integer, dimension(:), intent(in)                           :: Icomponents

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation

    ! Boundary region that is currently being processed.
    type(t_boundaryRegion), intent(in)                          :: rboundaryRegion

    ! The element number on the boundary which is currently being processed
    integer, intent(in)                                         :: ielement

    ! The type of information, the routine should calculate. One of the
    ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
    ! to return one or multiple information value in the result array.
    integer, intent(in)                                         :: cinfoNeeded

    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC :
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    integer, intent(in) :: iwhere

    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC :
    !   dwhere = parameter value of the point where the value should be computed
    real(DP), intent(in) :: dwhere
  !</input>

  !<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
  !</inputoutput>

  !<output>
    ! This array receives the calculated information. If the caller
    ! only needs one value, the computed quantity is put into Dvalues(1).
    ! The function may return SYS_INFINITY as a value, indicating
    ! 'do-nothing boundary conditions'.
    real(DP), dimension(:), intent(out) :: Dvalues
  !</output>

  !</subroutine>

    ! Local variables
    real(DP) :: dx,dy
    
    ! Get X/Y coordinates of the point
    call boundary_getCoords(rdiscretisation%p_rboundary, &
        rboundaryRegion%iboundCompIdx, dwhere, dx, dy)
    
    ! Boundary component?
    select case (Icomponents(1))
    
    ! 1st boundary component. Return a parabolic profile.
    ! This implies f=0 on the top and bottom edge.
    case (1)
      Dvalues(1) = mprim_getParabolicProfile (dy,0.41_DP,1.0_DP)
    
    ! 2nd boundary component, the circle. Return constant =1.
    case (2)
      Dvalues(1) = 1.0_DP

    end select

  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial015c

    ! Declare some variables.
    type(t_boundary) :: rboundary
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    type(t_matrixBlock) :: rmatrix
    type(t_vectorBlock) :: rrhs, rsolution
    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_discreteBC) :: rdiscreteBC
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_ucdExport) :: rexport

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 015c")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Read the underlying domain
    ! and the mesh
    ! =================================

    call boundary_read_prm(rboundary, "pre/bench1.prm")
    
    ! The mesh must always be in "standard" format to work with it.
    ! First read, then convert to standard.
    call tria_readTriFile2D (rtriangulation, "pre/bench1.tri", rboundary)
    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)

    ! =================================
    ! Discretise with Q1.
    !
    ! Create a structure rspatialDiscr
    ! which describes the discretisation.
    ! We generate a 1x1 system with Q1.
    ! =================================

    ! Create a spatial discretisation with Q1
    call spdiscr_initDiscr_simple (rspatialDiscr,EL_Q1_2D,rtriangulation,rboundary)
    
    ! Create a block discretisation with 1 block Q1.
    call spdiscr_initBlockDiscr (rblockDiscr,rtriangulation,rboundary)
    call spdiscr_appendBlockComponent (rblockDiscr,rspatialDiscr)
    call spdiscr_commitBlockDiscr (rblockDiscr)

    ! =================================
    ! Assemble a matrix, a RHS and create 
    ! an empty solution vector.
    ! =================================

    ! Use a Gauss 3x3 formula for the discretisation.
    call spdiscr_createDefCubStructure (rspatialDiscr,rcubatureInfo,CUB_GEN_AUTO_G3)

    ! -----------------------------------------------------
    ! Initialise a block martix and block vectors
    ! corresponding to the block discretisation
    call lsysbl_createMatrix (rblockDiscr,rmatrix)
    call lsysbl_createVector (rblockDiscr,rrhs)
    call lsysbl_createVector (rblockDiscr,rsolution)

    ! Create an empty matrix in block (1,1) in CSR format, 
    ! initialised with the FEM matrix structure.
    call bilf_createMatrixStructure (rspatialDiscr,LSYSSC_MATRIX9,rmatrix%RmatrixBlock(1,1))
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
    
    ! -----------------------------------------------------
    ! Clear the solution
    call lsysbl_clearVector (rsolution)
    
    ! -----------------------------------------------------
    ! Cubature done.
    call spdiscr_releaseCubStructure (rcubatureInfo)

    ! =================================
    ! Discretise boundary conditions
    ! =================================
    
    ! Initialise a boundary condition structure
    call bcasm_initDiscreteBC(rdiscreteBC)
    
    ! Dirichlet-BC for the bottom edge, 1st boundary component. The edge has segment-ID=1.
    call boundary_createRegion (rboundary, 1, 1, rboundaryRegion)
    call bcasm_newDirichletBConRealBd (rblockDiscr, &
        1, rboundaryRegion, rdiscreteBC, fgetBoundaryValues)

    ! Dirichlet-BC for the top edge, 1st boundary component. The edge has segment-ID=3.
    call boundary_createRegion (rboundary, 1, 3, rboundaryRegion)
    call bcasm_newDirichletBConRealBd (rblockDiscr, &
        1, rboundaryRegion, rdiscreteBC, fgetBoundaryValues)

    ! Dirichlet-BC for the left edge, 1st boundary component. The edge has segment-ID=4.
    call boundary_createRegion (rboundary, 1, 4, rboundaryRegion)
    call bcasm_newDirichletBConRealBd (rblockDiscr, &
        1, rboundaryRegion, rdiscreteBC, fgetBoundaryValues)

    ! Dirichlet-BC for the complete circle. This is boundary component 2.
    call boundary_createRegion (rboundary, 2, 0, rboundaryRegion)
    call bcasm_newDirichletBConRealBd (rblockDiscr, &
        1, rboundaryRegion, rdiscreteBC, fgetBoundaryValues)

    ! =================================
    ! Impose boundary conditions
    ! =================================

    ! Impose the BC into the matrix and the vectors.
    call matfil_discreteBC (rmatrix,rdiscreteBC)
    call vecfil_discreteBCrhs (rrhs,rdiscreteBC)
    call vecfil_discreteBCsol (rsolution,rdiscreteBC)

    ! =================================
    ! Create output files
    ! =================================

    call output_line ("Writing postprocessing files...")

    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeBlockMatrixHR (rmatrix, "matrix", .true., 0, &
        "post/tutorial015c_matrix.txt", "(E11.2)")

    ! Write the vector to a text file.
    call vecio_writeBlockVectorHR (rrhs, "rhs", .true., 0, &
        "post/tutorial015c_rhs.txt", "(E11.2)")

    ! Write the vector to a text file.
    call vecio_writeBlockVectorHR (rsolution, "vector", .true., 0, &
        "post/tutorial015c_sol.txt", "(E11.2)")

    ! Open / write / close; write the solution to a VTK file.
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       "post/tutorial015c.vtk")
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
    call lsysbl_releaseVector (rrhs)
    call lsysbl_releaseVector (rsolution)
    
    ! Release the discretisation
    call spdiscr_releaseBlockDiscr (rblockDiscr)
    call spdiscr_releaseDiscr (rspatialDiscr)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
    ! Release the boundary definition
    call boundary_release(rboundary)
    
  end subroutine

end module
