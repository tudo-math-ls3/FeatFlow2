!##############################################################################
!# Tutorial 015a: Impose Dirichlet boundary conditions, use mesh regions
!##############################################################################

module tutorial015a

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  
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
  
  use collection
  
  use matrixio
  use vectorio
  use ucd

  implicit none
  private
  
  public :: start_tutorial015a

contains

  ! ***************************************************************************
!<subroutine>

  subroutine fgetBoundaryValuesMR (Icomponents,rdiscretisation,rmeshRegion,&
                                    cinfoNeeded,Iwhere,Dwhere,Dcoords,Dvalues,&
                                    rcollection)

!<description>
  ! Auxiliary function that returns values on the boundary.
!</description>

!<input>
  ! Component specifier.
  ! For Dirichlet boundary:
  !   Icomponents(1) defines the number of the solution component. Here =1.
  integer, dimension(:), intent(in) :: Icomponents

  ! Discretisation structure of the underlying equation.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation

  ! Mesh region that is currently being processed.
  type(t_meshRegion), intent(in) :: rmeshRegion

  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  integer, intent(in) :: cinfoNeeded

  ! An array holding information about what type of DOF is currently processed.
  ! The information is build up as follows:
  ! Iwhere(1) = vertice number of the DOF, if the DOF is vertice-based, otherwise 0
  ! Iwhere(2) = edge number of the DOF, if the DOF is edge-based, otherwise 0
  ! Iwhere(3) = face number of the DOF, if the DOF is face-based, otherwise 0
  ! Iwhere(4) = currently processed element number.
  ! If Iwhere(1) = Iwhere(2) = Iwhere(3) = 0, then the DOF is element based.
  integer, dimension(4), intent(in) :: Iwhere

  ! The coordinates of the point which is currently processed, given in
  ! reference coordinates of the currently processed cell type (edge,face,element).
  ! If the DOF is vertice-based, then Dwhere is undefined.
  ! If the DOF is edge-based or element-based in 1D, then Dwhere has dimension 1.
  ! If the DOF is face-based or element-based in 2D, then Dwhere has dimension 2.
  ! IF the DOF is element-based in 3D, then Dwhere has dimension 3.
  real(DP), dimension(:), intent(in) :: Dwhere

  ! The coordinates of the point for which the boundary values are to be
  ! calculated.
  real(DP), dimension(:), intent(in) :: Dcoords
!</input>

!<inputoutput>
  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine. May point to NULL() if not defined.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1).
  ! If multiple values are needed, they are collected here (e.g. for
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  !
  ! The function may return SYS_INFINITY as a value. This indicates the
  ! framework to ignore the node and treat it as 'natural boundary condition'
  ! node.
  real(DP), dimension(:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    real(DP) :: dx,dy
    
    ! Get the coordinates
    dx = Dcoords(1)
    dy = Dcoords(2)
    
    ! We return the Dirichlet values: u(x,y) = y.
    ! This is =0 on the lower edge and =1 on the upper edge.
    Dvalues(1) = dy

  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial015a

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    type(t_matrixBlock) :: rmatrix
    type(t_vectorBlock) :: rrhs, rsolution
    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_discreteBC) :: rdiscreteBC
    type(t_meshRegion) :: rmeshRegion
    type(t_ucdExport) :: rexport

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 015a")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Create a brick mesh
    ! =================================

    ! The mesh must always be in "standard" format. 
    ! First create a 5x5-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (rtriangulation, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 4, 4)
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
    
    ! Get a mesh region for the bottommost edge and
    ! create Dirichlet-BC there. The values are returned
    ! by the above callback function.
    call mshreg_createFromExpression(rmeshRegion, rtriangulation, &
        MSHREG_IDX_ALL, .true., "Y = 0")
    
    call bcasm_newDirichletBConMR (rblockDiscr, 1, rdiscreteBC, &
        rmeshRegion, fgetBoundaryValuesMR)
    
    call mshreg_done(rmeshregion)
        
    ! Another one for the topmost edge
    call mshreg_createFromExpression(rmeshRegion, rtriangulation, &
        MSHREG_IDX_ALL, .true., "Y = 1")
    
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
    ! Create output files
    ! =================================

    call output_line ("Writing postprocessing files...")

    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeBlockMatrixHR (rmatrix, "matrix", .true., 0, &
        "post/tutorial015a_matrix.txt", "(E11.2)")

    ! Write the vector to a text file.
    call vecio_writeBlockVectorHR (rrhs, "rhs", .true., 0, &
        "post/tutorial015a_rhs.txt", "(E11.2)")

    ! Write the vector to a text file.
    call vecio_writeBlockVectorHR (rsolution, "vector", .true., 0, &
        "post/tutorial015a_sol.txt", "(E11.2)")

    ! Open / write / close; write the solution to a VTK file.
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       "post/tutorial015a.vtk")
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
    
  end subroutine

end module
