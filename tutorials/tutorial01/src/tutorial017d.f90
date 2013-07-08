!##############################################################################
!# Tutorial 017d: Solve simple linear system with MG, prol/rest matrices
!##############################################################################

module tutorial017d

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
  
  use filtersupport
  use multilevelprojection
  use multileveloperators
  use linearsolver
  use coarsegridcorrection
  use collection
  
  use ucd

  implicit none
  private
  
  public :: start_tutorial017d

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

  subroutine start_tutorial017d

    ! Declare some variables.
    integer, parameter :: NLMAX = 5
    
    type(t_triangulation), dimension(:), pointer :: p_Rtriangulations
    type(t_spatialDiscretisation), dimension(:), pointer :: p_RspatialDiscr
    type(t_blockDiscretisation), dimension(:), pointer :: p_RblockDiscr
    type(t_discreteBC), dimension(:), pointer :: p_RdiscreteBC
    type(t_filterChain), dimension(:,:), pointer :: p_RfilterChain
    type(t_scalarCubatureInfo), dimension(:), pointer :: p_RcubatureInfo
    integer, dimension(:), pointer :: Isize
    
    type(t_matrixBlock), dimension(:), pointer :: p_Rmatrices
    type(t_vectorBlock) :: rrhs, rsolution, rtemp
    
    type(t_matrixScalar), dimension(:), pointer :: p_RprolMatrices
    type(t_matrixScalar), dimension(:), pointer :: p_RrestMatrices
    type(t_interlevelProjectionBlock), dimension(:), pointer :: p_Rprojection
    
    type(t_meshRegion) :: rmeshRegion
    type(t_ucdExport) :: rexport
    
    type(t_linsolNode), pointer :: p_rsolverNode, p_rsmoother, p_rcoarsegridsolver
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    type(t_linsolMatrixSet) :: rmatrixSet
    integer :: ierror, ilevel

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 017d")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Allocate level structures
    ! =================================
    allocate (p_Rtriangulations(NLMAX))
    allocate (p_RspatialDiscr(NLMAX))
    allocate (p_RblockDiscr(NLMAX))
    allocate (p_RcubatureInfo(NLMAX))
    allocate (p_RdiscreteBC(NLMAX))
    allocate (p_Rmatrices(NLMAX))
    allocate (p_RfilterChain(1,NLMAX))
    allocate (Isize(NLMAX))
    
    allocate (p_RprolMatrices(NLMAX))
    allocate (p_RrestMatrices(NLMAX))
    allocate (p_Rprojection(NLMAX))
    
    ! =================================
    ! Create a brick mesh
    ! =================================

    ! The mesh must always be in "standard" format. 
    ! First create a 5x5-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (p_Rtriangulations(1), 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 4, 4)
    call tria_initStandardMeshFromRaw (p_Rtriangulations(1))
    
    ! =================================
    ! Create a mesh hierarchy
    ! =================================
    
    ! Refine until level NLMAX.
    do ilevel = 2,NLMAX
      call tria_refine2LevelOrdering(p_Rtriangulations(ilevel-1),p_Rtriangulations(ilevel))
      call tria_initStandardMeshFromRaw (p_Rtriangulations(ilevel))
    end do
    
    ! =================================
    ! Create a hierarchy of Q1 discretisations
    ! =================================

    ! On all levels, create a scalar and a block discretisation for a Q1 block.    
    do ilevel = 1,NLMAX
      ! Create a spatial discretisation with Q1
      call spdiscr_initDiscr_simple (&
          p_RspatialDiscr(ilevel),EL_Q1_2D,p_Rtriangulations(ilevel))
      
      ! Create a block discretisation with 1 block Q1.
      call spdiscr_initBlockDiscr (p_RblockDiscr(ilevel),p_Rtriangulations(ilevel))
      call spdiscr_appendBlockComponent (p_RblockDiscr(ilevel),p_RspatialDiscr(ilevel))
      call spdiscr_commitBlockDiscr (p_RblockDiscr(ilevel))
    end do

    ! =================================
    ! Create a hierarchy of cubature
    ! structures. Use Gauss 3x3 on all
    ! levels.
    ! =================================

    do ilevel = 1,NLMAX

      ! Use a Gauss 3x3 formula for the discretisation.
      call spdiscr_createDefCubStructure (&
          p_RspatialDiscr(ilevel),p_RcubatureInfo(ilevel),CUB_GEN_AUTO_G3)

    end do

    ! =================================
    ! Assemble matrices on all levels
    ! =================================

    do ilevel = 1,NLMAX

      ! Create a matrix on level ilevel
      call lsysbl_createMatrix (p_RblockDiscr(ilevel),p_Rmatrices(ilevel))

      ! Create an empty matrix in block (1,1) in CSR format
      call bilf_createMatrixStructure (p_Rmatrices(ilevel),1,1,LSYSSC_MATRIX9)
      call lsysbl_allocEmptyMatrix (p_Rmatrices(ilevel))

      ! -----------------------------------------------------
      ! Discretise a simple Laplace matrix into (1,1)
      call lsysbl_clearMatrix (p_Rmatrices(ilevel))
      call bma_buildMatrix (p_Rmatrices(ilevel),BMA_CALC_STANDARD,&
          bma_fcalc_laplace,rcubatureInfo=p_RcubatureInfo(ilevel))

    end do

    ! =================================
    ! Assemble a RHS and create 
    ! an empty solution vector.
    ! =================================

    ! -----------------------------------------------------
    ! Initialise RHS/temp/solution vectors on the topmost level
    call lsysbl_createVector (p_RblockDiscr(NLMAX),rrhs)
    call lsysbl_createVector (p_RblockDiscr(NLMAX),rsolution)
    call lsysbl_createVector (p_RblockDiscr(NLMAX),rtemp)

    ! Create a RHS vector into block (1)
    call lsysbl_clearVector (rrhs)
    call bma_buildVector (rrhs,BMA_CALC_STANDARD,&
        bma_fcalc_rhsBubble,rcubatureInfo=p_RcubatureInfo(NLMAX))

    ! =================================
    ! Discretise boundary conditions
    ! =================================
    
    do ilevel = 1,NLMAX
    
      ! Initialise a boundary condition structure
      call bcasm_initDiscreteBC(p_RdiscreteBC(ilevel))
      
      ! Get a mesh region for the complete boundary
      call mshreg_createFromNodalProp(rmeshRegion, &
          p_Rtriangulations(ilevel), MSHREG_IDX_ALL)
      
      ! Discretise Dirichlet boundary conditions
      call bcasm_newDirichletBConMR (p_RblockDiscr(ilevel), 1, &
          p_RdiscreteBC(ilevel), rmeshRegion, fgetBoundaryValuesMR)
          
      call mshreg_done(rmeshregion)

    end do

    ! =================================
    ! Impose boundary conditions
    ! =================================

    do ilevel = 1,NLMAX
      ! Impose the BC into the matrix
      call matfil_discreteBC (p_Rmatrices(ilevel),p_RdiscreteBC(ilevel))
    end do
    
    ! On the topmost level, impose to the RHS and the soluiton
    call vecfil_discreteBCrhs (rrhs,p_RdiscreteBC(NLMAX))
    call vecfil_discreteBCsol (rsolution,p_RdiscreteBC(NLMAX))
    
    ! =================================
    ! Create prolongation/restriction
    ! matrices for level change in MG.
    ! =================================
    
    ! Start with ilevel=2 which describes the solution exchange
    ! between level 1 and 2.
    do ilevel = 2,NLMAX
    
      ! Create the prolongation matrix structure
      call mlop_create2LvlMatrixStruct (p_RspatialDiscr(ilevel-1),&
          p_RspatialDiscr(ilevel), LSYSSC_MATRIX9, p_RprolMatrices(ilevel))
          
      ! and assemble the entries of the prolongation matrix.
      call mlop_build2LvlProlMatrix (&
          p_RspatialDiscr(ilevel-1),p_RspatialDiscr(ilevel),&
          .true., p_RprolMatrices(ilevel),&
          rcubatureInfoCoarse=p_RcubatureInfo(ilevel-1),&
          rcubatureInfoFine=p_RcubatureInfo(ilevel) )

      ! Transpose the matrix, this gives the restriction matrix
      call lsyssc_transposeMatrix (p_RprolMatrices(ilevel),p_RrestMatrices(ilevel))
    
    end do

    ! =================================
    ! Solve the system with Multigrid elimination
    ! =================================
    
    call output_line ("Solving linear system...")
    
    ! ----------------
    ! Solver preparation
    ! ----------------

    ! Initialise a multigrid solver, NLMAX levels
    call linsol_initMultigrid2 (p_rsolverNode,NLMAX)
    
    ! On level 1, add a Gauss elimination solver as coarse grid solver.

    ! On level 2..NLMAX, add Jacobi as pre- and postsmoother, 4 smoothing steps.
    do ilevel = 1, NLMAX
    
      ! Get the mutigrid level data
      call linsol_getMultigrid2Level (p_rsolverNode,ilevel,p_rlevelInfo)
    
      if (ilevel .eq. 1) then
      
        ! ---------------------------------------
        ! Create a coarse grid solver
        ! ---------------------------------------

        ! Create UMFPACK
        call linsol_initUMFPACK4 (p_rcoarsegridsolver)
        
        ! Set as coarse grid solver.
        p_rlevelInfo%p_rcoarseGridSolver => p_rcoarsegridsolver
      
      else
      
        ! ---------------------------------------
        ! Create a smoother
        ! ---------------------------------------

        ! Create Jacobi
        call linsol_initJacobi (p_rsmoother)
        
        ! Configure as smoother, 4 steps, damping parameter 0.7
        call linsol_convertToSmoother (p_rsmoother,4,0.7_DP)
        
        ! Set as pre- and postsmoother
        p_rlevelInfo%p_rpresmoother => p_rsmoother
        p_rlevelInfo%p_rpostsmoother => p_rsmoother
      
        ! Create and attach a filter chain that filters defect vectors on 
        ! this level to include boundary conditions.
        call filter_initFilterChain (p_RfilterChain(:,ilevel),Isize(ilevel))
        call filter_newFilterDiscBCDef (p_RfilterChain(:,ilevel),Isize(ilevel),p_RdiscreteBC(ilevel))
        
        p_rlevelInfo%p_RfilterChain => p_RfilterChain(:,ilevel)
      
        ! ---------------------------------------
        ! Prepare a matrix-based prolongation/
        ! restriction
        ! ---------------------------------------

        ! Create a projection structure that describes prolongation/restriction
        ! for all blocks.
        call mlprj_initProjectionDiscr (p_Rprojection(ilevel),p_RblockDiscr(ilevel))
        
        ! Initialise prolongation/restriction for the first block (our
        ! Q1 solution) according to the prolongation/restriction matrices.
        call mlprj_initMatrixProjection(p_Rprojection(ilevel)%RscalarProjection(1,1),&
            p_RprolMatrices(ilevel),rmatrixRest=p_RrestMatrices(ilevel))
            
        ! Tell MG about the projection we want to use.
        p_rlevelInfo%p_rprojection => p_Rprojection(ilevel)
        
      end if
    end do
    
    ! Attach the system matrices
    call linsol_newMatrixSet (rmatrixSet)
    do ilevel=1,NLMAX
      call linsol_addMatrix (rmatrixSet,p_Rmatrices(ilevel))
    end do
    
    call linsol_setMatrices (p_rsolverNode, rmatrixSet)
    
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

    ! ----------------
    ! Set multigrid specific settings
    ! ----------------
    
    p_rsolverNode%p_rsubnodeMultigrid2%icycle = 0   ! 0=F-cycle, 1=V-cycle, 2=W-cycle

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
    
    ! Matrix set
    call linsol_releaseMatrixSet (rmatrixSet)
    
    ! Remaining solver data
    call linsol_releaseSolver (p_rsolverNode)

    ! =================================
    ! Create output files
    ! =================================

    call output_line ("Writing postprocessing files...")

    ! Open / write / close; write the solution to a VTK file.
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_Rtriangulations(NLMAX),&
        "post/tutorial017d.vtk")
    call ucd_addVectorByVertex (rexport, "solution", &
        UCD_VAR_STANDARD, rsolution%RvectorBlock(1))
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release vectors
    call lsysbl_releaseVector (rtemp)
    call lsysbl_releaseVector (rrhs)
    call lsysbl_releaseVector (rsolution)
    
    ! Release the matrices/discretisation structures/BC/projection structurees
    do ilevel=1,NLMAX
      call mlprj_doneProjection(p_Rprojection(ilevel))
      call lsyssc_releaseMatrix (p_RrestMatrices(ilevel))
      call lsyssc_releaseMatrix (p_RprolMatrices(ilevel))

      call filter_doneFilterChain (p_RfilterChain(:,ilevel),Isize(ilevel))
      call lsysbl_releaseMatrix (p_Rmatrices(ilevel))

      call bcasm_releaseDiscreteBC (p_RdiscreteBC(ilevel))

      call spdiscr_releaseCubStructure (p_RcubatureInfo(ilevel))
      call spdiscr_releaseBlockDiscr (p_RblockDiscr(ilevel))
      call spdiscr_releaseDiscr (p_RspatialDiscr(ilevel))
      call tria_done (p_Rtriangulations(ilevel))
    end do

    deallocate (p_Rprojection)
    deallocate (p_RrestMatrices)
    deallocate (p_RprolMatrices)

    deallocate (Isize)
    deallocate (p_RfilterChain)
    deallocate (p_Rmatrices)
    deallocate (p_RdiscreteBC)
    deallocate (p_RcubatureInfo)
    deallocate (p_RblockDiscr)
    deallocate (p_RspatialDiscr)
    deallocate (p_Rtriangulations)
    
  end subroutine

end module
