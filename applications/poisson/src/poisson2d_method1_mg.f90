!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method1_mg </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!#
!# This module is based on poisson2d_method0_simple, but using a multi-grid
!# solver.
!# </purpose>
!##############################################################################

module poisson2d_method1_mg

  use fsystem
  use genoutput
  use storage
  use linearsolver
  use boundary
  use derivatives
  use element
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use filtersupport
  use linearsystemscalar
  use linearsystemblock
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use scalarpde
  use ucd
  use pprocerror
  use collection
    
  use poisson2d_callback
  
  implicit none

!<types>

!<typeblock description="Type block defining all information about one level">

  type t_level
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo) :: rcubatureInfo
    
    ! A system matrix for that specific level. The matrix will receive the
    ! discrete Laplace operator.
    type(t_matrixBlock) :: rmatrix

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC
  
    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1) :: RfilterChain
    
    ! Number of filters in the filter chain
    integer :: nfilters

  end type
  
!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine poisson2d_1_mg
  
!<description>
  ! This is an all-in-one poisson solver for directly solving a Poisson
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  ! 1.) Read in parametrisation
  ! 2.) Read in triangulation
  ! 3.) Set up RHS
  ! 4.) Set up matrix
  ! 5.) Create solver structure
  ! 6.) Solve the problem
  ! 7.) Write solution to VTK file
  ! 8.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let us see...
    !
    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels

    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_vectorBlock) :: rvecSol,rvecRhs,rvecTmp

    ! A variable that is used to specify a region on the boundary.
    type(t_boundaryRegion) :: rboundaryRegion

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rcoarseGridSolver,p_rsmoother

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices

    ! One level of multigrid
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo

    ! NLMIN receives the level of the coarse grid.
    integer :: NLMIN

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Error of FE function to reference function
    real(DP) :: derror
    
    ! Output block for UCD output to VTK file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    
    ! Some temporary variables
    integer :: i

    ! Ok, let us start.
    !
    ! We want to solve our Poisson problem on level...
    NLMIN = 1
    NLMAX = 7
    
    ! Allocate memory for all levels
    allocate(Rlevels(NLMIN:NLMAX))

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Read the domain, read the mesh, refine
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, trim(spredir)//"/QUAD.prm")
        
    ! Now read in the basic triangulation into our coarse level.
    call tria_readTriFile2D (Rlevels(NLMIN)%rtriangulation, &
                             trim(spredir)//"/QUAD.tri", rboundary)
    
    ! Refine it.
    call tria_quickRefine2LevelOrdering (NLMIN-1,&
        Rlevels(NLMIN)%rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs
    ! from a triangulation.
    call tria_initStandardMeshFromRaw (Rlevels(NLMIN)%rtriangulation,&
        rboundary)
    
    ! Now refine the grid for the fine levels.
    do i = NLMIN+1, NLMAX

      ! Refine the grid using the 2-Level-Ordering algorithm
      call tria_refine2LevelOrdering(Rlevels(i-1)%rtriangulation,&
          Rlevels(i)%rtriangulation,rboundary)
      
      ! Create a standard mesh
      call tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation,&
          rboundary)
    
    end do

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up discretisation structures which tells the code which
    ! finite element to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    ! Do this for all levels
    do i = NLMIN, NLMAX
      call spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation, 1, &
                                   Rlevels(i)%rtriangulation, rboundary)
    end do
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! for this solution component:
    do i = NLMIN, NLMAX
      call spdiscr_initDiscr_simple (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          EL_Q1,Rlevels(i)%rtriangulation, rboundary)
    end do

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up an cubature info structure to tell the code which cubature
    ! formula to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Create an assembly information structure on each level which tells the code
    ! the cubature formula to use. Standard: Gauss 3x3.
    do i = NLMIN, NLMAX
      call spdiscr_createDefCubStructure(&  
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),Rlevels(i)%rcubatureInfo,&
          CUB_GEN_AUTO_G3)
    end do

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create a 1x1 block matrix with the operator on every level
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    do i = NLMIN, NLMAX

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (&
          Rlevels(i)%rdiscretisation,Rlevels(i)%rmatrix)

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      call bilf_createMatrixStructure ( &
           Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
           LSYSSC_MATRIX9,Rlevels(i)%rmatrix%RmatrixBlock(1,1))
      
      ! And now to the entries of the matrix. For assembling of the entries,
      ! we need a bilinear form, which first has to be set up manually.
      ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
      ! scalar system matrix in 2D.
      rform%itermCount = 2
      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_DERIV_X
      rform%Idescriptors(1,2) = DER_DERIV_Y
      rform%Idescriptors(2,2) = DER_DERIV_Y

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%Dcoefficients(1)  = 1.0
      rform%Dcoefficients(2)  = 1.0

      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_Laplace for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = .FALSE. above,
      ! the framework will call the callback routine to get analytical
      ! data.
      call bilf_buildMatrixScalar (rform,.true.,&
           Rlevels(i)%rmatrix%RmatrixBlock(1,1),Rlevels(i)%rcubatureInfo)
    
    end do
      
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create RHS and solution vectors
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Next step: Create a RHS vector, a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rvecRhs,.true.)
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rvecSol,.true.)
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rvecTmp,.true.)
      
    ! The vector structure is ready but the entries are missing.
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to get a discrete version of it.
    ! Again we simply create a scalar vector based on the one and only
    ! discretisation structure.
    ! This scalar vector will later be used as the one and only first
    ! component in a block vector.
    call linf_buildVectorScalar (&
        rlinform,.true.,rvecRhs%RvectorBlock(1),Rlevels(NLMAX)%rcubatureInfo,coeff_RHS_2D)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Discretise the boundary conditions and apply them to the matrix/RHS/sol.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    do i = NLMIN, NLMAX
    
      ! Initialise the discrete BC structure
      call bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBC)

      ! On edge 1 of boundary component 1 add Dirichlet boundary conditions.
      call boundary_createRegion(rboundary,1,1,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
                               
      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rboundary,1,2,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
      
      ! Edge 4 of boundary component 1. That is it.
      call boundary_createRegion(rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                        rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                        getBoundaryValues_2D)
      
      ! Also implement the boundary conditions into the matrix.
      call matfil_discreteBC (Rlevels(i)%rmatrix,Rlevels(i)%rdiscreteBC)
      
      ! During the linear solver, the boundary conditions are also
      ! frequently imposed to the vectors. But as the linear solver
      ! does not work with the actual solution vectors but with
      ! defect vectors instead.
      ! So, set up a filter chain that filters the defect vector
      ! during the solution process to implement discrete boundary conditions.
      call filter_clearFilterChain (Rlevels(i)%RfilterChain,Rlevels(i)%nfilters)
      call filter_newFilterDiscBCDef (&
          Rlevels(i)%RfilterChain,Rlevels(i)%nfilters,Rlevels(i)%rdiscreteBC)

    end do

    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    call vecfil_discreteBCrhs (rvecRhs,Rlevels(NLMAX)%rdiscreteBC)
    call vecfil_discreteBCsol (rvecSol,Rlevels(NLMAX)%rdiscreteBC)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a linear solver
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Now we have to build up the level information for multigrid.
    !
    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    call linsol_initMultigrid2 (p_rsolverNode,NLMAX-NLMIN+1)
    
    ! Set up a coarse grid solver.
    ! The coarse grid in multigrid is always grid 1!
    call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
    call linsol_initUMFPACK4 (p_rcoarseGridSolver)
    
    p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
    p_rlevelInfo%p_RfilterChain => Rlevels(1)%RfilterChain
    
    ! Now set up the other levels...
    do i = NLMIN+1, NLMAX
    
      ! Create a Jacobi smoother
      !CALL linsol_initJacobi(p_rsmoother)
      
      ! Create an ILU(0) smoother
      call linsol_initMILUs1x1 (p_rsmoother,0,0.0_DP)
      
      ! We will use 4 smoothing steps with damping parameter 0.7
      call linsol_convertToSmoother(p_rsmoother, 4, 0.7_DP)
      
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (p_rsolverNode,i-NLMIN+1,p_rlevelInfo)
      p_rlevelInfo%p_rpresmoother => p_rsmoother
      p_rlevelInfo%p_rpostsmoother => p_rsmoother

      ! Attach the filter chain which imposes boundary conditions on that level.
      p_rlevelInfo%p_RfilterChain => Rlevels(i)%RfilterChain
      
    end do
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2
    
    ! Attach the system matrices to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array. Note that this does not
    ! allocate new memory, we create only "links" to existing matrices
    ! into Rmatrices(:)!
    allocate(Rmatrices(NLMIN:NLMAX))
    do i = NLMIN, NLMAX
      call lsysbl_duplicateMatrix (Rlevels(i)%rmatrix,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    
    call linsol_setMatrices(p_RsolverNode,Rmatrices(NLMIN:NLMAX))

    ! We can release Rmatrices immediately -- as long as we do not
    ! release Rlevels(i)%rmatrix!
    do i=NLMIN,NLMAX
      call lsysbl_releaseMatrix (Rmatrices(i))
    end do
    deallocate(Rmatrices)

    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    call linsol_initStructure (p_rsolverNode, ierror)
    
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Matrix structure invalid!",OU_CLASS_ERROR)
      call sys_halt()
    end if

    call linsol_initData (p_rsolverNode, ierror)
    
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Matrix singular!",OU_CLASS_ERROR)
      call sys_halt()
    end if
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Solve the system
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    call linsol_solveAdaptively (p_rsolverNode,rvecSol,rvecRhs,rvecTmp)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Postprocessing of the solution
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! That is it, rvecSol now contains our solution. We can now
    ! start the postprocessing.
    !
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,&
        Rlevels(NLMAX)%rtriangulation,trim(sucddir)//"/u2d_1_mg.vtk")
    
    ! Add the solution to the UCD exporter
    call ucd_addVectorByVertex (rexport, "sol", UCD_VAR_STANDARD, &
        rvecSol%RvectorBlock(1))
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Calculate the error to the reference function.
    call pperr_scalar (PPERR_L2ERROR,derror,rvecSol%RvectorBlock(1),&
        getReferenceFunction_2D, rcubatureInfo=Rlevels(NLMAX)%rcubatureInfo)
    call output_line ("L2-error: " // sys_sdEL(derror,10) )

    call pperr_scalar (PPERR_H1ERROR,derror,rvecSol%RvectorBlock(1),&
        getReferenceFunction_2D, rcubatureInfo=Rlevels(NLMAX)%rcubatureInfo)
    call output_line ("H1-error: " // sys_sdEL(derror,10) )
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Clean up
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the block matrices/vectors
    call lsysbl_releaseVector (rvecTmp)
    call lsysbl_releaseVector (rvecSol)
    call lsysbl_releaseVector (rvecRhs)
    do i = NLMAX, NLMIN, -1
      call lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
    end do

    ! Release our discrete version of the boundary conditions
    do i = NLMAX, NLMIN, -1
      call bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
    end do

    ! Release the cubature info structures
    do i=NLMAX,NLMIN,-1
      call spdiscr_releaseCubStructure(Rlevels(i)%rcubatureInfo)
    end do

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    do i = NLMAX, NLMIN, -1
      call spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
    end do
    
    ! Release the triangulation.
    do i = NLMAX, NLMIN, -1
      call tria_done (Rlevels(i)%rtriangulation)
    end do
    
    deallocate(Rlevels)

    ! Finally release the domain, that is it.
    call boundary_release (rboundary)

  end subroutine

end module
