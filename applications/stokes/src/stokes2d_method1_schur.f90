!##############################################################################
!# ****************************************************************************
!# <name> stokes2d_method1_schur </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a Stokes
!# problem on a simple domain.
!#
!# This module is based on stokes2d_method1_mg, but in contrast, this module
!# uses a Schur-complement preconditioner.
!#
!# The global 2d stokes system looks like follows:
!#
!#                   ( A B ) * ( u ) = ( f )
!#                   ( D 0 )   ( p )   ( g )
!#
!# Now the global system above is solved using BiCGStab and a Schur-complement
!# preconditioner of the type
!#
!#                P = ( A^-1  0    )
!#                    ( D    -S^-1 )
!#
!# Where S is an approximation of the Schur-complement of A. In this case,
!# S is chosen as the mass matrix of the pressure space.
!#
!# In this example A^-1 is approximated using a multigrid solver.
!#
!# Remark:
!# The documentation of this example is a bit more compressed than the other
!# examples to highlight the important steps to create a Schur-complement
!# preconditioner.
!# </purpose>
!##############################################################################

module stokes2d_method1_schur

  use fsystem
  use storage
  use genoutput
  use boundary
  use cubature
  use derivatives
  use matrixfilters
  use vectorfilters
  use linearalgebra
  use discretebc
  use bcassembly
  use triangulation
  use element
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use spdiscprojection
  use scalarpde
  use bilinearformevaluation
  use linearformevaluation
  use discretebc
  use filtersupport
  use coarsegridcorrection
  use linearsolver
  use ucd
  use stdoperators
  use collection, only: t_collection
  use pprocerror
  use genoutput
  
  use stokes2d_callback
  
  implicit none

!<types>

!<typeblock description="Type block defining all information about one level">

  type t_level
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo) :: rcubatureInfo

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A system matrix for that specific level. The matrix will receive the
    ! discrete Laplace operator.
    type(t_matrixBlock) :: rmatrix

    ! B1-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB1

    ! B2-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB2

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

  subroutine stokes2d_1_schur
  
!<description>
  ! This is an all-in-one stokes solver for directly solving a stokes
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  ! 1.) Read in parametrisation
  ! 2.) Read in triangulation
  ! 3.) Set up RHS
  ! 4.) Set up matrix
  ! 5.) Create solver structure
  ! 6.) Solve the problem
  ! 7.) Write solution to GMV file
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

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_vectorBlock) :: rvector,rrhs,rtempBlock
    
    ! A set of variables describing the analytic and discrete boundary
    ! conditions.
    type(t_boundaryRegion) :: rboundaryRegion

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsubsolverA, p_rsubsolverS, p_rsolverNode,&
                                   p_rsmoother, p_rsolverSchur

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices
    
    ! A discretisation for the Schur-complement matrix
    type(t_blockDiscretisation) :: rdiscrS
    
    ! An array for the Schur-complement matrices
    type(t_matrixBlock), dimension(1) :: RmatrixS

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain
    
    ! One level of multigrid
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    
    ! NLMIN receives the level of the coarse grid.
    integer :: NLMIN

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Viscosity parameter nu = 1/Re
    real(DP) :: dnu
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Du1,p_Du2, p_Dp

    ! A counter variable
    integer :: i

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! A collection structure for post-processing
    type(t_collection) :: rcollection

    ! Error data structures for post-processing
    type(t_errorScVec) :: rerrorU, rerrorP

    ! Error arrays for post-processing
    real(DP), dimension(2), target :: DerrorUL2, DerrorUH1
    real(DP), dimension(1), target :: DerrorPL2

    ! Ok, let us start.
    !
    ! We want to solve our Stokes problem on level...
    NLMIN = 2
    NLMAX = 7
    
    ! Viscosity parameter:
    dnu = 1.0_DP

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Read the domain, read the mesh, refine
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Allocate memory for all levels
    allocate(Rlevels(NLMIN:NLMAX))

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, trim(spredir)//"/QUAD.prm")
        
    ! Now read in the basic triangulation into our coarse level.
    call tria_readTriFile2D (Rlevels(NLMIN)%rtriangulation, &
                             trim(spredir)//"/QUAD.tri", rboundary)
    
    ! Refine the mesh up to the minimum level
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
    ! Set up a discretisation structure which tells the code which
    ! finite element to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies 3 blocks in the
    ! solution vector.
    do i = NLMIN, NLMAX
      call spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation, 3, &
                                   Rlevels(i)%rtriangulation, rboundary)
    end do

    ! Set up the block discretisation for (u1, u2, p)
    do i = NLMIN, NLMAX
      call spdiscr_initDiscr_simple (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          EL_EM30_NEW, Rlevels(i)%rtriangulation, rboundary)
      call spdiscr_duplicateDiscrSc (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(2))
      call spdiscr_deriveSimpleDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          EL_Q0, Rlevels(i)%rdiscretisation%RspatialDiscr(3))
    
    end do

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up an cubature info structure to tell the code which cubature
    ! formula to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
                 
    do i = NLMIN, NLMAX
      ! Create an cubature information structure which tells the code
      ! the cubature formula to use. Standard: Gauss 3x3.
      call spdiscr_createDefCubStructure(&  
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          Rlevels(i)%rcubatureInfo,CUB_GEN_AUTO_G3)
    end do

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create a 3x3 block matrix with the operator
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    do i = NLMIN, NLMAX
    
      ! Initialise a saddle-point block matrix
      call lsysbl_createMatBlockByDiscr (Rlevels(i)%rdiscretisation,&
                                         Rlevels(i)%rmatrix)
      Rlevels(i)%rmatrix%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT
      
      ! Create matrix structure for A
      call bilf_createMatrixStructure (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrix%RmatrixBlock(1,1))

      ! Create matrix structures for B1/B2
      call bilf_createMatrixStructure (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(3),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrixB1,&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1))

      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, Rlevels(i)%rmatrixB2,&
          LSYSSC_DUP_COPY, LSYSSC_DUP_REMOVE)
      
      ! Assemble a Laplace matrix for A1
      call stdop_assembleLaplaceMatrix(Rlevels(i)%rmatrix%RmatrixBlock(1,1),.true.,dnu,&
          rcubatureInfo=Rlevels(i)%rcubatureInfo)
      
      ! Create A2 as a shared copy of A1
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
          Rlevels(i)%rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
      ! Update discretisation of A2
      call lsyssc_assignDiscrDirectMat (Rlevels(i)%rmatrix%RmatrixBlock(2,2),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(2))
      
      ! Assemble B1/B2 matrices
      call stdop_assembleSimpleMatrix (Rlevels(i)%rmatrixB1,&
          DER_FUNC2D, DER_DERIV2D_X, -1.0_DP, rcubatureInfo=Rlevels(i)%rcubatureInfo)
      call stdop_assembleSimpleMatrix (Rlevels(i)%rmatrixB2,&
          DER_FUNC2D, DER_DERIV2D_Y, -1.0_DP, rcubatureInfo=Rlevels(i)%rcubatureInfo)
      
      ! Copy B1/B2 into the block matrix
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB1, &
          Rlevels(i)%rmatrix%RmatrixBlock(1,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrixB2, &
          Rlevels(i)%rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      
      ! Furthermore, put B1^T and B2^T to the block matrix.
      call lsyssc_transposeMatrix (Rlevels(i)%rmatrixB1, &
          Rlevels(i)%rmatrix%RmatrixBlock(3,1),LSYSSC_TR_VIRTUAL)

      call lsyssc_transposeMatrix (Rlevels(i)%rmatrixB2, &
          Rlevels(i)%rmatrix%RmatrixBlock(3,2),LSYSSC_TR_VIRTUAL)

    end do
    
    ! -------------------------------------------------------------------------
    ! Set up the Schur-complement approximation matrix
    ! -------------------------------------------------------------------------
    
    ! Derive a block discretisation which only consists of the pressure.
    ! We need to do this because the matrix S needs a block discretisation.
    call spdiscr_deriveBlockDiscr(Rlevels(NLMAX)%rdiscretisation, rdiscrS, 3, 3)
    
    ! Create a block matrix for the Schur-complement based on the previously
    ! created discretisation.
    call lsysbl_createMatBlockByDiscr(rdiscrS, rmatrixS(1))

    ! Create the matrix structure of the Schur-complement matrix
    call bilf_createMatrixStructure (rdiscrS%RspatialDiscr(1),&
        LSYSSC_MATRIX9, RmatrixS(1)%RmatrixBlock(1,1))
    
    ! Assemble a mass matrix as an approximation of the Schur-complement
    ! of the Laplace matrix blocks of the system matrix.
    call stdop_assembleSimpleMatrix (RmatrixS(1)%RmatrixBlock(1,1),&
        DER_FUNC2D, DER_FUNC2D, rcubatureInfo=Rlevels(NLMAX)%rcubatureInfo)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create RHS and solution vectors
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Next step: Create a RHS vector, a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rrhs,.true.)
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rvector,.true.)
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rtempBlock,.true.)

    ! Remark:
    ! As in this case the right-hand-side is zero everywhere, we will not call
    ! the linearformevaluation here to compress the code a bit.

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Assembly of matrices/vectors finished
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Discretise the boundary conditions and apply them to the matrix/RHS/sol.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Now implement the boundary conditions
    do i = NLMIN, NLMAX

      ! Create a BC structure
      call bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBC)
      
      ! Add boundary conditions for X- and Y-velocity
      call boundary_createRegion(rboundary,1,1,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
          rboundaryRegion,Rlevels(i)%rdiscreteBC,getBoundaryValues_2D)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,2,&
          rboundaryRegion,Rlevels(i)%rdiscreteBC,getBoundaryValues_2D)
      call boundary_createRegion(rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
          rboundaryRegion,Rlevels(i)%rdiscreteBC,getBoundaryValues_2D)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,2,&
          rboundaryRegion,Rlevels(i)%rdiscreteBC,getBoundaryValues_2D)
      call boundary_createRegion(rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
          rboundaryRegion,Rlevels(i)%rdiscreteBC,getBoundaryValues_2D)
      call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,2,&
          rboundaryRegion,Rlevels(i)%rdiscreteBC,getBoundaryValues_2D)

      ! Next step is to implement boundary conditions into the matrix.
      ! This is done using a matrix filter for discrete boundary conditions.
      ! The discrete boundary conditions are already attached to the
      ! matrix. Call the appropriate matrix filter that modifies the matrix
      ! according to the boundary conditions.
      call matfil_discreteBC (Rlevels(i)%rmatrix,Rlevels(i)%rdiscreteBC)
    
      ! Create a filter chain for the solver that implements boundary conditions
      ! during the solution process.
      call filter_clearFilterChain (Rlevels(i)%RfilterChain,Rlevels(i)%nfilters)
      call filter_newFilterDiscBCDef (&
          Rlevels(i)%RfilterChain,Rlevels(i)%nfilters,Rlevels(i)%rdiscreteBC)

    end do

    ! Also implement the discrete boundary conditions on the finest level
    ! onto our right-hand-side and solution vectors.

    call vecfil_discreteBCrhs (rrhs,Rlevels(NLMAX)%rdiscreteBC)
    call vecfil_discreteBCsol (rvector,Rlevels(NLMAX)%rdiscreteBC)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a linear solver
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Create a filter chain
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! -------------------------------------------------------------------------
    ! Set up subsolver for submatrix A
    ! -------------------------------------------------------------------------
    
    ! Create a Multigrid-solver for the A-submatrix.
    ! Please note that the filter chain is not passed to the multigrid solver,
    ! as the filter-chain belongs to the global system and not only to the
    ! A-submatrices!
    ! (Todo: Test whether multigrid survives without a filter chain in the
    !        case of "harder" problems, too...)
    call linsol_initMultigrid2 (p_rsubsolverA,NLMAX-NLMIN+1)
    
    ! Set up a coarse grid solver.
    ! The coarse grid in multigrid is always grid 1!
    call linsol_getMultigrid2Level (p_rsubsolverA,1,p_rlevelInfo)
    call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
    p_rlevelInfo%p_rfilterChain => Rlevels(NLMIN)%RfilterChain

    ! Now set up the other levels...
    do i = NLMIN+1, NLMAX
    
      ! Set up a Jacobi smoother.
      call linsol_initJacobi (p_rsmoother)
      
      ! We will use 4 smoothing steps with damping parameter 0.7
      call linsol_convertToSmoother(p_rsmoother, 4, 0.7_DP)

      ! Put the smoother into the level info structure as presmoother
      ! and postsmoother
      call linsol_getMultigrid2Level (p_rsubsolverA,i-NLMIN+1,p_rlevelInfo)
      p_rlevelInfo%p_rpresmoother => p_rsmoother
      p_rlevelInfo%p_rpostsmoother => p_rsmoother
      p_rlevelInfo%p_rfilterChain => Rlevels(i)%RfilterChain
      
    end do
    
    ! Make sure that multigrid always performs 3 iterations, disable
    ! the residual check and set its output level to -1 to avoid that
    ! multigrid compains about "accuracy warnings".
    p_rsubsolverA%nmaxiterations = 3
    p_rsubsolverA%irescheck = NO
    p_rsubsolverA%ioutputlevel = -1
    
    ! -------------------------------------------------------------------------
    ! Set up subsolver for Schur-complement matrix S
    ! -------------------------------------------------------------------------

    ! Create a jacobi solver for the Schur-complement matrix S.
    ! As S is a diagonal matrix, Jacobi will easily do the job.
    call linsol_initJacobi(p_rsubsolverS)
    
    ! -------------------------------------------------------------------------
    ! Set up system solver using a Schur-complement preconditioner
    ! -------------------------------------------------------------------------

    ! Create a Schur-complement preconditioner.
    ! This recieves the subsolvers for the A-submatrix and the Schur-complement
    ! matrix S and the Schur-complement matrix S itself.
    call linsol_initSchur(p_rsolverSchur, p_rsubsolverA, p_rsubsolverS,&
                          RmatrixS, LINSOL_SCHUR_TYPE_LTRI)
    
    ! Finally, create a BiCGStab solver for the global system.
    call linsol_initBiCGStab(p_rsolverNode, p_rsolverSchur, Rlevels(NLMAX)%RfilterChain)
   
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! Attach the system matrices to the solver.
    allocate(Rmatrices(NLMIN:NLMAX))
    do i = NLMIN, NLMAX
      call lsysbl_duplicateMatrix (Rlevels(i)%rmatrix,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    call linsol_setMatrices(p_rsolverNode,Rmatrices(NLMIN:NLMAX))
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

    ! Finally solve the system.
    call linsol_solveAdaptively (p_rsolverNode,rvector,rrhs,rtempBlock)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Postprocessing of the solution
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Project the solution sub-vectors for the output
    nullify(p_Du1)
    nullify(p_Du2)
    nullify(p_Dp)
    call spdp_projectToVertices(rvector%RvectorBlock(1), p_Du1)
    call spdp_projectToVertices(rvector%RvectorBlock(2), p_Du2)
    call spdp_projectToCells(rvector%RvectorBlock(3), p_Dp)

    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

    ! Write velocity and pressure to VTK file
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,&
        Rlevels(NLMAX)%rtriangulation,trim(sucddir)//"/u2d_1_schur.vtk")
    
    call ucd_addVarVertBasedVec(rexport, "velocity", p_Du1, p_Du2)
    call ucd_addVariableElementBased(rexport, "pressure", UCD_VAR_STANDARD, p_Dp)
    
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    deallocate(p_Dp)
    deallocate(p_Du2)
    deallocate(p_Du1)

    ! Store the viscosity parameter nu in the collection"s quick access array
    rcollection%DquickAccess(1) = dnu

    ! Set up the error structure for velocity
    rerrorU%p_RvecCoeff => rvector%RvectorBlock(1:2)
    rerrorU%p_DerrorL2 => DerrorUL2
    rerrorU%p_DerrorH1 => DerrorUH1

    ! Set up the error structure for pressure
    rerrorP%p_RvecCoeff => rvector%RvectorBlock(3:3)
    rerrorP%p_DerrorL2 => DerrorPL2

    ! Calculate errors of velocity and pressure against analytic solutions.
    call pperr_scalarVec(rerrorU, funcVelocity2D, rcollection, Rlevels(NLMAX)%rcubatureInfo)
    call pperr_scalarVec(rerrorP, funcPressure2D, rcollection, Rlevels(NLMAX)%rcubatureInfo)

    ! Print the errors.
    call output_lbrk()
    call output_line("|u - u_h|_L2 = " // trim(sys_sdEL(DerrorUL2(1), 10)) &
                                // " " // trim(sys_sdEL(DerrorUL2(2), 10)))
    call output_line("|u - u_h|_H1 = " // trim(sys_sdEL(DerrorUH1(1), 10)) &
                                // " " // trim(sys_sdEL(DerrorUH1(2), 10)))
    call output_line("|p - p_h|_L2 = " // trim(sys_sdEL(derrorPL2(1), 10)))

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Clean up
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    call linsol_releaseSolver (p_rsolverNode)
    call lsysbl_releaseVector (rtempBlock)
    call lsysbl_releaseVector (rvector)
    call lsysbl_releaseVector (rrhs)
    
    ! Release the Schur-complement matrix
    call lsysbl_releaseMatrix(RmatrixS(1))
    
    ! Release block discretisation of Schur-complement matrix
    call spdiscr_releaseBlockDiscr(rdiscrS)
    
    do i = NLMAX, NLMIN, -1
      call lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
      call lsyssc_releaseMatrix (Rlevels(i)%rmatrixB2)
      call lsyssc_releaseMatrix (Rlevels(i)%rmatrixB1)
      call spdiscr_releaseCubStructure(Rlevels(i)%rcubatureInfo)
      call bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
      call spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
      call tria_done (Rlevels(i)%rtriangulation)
    end do
    
    deallocate(Rlevels)
    
    call boundary_release (rboundary)

  end subroutine

end module
