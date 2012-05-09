!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method1_l2prj </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!#
!# This module is based on poisson2d_method1_mg, but using the EB30 element
!# (Q1~ with bubble) and L2-projection for multi-level operations and
!# post-processing, i.e. prolongation, restiction and the spatial projection
!# for the UCD output.
!#
!# </purpose>
!##############################################################################

module poisson2d_method1_l2prj

  use fsystem
  use genoutput
  use storage
  use linearalgebra
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use linearsystemscalar
  use linearsystemblock
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use ucd
  use scalarpde
  use pprocerror
  use genoutput
  use stdoperators
  use coarsegridcorrection
  use filtersupport
  use multileveloperators
  use multilevelprojection
    
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
    
    ! A mass matrix for that specific level. This matrix is needed by the
    ! L2-projection multi-level operators and is defined on each level except
    ! for the coarse-most one.
    type(t_matrixScalar) :: rmatrixMass
    
    ! A 2-Level-Mass matrix for that specific level. This matrix is needed by
    ! L2-projection multi-level operators and is defined on each level except
    ! for the coarse-most one.
    type(t_matrixScalar) :: rmatrix2Lvl

    ! An interlevel-projection structure for that specific level.
    type(t_interlevelProjectionBlock) :: rprojection

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC
  
  end type
  
!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine poisson2d_1_l2prj
  
!<description>
  ! This is an all-in-one poisson solver for directly solving a Poisson
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  !  1.) Read in parametrisation
  !  2.) Read in triangulation
  !  3.) Set up RHS
  !  4.) Set up system matrix
  !  5.) Set up 2-Level-projection operators for multi-grid
  !  6.) Create solver structure
  !  7.) Solve the problem
  !  8.) Perform L2-Projection of the solution to Q1 space
  !  9.) Write solution to VTK file
  ! 10.) Release all variables, finish
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
    type(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock

    ! A variable that is used to specify a region on the boundary.
    type(t_boundaryRegion) :: rboundaryRegion

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rsmoother

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices

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
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! One spatial discretisation for the L2-projection of the solution
    type(t_spatialDiscretisation) :: rdiscrQ1
    
    ! Two scalar matrices for the L2-projection
    type(t_matrixScalar) :: rmatrixMassPrj, rlumpedMassPrj
    
    ! Three scalar vectors for the L2-projection
    type(t_vectorScalar) :: rvecSolQ1, rvecRhsQ1, rvecDefQ1
    
    ! A real for the defect-correction-loop
    real(DP) :: dres
    
    ! Error of FE function to reference function
    real(DP) :: derror
    
    ! Output block for UCD output to VTK file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata

    ! A simple counter variable
    integer :: i
    
    ! Ok, let us start.
    !
    ! We want to solve our Poisson problem on level...
    NLMIN = 2
    NLMAX = 6
    
    ! Allocate memory for all levels
    allocate(Rlevels(NLMIN:NLMAX))
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Read the domain, read the mesh, refine
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = './pre'

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, trim(spredir)//'/QUAD.prm')
        
    ! Now read in the basic triangulation into our coarse level.
    call tria_readTriFile2D (Rlevels(NLMIN)%rtriangulation, &
                             trim(spredir)//'/QUAD.tri', rboundary)
    
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
          EL_EB30,Rlevels(i)%rtriangulation, rboundary)
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
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = 1.0
      rform%Dcoefficients(2)  = 1.0

      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_Laplace for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
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
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rrhsBlock,.true.)
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rvectorBlock,.true.)
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rtempBlock,.true.)
      
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
        rlinform,.true.,rrhsBlock%RvectorBlock(1),Rlevels(NLMAX)%rcubatureInfo,coeff_RHS_2D)
    
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
      
      ! Assign the BC`s to the matrix. That way, these
      ! boundary conditions are always connected to that matrix.
      call lsysbl_assignDiscreteBC(Rlevels(i)%rmatrix,Rlevels(i)%rdiscreteBC)
  
      ! Also implement the boundary conditions into the matrix.
      call matfil_discreteBC (Rlevels(i)%rmatrix)
      
    end do

    ! Our right-hand-side/solution/temp vectors also needs to 
    ! know the boundary conditions.
    call lsysbl_assignDiscreteBC(rrhsBlock,Rlevels(NLMAX)%rdiscreteBC)
    call lsysbl_assignDiscreteBC(rvectorBlock,Rlevels(NLMAX)%rdiscreteBC)
    call lsysbl_assignDiscreteBC(rtempBlock,Rlevels(NLMAX)%rdiscreteBC)

    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors/matrix. Call the appropriate vector/matrix filter that
    ! modifies the vectors/matrix according to the boundary conditions.
    call vecfil_discreteBCrhs (rrhsBlock)
    call vecfil_discreteBCsol (rvectorBlock)
    
    ! During the linear solver, the boundary conditions are also
    ! frequently imposed to the vectors. But as the linear solver
    ! does not work with the actual solution vectors but with
    ! defect vectors instead.
    ! So, set up a filter chain that filters the defect vector
    ! during the solution process to implement discrete boundary conditions.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up the L2-projection for Multigrid
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Now we need to set up the (2-Level) mass matrices for all levels except
    ! for the coarse-most one.
    do i = NLMIN+1, NLMAX
    
      ! Since the structure of the mass matrix is equal to the one of the
      ! Laplace matrix on the current level, we will simply create a shared
      ! copy of the structure.
      call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
          Rlevels(i)%rmatrixMass, LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      
      ! Assemble the mass matrix
      call stdop_assembleSimpleMatrix(Rlevels(i)%rmatrixMass,&
          DER_FUNC, DER_FUNC, rcubatureInfo=Rlevels(i)%rcubatureInfo)

      ! Now create the matrix structure of the 2-Level mass matrix.
      call mlop_create2LvlMatrixStruct(&
          Rlevels(i-1)%rdiscretisation%RspatialDiscr(1),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          LSYSSC_MATRIX9, Rlevels(i)%rmatrix2Lvl)
      
      ! And assemble the entries of the 2-Level mass matrix:
      call mlop_build2LvlMassMatrix (&
          Rlevels(i-1)%rdiscretisation%RspatialDiscr(1),&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
          .true., Rlevels(i)%rmatrix2Lvl,&
          Rlevels(i-1)%rcubatureInfo,Rlevels(i)%rcubatureInfo)
      
      ! Now set up an interlevel projecton structure for this level
      ! based on the Laplace matrix on this level.
      call mlprj_initProjectionMat (Rlevels(i)%rprojection,&
                                    Rlevels(i)%rmatrix)
      
      ! And initialise the L2-projection
      call mlprj_initL2Projection (&
          Rlevels(i)%rprojection%RscalarProjection(1,1),&
          Rlevels(i)%rmatrix2Lvl, Rlevels(i)%rmatrixMass)
      
      Rlevels(i)%rprojection%RscalarProjection(1,1)%depsL2 = 1e-20_DP
      
    end do

    ! And set up an interlevel projecton structure for the coarse-most level.
    call mlprj_initProjectionMat (Rlevels(NLMIN)%rprojection,&
                                  Rlevels(NLMIN)%rmatrix)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! L2-projection for Multigrid is set up now
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a linear solver
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Now we have to build up the level information for multigrid.
    !
    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    call linsol_initMultigrid2 (p_rsolverNode,NLMAX-NLMIN+1,RfilterChain)
    
    ! Set up a coarse grid solver.
    ! The coarse grid in multigrid is always grid 1!
    call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
    call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)

    ! Now set up the other levels...
    do i = NLMIN+1, NLMAX
    
      ! Create a Jacobi smoother
      call linsol_initJacobi(p_rsmoother)

      ! We will use 4 smoothing steps with damping parameter 0.7
      call linsol_convertToSmoother(p_rsmoother, 4, 0.7_DP)

      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (p_rsolverNode,i-NLMIN+1,p_rlevelInfo)
      p_rlevelInfo%p_rpresmoother => p_rsmoother
      p_rlevelInfo%p_rpostsmoother => p_rsmoother
      
      ! Attach our user-defined projection to the level.
      call linsol_initProjMultigrid2Level(p_rlevelInfo,Rlevels(i)%rprojection)
      
    end do
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2
    
    ! Set the coarse grid correction type to calculate an optimal
    ! damping parameter in repsect to the energy norm, as we apply
    ! the multigrid solver on a non-conforming element - this improves
    ! the multigrid convergence rate a bit.
    p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%ccorrectionType = &
      CGCOR_SCALARENERGYMIN
    
    ! Attach the system matrices to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array. Note that this does not
    ! allocate new memory, we create only 'links' to existing matrices
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
    call linsol_solveAdaptively (p_rsolverNode,rvectorBlock,rrhsBlock,rtempBlock)
      
    ! That is it, rvectorBlock now contains our solution.

    ! Calculate the error to the reference function.
    call pperr_scalar (PPERR_L2ERROR,derror,rvectorBlock%RvectorBlock(1),&
        getReferenceFunction_2D, rcubatureInfo=Rlevels(NLMAX)%rcubatureInfo)
    call output_line ('L2-error: ' // sys_sdEL(derror,10) )

    call pperr_scalar (PPERR_H1ERROR,derror,rvectorBlock%RvectorBlock(1),&
        getReferenceFunction_2D, rcubatureInfo=Rlevels(NLMAX)%rcubatureInfo)
    call output_line ('H1-error: ' // sys_sdEL(derror,10) )
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Perform L2-projection of solution into Q1 space
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    call output_lbrk()
    call output_line('Performing L2-projection of solution to Q1 space')
    call output_lbrk()

    ! We now have the solution vector, but unfortunately, it is a Q1~ solution
    ! vector and what we need are the function values in the vertices of the
    ! mesh. Instead of calling the spatial-projection module to interpolate
    ! the solution we will perform a "true" L2-projection of the Q1~ solution
    ! into the Q1 space, i.e. we want to "solve"
    !
    !                   v_h = u_h                                    (1)
    !
    ! with v_h in Q1 and u_h being our actual discrete Q1~ solution of the
    ! Poisson equation. So as a weak formulation of (1) we get
    !
    !          (v_h, phi_i) = (u_h, phi_i)  for all 1 <= i <= m      (2)
    !
    ! with (phi_1,...,phi_m) being the basis functions of Q1.
    ! Furthermore, we know that
    !
    !                          n
    !                   u_h = sum (x_k * psi_k)                      (3)
    !                         k=1
    !
    ! with (psi_1,...,psi_n) being the basis functions of Q1~ and
    ! x = (x_1,...,x_n) being our actual solution vector rvectorBlock.
    ! Now the discrete Q1 function v_h we search for also has a coefficient
    ! vector y = (y_1,...,y_m) and can be written as
    !
    !                          m
    !                   v_h = sum (y_j * phi_j)                      (4)
    !                         j=1
    !
    ! Inserting (3) and (4) into (2) and rewriting the result as a
    ! linear system we get
    !
    !                 M * y = N * x                                  (5)
    !
    !     for all 1 <= i,j <= m ; 1 <= k <= n :
    !
    !                 M_ij := (phi_j, phi_i)
    !                 N_ik := (psi_k, phi_i)
    !
    ! So there are 2 matrices we need to build: the Q1 mass matrix (M) and a
    ! mass matrix with Q1~ being its trial and Q1 being its test space (N).
    ! Afterwards we solve (5) and get the coefficient vector y of our
    ! Q1 solution v_h...
    !
    ! The first thing that we need is a Q1 discretisation on the fine mesh.
    call spdiscr_initDiscr_simple(rdiscrQ1, EL_Q1, &
                                  Rlevels(NLMAX)%rtriangulation, rboundary)

    ! Now create the the matrix structure of N.
    ! The trial space is EB30 and the test space is Q1:
    call bilf_createMatrixStructure (&
        Rlevels(NLMAX)%rdiscretisation%RspatialDiscr(1),&
        LSYSSC_MATRIX9, rmatrixMassPrj, rdiscrQ1)

    ! And assemble the mass matrix entries of N:
    call stdop_assembleSimpleMatrix(rmatrixMassPrj, DER_FUNC, DER_FUNC,&
        rcubatureInfo=Rlevels(NLMAX)%rcubatureInfo)

    ! We need to create the Q1 coefficient vector y:
    call lsyssc_createVecByDiscr (rdiscrQ1, rvecSolQ1, .true.)
    
    ! And we also need a Q1 rhs vector r which recieves N*x
    call lsyssc_createVecByDiscr (rdiscrQ1, rvecRhsQ1, .false.)
    
    ! Calculate r := N*x
    call lsyssc_scalarMatVec(rmatrixMassPrj, rvectorBlock%rvectorBlock(1),&
                             rvecRhsQ1, 1.0_DP, 0.0_DP)
    
    ! At this point we will not need the matrix N anymore, as we just needed
    ! it to get a rhs vector for our Q1 mass system - so we will release it now.
    call lsyssc_releaseMatrix(rmatrixMassPrj)

    ! As rmatrixMassPrj is free now, we will use it to store M
    call bilf_createMatrixStructure (rdiscrQ1,LSYSSC_MATRIX9,rmatrixMassPrj)
    call stdop_assembleSimpleMatrix(rmatrixMassPrj, DER_FUNC, DER_FUNC,&
        rcubatureInfo=Rlevels(NLMAX)%rcubatureInfo)
    
    ! We now have assembled the Q1 mass matrix M, created a Q1 solution
    ! vector y and calculated the right-hand-side vector r := N*x.
    ! Basically, we could now create a linear solver as CG or whatever, but
    ! this is not necessary as a system with a mass matrix is quite easy
    ! to solve (in contrast to a system with a Laplace matrix).
    ! We will solve the linear system by a Defect-Correction-Approach
    ! using the inverse "lumped" mass matrix as a preconditioner:
    !
    !                          -1
    !          y_k+1 := y_k + L   * (r - M * y_k)
    !
    !                    m
    !           L_ii := sum M_ij         L_ij := 0 for i /= j
    !                   j=1
    !
    ! This simple Defect-Correction-Loop converges quite fast (at least as
    ! long as L is regular, that is ^_^) so there is no need to seek for
    ! solvers as CG or even multigrid here...
    !
    ! Remark:
    ! Basically, we could also discretise the boundary conditions of our
    ! Poisson equation in Q1 and implement them into the mass matrix M.
    ! However, this is not crucial as the mass matrix M is regular even
    ! without any boundary conditions and the violation of the boundary
    ! conditions in our resulting Q1 solution will be tolerable...
    !
    ! So let us create a "lumped" mass matrix first...
    
    ! Create a copy of M
    call lsyssc_duplicateMatrix(rmatrixMassPrj, rlumpedMassPrj,&
                                LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
    
    ! And lump the copy to get L
    call lsyssc_lumpMatrixScalar(rlumpedMassPrj, LSYSSC_LUMP_DIAG)
    
    ! Furthermore we need a defect vector
    call lsyssc_createVecIndMat(rmatrixMassPrj, rvecDefQ1)
    
    ! And let us start the Defect-Correction-Loop
    do i = 1, 50
    
      ! Calculate current defect d_i := r - M * y_i
      call lsyssc_copyVector(rvecRhsQ1, rvecDefQ1)
      call lsyssc_scalarMatVec(rmatrixMassPrj, rvecSolQ1, rvecDefQ1,&
                               -1.0_DP, 1.0_DP)
      
      ! Calculate the current residual = || d_i ||_2
      dres = lsyssc_vectorNorm(rvecDefQ1, LINALG_NORML2)
      
      ! Print the residual to screen
      call output_line('L2prj: Iteration ' // trim(sys_siL(i,10)) // &
                       ',  !!RES!! = ' // trim(sys_sdEL(dres,15)))
      
      ! Is our L2-projection good enough?
      if (dres .le. 1e-7_DP) exit
      
      ! Otherwise multiply the defect by the inverse of L
      call lsyssc_invertedDiagMatVec(rlumpedMassPrj, rvecDefQ1, 1.0_DP,&
                                     rvecDefQ1)
      
      ! And add it onto y_i
      call lsyssc_vectorLinearComb(rvecDefQ1, rvecSolQ1, 1.0_DP, 1.0_DP)
    
    end do

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! L2-projection of solution into Q1 space done
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Postprocessing of the solution
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'

    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,&
        Rlevels(NLMAX)%rtriangulation,trim(sucddir)//'/u2d_1_l2prj.vtk')
    
    ! Add our Q1-solution to the UCD exporter:
    call lsyssc_getbase_double (rvecSolQ1,p_Ddata)
    call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Clean up
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.

    ! Release the three Q1 vectors
    call lsyssc_releaseVector(rvecDefQ1)
    call lsyssc_releaseVector(rvecRhsQ1)
    call lsyssc_releaseVector(rvecSolQ1)
    
    ! Release the Q1 mass matrix and its lumped version
    call lsyssc_releaseMatrix(rlumpedMassPrj)
    call lsyssc_releaseMatrix(rmatrixMassPrj)
    
    ! And release the Q1 discretisation
    call spdiscr_releaseDiscr(rdiscrQ1)
    
    ! That was all we used in the L2-projection of the solution and
    ! that has not been already released before.
    
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)

    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the L2-projection
    do i = NLMAX, NLMIN+1, -1
    
      ! Release the projection structure itself
      call mlprj_doneProjection(Rlevels(i)%rprojection)
      
      ! Release the 2-Level-Mass matrix
      call lsyssc_releaseMatrix (Rlevels(i)%rmatrix2Lvl)
      
      ! Release the mass matrix
      call lsyssc_releaseMatrix (Rlevels(i)%rmatrixMass)

    end do
    
    ! Release the projection structure on the coarse mesh
    call mlprj_doneProjection(Rlevels(NLMIN)%rprojection)
        
    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rtempBlock)
    call lsysbl_releaseVector (rvectorBlock)
    call lsysbl_releaseVector (rrhsBlock)
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
