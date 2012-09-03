!##############################################################################
!# ****************************************************************************
!# <name> codire_method5 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a simple
!# Convection-Diffusion-Reaction problem with constant coefficients
!# on a simple domain.
!#
!# The routine splits up the tasks of reading the domain, creating
!# triangulations, discretisation, solving, postprocessing and creanup into
!# different subroutines. The communication between these subroutines
!# is done using an application-specific structure saving problem data
!# as well as a collection structure for the communication with callback
!# routines.
!#
!# The routines here behave similar to codire_method3. In difference,
!# a multigrid-solver with ILU(0) smoother and BiCGStab coarse grid solver
!# is used and matrices/vectors are sorted for Cuthill-McKee.
!# </purpose>
!##############################################################################

module codire_method5

  use fsystem
  use genoutput
  use storage
  use boundary
  use cubature
  use derivatives
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use sortstrategybase
  use sortstrategy
  use triangulation
  use element
  use spatialdiscretisation
  use coarsegridcorrection
  use filtersupport
  use linearsystemscalar
  use linearsystemblock
  use scalarpde
  use bilinearformevaluation
  use linearformevaluation
  use multilevelprojection
  use linearsolver
  use discretebc
  use ucd
  
  use collection
  use paramlist
    
  use codire_callback
  
  implicit none
  
!<types>

!<typeblock description="Type block defining all information about one level">

  type t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo) :: rcubatureInfo
    
    ! A system matrix for that specific level. The matrix will receive the
    ! discrete Laplace operator.
    type(t_matrixBlock) :: rmatrix

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC
  
    ! Sorting strategy for resorting vectors/matrices.
    type(t_blockSortStrategy) :: rsortStrategy

    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    type(t_filterChain), dimension(1) :: RfilterChain
    
    ! Number of filters in the filter chain
    integer :: nfilters

  end type
  
!</typeblock>


!<typeblock description="Application-specific type block for codire problem">

  type t_problem
  
    ! Minimum refinement level; = Level i in RlevelInfo
    integer :: ilvmin
    
    ! Maximum refinement level
    integer :: ilvmax

    ! An object for saving the domain:
    type(t_boundary) :: rboundary

    ! A solution vector and a RHS vector on the finest level.
    type(t_vectorBlock) :: rvector,rrhs

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation.
    type(t_problem_lvl), dimension(:), pointer :: RlevelInfo
    
    ! A collection object that saves structural data and some
    ! problem-dependent information which is e.g. passed to
    ! callback routines.
    type(t_collection) :: rcollection

  end type

!</typeblock>

!</types>
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine cdrm5_initParamTriang (ilvmin,ilvmax,rproblem,sPRMfile,sTRIfile)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the parameter ilv.
!</description>

!<input>
  ! Minimum refinement level of the mesh; = coarse grid = level 1
  integer, intent(in) :: ilvmin
  
  ! Maximum refinement level
  integer, intent(in) :: ilvmax

  ! the PRM and TRI file of the mesh
  character(len=*), intent(in) :: sPRMfile,sTRIfile
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  
    ! Initialise the level in the problem structure
    rproblem%ilvmin = ilvmin
    rproblem%ilvmax = ilvmax
    allocate(rproblem%RlevelInfo(ilvmin:ilvmax))

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rproblem%rboundary, sPRMfile)
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation, &
        sTRIfile, rproblem%rboundary)
    
    ! Refine the mesh up to the minimum level
    call tria_quickRefine2LevelOrdering(rproblem%ilvmin-1,&
        rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation,rproblem%rboundary)
    
    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    call tria_initStandardMeshFromRaw (&
        rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation,rproblem%rboundary)
    
    ! Now, refine to level up to nlmax.
    do i=rproblem%ilvmin+1,rproblem%ilvmax
      call tria_refine2LevelOrdering (rproblem%RlevelInfo(i-1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)
      call tria_initStandardMeshFromRaw (rproblem%RlevelInfo(i)%rtriangulation,&
          rproblem%rboundary)
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cdrm5_initDiscretisation (rproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: I
  
    ! An object for saving the domain:
    type(t_boundary), pointer :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation

    ! An object for the block discretisation on one level
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! A cubature information structure on one level
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo
    
    do i=rproblem%ilvmin,rproblem%ilvmax

      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rboundary => rproblem%rboundary
      p_rtriangulation => rproblem%RlevelInfo(i)%rtriangulation
      
      ! Get the cubature information structure
      p_rcubatureInfo => rproblem%RlevelInfo(i)%rcubatureInfo
      
      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies the blocks in the
      ! solution vector. In this simple problem, we only have one block.
      allocate(p_rdiscretisation)
      call spdiscr_initBlockDiscr (p_rdiscretisation,1,&
                                   p_rtriangulation, p_rboundary)

      ! Save the discretisation structure to our local LevelInfo structure
      ! for later use.
      rproblem%RlevelInfo(i)%p_rdiscretisation => p_rdiscretisation

      ! p_rdiscretisation%Rdiscretisations is a list of scalar
      ! discretisation structures for every component of the solution vector.
      ! Initialise the first element of the list to specify the element
      ! and cubature rule for this solution component:
      call spdiscr_initDiscr_simple ( &
                  p_rdiscretisation%RspatialDiscr(1), &
                  EL_Q1,p_rtriangulation, p_rboundary)

      ! Set up an cubature info structure to tell the code which cubature
      ! formula to use
      ! Create an assembly information structure which tells the code
      ! the cubature formula to use. Standard: Gauss 3x3.
      call spdiscr_createDefCubStructure(&  
          p_rdiscretisation%RspatialDiscr(1),rproblem%RlevelInfo(i)%rcubatureInfo,&
          CUB_GEN_AUTO_G3)

    end do
                                   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cdrm5_initMatVec (rproblem,rparams)
  
!<description>
  ! Calculates the system matrix and RHS vector of the linear system
  ! by discretising the problem with the default discretisation structure
  ! in the problem structure.
  ! Sets up a solution vector for the linear system.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem

  ! A parameter list with informations from the DAT file.
  type(t_parlist), intent(in) :: rparams
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A pointer to the system matrix and the RHS/solution vectors.
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! A cubature information structure
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    ! Parameters from the DAT file
    real(DP) :: alpha11,alpha12,alpha21,alpha22,beta1,beta2,gamma
    character(LEN=10) :: Sstr
  
  
    ! And now to the entries of the matrix. For assembling of the entries,
    ! we need a bilinear form, which first has to be set up manually.
    ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
    ! scalar system matrix in 2D.
    
    rform%itermCount = 7
    
    ! alpha * Laplace(u)
    rform%Idescriptors(1,1) = DER_DERIV_X
    rform%Idescriptors(2,1) = DER_DERIV_X
    
    rform%Idescriptors(1,2) = DER_DERIV_Y
    rform%Idescriptors(2,2) = DER_DERIV_X
    
    rform%Idescriptors(1,3) = DER_DERIV_X
    rform%Idescriptors(2,3) = DER_DERIV_Y
    
    rform%Idescriptors(1,4) = DER_DERIV_Y
    rform%Idescriptors(2,4) = DER_DERIV_Y
    
    ! (beta1, beta2)^T * grad(u)
    rform%Idescriptors(1,5) = DER_DERIV_X
    rform%Idescriptors(2,5) = DER_FUNC
    
    rform%Idescriptors(1,6) = DER_DERIV_Y
    rform%Idescriptors(2,6) = DER_FUNC
    
    ! gamma * u
    rform%Idescriptors(1,7) = DER_FUNC
    rform%Idescriptors(2,7) = DER_FUNC

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff = .true.
    
    ! get the coefficients from the parameter list
    call parlst_getvalue_string (rparams, "EQUATION", "ALPHA11", Sstr, "1.0")
    read(Sstr,*) alpha11
    call parlst_getvalue_string (rparams, "EQUATION", "ALPHA12", Sstr, "0.0")
    read(Sstr,*) alpha12
    call parlst_getvalue_string (rparams, "EQUATION", "ALPHA21", Sstr, "0.0")
    read(Sstr,*) alpha21
    call parlst_getvalue_string (rparams, "EQUATION", "ALPHA22", Sstr, "1.0")
    read(Sstr,*) alpha22
    call parlst_getvalue_string (rparams, "EQUATION", "BETA1", Sstr, "0.0")
    read(Sstr,*) beta1
    call parlst_getvalue_string (rparams, "EQUATION", "BETA2", Sstr, "0.0")
    read(Sstr,*) beta2
    call parlst_getvalue_string (rparams, "EQUATION", "GAMMA", Sstr, "0.0")
    read(Sstr,*) gamma
    
    rform%Dcoefficients(1)  = alpha11
    rform%Dcoefficients(2)  = alpha12
    rform%Dcoefficients(3)  = alpha21
    rform%Dcoefficients(4)  = alpha22
    rform%Dcoefficients(5)  = beta1
    rform%Dcoefficients(6)  = beta2
    rform%Dcoefficients(7)  = gamma

    do i = rproblem%ilvmin, rproblem%ilvmax
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
      ! Get a pointer to the matrix
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! Get the cubature information structure
      p_rcubatureInfo => rproblem%RlevelInfo(i)%rcubatureInfo
            
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)
      
      ! Save matrix to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      call collct_setvalue_mat(rproblem%rcollection,"LAPLACE",p_rmatrix,.true.,i)

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
                p_rmatrix%RmatrixBlock(1,1))
    
      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_CoDiRe for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
      ! the framework will call the callback routine to get analytical data.
      !
      ! We pass our collection structure as well to this routine,
      ! so the callback routine has access to everything what is
      ! in the collection.
      call bilf_buildMatrixScalar (rform,.true.,&
           p_rmatrix%RmatrixBlock(1,1),p_rcubatureInfo,coeff_CoDiRe,rproblem%rcollection)
                                   
      ! Create a sort strategy structure for our discretisation
      call sstrat_initBlockSorting (rproblem%RlevelInfo(i)%rsortStrategy,p_rdiscretisation)

      ! Calculate the resorting strategy.
      call sstrat_initCuthillMcKee (rproblem%RlevelInfo(i)%rsortStrategy%p_Rstrategies(1),&
          p_rmatrix%RmatrixBlock(1,1))
      
      ! Attach the sorting strategy to the matrix. The matrix is not yet sorted.
      call lsysbl_setSortStrategy (p_rmatrix,&
          rproblem%RlevelInfo(i)%rsortStrategy,&
          rproblem%RlevelInfo(i)%rsortStrategy)
      
      ! Resort the matrix.
      call lsysbl_sortMatrix (p_rmatrix,.true.)
      
    end do

    ! (Only) on the finest level, we need to calculate a RHS vector
    ! and to allocate a solution vector.
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector

    ! Next step: Create a RHS vector and a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVectorBlock (p_rdiscretisation,p_rrhs,.true.)
    call lsysbl_createVectorBlock (p_rdiscretisation,p_rvector,.true.)

    ! Save the solution/RHS vector to the collection. Might be used
    ! later (e.g. in nonlinear problems)
    call collct_setvalue_vec(rproblem%rcollection,"RHS",p_rrhs,.true.)
    call collct_setvalue_vec(rproblem%rcollection,"SOLUTION",p_rvector,.true.)

    ! The vector structure is ready but the entries are missing.
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to the first subvector of
    ! the block vector using the discretisation structure of the
    ! first block.
    !
    ! We pass our collection structure as well to this routine,
    ! so the callback routine has access to everything what is
    ! in the collection.
    call linf_buildVectorScalar (&
              rlinform,.true.,p_rrhs%RvectorBlock(1),p_rcubatureInfo,&
              coeff_RHS,rproblem%rcollection)
                                
    ! Clear the solution vector on the finest level.
    call lsysbl_clearVector(rproblem%rvector)
    
    ! Install the resorting strategy in the RHS- and the solution
    ! vector, but do not resort them yet!
    ! We resort the vectors just before solving.
    call lsysbl_setSortStrategy (p_rrhs,&
        rproblem%RlevelInfo(rproblem%ilvmax)%rsortStrategy)
    call lsysbl_setSortStrategy (p_rvector,&
        rproblem%RlevelInfo(rproblem%ilvmax)%rsortStrategy)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cdrm5_initDiscreteBC (rproblem)
  
!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

  ! A pointer to the system matrix and the RHS vector as well as
  ! the discretisation
  type(t_matrixBlock), pointer :: p_rmatrix
  type(t_blockDiscretisation), pointer :: p_rdiscretisation

  type(t_boundaryRegion) :: rboundaryRegion

    do i=rproblem%ilvmin,rproblem%ilvmax
    
      ! Get our matrix from the problem structure.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! From the matrix or the RHS we have access to the discretisation
      ! boundary conditions.
      p_rdiscretisation => p_rmatrix%p_rblockDiscrTest
      
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions.
      call bcasm_initDiscreteBC(rproblem%RlevelInfo(i)%rdiscreteBC)
      !
      ! We "know" already (from the problem definition) that we have four boundary
      ! segments in the domain. Each of these, we want to use for enforcing
      ! some kind of boundary condition.
      !
      ! We ask the boundary routines to create a "boundary region" - which is
      ! simply a part of the boundary corresponding to a boundary segment.
      ! A boundary region roughly contains the type, the min/max parameter value
      ! and whether the endpoints are inside the region or not.
      call boundary_createRegion(rproblem%rboundary,1,1,rboundaryRegion)
      
      ! We use this boundary region and specify that we want to have Dirichlet
      ! boundary there. The following call does the following:
      ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
      !   We specify icomponent="1" to indicate that we set up the
      !   Dirichlet BC`s for the first (here: one and only) component in the
      !   solution vector.
      ! - Discretise the boundary condition so that the BC`s can be applied
      !   to matrices and vectors
      ! - Add the calculated discrete BC`s to rdiscreteBC for later use.
      call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
          rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
          getBoundaryValues)
                               
      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rproblem%rboundary,1,2,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
          rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
          getBoundaryValues)
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rproblem%rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
          rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
          getBoundaryValues)
      
      ! Edge 4 of boundary component 1. That is it.
      call boundary_createRegion(rproblem%rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
          rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
          getBoundaryValues)

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cdrm5_implementBC (rproblem)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,ilvmax
  type(t_discreteBC), pointer :: p_rdiscreteBC
  
    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector
    
    ! Get our the right hand side and solution from the problem structure
    ! on the finest level
    ilvmax = rproblem%ilvmax
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector

    ! Get the BCs    
    p_rdiscreteBC => rproblem%RlevelInfo(ilvmax)%rdiscreteBC

    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    call vecfil_discreteBCrhs (p_rrhs,p_rdiscreteBC)
    call vecfil_discreteBCsol (p_rvector,p_rdiscreteBC)

    ! Implement discrete boundary conditions into the matrices on all
    ! levels, too. Call the appropriate matrix filter to modify
    ! all matrices according to the attached discrete boundary conditions.
    do i=rproblem%ilvmin,rproblem%ilvmax
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      p_rdiscreteBC => rproblem%RlevelInfo(i)%rdiscreteBC
      call matfil_discreteBC (p_rmatrix,p_rdiscreteBC)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cdrm5_solve (rproblem)
  
!<description>
  ! Solves the given problem by applying a linear solver.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
    integer :: ilvmin,ilvmax
    integer :: i

    ! Error indicator during initialisation of the solver
    integer :: ierror
  
    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector
    type(t_vectorBlock), target :: rvecTmp

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rsmoother
    type(t_linsolNode), pointer :: p_rcoarseGridSolver,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices
    
    ! One level of multigrid
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    
    ilvmin = rproblem%ilvmin
    ilvmax = rproblem%ilvmax
    
    ! Get our right hand side / solution / matrix on the finest
    ! level from the problem structure.
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector
    p_rmatrix => rproblem%RlevelInfo(ilvmax)%rmatrix
    
    ! Create a temporary vector we need that for some preparation.
    call lsysbl_createVectorBlock (p_rrhs, rvecTmp, .false.)
    
    ! Resort the RHS and solution vector according to the resorting
    ! strategy given in the matrix.
    !
    ! Use the temporary vector from above to store intermediate data.
    ! The vectors are assumed to know how they are resorted (the strategy
    ! is already attached to them). So call the resorting routines
    ! to resort them as necessary!
    ! We use the first subvector of rvecTmp as temporary data; it is
    ! large enough, as we only have one block.
    call lsysbl_synchroniseSort(p_rmatrix,p_rrhs,rvecTmp%RvectorBlock(1))
    call lsysbl_synchroniseSort(p_rmatrix,p_rvector,rvecTmp%RvectorBlock(1))
    
    ! Now we have to build up the level information for multigrid.
    !
    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    call linsol_initMultigrid2 (p_rsolverNode,ilvmax-ilvmin+1)
    
    ! Then set up smoothers / coarse grid solver:
    do i=ilvmin,ilvmax
    
      ! Set up a filter chain for implementing boundary conditions on that level
      call filter_clearFilterChain (rproblem%RlevelInfo(i)%RfilterChain,&
          rproblem%RlevelInfo(i)%nfilters)
      call filter_newFilterDiscBCDef (rproblem%RlevelInfo(i)%RfilterChain,&
          rproblem%RlevelInfo(i)%nfilters,rproblem%RlevelInfo(i)%rdiscreteBC)
      
      ! On the coarsest grid, set up a coarse grid solver and no smoother
      ! On finer grids, set up a smoother but no coarse grid solver.
      nullify(p_rpreconditioner)
      nullify(p_rsmoother)
      nullify(p_rcoarseGridSolver)
      
      ! Get the level
      call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
      
      if (i .eq. ilvmin) then
        ! Set up a BiCGStab solver with ILU preconditioning as coarse grid solver
        ! would be:
        ! CALL linsol_initMILUs1x1 (p_rpreconditioner,0,0.0_DP)
        ! CALL linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,RfilterChain)
        
        ! Set up UMFPACK coarse grid solver.
        call linsol_initUMFPACK4 (p_rcoarseGridSolver)
        
      else
        ! Setting up Jacobi smoother for multigrid would be:
        ! CALL linsol_initJacobi (p_rsmoother)
        
        ! Setting up VANCA smoother for multigrid would be:
        ! CALL linsol_initVANKA (p_rsmoother)

        ! Set up an ILU smoother for multigrid with damping parameter 0.7,
        ! 4 smoothing steps:
        call linsol_initMILUs1x1 (p_rsmoother,0,0.0_DP)
        call linsol_convertToSmoother (p_rsmoother,4,0.7_DP)
      end if
    
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (p_rsolverNode,i-ilvmin+1,p_rlevelInfo)
      p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
      p_rlevelInfo%p_rpresmoother => p_rsmoother
      p_rlevelInfo%p_rpostsmoother => p_rsmoother

      ! Attach the filter chain which imposes boundary conditions on that level.
      p_rlevelInfo%p_RfilterChain => rproblem%RlevelInfo(i)%RfilterChain
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
    allocate(Rmatrices(ilvmin:ilvmax))
    do i=ilvmin,ilvmax
      call lsysbl_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrix,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    
    call linsol_setMatrices(p_RsolverNode,Rmatrices(ilvmin:ilvmax))
    
    ! We can release Rmatrices immediately -- as long as we do not
    ! release rproblem%RlevelInfo(i)%rmatrix!
    do i=ilvmin,ilvmax
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
    
    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    call linsol_solveAdaptively (p_rsolverNode,p_rvector,p_rrhs,rvecTmp)
    
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Unsort the vectors again in case they were resorted before calling
    ! the solver.
    ! We use the first subvector of rvecTmp as temporary data; it is
    ! large enough, as we only have one block.
    call lsysbl_sortVector (p_rrhs,.false.,rvecTmp%RvectorBlock(1))
    call lsysbl_sortVector (p_rvector,.false.,rvecTmp%RvectorBlock(1))
    
    ! Release the temporary vector
    call lsysbl_releaseVector (rvecTmp)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cdrm5_postprocessing (rproblem,sucddir)
  
!<description>
  ! Writes the solution into a VTK file.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem

  ! Path where to write output data to
  character(len=*), intent(in) :: sucddir
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! We need some more variables for postprocessing
    real(DP), dimension(:), pointer :: p_Ddata
    
    ! Output block for UCD output to VTK file
    type(t_ucdExport) :: rexport

    ! A pointer to the solution vector and to the triangulation.
    type(t_vectorBlock), pointer :: p_rvector
    type(t_triangulation), pointer :: p_rtriangulation

    ! Get the solution vector from the problem structure.
    p_rvector => rproblem%rvector
    
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      p_rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
    
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing.

    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
        trim(sucddir)//"/u5.vtk")
    
    call lsyssc_getbase_double (p_rvector%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,"sol",UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cdrm5_doneMatVec (rproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    integer :: i

    ! Delete solution/RHS vector
    call lsysbl_releaseVector (rproblem%rvector)
    call lsysbl_releaseVector (rproblem%rrhs)

    ! Delete the variables from the collection.
    call collct_deletevalue (rproblem%rcollection,"RHS")
    call collct_deletevalue (rproblem%rcollection,"SOLUTION")

    ! Release matrices and vectors on all levels
    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Delete the matrix
      call lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rmatrix)
      
      ! Delete the variables from the collection.
      call collct_deletevalue (rproblem%rcollection,"LAPLACE",i)

      ! Release the sorting strategy
      call sstrat_doneBlockSorting (rproblem%RlevelInfo(i)%rsortStrategy)

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cdrm5_doneBC (rproblem)
  
!<description>
  ! Releases discrete and analytic boundary conditions from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Release our discrete version of the boundary conditions
      call bcasm_releaseDiscreteBC (rproblem%RlevelInfo(i)%rdiscreteBC)
    end do
    
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine cdrm5_doneDiscretisation (rproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Release the cubature info structure.
      call spdiscr_releaseCubStructure(rproblem%RlevelInfo(i)%rcubatureInfo)

      ! Delete the block discretisation together with the associated
      ! scalar spatial discretisations....
      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%p_rdiscretisation)
      
      ! and remove the allocated block discretisation structure from the heap.
      deallocate(rproblem%RlevelInfo(i)%p_rdiscretisation)
    end do
    
  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine cdrm5_doneParamTriang (rproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Release the triangulation
      call tria_done (rproblem%RlevelInfo(i)%rtriangulation)
    end do
    
    deallocate(rproblem%RlevelInfo)
    
    ! Finally release the domain.
    call boundary_release (rproblem%rboundary)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine codire5
  
!<description>
  ! This is a "separated" CoDiRe solver for solving a convection-diffusion-
  ! reaction problem. The different tasks of the problem are separated into
  ! subroutines. The problem uses a problem-specific structure for the
  ! communication: All subroutines add their generated information to the
  ! structure, so that the other subroutines can work with them.
  ! (THis is somehow a cleaner implementation than using a collection!).
  ! For the communication to callback routines of black-box subroutines
  ! (matrix-assembly), a collection is used.
  !
  ! The following tasks are performed by the subroutines:
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

    ! A paramlist structure with parameters from the dat file
    type(t_parlist) :: rparams

    ! NLMIN receives the minimal level where to discretise for supporting
    ! the solution process.
    ! NLMAX receives the level where we want to solve
    integer :: NLMIN, NLMAX
    
    ! PRM/TRI file
    character(LEN=SYS_STRLEN) :: sfilePRM,sfileTRI

    ! A problem structure for our problem
    type(t_problem), pointer :: p_rproblem
    
    ! Strings that get data about the output directories for pre/postprocessing files
    character(len=SYS_STRLEN) :: sucddir,sstring,smaster
    
    integer :: i
    
    ! Allocate the problem structure on the heap -- it is rather large...
    ! too large for the stack!
    allocate(p_rproblem)
    
    ! Ok, let us start.
    ! Minimum level in Multigrid
    NLMIN = 1
    
    ! Initialise the parameter list
    call parlst_init(rparams)
    
    ! Read the parameters from disc and put a reference to it
    ! to the collection
    call sys_getcommandLineArg(1,smaster,sdefault="./data/codire.dat")
    call parlst_readfromfile (rparams, smaster)

    ! We want to solve our Laplace problem on level...
    call parlst_getvalue_int (rparams, "GENERAL", "NLMAX", NLMAX, 7)

    ! PRM file
    call parlst_getvalue_string (rparams, "GENERAL", &
                                 "sfilePRM", sstring)
    read(sstring,*) sfilePRM
                                 
    ! TRI file
    call parlst_getvalue_string (rparams, "GENERAL", &
                                 "sfileTRI", sstring)
    read(sstring,*) sfileTRI

    ! Initialise the collection
    call collct_init (p_rproblem%rcollection)
    do i=1,NLMAX
      call collct_addlevel (p_rproblem%rcollection)
    end do

    call collct_setvalue_parlst (p_rproblem%rcollection, "PARAMS", rparams, .true.)
    
    ! Get the path where to write VTK`s to.
    call parlst_getvalue_string (rparams, "GENERAL", &
                                 "sucddir", sstring)
    read(sstring,*) sucddir

    ! So now the different steps - one after the other.
    !
    ! Initialisation
    call cdrm5_initParamTriang (NLMIN,NLMAX,p_rproblem,sfilePRM,sfileTRI)
    call cdrm5_initDiscretisation (p_rproblem)
    call cdrm5_initMatVec (p_rproblem,rparams)
    call cdrm5_initDiscreteBC (p_rproblem)
    
    ! Implementation of boundary conditions
    call cdrm5_implementBC (p_rproblem)
    
    ! Solve the problem
    call cdrm5_solve (p_rproblem)
    
    ! Postprocessing
    call cdrm5_postprocessing (p_rproblem,sucddir)
    
    ! Cleanup
    call cdrm5_doneMatVec (p_rproblem)
    call cdrm5_doneBC (p_rproblem)
    call cdrm5_doneDiscretisation (p_rproblem)
    call cdrm5_doneParamTriang (p_rproblem)
    
    ! Release parameter list
    call collct_deletevalue (p_rproblem%rcollection,"PARAMS")
    call parlst_done (rparams)

    ! Print some statistical data about the collection - anything forgotten?
    print *
    print *,"Remaining collection statistics:"
    print *,"--------------------------------"
    print *
    call collct_printStatistics (p_rproblem%rcollection)
    
    ! Finally release the collection.
    call collct_done (p_rproblem%rcollection)
    deallocate(p_rproblem)
    
  end subroutine

end module
